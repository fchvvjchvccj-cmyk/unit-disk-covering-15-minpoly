"""受限展开九变量有理电路，实测有限域消元输入的规模。"""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from support import HERE, read, write, sha, bounded, binary_environment


class ExpansionStop(Exception):
    pass


class Rational:
    context = None
    prime = None
    term_limit = None
    pair_limit = None
    deadline = None
    metrics = None

    def __init__(self, n, d=None, reduce=True):
        n = self.context.constant(n % self.prime) if isinstance(n, int) else n
        d = self.context.constant(1) if d is None else d
        if d.is_zero():
            raise ZeroDivisionError("展开出现恒零分母")
        if reduce and not n.is_zero() and not d.is_one():
            common = n.gcd(d)
            if not common.is_one():
                n, d = n//common, d//common
        if n.is_zero():
            d = self.context.constant(1)
        leading = int(d.leading_coefficient())
        if leading != 1:
            inverse = pow(leading, -1, self.prime)
            n, d = n*inverse, d*inverse
        self.n, self.d = n, d
        count = max(len(n), len(d))
        self.metrics["largest_reduced_polynomial_terms"] = max(count, self.metrics["largest_reduced_polynomial_terms"])
        if count > self.term_limit:
            raise ExpansionStop("约分后单个多项式项数超过上限")
        if time.perf_counter() > self.deadline:
            raise ExpansionStop("在算术指令间达到展开墙钟上限")

    def cast(self, other):
        return other if isinstance(other, Rational) else Rational(other)

    @classmethod
    def product(cls, a, b):
        estimate = len(a)*len(b)
        cls.metrics["largest_product_pair_bound"] = max(estimate, cls.metrics["largest_product_pair_bound"])
        if estimate > cls.pair_limit:
            raise ExpansionStop("乘积的单项式对上界超过配置值，未执行该乘法")
        return a*b

    def __neg__(self):
        return Rational(-self.n, self.d, reduce=False)

    def __add__(self, other):
        other = self.cast(other)
        if self.d == other.d:
            return Rational(self.n+other.n, self.d)
        common = self.d.gcd(other.d)
        a, b = self.d//common, other.d//common
        numerator = self.product(self.n, b)+self.product(other.n, a)
        return Rational(numerator, self.product(a, other.d))

    __radd__ = __add__

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        g, h = self.n.gcd(other.d), other.n.gcd(self.d)
        return Rational(self.product(self.n//g, other.n//h), self.product(self.d//h, other.d//g), reduce=False)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        if other.n.is_zero():
            raise ZeroDivisionError("展开出现恒零除数")
        g, h = self.n.gcd(other.n), other.d.gcd(self.d)
        return Rational(self.product(self.n//g, other.d//h), self.product(self.d//h, other.n//g), reduce=False)

    def __rtruediv__(self, other):
        return self.cast(other)/self


def polynomial_data(p):
    return {"terms": len(p), "total_degree": int(p.total_degree()),
            "coefficients": [[list(e), int(c)] for e, c in p.to_dict().items()]}


def worker(args):
    from flint import nmod_mpoly_ctx
    target = Path(args.worker)
    source = HERE / "system-runs" / args.system / "circuit.json.gz"
    data = read(source)
    started = time.perf_counter()
    Rational.context = nmod_mpoly_ctx.get(data["variables"], args.prime, ordering="lex")
    Rational.prime, Rational.term_limit = args.prime, args.max_terms
    Rational.pair_limit = args.pair_bound
    Rational.deadline = started+args.seconds
    Rational.metrics = {"largest_reduced_polynomial_terms": 0, "largest_product_pair_bound": 0}
    equations = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    guards = sorted(set(data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+
                        data["simple_root_derivative_nodes"]+data["critical_locus_nonzero_guards"]))
    requested = set(equations+guards)
    active = set(requested)
    for i in range(len(data["nodes"])-1, -1, -1):
        if i in active:
            active.update(data["nodes"][i].get("args", []))
    uses = [0]*len(data["nodes"])
    for i in active:
        for j in data["nodes"][i].get("args", []):
            uses[j] += 1
    registers = [None]*len(data["nodes"])
    checkpoint, i = [], -1
    state = {"prime": args.prime, "variables": data["variables"], "circuit_sha256": sha(source),
             "active_nodes": len(active), "term_limit": args.max_terms, "product_pair_limit": args.pair_bound,
             "interpretation": "仅测量有理函数的直接稀疏展开；停止不说明最小多项式不存在或次数下界。"}
    try:
        for i, node in enumerate(data["nodes"]):
            if i not in active:
                continue
            op = node["op"]
            operands = [registers[j] for j in node.get("args", [])]
            if op == "var":
                value = Rational(Rational.context.gen(node["value"]))
            elif op == "const":
                value = Rational(node["value"])
            elif op == "add":
                value = operands[0]+operands[1]
            elif op == "neg":
                value = -operands[0]
            elif op == "mul":
                value = operands[0]*operands[1]
            elif op == "div":
                value = operands[0]/operands[1]
            else:
                raise ValueError("未知电路指令")
            registers[i] = value
            if i in requested:
                checkpoint.append({"node": i, "equation_index": equations.index(i) if i in equations else None,
                                   "numerator_terms": len(value.n), "denominator_terms": len(value.d),
                                   "numerator_degree": int(value.n.total_degree()), "denominator_degree": int(value.d.total_degree())})
                print(json.dumps(checkpoint[-1]), flush=True)
            if i % 50 == 0:
                print(json.dumps({"node": i, "seconds": time.perf_counter()-started, **Rational.metrics}), flush=True)
            for j in node.get("args", []):
                uses[j] -= 1
                if uses[j] == 0 and j not in requested:
                    registers[j] = None
        state.update(status="completed", last_node=i)
        export = {"prime": args.prime, "variables": data["variables"], "circuit_sha256": sha(source),
                  "equations": [{"numerator": polynomial_data(registers[n].n), "denominator": polynomial_data(registers[n].d)} for n in equations],
                  "nonzero_guards": [{"node": n, "numerator": polynomial_data(registers[n].n), "denominator": polynomial_data(registers[n].d)} for n in guards]}
        write(target / "polynomials.json.gz", export)
    except ExpansionStop as error:
        state.update(status="stopped", last_node=i, reason=str(error))
    state.update(checkpoints=checkpoint, seconds=time.perf_counter()-started, **Rational.metrics)
    write(target / "worker.json", state)
    print(json.dumps({k: v for k, v in state.items() if k not in ("checkpoints", "interpretation")}), flush=True)


def run(args):
    target = HERE / "expansion-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有展开试验")
    target.mkdir(parents=True)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--system", args.system, "--worker", str(target),
               "--prime", str(args.prime), "--seconds", str(args.seconds), "--max-terms", str(args.max_terms),
               "--pair-bound", str(args.pair_bound)]
    with bounded() as snapshot:
        tick = time.perf_counter()
        with (target / "worker.log").open("x", encoding="utf-8") as log:
            try:
                child = subprocess.run(command, env=binary_environment(), stdout=log, stderr=subprocess.STDOUT,
                                       timeout=args.seconds+5, check=False)
                status = "worker_completed" if child.returncode == 0 else "worker_error"
                returncode = child.returncode
            except subprocess.TimeoutExpired:
                status, returncode = "hard_timeout", None
        report = {"controller_status": status, "worker_returncode": returncode, "seconds": time.perf_counter()-tick,
                  "limits": snapshot(), "driver_sha256": sha(__file__), "log_sha256": sha(target / "worker.log")}
        if (target / "worker.json").exists():
            report["worker"] = read(target / "worker.json")
        write(target / "controller.json", report)
        print({k: v for k, v in report.items() if k not in ("worker", "limits")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--tag")
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--seconds", type=int, default=90)
    parser.add_argument("--max-terms", type=int, default=20000)
    parser.add_argument("--pair-bound", type=int, default=4_000_000)
    parser.add_argument("--worker")
    parsed = parser.parse_args()
    worker(parsed) if parsed.worker else run(parsed)
