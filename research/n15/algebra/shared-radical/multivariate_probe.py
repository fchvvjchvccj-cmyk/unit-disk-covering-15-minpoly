"""实测二次塔在多个符号参数下的规范形规模，禁止先做全局范数展开。"""

import argparse
from pathlib import Path
import subprocess
import sys
import time

from bridge_support import HERE, write, read, sha, bounded, binary_environment
import tower


class MultivariateRational(tower.Rational):
    context = None
    term_limit = None

    def __init__(self, numerator, denominator=None, cancel=True):
        n = self.context.constant(numerator%self.prime) if isinstance(numerator, int) else numerator
        d = self.one_poly if denominator is None else denominator
        if d.is_zero():
            raise ZeroDivisionError("多参数有理函数分母为零")
        if n.is_zero():
            d = self.one_poly
        elif cancel and not d.is_one():
            common = n.gcd(d)
            if not common.is_one():
                n, d = n//common, d//common
        leading = int(d.leading_coefficient())
        if leading != 1:
            inverse = pow(leading, -1, self.prime)
            n, d = n*inverse, d*inverse
        self.n, self.d = n, d
        degree, terms = int(max(n.total_degree(), d.total_degree())), max(len(n), len(d))
        self.metrics["max_coefficient_degree"] = max(self.metrics["max_coefficient_degree"], degree)
        self.metrics["max_coefficient_terms"] = max(self.metrics["max_coefficient_terms"], terms)
        if terms > self.term_limit:
            raise tower.TowerStop("单个系数的多项式项数超过上限")
        if time.perf_counter() > self.deadline:
            raise tower.TowerStop("多参数塔内约分达到墙钟上限")

    def encode(self):
        pack = lambda poly: [[list(map(int, exponent)), int(coefficient)] for exponent, coefficient in poly.to_dict().items()]
        return {"n": pack(self.n), "d": pack(self.d)}


def worker(args):
    from flint import nmod_mpoly_ctx
    target = Path(args.worker)
    tick = time.perf_counter()
    names = args.parameters
    context = nmod_mpoly_ctx.get(names, args.prime, ordering="lex")
    M = MultivariateRational
    M.context, M.prime = context, args.prime
    M.zero_poly, M.one_poly = context.constant(0), context.constant(1)
    M.term_limit, M.deadline = args.max_terms, tick+args.seconds
    M.metrics = {"max_coefficient_degree": 0, "max_coefficient_terms": 0, "rational_additions": 0, "rational_multiplications": 0}
    tower.Rational = M
    from quotient_probe import model
    algebra = tower.Algebra()
    events = []

    def checkpoint(stage):
        event = {"stage": stage, "levels": len(algebra.relations), "seconds": time.perf_counter()-tick, **M.metrics}
        events.append(event)
        print(event, flush=True)

    report = {"prime": args.prime, "symbolic_parameters": names, "derivatives_computed": args.derivatives,
              "driver_sha256": sha(__file__), "algebra_sha256": sha(HERE / "tower.py"),
              "model_sha256": sha(HERE / "quotient_probe.py")}
    try:
        fixed = {"u": 2, "v": 3, "w": 1, "r": 2}
        values = [M(context.gen(names.index(name))) if name in names else M(fixed[name]) for name in ("u", "v", "w", "r")]
        H, _ = model(algebra, values, args.derivatives, checkpoint)
        vectors = [H.value]+(list(H.gradient) if args.derivatives else [])
        total = sum(len(c.n)+len(c.d) for value in vectors for c in value.promote(5).coefficients)
        report.update(status="completed", output_coefficient_terms=total)
        if total <= args.save_terms:
            write(target / "quotient.json.gz", {**report,
                  "relations": [{"a": [c.encode() for c in a], "b": [c.encode() for c in b]} for a, b in algebra.relations],
                  "H": H.value.encode(5), "gradient_u_v_w_r": [g.encode(5) for g in H.gradient] if args.derivatives else None})
        else:
            report["coefficient_output_saved"] = False
    except (tower.TowerStop, ZeroDivisionError) as error:
        report.update(status="stopped", reason=str(error))
    report.update(events=events, seconds=time.perf_counter()-tick, **M.metrics)
    write(target / "worker.json", report)
    print({k: v for k, v in report.items() if k != "events"}, flush=True)


def run(args):
    target = HERE / "multivariate-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有多参数规范形探测")
    target.mkdir(parents=True)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(target), "--parameters", *args.parameters,
               "--prime", str(args.prime), "--seconds", str(args.seconds), "--max-terms", str(args.max_terms), "--save-terms", str(args.save_terms)]
    if args.derivatives:
        command.append("--derivatives")
    with bounded() as snapshot:
        tick = time.perf_counter()
        with (target / "worker.log").open("x", encoding="utf-8") as stream:
            try:
                child = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=binary_environment(), timeout=args.seconds+5, check=False)
                status, returncode = "worker_completed", child.returncode
            except subprocess.TimeoutExpired:
                status, returncode = "hard_timeout", None
        report = {"controller_status": status, "returncode": returncode, "seconds": time.perf_counter()-tick, "limits": snapshot()}
        if (target / "worker.json").exists():
            report["worker"] = read(target / "worker.json")
        write(target / "controller.json", report)
        print({k: v for k, v in report.items() if k not in ("limits", "worker")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag")
    parser.add_argument("--worker")
    parser.add_argument("--parameters", nargs="+", default=["u", "r"])
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--seconds", type=int, default=60)
    parser.add_argument("--max-terms", type=int, default=250000)
    parser.add_argument("--save-terms", type=int, default=200000)
    parser.add_argument("--derivatives", action="store_true")
    options = parser.parse_args()
    worker(options) if options.worker else run(options)
