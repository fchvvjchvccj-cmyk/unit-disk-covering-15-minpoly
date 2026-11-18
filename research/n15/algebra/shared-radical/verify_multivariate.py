"""独立回放Fp(u,r)闭环电路，并与已接受的Fp(r)截面精确比较。"""

import argparse
import time

from bridge_support import HERE, read, write, sha, bounded
from build_shared import construct
from verify_system import execute
from multivariate_probe import MultivariateRational as M
import tower


def prune(circuit, outputs):
    needed, pending = set(), list(outputs)
    while pending:
        index = pending.pop()
        if index in needed:
            continue
        if not 0 <= index < len(circuit["nodes"]):
            raise ValueError("电路输出越界")
        needed.add(index)
        args = circuit["nodes"][index].get("args", [])
        if any(not 0 <= j < index for j in args):
            raise ValueError("电路存在向前引用")
        pending.extend(args)
    indices = {old: new for new, old in enumerate(sorted(needed))}
    nodes = []
    for old in sorted(needed):
        node = dict(circuit["nodes"][old])
        if "args" in node:
            node["args"] = [indices[j] for j in node["args"]]
        nodes.append(node)
    return {"nodes": nodes}, [indices[index] for index in outputs]


def decode(encoded):
    poly = lambda rows: M.context.from_dict({tuple(exponent): coefficient for exponent, coefficient in rows})
    return M(poly(encoded["n"]), poly(encoded["d"]), cancel=False)


def specialize(encoded, fixed_u, prime):
    from flint import nmod_poly
    result = []
    for side in ("n", "d"):
        rows = encoded[side]
        coefficients = [0]*(1+max((exponent[1] for exponent, _ in rows), default=0))
        for exponent, coefficient in rows:
            coefficients[exponent[1]] = (coefficients[exponent[1]]+coefficient*pow(fixed_u, exponent[0], prime))%prime
        result.append(nmod_poly(coefficients, prime))
    if result[1].is_zero():
        raise ZeroDivisionError("双参数规范形在该截面分母恒零")
    return result


def check_specialization(data, line):
    from flint import nmod_poly
    if line["prime"] != data["prime"] or line["fixed_u_v_w"][1:] != [3, 1]:
        raise ValueError("截面与双参数问题不相容")
    fixed_u, prime = line["fixed_u_v_w"][0], line["prime"]
    pairs = list(zip(data["H"], line["H"]))
    for left, right in zip(data["relations"], line["relations"]):
        for side in ("a", "b"):
            if len(left[side]) != len(right[side]):
                raise ValueError("双参数关系与截面维数不符")
            pairs.extend(zip(left[side], right[side]))
    for left, right in pairs:
        n, d = specialize(left, fixed_u, prime)
        nn, dd = nmod_poly(right["n"], prime), nmod_poly(right["d"], prime)
        if n*dd != nn*d:
            raise ValueError("双参数规范形特化后与独立截面不一致")
    return {"verified": True, "fixed_u": fixed_u, "rational_function_equalities": len(pairs)}


def verify(data, circuit, seconds):
    from flint import nmod_mpoly_ctx
    if data["symbolic_parameters"] != ["u", "r"] or data["derivatives_computed"]:
        raise ValueError("本接受端只核验固定v=3,w=1的双参数闭环值")
    rebuilt = construct()
    for key in ("variables", "selected_sign", "nodes", "outputs"):
        if circuit[key] != rebuilt[key]:
            raise ValueError("双参数回放的电路绑定错误")
    M.context = nmod_mpoly_ctx.get(["u", "r"], data["prime"], ordering="lex")
    M.prime, M.term_limit, M.deadline = data["prime"], 250000, time.perf_counter()+seconds
    M.zero_poly, M.one_poly = M.context.constant(0), M.context.constant(1)
    M.metrics = {"max_coefficient_degree": 0, "max_coefficient_terms": 0, "rational_additions": 0, "rational_multiplications": 0}
    tower.Rational = M
    algebra = tower.Algebra()
    algebra.relations = [(tuple(decode(c) for c in row["a"]), tuple(decode(c) for c in row["b"])) for row in data["relations"]]
    if (len(algebra.relations) != 5 or len(data["H"]) != 32
            or any(len(a) != 2**j or len(b) != 2**j for j, (a, b) in enumerate(algebra.relations))):
        raise ValueError("双参数规范形的二次塔维数错误")
    generators = []
    for j in range(5):
        coefficients = [algebra.zero]*2**(j+1)
        coefficients[2**j] = algebra.one
        generators.append(tower.Element(algebra, j+1, tuple(coefficients)))
    parameters = [algebra.scalar(M(M.context.gen(0))), algebra.scalar(3), algebra.scalar(1), algebra.scalar(M(M.context.gen(1)))]
    outputs = circuit["outputs"]["finite_relations"]+[circuit["outputs"]["critical_equations"][0]]
    small, indices = prune(circuit, outputs)
    registers = execute(small, generators+parameters)
    if any(not registers[index].is_zero() for index in indices[:5]):
        raise ValueError("双参数二次塔没有满足电路关系")
    actual = registers[indices[5]].promote(5).coefficients
    expected = [decode(c) for c in data["H"]]
    if any(a.n != b.n or a.d != b.d for a, b in zip(actual, expected)):
        raise ValueError("双参数闭环规范形与独立电路不同")
    return {"verified": True, "symbolic_finite_relations": 5, "symbolic_closure_coefficients": 32,
            "pruned_circuit_nodes": len(small["nodes"]), "scope": "精确Fp(u,r)闭环规范形；v=3,w=1，未核验双参数全导数或求取双参数范数。"}


def run(args):
    target = HERE / "multivariate-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = read(target / "quotient.json.gz")
        circuit_path = HERE / "system-runs/20261015-selected/circuit.json.gz"
        symbolic = verify(data, read(circuit_path), args.seconds)
        lines = []
        for name in args.lines:
            path = HERE / "quotient-runs" / name / "quotient.json.gz"
            acceptance = read(path.parent / "verified-final.json.gz")
            if not acceptance["verified"] or acceptance["quotient_sha256"] != sha(path):
                raise ValueError("截面未通过绑定接受端")
            result = check_specialization(data, read(path))
            result.update(path=path.relative_to(HERE).as_posix(), sha256=sha(path))
            lines.append(result)
        report = {"verified": True, "symbolic": symbolic, "specializations": lines,
                  "quotient_sha256": sha(target / "quotient.json.gz"), "circuit_sha256": sha(circuit_path),
                  "verifier_sha256": sha(__file__), "seconds": time.perf_counter()-tick, "limits": snapshot()}
        write(target / args.output, report)
        print({"verified": True, "symbolic_coefficients": 32, "exact_specializations": len(lines), "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--seconds", type=int, default=60)
    parser.add_argument("--lines", nargs="+", default=["20261015-jets-p65521"])
    parser.add_argument("--output", default="verified.json.gz")
    run(parser.parse_args())
