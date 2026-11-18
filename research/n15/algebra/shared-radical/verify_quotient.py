"""精确回放Fp(r)电路，并以原几何纤维核验规范形、导数及范数。"""

import argparse
import time

from bridge_support import HERE, AXIS, FACTOR, read, write, sha, bounded
from tower import Rational, Algebra, Element
from verify_system import execute
from verify_shared import check_fiber
from fiber import compute


def setup(prime, seconds=300):
    from flint import nmod_poly
    Rational.prime = prime
    Rational.zero_poly, Rational.one_poly = nmod_poly([], prime), nmod_poly([1], prime)
    Rational.deadline, Rational.degree_limit = time.perf_counter()+seconds, 10000
    Rational.metrics = {"max_coefficient_degree": 0, "rational_additions": 0, "rational_multiplications": 0}


def decode_scalar(data, prime):
    from flint import nmod_poly
    return Rational(nmod_poly(data["n"], prime), nmod_poly(data["d"], prime))


def load(data):
    prime = data["prime"]
    algebra = Algebra()
    algebra.relations = [(tuple(decode_scalar(c, prime) for c in row["a"]),
                          tuple(decode_scalar(c, prime) for c in row["b"])) for row in data["relations"]]
    if len(algebra.relations) != 5 or any(len(a) != 2**j or len(b) != 2**j for j, (a, b) in enumerate(algebra.relations)):
        raise ValueError("二次塔关系的维数错误")
    return algebra


def vector(data, algebra):
    if len(data) != 32:
        raise ValueError("规范形的基长度不是32")
    return Element(algebra, 5, tuple(decode_scalar(c, Rational.prime) for c in data))


def equal(a, b):
    a, b = a.promote(5), b.promote(5)
    return all(x.n == y.n and x.d == y.d for x, y in zip(a.coefficients, b.coefficients))


def irreducible(poly, prime):
    from flint import nmod_poly
    degree = int(poly.degree())
    if degree < 1:
        return False
    divisors, remaining, trial = [], degree, 2
    while trial*trial <= remaining:
        if remaining%trial == 0:
            divisors.append(trial)
            while remaining%trial == 0:
                remaining //= trial
        trial += 1
    if remaining > 1:
        divisors.append(remaining)
    X = nmod_poly([0, 1], prime)
    # Rabin判据独立核验不可约性，不依赖产生端因子分解的标签。
    return X.pow_mod(prime**degree, poly) == X%poly and all(
        (X.pow_mod(prime**(degree//q), poly)-X).gcd(poly).is_one() for q in divisors)


def verify_symbolic(data, norm_data, circuit):
    from flint import nmod_poly
    from build_shared import construct
    rebuilt = construct()
    for key in ("variables", "selected_sign", "nodes", "outputs", "all_divisor_nodes",
                "tower_leading_coefficient_nodes", "simple_root_derivative_nodes", "critical_locus_nonzero_guards"):
        if circuit[key] != rebuilt[key]:
            raise ValueError("符号核验的共享根电路与构造器不一致")
    if data["source_system_sha256"] != sha(AXIS) or len(data["fixed_u_v_w"]) != 3:
        raise ValueError("规范形的输入绑定或固定参数个数错误")
    if data["derivatives_computed"] and len(data["gradient_u_v_w_r"]) != 4:
        raise ValueError("规范形必须包含四个全导数")
    prime = data["prime"]
    setup(prime)
    algebra = load(data)
    generators = []
    for j in range(5):
        coefficients = [algebra.zero]*(2**(j+1))
        coefficients[2**j] = algebra.one
        generators.append(Element(algebra, j+1, tuple(coefficients)))
    parameters = [algebra.scalar(value) for value in data["fixed_u_v_w"]]+[algebra.scalar(Rational(nmod_poly([0, 1], prime)))]
    registers = execute(circuit, generators+parameters)
    if any(not registers[n].is_zero() for n in circuit["outputs"]["finite_relations"]):
        raise ValueError("符号二次塔没有满足五条电路关系")
    expected_H = vector(data["H"], algebra)
    if not equal(registers[circuit["outputs"]["critical_equations"][0]], expected_H):
        raise ValueError("符号闭环规范形与电路不一致")
    if data["derivatives_computed"]:
        for n, row in zip(circuit["outputs"]["critical_equations"][1:]+[circuit["outputs"]["radius_derivative"]], data["gradient_u_v_w_r"]):
            if not equal(registers[n], vector(row, algebra)):
                raise ValueError("符号全导数规范形与电路不一致")
    for name, nodes in circuit["outputs"]["geometry"].items():
        if len(data["geometry"][name]) != len(nodes):
            raise ValueError("坐标规范形的维数错误")
        for node, row in zip(nodes, data["geometry"][name]):
            if not equal(registers[node], vector(row, algebra)):
                raise ValueError("符号坐标规范形与电路不一致")
    actual_norm = algebra.norm(expected_H.coefficients, 5)
    expected_norm = decode_scalar(norm_data["value"], prime)
    if actual_norm.n != expected_norm.n or actual_norm.d != expected_norm.d:
        raise ValueError("递归符号范数不一致")
    for side in ("numerator", "denominator"):
        product = nmod_poly([norm_data[side+"_unit"]], prime)
        for row in norm_data[side+"_factors"]:
            factor = nmod_poly(row["coefficients"], prime)
            if (int(factor.degree()) != row["degree"] or row["multiplicity"] < 1
                    or not irreducible(factor, prime)):
                raise ValueError("返回的范数因子并非不可约")
            product *= factor**row["multiplicity"]
        if product != (actual_norm.n if side == "numerator" else actual_norm.d):
            raise ValueError("范数因子没有精确乘回")
    return {"verified": True, "symbolic_finite_relations": 5, "symbolic_total_derivatives": 4 if data["derivatives_computed"] else 0,
            "norm_numerator_degree": int(actual_norm.n.degree()), "norm_denominator_degree": int(actual_norm.d.degree()),
            "irreducible_factorization_multiplied_back": True,
            "scope": "精确Fp(r)恒等式，u,v,w固定；未证明四参数全局恒等式或半径消元多项式。"}


def evaluate_vector(encoded, ports, radius, field, prime):
    from flint import nmod_poly
    weights = [field(1)]
    for port in ports:
        weights += [value*port for value in weights]
    result = field(0)
    for coefficient, weight in zip(encoded, weights):
        n, d = nmod_poly(coefficient["n"], prime), nmod_poly(coefficient["d"], prime)
        denominator = int(d(radius))
        if denominator == 0:
            raise ZeroDivisionError("规范形系数在测试半径有极点")
        result += field(int(n(radius))*pow(denominator, -1, prime)%prime)*weight
    return result


def verify_samples(data, norm_data, circuit, radii):
    from flint import fq_default_ctx, nmod_poly
    prime, fixed = data["prime"], data["fixed_u_v_w"]
    field = fq_default_ctx(prime, 8, "a", fq_type="FQ_NMOD")
    factor, system = read(FACTOR), read(AXIS)
    norm = norm_data["value"]
    N, D = nmod_poly(norm["n"], prime), nmod_poly(norm["d"], prime)
    reports, skipped = [], []
    for radius in radii:
        try:
            fiber = compute(field, fixed+[radius], factor)
            accepted = check_fiber(circuit, fiber, field, system)
            if int(D(radius)) == 0 or accepted["norm_h3"] != int(N(radius))*pow(int(D(radius)), -1, prime)%prime:
                raise ValueError("符号范数与独立原几何纤维乘积不一致")
            for i, j in accepted["selected_pairs"]:
                left, right = fiber["left"][i], fiber["right"][j]
                u, v, w, r = map(field, fixed+[radius])
                t, dx, L = field(left["t"]["value"]), v-u, (v-u)**2+w*w
                aa, bb = dx-w*t/L, w+dx*t/L
                ports = [bb/(2*r+aa), field(left["q0"][0]["value"]), -field(left["point"][1]["value"]),
                         field(right["c6"][1]["value"]), field(right["point"][0]["value"])]
                bridge = fiber["bridge_values"][8*i+j]
                expected = [field(bridge["value"]), *map(field, bridge["gradient"])]
                vectors = [data["H"]]+(data["gradient_u_v_w_r"] if data["derivatives_computed"] else [])
                actual = [evaluate_vector(row, ports, radius, field, prime) for row in vectors]
                if actual != expected[:len(vectors)]:
                    raise ValueError("规范形的某个分支值或全导数与原几何不一致")
            reports.append({"radius": radius, "verified_points": 32, "norm_h3": accepted["norm_h3"]})
        except ZeroDivisionError as error:
            skipped.append({"radius": radius, "reason": str(error)})
    if len(reports) < 2:
        raise ValueError("有效独立几何测试点不足")
    return {"verified": True, "samples": reports, "skipped": skipped,
            "field": {"prime": prime, "degree": 8, "modulus": list(map(int, field.modulus().coeffs()))}}


def run(args):
    target = HERE / "quotient-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, norm = read(target / "quotient.json.gz"), read(target / "norm.json.gz")
        circuit = read(HERE / "system-runs/20261015-selected/circuit.json.gz")
        symbolic = verify_symbolic(data, norm, circuit)
        samples = verify_samples(data, norm, circuit, args.radii)
        report = {"verified": True, "symbolic": symbolic, "samples": samples, "seconds": time.perf_counter()-tick,
                  "quotient_sha256": sha(target / "quotient.json.gz"), "norm_sha256": sha(target / "norm.json.gz"),
                  "circuit_sha256": sha(HERE / "system-runs/20261015-selected/circuit.json.gz"),
                  "verifier_sha256": sha(__file__), "limits": snapshot()}
        write(target / args.output, report)
        print({"verified": True, "symbolic_derivatives": symbolic["symbolic_total_derivatives"],
               "sample_points": sum(row["verified_points"] for row in samples["samples"]), "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--radii", type=int, nargs="+", default=[2, 3, 4, 5, 7, 11, 17, 29])
    parser.add_argument("--output", default="verified-final.json.gz")
    run(parser.parse_args())
