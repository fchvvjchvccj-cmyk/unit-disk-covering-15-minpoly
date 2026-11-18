"""用原几何纤维和严格区间独立核验32维共享根分量及九维目标单根。"""

import argparse
from fractions import Fraction as F
import time

from bridge_support import HERE, PARENT, TARGET, AXIS, read, write, sha, bounded
from build_shared import construct
from chain_fiber import Jet, prime_scalar
from verify_system import execute
from target_chart import Dual, image_box
from prove_shared import verify as verify_selection
import verify_fiber as parent_fiber


def check_fiber(data, fiber, field, system):
    parent_fiber.verify(fiber, field, system)
    params = [Jet(field(v), [field(int(i == j)) for j in range(4)]) for i, v in enumerate(fiber["parameters_u_v_w_r"])]
    u, v, w, r = params
    pairs, product = [], Jet(field(1))
    for i, left in enumerate(fiber["left"]):
        t = Jet(field(left["t"]["value"]), [field(a) for a in left["t"]["gradient"]])
        dx, L = v-u, (v-u)**2+w*w
        aa, bb = dx-w*t/L, w+dx*t/L
        s = bb/(2*r+aa)
        wanted = (1+r+v+dx-aa)/2
        matches = 0
        for j, right in enumerate(fiber["right"]):
            x = right["c6"][0]
            if field(x["value"]) != wanted.value:
                continue
            if tuple(map(field, x["gradient"])) != wanted.gradient:
                raise ValueError("共享根的右端四个微分不一致")
            values = [s.value, field(left["q0"][0]["value"]), -field(left["point"][1]["value"]),
                      field(right["c6"][1]["value"]), field(right["point"][0]["value"]), *[p.value for p in params]]
            registers = execute(data, values)
            if any(not registers[n].is_zero() for n in data["outputs"]["finite_relations"]):
                raise ValueError("共享根分量的五条载体方程不为零")
            for n in data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+data["simple_root_derivative_nodes"]:
                if registers[n].is_zero():
                    raise ValueError("共享根纤维落到坏主元")
            expected = {"c0": left["point"], "c3": left["c3"], "c4": left["c4"], "c5": right["point"], "c6": right["c6"]}
            for name, nodes in data["outputs"]["geometry"].items():
                if [registers[n] for n in nodes] != [field(a["value"]) for a in expected[name]]:
                    raise ValueError("共享根图表重建坐标不一致")
            bridge = fiber["bridge_values"][8*i+j]
            expected_h = [field(bridge["value"]), *map(field, bridge["gradient"])]
            indices = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
            if [registers[n] for n in indices] != expected_h:
                raise ValueError("共享根临界式的闭环值或全导数错误")
            product *= Jet(expected_h[0], expected_h[1:])
            pairs.append([i, j])
            matches += 1
        if matches != 4:
            raise ValueError("每个左端根没有恰好四个同号右端根")
    if len(pairs) != 32:
        raise ValueError("选出的共享根分量秩不是32")
    return {"verified": True, "rank": 32, "selected_pairs": pairs,
            "norm_h3": prime_scalar(product.value), "norm_gradient_u_v_w_r": [prime_scalar(a) for a in product.gradient]}


def box_system(data, root, bits):
    arithmetic, old_ports = image_box(root, bits, ["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"])
    ports = old_ports[:3]+old_ports[4:]
    registers = execute(data, [Dual(v, [v*0+int(i == j) for j in range(9)]) for i, v in enumerate(ports)])
    equations = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    if len(equations) != 9:
        raise ValueError("共享根临界系统的方程数错误")
    for n in equations:
        if not registers[n].value.value[0] <= 0 <= registers[n].value.value[1]:
            raise ValueError("共享根方程与目标根不相容")
    guards = sorted(set(data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+
                        data["simple_root_derivative_nodes"]+data["critical_locus_nonzero_guards"]))
    for n in guards:
        lo, hi = registers[n].value.value
        if lo <= 0 <= hi:
            raise ValueError("目标共享根图表的主元包含零")
    return arithmetic, ports, registers, equations, guards


def verify_target(certificate, data, root):
    rebuilt = construct()
    for key in ("variables", "selected_sign", "nodes", "outputs", "all_divisor_nodes", "tower_leading_coefficient_nodes",
                "simple_root_derivative_nodes", "critical_locus_nonzero_guards"):
        if data[key] != rebuilt[key]:
            raise ValueError("共享根电路与构造器不一致")
    if certificate["root_sha256"] != sha(TARGET / "refined-root.json") or certificate["source_system_sha256"] != sha(AXIS):
        raise ValueError("共享根目标证书的输入绑定错误")
    if root != read(TARGET / "refined-root.json"):
        raise ValueError("内存中的根盒与已认证输入不一致")
    verify_selection(root)
    arithmetic, ports, registers, equations, guards = box_system(data, root, certificate["bits"])
    numerators, Q = certificate["inverse_numerators"], int(certificate["denominator"])
    if Q <= 0 or len(numerators) != 9 or any(len(row) != 9 for row in numerators):
        raise ValueError("目标逆矩阵的维数错误")
    P = [[arithmetic.rational(F(int(v), Q)) for v in row] for row in numerators]
    J = [[v.value for v in registers[n].gradient] for n in equations]
    row_sums = []
    for i in range(9):
        total = 0
        for j in range(9):
            value = arithmetic.zero
            for k in range(9):
                value = arithmetic.add(value, arithmetic.mul(P[i][k], J[k][j]))
            total += arithmetic.absolute(arithmetic.sub(arithmetic.one if i == j else arithmetic.zero, value))
        row_sums.append(total)
    defect = F(max(row_sums), arithmetic.scale)
    if defect >= 1:
        raise ValueError("共享根九维雅可比未能认证可逆")
    encode = lambda value: [str(F(a, arithmetic.scale)) for a in value]
    return {"verified": True, "jacobian_dimension": 9, "jacobian_invertible": True,
            "inverse_defect_upper": str(defect), "target_local_length": 1,
            "nonzero_guards": {str(n): encode(registers[n].value.value) for n in guards},
            "target_ports": [encode(v.value) for v in ports],
            "radius_derivative_interval": encode(registers[data["outputs"]["radius_derivative"]].value.value),
            "exact_membership_argument": "沿已接受的原临界图表选择共享根同号分量，另一符号因子为单位，故删除x方程和坐标保持目标局部环。其余五级单根关系和原KKT拉回给出精确九条零式；区间只认证相容性、非零主元和雅可比可逆。",
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    import mpmath as mp
    source = HERE / "system-runs" / args.system
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, root, system = read(source / "circuit.json.gz"), read(TARGET / "refined-root.json"), read(AXIS)
        arithmetic, _, registers, equations, _ = box_system(data, root, 1024)
        with mp.workdps(150):
            J = mp.matrix([[mp.mpf(g.value[0]+g.value[1])/(2*arithmetic.scale) for g in registers[n].gradient] for n in equations])
            inverse, Q = J**-1, 10**100
            certificate = {"bits": 1024, "denominator": str(Q),
                           "inverse_numerators": [[str(int(mp.nint(inverse[i, j]*Q))) for j in range(9)] for i in range(9)],
                           "root_sha256": sha(TARGET / "refined-root.json"), "source_system_sha256": sha(AXIS),
                           "circuit_sha256": sha(source / "circuit.json.gz"), "verifier_sha256": sha(__file__)}
        report = verify_target(certificate, data, root)
        fibers = []
        for prime in (65521, 65519):
            folder = PARENT / "fiber-runs" / f"20261015-p{prime}"
            metadata = read(folder / "field.json")
            field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(prime)(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
            for path in sorted(folder.glob("fiber-*.json.gz")):
                accepted = check_fiber(data, read(path), field, system)
                accepted.update(prime=prime, path=path.relative_to(HERE.parents[2]).as_posix(), sha256=sha(path))
                fibers.append(accepted)
        report.update(fibers=fibers, seconds=time.perf_counter()-tick, limits=snapshot())
        write(source / (args.prefix+"target-certificate.json.gz"), certificate)
        write(source / (args.prefix+"verified.json.gz"), report)
        print({"verified": True, "selected_points": sum(r["rank"] for r in fibers),
               "jacobian_defect_float": float(F(report["inverse_defect_upper"])), "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--prefix", default="final-")
    run(parser.parse_args())
