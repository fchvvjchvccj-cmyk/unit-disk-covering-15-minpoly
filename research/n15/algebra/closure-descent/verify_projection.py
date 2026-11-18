"""用原认证根、恢复图表和严格区间验证四参数临界系统的目标单根。"""

import argparse
from fractions import Fraction as F
import hashlib
import json
import time

from descent_support import HERE, SHARED, TARGET, SOURCE, AXIS, read, write, sha, bounded
from build_projected import construct
from verify_shared import verify_target as parent_target
from target_chart import image_box, Dual
from verify_system import execute
from prove_relative import verify as relative_identities


def source_acceptance(root):
    folder = SOURCE.parent
    certificate = read(folder / "final-target-certificate.json.gz")
    report = parent_target(certificate, read(SOURCE), root)
    if not report["verified"] or not report["jacobian_invertible"]:
        raise ValueError("原九变量目标系统未通过接受端")
    return {"verified": True, "certificate_sha256": sha(folder / "final-target-certificate.json.gz"),
            "source_circuit_sha256": sha(SOURCE), "source_local_length": 1}


def content_sha(data):
    return hashlib.sha256(json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def box(data, root, bits, derivatives=True):
    arithmetic, old = image_box(root, bits, ["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"])
    k, n = data["keep_left_roots"], len(data["variables"])
    ports = old[:k]+old[6:]
    values = [Dual(v, [v*0+int(i == j) for j in range(n)]) for i, v in enumerate(ports)] if derivatives else ports
    registers = execute(data, values)
    interval = lambda index: registers[index].value.value if derivatives else registers[index].value
    equations = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    if len(equations) != n:
        raise ValueError("投影系统的方程个数错误")
    if any(not interval(index)[0] <= 0 <= interval(index)[1] for index in equations):
        raise ValueError("投影零式与认证目标根不相容")
    guards = sorted(set(sum([data[name] for name in ["all_divisor_nodes", "tower_leading_coefficient_nodes",
                        "simple_root_derivative_nodes", "closed_chart_nonzero_guards", "critical_locus_nonzero_guards"]], [])))
    if any(interval(index)[0] <= 0 <= interval(index)[1] for index in guards):
        raise ValueError("投影图表的单位或导数乘子没有排除零")
    expected = {name: old[i] for i, name in enumerate(["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"])}
    for name, index in data["outputs"]["recovered_left_and_y"].items():
        difference = arithmetic.sub(interval(index), expected[name].value)
        if not difference[0] <= 0 <= difference[1]:
            raise ValueError("被消去生成元的恢复值与目标分支不一致")
    difference = arithmetic.sub(interval(data["outputs"]["recovered_d"]), expected["d"].value)
    if not difference[0] <= 0 <= difference[1]:
        raise ValueError("最后二次根的恢复值与目标分支不一致")
    return arithmetic, ports, registers, equations, guards


def verify(certificate, data, root):
    rebuilt = construct(data["keep_left_roots"])
    for key in ("variables", "nodes", "outputs", "all_divisor_nodes", "tower_leading_coefficient_nodes",
                "simple_root_derivative_nodes", "closed_chart_nonzero_guards", "critical_locus_nonzero_guards"):
        if data[key] != rebuilt[key]:
            raise ValueError("投影电路与构造器不一致")
    if (certificate["root_sha256"] != sha(TARGET / "refined-root.json") or root != read(TARGET / "refined-root.json")
            or certificate["source_circuit_sha256"] != sha(SOURCE) or certificate["source_system_sha256"] != sha(AXIS)):
        raise ValueError("投影证书没有绑定原根和方程")
    if certificate["circuit_content_sha256"] != content_sha(data):
        raise ValueError("投影证书与实际电路内容失去绑定")
    parents = source_acceptance(root)
    universal = relative_identities()
    arithmetic, ports, registers, equations, guards = box(data, root, certificate["bits"])
    n, Q = len(equations), int(certificate["denominator"])
    numerators = certificate["inverse_numerators"]
    if Q <= 0 or len(numerators) != n or any(len(row) != n for row in numerators):
        raise ValueError("投影雅可比的有理预条件器维数错误")
    P = [[arithmetic.rational(F(int(value), Q)) for value in row] for row in numerators]
    J = [[entry.value for entry in registers[index].gradient] for index in equations]
    sums = []
    for i in range(n):
        total = 0
        for j in range(n):
            product = arithmetic.zero
            for k in range(n):
                product = arithmetic.add(product, arithmetic.mul(P[i][k], J[k][j]))
            total += arithmetic.absolute(arithmetic.sub(arithmetic.one if i == j else arithmetic.zero, product))
        sums.append(total)
    defect = F(max(sums), arithmetic.scale)
    if defect >= 1:
        raise ValueError("投影雅可比未能以严格区间认证可逆")
    encode = lambda value: [str(F(v, arithmetic.scale)) for v in value]
    report = {"verified": True, "parent_acceptance": parents, "universal_identities": universal,
              "variables": data["variables"], "jacobian_dimension": n, "jacobian_invertible": True,
              "inverse_defect_upper": str(defect), "target_local_length": 1,
              "recovered": {name: encode(registers[index].value.value) for name, index in data["outputs"]["recovered_left_and_y"].items()},
              "recovered_d": encode(registers[data["outputs"]["recovered_d"]].value.value),
              "critical_multiplier_interval": encode(registers[data["outputs"]["critical_multiplier"]].value.value),
              "radius_derivative_interval": encode(registers[data["outputs"]["radius_derivative"]].value.value),
              "nonzero_guards": {str(index): encode(registers[index].value.value) for index in guards},
              "exact_membership_argument": "原九变量认证根满足精确闭环临界式。通用整数恒等式把H约成J0+J1*d；逐级范数在闭环零式上只乘共轭单位。各恢复一次项、根简单导数和全部图表分母均严格非零，故删除生成元保持目标局部临界环；区间只用于非零、相容性和雅可比可逆认证。",
              "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}
    return report


def run(args):
    import mpmath as mp
    source = HERE / "system-runs" / args.system
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, root = read(source / "circuit.json.gz"), read(TARGET / "refined-root.json")
        prepared = box(data, root, args.bits)
        arithmetic, _, registers, equations, _ = prepared
        n = len(equations)
        with mp.workdps(220):
            J = mp.matrix([[mp.mpf(v.value[0]+v.value[1])/(2*arithmetic.scale) for v in registers[index].gradient] for index in equations])
            inverse, Q = J**-1, 10**140
            certificate = {"bits": args.bits, "denominator": str(Q),
                           "inverse_numerators": [[str(int(mp.nint(inverse[i, j]*Q))) for j in range(n)] for i in range(n)],
                           "root_sha256": sha(TARGET / "refined-root.json"), "source_circuit_sha256": sha(SOURCE),
                           "source_system_sha256": sha(AXIS), "circuit_sha256": sha(source / "circuit.json.gz"),
                           "circuit_content_sha256": content_sha(data),
                           "verifier_sha256": sha(__file__)}
        report = verify(certificate, data, root)
        report.update(seconds=time.perf_counter()-tick, limits=snapshot())
        write(source / args.certificate, certificate)
        write(source / args.output, report)
        print({"verified": True, "dimension": n, "inverse_defect_float": float(F(report["inverse_defect_upper"])),
               "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--bits", type=int, default=2048)
    parser.add_argument("--certificate", default="target-certificate.json.gz")
    parser.add_argument("--output", default="target-verified.json.gz")
    run(parser.parse_args())
