"""在原认证根的像盒上，用纯整数区间认证临界系统雅可比可逆。"""

import argparse
from fractions import Fraction as F
import hashlib
import json
import time

from support import HERE, TARGET, AXIS, FACTOR, read, write, sha, bounded
from prove_branch import Interval
from geometry import circumcenter
from verify_system import execute


class Dual:
    def __init__(self, value, gradient):
        self.value, self.gradient = value, list(gradient)

    def cast(self, other):
        return other if isinstance(other, Dual) else Dual(self.value*0+other, [self.value*0 for _ in self.gradient])

    def __add__(self, other):
        other = self.cast(other)
        return Dual(self.value+other.value, [a+b for a, b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-a for a in self.gradient])

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        return Dual(self.value*other.value, [a*other.value+self.value*b for a, b in zip(self.gradient, other.gradient)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        value = self.value/other.value
        return Dual(value, [(a-value*b)/other.value for a, b in zip(self.gradient, other.gradient)])

    def __rtruediv__(self, other):
        return self.cast(other)/self


DEPENDENCIES = {
    "original_root_and_stationarity": HERE.parent / "chain-radius/critical-runs/20261015-target/verified.json.gz",
    "rational_contact_branch": HERE / "proof-runs/20261015-complementary/verified.json.gz",
    "left_factor_selection": FACTOR,
    "right_B2_quotient": HERE.parent / "three-parameter-elimination/localized-quotient/quotient-final-verified.json",
    "right_B2_selector": HERE.parent / "three-parameter-elimination/localized-quotient/target-selector-final.json",
}


def dependencies(data=None):
    reports = {name: read(path) for name, path in DEPENDENCIES.items()}
    if not all(reports[k]["verified"] for k in reports if k != "right_B2_selector"):
        raise ValueError("父证明未通过")
    parent = reports["original_root_and_stationarity"]["root_replay"]
    if parent["refined_root_sha256"] != sha(TARGET / "refined-root.json") or parent["system_sha256"] != sha(AXIS):
        raise ValueError("原认证根的绑定改变")
    branch = reports["rational_contact_branch"]
    if branch["root_sha256"] != sha(TARGET / "refined-root.json") or branch["source_system_sha256"] != sha(AXIS):
        raise ValueError("互补交点分支证明的绑定改变")
    if not reports["left_factor_selection"]["factorization_checked_over_characteristic_zero"]:
        raise ValueError("左端因子没有特征零乘回证明")
    if reports["right_B2_selector"]["selected"] != "B2" or reports["right_B2_quotient"]["target_carrier"] != "B2":
        raise ValueError("右端分量选择不一致")
    result = {name: {"path": path.relative_to(HERE.parents[2]).as_posix(), "sha256": sha(path)}
              for name, path in DEPENDENCIES.items()}
    if data and "z" in data["variables"]:
        from quadratic_system import prove_left_identity
        prove_left_identity()
        path = HERE / "quadratic_system.py"
        result["common_height_identity_recomputed"] = {"path": path.relative_to(HERE.parents[2]).as_posix(), "sha256": sha(path)}
    return result


def bind_circuit(data):
    if data["variables"] == ["s", "k", "x", "y", "d", "u", "v", "w", "r"]:
        from build_system import construct
    elif data["variables"] == ["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"]:
        from quadratic_system import construct
    else:
        raise ValueError("不认识目标临界系统的端口")
    rebuilt = construct()
    for key in ("variables", "nodes", "outputs", "all_divisor_nodes", "tower_leading_coefficient_nodes",
                "simple_root_derivative_nodes", "critical_locus_nonzero_guards", "source_system_sha256", "left_factor_sha256"):
        if data[key] != rebuilt[key]:
            raise ValueError("保存电路与精确构造器不一致："+key)
    return hashlib.sha256(json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def image_box(root, bits, names=None):
    from intervals import Dyadic
    arithmetic = Dyadic(bits)
    Interval.arithmetic = arithmetic
    rho = F(root["rho"])
    z = [Interval((arithmetic.rational(F(int(v), int(root["Qx"]))-rho)[0],
                   arithmetic.rational(F(int(v), int(root["Qx"]))+rho)[1])) for v in root["xnum"][:14]]
    a, u, b, r = [z[0], -z[1]], z[10], z[11:13], z[13]
    p = circumcenter(a, [u, 0*u], b)
    s, k = p[1]/(r+p[0]-u), z[1]/(z[0]-u)
    if names and "z" in names:
        from geometry import previous_boundary
        q = previous_boundary(a, u, r)
        return arithmetic, [s, q[0], a[1], z[8], z[9], z[6], u, b[0], b[1], r]
    return arithmetic, [s, k, z[8], z[9], z[6], u, b[0], b[1], r]


def interval_system(data, root, bits):
    arithmetic, box = image_box(root, bits, data["variables"])
    dimension = len(box)
    values = [Dual(v, [v*0+int(i == j) for j in range(dimension)]) for i, v in enumerate(box)]
    registers = execute(data, values)
    equations = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    if len(equations) != dimension:
        raise ValueError("临界系统的方程数或维数改变")
    residual = [registers[n].value.value for n in equations]
    if any(not lo <= 0 <= hi for lo, hi in residual):
        raise ValueError("九变量方程与原认证根不相容")
    guards = sorted(set(data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+
                        data["simple_root_derivative_nodes"]+data["critical_locus_nonzero_guards"]))
    for node in guards:
        lo, hi = registers[node].value.value
        if lo <= 0 <= hi:
            raise ValueError("目标根盒不能排除零主元："+str(node))
    return arithmetic, box, registers, equations, residual, guards


def verify(certificate, data, root):
    if certificate["root_sha256"] != sha(TARGET / "refined-root.json") or certificate["source_system_sha256"] != sha(AXIS):
        raise ValueError("雅可比证书输入绑定错误")
    if certificate["dependencies"] != dependencies(data):
        raise ValueError("雅可比证书的父证明绑定改变")
    if certificate["circuit_payload_sha256"] != bind_circuit(data):
        raise ValueError("雅可比证书的电路绑定错误")
    arithmetic, box, registers, equations, residual, guards = interval_system(data, root, certificate["interval_bits"])
    dimension = len(box)
    Q = int(certificate["denominator"])
    entries = certificate["inverse_numerators"]
    if len(entries) != dimension or any(len(row) != dimension for row in entries):
        raise ValueError("预条件器的维数错误")
    P = [[arithmetic.rational(F(int(v), Q)) for v in row] for row in entries]
    J = [[v.value for v in registers[n].gradient] for n in equations]
    errors = []
    for i in range(dimension):
        row = []
        for j in range(dimension):
            product = arithmetic.zero
            for k in range(dimension):
                product = arithmetic.add(product, arithmetic.mul(P[i][k], J[k][j]))
            row.append(arithmetic.sub(arithmetic.one if i == j else arithmetic.zero, product))
        errors.append(row)
    defect = F(max(sum(arithmetic.absolute(v) for v in row) for row in errors), arithmetic.scale)
    if defect >= 1:
        raise ValueError("不能认证临界系统雅可比可逆")
    encode = lambda interval: [str(F(v, arithmetic.scale)) for v in interval]
    return {"verified": True, "jacobian_dimension": dimension, "jacobian_invertible": True,
            "inverse_defect_infinity_norm_upper": str(defect),
            "target_ports": dict(zip(data["variables"], [encode(v.value) for v in box])),
            "nonzero_guards": {str(n): encode(registers[n].value.value) for n in guards},
            "radius_derivative_interval": encode(registers[data["outputs"]["radius_derivative"]].value.value),
            "residual_intervals_consistency_only": list(map(encode, residual)),
            "target_local_length": 1,
            "exact_root_membership_dependencies": certificate["dependencies"],
            "argument": "原精确KKT根经已证明的互补交点及四次/B2分支映射进入临界系统。前面的三角方程定义单根有限载体；原拉格朗日恒等式在该图表给出H=Hu=Hv=Hw=0。区间仅检查映射相容性和非零主元；精确零依赖上述恒等式和父证明。纯整数向外区间证明||I-PJ||_infinity<1，故目标解为本系统的局部孤立单根。",
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    import mpmath as mp
    source = HERE / "system-runs" / args.system
    target = source / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有目标图表")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, root = read(source / "circuit.json.gz"), read(TARGET / "refined-root.json")
        parents = dependencies(data)
        arithmetic, _, registers, equations, _, _ = interval_system(data, root, 1024)
        dimension = len(data["variables"])
        with mp.workdps(150):
            middle = lambda a: mp.mpf(a.value[0]+a.value[1])/(2*arithmetic.scale)
            matrix = mp.matrix([[middle(v) for v in registers[n].gradient] for n in equations])
            inverse = matrix**-1
            Q = 10**100
            certificate = {"interval_bits": 1024, "denominator": str(Q),
                           "inverse_numerators": [[str(int(mp.nint(inverse[i, j]*Q))) for j in range(dimension)] for i in range(dimension)],
                           "root_sha256": sha(TARGET / "refined-root.json"), "source_system_sha256": sha(AXIS),
                           "circuit_sha256": sha(source / "circuit.json.gz"), "dependencies": parents,
                           "circuit_payload_sha256": bind_circuit(data),
                           "verifier_sha256": sha(__file__)}
        report = verify(certificate, data, root)
        report.update(seconds=time.perf_counter()-tick, limits=snapshot())
        write(target / "certificate.json.gz", certificate)
        write(target / "verified.json.gz", report)
        print({"verified": True, "local_length": 1, "nonzero_guards": len(report["nonzero_guards"]),
               "defect_upper_float": float(F(report["inverse_defect_infinity_norm_upper"])),
               "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
