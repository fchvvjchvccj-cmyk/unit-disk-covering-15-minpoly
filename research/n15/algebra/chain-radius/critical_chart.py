"""在认证目标根盒上消去驻定乘子，核验十四坐标的局部临界点图表。"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import time

from common import HERE, AXIS, TARGET, bounded, read, write, sha


RETAIN = [i for i in range(11) if i != 3]


def terms(system):
    return [[(t["coefficient"], tuple(map(tuple, t["powers"]))) for t in row]
            for row in system["equations"]]


def check_stationarity(system):
    from intervals import derivatives
    rows = terms(system)
    jacobian = derivatives(rows[:11], 14)
    for j in range(14):
        expected = defaultdict(int)
        if j == 13:
            expected[()] = 1
        for i in range(11):
            for c, e in jacobian[i][j]:
                expected[tuple(sorted(e+((14+i, 1),)))] += c
        actual = defaultdict(int)
        for c, e in rows[11+j]:
            actual[tuple(sorted(e))] += c
        if {e: c for e, c in expected.items() if c} != {e: c for e, c in actual.items() if c}:
            raise ValueError("原系统的驻定行不等于dr+Σtheta_i dh_i的系数")
    return jacobian


def interval_jacobian(system, root, bits):
    from intervals import Dyadic, evaluate_rows
    arithmetic = Dyadic(bits)
    rho = F(root["rho"])
    box = [(arithmetic.rational(F(int(v), int(root["Qx"]))-rho)[0],
            arithmetic.rational(F(int(v), int(root["Qx"]))+rho)[1]) for v in root["xnum"]]
    jacobian = check_stationarity(system)
    flat = evaluate_rows(arithmetic, [c for row in jacobian for c in row], box)
    return arithmetic, box, [flat[i*14:(i+1)*14] for i in range(11)]


def verify(certificate, root, system):
    if certificate["root_sha256"] != sha(TARGET / "refined-root.json") or certificate["system_sha256"] != sha(AXIS):
        raise ValueError("临界图表的认证根或整数系统摘要不一致")
    arithmetic, box, J = interval_jacobian(system, root, certificate["interval_bits"])
    if certificate["retained_constraint_indices"] != RETAIN:
        raise ValueError("被消去的闭环约束索引改变")
    Q = int(certificate["denominator"])
    P = [[arithmetic.rational(F(int(v), Q)) for v in row] for row in certificate["inverse_numerators"]]
    V0 = [[arithmetic.rational(F(int(v), Q)) for v in row] for row in certificate["tangent_numerators"]]
    if len(P) != 10 or any(len(row) != 10 for row in P) or len(V0) != 10 or any(len(row) != 4 for row in V0):
        raise ValueError("预条件器或切向量的维数错误")
    A = [[J[i][j] for j in range(10)] for i in RETAIN]
    B = [[J[i][j] for j in range(10, 14)] for i in RETAIN]

    def dot(a, b):
        value = arithmetic.zero
        for x, y in zip(a, b):
            value = arithmetic.add(value, arithmetic.mul(x, y))
        return value

    error = [[arithmetic.sub(arithmetic.one if i == j else arithmetic.zero,
                             dot(P[i], [A[k][j] for k in range(10)])) for j in range(10)] for i in range(10)]
    q = F(max(sum(arithmetic.absolute(v) for v in row) for row in error), arithmetic.scale)
    if q >= 1:
        raise ValueError("无法认证十坐标隐函数雅可比可逆")
    # A V0+B的预条件残差给出真实切向量的统一无穷范数误差。
    residual = [[arithmetic.add(dot(A[i], [V0[k][j] for k in range(10)]), B[i][j])
                 for j in range(4)] for i in range(10)]
    radii = [F(max(arithmetic.absolute(dot(P[i], [residual[k][j] for k in range(10)]))
                   for i in range(10)), arithmetic.scale)/(1-q) for j in range(4)]
    V = [[arithmetic.add(V0[i][j], (arithmetic.rational(-radii[j])[0], arithmetic.rational(radii[j])[1]))
          for j in range(4)] for i in range(10)]
    reduced = [arithmetic.add(J[3][10+j], dot(J[3][:10], [V[k][j] for k in range(10)])) for j in range(4)]
    theta = box[17]
    if theta[0] <= 0 <= theta[1] or reduced[3][0] <= 0 <= reduced[3][1]:
        raise ValueError("闭环乘子或半径方向导数不能认证非零")
    if any(not a <= 0 <= b for a, b in reduced[:3]):
        raise ValueError("三个局部驻定方向与认证根不相容")
    identity = arithmetic.add(arithmetic.one, arithmetic.mul(theta, reduced[3]))
    if not identity[0] <= 0 <= identity[1]:
        raise ValueError("半径方向驻定恒等式与认证根不相容")
    encode = lambda interval: [str(F(v, arithmetic.scale)) for v in interval]
    return {"verified": True, "stationarity_polynomial_identities_checked": 14,
            "inverse_defect_infinity_norm_upper": str(q), "ten_coordinate_jacobian_invertible": True,
            "theta3_interval": encode(theta), "reduced_h3_gradient_u_v_w_r_intervals": list(map(encode, reduced)),
            "h3_radius_derivative_nonzero": True,
            "exact_zero_directions": ["u", "v", "w"],
            "argument": "已认证根满足精确KKT多项式。十个G约束的雅可比可逆，故有局部隐函数y(q)。将dr+Σtheta_i dh_i=0拉回G=0得dr+theta3 dH=0；theta3非零，故Hu=Hv=Hw=0且Hr=-1/theta3。",
            "equivalent_local_square_system": "G0=...=G9=0, h3=0, det([[Gy,Gqj],[h3y,h3qj]])=0 (qj=u,v,w)",
            "scope": "只认证目标根邻域的临界点图表；不声称全局分量、最优性或半径最小多项式"}


def run(args):
    import mpmath as mp
    from verify_root import verify as verify_root
    source = HERE / "critical-runs" / args.tag
    if source.exists():
        raise FileExistsError("拒绝覆盖临界图表实验")
    source.mkdir(parents=True)
    root, system = read(TARGET / "refined-root.json"), read(AXIS)
    with bounded() as snapshot:
        tick = time.perf_counter()
        root_report = verify_root(root)
        arithmetic, box, J = interval_jacobian(system, root, 1024)
        with mp.workdps(150):
            middle = lambda a: mp.mpf(a[0]+a[1])/(2*arithmetic.scale)
            A = mp.matrix([[middle(J[i][j]) for j in range(10)] for i in RETAIN])
            B = mp.matrix([[middle(J[i][j]) for j in range(10, 14)] for i in RETAIN])
            P = A**-1
            V = -P*B
            Q = 10**100
            certificate = {"retained_constraint_indices": RETAIN, "interval_bits": 1024, "denominator": str(Q),
                           "inverse_numerators": [[str(int(mp.nint(P[i, j]*Q))) for j in range(10)] for i in range(10)],
                           "tangent_numerators": [[str(int(mp.nint(V[i, j]*Q))) for j in range(4)] for i in range(10)],
                           "root_sha256": sha(TARGET / "refined-root.json"), "system_sha256": sha(AXIS),
                           "driver_sha256": sha(__file__)}
        report = verify(certificate, root, system)
        report.update(root_replay=root_report, seconds=time.perf_counter()-tick, limits=snapshot())
        write(source / "certificate.json.gz", certificate)
        write(source / "verified.json.gz", report)
        print({k: report[k] for k in ("verified", "seconds", "ten_coordinate_jacobian_invertible", "h3_radius_derivative_nonzero", "exact_zero_directions")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
