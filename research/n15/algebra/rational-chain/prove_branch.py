"""核验整数几何恒等式，并在认证根盒内严格选择两个互补交点分支。"""

import argparse
from fractions import Fraction as F
import time

from support import HERE, TARGET, AXIS, read, write, sha, bounded
from geometry import dot, subtract, common_boundary, circumcenter, previous_boundary, reconstruct


class Interval:
    arithmetic = None

    def __init__(self, value):
        self.value = value

    def cast(self, other):
        return other if isinstance(other, Interval) else Interval(self.arithmetic.rational(other))

    def __add__(self, other):
        return Interval(self.arithmetic.add(self.value, self.cast(other).value))

    __radd__ = __add__

    def __neg__(self):
        return Interval(self.arithmetic.neg(self.value))

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        return Interval(self.arithmetic.mul(self.value, self.cast(other).value))

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Interval(self.arithmetic.div(self.value, self.cast(other).value))

    def __rtruediv__(self, other):
        return self.cast(other)/self


def source_identities(system):
    import sympy as S
    z = S.symbols("x:14")
    a, c3, c4, b = [z[0], -z[1]], list(z[2:4]), list(z[4:6]), list(z[11:13])
    u, r = z[10], z[13]

    def triangle(a, b, c):
        d, e = subtract(b, a), subtract(c, a)
        area = d[0]*e[1]-d[1]*e[0]
        return dot(d, d)*dot(e, e)*dot(subtract(b, c), subtract(b, c))-4*r*r*area*area

    from geometry import boundary
    t, T, beta = 2*z[0]-u, z[0]**2+z[1]**2, 1-r*r
    left = z[1]**2*(t*t+beta)**2+(t*(T+beta)-z[0]*(t*t+beta))**2-4*t*t*z[1]**2
    expected = {1: boundary(a, c3, r), 2: boundary(c3, c4, r), 5: triangle(a, c3, b),
                6: triangle(a, [u, 0], b), 7: triangle(c3, c4, b), 10: left}
    for index, expression in expected.items():
        original = sum(term["coefficient"]*S.prod(z[j]**e for j, e in term["powers"])
                       for term in system["equations"][index])
        if not S.Poly(expression-original, *z).is_zero:
            raise ValueError("原整数几何恒等式失败：h"+str(index))
    return sorted(expected)


def verify(root, system):
    from intervals import Dyadic
    identities = source_identities(system)
    arithmetic = Dyadic(1024)
    Interval.arithmetic = arithmetic
    rho = F(root["rho"])
    values = [Interval((arithmetic.rational(F(int(v), int(root["Qx"]))-rho)[0],
                        arithmetic.rational(F(int(v), int(root["Qx"]))+rho)[1])) for v in root["xnum"][:14]]
    a, c3, c4, b = [values[0], -values[1]], values[2:4], values[4:6], values[11:13]
    u, r = values[10], values[13]
    q0, p0 = previous_boundary(a, u, r), circumcenter(a, [u, 0*u], b)
    q1, p1 = common_boundary(a, c3, r), circumcenter(a, c3, b)
    q2, p2 = common_boundary(c3, c4, r), circumcenter(c3, c4, b)
    guards = {}

    def positive(label, value):
        if value.value[0] <= 0:
            raise ValueError("不能认证非退化条件："+label)
        guards[label] = [str(F(v, arithmetic.scale)) for v in value.value]

    for label, pair in (("q1_ne_q0", subtract(q1, q0)), ("p1_ne_p0", subtract(p1, p0)),
                        ("q2_ne_q1", subtract(q2, q1)), ("p2_ne_p1", subtract(p2, p1)),
                        ("a_ne_c3", subtract(a, c3)), ("c3_ne_c4", subtract(c3, c4)),
                        ("q1_ne_p1", subtract(q1, p1)), ("q2_ne_p2", subtract(q2, p2)),
                        ("a_nonzero", a), ("c3_nonzero", c3),
                        ("a_ne_b", subtract(a, b)), ("c3_ne_b", subtract(c3, b))):
        positive(label+"_squared_distance", dot(pair, pair))
    predicted = reconstruct(values[:2], u, b, r)
    consistency = {}
    for name, actual in (("c3", c3), ("c4", c4)):
        errors = [arithmetic.sub(v.value, w.value) for v, w in zip(predicted[name], actual)]
        if any(not lo <= 0 <= hi for lo, hi in errors):
            raise ValueError("目标坐标不相容于有理重建")
        consistency[name] = [[str(F(v, arithmetic.scale)) for v in interval] for interval in errors]
    return {"verified": True, "integer_source_identities_checked": identities,
            "strict_branch_guards": guards, "interval_consistency_only": consistency,
            "exact_branch_selection_argument": "h10给出旧边界交点q0，h6给出旧三角形外心p0。原h1,h5给出q1,p1。严格区分两个交点后，韦达和为q0+q1=Ka*a/(a.a)，p0+p1=a+b；等半径圆的两个交点和为a+c3=q1+p1。同理h2,h7及严格分支区分给出c4=b+(1-r^2)c3/(c3.c3)-a。精确相等由这些恒等式及原方程推出，不由区间包含零推出。",
            "eliminated_coordinate_indices": [2, 3, 4, 5], "remaining_coordinate_count": 10,
            "remaining_geometric_constraint_indices": [0, 3, 4, 6, 8, 9, 10],
            "remaining_geometric_constraint_count": 7, "local_critical_square_system_coordinates": 10,
            "generic_selected_carrier_rank_before_closure": 64,
            "scope": "目标根邻域内的有理分支压缩，不能直接将全局原理想替换为此分支"}


def run(args):
    source = HERE / "proof-runs" / args.tag
    if source.exists():
        raise FileExistsError("拒绝覆盖已有分支证明")
    source.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        root, system = read(TARGET / "refined-root.json"), read(AXIS)
        prior_path = HERE.parent / "chain-radius/critical-runs/20261015-target/verified.json.gz"
        prior = read(prior_path)["root_replay"]
        if not prior["verified"] or prior["refined_root_sha256"] != sha(TARGET / "refined-root.json") or prior["system_sha256"] != sha(AXIS):
            raise ValueError("先前认证根证书的绑定改变")
        report = verify(root, system)
        report.update(source_system_sha256=sha(AXIS), root_sha256=sha(TARGET / "refined-root.json"),
                      root_proof_reference=prior_path.relative_to(HERE.parents[2]).as_posix(),
                      root_proof_sha256=sha(prior_path), verifier_sha256=sha(__file__),
                      seconds=time.perf_counter()-tick, limits=snapshot())
        write(source / "verified.json.gz", report)
        print({k: report[k] for k in ("verified", "remaining_coordinate_count", "generic_selected_carrier_rank_before_closure", "seconds")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
