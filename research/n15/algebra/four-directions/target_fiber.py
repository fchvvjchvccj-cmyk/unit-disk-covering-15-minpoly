"""在实际半径的认证盒上，检查 104 次多项式的两个方向分支。"""

import argparse
from fractions import Fraction as Q
import mpmath as mp

from chart import HERE, direction, forward, geometric_critical, read
from context import write
from binary_projection import polynomial_value
from verify import independent_variables, require
from verify_chart import certify_system


def refine():
    old = read(HERE / "target.json")
    old_center = list(map(Q, old["center"]))
    with mp.workdps(430):
        mid = mp.matrix([mp.mpf(v.numerator)/v.denominator for v in old_center])
        for _ in range(4):
            rows = geometric_critical(independent_variables(list(mid)))
            mid -= mp.lu_solve(mp.matrix([row.gradient for row in rows]), mp.matrix([row.value for row in rows]))
        scale = 10**400
        center = [Q(int(mp.nint(v*scale)), scale) for v in mid]
    width = Q(1, 10**350)
    require(max(abs(a-b) for a, b in zip(old_center, center))+width < Q(old["half_width"]),
            "细化盒未包含在已绑定的四变量盒内")
    return center, width, certify_system(geometric_critical, center, width, 1536, 380)


def coefficients(e, a, nodes):
    from flint import acb
    values = [polynomial_value(e, a, w) for w in nodes]
    require(all(v.is_finite() for v in values), "单位根取值遇到图表极点")
    n = len(nodes)
    return [sum((value*nodes[-(j*k) % n] for k, value in enumerate(values)), acb(0))/n
            for j in range(105)]


def divide_at_root(coeffs, root):
    quotient = [coeffs[-1]]
    for c in reversed(coeffs[1:-1]):
        quotient.append(c+root*quotient[-1])
    require(0 in coeffs[0]+root*quotient[-1], "已知根的除法余数不相容")
    return quotient[::-1]


def analyze():
    from flint import arb, acb, acb_poly, ctx, fmpq
    center, width, certificate = refine()
    ctx.prec = 1450
    params = [arb(fmpq(m.numerator, m.denominator)) + arb(0, str(width)) for m in center]
    a, _, _, d = [acb(*direction(s)) for s in params]
    _, z, _, _, _, _, _, r = forward(params)
    e = acb(z, r*a.imag)
    nodes = [acb(fmpq(2*k, 128)).exp_pi_i() for k in range(128)]
    branches, isolated = [], []
    require(0 not in e+1/e, "椭圆曲线两支未分离")
    for index, branch in enumerate((e, -1/e)):
        coeffs = coefficients(branch, a, nodes)
        require(0 not in coeffs[-1], "实际纤维的 104 次首项未认证")
        removed = 0
        if index == 0:
            coeffs = divide_at_root(divide_at_root(coeffs, d), d)
            removed = 2
        poly = acb_poly(coeffs)
        require(0 not in poly(d), "目标根仍出现在剩余因子中")
        roots = poly.roots(tol=1e-18, maxprec=1450)
        require(len(roots) == 104-removed and all(
            not p.overlaps(q) for i, p in enumerate(roots) for q in roots[:i]), "复根隔离不完整")
        isolated.append(roots)
        branches.append({"degree_before_deflation": 104, "known_target_multiplicity_removed": removed,
                         "isolated_simple_roots": len(roots), "leading_abs_lower": str(abs(coeffs[-1]).lower()),
                         "deflated_at_D_abs_lower": str(abs(poly(d)).lower())})
    require(all(not p.overlaps(q) for p in isolated[0]+[d] for q in isolated[1]), "两个方向分支有重合根盒")
    return {"precision_bits": ctx.prec, "refined_center": list(map(str, center)), "width": str(width),
            "distinct_branches_have_disjoint_roots": True,
            "refined_root_check": certificate, "branches": branches, "target_in_multiplicity_one_layer": True,
            "radius_minimal_polynomial_obtained": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        report = analyze()
        write(args.output, report)
        print([item["isolated_simple_roots"] for item in report["branches"]])
    finally:
        close()
