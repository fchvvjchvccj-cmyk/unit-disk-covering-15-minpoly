"""重验候选根及两种四变量临界方程，按有理外包检查全部图表选择。"""

import argparse
from fractions import Fraction as Q

import mpmath as mp

from chart import (RESULTS, circle_closure, critical, forward, geometric_critical,
                   inverse, old_critical, read)
from context import execute, prune, sha, write
from verify import (away, box_values, contraction, encode, independent_variables,
                    require, transfer_audit, universal_identities, verify_new, verify_sources)


def certify_system(function, center, width, bits=512, inverse_digits=75):
    box = box_values(center, width, bits)
    exact = box_values(center, Q(0), bits)
    jacobian = [v.gradient for v in function(independent_variables(box))]
    # 浮点数只用来提出有理预条件器；随后所有接受步骤均向外舍入。
    with mp.workdps(inverse_digits+15):
        mid = [mp.mpf(x.numerator)/x.denominator for x in center]
        matrix = mp.matrix([v.gradient for v in function(independent_variables(mid))])
        estimate = matrix**-1
        scale = 10**inverse_digits
        preconditioner = [[Q(int(mp.nint(estimate[i, j]*scale)), scale) for j in range(estimate.cols)]
                          for i in range(estimate.rows)]
    return contraction(function(exact), jacobian, preconditioner, width, bits)


def prove():
    _, source, center, _ = verify_sources()
    universal_identities()
    transfer_audit(center, Q(source["rho"]))
    circuit_file = RESULTS / "primitive16.json.gz"
    source_file = RESULTS / "target-final-certificate.json.gz"
    circuit, certificate = read(circuit_file), read(source_file)
    verify_new(circuit, center, Q(source["rho"]), certificate)
    carrier_center = list(map(Q, certificate["center"]))
    carrier_box = box_values(carrier_center, Q(certificate["radius"]))
    param_box = inverse(carrier_box)
    scale = 10**60
    param_center = [Q(round(x*scale), scale) for x in inverse(carrier_center)]
    width = Q(1, 10**42)
    for image, mid in zip(param_box, param_center):
        lo, hi = map(Q, encode(image))
        require(mid-width < lo <= hi < mid+width, "原根像未进入四变量盒")

    params = box_values(param_center, width)
    lifted = forward(params)
    r = lifted[-1]
    radius_pivot = forward(independent_variables(params))[-1].gradient[1]
    closure, geometry = circle_closure(independent_variables(params))
    left, right, inner, area, previous = geometry
    reduced, outputs = prune(circuit, circuit["outputs"]["geometry"]["c5"])
    evaluated = execute(reduced, carrier_box)
    c5 = [evaluated[j] for j in outputs]

    def squared(point):
        return sum((coordinate-center)**2 for coordinate, center in zip(point, c5))

    gap_previous = squared([x.value for x in previous])-r*r
    gap_fixed = squared([r*0+1, r*0])-r*r
    require(gap_previous.value[0] > 0 and gap_fixed.value[0] > 0, "圆心未分离另一个边界接点")
    for point in (left, right, inner):
        residual = squared([x.value for x in point])-r*r
        require(residual.value[0] <= 0 <= residual.value[1], "三接点与目标盒不相容")
    old = old_critical()
    return {
        "parameters": ["s_A", "s_B", "s_C", "s_D"],
        "center": list(map(str, param_center)), "half_width": str(width),
        "carrier_certificate_sha256": sha(source_file), "carrier_circuit_sha256": sha(circuit_file),
        "source_image_inside": True, "radius_partial_B": away(radius_pivot, "半径偏导"),
        "closure_partial_B": away(closure.gradient[1], "闭环偏导"),
        "twice_signed_area": away(area.value, "三点共线"),
        "other_anchor_squared_distance_excess": [encode(gap_previous), encode(gap_fixed)],
        "pulled_back_system": certify_system(lambda p: critical(p, old), param_center, width),
        "geometric_system": certify_system(geometric_critical, param_center, width),
        "univariate_radius_polynomial_obtained": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        report = prove()
        write(args.output, report)
        print("两套四变量临界系统的根盒均通过严格压缩检验。")
    finally:
        close()
