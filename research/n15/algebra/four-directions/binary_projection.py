"""从三方向关系构造关于 E、A 的判别式，并验证清分母及次数界。"""

import argparse
from fractions import Fraction as Q

from chart import HERE, read
from context import write
from eliminate import (anchor_quadratic, discriminant27, endpoint, pencil_discriminant,
                       primitive_resultant)
from verify import box_values, require
from verify_elimination import base_coords, norm_guard


def polynomial_value(e, a, d):
    r = (e-1/e)/(a-1/a)
    return d**36*(d+r)**16*primitive_resultant(e, a, d)


def identity_and_orders():
    from flint import fmpz_mpoly_ctx
    ring = fmpz_mpoly_ctx.get(("x", "y", "z", "w", "t"))
    x, y, z, w, t = ring.gens()
    coeffs = [y*(y+z), y*(w+z*t)+w*(y+z),
              y*(x+z)+w*(w+z*t)+x*(y+z)+z*z,
              x*(w+z*t)+w*(x+z), x*(x+z)]
    require(discriminant27(coeffs) == 27*z**4*pencil_discriminant(x, y, z, w, t),
            "四次判别式与三次铅笔的整数恒等式失败")
    # v=1/T；(j:k) 是左接点的齐次坐标。以下运算不对系数点采样。
    ring = fmpz_mpoly_ctx.get(("v", "r", "e", "j", "k", "s"))
    v, r, e, j, k, s = ring.gens()
    def scaled_disc(j, k):
        unit = 1-r*r
        a, b = -r*r*e*k, r*r*j
        c, d = -v*r*unit*e*j, v*r*unit*k
        base = j*(v*unit-r*e)+k*(r-v*unit*e)
        x, y, z, w = v*(b+c), v*(b+d), v*(a-b), v*base+b
        c0 = z*z
        c1 = w*z+v*(w*w-4*x*y-2*(x+y)*z+z*z)
        c2 = x*y-v*w*(x+y)+v*v*(w*w+(x-y)**2-4*x*y-2*(x+y)*z)
        c3 = (x-y)**2
        return (c1*c1*c2*c2-4*v**3*c3*c1**3-4*c0*c2**3
                -27*v**6*c0*c0*c3*c3+18*v**3*c0*c1*c2*c3)
    orders = []
    for pair in ((j, k), (r, -v+v*v*s), (-v+v*v*s, r)):
        expanded = scaled_disc(*pair)
        orders.append(int(min(monomial[0] for monomial in expanded.monoms())))
    require(orders == [10, 14, 14], "齐次判别式在三个坐标图的阶数错误")
    return orders


def verify():
    from sympy import QQ
    from sympy.polys.fields import field
    from flint import fmpz_poly
    orders = identity_and_orders()
    field_d, d = field("D", QQ)
    e, a = field_d(2), field_d(3)
    f = polynomial_value(e, a, d)
    require(f.denom.is_ground and f.numer.degree() == 104, "特征零截面未成为 104 次多项式")
    from math import gcd, lcm
    coefficients = [Q(str(f.numer.get((i,), 0)))/Q(str(f.denom[(0,)])) for i in range(105)]
    common = lcm(*(c.denominator for c in coefficients))
    integers = [int(c*common) for c in coefficients]
    content = gcd(*integers)
    if integers[-1] < 0:
        content = -content
    poly = fmpz_poly([n//content for n in integers])
    resultant = int(poly.resultant(poly.derivative()))
    require(resultant != 0, "截面判别式为零")
    # 将独立的直接数值与所得截面整系数表示的精确值核对。
    for point in (Q(4), Q(7, 2), Q(-3)):
        val = Q(0)
        for c in reversed(coefficients):
            val = val*point+c
        require(val == polynomial_value(Q(2), Q(3), point), "截面与原始圆杆方程不符")
    target = read(HERE / "target.json")
    params = box_values(list(map(Q, target["center"])), Q(target["half_width"]))
    e, a, _, _, d = base_coords(params)
    r, t, right, inner, inner_bar = endpoint(e, a, d)
    h0, h1, h2 = anchor_quadratic(right, inner, inner_bar, r)
    guards = {name: norm_guard(z, name) for name, z in {
        "D": d, "D_plus_r": d+r, "E": e, "A": a, "A_minus_inverse": a-1/a,
        "T": t, "H2": h2, "H_at_minus_E": h0-e*h1+e*e*h2}.items()}
    import hashlib
    return {"valuation_orders_after_scaling_by_v22": orders, "target_units_squared": guards,
            "D_degree_bound": 104, "specialization": {"E": 2, "A": 3},
            "specialized_primitive_coefficients_ascending": list(map(int, poly.coeffs())),
            "specialized_scalar": str(Q(content, common)),
            "specialized_integer_resultant_decimal_digits": len(str(abs(resultant))),
            "specialized_integer_resultant_sha256": hashlib.sha256(str(resultant).encode("ascii")).hexdigest(),
            "specialized_integer_resultant_mod_65521": resultant % 65521,
            "radius_minimal_polynomial_obtained": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        write(args.output, verify())
        print("已验证二参数结式的非零特化与目标分支。")
    finally:
        close()
