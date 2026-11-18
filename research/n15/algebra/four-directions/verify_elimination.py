"""用有理区间验证四次重根、去掉辅助接点后的局部坐标及一条特征零截面。"""

import argparse
from fractions import Fraction as Q

from chart import HERE, circle_closure, direction, forward, read
from context import write
from eliminate import anchor_quadratic, discriminant27, eliminant, endpoint, rod_quartic
from verify import away, box_values, independent_variables, require
from verify_chart import certify_system


class Complex:
    """实部和虚部可取有理区间或带导数的区间，避免使用浮点复数判零。"""

    def __init__(self, x, y):
        self.x, self.y = x, y

    def cast(self, other):
        return other if isinstance(other, Complex) else Complex(self.x*0+other, self.x*0)

    def __add__(self, other):
        other = self.cast(other)
        return Complex(self.x+other.x, self.y+other.y)

    __radd__ = __add__

    def __neg__(self):
        return Complex(-self.x, -self.y)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = self.cast(other)
        return Complex(self.x*other.x-self.y*other.y, self.x*other.y+self.y*other.x)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        denominator = other.x*other.x+other.y*other.y
        return Complex((self.x*other.x+self.y*other.y)/denominator,
                       (self.y*other.x-self.x*other.y)/denominator)

    def __rtruediv__(self, other):
        return self.cast(other)/self


def value(polynomial, point):
    result = point*0
    for coefficient in reversed(polynomial):
        result = point*result+coefficient
    return result


def derivative(poly):
    return [i*coefficient for i, coefficient in enumerate(poly) if i]


def norm_guard(z, label):
    return away(z.x*z.x+z.y*z.y, label)


def base_coords(parameters):
    a, b, c, d = [Complex(*direction(s)) for s in parameters]
    _, z, _, _, _, _, _, r = forward(parameters)
    e = Complex(z, r*a.y)
    return e, a, b, c, d


def three_system(parameters):
    from verify import independent_variables
    e, a, d = [Complex(*direction(s)) for s in independent_variables(parameters)]
    radius = endpoint(e, a, d)[0].x
    k = eliminant(e, a, d).x
    return [k.value, k.gradient[2], k.gradient[0]*radius.gradient[1]-k.gradient[1]*radius.gradient[0]]


def verify_projection():
    target = read(HERE / "target.json")
    center = list(map(Q, target["center"]))
    width = Q(target["half_width"])
    params = box_values(center, width)
    e, a, b, c, d = base_coords(params)
    r, transfer, right, inner, inner_bar = endpoint(e, a, d)
    left = Complex(*circle_closure(params)[1][0])
    h = anchor_quadratic(right, inner, inner_bar, r)
    q = rod_quartic(left, e, r, transfer)
    other_left = -h[1]/h[2]-left
    scales = {
        "r": r, "C_minus_1_over_C": c-1/c, "r2_T_times_E_plus_left": r*r*transfer*(e+left),
        "H2": h[2], "H_left": value(derivative(h), left),
        "quartic_leading": q[4], "quartic_BB": value(derivative(derivative(q)), b),
        "discriminant_at_other_anchor": discriminant27(rod_quartic(other_left, e, r, transfer))}
    # 用精确多项式的一阶微分检验 Q_L 与判别式对 L 的导数。
    from verify import Dual
    def variable(z, moving):
        return Complex(Dual(z.x, [z.x*0+int(moving)]), Dual(z.y, [z.y*0]))
    def tangent(z):
        return Complex(z.x.gradient[0], z.y.gradient[0])
    qj = rod_quartic(variable(left, True), variable(e, False), variable(r, False), variable(transfer, False))
    scales["quartic_L"] = tangent(value(qj, variable(b, False)))
    scales["discriminant_L"] = tangent(discriminant27(qj))
    guards = {name: norm_guard(z, name) for name, z in scales.items()}

    def half_e(p):
        e = base_coords(p)[0]
        return e.y/(1+e.x)
    mid = [half_e(center), center[0], center[3]]
    scale = 10**55
    mid = [Q(round(x*scale), scale) for x in mid]
    box3 = Q(1, 10**37)
    from verify import encode
    for x, m in zip([half_e(params), params[0], params[3]], mid):
        lo, hi = map(Q, encode(x))
        require(m-box3 < lo <= hi < m+box3, "四变量根未映入三变量盒")
    p = independent_variables(box_values(mid, box3))
    ep, ap, dp = [Complex(*direction(s)) for s in p]
    k = eliminant(ep, ap, dp).x
    gradient_e = away(k.gradient[0], "实部对第一参数的偏导")
    p2 = independent_variables(p)
    ep, ap, dp = [Complex(*direction(s)) for s in p2]
    kdd = away(eliminant(ep, ap, dp).x.gradient[2].gradient[2], "实部对第三参数的二阶偏导")
    return {"nonzero_squared_norms": guards, "parameters": ["s_E", "s_A", "s_D"],
            "center": list(map(str, mid)), "half_width": str(box3), "source_image_inside": True,
            "real_K_E": gradient_e, "real_K_DD": kdd,
            "critical_root": certify_system(three_system, mid, box3)}


def characteristic_zero_slice():
    from sympy import QQ
    from sympy.polys.fields import field
    ring, d = field("D", QQ)
    e, a = ring(2), ring(3)
    k = eliminant(e, a, d)
    # 特化在复化参数空间中进行；固定的两个有理数不代表目标方向。
    numerator, denominator = k.numer, k.denom
    content, layers = numerator.sqf_list()
    product = numerator.ring(content)
    for factor, multiplicity in layers:
        require(factor.gcd(factor.diff(0)).degree() == 0, "截面因子含重根")
        product *= factor**multiplicity
    require(product == numerator and numerator.gcd(denominator).degree() == 0, "截面回乘失败")
    linear_layer = next(f for f, m in layers if m == 1)
    return {"specialization": {"E": 2, "A": 3}, "numerator_degree": numerator.degree(),
            "denominator_degree": denominator.degree(), "squarefree_layers": [[f.degree(), m] for f, m in layers],
            "multiplicity_one_coefficients_ascending": [str(linear_layer.get((i,), 0)) for i in range(linear_layer.degree()+1)],
            "radius_minimal_polynomial_obtained": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        result = verify_projection()
        result["slice"] = characteristic_zero_slice()
        write(args.output, result)
        print("三变量消元、重根图表与特征零截面已验证。")
    finally:
        close()
