"""独立核验左右端口因式分解，显式说明目标分量选择所需的单位。"""

from context import INPUTS, read


def verify():
    import sympy as S
    from geometry import boundary, dot, subtract
    x, y, u, px, py, r = S.symbols("x y u px py r")
    axis = read(INPUTS / "axis25-system.json.gz")
    ports = S.symbols("c0x c0y c3x c3y c4x c4y c5x c5y c6x c6y u v w r")
    c4, c5, c6, b = ports[4:6], ports[6:8], ports[8:10], ports[11:13]
    rr, uu, vv, ww = ports[13], ports[10], ports[11], ports[12]
    zz = 2*c6[0]-1-rr

    def triangle(a, b, c):
        d, e = subtract(b, a), subtract(c, a)
        cross = d[0]*e[1]-d[1]*e[0]
        return dot(d, d)*dot(e, e)*dot(subtract(b, c), subtract(b, c))-4*rr*rr*cross*cross

    norm6 = dot(c6, c6)
    expected_rows = {0: (c6[0]-1)**2+c6[1]**2-rr**2,
                     3: boundary(c4, c5, rr),
                     4: norm6*(dot(c5, c5)-norm6+2*c6[0])-2*c5[0]*(c6[0]**2-c6[1]**2)-4*c5[1]*c6[0]*c6[1],
                     8: triangle(c5, b, [zz, 0]),
                     9: ((uu-vv)**2+ww**2)*((vv-zz)**2+ww**2)-4*rr*rr*ww*ww}
    for index, expression in expected_rows.items():
        original = sum(term["coefficient"]*S.prod(ports[i]**e for i, e in term["powers"]) for term in axis["equations"][index])
        if S.Poly(expression-original, *ports) != 0:
            raise ValueError("剩余原几何约束不匹配："+str(index))
    row = axis["equations"][10]
    H = sum(term["coefficient"]*S.prod({0: x, 1: y, 10: u, 13: r}[i]**e for i, e in term["powers"]) for term in row)
    radical = (px-u)**2+py*py
    terms = S.Poly(H, r)
    if any(e[0] % 2 for e, _ in terms.terms()):
        raise ValueError("左端原约束含未处理的奇次半径")
    H = sum(c*radical**(e[0]//2) for e, c in terms.terms())
    circle = y*y-2*py*y+x*x-2*px*x+2*u*px-u*u
    quotient, remainder = S.div(H, circle, y)
    G = (2*x-px)**2+py*py-1
    C1 = px*u-px*x-py*y-u*u+u*x
    C2 = -px*px+4*px*u-4*px*x-py*py-2*u*u+2*u*x+1
    if S.expand(H-quotient*circle-2*G*C1*C2) != 0:
        raise ValueError("左端八分量图表的整数因式分解失败")
    d, e, v, w, Z = S.symbols("d e v w Z")
    edge = (v-Z)**2+w*w
    E = (u-v)**2+w*w
    cross = (v-d)*(-e)-(w-e)*(Z-d)
    triangle = ((v-d)**2+(w-e)**2)*edge*((d-Z)**2+e*e)-4*r*r*cross*cross
    plus = w*(d*d+e*e-(Z+2*v-u)*d+Z*(2*v-u))-(w*w-(v-Z)*(v-u))*e
    minus = w*(d*d+e*e-(Z+u)*d+Z*u)-(w*w+(v-Z)*(v-u))*e
    axis_circle = E*edge-4*r*r*w*w
    if S.expand(w*w*triangle-edge*plus*minus-axis_circle*cross*cross) != 0:
        raise ValueError("右端两圆分量的整数因式分解失败")
    aa, bb, dx, t = S.symbols("aa bb dx t")
    L = dx*dx+w*w
    line = dx*aa+w*bb-L
    circ = aa*aa+bb*bb-4*r*r
    if S.expand((dx*bb-w*aa)**2-L*(4*r*r-L)-L*circ+line*(line+2*L)) != 0:
        raise ValueError("交点叉积基恒等式失败")
    shared = L*(dx-aa)**2-w*w*(4*r*r-L)
    if S.expand(shared-w*w*circ-line*(dx*aa-w*bb-L)) != 0:
        raise ValueError("左右共享二次根的无半角恒等式失败")
    X, Y, D, E = S.symbols("X Y D E")
    N, beta = X*X+Y*Y, r*r-1
    Sx = 2*X+beta
    H0 = N-Sx
    H4 = N*(D*D+E*E-N+2*X)-2*D*(X*X-Y*Y)-4*E*X*Y
    circle_right = Sx*(D*D+E*E+2*D-beta)-4*D*X*X-4*E*X*Y
    if S.expand(circle_right-H4-H0*(N+beta-D*D-E*E-2*D)) != 0:
        raise ValueError("右端原第四行与圆方程的理想身份失败")
    return {"left_factorization": "H10|r^2=(px-u)^2+py^2 = quotient*circle + 2*G*C1*C2",
            "remaining_source_geometric_rows_checked": sorted(expected_rows),
            "right_circle_ideal_identity": "Sx*Cright=H4+H0*(N+beta-d^2-e^2-2*d)",
            "left_circle": "circle=y^2-2*py*y+x^2-2*px*x+2*u*px-u^2",
            "left_G": "G=(2*x-px)^2+py^2-1", "left_C1": str(C1), "left_C2": str(C2),
            "left_quotient_variables": ["x", "y", "u", "px", "py"],
            "left_quotient_terms": [[list(exponent), str(c)] for exponent, c in S.Poly(quotient, x, y, u, px, py).terms()],
            "right_factorization": "w^2*T(c5,b,c13)=edge*plus*minus+H9*cross^2",
            "right_plus": str(plus), "right_minus": str(minus),
            "cross_product_basis_identity": True, "shared_root_identity_without_half_angle": True,
            "scope": "通用特征零多项式恒等式；因子选择另由原根整盒的单位界支持。"}
