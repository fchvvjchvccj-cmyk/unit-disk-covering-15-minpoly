"""把左端四次条件拆成两个几何二次层，构造十变量全二次载体。"""

import argparse
import time

from support import HERE, FACTOR, AXIS, read, write, sha, bounded
from circuit import Circuit
from chain_fiber import Jet
from geometry import dot, boundary


NAMES = ["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"]


def prove_left_identity():
    import sympy as S
    s, k, u, r = S.symbols("s k u r")
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    X = (aa-bb*k)/(1+k*k)
    ax, ay = u+X, -k*X
    t, T = 2*ax-u, ax*ax+ay*ay
    qx = (t*t+1-r*r)/(2*t)
    expression = S.together((1+T-r*r-2*ax*qx)/(2*ay)-bb/2)
    numerator, denominator = expression.as_numer_denom()
    factor = read(FACTOR)
    ports = {0: k, 1: s, 2: u, 3: r}
    f4 = sum(row["coefficient"]*S.prod(ports[j]**e for j, e in row["powers"])
             for row in factor["factor_polynomials"][factor["target_selected_factor_index"]])
    expected = -2*r*(2*k*s+s*s-1)*f4
    if not S.Poly(numerator-expected, s, k, u, r).is_zero:
        raise ValueError("四次因子与共同纵坐标的整数恒等式失败")
    encode = lambda expression: [{"coefficient": int(c), "powers": [[i, e] for i, e in enumerate(exponents) if e]}
                                  for exponents, c in S.Poly(expression, s, k, u, r).terms()]
    return {"verified": True, "variables": ["s", "k", "u", "r"],
            "identity": "q0y-bb/2=-2*r*(2*k*s+s^2-1)*f4/denominator",
            "denominator": encode(denominator), "left_factor_sha256": sha(FACTOR),
            "argument": "q0为单位圆与旧圆的交点，p0=(u+aa/2,bb/2)。四次目标分量给出q0y=p0y。等半径圆的公共弦中垂线随后给出a_x=(q0x+p0x)/2；a_y由到p0距离r的二次式给出。反向映射在原四次图表的单位非零处恢复f4=0。",
            "left_tower_degrees": [2, 2, 2], "left_rank_before_closure": 8}


def implicit(value, residual, derivative):
    return Jet(value, [-a/derivative for a in residual.gradient])


def construct():
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(10)]
    s0, z0, h0, x0, y0, d0 = [Jet(value) for value in values[:6]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)]) for i, value in enumerate(values[6:])]
    dx, L = v-u, (v-u)**2+w*w
    g0 = (L+2*r*dx)*s0*s0-4*r*w*s0+L-2*r*dx
    gs = 2*(L+2*r*dx)*s0-4*r*w
    s = implicit(values[0], g0, gs.value)
    px, py = u+r*(1-s*s)/(1+s*s), 2*r*s/(1+s*s)
    g1 = z0*z0+py*py-1
    z = implicit(values[1], g1, 2*values[1])
    ax = (px+z)/2
    g2 = (h0-py)**2+(ax-px)**2-r*r
    gh = 2*(h0-py)
    h = implicit(values[2], g2, gh.value)
    a = [ax, h]
    T = dot(a, a)
    K = 1+T-r*r
    c3 = [v+K*ax/T-z-px, w+K*h/T-2*py]
    TT = dot(c3, c3)
    c4 = [v+(1-r*r)*c3[0]/TT-ax, w+(1-r*r)*c3[1]/TT-h]
    E = (u-v)**2+w*w
    hright = 2*x0-1-r-v
    g3 = E*hright*hright-w*w*(4*r*r-E)
    gx = 4*E*hright
    x = implicit(values[3], g3, gx.value)
    g4 = (x-1)**2+y0*y0-r*r
    y = implicit(values[4], g4, 2*values[4])
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    g5 = A*d0*d0+B*d0+C
    gd = 2*A*d0+B
    d = implicit(values[5], g5, gd.value)
    e = e0+e1*d
    H = boundary(c4, [d, e], r)
    return {"schema": "n15-quadratic-chain-critical10-v1", "variables": NAMES, "nodes": D.nodes,
            "outputs": {"finite_relations": [p.value.index for p in (g0, g1, g2, g3, g4, g5)],
                        "critical_equations": [H.value.index]+[a.index for a in H.gradient[:3]],
                        "radius_derivative": H.gradient[3].index,
                        "geometry": {"c0": [ax.value.index, (-h).value.index],
                                     "c3": [a.value.index for a in c3], "c4": [a.value.index for a in c4],
                                     "c5": [d.value.index, e.value.index], "c6": [x.value.index, y.value.index]}},
            "all_divisor_nodes": sorted(D.divisors),
            "tower_leading_coefficient_nodes": [(L+2*r*dx).value.index, D.constant(1).index, D.constant(1).index,
                                                (4*E).value.index, D.constant(1).index, A.value.index],
            "simple_root_derivative_nodes": [gs.value.index, (2*values[1]).index, gh.value.index,
                                             gx.value.index, (2*values[4]).index, gd.value.index],
            "critical_locus_nonzero_guards": [H.gradient[3].index],
            "generic_tower_degrees": [2]*6, "generic_rank_before_closure": 64,
            "equation_count": 10, "equation_semantics": "G0..G5=0,H=Hu=Hv=Hw=0，偏导沿六级二次塔取全导数。",
            "source_system_sha256": sha(AXIS), "left_factor_sha256": sha(FACTOR),
            "generator_sha256": sha(__file__), "geometry_sha256": sha(HERE / "geometry.py"),
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "system-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有全二次系统")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        proof = prove_left_identity()
        data = construct()
        write(target / "left-identity.json", proof)
        write(target / "circuit.json.gz", data)
        write(target / "construction.json", {"node_count": len(data["nodes"]), "variables": data["variables"],
                                             "equations": 10, "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": 10, "equations": 10}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
