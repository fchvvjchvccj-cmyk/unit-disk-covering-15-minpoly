"""继续消去右端的两个二次生成元，保留全参数七变量临界电路。"""

import argparse
import time

from descent_support import HERE, SOURCE, AXIS, write, sha, bounded
from build_descent import implicit, closure_coefficients
from quadratic_pair import Pair
from circuit import Circuit
from chain_fiber import Jet
from geometry import dot


NAMES = ["s", "z", "h", "u", "v", "w", "r"]


def right_pairs(c4, x, u, v, w, r):
    V = r*r-(x-1)**2
    Y = Pair(x*0, x*0+1, V)
    Q = Y.cast
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = Q(2-4*x*x/Sx), Q(-4*x/Sx)*Y
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = Q(bc)-ac, Q(be)-ec, Q(bc0+beta)
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-Q(beta)
    J0, J1, F = closure_coefficients([Q(a) for a in c4], e0, e1, ac, ec, A, B, C, Q(r))
    return V, e0, e1, A, B, C, J0, J1, F


def construct():
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(7)]
    s0, z0, h0 = [Jet(value) for value in values[:3]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)]) for i, value in enumerate(values[3:])]
    dx, L = v-u, (v-u)**2+w*w
    g0 = (L+2*r*dx)*s0*s0-4*r*w*s0+L-2*r*dx
    gs = 2*(L+2*r*dx)*s0-4*r*w
    s = implicit(values[0], g0, gs.value)
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    px, py = u+aa/2, bb/2
    g1 = z0*z0+py*py-1
    z = implicit(values[1], g1, 2*values[1])
    ax = (px+z)/2
    g2 = (h0-py)**2+(ax-px)**2-r*r
    gh = 2*(h0-py)
    h = implicit(values[2], g2, gh.value)
    a = [ax, h]
    T, K = dot(a, a), 1+dot(a, a)-r*r
    c3 = [v+K*ax/T-z-px, w+K*h/T-2*py]
    TT = dot(c3, c3)
    c4 = [v+(1-r*r)*c3[0]/TT-ax, w+(1-r*r)*c3[1]/TT-h]
    x = (1+r+v+dx-aa)/2
    V, e0, e1, A, B, C, J0, J1, F = right_pairs(c4, x, u, v, w, r)
    F0, F1 = F.a, F.b
    G = F0*F0-V*F1*F1
    y = -F0/F1
    J0r, J1r = J0.evaluate(y), J1.evaluate(y)
    Ar, Br, Cr = A.evaluate(y), B.evaluate(y), C.evaluate(y)
    d = -J0r/J1r
    gd = 2*Ar*d+Br
    e = e0.evaluate(y)+e1.evaluate(y)*d
    multiplier_y = -2*y*F1
    multiplier = multiplier_y*(-J1r*gd)
    pack = lambda q: [q.a.value.index, q.b.value.index]
    return {"schema": "n15-closure-descent-critical7-v1", "variables": NAMES, "nodes": D.nodes,
            "outputs": {"finite_relations": [g.value.index for g in (g0, g1, g2)],
                        "critical_equations": [G.value.index]+[a.index for a in G.gradient[:3]],
                        "radius_derivative": G.gradient[3].index,
                        "first_norm_coefficients": [F0.value.index, F1.value.index],
                        "right_pairs": {name: pack(q) for name, q in zip(["e0", "e1", "A", "B", "C", "J0", "J1"], [e0, e1, A, B, C, J0, J1])},
                        "right_radicand": V.value.index,
                        "recovered_y": y.value.index, "recovered_d": d.value.index,
                        "right_simple_derivative": gd.value.index,
                        "critical_multiplier_y": multiplier_y.value.index, "critical_multiplier": multiplier.value.index,
                        "geometry": {"c0": [ax.value.index, (-h).value.index], "c3": [a.value.index for a in c3],
                                     "c4": [a.value.index for a in c4], "c5": [d.value.index, e.value.index],
                                     "c6": [x.value.index, y.value.index]}},
            "all_divisor_nodes": sorted(D.divisors),
            "tower_leading_coefficient_nodes": [(L+2*r*dx).value.index, D.constant(1).index, D.constant(1).index],
            "simple_root_derivative_nodes": [gs.value.index, (2*values[1]).index, gh.value.index],
            "closed_chart_nonzero_guards": [V.value.index, F1.value.index, Ar.value.index, J1r.value.index, gd.value.index],
            "critical_locus_nonzero_guards": [G.gradient[3].index],
            "generic_tower_degrees_before_closure": [2]*3, "generic_rank_before_closure": 8,
            "equation_count": 7, "source_circuit_sha256": sha(SOURCE), "source_system_sha256": sha(AXIS),
            "generator_sha256": sha(__file__), "coefficient_formula_sha256": sha(HERE / "build_descent.py"),
            "pair_arithmetic_sha256": sha(HERE / "quadratic_pair.py"),
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "system-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有二重闭环消元电路")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = construct()
        write(target / "circuit.json.gz", data)
        write(target / "construction.json", {"nodes": len(data["nodes"]), "variables": 7,
                                             "carrier_rank_before_closure": 8,
                                             "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": 7, "carrier_rank_before_closure": 8}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
