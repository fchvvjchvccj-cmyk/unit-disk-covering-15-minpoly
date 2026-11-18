"""消去重复二次根，构造九变量、秩32的目标临界电路。"""

import argparse
import time

from bridge_support import HERE, AXIS, FACTOR, read, write, sha, bounded
from circuit import Circuit
from chain_fiber import Jet
from geometry import dot, boundary


NAMES = ["s", "z", "h", "y", "d", "u", "v", "w", "r"]


def implicit(value, residual, derivative):
    return Jet(value, [-a/derivative for a in residual.gradient])


def construct(sign=1):
    if sign not in (-1, 1):
        raise ValueError("共享根符号必须为正负一")
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(9)]
    s0, z0, h0, y0, d0 = [Jet(value) for value in values[:5]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)]) for i, value in enumerate(values[5:])]
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
    x = (1+r+v+sign*(dx-aa))/2
    g3 = (x-1)**2+y0*y0-r*r
    y = implicit(values[3], g3, 2*values[3])
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    g4 = A*d0*d0+B*d0+C
    gd = 2*A*d0+B
    d = implicit(values[4], g4, gd.value)
    e = e0+e1*d
    H = boundary(c4, [d, e], r)
    return {"schema": "n15-shared-root-critical9-v1", "variables": NAMES, "selected_sign": sign, "nodes": D.nodes,
            "outputs": {"finite_relations": [p.value.index for p in (g0, g1, g2, g3, g4)],
                        "critical_equations": [H.value.index]+[a.index for a in H.gradient[:3]],
                        "radius_derivative": H.gradient[3].index,
                        "geometry": {"c0": [ax.value.index, (-h).value.index], "c3": [a.value.index for a in c3],
                                     "c4": [a.value.index for a in c4], "c5": [d.value.index, e.value.index],
                                     "c6": [x.value.index, y.value.index]}},
            "all_divisor_nodes": sorted(D.divisors),
            "tower_leading_coefficient_nodes": [(L+2*r*dx).value.index, D.constant(1).index, D.constant(1).index,
                                                D.constant(1).index, A.value.index],
            "simple_root_derivative_nodes": [gs.value.index, (2*values[1]).index, gh.value.index,
                                             (2*values[3]).index, gd.value.index],
            "critical_locus_nonzero_guards": [H.gradient[3].index],
            "generic_tower_degrees": [2]*5, "generic_rank_before_closure": 32,
            "equation_count": 9, "source_system_sha256": sha(AXIS), "left_factor_sha256": sha(FACTOR),
            "generator_sha256": sha(__file__), "global_zero_dimensionality_proved": False,
            "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "system-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有共享根电路")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = construct(args.sign)
        write(target / "circuit.json.gz", data)
        write(target / "construction.json", {"nodes": len(data["nodes"]), "variables": 9, "equations": 9,
                                             "rank": 32, "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": 9, "rank": 32}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--sign", type=int, default=1)
    run(parser.parse_args())
