"""在全部参数可变时消去最后二次生成元，构造八变量临界系统。"""

import argparse
import time

from descent_support import HERE, SOURCE, AXIS, read, write, sha, bounded
from circuit import Circuit
from chain_fiber import Jet
from geometry import dot


NAMES = ["s", "z", "h", "y", "u", "v", "w", "r"]


def implicit(value, residual, derivative):
    return Jet(value, [-a/derivative for a in residual.gradient])


def closure_coefficients(a, e0, e1, ac, ec, A, B, C, r):
    ka = 1+dot(a, a)-r*r
    kb0, kb1 = -ec*e0, -ac-ec*e1
    nx0, nx1 = ka*e0-kb0*a[1], ka*e1-kb1*a[1]
    ny0, ny1 = a[0]*kb0, a[0]*kb1-ka
    t0, t1 = a[0]*e0, a[0]*e1-a[1]
    U0 = nx0*nx0+ny0*ny0-4*t0*t0
    U1 = 2*(nx0*nx1+ny0*ny1)-8*t0*t1
    U2 = nx1*nx1+ny1*ny1-4*t1*t1
    J0, J1 = U0-U2*C/A, U1-U2*B/A
    F = A*J0*J0-B*J0*J1+C*J1*J1
    return J0, J1, F


def construct():
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(8)]
    s0, z0, h0, y0 = [Jet(value) for value in values[:4]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)]) for i, value in enumerate(values[4:])]
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
    g3 = (x-1)**2+y0*y0-r*r
    y = implicit(values[3], g3, 2*values[3])
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    J0, J1, F = closure_coefficients(c4, e0, e1, ac, ec, A, B, C, r)
    d = -J0/J1
    gd = 2*A*d+B
    multiplier = -J1*gd
    e = e0+e1*d
    return {"schema": "n15-closure-descent-critical8-v1", "variables": NAMES, "nodes": D.nodes,
            "outputs": {"finite_relations": [g.value.index for g in (g0, g1, g2, g3)],
                        "critical_equations": [F.value.index]+[a.index for a in F.gradient[:3]],
                        "radius_derivative": F.gradient[3].index,
                        "closure_coefficients": [J0.value.index, J1.value.index],
                        "right_quadratic": [A.value.index, B.value.index, C.value.index],
                        "recovered_d": d.value.index, "right_simple_derivative": gd.value.index,
                        "critical_multiplier": multiplier.value.index,
                        "geometry": {"c0": [ax.value.index, (-h).value.index], "c3": [a.value.index for a in c3],
                                     "c4": [a.value.index for a in c4], "c5": [d.value.index, e.value.index],
                                     "c6": [x.value.index, y.value.index]}},
            "all_divisor_nodes": sorted(D.divisors),
            "tower_leading_coefficient_nodes": [(L+2*r*dx).value.index, D.constant(1).index, D.constant(1).index, D.constant(1).index],
            "simple_root_derivative_nodes": [gs.value.index, (2*values[1]).index, gh.value.index, (2*values[3]).index],
            "closed_chart_nonzero_guards": [A.value.index, J1.value.index, gd.value.index],
            "critical_locus_nonzero_guards": [F.gradient[3].index],
            "generic_tower_degrees_before_closure": [2]*4, "generic_rank_before_closure": 16,
            "equation_count": 8, "source_circuit_sha256": sha(SOURCE), "source_system_sha256": sha(AXIS),
            "generator_sha256": sha(__file__), "global_zero_dimensionality_proved": False,
            "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "system-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有相对闭环电路")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = construct()
        write(target / "circuit.json.gz", data)
        write(target / "construction.json", {"nodes": len(data["nodes"]), "variables": 8,
                                             "carrier_rank_before_closure": 16,
                                             "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": 8, "carrier_rank_before_closure": 16}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
