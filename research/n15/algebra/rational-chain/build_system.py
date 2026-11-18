"""构造五级有限载体及四条临界条件，得到九变量九方程有理电路。"""

import argparse
import time

from support import HERE, AXIS, FACTOR, read, write, sha, bounded
from circuit import Circuit
from chain_fiber import Jet
from geometry import reconstruct, boundary


NAMES = ["s", "k", "x", "y", "d", "u", "v", "w", "r"]


def polynomial(rows, ports):
    value = ports[0].cast(0)
    for row in rows:
        term = ports[0].cast(row["coefficient"])
        for index, exponent in row["powers"]:
            term *= ports[index]**exponent
        value += term
    return value


def construct():
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(9)]
    s0, k0, x0, y0, d0 = [Jet(v) for v in values[:5]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)]) for i, value in enumerate(values[5:])]
    L, dx = (v-u)**2+w*w, v-u
    g0 = (L+2*r*dx)*s0*s0-4*r*w*s0+L-2*r*dx
    gs = 2*(L+2*r*dx)*s0-4*r*w
    s = Jet(values[0], [-a/gs.value for a in g0.gradient])
    factor = read(FACTOR)
    rows = factor["factor_polynomials"][factor["target_selected_factor_index"]]
    g1 = polynomial(rows, {0: k0, 1: s, 2: u, 3: r})
    differentiated = [dict(row, coefficient=row["coefficient"]*dict(row["powers"])[0],
                           powers=[[j, e-int(j == 0)] for j, e in row["powers"] if e-int(j == 0)])
                      for row in rows if dict(row["powers"]).get(0, 0)]
    gk = polynomial(differentiated, {0: k0, 1: s, 2: u, 3: r})
    k = Jet(values[1], [-a/gk.value for a in g1.gradient])
    E = (u-v)**2+w*w
    h = 2*x0-1-r-v
    g2 = E*h*h-w*w*(4*r*r-E)
    gx = 4*E*h
    x = Jet(values[2], [-a/gx.value for a in g2.gradient])
    g3 = (x-1)**2+y0*y0-r*r
    y = Jet(values[3], [-a/(2*values[3]) for a in g3.gradient])
    beta = r*r-1
    Sx = 2*x+beta
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z = 2*x-1-r
    Pz = (v-Z)*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    g4 = A*d0*d0+B*d0+C
    gd = 2*A*d0+B
    d = Jet(values[4], [-a/gd.value for a in g4.gradient])
    e = e0+e1*d
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    X = (aa-bb*k)/(1+k*k)
    c0 = [u+X, k*X]
    geometry = reconstruct(c0, u, [v, w], r)
    H = boundary(geometry["c4"], [d, e], r)
    leading_k = [dict(row, powers=[[j, e] for j, e in row["powers"] if j != 0])
                 for row in rows if dict(row["powers"]).get(0, 0) == 4]
    lk = polynomial(leading_k, {0: k0, 1: s, 2: u, 3: r})
    outputs = {"finite_relations": [p.value.index for p in (g0, g1, g2, g3, g4)],
               "critical_equations": [H.value.index]+[a.index for a in H.gradient[:3]],
               "radius_derivative": H.gradient[3].index,
               "geometry": {"c0": [a.value.index for a in c0],
                            "c3": [a.value.index for a in geometry["c3"]],
                            "c4": [a.value.index for a in geometry["c4"]],
                            "c5": [d.value.index, e.value.index], "c6": [x.value.index, y.value.index]}}
    return {"schema": "n15-rational-chain-critical9-v1", "variables": NAMES, "nodes": D.nodes,
            "outputs": outputs, "all_divisor_nodes": sorted(D.divisors),
            "tower_leading_coefficient_nodes": [(L+2*r*dx).value.index, lk.value.index, (4*E).value.index, D.constant(1).index, A.value.index],
            "simple_root_derivative_nodes": [gs.value.index, gk.value.index, gx.value.index, (2*values[3]).index, gd.value.index],
            "critical_locus_nonzero_guards": [H.gradient[3].index],
            "generic_tower_degrees": [2, 4, 2, 2, 2], "generic_rank_before_closure": 64,
            "equation_count": 9, "equation_semantics": "G0..G4=0,H=Hu=Hv=Hw=0；各偏导沿五级隐函数有限载体取全导数，保持r为独立参数。",
            "localization": "电路所有分母、五个单根导数与各载体首项非零；目标分支另依赖已认证互补交点图表。",
            "source_system_sha256": sha(AXIS), "left_factor_sha256": sha(FACTOR),
            "generator_sha256": sha(__file__), "geometry_sha256": sha(HERE / "geometry.py"),
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    source = HERE / "system-runs" / args.tag
    if source.exists():
        raise FileExistsError("拒绝覆盖已有九变量系统")
    source.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = construct()
        write(source / "circuit.json.gz", data)
        write(source / "construction.json", {"node_count": len(data["nodes"]),
                                             "division_guard_count": len(data["all_divisor_nodes"]),
                                             "variables": data["variables"], "equations": data["equation_count"],
                                             "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": 9, "equations": 9}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
