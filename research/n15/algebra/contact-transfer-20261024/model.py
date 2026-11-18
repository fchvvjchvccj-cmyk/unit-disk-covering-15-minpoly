"""以交点叉积基和本原闭环构造八变量、秩十六的临界系统。"""

import hashlib

from context import INPUTS, RESULTS, payload, sha, write, bounded
from circuit import Circuit
from chain_fiber import Jet

NAMES = ["t", "z", "h", "y", "u", "v", "w", "r"]


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]


def implicit(value, equation, partial):
    return Jet(value, [-a/partial for a in equation.gradient])


def construct():
    D = Circuit(NAMES)
    values = [D.variable(i) for i in range(8)]
    t0, z0, h0, y0 = [Jet(a) for a in values[:4]]
    u, v, w, r = [Jet(value, [D.constant(int(i == j)) for j in range(4)])
                  for i, value in enumerate(values[4:])]
    dx, L = v-u, (v-u)**2+w*w
    g0 = t0*t0-L*(4*r*r-L)
    t = implicit(values[0], g0, 2*values[0])
    aa, bb = dx-w*t/L, w+dx*t/L
    px, py = u+aa/2, bb/2
    g1 = z0*z0+py*py-1
    z = implicit(values[1], g1, 2*values[1])
    ax = (px+z)/2
    g2 = (h0-py)**2+(ax-px)**2-r*r
    h = implicit(values[2], g2, 2*(values[2]-py.value))
    a = [ax, h]
    S, K = dot(a, a), 1+dot(a, a)-r*r
    c3 = [v+K*ax/S-z-px, w+K*h/S-2*py]
    T = dot(c3, c3)
    c4 = [v+(1-r*r)*c3[0]/T-ax, w+(1-r*r)*c3[1]/T-h]
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
    ka = dot(c4, c4)-beta
    aq = -(ac*c4[0]+ec*c4[1])/2
    q0, q1 = -ec*e0/2, -(ac+ec*e1)/2
    ab0, ab1 = c4[1]*e0, c4[0]+c4[1]*e1
    U0 = (2*q0+beta)*(ka+2*aq)-4*q0*ab0
    U1 = 2*q1*(ka+2*aq)-4*(q0*ab1+q1*ab0)
    U2 = -4*q1*ab1
    J0, J1 = U0-U2*C/A, U1-U2*B/A
    F = A*J0*J0-B*J0*J1+C*J1*J1
    d = -J0/J1
    e = e0+e1*d
    gd = 2*A*d+B
    removed = ka-2*aq
    multiplier = -J1*gd/removed
    return {"schema": "n15-transfer-primitive8-v1", "variables": NAMES, "nodes": D.nodes,
            "outputs": {"finite_relations": [g.value.index for g in (g0, g1, g2, g3)],
                        "critical_equations": [F.value.index]+[a.index for a in F.gradient[:3]],
                        "radius_derivative": F.gradient[3].index,
                        "right_quadratic": [q.value.index for q in (A, B, C)],
                        "closure_coefficients": [J0.value.index, J1.value.index],
                        "removed_factor": removed.value.index, "critical_multiplier": multiplier.value.index,
                        "right_simple_derivative": gd.value.index,
                        "geometry": {"c0": [ax.value.index, (-h).value.index],
                                     "c3": [q.value.index for q in c3], "c4": [q.value.index for q in c4],
                                     "c5": [d.value.index, e.value.index], "c6": [x.value.index, y.value.index]},
                        "intersection": [px.value.index, py.value.index]},
            "all_divisor_nodes": sorted(D.divisors),
            "simple_root_derivative_nodes": [2*values[0], 2*values[1], 2*(values[2]-py.value), 2*values[3]],
            "closed_chart_nonzero_guards": [values[7].index, A.value.index, J1.value.index, gd.value.index, removed.value.index],
            "critical_locus_nonzero_guards": [F.gradient[3].index],
            "generic_tower_degrees_before_closure": [2]*4, "generic_rank_before_closure": 16,
            "equation_count": 8, "source_system_payload_sha256": hashlib.sha256(payload(INPUTS / "axis25-system.json.gz")).hexdigest(),
            "generator_sha256": sha(__file__), "global_zero_dimensionality_proved": False,
            "positive_radius_required_for_target": True,
            "radius_minimal_polynomial_obtained": False}


def build():
    data = construct()
    data["simple_root_derivative_nodes"] = [a.index for a in data["simple_root_derivative_nodes"]]
    return data


if __name__ == "__main__":
    with bounded():
        data = build()
        write(RESULTS / "primitive16.json.gz", data)
        print({"variables": 8, "carrier_rank": 16, "nodes": len(data["nodes"])})
