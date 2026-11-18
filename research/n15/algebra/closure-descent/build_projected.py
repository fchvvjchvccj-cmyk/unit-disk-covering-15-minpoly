"""以嵌套系数对逐级求相对范数，构造四至七变量的全参数投影电路。"""

import argparse
import time

from descent_support import HERE, SOURCE, AXIS, write, sha, bounded
from build_descent import implicit, closure_coefficients
from quadratic_pair import Pair
from circuit import Circuit
from chain_fiber import Jet
from geometry import dot


def evaluate_known(value, roots):
    if not isinstance(value, Pair):
        return value
    root = roots[id(value.radicand)]
    return evaluate_known(value.a, roots)+evaluate_known(value.b, roots)*root


def construct(keep):
    if keep not in (0, 1, 2, 3):
        raise ValueError("只接受保留零至三个左端生成元")
    names = ["s", "z", "h"][:keep]+["u", "v", "w", "r"]
    D = Circuit(names)
    values = [D.variable(i) for i in range(len(names))]
    p = {name: Jet(value, [D.constant(int(i == j)) for j in range(4)])
         for i, (name, value) in enumerate(zip(["u", "v", "w", "r"], values[keep:]))}
    relations, leading, simple, formal = [], [], [], []

    def lift(root):
        # 每加一层，旧变量都作为新层的常数，避免混淆不同根的系数域。
        for key in p:
            p[key] = root.cast(p[key])

    def root(name, offset, radicand, index, residual=None, partial=None, lead=None):
        if index < keep:
            relations.append(residual.value.index)
            leading.append(lead.value.index)
            simple.append(partial.value.index)
            result = implicit(values[index], residual, partial.value)
        else:
            result = Pair(offset, offset*0+1, radicand)
            formal.append({"name": name, "offset": offset, "radicand": radicand})
            lift(result)
        p[name] = result

    u, v, w, r = [p[name] for name in ("u", "v", "w", "r")]
    dx, L = v-u, (v-u)**2+w*w
    lead = L+2*r*dx
    a_s, b_s = -(L-2*r*dx)/lead, 4*r*w/lead
    if keep:
        s0 = Jet(values[0])
        residual = lead*s0*s0-4*r*w*s0+L-2*r*dx
        partial = 2*lead*s0-4*r*w
    else:
        residual = partial = None
    root("s", b_s/2, a_s+b_s*b_s/4, 0, residual, partial, lead)
    u, r, s = p["u"], p["r"], p["s"]
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    p.update(px=u+aa/2, py=bb/2, aa=aa)
    py = p["py"]
    if keep > 1:
        z0 = Jet(values[1])
        residual, partial = z0*z0+py*py-1, 2*z0
    else:
        residual = partial = None
    root("z", py*0, 1-py*py, 1, residual, partial, py*0+1)
    p["ax"] = (p["px"]+p["z"])/2
    ax, px, py, r = [p[name] for name in ("ax", "px", "py", "r")]
    if keep > 2:
        h0 = Jet(values[2])
        residual, partial = (h0-py)**2+(ax-px)**2-r*r, 2*(h0-py)
    else:
        residual = partial = None
    root("h", py, r*r-(ax-px)**2, 2, residual, partial, py*0+1)
    u, v, w, r, ax, h, px, py, z, aa = [p[name] for name in ("u", "v", "w", "r", "ax", "h", "px", "py", "z", "aa")]
    a = [ax, h]
    T, K = dot(a, a), 1+dot(a, a)-r*r
    c3 = [v+K*ax/T-z-px, w+K*h/T-2*py]
    TT = dot(c3, c3)
    c4 = [v+(1-r*r)*c3[0]/TT-ax, w+(1-r*r)*c3[1]/TT-h]
    x = (1+r+v+(v-u)-aa)/2
    p.update(c0x=ax, c0y=-h, c3x=c3[0], c3y=c3[1], c4x=c4[0], c4y=c4[1], x=x)
    root("y", x*0, r*r-(x-1)**2, 3)
    u, v, w, r, x, y = [p[name] for name in ("u", "v", "w", "r", "x", "y")]
    c4 = [p["c4x"], p["c4y"]]
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    J0, J1, current = closure_coefficients(c4, e0, e1, ac, ec, A, B, C, r)
    operands = []
    while isinstance(current, Pair):
        operands.append(current)
        current = current.norm()
    G = current
    if len(operands) != 4-keep:
        raise ValueError("相对范数的层数不符")
    roots, recovered, coeffs, nonzero, multipliers = {}, {}, {}, [], []
    multiplier = G.cast(1)
    for metadata, operand in zip(formal, reversed(operands)):
        if operand.radicand is not metadata["radicand"]:
            raise ValueError("恢复次序与范数次序不一致")
        c0, c1 = evaluate_known(operand.a, roots), evaluate_known(operand.b, roots)
        t = -c0/c1
        roots[id(metadata["radicand"])] = t
        recovered[metadata["name"]] = evaluate_known(metadata["offset"], roots)+t
        coeffs[metadata["name"]] = [c0.value.index, c1.value.index]
        nonzero.extend([c1.value.index, t.value.index])
        local = -2*c1*t
        multipliers.append(local.value.index)
        multiplier *= local
    at = lambda value: evaluate_known(value, roots)
    Ar, Br, Cr, J0r, J1r = map(at, [A, B, C, J0, J1])
    d = -J0r/J1r
    gd = 2*Ar*d+Br
    multiplier *= -J1r*gd
    e = at(e0)+at(e1)*d
    geometry = {"c0": [at(p["c0x"]), at(p["c0y"])], "c3": [at(p["c3x"]), at(p["c3y"])],
                "c4": [at(p["c4x"]), at(p["c4y"])], "c5": [d, e], "c6": [at(x), recovered["y"]]}
    return {"schema": "n15-relative-norm-projected-v1", "keep_left_roots": keep, "variables": names, "nodes": D.nodes,
            "outputs": {"finite_relations": relations,
                        "critical_equations": [G.value.index]+[a.index for a in G.gradient[:3]],
                        "radius_derivative": G.gradient[3].index,
                        "recovery_coefficients": coeffs,
                        "recovered_left_and_y": {name: a.value.index for name, a in recovered.items()},
                        "recovered_d": d.value.index, "critical_multiplier": multiplier.value.index,
                        "relative_multipliers": multipliers,
                        "geometry": {name: [a.value.index for a in row] for name, row in geometry.items()}},
            "all_divisor_nodes": sorted(D.divisors), "tower_leading_coefficient_nodes": leading,
            "simple_root_derivative_nodes": simple,
            "closed_chart_nonzero_guards": nonzero+[Ar.value.index, J1r.value.index, gd.value.index],
            "critical_locus_nonzero_guards": [G.gradient[3].index],
            "generic_rank_before_closure": 2**keep, "equation_count": keep+4,
            "source_circuit_sha256": sha(SOURCE), "source_system_sha256": sha(AXIS),
            "generator_sha256": sha(__file__), "coefficient_formula_sha256": sha(HERE / "build_descent.py"),
            "pair_arithmetic_sha256": sha(HERE / "quadratic_pair.py"),
            "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "system-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有相对范数投影")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data = construct(args.keep)
        write(target / "circuit.json.gz", data)
        write(target / "construction.json", {"nodes": len(data["nodes"]), "variables": args.keep+4,
                                             "carrier_rank_before_closure": 2**args.keep,
                                             "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"nodes": len(data["nodes"]), "variables": args.keep+4, "seconds": time.perf_counter()-tick}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--keep", type=int, choices=[0, 1, 2, 3], required=True)
    run(parser.parse_args())
