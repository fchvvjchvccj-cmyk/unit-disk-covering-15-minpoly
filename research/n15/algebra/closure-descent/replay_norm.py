"""以旧32维规范形的范数和迹独立重验四参数电路的符号截面及四个导数。"""

import argparse
import time

from descent_support import HERE, SHARED, SOURCE, read, write, sha, bounded
from verify_quotient import setup, load, vector, decode_scalar
from verify_multivariate import prune
from verify_system import execute
from tower import Rational, Element
from chain_fiber import Jet


def generators(algebra, count):
    result = []
    for j in range(count):
        coefficients = [algebra.zero]*2**(j+1)
        coefficients[2**j] = algebra.one
        result.append(Element(algebra, j+1, tuple(coefficients)))
    return result


def right_leading(algebra, parameters):
    s0, z0, h0, y0 = map(Jet, generators(algebra, 4))
    u, v, w, r = [Jet(algebra.scalar(value), [algebra.scalar(int(i == j)) for j in range(4)])
                   for i, value in enumerate(parameters)]
    dx, L = v-u, (v-u)**2+w*w
    g = (L+2*r*dx)*s0*s0-4*r*w*s0+L-2*r*dx
    gs = 2*(L+2*r*dx)*s0-4*r*w
    s = Jet(s0.value, [-a/gs.value for a in g.gradient])
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    px, py = u+aa/2, bb/2
    g = z0*z0+py*py-1
    z = Jet(z0.value, [-a/(2*z0.value) for a in g.gradient])
    ax = (px+z)/2
    g = (h0-py)**2+(ax-px)**2-r*r
    h = Jet(h0.value, [-a/(2*(h0-py).value) for a in g.gradient])
    x = (1+r+v+dx-aa)/2
    g = (x-1)**2+y0*y0-r*r
    y = Jet(y0.value, [-a/(2*y0.value) for a in g.gradient])
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be = -(Z+2*v-u), -(w*w-Pz)/w
    ld, le = bc-ac, be-ec
    e1 = -ld/le
    return 1+e1*e1


def trace(algebra, coefficients, level):
    if level == 0:
        return coefficients[0]
    half = len(coefficients)//2
    a, b = coefficients[:half], coefficients[half:]
    relation_b = algebra.relations[level-1][1]
    lower = algebra.add(algebra.add(a, a), algebra.multiply(b, relation_b, level-1))
    return trace(algebra, lower, level-1)


def same(a, b):
    return a.n == b.n and a.d == b.d


def compare(actual, expected):
    if len(actual) != 5 or len(expected) != 5:
        raise ValueError("闭环值及四个导数的个数错误")
    for name, a, b in zip(["G", "Gu", "Gv", "Gw", "Gr"], actual, expected):
        if not same(a, b):
            raise ValueError("投影电路与独立范数迹公式不一致："+name)


def run(args):
    from flint import nmod_poly
    source = HERE / "system-runs" / args.system
    line = SHARED / "quotient-runs" / args.line
    destination = HERE / "norm-runs" / args.tag
    if destination.exists():
        raise FileExistsError("拒绝覆盖已有范数回放")
    destination.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, norm = read(line / "quotient.json.gz"), read(line / "norm.json.gz")
        accepted = read(line / "verified-final.json.gz")
        if (not accepted["verified"] or accepted["quotient_sha256"] != sha(line / "quotient.json.gz")
                or accepted["norm_sha256"] != sha(line / "norm.json.gz")):
            raise ValueError("旧规范形没有绑定接受结果")
        circuit = read(source / "circuit.json.gz")
        if circuit["keep_left_roots"] != 0:
            raise ValueError("只接受四参数投影电路")
        setup(data["prime"], args.seconds)
        parameters = [Rational(value) for value in data["fixed_u_v_w"]]+[Rational(nmod_poly([0, 1], data["prime"]))]
        outputs = circuit["outputs"]["critical_equations"]+[circuit["outputs"]["radius_derivative"]]
        small, indices = prune(circuit, outputs)
        registers = execute(small, parameters)
        actual = [registers[index] for index in indices]
        print({"stage": "projected_circuit", "seconds": time.perf_counter()-tick, "pruned_nodes": len(small["nodes"])}, flush=True)
        algebra = load(data)
        H = vector(data["H"], algebra)
        Hgrad = [vector(row, algebra) for row in data["gradient_u_v_w_r"]]
        A = right_leading(algebra, parameters)
        NA = algebra.norm(A.value.promote(4).coefficients, 4)
        NH = algebra.norm(H.coefficients, 5)
        if not same(NH, decode_scalar(norm["value"], data["prime"])):
            raise ValueError("旧闭环范数未能从32维规范形重算")
        G = NA*NH
        if not same(actual[0], G):
            raise ValueError("四参数闭环不等于Norm(A)*Norm(H)")
        print({"stage": "independent_norms", "seconds": time.perf_counter()-tick}, flush=True)
        Hinv, Ainv = H.inverse(), A.value.inverse()
        expected = [G]
        for j in range(4):
            relative_H, relative_A = Hinv*Hgrad[j], Ainv*A.gradient[j]
            logarithmic = trace(algebra, relative_H.coefficients, 5)+trace(algebra, relative_A.promote(4).coefficients, 4)
            expected.append(G*logarithmic)
            print({"stage": "independent_trace_"+str(j), "seconds": time.perf_counter()-tick}, flush=True)
        compare(actual, expected)
        direct_r = Rational(G.n.derivative()*G.d-G.n*G.d.derivative(), G.d*G.d)
        if not same(actual[4], direct_r):
            raise ValueError("半径导数与单变量有理函数的直接微分不一致")
        write(destination / "functions.json.gz", {"prime": data["prime"], "fixed_u_v_w": data["fixed_u_v_w"],
              "G_and_gradient_u_v_w_r": [value.encode() for value in actual], "norm_A": NA.encode(), "norm_H": NH.encode()})
        report = {"verified": True, "symbolic_total_derivatives": 4, "identity": "G=Norm_16(A)*Norm_32(H)",
                  "derivative_identity": "G_theta=G*[Trace_16(A^-1*A_theta)+Trace_32(H^-1*H_theta)]",
                  "G_numerator_degree": int(G.n.degree()), "G_denominator_degree": int(G.d.degree()),
                  "norm_A_degrees": [int(NA.n.degree()), int(NA.d.degree())],
                  "prime": data["prime"], "fixed_u_v_w": data["fixed_u_v_w"],
                  "source_quotient_sha256": sha(line / "quotient.json.gz"), "source_norm_sha256": sha(line / "norm.json.gz"),
                  "circuit_sha256": sha(source / "circuit.json.gz"), "driver_sha256": sha(__file__),
                  "functions_sha256": sha(destination / "functions.json.gz"), "seconds": time.perf_counter()-tick,
                  "limits": snapshot(), "scope": "精确Fp(r)截面及四个全参数导数；完整参数等价性另由通用恒等式和图表单位证明，不以有限截面推出特征零恒等式。"}
        write(destination / "verified.json.gz", report)
        print({"verified": True, "degrees": [report["G_numerator_degree"], report["G_denominator_degree"]], "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", default="20261015-projected4-final")
    parser.add_argument("--line", default="20261015-jets-p65521")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--seconds", type=int, default=120)
    run(parser.parse_args())
