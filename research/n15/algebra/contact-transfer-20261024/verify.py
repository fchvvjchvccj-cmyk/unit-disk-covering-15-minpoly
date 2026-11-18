"""重验原根、几何传递边界、本原恒等式及新八变量系统的压缩证书。"""

from fractions import Fraction as Q
import argparse
import hashlib
import time

from context import ROOT, INPUTS, RESULTS, read, write, payload, sha, bounded, execute
from intervals import Dyadic, derivatives, evaluate_rows
from model import build


class Interval:
    arithmetic = None

    def __init__(self, value):
        self.value = value

    def cast(self, other):
        return other if isinstance(other, Interval) else Interval(self.arithmetic.rational(other))

    def __add__(self, other):
        return Interval(self.arithmetic.add(self.value, self.cast(other).value))

    __radd__ = __add__

    def __neg__(self):
        return Interval(self.arithmetic.neg(self.value))

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        return Interval(self.arithmetic.mul(self.value, self.cast(other).value))

    __rmul__ = __mul__

    def __truediv__(self, other):
        return Interval(self.arithmetic.div(self.value, self.cast(other).value))

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def __pow__(self, exponent):
        if exponent < 0:
            return 1/self**(-exponent)
        result = self.cast(1)
        for _ in range(exponent):
            result *= self
        return result


class Dual:
    """对全部独立端口求导，与生产端的四个隐函数方向分开实现。"""

    def __init__(self, value, gradient):
        self.value, self.gradient = value, tuple(gradient)

    def cast(self, other):
        return other if isinstance(other, Dual) else Dual(self.value*0+other, [a*0 for a in self.gradient])

    def __add__(self, other):
        other = self.cast(other)
        return Dual(self.value+other.value, [a+b for a, b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-a for a in self.gradient])

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        return Dual(self.value*other.value, [a*other.value+self.value*b for a, b in zip(self.gradient, other.gradient)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        value = self.value/other.value
        return Dual(value, [(a-value*b)/other.value for a, b in zip(self.gradient, other.gradient)])

    def __rtruediv__(self, other):
        return self.cast(other)/self


def require(condition, message):
    if not condition:
        raise ValueError(message)


def encode(value):
    return [str(Q(a, Interval.arithmetic.scale)) for a in value.value]


def away(value, label):
    lo, hi = value.value
    require(lo > 0 or hi < 0, "非零条件失败："+label)
    return encode(value)


def box_values(center, radius, bits=512):
    Interval.arithmetic = Dyadic(bits)
    A = Interval.arithmetic
    return [Interval((A.rational(x-radius)[0], A.rational(x+radius)[1])) for x in center]


def independent_variables(values):
    n = len(values)
    return [Dual(x, [x*0+int(i == j) for j in range(n)]) for i, x in enumerate(values)]


def contraction(values, jacobian, inverse, radius, bits):
    A, n = Interval.arithmetic, len(values)
    inverse = [[Interval(A.rational(x)) for x in row] for row in inverse]
    rows, residual = [], []
    for i in range(n):
        error = 0
        for j in range(n):
            product = sum((inverse[i][k]*jacobian[k][j] for k in range(n)), Interval(A.zero))
            error += A.absolute((int(i == j)-product).value)
        defect = sum((inverse[i][k]*values[k] for k in range(n)), Interval(A.zero))
        rows.append(Q(error, A.scale))
        residual.append(Q(A.absolute(defect.value), A.scale))
    q, inclusion = max(rows), max(e/radius+c for e, c in zip(residual, rows))
    require(q < 1 and inclusion < 1, "严格自映与压缩失败")
    return {"contraction_upper": str(q), "inclusion_ratio_upper": str(inclusion), "bits": bits,
            "radius": str(radius), "unique_root_in_box": True}


def verify_sources():
    manifest = read(INPUTS / "manifest.json")
    for item in manifest["files"]:
        path = INPUTS / item["path"]
        require(sha(path) == item["stored_sha256"], "输入存储摘要错误")
        require(hashlib.sha256(payload(path)).hexdigest() == item["payload_sha256"], "输入内容摘要错误")
    data, root = read(INPUTS / "axis25-system.json.gz"), read(INPUTS / "axis25-root.json.gz")
    require(root["system_sha256"] == hashlib.sha256(payload(INPUTS / "axis25-system.json.gz")).hexdigest(), "原根未绑定系统")
    center = [Q(int(v), int(root["Qx"])) for v in root["xnum"]]
    radius = Q(root["rho"])
    n = len(center)
    require(n == 25 and root["variables"] == data["variables"], "原系统布局错误")
    terms = [[(t["coefficient"], tuple(map(tuple, t["powers"]))) for t in row] for row in data["equations"]]
    parsed = [dict((powers, coefficient) for coefficient, powers in row) for row in terms]
    for j in range(14):
        expected = {(): 1} if j == 13 else {}
        for k, polynomial in enumerate(terms[:11]):
            for coefficient, powers in polynomial:
                exponent = dict(powers)
                degree = exponent.get(j, 0)
                if degree:
                    exponent[j] -= 1
                    exponent[14+k] = 1
                    key = tuple(sorted((i, e) for i, e in exponent.items() if e))
                    expected[key] = expected.get(key, 0)+coefficient*degree
        require(parsed[11+j] == {k: v for k, v in expected.items() if v}, "原驻定行不等于实际拉格朗日导数")
    bits = 512
    box_values(center, radius, bits)
    A = Interval.arithmetic
    exact = [A.rational(x) for x in center]
    box = [(A.rational(x-radius)[0], A.rational(x+radius)[1]) for x in center]
    values = [Interval(x) for x in evaluate_rows(A, terms, exact)]
    jacobian = [[Interval(x) for x in evaluate_rows(A, row, box)] for row in derivatives(terms, n)]
    inverse = [[Q(int(v), int(root["Qy"])) for v in row] for row in root["Ynum"]]
    result = contraction(values, jacobian, inverse, radius, bits)
    from prove_branch import source_identities
    result["integer_geometry_rows_checked"] = source_identities(data)
    result["kkt_gradient_rows_checked"] = 14
    result["scope"] = "独立认证恢复的25元轴系统目标根；本轮未重放原148元来源链及全局证明。"
    return data, root, center, result


def universal_identities():
    import sympy as S
    from prove_shared import identity
    require(identity()["verified"], "共享根恒等式失败")
    ax, ay, cx, cy, bx, by, qx, qy, k = S.symbols("ax ay cx cy bx by qx qy k")
    norm = cx*cx+cy*cy
    p = [ax+cx-qx, ay+cy-qy]
    qp = [(norm+k)*cx/norm-qx, (norm+k)*cy/norm-qy]
    pp = [cx+bx-p[0], cy+by-p[1]]
    transfer = [bx-ax+k*cx/norm, by-ay+k*cy/norm]
    require(all(S.cancel(qp[j]+pp[j]-[cx, cy][j]-transfer[j]) == 0 for j in (0, 1)), "二阶传递恒等式失败")
    aa, ab = ax*ax+ay*ay, ax*bx+ay*by
    ba, qa, qb = bx*bx+by*by, ax*qx+ay*qy, bx*qx+by*qy
    ka, kb, area = aa+k, ba+k, ax*by-ay*bx
    H = (ka*by-kb*ay)**2+(ax*kb-bx*ka)**2-4*area*area
    f = ka-2*qa
    P = ba*(ka+2*qa)-4*qb*ab
    R = kb-2*qb
    quotient = aa*(kb+2*qb)-2*ka*ab
    require(S.expand(H-f*P-R*quotient-4*(qx*qx+qy*qy-1)*area*area) == 0, "本原闭环几何因式分解失败")
    a, b, c, j0, j1, d = S.symbols("A B C J0 J1 d")
    g, F = a*d*d+b*d+c, a*j0*j0-b*j0*j1+c*j1*j1
    require(S.cancel(F-a*(j0+j1*d)*(j0+j1*(-b/a-d))-j1*j1*g) == 0, "相对二次范数恒等式失败")
    # 此微分恒等式在共同零集上验证，避免只比较浮点导数。
    da, db, dc, dj0, dj1 = S.symbols("dA dB dC dJ0 dJ1")
    dd = -(da*d*d+db*d+dc)/(2*a*d+b)
    dF = da*j0*j0+2*a*j0*dj0-db*j0*j1-b*dj0*j1-b*j0*dj1+dc*j1*j1+2*c*j1*dj1
    require(S.cancel((dF+j1*(2*a*d+b)*(dj0+dj1*d+j1*dd)).subs({j0: -j1*d, c: -a*d*d-b*d})) == 0, "相对临界微分恒等式失败")
    from identities import verify as endpoint_identities
    return {"endpoint_factorizations": endpoint_identities(),
            "shared_discriminant_identity": True, "second_order_transfer_identity": True,
            "primitive_boundary_factor_identity": True, "relative_norm_identity": True,
            "relative_critical_differential_identity": True}


def target_map(values):
    from geometry import circumcenter
    a, U, b = [values[0], -values[1]], [values[10], values[10]*0], values[11:13]
    px, py = circumcenter(a, U, b)
    u, v, w, r = values[10:14]
    t = 2*((v-u)*py-w*(px-u))
    return [t, 2*a[0]-px, a[1], values[9], u, v, w, r]


def transfer_audit(center, radius):
    from geometry import dot, subtract, reconstruct, previous_boundary, common_boundary, circumcenter
    values = box_values(center, radius)
    r, b, a = values[13], values[11:13], [values[0], -values[1]]
    c3, c4, c5 = values[2:4], values[4:6], values[6:8]
    u, v, w = values[10:13]
    px, py = circumcenter(a, [u, u*0], b)
    C1 = px*u-px*a[0]-py*a[1]-u*u+u*a[0]
    C2 = -px*px+4*px*u-4*px*a[0]-py*py-2*u*u+2*u*a[0]+1
    Z = 2*values[8]-1-r
    d0, e0 = c5
    minus = w*(d0*d0+e0*e0-(Z+u)*d0+Z*u)-(w*w+(v-Z)*(v-u))*e0
    aa = 2*(px-u)
    opposite = (2*values[8]-1-r-v)+(v-u-aa)
    units = {"left_C1": away(C1, "左端C1"), "left_C2": away(C2, "左端C2"),
             "right_minus": away(minus, "右端补圆因子"), "opposite_shared_root_factor": away(opposite, "共享根另一因子")}
    q0, p0 = previous_boundary(a, u, r), [px, py]
    q1, p1 = common_boundary(a, c3, r), circumcenter(a, c3, b)
    q2, p2 = common_boundary(c3, c4, r), circumcenter(c3, c4, b)
    for label, pair in (("q1_ne_q0", subtract(q1, q0)), ("p1_ne_p0", subtract(p1, p0)),
                        ("q2_ne_q1", subtract(q2, q1)), ("p2_ne_p1", subtract(p2, p1)),
                        ("a_ne_c3", subtract(a, c3)), ("c3_ne_c4", subtract(c3, c4)),
                        ("q1_ne_p1", subtract(q1, p1)), ("q2_ne_p2", subtract(q2, p2))):
        squared = dot(pair, pair)
        require(squared.value[0] > 0, "互补交点未严格分离："+label)
        units[label] = encode(squared)
    predicted = reconstruct(values[:2], values[10], b, r)
    for name, actual in (("c3", c3), ("c4", c4)):
        require(all((x-y).value[0] <= 0 <= (x-y).value[1] for x, y in zip(predicted[name], actual)), "有理传递与原根不相容")
    d, e, f = subtract(c5, c4), subtract(b, c4), subtract(b, c5)
    cross = d[0]*e[1]-d[1]*e[0]
    inactive = dot(d, d)*dot(e, e)*dot(f, f)-4*r*r*cross*cross
    require(inactive.value[1] < 0, "未能认证非活跃面有严格余量")
    wrong = [b[i]-c3[i]+(1-r*r)*c4[i]/dot(c4, c4) for i in (0, 1)]
    error = subtract(wrong, c5)
    distance = dot(error, error)
    require(distance.value[0] > 0, "错误延长传递未被排除")
    from importlib.util import spec_from_file_location, module_from_spec
    path = ROOT / "references/diskcoveringsolve-089584973f41e6374ce110f6e5b5ee456df19352/cover15/proof_bundle/core/system15.py"
    spec = spec_from_file_location("transfer_source_graph", path)
    graph = module_from_spec(spec)
    spec.loader.exec_module(graph)
    require((2, 3, 12) in graph.ACTIVE and (3, 4, 12) in graph.ACTIVE and (4, 5, 12) not in graph.ACTIVE, "活跃面表不符")
    return {"active_fan": [[2, 3, 12], [3, 4, 12]], "inactive_next_face": [4, 5, 12],
            "endpoint_and_orientation_units": units,
            "inactive_triangle_polynomial_interval": encode(inactive),
            "false_extension_squared_error_interval": encode(distance), "source_graph_sha256": sha(path),
            "conclusion": "固定内圆心时是二阶有理传递；前一圆心随步变化，且下一面非活跃，不能将整条外链当作固定系数Möbius复合。"}


def verify_new(data, source_center, source_radius, certificate=None):
    import mpmath as mp
    require(data == build(), "新电路与当前构造器不同")
    n = 8
    center = target_map(source_center)
    indices = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    radius = Q(1, 10**50)
    if certificate is None:
        mp.mp.dps = 100
        numeric = [mp.mpf(x.numerator)/x.denominator for x in center]
        registers = execute(data, independent_variables(numeric))
        matrix = mp.matrix([registers[i].gradient for i in indices])
        inv = matrix**-1
        scale = 10**90
        inverse = [[Q(int(mp.nint(inv[i, j]*scale)), scale) for j in range(n)] for i in range(n)]
    else:
        require(certificate["center"] == [str(x) for x in center] and Q(certificate["radius"]) == radius,
                "新证书未绑定实际原根像或既定根盒")
        require(certificate["circuit_sha256"] == sha(RESULTS / "primitive16.json.gz") and certificate["bits"] == 512,
                "新证书的电路或区间精度不匹配")
        inverse = [[Q(x) for x in row] for row in certificate["inverse"]]
        require(len(inverse) == n and all(len(row) == n for row in inverse), "预条件器维数错误")
    box = box_values(center, radius)
    registers_box = execute(data, independent_variables(box))
    jacobian = [registers_box[i].gradient for i in indices]
    exact = box_values(center, Q(0))
    registers_exact = execute(data, exact)
    values = [registers_exact[i] for i in indices]
    result = contraction(values, jacobian, inverse, radius, 512)
    box_values(center, radius)
    guards = data["all_divisor_nodes"]+data["simple_root_derivative_nodes"]+data["closed_chart_nonzero_guards"]+data["critical_locus_nonzero_guards"]
    guard_values = {str(i): away(registers_box[i].value, "新系统门"+str(i)) for i in sorted(set(guards))}
    original_box = box_values(source_center, source_radius)
    mapped = target_map(original_box)
    distances = [max(abs(Q(a, Interval.arithmetic.scale)-c) for a in x.value) for c, x in zip(center, mapped)]
    require(max(distances) < radius, "原根像未包含于新根盒")
    transformed = execute(data, mapped)
    source_geometry = {"c0": original_box[:2], "c3": original_box[2:4], "c4": original_box[4:6],
                       "c5": original_box[6:8], "c6": original_box[8:10]}
    for name, nodes in data["outputs"]["geometry"].items():
        for index, expected in zip(nodes, source_geometry[name]):
            difference = transformed[index]-expected
            require(difference.value[0] <= 0 <= difference.value[1], "回代与原根不相容："+name)
    result.update(variables=data["variables"], equations=8, carrier_rank=16,
                  source_image_max_distance_upper=str(max(distances)), source_image_in_new_box=True,
                  nonzero_guards=guard_values, recovered_coordinates_consistent=10,
                  target_connection="通用传递/共享根/本原因式分解及局部可逆回代提供精确连接；区间残差只作相容性检查。",
                  critical_multiplier_interval=encode(transformed[data["outputs"]["critical_multiplier"]]),
                  removed_factor_interval=encode(transformed[data["outputs"]["removed_factor"]]),
                  radius_derivative_interval=encode(transformed[data["outputs"]["radius_derivative"]]))
    certificate = {"center": [str(x) for x in center], "radius": str(radius),
                   "inverse": [[str(x) for x in row] for row in inverse], "bits": 512,
                   "circuit_sha256": sha(RESULTS / "primitive16.json.gz")}
    return result, certificate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default="initial")
    args = parser.parse_args()
    tick = time.perf_counter()
    with bounded() as limits:
        _, root, center, source = verify_sources()
        print({"stage": "source_root", "verified": True}, flush=True)
        identities = universal_identities()
        audit = transfer_audit(center, Q(root["rho"]))
        print({"stage": "transfer_boundary", "verified": True}, flush=True)
        candidate, certificate = verify_new(read(RESULTS / "primitive16.json.gz"), center, Q(root["rho"]))
        certificate_path = RESULTS / ("target-"+args.tag+"-certificate.json.gz")
        if certificate_path.exists():
            require(read(certificate_path) == certificate, "已有目标证书与本轮生成数据不同")
        else:
            write(certificate_path, certificate)
        report = {"verified": True, "source_root": source, "universal_identities": identities,
                  "transfer_audit": audit, "candidate": candidate,
                  "certificate_sha256": sha(certificate_path),
                  "verifier_sha256": sha(__file__), "limits": limits(), "seconds": time.perf_counter()-tick,
                  "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False,
                  "global_covering_proof_replayed": False}
        write(RESULTS / ("verified-"+args.tag+".json"), report)
        print({"verified": True, "seconds": report["seconds"], "equations": 8}, flush=True)


if __name__ == "__main__":
    main()
