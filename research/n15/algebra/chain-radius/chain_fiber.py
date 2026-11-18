"""在共同有限扩域中双向重建接触链，直接计算闭环范数及四个偏导。"""

import argparse
import json
import time

from common import HERE, AXIS, bounded, read, write, sha

FACTOR = HERE.parent / "three-parameter-elimination/c0-critical/factor-branch-verified.json"


class Jet:
    """值和四个参数方向的一阶微分；所有算术都在同一个精确有限域内。"""

    sqrt_cache = {}
    sqrt_calls = 0
    sqrt_hits = 0

    def __init__(self, value, gradient=None):
        self.value = value
        self.gradient = tuple(gradient) if gradient is not None else (value*0,)*4

    def cast(self, other):
        return other if isinstance(other, Jet) else Jet(self.value*0+other)

    def __add__(self, other):
        other = self.cast(other)
        return Jet(self.value+other.value, [a+b for a, b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, [-a for a in self.gradient])

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        return Jet(self.value*other.value, [a*other.value+self.value*b for a, b in zip(self.gradient, other.gradient)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        if other.value.is_zero():
            raise ZeroDivisionError("当前纤维落入图表分母零点")
        inv = other.value.inverse()
        value = self.value*inv
        return Jet(value, [(a-value*b)*inv for a, b in zip(self.gradient, other.gradient)])

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def __pow__(self, exponent):
        if exponent < 0:
            return 1/(self**(-exponent))
        if exponent == 0:
            return self.cast(1)
        value = self.value**exponent
        factor = exponent*self.value**(exponent-1)
        return Jet(value, [factor*a for a in self.gradient])

    def sqrt_pair(self):
        if self.value.is_zero():
            raise ValueError("纤维出现二次重根")
        key = tuple(map(int, self.value.to_list()))
        if key in self.sqrt_cache:
            value = self.sqrt_cache[key]
            Jet.sqrt_hits += 1
        else:
            if not self.value.is_square():
                raise ValueError("共同扩域未分裂该二次关系")
            value = self.value.sqrt()
            Jet.sqrt_calls += 1
            # 同一Frobenius轨道的根可直接复用；两种符号仍全部保留。
            rad, root = self.value, value
            while True:
                self.sqrt_cache[tuple(map(int, rad.to_list()))] = root
                rad, root = rad.frobenius(), root.frobenius()
                if rad == self.value:
                    break
        result = Jet(value, [a/(2*value) for a in self.gradient])
        return result, -result

    def encode(self):
        return {"value": list(map(int, self.value.to_list())),
                "gradient": [list(map(int, a.to_list())) for a in self.gradient]}


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]


def add(a, b):
    return [a[i]+b[i] for i in range(2)]


def sub(a, b):
    return [a[i]-b[i] for i in range(2)]


def rotate(a):
    return [-a[1], a[0]]


def scale(a, k):
    return [x*k for x in a]


def boundary(a, b, r):
    ka, kb = 1+dot(a, a)-r*r, 1+dot(b, b)-r*r
    nx, ny = ka*b[1]-kb*a[1], a[0]*kb-b[0]*ka
    area = a[0]*b[1]-a[1]*b[0]
    return nx*nx+ny*ny-4*area*area


def circle_step(a, b, r):
    """两个接触根各保留正负，返回完整的四个非重合分支。"""
    T = dot(a, a)
    K = 1+T-r*r
    d = sub(b, a)
    L = dot(d, d)
    delta, eta = 4*T-K*K, L*(4*r*r-L)
    output = []
    for s in delta.sqrt_pair():
        q = scale(add(scale(a, K), scale(rotate(a), s)), 1/(2*T))
        for t in eta.sqrt_pair():
            p = add(scale(add(a, b), r.cast(1)/2), scale(rotate(d), t/(2*L)))
            output.append({"point": sub(add(q, p), a), "s": s, "t": t})
    return output


def polynomial(coefficients, x):
    value = x.cast(0)
    for c in reversed(coefficients):
        value = value*x+c
    return value


def left_quartic(params, t, rows):
    u, v, w, r = params
    dx, L = v-u, (v-u)**2+w*w
    a, b = dx-w*t/L, w+dx*t/L
    s = b/(2*r+a)
    coefficients = [r.cast(0) for _ in range(5)]
    for row in rows:
        powers = dict(row["powers"])
        e = powers.pop(0, 0)
        term = r.cast(row["coefficient"])
        for j, d in powers.items():
            term *= {1: s, 2: u, 3: r}[j]**d
        coefficients[e] += term
    return a, b, coefficients


def left_endpoint(params, factor, field):
    from flint import fq_default_poly_ctx
    u, v, w, r = params
    L = (v-u)**2+w*w
    eta = L*(4*r*r-L)
    rows = factor["factor_polynomials"][factor["target_selected_factor_index"]]
    output = []
    for t in eta.sqrt_pair():
        a, b, coefficients = left_quartic(params, t, rows)
        fp = fq_default_poly_ctx(field)([c.value for c in coefficients])
        roots = fp.roots()
        if len(roots) != 4 or any(m != 1 for k, m in roots):
            raise ValueError("四次载体在共同扩域中没有四个单根")
        for value, multiplicity in roots:
            partial = [value*0]*4
            for c in reversed(coefficients):
                partial = [a*value+b for a, b in zip(partial, c.gradient)]
            derivative = fp.derivative()(value)
            if derivative.is_zero():
                raise ValueError("四次载体出现重根")
            k = Jet(value, [-a/derivative for a in partial])
            X = (a-b*k)/(1+k*k)
            output.append({"point": [u+X, k*X], "t": t, "k": k})
    return output


def right_endpoint(params):
    u, v, w, r = params
    E = (u-v)**2+w*w
    beta = r*r-1
    output = []
    for root in (4*w*w*r*r/E-w*w).sqrt_pair():
        x = (1+r+v+root)/2
        for y in (r*r-(x-1)**2).sqrt_pair():
            S = 2*x+beta
            ac, ec = 2-4*x*x/S, -4*x*y/S
            Z = 2*x-1-r
            Pz = (v-Z)*(v-u)
            bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
            ld, le, l0 = bc-ac, be-ec, bc0+beta
            e0, e1 = -l0/le, -ld/le
            A = 1+e1*e1
            B = 2*e0*e1+ac+ec*e1
            C = e0*e0+ec*e0-beta
            for disc in (B*B-4*A*C).sqrt_pair():
                d = (-B+disc)/(2*A)
                e = e0+e1*d
                output.append({"c6": [x, y], "point": [d, e]})
    return output


def prime_scalar(a):
    coefficients = list(map(int, a.to_list()))
    if any(coefficients[1:]):
        raise ValueError("完整分支乘积没有下降到素域")
    return coefficients[0] if coefficients else 0


def fiber(field, values, factor, save_points=True):
    Jet.sqrt_cache = {}
    Jet.sqrt_calls = Jet.sqrt_hits = 0
    params = [Jet(field(v), [field(int(i == j)) for j in range(4)]) for i, v in enumerate(values)]
    u, v, w, r = params
    c12 = [v, w]
    tick = time.perf_counter()
    left0 = left_endpoint(params, factor, field)
    left3, left4 = [], []
    for i, root in enumerate(left0):
        a = [root["point"][0], -root["point"][1]]
        for root3 in circle_step(a, c12, r):
            root3["parent"] = i
            j = len(left3)
            left3.append(root3)
            for root4 in circle_step(root3["point"], c12, r):
                root4["parent"] = j
                left4.append(root4)
    right = right_endpoint(params)
    if (len(left0), len(left3), len(left4), len(right)) != (8, 32, 128, 8):
        raise ValueError("链载体分支数不符合单首项扩张的维数")
    reconstructed = time.perf_counter()
    product = r.cast(1)
    values_at_bridge = []
    zeros = []
    for i, root4 in enumerate(left4):
        for j, root5 in enumerate(right):
            contact = boundary(root4["point"], root5["point"], r)
            product *= contact
            if contact.value.is_zero():
                zeros.append([i, j])
            if save_points:
                values_at_bridge.append(contact.encode())
    result = {"parameters_u_v_w_r": values,
              "dimensions_at_fixed_parameters": {"c0": 8, "c0_c3": 32, "c0_c3_c4": 128, "c5_c6_B2": 8, "before_closure": 1024},
              "norm_h3_mod_p": prime_scalar(product.value),
              "norm_gradient_u_v_w_r_mod_p": [prime_scalar(a) for a in product.gradient],
              "vanishing_bridge_pairs": zeros,
              "native_square_roots": Jet.sqrt_calls, "square_root_cache_hits": Jet.sqrt_hits,
              "reconstruction_seconds": reconstructed-tick, "bridge_product_seconds": time.perf_counter()-reconstructed}
    if save_points:
        def encode_root(root):
            return {k: [a.encode() for a in value] if isinstance(value, list)
                    else value.encode() if isinstance(value, Jet) else value for k, value in root.items()}
        result.update(left0=[encode_root(a) for a in left0], left3=[encode_root(a) for a in left3],
                      left4=[encode_root(a) for a in left4], right=[encode_root(a) for a in right],
                      bridge_values=values_at_bridge)
    return result


def run(args):
    from flint import fq_default_ctx
    out = HERE / "fiber-runs" / args.tag
    if out.exists():
        raise FileExistsError("拒绝覆盖已有实验")
    out.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        factor = read(FACTOR)
        if not factor["verified"] or factor["target_selected_degree"] != 4:
            raise ValueError("左端目标载体证书发生变化")
        field = fq_default_ctx(args.prime, args.degree, "a", fq_type="FQ_NMOD")
        write(out / "field.json", {"prime": args.prime, "degree": args.degree,
                                   "modulus_low_to_high": list(map(int, field.modulus().coeffs())),
                                   "source_system_sha256": sha(AXIS), "left_factor_sha256": sha(FACTOR),
                                   "generator_sha256": sha(__file__),
                                   "scope": "参数专门化的精确有限域纤维与偏导，不是半径消元多项式"})
        records, failures = [], []
        for i, values in enumerate(([2, 3, 1, 2], [3, 5, 2, 3], [2, 5, 3, 1], [7, 11, 3, 2])):
            try:
                data = fiber(field, values, factor)
            except (ValueError, ZeroDivisionError) as error:
                failure = {"parameters": values, "error": str(error)}
                failures.append(failure)
                write(out / f"rejected-fiber-{i}.json", failure)
                print(json.dumps(failure, ensure_ascii=False), flush=True)
                continue
            write(out / f"fiber-{i}.json.gz", data)
            record = {k: v for k, v in data.items() if k not in ("left0", "left3", "left4", "right", "bridge_values")}
            records.append(record)
            print(json.dumps(record), flush=True)
            if len(records) == 2:
                break
        write(out / "summary.json", {"fibers": records, "rejected_fibers": failures, "seconds": time.perf_counter()-tick,
                                      "limits": snapshot(), "minimal_polynomial_obtained": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--degree", type=int, default=96)
    run(parser.parse_args())
