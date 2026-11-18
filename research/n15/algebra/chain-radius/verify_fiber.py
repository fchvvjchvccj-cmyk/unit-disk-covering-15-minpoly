"""独立重放原接触方程、切向微分和完整1024项闭环乘积。"""

import argparse
import json
import time

from common import HERE, AXIS, bounded, read, write, sha

FACTOR = HERE.parent / "three-parameter-elimination/c0-critical/factor-branch-verified.json"


def decode(data, field):
    def element(coefficients):
        return field(coefficients)
    return element(data["value"]), tuple(element(v) for v in data["gradient"])


def add(a, b):
    return a[0]+b[0], tuple(x+y for x, y in zip(a[1], b[1]))


def mul(a, b):
    return a[0]*b[0], tuple(x*b[0]+a[0]*y for x, y in zip(a[1], b[1]))


def times(a, n):
    return a[0]*n, tuple(x*n for x in a[1])


def power(a, n):
    one = a[0]*0+1
    if not n:
        return one, (one*0,)*4
    return a[0]**n, tuple(n*a[0]**(n-1)*x for x in a[1])


def subtract(a, b):
    return add(a, times(b, -1))


def divide(a, b):
    if b[0].is_zero():
        raise ValueError("重放公式出现零分母")
    v = a[0]/b[0]
    return v, tuple((x-v*y)/b[0] for x, y in zip(a[1], b[1]))


def key(value):
    return tuple(map(int, value[0].to_list()))


def point_key(point):
    return tuple(key(a) for a in point)


def row_value(row, values, field):
    result = field(0), (field(0),)*4
    for term in row:
        value = field(term["coefficient"]), (field(0),)*4
        for j, e in term["powers"]:
            value = mul(value, power(values[j], e))
        result = add(result, value)
    return result


def nonzero_root_assertion(value):
    if not value[0].is_zero() or any(not a.is_zero() for a in value[1]):
        raise ValueError("原方程或沿参数方向的微分不为零")


def close_boundary(a, b, radius, field):
    one = field(1), (field(0),)*4
    ka = subtract(add(one, add(power(a[0], 2), power(a[1], 2))), power(radius, 2))
    kb = subtract(add(one, add(power(b[0], 2), power(b[1], 2))), power(radius, 2))
    nx = subtract(mul(ka, b[1]), mul(kb, a[1]))
    ny = subtract(mul(a[0], kb), mul(b[0], ka))
    det = subtract(mul(a[0], b[1]), mul(a[1], b[0]))
    return subtract(add(power(nx, 2), power(ny, 2)), times(power(det, 2), 4))


def check_boundary_identity(system):
    import sympy as S
    symbols = S.symbols("x:14")
    ax, ay, bx, by, r = (symbols[i] for i in (4, 5, 6, 7, 13))
    ka, kb = 1+ax*ax+ay*ay-r*r, 1+bx*bx+by*by-r*r
    B = (ka*by-kb*ay)**2+(ax*kb-bx*ka)**2-4*(ax*by-ay*bx)**2
    source = sum(t["coefficient"]*S.prod(symbols[j]**e for j, e in t["powers"]) for t in system["equations"][3])
    if S.Poly(B-source, *symbols):
        raise ValueError("快速闭环公式与原h3整数系数不同")


def check_left_completeness(data, decoded, base, field, factor):
    zero = field(0), (field(0),)*4
    one = field(1), (field(0),)*4
    u, v, w, r = [base[i] for i in (10, 11, 12, 13)]
    dx = subtract(v, u)
    L = add(power(dx, 2), power(w, 2))
    eta = mul(L, subtract(times(power(r, 2), 4), L))
    groups = {}
    rows = factor["factor_polynomials"][factor["target_selected_factor_index"]]
    for entry, point in zip(data["left0"], decoded["left0"]):
        t, k = [decode(entry[name], field) for name in ("t", "k")]
        nonzero_root_assertion(subtract(power(t, 2), eta))
        a = subtract(dx, divide(mul(w, t), L))
        b = add(w, divide(mul(dx, t), L))
        s = divide(b, add(times(r, 2), a))
        leading = [dict(row, powers=[[j, e] for j, e in row["powers"] if j != 0])
                   for row in rows if dict(row["powers"]).get(0, 0) == 4]
        if row_value(leading, {1: s, 2: u, 3: r}, field)[0].is_zero():
            raise ValueError("左端四次首项为零，不能据四个根判断完备性")
        nonzero_root_assertion(row_value(rows, {0: k, 1: s, 2: u, 3: r}, field))
        X = divide(subtract(a, mul(b, k)), add(one, power(k, 2)))
        for value in (subtract(point[0], add(u, X)), subtract(point[1], mul(k, X))):
            nonzero_root_assertion(value)
        groups.setdefault(key(t), set()).add(key(k))
    if len(groups) != 2 or any(len(roots) != 4 for roots in groups.values()):
        raise ValueError("左端没有覆盖两个二次方向上的全部四次根")
    if len({point_key(p) for p in decoded["left0"]}) != 8:
        raise ValueError("左端存在重复根")
    for name, parent_name in (("left3", "left0"), ("left4", "left3")):
        groups = {}
        for entry, point in zip(data[name], decoded[name]):
            a = decoded[parent_name][entry["parent"]]
            if name == "left3":
                a = [a[0], times(a[1], -1)]
            b = [v, w]
            T = add(power(a[0], 2), power(a[1], 2))
            K = subtract(add(one, T), power(r, 2))
            d = [subtract(b[i], a[i]) for i in range(2)]
            LL = add(power(d[0], 2), power(d[1], 2))
            s, t = [decode(entry[n], field) for n in ("s", "t")]
            nonzero_root_assertion(subtract(power(s, 2), subtract(times(T, 4), power(K, 2))))
            nonzero_root_assertion(subtract(power(t, 2), mul(LL, subtract(times(power(r, 2), 4), LL))))
            if any(q[0].is_zero() for q in (T, LL, s, t)):
                raise ValueError("几何四分支的非退化条件失败")
            ja, jd = [times(a[1], -1), a[0]], [times(d[1], -1), d[0]]
            for i in range(2):
                qn = add(mul(K, a[i]), mul(s, ja[i]))
                pn = add(mul(LL, add(a[i], b[i])), mul(t, jd[i]))
                numerator = subtract(add(mul(LL, qn), mul(T, pn)), times(mul(mul(T, LL), a[i]), 2))
                nonzero_root_assertion(subtract(numerator, times(mul(mul(T, LL), point[i]), 2)))
            groups.setdefault(entry["parent"], set()).add((key(s), key(t)))
        if set(groups) != set(range(len(decoded[parent_name]))) or any(len(a) != 4 for a in groups.values()):
            raise ValueError("接触链四分支未完整覆盖或出现重复")


def check_right_completeness(data, decoded, base, field):
    one = field(1), (field(0),)*4
    u, v, w, r = [base[i] for i in (10, 11, 12, 13)]
    beta = subtract(power(r, 2), one)
    groups, c6groups = {}, {}
    for entry, point in zip(data["right"], decoded["right"]):
        x, y = [decode(a, field) for a in entry["c6"]]
        d, e = point
        E = add(power(subtract(u, v), 2), power(w, 2))
        root = subtract(subtract(subtract(times(x, 2), one), r), v)
        nonzero_root_assertion(subtract(power(root, 2), subtract(divide(times(mul(power(w, 2), power(r, 2)), 4), E), power(w, 2))))
        nonzero_root_assertion(subtract(power(y, 2), subtract(power(r, 2), power(subtract(x, one), 2))))
        if root[0].is_zero() or y[0].is_zero():
            raise ValueError("右端c6存在二次重根")
        Sx = add(times(x, 2), beta)
        ac = subtract(times(one, 2), divide(times(power(x, 2), 4), Sx))
        ec = divide(times(mul(x, y), -4), Sx)
        Z = subtract(subtract(times(x, 2), one), r)
        Pz = mul(subtract(v, Z), subtract(v, u))
        bc, be, bc0 = times(subtract(add(Z, times(v, 2)), u), -1), times(divide(subtract(power(w, 2), Pz), w), -1), mul(Z, subtract(times(v, 2), u))
        Qa = subtract(add(add(power(d, 2), power(e, 2)), add(mul(ac, d), mul(ec, e))), beta)
        B2 = add(add(power(d, 2), power(e, 2)), add(add(mul(bc, d), mul(be, e)), bc0))
        nonzero_root_assertion(Qa)
        nonzero_root_assertion(B2)
        ld, le = subtract(bc, ac), subtract(be, ec)
        if le[0].is_zero() or add(one, power(divide(ld, le), 2))[0].is_zero():
            raise ValueError("右端线性主元或二次首项为零")
        groups.setdefault(point_key([x, y]), set()).add(point_key(point))
        c6groups.setdefault(key(x), set()).add(key(y))
    if len(groups) != 4 or any(len(a) != 2 for a in groups.values()):
        raise ValueError("右端没有覆盖每个c6上的两个B2分支")
    if len(c6groups) != 2 or any(len(a) != 2 for a in c6groups.values()):
        raise ValueError("右端没有覆盖c6的两个二次方向")


def verify(data, field, system):
    zero = field(0)
    base = {j: (field(v), tuple(field(int(i == k)) for k in range(4)))
            for i, (j, v) in enumerate(zip((10, 11, 12, 13), data["parameters_u_v_w_r"]))}
    for name, count in (("left0", 8), ("left3", 32), ("left4", 128), ("right", 8)):
        if len(data[name]) != count:
            raise ValueError("纤维分支数改变")
    decoded = {}
    checks = 0
    for name, indices, own_indices in (("left0", (6, 10), (0, 1)),
                                        ("left3", (1, 5), (2, 3)),
                                        ("left4", (2, 7), (4, 5))):
        decoded[name] = []
        for entry in data[name]:
            point = [decode(v, field) for v in entry["point"]]
            decoded[name].append(point)
            values = dict(base)
            if name == "left3":
                parent = decoded["left0"][entry["parent"]]
                values.update(zip((0, 1), parent))
            if name == "left4":
                parent = decoded["left3"][entry["parent"]]
                values.update(zip((2, 3), parent))
            values.update(zip(own_indices, point))
            for index in indices:
                nonzero_root_assertion(row_value(system["equations"][index], values, field))
                checks += 1
    decoded["right"] = []
    for entry in data["right"]:
        point, c6 = [[decode(v, field) for v in entry[key]] for key in ("point", "c6")]
        decoded["right"].append(point)
        values = {**base, **dict(zip((6, 7), point)), **dict(zip((8, 9), c6))}
        for index in (0, 4, 8, 9):
            nonzero_root_assertion(row_value(system["equations"][index], values, field))
            checks += 1
    factor = read(FACTOR)
    check_left_completeness(data, decoded, base, field, factor)
    check_right_completeness(data, decoded, base, field)
    if len(data["bridge_values"]) != 1024:
        raise ValueError("闭环乘积不是完整1024项")
    product = field(1), (zero,)*4
    zeros = []
    for i, a in enumerate(decoded["left4"]):
        for j, b in enumerate(decoded["right"]):
            value = close_boundary(a, b, base[13], field)
            if value != decode(data["bridge_values"][8*i+j], field):
                raise ValueError("闭环值或偏导与独立重算不符")
            product = mul(product, value)
            if value[0].is_zero():
                zeros.append([i, j])
    if product[0] != field(data["norm_h3_mod_p"]) or product[1] != tuple(field(v) for v in data["norm_gradient_u_v_w_r_mod_p"]):
        raise ValueError("闭环范数或参数偏导错误")
    if zeros != data["vanishing_bridge_pairs"]:
        raise ValueError("零闭环分支清单错误")
    return {"verified": True, "original_constraint_rows_with_tangents_checked": checks,
            "all_endpoint_and_transfer_branches_verified": True,
            "scalar_equalities_checked": checks*5+1024*5+5,
            "bridge_factors_checked": 1024, "norm_h3_mod_p": data["norm_h3_mod_p"],
            "norm_gradient_mod_p": data["norm_gradient_u_v_w_r_mod_p"],
            "scope": "精确有限域纤维及一阶微分；不认证半径最小多项式"}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    source = HERE / "fiber-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        metadata = read(source / "field.json")
        if metadata["source_system_sha256"] != sha(AXIS) or metadata["left_factor_sha256"] != sha(FACTOR):
            raise ValueError("原整数系统改变")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        system = read(AXIS)
        check_boundary_identity(system)
        reports = []
        for path in sorted(source.glob("fiber-*.json.gz")):
            report = verify(read(path), field, system)
            report.update(source=path.name, source_sha256=sha(path))
            reports.append(report)
            print(json.dumps(report), flush=True)
        if not reports:
            raise ValueError("没有可核验纤维")
        write(source / args.output, {"verified": True, "fibers": reports,
                                          "h3_formula_verified_over_integers": True,
                                          "seconds": time.perf_counter()-tick, "limits": snapshot()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--output", default="verified-complete.json")
    run(parser.parse_args())
