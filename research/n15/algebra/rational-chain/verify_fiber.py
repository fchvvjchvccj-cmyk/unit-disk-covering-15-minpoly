"""独立核验有理链纤维、原几何方程及六十四项闭环微分。"""

import argparse
import time

from support import HERE, AXIS, FACTOR, read, write, sha, bounded
import importlib.util

spec = importlib.util.spec_from_file_location("old_fiber_acceptor", HERE.parent / "chain-radius/verify_fiber.py")
OLD = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OLD)
decode, add, mul, times, power, subtract, divide = OLD.decode, OLD.add, OLD.mul, OLD.times, OLD.power, OLD.subtract, OLD.divide
zero_assert = OLD.nonzero_root_assertion


def verify(data, field, system):
    zero, one = (field(0), (field(0),)*4), (field(1), (field(0),)*4)
    base = {j: (field(v), tuple(field(int(i == k)) for k in range(4)))
            for i, (j, v) in enumerate(zip((10, 11, 12, 13), data["parameters_u_v_w_r"]))}
    u, v, w, r = [base[i] for i in (10, 11, 12, 13)]
    if len(data["left"]) != 8 or len(data["right"]) != 8 or len(data["bridge_values"]) != 64:
        raise ValueError("有理链载体维数错误")
    left, right = [], []
    left0 = []
    rows = 0
    for entry in data["left"]:
        c0, c3, c4, q0, p0 = [[decode(v, field) for v in entry[name]] for name in ("point", "c3", "c4", "q0", "p0")]
        left0.append(c0)
        a = [c0[0], times(c0[1], -1)]
        T = add(power(a[0], 2), power(a[1], 2))
        K = subtract(add(one, T), power(r, 2))
        t = subtract(times(a[0], 2), u)
        zero_assert(subtract(times(mul(t, q0[0]), 2), subtract(add(power(t, 2), one), power(r, 2))))
        zero_assert(subtract(times(mul(a[1], q0[1]), 2), subtract(K, times(mul(a[0], q0[0]), 2))))
        for p in (a, [u, zero], [v, w]):
            distance = add(power(subtract(p0[0], p[0]), 2), power(subtract(p0[1], p[1]), 2))
            zero_assert(subtract(distance, power(r, 2)))
        if any(z[0].is_zero() for z in (t, a[1], T)):
            raise ValueError("旧边界或第一次重建分母为零")
        area = subtract(mul(subtract(u, a[0]), subtract(w, a[1])), mul(times(a[1], -1), subtract(v, a[0])))
        if area[0].is_zero():
            raise ValueError("旧三角形外心不唯一")
        for i, b in enumerate((v, w)):
            expected = subtract(subtract(add(b, divide(mul(K, a[i]), T)), q0[i]), p0[i])
            zero_assert(subtract(c3[i], expected))
        TT = add(power(c3[0], 2), power(c3[1], 2))
        for i, b in enumerate((v, w)):
            expected = subtract(add(b, divide(mul(subtract(one, power(r, 2)), c3[i]), TT)), a[i])
            zero_assert(subtract(c4[i], expected))
        values = {**base, **dict(zip(range(6), c0+c3+c4))}
        for index in (1, 2, 5, 6, 7, 10):
            zero_assert(OLD.row_value(system["equations"][index], values, field))
            rows += 1
        left.append(c4)
    # 以原四次载体核验左端两个二次方向、全部四个单根及坐标映射。
    factor = read(FACTOR)
    quartic = factor["factor_polynomials"][factor["target_selected_factor_index"]]
    dx, L = subtract(v, u), add(power(subtract(v, u), 2), power(w, 2))
    eta = mul(L, subtract(times(power(r, 2), 4), L))
    groups = {}
    for entry, c0 in zip(data["left"], left0):
        t, k = [decode(entry[n], field) for n in ("t", "k")]
        zero_assert(subtract(power(t, 2), eta))
        aa, bb = subtract(dx, divide(mul(w, t), L)), add(w, divide(mul(dx, t), L))
        s = divide(bb, add(times(r, 2), aa))
        zero_assert(OLD.row_value(quartic, {0: k, 1: s, 2: u, 3: r}, field))
        lead = [dict(row, powers=[[j, e] for j, e in row["powers"] if j != 0]) for row in quartic if dict(row["powers"]).get(0, 0) == 4]
        if OLD.row_value(lead, {1: s, 2: u, 3: r}, field)[0].is_zero():
            raise ValueError("左端四次首项消失")
        X = divide(subtract(aa, mul(bb, k)), add(one, power(k, 2)))
        zero_assert(subtract(c0[0], add(u, X)))
        zero_assert(subtract(c0[1], mul(k, X)))
        groups.setdefault(OLD.key(t), set()).add(OLD.key(k))
    if len(groups) != 2 or any(len(roots) != 4 for roots in groups.values()):
        raise ValueError("左端载体根不完整或重复")
    for entry in data["right"]:
        point, c6 = [[decode(v, field) for v in entry[key]] for key in ("point", "c6")]
        right.append(point)
        values = {**base, **dict(zip((6, 7), point)), **dict(zip((8, 9), c6))}
        for index in (0, 4, 8, 9):
            zero_assert(OLD.row_value(system["equations"][index], values, field))
            rows += 1
    OLD.check_right_completeness(data, {"right": right}, base, field)
    product, zeros = one, []
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            value = OLD.close_boundary(a, b, r, field)
            if value != decode(data["bridge_values"][8*i+j], field):
                raise ValueError("有理链闭环值或微分错误")
            product = mul(product, value)
            if value[0].is_zero():
                zeros.append([i, j])
    if product != (field(data["norm_h3_mod_p"]), tuple(field(v) for v in data["norm_gradient_u_v_w_r_mod_p"])) or zeros != data["vanishing_bridge_pairs"]:
        raise ValueError("有理链范数或微分错误")
    return {"verified": True, "original_constraint_rows_with_tangents": rows, "complete_selected_rank": 64,
            "bridge_values_with_tangents_checked": 64, "norm_h3_mod_p": data["norm_h3_mod_p"],
            "norm_gradient_u_v_w_r_mod_p": data["norm_gradient_u_v_w_r_mod_p"],
            "scope": "目标分支所对应的精确有限域载体，不是半径最小多项式"}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    source = HERE / "fiber-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        metadata = read(source / "field.json")
        if metadata["source_system_sha256"] != sha(AXIS) or metadata["left_factor_sha256"] != sha(FACTOR):
            raise ValueError("有理链原始输入摘要改变")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        system = read(AXIS)
        OLD.check_boundary_identity(system)
        reports = []
        for path in sorted(source.glob("fiber-*.json.gz")):
            report = verify(read(path), field, system)
            report.update(source=path.name, sha256=sha(path))
            reports.append(report)
            print(report, flush=True)
        if not reports:
            raise ValueError("没有可重放的有理链纤维")
        write(source / args.output, {"verified": True, "reports": reports, "seconds": time.perf_counter()-tick,
                                    "limits": snapshot(), "verifier_sha256": sha(__file__)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--output", default="verified.json")
    run(parser.parse_args())
