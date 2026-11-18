"""由原几何纤维独立核验四参数范数，并在闭环零点核验全部有理恢复。"""

import argparse
import time

from descent_support import HERE, SHARED, FACTOR, AXIS, read, write, sha, bounded
from verify_multivariate import prune
from verify_system import execute
from verify_shared import check_fiber
from chain_fiber import Jet, prime_scalar
from fiber import compute


def jet(encoded, field):
    return Jet(field(encoded["value"]), [field(c) for c in encoded["gradient"]])


def leading_A(params, right, field):
    u, v, w, r = params
    x, y = [jet(row, field) for row in right["c6"]]
    Sx = 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be = -(Z+2*v-u), -(w*w-Pz)/w
    e1 = -(bc-ac)/(be-ec)
    return 1+e1*e1


def product_check(data, fiber, field, axis, reduced=None):
    source = read(SHARED / "system-runs/20261015-selected/circuit.json.gz")
    accepted = check_fiber(source, fiber, field, axis)
    params = [Jet(field(value), [field(int(i == j)) for j in range(4)])
              for i, value in enumerate(fiber["parameters_u_v_w_r"])]
    NA = Jet(field(1))
    count = 0
    for i in range(8):
        by_y = {}
        for left, right in accepted["selected_pairs"]:
            if left != i:
                continue
            row = fiber["right"][right]
            key = tuple(row["c6"][1]["value"])
            A = leading_A(params, row, field)
            if key in by_y:
                old = by_y[key]
                if (A.value, A.gradient) != (old.value, old.gradient):
                    raise ValueError("同一y根的两个d分支没有共享A及其导数")
            else:
                by_y[key] = A
        if len(by_y) != 2:
            raise ValueError("每个左端根没有恰好两个y根")
        for A in by_y.values():
            NA *= A
            count += 1
    NH = Jet(field(accepted["norm_h3"]), list(map(field, accepted["norm_gradient_u_v_w_r"])))
    expected = NH*NA
    outputs = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
    reduced, indices = reduced or prune(data, outputs)
    registers = execute(reduced, [p.value for p in params])
    if [registers[index] for index in indices] != [expected.value, *expected.gradient]:
        raise ValueError("四参数范数或全导数与原几何分支乘积不一致")
    closed = [(i, j) for i, j in accepted["selected_pairs"] if field(fiber["bridge_values"][8*i+j]["value"]).is_zero()]
    return {"verified": True, "original_geometry_points": 32, "A_factors": count,
            "G": prime_scalar(expected.value), "gradient_u_v_w_r": [prime_scalar(a) for a in expected.gradient],
            "closed_pairs": closed}


def closed_check(data, fiber, field, axis):
    product = product_check(data, fiber, field, axis)
    if len(product["closed_pairs"]) != 1:
        raise ValueError("闭环测试未取得唯一零分支")
    i, j = product["closed_pairs"][0]
    left, right = fiber["left"][i], fiber["right"][j]
    params = list(map(field, fiber["parameters_u_v_w_r"]))
    registers = execute(data, params)
    groups = ["all_divisor_nodes", "tower_leading_coefficient_nodes", "simple_root_derivative_nodes",
              "closed_chart_nonzero_guards", "critical_locus_nonzero_guards"]
    if any(registers[index].is_zero() for name in groups for index in data[name]):
        raise ValueError("闭环测试落到投影图表的坏主元")
    u, v, w, r = params
    t, dx, L = field(left["t"]["value"]), v-u, (v-u)**2+w*w
    aa, bb = dx-w*t/L, w+dx*t/L
    expected = {"s": bb/(2*r+aa), "z": field(left["q0"][0]["value"]),
                "h": -field(left["point"][1]["value"]), "y": field(right["c6"][1]["value"])}
    for name, index in data["outputs"]["recovered_left_and_y"].items():
        if registers[index] != expected[name]:
            raise ValueError("四参数图表的生成元恢复错误："+name)
    if registers[data["outputs"]["recovered_d"]] != field(right["point"][0]["value"]):
        raise ValueError("四参数图表的d恢复错误")
    geometry = {"c0": left["point"], "c3": left["c3"], "c4": left["c4"], "c5": right["point"], "c6": right["c6"]}
    for name, indices in data["outputs"]["geometry"].items():
        if [registers[index] for index in indices] != [field(a["value"]) for a in geometry[name]]:
            raise ValueError("四参数图表的原几何坐标恢复错误")
    bridge = fiber["bridge_values"][8*i+j]
    original_gradient = list(map(field, bridge["gradient"]))
    multiplier = registers[data["outputs"]["critical_multiplier"]]
    if [registers[index] for index in data["outputs"]["critical_equations"][1:]+[data["outputs"]["radius_derivative"]]] != [multiplier*a for a in original_gradient]:
        raise ValueError("闭环处四个临界导数的单位乘子错误")
    if not registers[data["outputs"]["critical_equations"][0]].is_zero():
        raise ValueError("闭环恢复后投影零式不为零")
    return {"verified": True, "parameters_u_v_w_r": fiber["parameters_u_v_w_r"],
            "unique_closed_pair": [i, j], "recovered_generators": 5, "recovered_geometry_coordinates": 10,
            "critical_total_derivatives_checked": 4, "nonzero_multiplier": list(map(int, multiplier.to_list()))}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    source = HERE / "system-runs" / args.system
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, axis, factor = read(source / "circuit.json.gz"), read(AXIS), read(FACTOR)
        reports, closed_reports, inputs = [], [], []
        for prime in (65521, 65519):
            folder = SHARED.parent / "rational-chain/fiber-runs" / f"20261015-p{prime}"
            metadata = read(folder / "field.json")
            field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(prime)(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
            outputs = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
            reduced = prune(data, outputs)
            for path in sorted(folder.glob("fiber-*.json.gz")):
                report = product_check(data, read(path), field, axis, reduced)
                report.update(prime=prime, fiber_sha256=sha(path))
                reports.append(report)
                inputs.append({"path": path.relative_to(HERE.parents[2]).as_posix(), "sha256": sha(path)})
            norm = read(SHARED / "quotient-runs" / f"20261015-jets-p{prime}" / "norm.json.gz")
            linear = next(row for row in norm["numerator_factors"] if row["degree"] == 1 and row["multiplicity"] == 1)
            radius = -linear["coefficients"][0]*pow(linear["coefficients"][1], -1, prime)%prime
            closed_field = fq_default_ctx(prime, 8, "a", fq_type="FQ_NMOD")
            fiber = compute(closed_field, [2, 3, 1, radius], factor)
            accepted = closed_check(data, fiber, closed_field, axis)
            accepted.update(prime=prime, field_modulus=list(map(int, closed_field.modulus().coeffs())))
            closed_reports.append(accepted)
            print({"prime": prime, "generic_fibers": 2, "closed_radius": radius, "closed_recovery_verified": True}, flush=True)
        report = {"verified": True, "generic_fibers": reports, "closed_recovery_tests": closed_reports,
                  "parent_inputs": inputs, "circuit_sha256": sha(source / "circuit.json.gz"), "verifier_sha256": sha(__file__),
                  "seconds": time.perf_counter()-tick, "limits": snapshot(),
                  "scope": "原几何有限域独立回放和闭环有理恢复测试；特征零等价性依赖通用恒等式，目标单根由独立区间证书验证。"}
        write(HERE / args.output, report)
        print({"verified": True, "generic_points": sum(row["original_geometry_points"] for row in reports),
               "closed_recovery_tests": len(closed_reports), "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", default="20261015-projected4-final")
    parser.add_argument("--output", default="geometry-verified.json.gz")
    run(parser.parse_args())
