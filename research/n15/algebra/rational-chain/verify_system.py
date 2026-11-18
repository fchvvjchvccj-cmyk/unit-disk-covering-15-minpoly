"""用已接受的几何纤维逐点重验九变量电路的关系、坐标和全导数。"""

import argparse
import time

from support import HERE, AXIS, FACTOR, read, write, sha, bounded


def execute(data, values):
    registers = []
    for i, node in enumerate(data["nodes"]):
        indices = node.get("args", [])
        if any(not isinstance(j, int) or not 0 <= j < i for j in indices):
            raise ValueError("电路不是顺序无环依赖")
        args = [registers[j] for j in indices]
        op = node["op"]
        if op == "var":
            result = values[node["value"]]
        elif op == "const":
            result = values[0]*0+node["value"]
        elif op == "neg" and len(args) == 1:
            result = -args[0]
        elif op == "add" and len(args) == 2:
            result = args[0]+args[1]
        elif op == "mul" and len(args) == 2:
            result = args[0]*args[1]
        elif op == "div" and len(args) == 2:
            result = args[0]/args[1]
        else:
            raise ValueError("电路指令或元数非法")
        registers.append(result)
    return registers


def verify(data, fiber, field):
    u, v, w, r = map(field, fiber["parameters_u_v_w_r"])
    L, dx = (v-u)**2+w*w, v-u
    points, checks = 0, 0
    for i, left in enumerate(fiber["left"]):
        t, k = field(left["t"]["value"]), field(left["k"]["value"])
        aa, bb = dx-w*t/L, w+dx*t/L
        s = bb/(2*r+aa)
        for j, right in enumerate(fiber["right"]):
            x, y = [field(a["value"]) for a in right["c6"]]
            d = field(right["point"][0]["value"])
            if data["variables"] == ["s", "k", "x", "y", "d", "u", "v", "w", "r"]:
                values = [s, k, x, y, d, u, v, w, r]
            elif data["variables"] == ["s", "z", "h", "x", "y", "d", "u", "v", "w", "r"]:
                values = [s, field(left["q0"][0]["value"]), -field(left["point"][1]["value"]), x, y, d, u, v, w, r]
            else:
                raise ValueError("未知的接触链端口次序")
            registers = execute(data, values)
            guards = data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+data["simple_root_derivative_nodes"]
            if any(registers[n].is_zero() for n in guards):
                raise ValueError("电路纤维进入分母、首项或重根坏点")
            for n in data["outputs"]["finite_relations"]:
                if not registers[n].is_zero():
                    raise ValueError("五级载体关系在几何纤维上不为零")
                checks += 1
            expected = {"c0": left["point"], "c3": left["c3"], "c4": left["c4"],
                        "c5": right["point"], "c6": right["c6"]}
            for name, nodes in data["outputs"]["geometry"].items():
                for n, value in zip(nodes, expected[name]):
                    if registers[n] != field(value["value"]):
                        raise ValueError("九变量重建坐标与独立纤维不同")
                    checks += 1
            bridge = fiber["bridge_values"][8*i+j]
            target = [field(bridge["value"])]+[field(a) for a in bridge["gradient"]]
            nodes = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
            if [registers[n] for n in nodes] != target:
                raise ValueError("九变量电路的闭环全导数与原几何纤维不同")
            checks += 5
            points += 1
    return {"verified": True, "tower_points_checked": points, "scalar_equalities_checked": checks,
            "scope": "验证指定有限域纤维上的电路；不将有限样本比较当作特征零恒等式证明"}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    source = HERE / "system-runs" / args.system
    fiber_folder = HERE / "fiber-runs" / args.fiber
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, metadata = read(source / "circuit.json.gz"), read(fiber_folder / "field.json")
        if data["source_system_sha256"] != sha(AXIS) or data["left_factor_sha256"] != sha(FACTOR):
            raise ValueError("九变量电路输入摘要不符")
        accepted = read(fiber_folder / "verified.json")
        if not accepted["verified"]:
            raise ValueError("几何纤维未通过独立接受端")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        reports = []
        for path in sorted(fiber_folder.glob("fiber-*.json.gz")):
            report = verify(data, read(path), field)
            report.update(fiber=path.name, sha256=sha(path))
            reports.append(report)
        write(source / args.output, {"verified": True, "reports": reports,
                                    "circuit_sha256": sha(source / "circuit.json.gz"),
                                    "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"verified": True, "points": sum(r["tower_points_checked"] for r in reports),
               "seconds": time.perf_counter()-tick}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--fiber", required=True)
    parser.add_argument("--output", default="fibers-verified.json")
    run(parser.parse_args())
