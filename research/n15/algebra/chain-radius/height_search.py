"""扩大次数和高度的二维搜索；精确接受端独立检查整数格与真实根。"""

import argparse
from fractions import Fraction
import json
import subprocess
import sys
import time

from common import HERE, TARGET, ROOT, bounded, read, write, sha, payload, binary_environment


CASES = [
    ("r", 32, 100), ("r2", 32, 100),
    ("r", 64, 60), ("r2", 64, 60),
    ("r", 96, 40), ("r", 128, 24),
]


def verify(input_path, data_path, root, config):
    from intervals import Dyadic
    from verify_lattices import gs_lower_bounds, polynomial_interval
    entries = list(map(int, payload(input_path).decode("ascii").split()))
    size, entries = entries[0], entries[1:]
    if size != config["degree"]+1 or len(entries) != size:
        raise ValueError("格输入的维数错误")
    proposal = read(data_path)
    basis = [list(map(int, row)) for row in proposal["basis"]]
    inverse = [list(map(int, row)) for row in proposal["inverse_transform"]] if "inverse_transform" in proposal else None
    if len(basis) != size or any(len(row) != size+1 for row in basis):
        raise ValueError("返回基维数错误")
    if inverse is not None and (len(inverse) != size or any(len(row) != size for row in inverse)):
        raise ValueError("返回逆变换维数错误")
    for i, row in enumerate(basis):
        if sum(a*b for a, b in zip(row[:-1], entries)) != row[-1]:
            raise ValueError("返回行不属于输入整数格")
        if inverse is not None:
            for j in range(size):
                if sum(row[k]*inverse[k][j] for k in range(size)) != int(i == j):
                    raise ValueError("返回变换不是整数幺模变换")
    if inverse is None:
        from flint import fmpz_mat
        if abs(int(fmpz_mat([row[:-1] for row in basis]).det())) != 1:
            raise ValueError("紧凑格证书的整数变换行列式不是正负1")
    arithmetic = Dyadic(root["interval_bits"])
    center = Fraction(int(root["xnum"][13]), int(root["Qx"]))
    rho = Fraction(root["rho"])
    scalar = (arithmetic.rational(center-rho)[0], arithmetic.rational(center+rho)[1])
    if config["seed"] == "r2":
        scalar = arithmetic.power(scalar, 2)
    elif config["seed"] != "r":
        raise ValueError("未知搜索元")
    power, max_error = arithmetic.one, 0
    for entry in entries:
        error = arithmetic.absolute(arithmetic.sub(arithmetic.integer(power, 10**config["scale_digits"]), arithmetic.rational(entry)))
        max_error = max(max_error, error)
        if error > arithmetic.scale:
            raise ValueError("输入幂的认证误差超过1")
        power = arithmetic.mul(power, scalar)
    diagonals, scale = gs_lower_bounds(basis, 1536)
    minimum = min(a for a, b in diagonals)
    height = 10**config["height_digits"]
    threshold = (size+size**2)*height**2
    possible, rejected = [], []
    for row in sorted(basis, key=lambda a: sum(v*v for v in a)):
        coefficients = row[:-1]
        value = polynomial_interval(coefficients, scalar, arithmetic)
        record = {"coefficients_low_to_high": list(map(str, coefficients)),
                  "height": str(max(map(abs, coefficients))),
                  "inside_height_budget": max(map(abs, coefficients)) <= height,
                  "value_interval": [str(Fraction(v, arithmetic.scale)) for v in value]}
        if value[0] <= 0 <= value[1]:
            possible.append(record)
        elif not rejected:
            rejected.append(record)
    return {"verified": True, "basis_and_unimodularity_verified": True,
            "unimodularity_method": "integer_inverse" if inverse is not None else "exact_integer_determinant",
            "power_rounding_error_upper": str(Fraction(max_error, arithmetic.scale)),
            "gs_squared_norm_lower": str(Fraction(minimum, scale)),
            "relation_squared_norm_upper_if_zero": str(threshold),
            "no_nonzero_integer_polynomial_in_budget_proved": minimum > threshold*scale,
            "candidates_needing_exact_algebraic_validation": possible,
            "shortest_rejected_basis_row": rejected[0] if rejected else None,
            "flint_version": proposal["flint_version"]}


def run(args):
    import mpmath as mp
    from verify_root import verify as verify_root
    out = HERE / "height-runs" / args.tag
    if out.exists():
        raise FileExistsError("运行目录已经存在，拒绝覆盖："+str(out))
    out.mkdir(parents=True)
    root = read(TARGET / "refined-root.json")
    binary = ROOT / ".tools/target-relation-lll.exe"
    with bounded() as snapshot:
        tick = time.perf_counter()
        root_report = verify_root(root)
        write(out / "root-replayed.json", root_report)
        cases = [{"name": f"{s}-d{d}-h{h}", "seed": s, "degree": d, "height_digits": h,
                  "scale_digits": (d+1)*(h+8)+50} for s, d, h in CASES]
        config = {"schema": "chain-radius-height-search-v1", "cases": cases,
                  "per_case_wall_seconds": args.seconds, "root_sha256": sha(TARGET / "refined-root.json"),
                  "source_system_sha256": sha(ROOT / "n15/algebra/minpoly-final/axis25-system.json"),
                  "binary_sha256": sha(binary), "binary_source_sha256": sha(TARGET / "relation_lll.c"),
                  "driver_sha256": sha(__file__), "coefficient_height_inclusive": True,
                  "numerical_relation_is_exact_proof": False}
        write(out / "budget.json", config)
        records = []
        for case in cases:
            begin = time.perf_counter()
            name = case["name"]
            incoming, outgoing = out / (name+".input.txt"), out / (name+".basis.json")
            with mp.workdps(5760):
                value = mp.mpf(root["xnum"][13])/mp.mpf(root["Qx"])
                if case["seed"] == "r2":
                    value *= value
                scale = mp.mpf(10)**case["scale_digits"]
                entries = [int(mp.nint(scale*value**i)) for i in range(case["degree"]+1)]
            incoming.write_text(str(len(entries))+"\n"+"\n".join(map(str, entries))+"\n", encoding="ascii")
            print(json.dumps({"case": name, "status": "started"}), flush=True)
            with incoming.open("rb") as src, outgoing.open("xb") as dst:
                try:
                    result = subprocess.run([str(binary)], stdin=src, stdout=dst, stderr=subprocess.PIPE,
                                            env=binary_environment(), timeout=args.seconds)
                    status = "completed" if result.returncode == 0 else "worker_failed"
                    error = result.stderr.decode("utf-8", errors="replace")
                except subprocess.TimeoutExpired:
                    status, error = "wall_timeout", "单病例达到墙钟上限；未完成范围不作排除"
            record = {**case, "status": status, "worker_seconds": time.perf_counter()-begin, "error": error,
                      "input_sha256": sha(incoming), "output_sha256": sha(outgoing)}
            if status == "completed":
                validation = verify(incoming, outgoing, root, case)
                write(out / (name+".verified.json.gz"), validation)
                record.update(verified=True,
                              excluded=validation["no_nonzero_integer_polynomial_in_budget_proved"],
                              possible_candidates=len(validation["candidates_needing_exact_algebraic_validation"]))
            record["total_seconds"] = time.perf_counter()-begin
            write(out / (name+".run.json"), record)
            records.append(record)
            print(json.dumps(record, ensure_ascii=False), flush=True)
        write(out / "summary.json", {"schema": "chain-radius-height-summary-v1", "runs": records,
                                      "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                      "minimal_polynomial_obtained": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--seconds", type=float, default=180)
    run(parser.parse_args())
