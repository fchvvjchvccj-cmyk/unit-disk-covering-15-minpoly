"""从父25元盒迭代并产生更高精度的候选；认证另由整数区间接受器完成。"""

import hashlib
import gzip
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from windows_compute_limits import apply_limits


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def artifact_bytes(path):
    if path.exists():
        return path.read_bytes()
    compressed = path.with_name(path.name+".gz")
    with gzip.open(compressed, "rb") as stream:
        return stream.read()


def artifact_sha(path):
    return hashlib.sha256(artifact_bytes(path)).hexdigest()


def write_new(path, data):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def main():
    import mpmath as mp
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        system_path = HERE.parent / "minpoly-final/axis25-system.json"
        parent_path = HERE.parent / "minpoly-final/axis25-root.json"
        system = json.loads(system_path.read_text(encoding="utf-8"))
        parent = json.loads(parent_path.read_text(encoding="utf-8"))
        symbols = system["variables"]
        terms = [[(t["coefficient"], t["powers"]) for t in row]
                 for row in system["equations"]]
        records = []
        with mp.workdps(100):
            x = [mp.mpf(v)/mp.mpf(parent["Qx"]) for v in parent["xnum"]]

        def evaluate(values):
            powers = [[mp.mpf(1)] + [v**k for k in range(1, 7)] for v in values]
            f = []
            jacobian = [[mp.mpf(0) for _ in values] for _ in values]
            for j, polynomial in enumerate(terms):
                total = mp.mpf(0)
                for coefficient, exponent in polynomial:
                    term = mp.mpf(coefficient)
                    for i, e in exponent:
                        term *= powers[i][e]
                    total += term
                    for column, degree in exponent:
                        term = mp.mpf(coefficient*degree)
                        for i, e in exponent:
                            term *= powers[i][e-(i == column)]
                        jacobian[j][column] += term
                f.append(total)
            return mp.matrix(f), mp.matrix(jacobian)

        for dps in (180, 360, 720, 1440, 2880, 5760):
            with mp.workdps(dps):
                for iteration in range(3):
                    f, jacobian = evaluate(x)
                    correction = mp.lu_solve(jacobian, f)
                    x = [v-d for v, d in zip(x, correction)]
                    record = {"dps": dps, "iteration": iteration+1,
                              "residual_before": mp.nstr(max(abs(v) for v in f), 12),
                              "correction": mp.nstr(max(abs(v) for v in correction), 12)}
                    records.append(record)
                print(json.dumps(records[-1]), flush=True)

        with mp.workdps(5760):
            scale = 10**5500
            center = [int(mp.nint(v*scale)) for v in x]
            values = [mp.mpf(v)/scale for v in center]
            f, jacobian = evaluate(values)
            inverse = jacobian**-1
            inverse_scale = 10**150
            root = {"schema": "target-reconstruction-root-v1",
                    "system_sha256": sha(system_path),
                    "parent_root_sha256": sha(parent_path),
                    "variables": symbols, "Qx": str(scale),
                    "xnum": list(map(str, center)), "Qy": str(inverse_scale),
                    "Ynum": [[str(int(mp.nint(v*inverse_scale))) for v in row]
                             for row in inverse.tolist()],
                    "rho": "1/"+str(10**5400), "interval_bits": 19000,
                    "working_dps": 5760,
                    "numerical_construction_is_proof": False}
            write_new(HERE / "refined-root.json", root)
            write_new(HERE / "numerical-refinement.json",
                      {"schema": "target-refinement-run-v1", "iterations": records,
                       "refined_root_sha256": sha(HERE / "refined-root.json"),
                       "residual_at_rational_center": mp.nstr(max(abs(v) for v in f), 30),
                       "radius_decimal": mp.nstr(values[13], 5501),
                       "mpmath_version": mp.__version__, "seconds": time.perf_counter()-started,
                       "windows_job_limits": limits()})
        print(json.dumps({"completed": True, "seconds": time.perf_counter()-started}), flush=True)
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
