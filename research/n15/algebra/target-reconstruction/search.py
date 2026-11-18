"""记录显式次数、系数高度、精度预算的PSLQ与整数格搜索。"""

import contextlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from refine import HERE, ROOT, sha, write_new
from windows_compute_limits import apply_limits

SEEDS = {
    "r": {"kind": "coordinate", "index": 13},
    "r2": {"kind": "square", "index": 13},
    "c0x": {"kind": "coordinate", "index": 0},
    "c12x": {"kind": "coordinate", "index": 11},
    "tangent6": {"kind": "ratio", "numerator": {"9": 1}, "denominator": {"constant": 1, "8": -1}},
    "tangent0": {"kind": "ratio", "numerator": {"1": 1}, "denominator": {"constant": 1, "0": 1}},
    "mixed": {"kind": "linear", "terms": {"13": 1, "0": 1, "1": 2}},
    "mobius_r": {"kind": "ratio", "numerator": {"constant": 1, "13": -1}, "denominator": {"constant": 1, "13": 1}},
    "affine_r": {"kind": "linear", "terms": {"constant": -3, "13": 8}},
}


def scalar_numeric(definition, values):
    def linear(terms):
        return sum(v if i == "constant" else v*values[int(i)] for i, v in terms.items())
    kind = definition["kind"]
    if kind == "coordinate":
        return values[definition["index"]]
    if kind == "square":
        return values[definition["index"]]**2
    if kind == "linear":
        return linear(definition["terms"])
    return linear(definition["numerator"])/linear(definition["denominator"])


def main():
    import mpmath as mp
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        root = json.loads((HERE / "refined-root.json").read_text(encoding="utf-8"))
        accepted = json.loads((HERE / "refined-root-verified.json").read_text(encoding="utf-8"))
        if not accepted["verified"] or accepted["refined_root_sha256"] != sha(HERE / "refined-root.json"):
            raise ValueError("尚未接受高精度目标根")
        binary = ROOT / ".tools/target-relation-lll.exe"
        chain = ROOT / ".tools/msolve-paged/chain/ucrt64"
        env = os.environ.copy()
        env["PATH"] = str(chain / "bin")+os.pathsep+env["PATH"]
        cases = [("r", 16), ("r", 32), ("r2", 32)]
        cases += [(seed, 64) for seed in SEEDS]
        cases += [("r", 128), ("r2", 128)]
        config = {"schema": "target-search-budget-v1", "seeds": SEEDS,
                  "coefficient_bound_inclusive": str(10**24),
                  "lll_delta": 0.99, "lll_eta": 0.51,
                  "lll_scale_digits_rule": "(degree+1)*(24+14)+50",
                  "lll_working_dps": 5760, "gram_schmidt_interval_bits": 1024,
                  "per_lll_case_wall_seconds": 180,
                  "cases": [{"seed": s, "degree": d, "scale_digits": (d+1)*38+50} for s, d in cases],
                  "pslq": [{"seed": "r", "degree": d, "dps": p,
                            "maxcoeff_exclusive": str(10**24+1), "tol_digits": p-60,
                            "maxsteps": 4000} for d, p in ((8, 360), (16, 720), (32, 1440))],
                  "root_sha256": sha(HERE / "refined-root.json"),
                  "binary_sha256": sha(binary), "source_sha256": sha(HERE / "relation_lll.c"),
                  "flint_dll_sha256": sha(chain / "bin/libflint-24.dll"),
                  "compiler_sha256": sha(chain / "bin/gcc.exe"),
                  "mpmath_pslq_source_sha256": sha(Path(mp.pslq.__code__.co_filename)),
                  "complete_groebner_basis_requested": False,
                  "primitive_element_claimed": False}
        write_new(HERE / "search-budget.json", config)
        output = HERE / "lattices"
        output.mkdir(exist_ok=True)
        runs = []
        for seed, degree in cases:
            name = f"{seed}-d{degree}"
            begin = time.perf_counter()
            digits = (degree+1)*38+50
            with mp.workdps(5760):
                values = [mp.mpf(v)/mp.mpf(root["Qx"]) for v in root["xnum"]]
                scalar = scalar_numeric(SEEDS[seed], values)
                scale = mp.mpf(10)**digits
                entries = [int(mp.nint(scale*scalar**i)) for i in range(degree+1)]
                scalar_text = mp.nstr(scalar, 50)
            input_path = output / (name+".input.txt")
            with input_path.open("x", encoding="ascii", newline="\n") as stream:
                stream.write(str(degree+1)+"\n"+"\n".join(map(str, entries))+"\n")
            data_path = output / (name+".json")
            with input_path.open("rb") as incoming, data_path.open("xb") as outgoing:
                try:
                    result = subprocess.run([str(binary)], stdin=incoming, stdout=outgoing,
                                            stderr=subprocess.PIPE, env=env, timeout=180)
                    status = "completed" if result.returncode == 0 else "failed"
                    error = result.stderr.decode("utf-8", errors="replace")
                except subprocess.TimeoutExpired:
                    status, error = "wall_timeout", "180 second cap"
            record = {"name": name, "seed": seed, "degree": degree, "scalar_decimal": scalar_text,
                      "coefficient_bound_inclusive": str(10**24), "scale_digits": digits,
                      "status": status, "seconds": time.perf_counter()-begin,
                      "input_sha256": sha(input_path), "data_sha256": sha(data_path), "error": error}
            write_new(output / (name+".run.json"), record)
            runs.append(record)
            print(json.dumps({k: record[k] for k in ("name", "status", "seconds")}), flush=True)
        pslq_runs = []
        for case in config["pslq"]:
            begin = time.perf_counter()
            with mp.workdps(case["dps"]):
                value = mp.mpf(root["xnum"][13])/mp.mpf(root["Qx"])
                vector = [value**i for i in range(case["degree"]+1)]
                captured = io.StringIO()
                with contextlib.redirect_stdout(captured):
                    candidate = mp.pslq(vector, tol=mp.mpf(10)**(-case["tol_digits"]),
                                        maxcoeff=int(case["maxcoeff_exclusive"]), maxsteps=case["maxsteps"], verbose=True)
                lines = captured.getvalue().splitlines()
                iterations = [int(m.group(1)) for line in lines if (m := re.match(r"(\d+)/\d+:.*Error", line))]
                record = {**case, "candidate_low_to_high": candidate,
                          "iterations_observed": (max(iterations)+1) if iterations else None,
                          "last_diagnostics": lines[-3:], "seconds": time.perf_counter()-begin,
                          "absence_is_a_proof": False}
                pslq_runs.append(record)
            print(json.dumps(record), flush=True)
        write_new(HERE / "search-runs.json",
                  {"schema": "target-search-runs-v1", "lll": runs, "pslq": pslq_runs,
                   "seconds": time.perf_counter()-started, "windows_job_limits": limits()})
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
