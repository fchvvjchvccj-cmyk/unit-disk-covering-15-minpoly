"""逐级放大格精度并沿用整数变换，避免高维搜索每次从单位基开始。"""

import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import time

from common import HERE, TARGET, bounded, read, write, sha, payload


def worker(args):
    from flint import fmpz_mat, __version__
    incoming = list(map(int, payload(args.input).decode("ascii").split()))
    n, entries = incoming[0], incoming[1:]
    if n != len(entries):
        raise ValueError("格输入维数不符")
    digits = len(str(entries[0]))-1
    U = fmpz_mat([[int(i == j) for j in range(n)] for i in range(n)])
    records = []
    for d in sorted(set([min(digits, a) for a in (1000, 2000, 3000, 4000, digits)])):
        tick = time.perf_counter()
        truncated = [v//10**(digits-d) for v in entries]
        original = fmpz_mat([[int(i == j) for j in range(n)]+[truncated[i]] for i in range(n)])
        basis, transform = (U*original).lll(transform=True, delta=0.99, eta=0.51)
        U = transform*U
        # 每一阶段的变换都与原单位基关联；最终仍由整数接受器重验。
        if U*original != basis:
            raise ValueError("分级格变换没有正确复合")
        record = {"digits": d, "seconds": time.perf_counter()-tick}
        records.append(record)
        write(Path(args.output).parent / f"stage-{d}.json", record)
        write(Path(args.output).parent / f"transform-{d}.json.gz",
              {"digits": d, "transform": [[str(U[i, j]) for j in range(n)] for i in range(n)]})
        print(record, flush=True)
    inverse = U.inv()
    write(args.output, {"flint_version": "python-flint "+__version__,
                        "basis": [[str(basis[i, j]) for j in range(n+1)] for i in range(n)],
                        "inverse_transform": [[str(inverse[i, j]) for j in range(n)] for i in range(n)],
                        "stages": records})


def run(args):
    from height_search import verify
    parent = HERE / "height-runs" / args.parent
    root = read(TARGET / "refined-root.json")
    config = read(parent / "budget.json")
    cases = [case for case in config["cases"] if case["degree"] in ((args.degree,) if args.degree else (96, 128))]
    if not cases:
        raise ValueError("指定次数不在母实验预算中")
    out = HERE / "height-runs" / args.tag
    if out.exists():
        raise FileExistsError("拒绝覆盖已有分级搜索")
    out.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        records = []
        for case in cases:
            case_dir = out / case["name"]
            case_dir.mkdir()
            incoming = parent / (case["name"]+".input.txt")
            result_path = case_dir / "basis.json"
            begin = time.perf_counter()
            with (case_dir / "stdout.txt").open("xb") as log:
                try:
                    process = subprocess.run([sys.executable, "-B", __file__, "--worker", "--input", str(incoming),
                                              "--output", str(result_path)], stdout=log, stderr=subprocess.STDOUT,
                                             timeout=args.seconds)
                    status = "completed" if process.returncode == 0 else "worker_failed"
                except subprocess.TimeoutExpired:
                    status = "wall_timeout"
            record = {**case, "status": status, "worker_seconds": time.perf_counter()-begin,
                      "input_sha256": hashlib.sha256(payload(incoming)).hexdigest()}
            if status == "completed":
                accepted = verify(incoming, result_path, root, case)
                write(case_dir / "verified.json.gz", accepted)
                record.update(verified=True, excluded=accepted["no_nonzero_integer_polynomial_in_budget_proved"],
                              possible_candidates=len(accepted["candidates_needing_exact_algebraic_validation"]),
                              output_sha256=sha(result_path))
            record["total_seconds"] = time.perf_counter()-begin
            records.append(record)
            write(case_dir / "run.json", record)
            print(record, flush=True)
        write(out / "summary.json", {"runs": records, "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                      "driver_sha256": sha(__file__), "source_budget_sha256": sha(parent / "budget.json"),
                                      "minimal_polynomial_obtained": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--parent", default="20261015-wide")
    parser.add_argument("--tag")
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--degree", type=int, choices=(96, 128))
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        if not args.tag:
            parser.error("需要新的运行标签")
        run(args)
