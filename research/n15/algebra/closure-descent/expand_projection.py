"""在受限资源内实测四参数标量闭环的平坦多项式展开，不把停止报告当作消元结果。"""

import argparse
from pathlib import Path
import subprocess
import sys
import time

from descent_support import HERE, read, write, sha, bounded, binary_environment
from multivariate_probe import MultivariateRational as M
from verify_multivariate import prune
from verify_system import execute
import tower


def worker(args):
    from flint import nmod_mpoly_ctx
    target, tick = Path(args.worker), time.perf_counter()
    path = HERE / "system-runs" / args.system / "circuit.json.gz"
    data = read(path)
    M.context = nmod_mpoly_ctx.get(args.parameters, args.prime, ordering="lex")
    M.prime, M.term_limit, M.deadline = args.prime, args.max_terms, tick+args.seconds
    M.zero_poly, M.one_poly = M.context.constant(0), M.context.constant(1)
    M.metrics = {"max_coefficient_degree": 0, "max_coefficient_terms": 0, "rational_additions": 0, "rational_multiplications": 0}
    tower.Rational = M
    fixed = {"u": 2, "v": 3, "w": 1, "r": 2}
    values = [M(M.context.gen(args.parameters.index(name))) if name in args.parameters else M(fixed[name])
              for name in data["variables"]]
    small, indices = prune(data, [data["outputs"]["critical_equations"][0]])
    report = {"prime": args.prime, "symbolic_parameters": args.parameters, "fixed_parameters": fixed,
              "pruned_nodes": len(small["nodes"]), "circuit_sha256": sha(path), "driver_sha256": sha(__file__),
              "radius_minimal_polynomial_obtained": False}
    print({"stage": "start", "pruned_nodes": len(small["nodes"])}, flush=True)
    try:
        result = execute(small, values)[indices[0]]
        write(target / "function.json.gz", {"symbolic_parameters": args.parameters, "prime": args.prime,
              "value": result.encode(), "radius_minimal_polynomial_obtained": False})
        report.update(status="completed", numerator_terms=len(result.n), denominator_terms=len(result.d),
                      numerator_total_degree=int(result.n.total_degree()), denominator_total_degree=int(result.d.total_degree()))
    except (tower.TowerStop, ZeroDivisionError) as error:
        report.update(status="stopped", reason=str(error))
    report.update(seconds=time.perf_counter()-tick, **M.metrics)
    write(target / "worker.json", report)
    print(report, flush=True)


def run(args):
    target = HERE / "expansion-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有全参数展开试验")
    target.mkdir(parents=True)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(target), "--system", args.system,
               "--prime", str(args.prime), "--parameters", *args.parameters,
               "--seconds", str(args.seconds), "--max-terms", str(args.max_terms)]
    with bounded() as snapshot:
        tick = time.perf_counter()
        with (target / "worker.log").open("x", encoding="utf-8") as stream:
            try:
                child = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                       env=binary_environment(), timeout=args.seconds+5, check=False)
                status, returncode = "worker_completed", child.returncode
            except subprocess.TimeoutExpired:
                status, returncode = "hard_timeout", None
        report = {"controller_status": status, "returncode": returncode, "seconds": time.perf_counter()-tick,
                  "limits": snapshot(), "log_sha256": sha(target / "worker.log")}
        if (target / "worker.json").exists():
            report["worker"] = read(target / "worker.json")
        write(target / "controller.json", report)
        print({key: value for key, value in report.items() if key not in ("limits", "worker")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", default="20261015-projected4")
    parser.add_argument("--parameters", nargs="+", default=["u", "v", "w", "r"])
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--seconds", type=int, default=60)
    parser.add_argument("--max-terms", type=int, default=100000)
    parser.add_argument("--tag")
    parser.add_argument("--worker")
    args = parser.parse_args()
    worker(args) if args.worker else run(args)
