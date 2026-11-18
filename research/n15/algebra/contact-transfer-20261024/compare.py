"""在相同素域参数线上完整比较有限载体与四参数范数，不运行通用求基。"""

import argparse
from collections import Counter
import json
import subprocess
import sys
import time

from context import HERE, INPUTS, RESULTS, read, write, sha, bounded, prune, execute

CASES = {"primitive16": RESULTS / "primitive16.json.gz", "descent16": INPUTS / "descent16.json.gz",
         "primitive4": INPUTS / "primitive4.json.gz"}


def statistics(data):
    outputs = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
    value, _ = prune(data, outputs[:1])
    critical, _ = prune(data, outputs)
    return {"variables": len(data["variables"]), "stored_nodes": len(data["nodes"]),
            "closure_nodes": len(value["nodes"]), "closure_operations": dict(Counter(n["op"] for n in value["nodes"])),
            "critical_nodes": len(critical["nodes"]), "critical_operations": dict(Counter(n["op"] for n in critical["nodes"]))}


def initialize(prime, seconds):
    import tower
    from flint import nmod_poly
    R = tower.Rational
    R.prime = prime
    R.zero_poly, R.one_poly = nmod_poly([], prime), nmod_poly([1], prime)
    R.deadline, R.degree_limit = time.perf_counter()+seconds, 4096
    R.metrics = {"max_coefficient_degree": 0, "rational_additions": 0, "rational_multiplications": 0}
    return tower, R


def carrier(case, prime, seconds):
    from flint import nmod_poly
    tower, R = initialize(prime, seconds)
    params = [R(2), R(3), R(1), R(nmod_poly([0, 1], prime))]
    if case == "primitive4":
        return None, params
    A = tower.Algebra()
    u, v, w, r = [A.scalar(x) for x in params]
    dx, L = v-u, (v-u)**2+w*w
    if case == "primitive16":
        first = A.add_level(L*(4*r*r-L), 0)
        aa, bb = dx-w*first/L, w+dx*first/L
    else:
        lead = L+2*r*dx
        first = A.add_level(-(L-2*r*dx)/lead, 4*r*w/lead)
        aa, bb = 2*r*(1-first*first)/(1+first*first), 4*r*first/(1+first*first)
    px, py = u+aa/2, bb/2
    z = A.add_level(1-py*py, 0)
    ax = (px+z)/2
    h = A.add_level(r*r-(ax-px)**2-py*py, 2*py)
    x = (1+r+v+dx-aa)/2
    y = A.add_level(r*r-(x-1)**2, 0)
    return A, [first, z, h, y, u, v, w, r]


def trace(A, coefficients, level):
    if level == 0:
        return coefficients[0]
    half = len(coefficients)//2
    a, b = coefficients[:half], coefficients[half:]
    relation_b = A.relations[level-1][1]
    lower = A.add(tuple(2*x for x in a), A.multiply(b, relation_b, level-1))
    return trace(A, lower, level-1)


def worker(args):
    tick = time.perf_counter()
    path = RESULTS / (args.tag+"-"+args.case+"-p"+str(args.prime)+".json.gz")
    data = read(CASES[args.case])
    outputs = data["outputs"]["critical_equations"]+[data["outputs"]["radius_derivative"]]
    small, outputs = prune(data, outputs)
    report = {"case": args.case, "prime": args.prime, "fixed_u_v_w": [2, 3, 1],
              "source_sha256": sha(CASES[args.case]), "statistics": statistics(data),
              "scope": "精确Fp(r)参数线；不是目标半径的多项式或特征零次数证书。"}
    stage = "carrier"
    try:
        A, values = carrier(args.case, args.prime, args.seconds)
        stage = "closure_and_total_derivatives"
        registers = execute(small, values)
        actual = [registers[i] for i in outputs]
        report["circuit_seconds"] = time.perf_counter()-tick
        stage = "complete_norm"
        if A is None:
            N = actual[0]
            norm_gradients = actual[1:]
        else:
            N = A.norm(actual[0].promote(4).coefficients, 4)
            report["carrier_rank"] = 16
            report["closure_vector"] = [q.encode() for q in actual[0].promote(4).coefficients]
            report["total_derivative_vectors"] = [[q.encode() for q in value.promote(4).coefficients] for value in actual[1:]]
            report["relations"] = [{"a": [q.encode() for q in a], "b": [q.encode() for q in b]} for a, b in A.relations]
            stage = "norm_total_derivatives"
            inverse = actual[0].inverse()
            norm_gradients = [N*trace(A, (inverse*value).promote(4).coefficients, 4) for value in actual[1:]]
            if args.case == "primitive16":
                # 因子只依赖四层载体，范数恒等式可直接与旧闭环比较。
                factor_data, factor_outputs = prune(data, [data["outputs"]["removed_factor"]])
                factor = execute(factor_data, values)[factor_outputs[0]].promote(4)
                report["removed_factor_norm"] = A.norm(factor.coefficients, 4).encode()
            report["recursive_products"], report["recursive_inversions"] = A.products, A.inversions
        report.update(status="completed", norm=N.encode(), norm_total_derivatives=[q.encode() for q in norm_gradients],
                      numerator_degree=int(N.n.degree()), denominator_degree=int(N.d.degree()))
    except Exception as error:
        from tower import TowerStop
        if not isinstance(error, (TowerStop, ZeroDivisionError)):
            raise
        report.update(status="stopped", stopped_stage=stage, reason=str(error))
    from tower import Rational
    report.update(seconds=time.perf_counter()-tick, metrics=Rational.metrics)
    write(path, report)
    print({key: report[key] for key in ("case", "status", "seconds")}, flush=True)


def controller(args):
    prefix = RESULTS / args.tag
    if prefix.with_suffix(".json").exists():
        raise FileExistsError("拒绝覆盖比较记录")
    records = []
    with bounded() as limits:
        for case in CASES:
            log = RESULTS / (args.tag+"-"+case+"-p"+str(args.prime)+".log")
            command = [sys.executable, "-B", __file__, "--worker", "--case", case, "--tag", args.tag,
                       "--prime", str(args.prime), "--seconds", str(args.seconds)]
            tick = time.perf_counter()
            with log.open("xb") as stream:
                try:
                    child = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=args.seconds+5, check=False)
                    status = "returned"
                    code = child.returncode
                except subprocess.TimeoutExpired:
                    status, code = "hard_timeout", None
            output = RESULTS / (args.tag+"-"+case+"-p"+str(args.prime)+".json.gz")
            record = {"case": case, "controller_status": status, "exit_code": code,
                      "wall_seconds": time.perf_counter()-tick, "log_sha256": sha(log)}
            if output.exists():
                result = read(output)
                record.update(status=result["status"], result_sha256=sha(output))
            records.append(record)
            print(record, flush=True)
        report = {"cases": records, "prime": args.prime, "per_case_seconds": args.seconds,
                  "limits": limits(), "driver_sha256": sha(__file__),
                  "minimal_polynomial_obtained": False, "timing_is_single_sample": True}
        write(prefix.with_suffix(".json"), report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--case", choices=list(CASES))
    parser.add_argument("--tag", required=True)
    parser.add_argument("--prime", type=int, choices=[65521, 65537], default=65521)
    parser.add_argument("--seconds", type=int, choices=range(5, 121), default=45)
    arguments = parser.parse_args()
    (worker if arguments.worker else controller)(arguments)
