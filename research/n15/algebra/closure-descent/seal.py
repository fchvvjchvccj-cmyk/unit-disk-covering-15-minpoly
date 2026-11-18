"""封存四参数投影的源码、活动证书、诊断记录及父系统指纹。"""

import argparse
from fractions import Fraction
import platform
import sys

from descent_support import HERE, SHARED, ROOT, SOURCE, read, write, sha
from verify_projection import content_sha


NORM_RUNS = ["20261015-p65521-final", "20261015-p65519-final", "20261015-second-line-final"]


def fingerprint(path):
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}


def parents():
    path = SHARED / "provenance.json"
    manifest = read(path)
    return sorted({path}|{ROOT / row["path"] for row in manifest["files"]+manifest["parent_inputs"]})


def verify(manifest):
    paths = set()
    for row in manifest["files"]+manifest["parent_inputs"]:
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != row["bytes"] or sha(path) != row["sha256"]:
            raise ValueError("四参数投影文件指纹改变："+row["path"])
        if row["path"] in paths:
            raise ValueError("封存清单存在重复路径")
        paths.add(row["path"])
    for path in HERE.glob("*.py"):
        compile(path.read_text(encoding="utf-8"), str(path), "exec")
    source = HERE / "system-runs/20261015-projected4-final"
    circuit = read(source / "circuit.json.gz")
    if (circuit["variables"] != ["u", "v", "w", "r"] or circuit["keep_left_roots"] != 0
            or circuit["generator_sha256"] != sha(HERE / "build_projected.py")
            or circuit["coefficient_formula_sha256"] != sha(HERE / "build_descent.py")
            or circuit["pair_arithmetic_sha256"] != sha(HERE / "quadratic_pair.py")
            or circuit["source_circuit_sha256"] != sha(SOURCE)):
        raise ValueError("四参数投影构造与源代码未绑定")
    proof = read(HERE / "relative-identities-accepted.json.gz")
    if (not proof["verified"] or proof["characteristic_zero_universal_identities"] != 5
            or proof["verifier_sha256"] != sha(HERE / "prove_relative.py")):
        raise ValueError("通用特征零恒等式未通过接受端")
    certificate, target = read(source / "target-final-certificate.json.gz"), read(source / "target-final-verified.json.gz")
    if (certificate["verifier_sha256"] != sha(HERE / "verify_projection.py")
            or certificate["circuit_sha256"] != sha(source / "circuit.json.gz")
            or certificate["circuit_content_sha256"] != content_sha(circuit)
            or not target["verified"] or not target["jacobian_invertible"]
            or target["jacobian_dimension"] != 4 or target["target_local_length"] != 1
            or Fraction(target["inverse_defect_upper"]) >= 1):
        raise ValueError("四维目标单根证书没有保持绑定")
    for name in NORM_RUNS:
        folder = HERE / "norm-runs" / name
        accepted, functions = read(folder / "verified.json.gz"), read(folder / "functions.json.gz")
        if (not accepted["verified"] or accepted["symbolic_total_derivatives"] != 4
                or accepted["circuit_sha256"] != sha(source / "circuit.json.gz")
                or accepted["functions_sha256"] != sha(folder / "functions.json.gz")
                or accepted["driver_sha256"] != sha(HERE / "replay_norm.py")
                or len(functions["G_and_gradient_u_v_w_r"]) != 5):
            raise ValueError("符号范数与四个导数没有完整回放："+name)
    geometry = read(HERE / "geometry-final-verified.json.gz")
    if (not geometry["verified"] or geometry["verifier_sha256"] != sha(HERE / "verify_geometry.py")
            or geometry["circuit_sha256"] != sha(source / "circuit.json.gz")
            or sum(row["original_geometry_points"] for row in geometry["generic_fibers"]) != 128
            or len(geometry["closed_recovery_tests"]) != 2
            or any(not row["verified"] for row in geometry["closed_recovery_tests"])):
        raise ValueError("原几何及闭环恢复回放不完整")
    checks = read(HERE / "checks-final-verified.json.gz")
    if (not checks["verified"] or checks["positive_replays"] != 5 or len(checks["negative_controls_rejected"]) != 12
            or checks["driver_sha256"] != sha(HERE / "checks.py")
            or checks["circuit_sha256"] != sha(source / "circuit.json.gz")
            or checks["certificate_sha256"] != sha(source / "target-final-certificate.json.gz")):
        raise ValueError("正例与错误材料拒绝检查未通过")
    experiment = read(HERE / "expansion-runs/20261015-four-parameter/controller.json")
    experiment_source = HERE / "system-runs/20261015-projected4/circuit.json.gz"
    experiment_circuit = read(experiment_source)
    if (experiment["worker"]["driver_sha256"] != sha(HERE / "expand_projection.py")
            or experiment["worker"]["circuit_sha256"] != sha(experiment_source)
            or any(experiment_circuit[key] != circuit[key] for key in ("variables", "nodes", "outputs", "all_divisor_nodes"))
            or experiment["worker"]["status"] != "stopped"):
        raise ValueError("平坦展开的停止诊断与活动电路不符")
    return {"verified": True, "artifact_files": len(manifest["files"]), "parent_inputs": len(manifest["parent_inputs"]),
            "python_sources": len(list(HERE.glob("*.py"))), "active_variables": circuit["variables"],
            "target_local_length": 1, "symbolic_slices_with_four_derivatives": 3,
            "closed_recovery_tests": 2, "radius_minimal_polynomial_obtained": False}


def run(args):
    path = HERE / "provenance.json"
    if args.create:
        import flint
        if path.exists():
            raise FileExistsError("拒绝覆盖已有四参数投影封存")
        files = [entry for entry in sorted(HERE.rglob("*")) if entry.is_file() and entry != path]
        manifest = {"schema": "n15-closure-descent-provenance-v1", "date": "2026-10-15",
                    "files": [fingerprint(entry) for entry in files], "parent_inputs": [fingerprint(entry) for entry in parents()],
                    "runtime": {"python": sys.version, "python_executable": sys.executable,
                                "flint": flint.__version__, "platform": platform.platform()},
                    "active_system": "20261015-projected4-final", "active_target_certificate": "target-final-certificate.json.gz",
                    "active_target_report": "target-final-verified.json.gz", "active_norm_runs": NORM_RUNS,
                    "historical_notes": "保留原型、中间投影、初次接受和初始封存清单。三处源码尾部空行规范化后重新生成活动四变量电路及全部接受报告；旧平坦展开电路的数学节点及输出与活动电路完全相同，因项数上限停止。",
                    "target_local_length": 1, "global_zero_dimensionality_proved": False,
                    "radius_minimal_polynomial_obtained": False}
        result = verify(manifest)
        write(path, manifest)
    else:
        result = verify(read(path))
    print(result, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--create", action="store_true")
    run(parser.parse_args())
