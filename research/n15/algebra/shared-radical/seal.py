"""封存共享根实验并只读核验活动证书、源代码及父输入的完整性。"""

import argparse
from fractions import Fraction
import platform
import sys

from bridge_support import HERE, PARENT, ROOT, read, write, sha


LINES = ["20261015-jets-p65521", "20261015-jets-p65519", "20261015-jets-second-line"]


def fingerprint(path):
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}


def parents():
    path = PARENT / "provenance.json"
    sealed = read(path)
    return sorted({path}|{ROOT / row["path"] for row in sealed["files"]+sealed["parent_inputs"]})


def verify(manifest):
    paths = set()
    for entry in manifest["files"]+manifest["parent_inputs"]:
        path = ROOT / entry["path"]
        if not path.is_file() or path.stat().st_size != entry["bytes"] or sha(path) != entry["sha256"]:
            raise ValueError("共享根封存指纹改变："+entry["path"])
        if entry["path"] in paths:
            raise ValueError("封存清单中有重复路径")
        paths.add(entry["path"])
    for path in HERE.glob("*.py"):
        compile(path.read_text(encoding="utf-8"), str(path), "exec")
    proof = read(HERE / "proof-runs/20261015-target/verified.json.gz")
    if not proof["verified"] or proof["verifier_sha256"] != sha(HERE / "prove_shared.py"):
        raise ValueError("共享根选择证明与源码不符")
    source = HERE / "system-runs/20261015-selected"
    circuit = read(source / "circuit.json.gz")
    certificate = read(source / "final-target-certificate.json.gz")
    target = read(source / "final-verified.json.gz")
    if circuit["generator_sha256"] != sha(HERE / "build_shared.py") or circuit["generic_rank_before_closure"] != 32:
        raise ValueError("九变量电路构造绑定改变")
    if certificate["verifier_sha256"] != sha(HERE / "verify_shared.py") or certificate["circuit_sha256"] != sha(source / "circuit.json.gz"):
        raise ValueError("目标证书与电路或接受端失去绑定")
    if (not target["verified"] or not target["jacobian_invertible"] or target["target_local_length"] != 1
            or target["jacobian_dimension"] != 9 or Fraction(target["inverse_defect_upper"]) >= 1):
        raise ValueError("九维目标单根证书未通过")
    for name in LINES:
        folder = HERE / "quotient-runs" / name
        data, accepted = read(folder / "quotient.json.gz"), read(folder / "verified-final.json.gz")
        if (data["driver_sha256"] != sha(HERE / "quotient_probe.py") or data["algebra_sha256"] != sha(HERE / "tower.py")
                or not accepted["verified"] or accepted["quotient_sha256"] != sha(folder / "quotient.json.gz")
                or accepted["norm_sha256"] != sha(folder / "norm.json.gz")
                or accepted["circuit_sha256"] != sha(source / "circuit.json.gz")
                or accepted["verifier_sha256"] != sha(HERE / "verify_quotient.py")
                or accepted["symbolic"]["symbolic_total_derivatives"] != 4
                or sum(row["verified_points"] for row in accepted["samples"]["samples"]) != 256):
            raise ValueError("单参数规范形没有完整绑定接受结果："+name)
        if name in LINES[:2]:
            analysis = read(folder / "factor-analysis-final.json.gz")
            if (not analysis["verified"] or not analysis["all_factors_visited"]
                    or analysis["verifier_sha256"] != sha(HERE / "analyze_factors.py")
                    or analysis["quotient_sha256"] != sha(folder / "quotient.json.gz")
                    or analysis["norm_sha256"] != sha(folder / "norm.json.gz")):
                raise ValueError("剩余域分析未通过输入和源码绑定")
            classified = [row for row in analysis["reports"] if row["status"] == "rank_computed"]
            if any(row["saturation_power"] != 32 for row in classified):
                raise ValueError("剩余域分析未采用稳定饱和")
    folder = HERE / "multivariate-runs/20261015-two-parameter-final"
    data, accepted = read(folder / "quotient.json.gz"), read(folder / "verified.json.gz")
    if (not accepted["verified"] or accepted["quotient_sha256"] != sha(folder / "quotient.json.gz")
            or accepted["circuit_sha256"] != sha(source / "circuit.json.gz")
            or accepted["verifier_sha256"] != sha(HERE / "verify_multivariate.py")
            or data["driver_sha256"] != sha(HERE / "multivariate_probe.py")
            or data["model_sha256"] != sha(HERE / "quotient_probe.py") or data["algebra_sha256"] != sha(HERE / "tower.py")):
        raise ValueError("双参数规范形未通过绑定接受端")
    controls = read(HERE / "checks-verified.json")
    if (not controls["verified"] or controls["positive_replays"] != 6 or len(controls["negative_controls_rejected"]) != 15
            or controls["driver_sha256"] != sha(HERE / "checks.py")
            or controls["finite_algebra_controls"]["nilpotent_stable_rank"] != 0
            or controls["finite_algebra_controls"]["split_joint_zero_dimension"] != 0):
        raise ValueError("错误材料或有限代数反例检查未完整通过")
    return {"verified": True, "artifact_files": len(manifest["files"]), "parent_inputs": len(manifest["parent_inputs"]),
            "python_sources": len(list(HERE.glob("*.py"))), "carrier_rank_before_closure": 32,
            "target_local_length": 1, "symbolic_one_parameter_slices": 3, "symbolic_two_parameter_closures": 1,
            "radius_minimal_polynomial_obtained": False}


def run(args):
    path = HERE / "provenance.json"
    if args.create:
        import flint
        if path.exists():
            raise FileExistsError("拒绝覆盖已有共享根封存")
        files = [entry for entry in sorted(HERE.rglob("*")) if entry.is_file() and entry != path]
        manifest = {"schema": "n15-shared-radical-provenance-v1", "date": "2026-10-15",
                    "files": [fingerprint(entry) for entry in files], "parent_inputs": [fingerprint(entry) for entry in parents()],
                    "runtime": {"python": sys.version, "python_executable": sys.executable,
                                "flint": flint.__version__, "platform": platform.platform()},
                    "active_one_parameter_runs": LINES, "active_two_parameter_run": "20261015-two-parameter-final",
                    "active_target_certificate": "system-runs/20261015-selected/final-target-certificate.json.gz",
                    "historical_notes": "保留旧接受报告、单次像分析、限额停止和首次序列化失败的日志及未完成文件；活动路径在本清单中明确列出。",
                    "carrier_rank_before_closure": 32, "target_local_length": 1,
                    "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False}
        result = verify(manifest)
        write(path, manifest)
    else:
        result = verify(read(path))
    print(result, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--create", action="store_true")
    run(parser.parse_args())
