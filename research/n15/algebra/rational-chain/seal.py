"""封存本轮输入与产物指纹，并提供只读的完整性及语法检查。"""

import argparse
import platform
import sys

from support import HERE, ROOT, TARGET, AXIS, FACTOR, read, write, sha


def fingerprint(path):
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}


def parent_paths():
    from target_chart import DEPENDENCIES
    candidates = set([AXIS, FACTOR, TARGET / "refined-root.json", TARGET / "intervals.py",
                       HERE.parent / "chain-radius/common.py", HERE.parent / "chain-radius/chain_fiber.py",
                       HERE.parent / "chain-radius/verify_fiber.py", ROOT / "scripts/windows_compute_limits.py",
                      ROOT / "scripts/compact_computation_storage.py"]+list(DEPENDENCIES.values()))
    return sorted(path for path in candidates if HERE not in path.parents)


def verify(manifest):
    paths = set()
    for entry in manifest["files"]+manifest["parent_inputs"]:
        path = ROOT / entry["path"]
        if not path.is_file() or path.stat().st_size != entry["bytes"] or sha(path) != entry["sha256"]:
            raise ValueError("文件指纹不符："+entry["path"])
        if entry["path"] in paths:
            raise ValueError("清单存在重复路径")
        paths.add(entry["path"])
    for path in sorted(HERE.glob("*.py")):
        compile(path.read_text(encoding="utf-8"), str(path), "exec")
    checks = read(HERE / "checks-final-verified.json")
    if not checks["verified"] or len(checks["negative_controls_rejected"]) != 11 or checks["positive_replays"] != 7:
        raise ValueError("本轮检查结果不完整")
    if checks["driver_sha256"] != sha(HERE / "checks.py"):
        raise ValueError("检查脚本与接受结果绑定改变")
    for name in ("20261015-critical9-final", "20261015-critical10-final"):
        source = HERE / "system-runs" / name
        cert = read(source / "20261015-target-bound/certificate.json.gz")
        result = read(source / "20261015-target-bound/verified.json.gz")
        if cert["verifier_sha256"] != sha(HERE / "target_chart.py") or cert["circuit_sha256"] != sha(source / "circuit.json.gz"):
            raise ValueError("目标证书与当前验证器或电路不一致")
        if not result["verified"] or not result["jacobian_invertible"] or result["target_local_length"] != 1:
            raise ValueError("目标单根证书未通过")
    return {"verified": True, "artifact_files": len(manifest["files"]),
            "parent_inputs": len(manifest["parent_inputs"]), "python_sources": len(list(HERE.glob("*.py"))),
            "radius_minimal_polynomial_obtained": False}


def run(args):
    target = HERE / "provenance.json"
    if args.create:
        import flint
        if target.exists():
            raise FileExistsError("拒绝覆盖已有封存清单")
        files = [path for path in sorted(HERE.rglob("*")) if path.is_file() and path != target]
        manifest = {"schema": "n15-rational-chain-provenance-v1", "date": "2026-10-15",
                    "files": [fingerprint(path) for path in files], "parent_inputs": [fingerprint(path) for path in parent_paths()],
                    "runtime": {"python": sys.version, "python_executable": sys.executable, "flint": flint.__version__,
                                "platform": platform.platform()},
                    "active_systems": ["20261015-critical9-final", "20261015-critical10-final"],
                    "active_target_certificate_tag": "20261015-target-bound",
                    "generic_carrier_rank_before_closure": 64, "target_local_length": 1,
                    "global_zero_dimensionality_proved": False, "radius_minimal_polynomial_obtained": False,
                    "history_note": "包括原型、旧目标报告及失败展开试验的原始字节；活动接受端只使用明确指定的系统和目标证书。"}
        report = verify(manifest)
        write(target, manifest)
    else:
        report = verify(read(target))
    print(report, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--create", action="store_true")
    run(parser.parse_args())
