"""从固定 Git 版本取出枚举输入，在独立临时目录重新执行全部数学检查。"""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REVISION = "ecadc602d1f93c74d528d222dae24178ac255e6c"
sys.path.insert(0, str(ROOT / "scripts"))
from windows_compute_limits import apply_limits


def archived_blob(relative):
    key = f"{REVISION}:n15/proofs/global-final/{relative}"
    return subprocess.check_output(["git", "cat-file", "blob", key], cwd=ROOT)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def replay(bundle):
    manifest = json.loads((HERE / "input-manifest.json").read_bytes())
    for item in manifest["files"]:
        destination = bundle / item["path"]
        destination.parent.mkdir(exist_ok=True, parents=True)
        raw = archived_blob("source/" + item["path"])
        if hashlib.sha256(raw).hexdigest() != item.get("source_sha256", item["sha256"]):
            raise ValueError("来源摘要不匹配：" + item["path"])
        destination.write_bytes(raw)

    workspace = bundle / "replay_workspace"
    workspace.mkdir()
    certificate = json.loads((bundle / "enumeration/metric15_certificate.json").read_bytes())
    for family, data in certificate.items():
        name = family + ".txt"
        archive = workspace / (name + ".gz")
        archive.write_bytes(archived_blob("generated-inputs/" + archive.name))
        with gzip.open(archive, "rb") as source, (workspace / name).open("xb") as destination:
            shutil.copyfileobj(source, destination)
        if digest(workspace / name) != data["source_txt_sha256"]:
            raise ValueError("枚举表摘要不匹配：" + family)

    config = json.loads((HERE / "toolchain.json").read_bytes())
    compiler = Path(config["compiler"])
    if digest(compiler) != config["compiler_sha256"]:
        raise ValueError("C++ 编译器摘要不匹配")
    for name, expected in config["headers"].items():
        if digest(ROOT / ".tools/global-final/boost" / name) != expected:
            raise ValueError("Boost 头文件摘要不匹配：" + name)
    env = {**os.environ, "CXX": str(compiler), "CPLUS_INCLUDE_PATH": config["boost_include"],
           "PATH": str(compiler.parent) + os.pathsep + os.environ["PATH"], "PYTHONUTF8": "1"}
    env.pop("PYTHONOPTIMIZE", None)
    subprocess.run([sys.executable, "-S", "-B", "verify_all.py", "--jobs", "1"],
                   cwd=bundle, env=env, check=True)
    report = json.loads((bundle / "reports/MASTER_VERIFIED.json").read_bytes())
    if report["global_optimality_certificate_chain_closed"] is not True or report["unresolved_leaves"] != 0:
        raise ValueError("全局下界未闭合")
    return report


def replay_axis(workspace):
    manifest = json.loads(archived_blob("final-review/inputs-reviewed.json"))
    needed = {"root-certificate.json", "source-manifest.json", "root15.json", "axis-certificate.json",
              "axis25-inputs.json", "axis25-root.json", "axis25-system.json", "elimination-certificate.json",
              "multiplier-free-system.json", "reduced-kkt-system.json", "radius-kkt-system.json",
              "symmetric-kkt-system.json"}
    for name, expected in manifest.items():
        path = workspace / name
        if not (path.name in needed or path.suffix == ".py" or
                path.suffix == ".ms" and path.name.startswith("axis25-")):
            continue
        raw = subprocess.check_output(["git", "cat-file", "blob", f"{REVISION}:{name}"], cwd=ROOT)
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError("根链输入摘要不符：" + name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    axis = workspace / "n15/algebra/minpoly-final"
    subprocess.run([sys.executable, "-S", "-B", "verify_axis.py", "--output", "new-axis.json"],
                   cwd=axis, check=True)
    report = json.loads((axis / "new-axis.json").read_bytes())
    archive = ROOT / "n15/algebra/contact-transfer-20261024/inputs"
    for kind in ("root", "system"):
        raw = gzip.decompress((archive / f"axis25-{kind}.json.gz").read_bytes())
        if report[f"{kind}_sha256"] != hashlib.sha256(raw).hexdigest():
            raise ValueError("轴根未与接触链认证输入绑定")
    bound_root = next(digest(workspace / name) for name in manifest if name.endswith("/root15.json"))
    if bound_root != digest(HERE / "replay/core/root15.json") or not (report["verified"] and report["source_replayed"]):
        raise ValueError("全局根到轴根的连接未闭合")
    report["global_root_sha256"] = bound_root
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--axis-chain", action="store_true", help="重放全局根到轴根的精确映射")
    args = parser.parse_args()
    output = args.output
    if output.exists():
        raise FileExistsError(output)
    snapshot, close = apply_limits()
    try:
        with tempfile.TemporaryDirectory(prefix="cover15-full-", dir=ROOT / ".tools") as temporary:
            report = replay_axis(Path(temporary)) if args.axis_chain else replay(Path(temporary))
        report.update(archive_revision=REVISION, resource_limits=snapshot())
        with output.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    finally:
        close()


if __name__ == "__main__":
    main()
