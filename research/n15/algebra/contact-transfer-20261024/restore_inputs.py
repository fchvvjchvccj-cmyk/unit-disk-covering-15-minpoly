"""从固定 Git 提交按原始字节恢复本轮必要输入，不复原旧实验目录。"""

import gzip
import hashlib
import subprocess

from context import ROOT, INPUTS, write

SOURCE_REVISION = "ecadc60"
SOURCES = {
    "axis25-system.json.gz": "n15/algebra/minpoly-final/axis25-system.json",
    "axis25-root.json.gz": "n15/algebra/minpoly-final/axis25-root.json",
    "shared32.json.gz": "n15/algebra/shared-radical/system-runs/20261015-selected/circuit.json.gz",
    "descent16.json.gz": "n15/algebra/closure-descent/system-runs/20261015-critical8/circuit.json.gz",
    "primitive4.json.gz": "n15/algebra/primitive-closure/system-runs/20261015-linear4-final/circuit.json.gz",
}


def restore():
    revision = subprocess.check_output(["git", "rev-parse", SOURCE_REVISION], cwd=ROOT, text=True).strip()
    records = []
    INPUTS.mkdir(parents=True, exist_ok=True)
    for filename, source in SOURCES.items():
        raw = subprocess.check_output(["git", "show", f"{revision}:{source}"], cwd=ROOT)
        stored = raw if source.endswith(".gz") else gzip.compress(raw, mtime=0)
        target = INPUTS / filename
        if target.exists():
            if target.read_bytes() != stored:
                raise ValueError("已有恢复输入与固定提交不同："+filename)
        else:
            with target.open("xb") as stream:
                stream.write(stored)
        records.append({"path": filename, "git_path": source,
                        "git_blob": subprocess.check_output(["git", "rev-parse", f"{revision}:{source}"], cwd=ROOT, text=True).strip(),
                        "git_bytes": len(raw), "stored_bytes": len(stored),
                        "stored_sha256": hashlib.sha256(stored).hexdigest(),
                        "payload_sha256": hashlib.sha256(gzip.decompress(stored)).hexdigest()})
    return {"revision": revision, "files": records,
            "scope": "恢复数据作为输入；历史成功标志未作为本轮核验证据。"}


if __name__ == "__main__":
    result = restore()
    write(INPUTS / "manifest.json", result)
    print({"restored_files": len(result["files"]), "stored_bytes": sum(x["stored_bytes"] for x in result["files"])})
