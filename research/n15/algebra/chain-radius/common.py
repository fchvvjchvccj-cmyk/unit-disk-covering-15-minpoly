"""接触链实验的固定输入、资源限制和只追加产物接口。"""

from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = HERE.parent / "target-reconstruction"
AXIS = HERE.parent / "minpoly-final/axis25-system.json"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "native/direct-rur-deps/runtime"))
sys.path.insert(0, str(TARGET))
from windows_compute_limits import apply_limits
from compact_computation_storage import lease, KERNEL


def payload(path):
    """归档文件仍按原始字节参与输入哈希与接受端重放。"""
    path = Path(path)
    if not path.exists() and path.suffix != ".gz":
        path = path.with_name(path.name+".gz")
    return gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()


def read(path):
    return json.loads(payload(path).decode("utf-8"))


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("xb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
                stream.write((json.dumps(data, ensure_ascii=False, separators=(",", ":"))+"\n").encode("utf-8"))
    else:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write("\n")


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


@contextmanager
def bounded():
    handle = lease()
    close = None
    try:
        snapshot, close = apply_limits(memory_bytes=8_000_000_000, cpu_percent=50)
        yield snapshot
    finally:
        if close is not None:
            close()
        KERNEL.ReleaseMutex(handle)
        KERNEL.CloseHandle(handle)


def binary_environment():
    env = os.environ.copy()
    env["PATH"] = str(ROOT / ".tools/msolve-paged/chain/ucrt64/bin")+os.pathsep+env["PATH"]
    return env


sys.set_int_max_str_digits(0)
