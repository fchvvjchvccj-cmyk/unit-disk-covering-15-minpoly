"""本轮接触链研究的固定输入、资源配额和只追加产物接口。"""

from contextlib import contextmanager
import gzip
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUTS = HERE / "inputs"
RESULTS = HERE / "results"
for directory in ("native/direct-rur-deps/runtime", "scripts", "n15/algebra/chain-radius",
                  "n15/algebra/rational-chain", "n15/algebra/shared-radical",
                  "n15/algebra/closure-descent", "n15/algebra/target-reconstruction"):
    sys.path.append(str(ROOT / directory))
sys.set_int_max_str_digits(0)


def payload(path):
    path = Path(path)
    raw = path.read_bytes()
    return gzip.decompress(raw) if path.suffix == ".gz" else raw


def read(path):
    return json.loads(payload(path).decode("utf-8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))+"\n").encode("utf-8")
    if path.suffix == ".gz":
        raw = gzip.compress(raw, mtime=0)
    with path.open("xb") as stream:
        stream.write(raw)


@contextmanager
def bounded():
    from compact_computation_storage import lease, KERNEL
    from windows_compute_limits import apply_limits
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


def prune(data, outputs):
    needed, todo = set(), list(outputs)
    while todo:
        index = todo.pop()
        if index not in needed:
            needed.add(index)
            todo.extend(data["nodes"][index].get("args", []))
    order = sorted(needed)
    remap = {old: new for new, old in enumerate(order)}
    nodes = []
    for old in order:
        node = dict(data["nodes"][old])
        if "args" in node:
            node["args"] = [remap[index] for index in node["args"]]
        nodes.append(node)
    return {"variables": data["variables"], "nodes": nodes}, [remap[i] for i in outputs]


def execute(data, values):
    registers = []
    for i, node in enumerate(data["nodes"]):
        indices = node.get("args", [])
        if any(type(j) is not int or not 0 <= j < i for j in indices):
            raise ValueError("电路不是顺序无环依赖")
        args = [registers[j] for j in indices]
        op = node["op"]
        if op == "const":
            value = values[0]*0+node["value"]
        elif op == "var":
            value = values[node["value"]]
        elif op == "neg":
            value = -args[0]
        elif op == "add":
            value = args[0]+args[1]
        elif op == "mul":
            value = args[0]*args[1]
        elif op == "div":
            value = args[0]/args[1]
        else:
            raise ValueError("未知电路指令")
        registers.append(value)
    return registers
