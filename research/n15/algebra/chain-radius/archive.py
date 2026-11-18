"""压缩大整数格的原始字节并登记哈希，工作副本留在本地。"""

import gzip
import hashlib
from pathlib import Path

from common import HERE, write


def run():
    records = []
    for path in sorted((HERE / "height-runs").rglob("*")):
        if not path.is_file() or path.suffix == ".gz":
            continue
        if not (path.name.endswith(".input.txt") or path.name.endswith("basis.json")):
            continue
        raw = path.read_bytes()
        compressed = path.with_name(path.name+".gz")
        if compressed.exists():
            if gzip.decompress(compressed.read_bytes()) != raw:
                raise ValueError("同名归档与原数据不同")
        else:
            with compressed.open("xb") as stream:
                with gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as output:
                    output.write(raw)
        records.append({"original": path.relative_to(HERE).as_posix(),
                        "archive": compressed.relative_to(HERE).as_posix(),
                        "original_sha256": hashlib.sha256(raw).hexdigest(), "original_bytes": len(raw),
                        "archive_bytes": compressed.stat().st_size})
    write(HERE / "archive-manifest.json", {"archives": records,
                                           "raw_working_copies_intentionally_ignored": True})


if __name__ == "__main__":
    run()
