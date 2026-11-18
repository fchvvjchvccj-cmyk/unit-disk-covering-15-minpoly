"""只无损压缩本目录的格原始产物；删除原件前逐字节检查解压结果。"""

import gzip
import hashlib
import json
import time

from refine import HERE, sha, write_new
from windows_compute_limits import apply_limits


def main():
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        records = []
        selected = list((HERE/"lattices").glob("*.input.txt"))
        selected += [path for path in (HERE/"lattices").glob("*.json")
                     if not path.name.endswith((".run.json", ".verified.json"))]
        selected += list((HERE/"controls").glob("*.input.txt"))
        selected += list((HERE/"controls").glob("*.lattice.json"))
        for path in sorted(selected):
            if path.stat().st_size == 0:
                continue
            if not path.resolve().is_relative_to(HERE.resolve()):
                raise ValueError("压缩对象越出本任务目录")
            original = path.read_bytes()
            compressed = gzip.compress(original, compresslevel=9, mtime=0)
            target = path.with_name(path.name+".gz")
            target.open("xb").write(compressed)
            if gzip.decompress(target.read_bytes()) != original:
                raise ValueError("压缩回读逐字节不一致")
            records.append({"original_path": path.relative_to(HERE).as_posix(),
                            "original_bytes": len(original),
                            "original_sha256": hashlib.sha256(original).hexdigest(),
                            "archive_path": target.relative_to(HERE).as_posix(),
                            "archive_bytes": len(compressed), "archive_sha256": sha(target),
                            "decompressed_bytes_identical": True})
            path.unlink()
        write_new(HERE / "artifact-compression.json",
                  {"verified": True, "records": records,
                   "original_bytes": sum(v["original_bytes"] for v in records),
                   "archive_bytes": sum(v["archive_bytes"] for v in records),
                   "seconds": time.perf_counter()-started, "windows_job_limits": limits()})
        print(json.dumps({"compressed_files": len(records), "seconds": time.perf_counter()-started}), flush=True)
    finally:
        close()


if __name__ == "__main__":
    main()
