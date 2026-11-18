"""核对硬链接占用并无损压缩已结束计算的大文件，保持原路径与全部字节。"""

import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from windows_compute_limits import apply_limits

ROOT = Path(__file__).resolve().parents[1]
KERNEL = ctypes.WinDLL("kernel32", use_last_error=True) if os.name == "nt" else None


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def stored_bytes(path):
    # 对透明压缩文件查询实际存储字节；普通文件返回数据长度，不含簇尾碎片。
    KERNEL.GetCompressedFileSizeW.argtypes = [ctypes.c_wchar_p, ctypes.POINTER(ctypes.c_ulong)]
    KERNEL.GetCompressedFileSizeW.restype = ctypes.c_ulong
    high = ctypes.c_ulong()
    ctypes.set_last_error(0)
    low = KERNEL.GetCompressedFileSizeW(str(path), ctypes.byref(high))
    if low == 0xFFFFFFFF and ctypes.get_last_error():
        raise ctypes.WinError(ctypes.get_last_error())
    return (high.value << 32) | low


def inventory():
    files = {}
    count = logical = 0
    for base, directories, names in os.walk(ROOT, followlinks=False):
        directories[:] = [name for name in directories if not (Path(base) / name).is_symlink()]
        for name in names:
            path = Path(base) / name
            if path.is_symlink():
                continue
            stat = path.stat()
            count += 1
            logical += stat.st_size
            identity = (stat.st_dev, stat.st_ino)
            if identity in files:
                files[identity]["paths"].append(path)
                continue
            files[identity] = {"paths": [path], "bytes": stat.st_size,
                               "stored_bytes": stored_bytes(path), "links": stat.st_nlink}
    summary = {"files": count, "logical_bytes": logical, "unique_files": len(files),
               "unique_logical_bytes": sum(x["bytes"] for x in files.values()),
               "stored_bytes": sum(x["stored_bytes"] for x in files.values())}
    return summary, files


def candidate(path, size):
    if size < 1024 * 1024:
        return False
    relative = path.relative_to(ROOT)
    if relative.parts[:2] == (".tools", "msolve-paged"):
        return path.name.endswith((".f4ck", ".f4ck.prev", ".ck"))
    if (relative.parts[:2] == (".tools", "radius-macaulay")
            and len(relative.parts) > 2 and relative.parts[2].startswith("verify-")):
        return path.suffix == ".bin"
    if relative.parts[:2] == ("n15", "algebra"):
        if path.name.endswith((".f4ck", ".f4ck.prev")):
            return True
        # 新路线的原始矩阵与行操作证明也保留原路径和全部字节。
        return path.suffix == ".bin" and "runs" in relative.parts
    return (relative.parts[:5] == ("n15", "proofs", "global-final", "replay", "replay_workspace")
            and path.suffix in (".txt", ".bin"))


def known_snapshots():
    def snapshots(value):
        # 同时读取旧版直接断点和新版snapshot/previous_snapshot内的断点。
        if isinstance(value, dict):
            if isinstance(value.get("file"), str) and isinstance(value.get("sha256"), str):
                yield value
            for child in value.values():
                yield from snapshots(child)
        elif isinstance(value, list):
            for child in value:
                yield from snapshots(child)

    expected = {}
    for path in (ROOT / "n15").rglob("record.json"):
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("status") == "running":
            raise RuntimeError("存在尚未结束的计算记录：" + str(path))
        for snapshot in snapshots(record):
            source = (ROOT / snapshot["file"]).resolve()
            if not source.is_relative_to(ROOT):
                raise RuntimeError("断点引用超出工作区：" + str(source))
            prior = expected.setdefault(source, snapshot["sha256"])
            if prior != snapshot["sha256"]:
                raise RuntimeError("同一断点路径存在冲突摘要：" + str(source))
    # Macaulay原始矩阵绑定到已提交的行操作证明摘要。
    for path in (ROOT / "n15/algebra/macaulay-radius/runs").glob("*/proof-manifest.json"):
        manifest = json.loads(path.read_text(encoding="utf-8"))
        for name, item in manifest["raw_files"].items():
            source = (path.parent / name).resolve()
            cache = ROOT / ".tools/radius-macaulay" / ("verify-" + path.parent.name) / name
            for target in (source, cache.resolve()):
                if not target.is_relative_to(ROOT):
                    raise RuntimeError("矩阵引用超出工作区：" + str(target))
                prior = expected.setdefault(target, item["sha256"])
                if prior != item["sha256"]:
                    raise RuntimeError("同一矩阵路径存在冲突摘要：" + str(target))
    return expected


def lease():
    # 立即尝试共享租约，防止与求解器同时修改断点；忙时直接退出。
    KERNEL.CreateMutexW.argtypes = [ctypes.c_void_p, ctypes.c_bool, ctypes.c_wchar_p]
    KERNEL.CreateMutexW.restype = ctypes.c_void_p
    KERNEL.WaitForSingleObject.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
    KERNEL.ReleaseMutex.argtypes = [ctypes.c_void_p]
    KERNEL.CloseHandle.argtypes = [ctypes.c_void_p]
    handle = KERNEL.CreateMutexW(None, False, "Local\\CodexSxypComputeLease")
    if not handle:
        raise ctypes.WinError(ctypes.get_last_error())
    if KERNEL.WaitForSingleObject(handle, 0) not in (0, 0x80):
        KERNEL.CloseHandle(handle)
        raise RuntimeError("计算租约正在使用，请在计算结束后运行存储维护")
    return handle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="实际执行透明压缩；默认只统计")
    parser.add_argument("--algorithm", choices=("XPRESS16K", "LZX"), default="XPRESS16K",
                        help="默认采用速度较快的XPRESS16K；LZX进一步压缩但耗时更长")
    parser.add_argument("--output", help="新报告路径，拒绝覆盖已有报告")
    parser.add_argument("--baseline", help="可选的此前占用调查JSON，附入报告供比较")
    parser.add_argument("--only-uncompressed", action="store_true",
                        help="跳过已有压缩或稀疏存储的文件，只处理新增数据")
    args = parser.parse_args()
    if os.name != "nt":
        parser.error("此维护工具需要Windows透明文件压缩")
    if args.apply and not args.output:
        parser.error("实际压缩必须指定新的报告路径")
    output = (ROOT / args.output).resolve() if args.output else None
    if output and (not output.is_relative_to(ROOT) or output.exists()):
        parser.error("报告必须使用工作区内未占用的路径")
    baseline = json.loads(Path(args.baseline).read_text(encoding="utf-8")) if args.baseline else None
    handle = lease()
    close = None
    try:
        limits, close = apply_limits(memory_bytes=8_000_000_000)
        expected = known_snapshots()
        before, files = inventory()
        selected = [(identity, item) for identity, item in files.items()
                    if any(candidate(path, item["bytes"]) for path in item["paths"])
                    and not (args.only_uncompressed and
                             (item["stored_bytes"] < item["bytes"]
                              or getattr(item["paths"][0].stat(), "st_reparse_tag", 0) == 0x80000017
                              or item["paths"][0].stat().st_file_attributes & 0x800))]
        selected.sort(key=lambda entry: -entry[1]["bytes"])
        print(json.dumps({"before": before, "candidate_unique_files": len(selected),
                          "candidate_bytes": sum(x["bytes"] for _, x in selected)}, indent=2), flush=True)
        if not args.apply:
            return
        journal = ROOT / ".tools/storage-cleanup-progress.jsonl"
        journal.parent.mkdir(exist_ok=True)
        rows = []
        start = time.monotonic()
        for index, (identity, item) in enumerate(selected, 1):
            path = item["paths"][0].resolve()
            if not path.is_relative_to(ROOT) or path.is_symlink():
                raise RuntimeError("压缩目标不在工作区或为符号链接：" + str(path))
            original = sha(path)
            for alias in item["paths"]:
                digest = expected.get(alias.resolve())
                if digest and original != digest:
                    raise RuntimeError("既有断点SHA256不符：" + str(alias))
            # 早期试验使用NTFS压缩；转换前先透明解压，路径与硬链接均保持。
            if path.stat().st_file_attributes & 0x800:
                subprocess.run(["compact.exe", "/u", str(path)], check=True, capture_output=True)
            result = subprocess.run(["compact.exe", "/c", "/exe:" + args.algorithm, str(path)],
                                    check=True, capture_output=True)
            if sha(path) != original:
                raise RuntimeError("压缩后的SHA256发生变化：" + str(path))
            for alias in item["paths"]:
                stat = alias.stat()
                if (stat.st_dev, stat.st_ino) != identity or stat.st_size != item["bytes"]:
                    raise RuntimeError("压缩后的文件身份或长度发生变化：" + str(alias))
            row = {"paths": [p.relative_to(ROOT).as_posix() for p in item["paths"]],
                   "bytes": item["bytes"], "before_stored_bytes": item["stored_bytes"],
                   "after_stored_bytes": stored_bytes(path), "sha256": original,
                   "sha256_unchanged": True, "hardlinks_preserved": True,
                   "sealed_digest_checked": any(p.resolve() in expected for p in item["paths"]),
                   "compact_exit_code": result.returncode}
            row["requested_algorithm"] = args.algorithm
            rows.append(row)
            with journal.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(row, ensure_ascii=False) + "\n")
            print(f"[{index}/{len(selected)}] {row['paths'][0]}: "
                  f"{row['bytes']} -> {row['after_stored_bytes']} bytes; SHA256 OK", flush=True)
        after, _ = inventory()
        report = {"schema": "sxyp-storage-compaction-v1", "verified": True,
                  "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                  "method": "Windows transparent compression of completed computation files",
                  "requested_algorithm": args.algorithm,
                  "only_uncompressed": args.only_uncompressed,
                  "before": before, "after": after,
                  "saved_stored_bytes": before["stored_bytes"] - after["stored_bytes"],
                  "accounting": "hardlinks counted once; GetCompressedFileSizeW; cluster slack and filesystem metadata excluded",
                  "files_deleted": 0, "content_changed": 0, "compressed_files": rows,
                  "seconds": time.monotonic() - start, "windows_job_limits": limits(),
                  "maintenance_script_sha256": sha(Path(__file__)),
                  "baseline_observation": baseline}
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print(json.dumps({"verified": True, "after": after,
                          "saved_stored_bytes": report["saved_stored_bytes"],
                          "report": str(output)}, indent=2), flush=True)
    finally:
        if close:
            close()
        KERNEL.ReleaseMutex(handle)
        KERNEL.CloseHandle(handle)


if __name__ == "__main__":
    main()
