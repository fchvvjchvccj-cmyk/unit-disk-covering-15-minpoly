"""封存并受限重放本地旧实验；区分找到的脚本与无法确认的截图出处。"""

import json
from pathlib import Path
import subprocess
import sys
import time

from refine import HERE, ROOT, sha, write_new
from windows_compute_limits import apply_limits


def main():
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        output = HERE / "prior-experiment"
        output.mkdir(exist_ok=True)
        sources = []
        for name in ("pslq_alpha.py", "pslq_selftest.py", "lll_relation.py", "wide_search.py"):
            source = ROOT / ".tools" / name
            data = source.read_bytes()
            snapshot = output / name
            snapshot.open("xb").write(data)
            sources.append({"path": str(source), "sha256": sha(source),
                            "snapshot": snapshot.relative_to(HERE).as_posix(),
                            "snapshot_sha256": sha(snapshot), "bytes": len(data)})
        root_path = ROOT / "n15/proofs/global-final/source/core/root15.json"
        data = json.loads(root_path.read_text(encoding="utf-8"))
        source_input = {"path": str(root_path), "sha256": sha(root_path),
                        "Qx": data["Qx"], "xnum76": data["xnum"][76], "rho": data["rho"]}
        # 此处执行封存原件，没有把本轮算法放进旧实验。
        replays = []
        for name, timeout in (("pslq_selftest.py", 15), ("wide_search.py", 45)):
            begin = time.perf_counter()
            log_path = output / (name+".replay.log")
            with log_path.open("xb") as log:
                try:
                    result = subprocess.run([sys.executable, "-B", "-u", str(output/name)],
                                            stdout=log, stderr=subprocess.STDOUT,
                                            cwd=ROOT, timeout=timeout)
                    status, code = "completed", result.returncode
                except subprocess.TimeoutExpired:
                    status, code = "timeout", None
            replays.append({"script": name, "wall_limit_seconds": timeout, "status": status,
                            "exit_code": code, "seconds": time.perf_counter()-begin,
                            "log": log_path.relative_to(HERE).as_posix(), "log_sha256": sha(log_path)})
        report = {"schema": "prior-integer-relation-audit-v1", "sources": sources,
                  "original_input": source_input, "replays": replays,
                  "screenshot_path": "C:/Users/Administrator/Pictures/Screenshots/屏幕截图 2026-10-10 235506.png",
                  "screenshot_sha256": sha(Path("C:/Users/Administrator/Pictures/Screenshots/屏幕截图 2026-10-10 235506.png")),
                  "scripts_are_proven_to_be_screenshot_run": False,
                  "original_frozen_input_basis_transform_and_run_logs_found": False,
                  "static_findings": [
                      "旧源根证书半宽为10^-60；90位中点数字本身不是90位误差证书。",
                      "旧wide_search实际格缩放为10^(88+20+12)=10^120；打印文案10^118不一致。",
                      "lll_dec将整数基转为Decimal；未输出幺模变换或经精确接受的长度下界。",
                      "Gram-Schmidt只重算交换附近行，未重算后续行；未核验最终约化条件。",
                      "旧residual_t对t中点开平方，所评估区间实际围住sqrt(t_mid)而非t；也未包括认证盒误差。",
                      "旧脚本把区间含零标成EXACT ZERO；区间含零不能证明代数值为零。",
                      "pslq_selftest从单位矩阵开始，H全部为1，m=0即在第一轮退出。",
                      "搜索t=r^2失败不能直接推出r的相同次数、相同系数高度下界。"],
                  "screenshot_exclusion_claim_verified": False,
                  "new_experiment_independent": True,
                  "seconds": time.perf_counter()-started, "windows_job_limits": limits()}
        write_new(HERE / "prior-experiment-audit.json", report)
        print(json.dumps({"audit_completed": True, "replays": replays}), flush=True)
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
