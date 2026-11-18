"""重放认证父根、新根与全部已完成格证书，不重新启动数值搜索。"""

import argparse
import json
import sys
import time

from refine import HERE, ROOT, sha, write_new
from verify_root import verify as verify_refined_root
from verify_lattices import verify_case
from windows_compute_limits import apply_limits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="delivery-replayed.json")
    args = parser.parse_args()
    if (HERE/args.output).exists():
        raise ValueError("重放输出已存在；请指定新的文件名")
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        sys.path.insert(0, str(HERE.parent / "minpoly-final"))
        from verify_axis import verify as verify_parent
        parent = verify_parent()
        print(json.dumps({"parent_replayed": parent["verified"]}), flush=True)
        root = json.loads((HERE / "refined-root.json").read_text(encoding="utf-8"))
        refined = verify_refined_root(root)
        print(json.dumps({"refined_root_replayed": refined["verified"]}), flush=True)
        config = json.loads((HERE / "search-budget.json").read_text(encoding="utf-8"))
        runs = json.loads((HERE / "search-runs.json").read_text(encoding="utf-8"))
        from search import SEEDS
        if config["seeds"] != SEEDS or config["root_sha256"] != sha(HERE / "refined-root.json"):
            raise ValueError("搜索种子或根摘要不符")
        expected = [(v["seed"], v["degree"], v["scale_digits"]) for v in config["cases"]]
        if [(v["seed"], v["degree"], v["scale_digits"]) for v in runs["lll"]] != expected:
            raise ValueError("实际格病例与预算不符")
        reports = []
        for case in runs["lll"]:
            if case["status"] != "completed":
                if case["status"] != "wall_timeout":
                    raise ValueError("存在未说明的搜索故障")
                reports.append({"name": case["name"], "status": case["status"], "exclusion_proved": False})
                continue
            report = verify_case(case, root, config)
            if not report["no_nonzero_integer_polynomial_in_budget_proved"]:
                raise ValueError("本轮声称的有限排除未重放成功")
            if report["candidates_within_height"] or not report["shortest_basis_row"]["exact_interval_excludes_zero"]:
                raise ValueError("候选拒绝结果与本轮交付不同")
            reports.append({"name": case["name"], "exclusion_proved": True})
            print(json.dumps(reports[-1]), flush=True)
        manifest = json.loads((HERE / "SHA256SUMS.json").read_text(encoding="utf-8"))
        for item in manifest["files"]:
            if sha(HERE / item["path"]) != item["sha256"]:
                raise ValueError("交付文件摘要改变："+item["path"])
        write_new(HERE / args.output,
                  {"verified": True, "parent_replayed": True, "refined_root_replayed": True,
                   "cases": reports, "manifest_files_checked": len(manifest["files"]),
                   "manifest_sha256": sha(HERE / "SHA256SUMS.json"),
                   "exact_algebraic_representation_obtained": False,
                   "minimal_polynomial_obtained": False,
                   "seconds": time.perf_counter()-started, "windows_job_limits": limits()})
        print(json.dumps({"delivery_replayed": True, "seconds": time.perf_counter()-started}), flush=True)
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
