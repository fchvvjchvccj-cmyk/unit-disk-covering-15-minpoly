"""从压缩的原格和变换证书重放有界关系排除；超时记录不转为结论。"""

import argparse
import hashlib
import time

from common import HERE, TARGET, bounded, payload, read, write, sha
from height_search import verify


def run(args, root_report=None):
    from verify_root import verify as verify_root
    source = HERE / "height-runs" / args.tag
    parent = source if (source / "budget.json").exists() else HERE / "height-runs" / args.parent
    budget = read(parent / "budget.json")
    summary = read(source / "summary.json")
    root_path = TARGET / "refined-root.json"
    if sha(root_path) != budget["root_sha256"]:
        raise ValueError("认证实根与格实验输入不一致")
    with bounded() as snapshot:
        tick = time.perf_counter()
        root = read(root_path)
        if root_report is None:
            root_report = verify_root(root)
        reports = []
        for case in summary["runs"]:
            if case["status"] != "completed":
                reports.append({"name": case["name"], "status": case["status"], "excluded": False})
                continue
            incoming = parent / (case["name"]+".input.txt")
            outgoing = source / (case["name"]+".basis.json") if source == parent else source / case["name"] / "basis.json"
            # 优先核验即将提交的压缩归档，比较的是解压后的原始哈希。
            incoming = incoming.with_name(incoming.name+".gz") if incoming.with_name(incoming.name+".gz").exists() else incoming
            compact = outgoing.with_name(outgoing.stem+".compact.json.gz")
            if compact.exists():
                artifact = read(compact)
                if artifact["original_payload_sha256"] != case["output_sha256"] or artifact["input_sha256"] != case["input_sha256"]:
                    raise ValueError("紧凑格证书与原始运行的输入或输出标识不同")
                outgoing = compact
                paths_to_hash = ((incoming, "input_sha256"),)
            else:
                outgoing = outgoing.with_name(outgoing.name+".gz") if outgoing.with_name(outgoing.name+".gz").exists() else outgoing
                paths_to_hash = ((incoming, "input_sha256"), (outgoing, "output_sha256"))
            for path, key in paths_to_hash:
                if hashlib.sha256(payload(path)).hexdigest() != case[key]:
                    raise ValueError("解压后的格数据与原始哈希不同")
            report = verify(incoming, outgoing, root, case)
            if report["no_nonzero_integer_polynomial_in_budget_proved"] != case["excluded"]:
                raise ValueError("重新接受的结论与原报告不同")
            reports.append({"name": case["name"], "verified": True, "excluded": case["excluded"],
                            "possible_candidates": len(report["candidates_needing_exact_algebraic_validation"])})
            print(reports[-1], flush=True)
        write(source / args.output, {"verified": True, "root_replay": root_report, "runs": reports,
                                    "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                    "verifier_sha256": sha(__file__)})
        return root_report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True, nargs="+")
    parser.add_argument("--parent", default="20261015-wide")
    parser.add_argument("--output", default="replayed.json")
    args = parser.parse_args()
    tags, root_report = args.tag, None
    for tag in tags:
        args.tag = tag
        root_report = run(args, root_report)
