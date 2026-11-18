"""核对封存字节并重新执行本轮严格接受检查，拒绝仅复用成功标志。"""

import argparse
from fractions import Fraction as Q
import time

from context import ROOT, HERE, RESULTS, read, write, sha, bounded
from verify import verify_sources, verify_new, universal_identities, transfer_audit, require
from checks import run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="replayed.json")
    args = parser.parse_args()
    with bounded() as limits:
        tick = time.perf_counter()
        manifest = read(HERE / "manifest.json")
        for item in manifest["files"]:
            path = (ROOT / item["path"]).resolve()
            require(path.is_relative_to(ROOT) and path.stat().st_size == item["bytes"] and sha(path) == item["sha256"],
                    "封存文件摘要不匹配："+item["path"])
        _, root, center, source = verify_sources()
        universal_identities()
        transfer_audit(center, Q(root["rho"]))
        result, _ = verify_new(read(RESULTS / "primitive16.json.gz"), center, Q(root["rho"]),
                               read(RESULTS / "target-final-certificate.json.gz"))
        controls = run()
        report = {"verified": True, "files_checked": len(manifest["files"]), "source_root_replayed": source["unique_root_in_box"],
                  "new_root_replayed": result["unique_root_in_box"], "controls": controls,
                  "seconds": time.perf_counter()-tick, "limits": limits(), "manifest_sha256": sha(HERE / "manifest.json"),
                  "radius_minimal_polynomial_obtained": False, "global_covering_proof_replayed": False}
        write(RESULTS / args.output, report)
        print({"verified": True, "files": len(manifest["files"]), "seconds": report["seconds"]})


if __name__ == "__main__":
    main()
