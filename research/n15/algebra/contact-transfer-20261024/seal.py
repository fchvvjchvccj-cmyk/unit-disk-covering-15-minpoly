"""封存本轮源代码、最小输入、证书和正式比较结果，不包含探索性旧输出。"""

from context import ROOT, HERE, INPUTS, RESULTS, sha, write

DEPENDENCIES = [
    "scripts/windows_compute_limits.py", "scripts/compact_computation_storage.py",
    "n15/algebra/chain-radius/common.py", "n15/algebra/chain-radius/chain_fiber.py",
    "n15/algebra/rational-chain/circuit.py", "n15/algebra/rational-chain/geometry.py",
    "n15/algebra/rational-chain/prove_branch.py", "n15/algebra/rational-chain/support.py",
    "n15/algebra/shared-radical/tower.py", "n15/algebra/shared-radical/prove_shared.py",
    "n15/algebra/shared-radical/verify_quotient.py", "n15/algebra/shared-radical/bridge_support.py",
    "n15/algebra/target-reconstruction/intervals.py",
    "references/diskcoveringsolve-089584973f41e6374ce110f6e5b5ee456df19352/cover15/proof_bundle/core/system15.py",
    "docs/n15/algebra/contact-transfer-20261024/proof.zh.md",
]


def main():
    paths = list(HERE.glob("*.py"))+list(INPUTS.glob("*"))
    paths += [RESULTS / name for name in ("primitive16.json.gz", "target-final-certificate.json.gz",
                                       "verified-final.json", "comparison-final.json", "checks.json")]
    for tag, prime in (("exact65521", 65521), ("exact65537", 65537)):
        paths.append(RESULTS / (tag+".json"))
        for case in ("primitive16", "descent16", "primitive4"):
            paths += [RESULTS / (tag+"-"+case+"-p"+str(prime)+suffix) for suffix in (".json.gz", ".log")]
    paths += [ROOT / path for path in DEPENDENCIES]
    records = [{"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}
               for path in sorted(set(paths))]
    write(HERE / "manifest.json", {"files": records, "scope": "候选局部结构与参数线比较；未取得半径最小多项式或重放全局证明。"})
    print({"sealed_files": len(records)})


if __name__ == "__main__":
    main()
