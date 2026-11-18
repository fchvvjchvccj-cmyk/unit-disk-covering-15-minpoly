"""省去大逆矩阵，以小整数格基和精确行列式接受器保存同一排除证书。"""

import hashlib
import time

from common import HERE, bounded, payload, read, write, sha


def run():
    from flint import fmpz_mat
    with bounded() as snapshot:
        tick = time.perf_counter()
        records = []
        for summary_path in sorted((HERE / "height-runs").glob("*/summary.json")):
            source = summary_path.parent
            parent = source if (source / "budget.json").exists() else HERE / "height-runs/20261015-wide"
            for case in read(summary_path)["runs"]:
                if case["status"] != "completed":
                    continue
                path = source / (case["name"]+".basis.json") if source == parent else source / case["name"] / "basis.json"
                data = read(path)
                if hashlib.sha256(payload(path)).hexdigest() != case["output_sha256"]:
                    raise ValueError("原始格产物哈希不一致")
                determinant = int(fmpz_mat([[int(c) for c in row[:-1]] for row in data["basis"]]).det())
                if abs(determinant) != 1:
                    raise ValueError("原始格的整数变换不幺模")
                compact = {"schema": "compact-integer-lattice-v1", "basis": data["basis"],
                           "flint_version": data["flint_version"], "integer_transform_determinant": determinant,
                           "original_payload_sha256": case["output_sha256"], "input_sha256": case["input_sha256"],
                           "acceptor_recomputes_determinant": True}
                target = path.with_name(path.stem+".compact.json.gz")
                write(target, compact)
                record = {"path": target.relative_to(HERE).as_posix(), "sha256": sha(target),
                          "original_payload_sha256": case["output_sha256"], "bytes": target.stat().st_size}
                records.append(record)
                print({"path": record["path"], "bytes": record["bytes"]}, flush=True)
        write(HERE / "compact-lattices-manifest.json", {"verified": True, "lattices": records,
                                                       "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                                       "huge_inverse_archives_retained_locally_only": True})


if __name__ == "__main__":
    run()
