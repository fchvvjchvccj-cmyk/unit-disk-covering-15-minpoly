"""对证明接受端做损坏拒绝检查，并核对区间舍入在正负边界的包含性。"""

import copy
from fractions import Fraction as F
import hashlib
import json
import sys

from intervals import Dyadic
from refine import HERE, artifact_bytes, write_new
import verify_lattices as lattice
from verify_root import verify as verify_root
from windows_compute_limits import apply_limits


def main():
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    try:
        arithmetic = Dyadic(16)
        intervals = [(F(-9, 7), F(-8, 7)), (F(-1, 11), F(1, 13)),
                     (F(2, 17), F(3, 17)), (F(3, 2), F(7, 3))]
        checked = 0
        for a in intervals:
            ai = (arithmetic.rational(a[0])[0], arithmetic.rational(a[1])[1])
            for b in intervals:
                bi = (arithmetic.rational(b[0])[0], arithmetic.rational(b[1])[1])
                product = arithmetic.mul(ai, bi)
                for x in a:
                    for y in b:
                        if not F(product[0], arithmetic.scale) <= x*y <= F(product[1], arithmetic.scale):
                            raise ValueError("区间乘法不包含精确端点积")
                        checked += 1
                if not b[0] <= 0 <= b[1]:
                    quotient = arithmetic.div(ai, bi)
                    for x in a:
                        for y in b:
                            if not F(quotient[0], arithmetic.scale) <= x/y <= F(quotient[1], arithmetic.scale):
                                raise ValueError("区间除法不包含精确端点商")
                            checked += 1
        root = json.loads((HERE / "refined-root.json").read_text(encoding="utf-8"))
        config = json.loads((HERE / "search-budget.json").read_text(encoding="utf-8"))
        case = json.loads((HERE / "lattices/r-d16.run.json").read_text(encoding="utf-8"))
        basis_bytes = artifact_bytes(HERE / "lattices/r-d16.json")
        input_bytes = artifact_bytes(HERE / "lattices/r-d16.input.txt")
        proposal = json.loads(basis_bytes)
        originals = lattice.artifact_bytes, lattice.artifact_sha
        rejected = []
        try:
            for label in ("last_column_changed", "integer_inverse_changed", "power_error_outside_bound"):
                damaged = copy.deepcopy(proposal)
                raw_input = input_bytes
                if label == "last_column_changed":
                    damaged["basis"][0][-1] = str(int(damaged["basis"][0][-1])+1)
                elif label == "integer_inverse_changed":
                    damaged["inverse_transform"][0][0] = str(int(damaged["inverse_transform"][0][0])+1)
                else:
                    entries = list(map(int, input_bytes.decode().split()))
                    size = entries[0]
                    entries[2] += 100
                    raw_input = ("\n".join(map(str, entries))+"\n").encode()
                    damaged["basis"] = [[str(int(i == j)) for j in range(size)]+[str(entries[1+i])]
                                        for i in range(size)]
                    damaged["inverse_transform"] = [[str(int(i == j)) for j in range(size)] for i in range(size)]
                raw_proposal = json.dumps(damaged).encode()
                fake = {"r-d16.input.txt": raw_input, "r-d16.json": raw_proposal}
                lattice.artifact_bytes = lambda path: fake[path.name]
                lattice.artifact_sha = lambda path: hashlib.sha256(fake[path.name]).hexdigest()
                test_case = {**case, "input_sha256": hashlib.sha256(raw_input).hexdigest(),
                             "data_sha256": hashlib.sha256(raw_proposal).hexdigest()}
                try:
                    lattice.verify_case(test_case, root, config)
                except ValueError as error:
                    rejected.append({"name": label, "rejected": True, "reason": str(error)})
                else:
                    raise AssertionError("损坏格证书未被拒绝")
        finally:
            lattice.artifact_bytes, lattice.artifact_sha = originals
        damaged_root = copy.deepcopy(root)
        damaged_root["xnum"][0] = str(int(damaged_root["xnum"][0])+10**5500)
        try:
            verify_root(damaged_root)
        except ValueError as error:
            rejected.append({"name": "root_outside_parent_box", "rejected": True, "reason": str(error)})
        else:
            raise AssertionError("错误根盒未被拒绝")
        report = {"verified": True, "exact_endpoint_checks": checked,
                  "damaged_certificates": rejected, "windows_job_limits": limits()}
        write_new(HERE / "acceptor-rejection-tests.json", report)
        print(json.dumps(report, ensure_ascii=False), flush=True)
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
