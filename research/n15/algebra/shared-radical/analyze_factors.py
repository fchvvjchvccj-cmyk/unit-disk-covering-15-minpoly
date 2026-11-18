"""在范数因子的剩余域中用32维乘法矩阵区分闭环零点与临界零点。"""

import argparse
import time

from bridge_support import HERE, read, write, sha, bounded
from tower import Algebra


class FiniteAlgebra(Algebra):
    def __init__(self, field, relations):
        self.field, self.relations = field, relations
        self.zero, self.one = field(0), field(1)
        self.products = self.inversions = 0


def evaluate(coefficients, value, field):
    result = field(0)
    for c in reversed(coefficients):
        result = result*value+c
    return result


def scalar(data, value, field):
    numerator, denominator = evaluate(data["n"], value, field), evaluate(data["d"], value, field)
    if denominator.is_zero():
        raise ZeroDivisionError("规范形系数在该剩余域有极点")
    return numerator/denominator


def columns(algebra, coefficients):
    result = []
    for j in range(32):
        basis = [algebra.zero]*32
        basis[j] = algebra.one
        result.append(algebra.multiply(tuple(coefficients), tuple(basis), 5))
    return result


def append_span(basis, vectors):
    for vector in vectors:
        row = list(vector)
        for i in range(32):
            if row[i].is_zero():
                continue
            if basis[i] is None:
                inverse = row[i].inverse()
                basis[i] = tuple(value*inverse for value in row)
                break
            factor = row[i]
            for j in range(i, 32):
                row[j] -= factor*basis[i][j]
    return sum(row is not None for row in basis)


def stable_power(algebra, coefficients):
    # 任意不超过32维的商代数中，32次幂已杀去乘法算子的零特征值幂零部分。
    power = tuple(coefficients)
    for _ in range(5):
        power = algebra.multiply(power, power, 5)
    return power


def check_factor(data, factor):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    prime = data["prime"]
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(prime)(factor["coefficients"]), fq_type="FQ_NMOD")
    alpha = field.gen()
    relations = [(tuple(scalar(c, alpha, field) for c in row["a"]),
                  tuple(scalar(c, alpha, field) for c in row["b"])) for row in data["relations"]]
    algebra = FiniteAlgebra(field, relations)
    H = tuple(scalar(c, alpha, field) for c in data["H"])
    gradients = [tuple(scalar(c, alpha, field) for c in row) for row in data["gradient_u_v_w_r"]]
    MH = columns(algebra, H)
    Mr = columns(algebra, gradients[3])
    Mr_stable = columns(algebra, stable_power(algebra, gradients[3]))
    basis_H = [None]*32
    rank_H = append_span(basis_H, MH)
    basis_good_H = list(basis_H)
    rank_H_Hr = append_span(list(basis_good_H), Mr)
    rank_H_stable = append_span(basis_good_H, Mr_stable)
    basis_critical = list(basis_H)
    for gradient in gradients[:3]:
        append_span(basis_critical, columns(algebra, gradient))
    rank_critical = sum(row is not None for row in basis_critical)
    rank_critical_Hr = append_span(list(basis_critical), Mr)
    rank_critical_stable = append_span(list(basis_critical), Mr_stable)
    return {"verified": True, "status": "rank_computed", "rank_H": rank_H, "closed_quotient_dimension": 32-rank_H,
            "closed_raw_Hr_image_dimension": rank_H_Hr-rank_H,
            "closed_dimension_where_Hr_nonzero": rank_H_stable-rank_H,
            "critical_quotient_dimension": 32-rank_critical,
            "critical_raw_Hr_image_dimension": rank_critical_Hr-rank_critical,
            "critical_dimension_where_Hr_nonzero": rank_critical_stable-rank_critical,
            "closed_dimension_where_r_and_Hr_nonzero": 0 if alpha.is_zero() else rank_H_stable-rank_H,
            "critical_dimension_where_r_and_Hr_nonzero": 0 if alpha.is_zero() else rank_critical_stable-rank_critical,
            "saturation_power": 32,
            "interpretation": "乘法矩阵列张成(H,Hu,Hv,Hw)理想，再以Hr^32的稳定像取Hr可逆部分；适用于含幂零元的有限代数。尚未施加其余图册非零条件，未以各范数的公因子替代同一分支条件。"}


def run(args):
    target = HERE / "quotient-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        data, norm = read(target / "quotient.json.gz"), read(target / "norm.json.gz")
        accepted = read(target / "verified-final.json.gz")
        if (not accepted["verified"] or accepted["quotient_sha256"] != sha(target / "quotient.json.gz")
                or accepted["norm_sha256"] != sha(target / "norm.json.gz")
                or accepted["verifier_sha256"] != sha(HERE / "verify_quotient.py")):
            raise ValueError("规范形未通过接受端绑定")
        reports = []
        for i, factor in enumerate(norm["numerator_factors"]):
            started = time.perf_counter()
            if time.perf_counter()-tick > args.seconds:
                break
            try:
                report = check_factor(data, factor)
            except ZeroDivisionError as error:
                report = {"status": "coefficient_pole", "reason": str(error)}
            report.update(index=i, degree=factor["degree"], multiplicity=factor["multiplicity"], seconds=time.perf_counter()-started)
            reports.append(report)
            print({k: report[k] for k in ("index", "degree", "status", "seconds")}, flush=True)
        result = {"verified": True, "all_factors_visited": len(reports) == len(norm["numerator_factors"]),
                  "reports": reports, "seconds": time.perf_counter()-tick, "limits": snapshot(),
                  "quotient_sha256": sha(target / "quotient.json.gz"), "norm_sha256": sha(target / "norm.json.gz"),
                  "verifier_sha256": sha(__file__),
                  "scope": "固定三个参数的单变量截面剩余域分析；系数有极点的因子未分类，不推断全局目标分量或半径次数。"}
        write(target / args.output, result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--seconds", type=int, default=120)
    parser.add_argument("--output", default="factor-analysis-final.json.gz")
    run(parser.parse_args())
