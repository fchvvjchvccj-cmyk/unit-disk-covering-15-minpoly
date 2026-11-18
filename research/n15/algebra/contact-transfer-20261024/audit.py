"""独立核验参数线的范数、全导数和分支公零点，输出下一轮路线依据。"""

import argparse
import time

from context import HERE, INPUTS, RESULTS, read, write, sha, bounded, execute, prune
from compare import CASES, initialize, statistics
from tower import Algebra, Element


def require(condition, message):
    if not condition:
        raise ValueError(message)


def columns(A, coefficients):
    result = []
    for j in range(16):
        vector = [A.zero]*16
        vector[j] = A.one
        result.append(A.multiply(coefficients, tuple(vector), 4))
    return result


def span(basis, vectors):
    for vector in vectors:
        row = list(vector)
        for j in range(16):
            if row[j].is_zero():
                continue
            if basis[j] is None:
                inv = row[j].inverse()
                basis[j] = tuple(x*inv for x in row)
                break
            scale = row[j]
            for k in range(j, 16):
                row[k] -= scale*basis[j][k]
    return sum(row is not None for row in basis)


class FiniteAlgebra(Algebra):
    def __init__(self, field, relations):
        self.field, self.relations = field, relations
        self.zero, self.one = field(0), field(1)
        self.products = self.inversions = 0

    def scalar(self, value):
        return Element(self, 0, (value if isinstance(value, type(self.one)) else self.field(value),))

    def inverse(self, coefficients, level):
        if level == 0 and coefficients[0].is_zero():
            raise ZeroDivisionError("该有限代数元素不是单位")
        return super().inverse(coefficients, level)


def classify(data, factor, prime):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(prime)(list(map(int, factor.coeffs()))), fq_type="FQ_NMOD")
    root = field.gen()

    def evaluate(coeffs):
        value = field(0)
        for coefficient in reversed(coeffs):
            value = value*root+coefficient
        return value

    def decode(q):
        n, d = evaluate(q["n"]), evaluate(q["d"])
        if d.is_zero():
            raise ZeroDivisionError("规范形系数在此因子上有极点")
        return n/d

    report = {"factor": list(map(int, factor.coeffs())), "degree": int(factor.degree())}
    try:
        relations = [(tuple(decode(q) for q in row["a"]), tuple(decode(q) for q in row["b"])) for row in data["relations"]]
        A = FiniteAlgebra(field, relations)
        H = tuple(decode(q) for q in data["closure_vector"])
        gradients = [tuple(decode(q) for q in row) for row in data["total_derivative_vectors"]]
        basis_H = [None]*16
        rank_H = span(basis_H, columns(A, H))
        basis_critical = list(basis_H)
        for gradient in gradients[:3]:
            span(basis_critical, columns(A, gradient))
        rank_critical = sum(x is not None for x in basis_critical)
        report.update(status="classified", closed_dimension=16-rank_H, common_critical_dimension=16-rank_critical)
        if rank_critical == 16:
            report.update(critical_dimension_on_required_open_chart=0, exclusion="共同临界理想已为单位理想")
            return report
        if root.is_zero():
            report.update(critical_dimension_on_required_open_chart=0, exclusion="r=0，不属于目标正半径图表")
            return report
        generators = []
        for j in range(4):
            vector = [A.zero]*2**(j+1)
            vector[2**j] = A.one
            generators.append(Element(A, j+1, tuple(vector)))
        circuit = read(CASES["primitive16"])
        guard_indices = sorted(set(circuit["all_divisor_nodes"]+circuit["simple_root_derivative_nodes"]+
                                   circuit["outputs"]["right_quadratic"]+[
                                       circuit["outputs"]["closure_coefficients"][1], circuit["outputs"]["removed_factor"],
                                       circuit["outputs"]["radius_derivative"]]))
        # 二次式常数项与一次项不要求非零；只加入首项、判别式和必要主元。
        right = circuit["outputs"]["right_quadratic"]
        guard_indices = [i for i in guard_indices if i not in right[1:]]
        needed = sorted(set(guard_indices+right))
        small, outputs = prune(circuit, needed)
        values = generators+[A.scalar(2), A.scalar(3), A.scalar(1), A.scalar(root)]
        registers = execute(small, values)
        mapping = {old: registers[new] for old, new in zip(needed, outputs)}
        guard = A.scalar(1)
        for index in guard_indices:
            guard *= mapping[index]
        aa, bb, cc = [mapping[index] for index in right]
        guard *= bb*bb-4*aa*cc
        stable = (guard**16).promote(4)
        saturated_rank = span(list(basis_critical), columns(A, stable.coefficients))
        report["critical_dimension_on_required_open_chart"] = saturated_rank-rank_critical
        report["saturation_power"] = 16
        report["individually_excluding_guard_nodes"] = [index for index in guard_indices
            if span(list(basis_critical), columns(A, (mapping[index]**16).promote(4).coefficients)) == rank_critical]
    except (ZeroDivisionError, ValueError) as error:
        report.update(status="coefficient_or_chart_pole", reason=str(error))
    return report


def verify_line(tag, prime):
    from flint import nmod_poly
    from verify_quotient import irreducible
    _, R = initialize(prime, 90)
    decode = lambda q: R(nmod_poly(q["n"], prime), nmod_poly(q["d"], prime))
    controller = read(RESULTS / (tag+".json"))
    require(controller["driver_sha256"] == sha(HERE / "compare.py"), "比较驱动版本不匹配")
    rows = {}
    for record in controller["cases"]:
        case = record["case"]
        path = RESULTS / (tag+"-"+case+"-p"+str(prime)+".json.gz")
        require(record["controller_status"] == "returned" and record["exit_code"] == 0, "存在未完成比较")
        require(record["result_sha256"] == sha(path), "比较结果摘要不匹配")
        row = read(path)
        require(row["status"] == "completed" and row["source_sha256"] == sha(CASES[case]), "比较未完成或未绑定输入")
        require(row["statistics"] == statistics(read(CASES[case])), "静态电路统计不匹配")
        rows[case] = row
    new, old, scalar = [rows[k] for k in ("primitive16", "descent16", "primitive4")]
    N, previous, G = [decode(row["norm"]) for row in (new, old, scalar)]
    factor_norm = decode(new["removed_factor_norm"])
    require((N-G).is_zero(), "有限载体范数与四参数电路不同")
    require((previous-factor_norm*factor_norm*N).is_zero(), "本原剥离因子不匹配旧闭环")
    require(all((decode(a)-decode(b)).is_zero() for a, b in zip(new["norm_total_derivatives"], scalar["norm_total_derivatives"])), "四个全参数导数不匹配")
    direct = R(N.n.derivative()*N.d-N.n*N.d.derivative(), N.d*N.d)
    require((direct-decode(new["norm_total_derivatives"][3])).is_zero(), "半径导数不等于直接有理函数微分")
    common = N.n
    for gradient in new["norm_total_derivatives"][:3]:
        common = common.gcd(decode(gradient).n)
    unit, factors = common.factor()
    reconstructed = nmod_poly([int(unit)], prime)
    classifications = []
    for factor, multiplicity in factors:
        require(irreducible(factor, prime), "因子未通过独立不可约判据")
        reconstructed *= factor**multiplicity
        result = classify(new, factor, prime)
        result["multiplicity_in_scalar_common_factor"] = multiplicity
        classifications.append(result)
    require(reconstructed == common, "标量公因子未完整乘回")
    return {"prime": prime, "fixed_u_v_w": [2, 3, 1], "all_cases_completed": True,
            "norms_identical": True, "four_total_derivatives_identical": True,
            "radius_derivative_independently_checked": True, "removed_square_factor_identity": True,
            "old_norm_degrees": [int(previous.n.degree()), int(previous.d.degree())],
            "primitive_norm_degrees": [int(N.n.degree()), int(N.d.degree())],
            "scalar_critical_common_factor_degree": int(common.degree()), "common_factor_classifications": classifications,
            "statistics": {case: row["statistics"] for case, row in rows.items()},
            "seconds_single_sample": {case: row["seconds"] for case, row in rows.items()},
            "metrics": {case: row["metrics"] for case, row in rows.items()},
            "scope": "两个固定端口素域参数线的完整比较；公因子与分支诊断不证明目标半径次数、特征零全局维数或求解加速。"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="comparison-audited.json")
    args = parser.parse_args()
    with bounded() as limits:
        tick = time.perf_counter()
        report = {"verified": True, "lines": [verify_line("exact65521", 65521), verify_line("exact65537", 65537)],
                  "limits": limits(), "seconds": time.perf_counter()-tick, "auditor_sha256": sha(__file__)}
        write(RESULTS / args.output, report)
        print({"verified": True, "seconds": report["seconds"]}, flush=True)
