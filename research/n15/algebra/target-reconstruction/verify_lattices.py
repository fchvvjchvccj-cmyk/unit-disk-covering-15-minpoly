"""不信任LLL；重验整数基、逆矩阵、向外舍入Gram-Schmidt及有限排除范围。"""

from fractions import Fraction as F
import json
import sys
import time

from intervals import Dyadic
from refine import HERE, sha, write_new, artifact_bytes, artifact_sha
from search import SEEDS
from windows_compute_limits import apply_limits


def scalar_interval(definition, values, arithmetic):
    def linear(terms):
        total = arithmetic.zero
        for index, coefficient in terms.items():
            value = arithmetic.one if index == "constant" else values[int(index)]
            total = arithmetic.add(total, arithmetic.integer(value, coefficient))
        return total
    kind = definition["kind"]
    if kind == "coordinate":
        return values[definition["index"]], None
    if kind == "square":
        return arithmetic.power(values[definition["index"]], 2), None
    if kind == "linear":
        return linear(definition["terms"]), None
    denominator = linear(definition["denominator"])
    return arithmetic.div(linear(definition["numerator"]), denominator), denominator


def gs_lower_bounds(basis, bits):
    """Gram矩阵的区间LDL分解；对角元为正交化向量长度平方。"""
    arithmetic = Dyadic(bits)
    size = len(basis)
    gram = [[sum(x*y for x, y in zip(basis[i], basis[j])) for j in range(i+1)] for i in range(size)]
    lower = [[arithmetic.zero]*size for _ in range(size)]
    diagonal = []
    for i in range(size):
        for j in range(i):
            value = arithmetic.rational(gram[i][j])
            for k in range(j):
                value = arithmetic.sub(value, arithmetic.mul(arithmetic.mul(lower[i][k], lower[j][k]), diagonal[k]))
            lower[i][j] = arithmetic.div(value, diagonal[j])
        value = arithmetic.rational(gram[i][i])
        for k in range(i):
            value = arithmetic.sub(value, arithmetic.mul(arithmetic.power(lower[i][k], 2), diagonal[k]))
        if value[0] <= 0:
            raise ValueError("正交长度下界不为正，须提高验证精度")
        diagonal.append(value)
    return diagonal, arithmetic.scale


def polynomial_interval(coefficients, value, arithmetic):
    total = arithmetic.zero
    for coefficient in reversed(coefficients):
        total = arithmetic.add(arithmetic.mul(total, value), arithmetic.rational(coefficient))
    return total


def verify_case(case, root, config):
    started = time.perf_counter()
    name = case["name"]
    degree = case["degree"]
    size = degree+1
    input_path = HERE / "lattices" / (name+".input.txt")
    data_path = HERE / "lattices" / (name+".json")
    if artifact_sha(input_path) != case["input_sha256"] or artifact_sha(data_path) != case["data_sha256"]:
        raise ValueError("格搜索产物摘要改变")
    entries = list(map(int, artifact_bytes(input_path).decode("ascii").split()))
    if entries[0] != size or len(entries) != size+1:
        raise ValueError("格输入维数错误")
    entries = entries[1:]
    proposal = json.loads(artifact_bytes(data_path).decode("utf-8"))
    basis = [list(map(int, row)) for row in proposal["basis"]]
    inverse = [list(map(int, row)) for row in proposal["inverse_transform"]]
    if len(basis) != size or any(len(row) != size+1 for row in basis):
        raise ValueError("整数基维数错误")
    if len(inverse) != size or any(len(row) != size for row in inverse):
        raise ValueError("逆变换维数错误")
    # 基的前n列就是整数变换U；最后一列必须逐项等于U乘输入整数向量。
    for row in basis:
        if sum(x*y for x, y in zip(row[:-1], entries)) != row[-1]:
            raise ValueError("基不属于原整数格")
    columns = list(zip(*inverse))
    for i, row in enumerate(basis):
        for j, column in enumerate(columns):
            if sum(x*y for x, y in zip(row[:-1], column)) != (i == j):
                raise ValueError("整数变换不为幺模矩阵")
    arithmetic = Dyadic(root["interval_bits"])
    center = [F(int(v), int(root["Qx"])) for v in root["xnum"]]
    radius = F(root["rho"])
    values = [(arithmetic.rational(v-radius)[0], arithmetic.rational(v+radius)[1]) for v in center]
    scalar, denominator = scalar_interval(SEEDS[case["seed"]], values, arithmetic)
    scale = 10**case["scale_digits"]
    powers = [arithmetic.one]
    errors = []
    for i, entry in enumerate(entries):
        if i:
            powers.append(arithmetic.mul(powers[-1], scalar))
        scaled = arithmetic.integer(powers[-1], scale)
        error = arithmetic.absolute(arithmetic.sub(scaled, arithmetic.rational(entry)))
        if error > arithmetic.scale:
            raise ValueError("认证根的幂与格输入误差大于1")
        errors.append(error)
    diagonal, gs_scale = gs_lower_bounds(basis, config["gram_schmidt_interval_bits"])
    minimum = min(row[0] for row in diagonal)
    height = int(config["coefficient_bound_inclusive"])
    threshold = (size+size**2)*height**2
    excluded = minimum > threshold*gs_scale
    candidates = []
    for row in basis:
        coefficients = row[:-1]
        if max(map(abs, coefficients)) <= height:
            interval = polynomial_interval(coefficients, scalar, arithmetic)
            candidates.append({"coefficients_low_to_high": list(map(str, coefficients)),
                               "interval_contains_zero": interval[0] <= 0 <= interval[1],
                               "exact_value_interval": [str(F(v, arithmetic.scale)) for v in interval]})
    best = min(basis, key=lambda row: sum(v*v for v in row))
    interval = polynomial_interval(best[:-1], scalar, arithmetic)
    return {"name": name, "verified": True, "degree_budget": degree,
            "height_budget_inclusive": str(height), "basis_and_unimodularity_verified": True,
            "power_rounding_error_upper": str(F(max(errors), arithmetic.scale)),
            "scalar_interval": [str(F(v, arithmetic.scale)) for v in scalar],
            "denominator_nonzero_interval": None if denominator is None else [str(F(v, arithmetic.scale)) for v in denominator],
            "gs_squared_norm_lower": str(F(minimum, gs_scale)),
            "relation_squared_norm_upper_if_zero": str(threshold),
            "no_nonzero_integer_polynomial_in_budget_proved": excluded,
            "candidates_within_height": candidates,
            "shortest_basis_row": {"coefficients_low_to_high": list(map(str, best[:-1])),
                                   "height": str(max(map(abs, best[:-1]))),
                                   "last_lattice_coordinate": str(best[-1]),
                                   "exact_value_interval": [str(F(v, arithmetic.scale)) for v in interval],
                                   "exact_interval_excludes_zero": not interval[0] <= 0 <= interval[1]},
            "flint_version": proposal["flint_version"], "seconds": time.perf_counter()-started}


def main():
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        root = json.loads((HERE / "refined-root.json").read_text(encoding="utf-8"))
        accepted = json.loads((HERE / "refined-root-verified.json").read_text(encoding="utf-8"))
        if not accepted["verified"] or accepted["refined_root_sha256"] != sha(HERE / "refined-root.json"):
            raise ValueError("高精度根未接受或摘要不符")
        config = json.loads((HERE / "search-budget.json").read_text(encoding="utf-8"))
        if config["seeds"] != SEEDS or config["root_sha256"] != sha(HERE / "refined-root.json"):
            raise ValueError("搜索配置或根改变")
        runs = json.loads((HERE / "search-runs.json").read_text(encoding="utf-8"))
        expected = [(v["seed"], v["degree"], v["scale_digits"]) for v in config["cases"]]
        actual = [(v["seed"], v["degree"], v["scale_digits"]) for v in runs["lll"]]
        if actual != expected or len({v["name"] for v in runs["lll"]}) != len(expected):
            raise ValueError("实际搜索病例与预算清单不同")
        reports = []
        for case in runs["lll"]:
            if case["status"] != "completed":
                reports.append({"name": case["name"], "verified": False, "status": case["status"]})
                continue
            report = verify_case(case, root, config)
            write_new(HERE / "lattices" / (case["name"]+".verified.json"), report)
            reports.append({k: report[k] for k in ("name", "verified", "no_nonzero_integer_polynomial_in_budget_proved", "seconds")})
            print(json.dumps(reports[-1]), flush=True)
        write_new(HERE / "lattices-verified.json",
                  {"schema": "target-lattice-verification-v1", "reports": reports,
                   "all_completed_cases_verified": all(v["verified"] for v in reports if v.get("status") != "wall_timeout"),
                   "completed_cases": sum(v["verified"] for v in reports),
                   "unfinished_cases": [v["name"] for v in reports if not v["verified"]],
                   "all_budgets_excluded": all(v.get("no_nonzero_integer_polynomial_in_budget_proved", False) for v in reports),
                   "search_runs_sha256": sha(HERE / "search-runs.json"),
                   "search_budget_sha256": sha(HERE / "search-budget.json"),
                   "root_sha256": sha(HERE / "refined-root.json"),
                   "seconds": time.perf_counter()-started, "windows_job_limits": limits(),
                   "minimal_polynomial_obtained": False})
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
