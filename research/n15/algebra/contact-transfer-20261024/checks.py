"""以损坏证书、分裂分支和幂零元检查本轮接受器是否会误接收。"""

from copy import deepcopy
from fractions import Fraction as Q

from context import INPUTS, RESULTS, read, write, bounded, sha, execute
from verify import verify_sources, verify_new, require
from audit import FiniteAlgebra, columns, span, verify_line


def reject(action, label, reports):
    try:
        action()
    except (ValueError, ZeroDivisionError):
        reports.append(label)
        return
    raise ValueError("损坏数据未被拒绝："+label)


def algebra_controls():
    from flint import fq_default_ctx
    field = fq_default_ctx(65521, 1, "a", fq_type="FQ_NMOD")
    zero, one = field(0), field(1)
    relations = []
    for j in range(4):
        a, b = [zero]*2**j, [zero]*2**j
        b[0] = one
        relations.append((tuple(a), tuple(b)))
    A = FiniteAlgebra(field, relations)
    t, complement = tuple([zero, one]+[zero]*14), tuple([one, -one]+[zero]*14)
    rank_t = span([None]*16, columns(A, t))
    rank_complement = span([None]*16, columns(A, complement))
    combined = [None]*16
    rank_joint = span(combined, columns(A, t)+columns(A, complement))
    require((rank_t, rank_complement, rank_joint) == (8, 8, 16), "不同分支零点被误认为共同零点")
    relations[0] = ((zero,), (zero,))
    A = FiniteAlgebra(field, relations)
    raw = span([None]*16, columns(A, t))
    stable = t
    for _ in range(4):
        stable = A.multiply(stable, stable, 4)
    stable_rank = span([None]*16, columns(A, stable))
    require((raw, stable_rank) == (8, 0), "稳定幂未清除幂零部分")
    rejected = []
    reject(lambda: A.inverse((zero,), 0), "zero_field_inverse", rejected)
    return {"split_individual_zero_dimensions": [8, 8], "split_common_zero_dimension": 0,
            "nilpotent_raw_rank": raw, "nilpotent_stable_rank": stable_rank, "zero_inverse_rejected": bool(rejected)}


def run():
    _, root, center, _ = verify_sources()
    data = read(RESULTS / "primitive16.json.gz")
    certificate = read(RESULTS / "target-final-certificate.json.gz")
    positive, _ = verify_new(data, center, Q(root["rho"]), certificate)
    require(positive["unique_root_in_box"], "正证书重放失败")
    rejected = []
    bad = deepcopy(certificate)
    bad["inverse"] = [["0"]*8 for _ in range(8)]
    reject(lambda: verify_new(data, center, Q(root["rho"]), bad), "zero_preconditioner", rejected)
    bad_center = deepcopy(certificate)
    bad_center["center"][0] = str(Q(bad_center["center"][0])+1)
    reject(lambda: verify_new(data, center, Q(root["rho"]), bad_center), "unbound_root_center", rejected)
    bad_radius = deepcopy(certificate)
    bad_radius["radius"] = "-1"
    reject(lambda: verify_new(data, center, Q(root["rho"]), bad_radius), "negative_radius", rejected)
    bad_circuit = deepcopy(data)
    bad_circuit["nodes"][0]["value"] = 7
    reject(lambda: verify_new(bad_circuit, center, Q(root["rho"]), certificate), "changed_circuit_port", rejected)
    reject(lambda: execute({"variables": ["x"], "nodes": [{"op": "neg", "args": [0]}]}, [Q(1)]),
           "forward_circuit_reference", rejected)
    lines = [verify_line("exact65521", 65521), verify_line("exact65537", 65537)]
    require(all(all(q.get("critical_dimension_on_required_open_chart") == 0 for q in line["common_factor_classifications"])
                for line in lines), "当前截面的公因子仍有未说明分支")
    return {"verified": True, "root_certificate_replayed": True, "rejected": rejected,
            "algebra_controls": algebra_controls(), "both_exact_parameter_lines_replayed": True,
            "checks_sha256": sha(__file__)}


if __name__ == "__main__":
    with bounded():
        result = run()
        write(RESULTS / "checks.json", result)
        print({"verified": True, "rejections": len(result["rejected"]), "algebra_controls": 3})
