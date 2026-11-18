"""汇总真实状态，保留实际拒绝候选，生成交付摘要；不宣称取得最小多项式。"""

from fractions import Fraction as F
import json
import sys

from refine import HERE, sha, write_new


def main():
    budgets = json.loads((HERE / "search-budget.json").read_text(encoding="utf-8"))
    runs = json.loads((HERE / "search-runs.json").read_text(encoding="utf-8"))
    accepted = json.loads((HERE / "lattices-verified.json").read_text(encoding="utf-8"))
    controls = json.loads((HERE / "controls-verified.json").read_text(encoding="utf-8"))
    rejection = json.loads((HERE / "acceptor-rejection-tests.json").read_text(encoding="utf-8"))
    root = json.loads((HERE / "refined-root-verified.json").read_text(encoding="utf-8"))
    compression = json.loads((HERE / "artifact-compression.json").read_text(encoding="utf-8"))
    cases = []
    for run in runs["lll"]:
        if run["status"] != "completed":
            cases.append({"name": run["name"], "degree_attempted": run["degree"],
                          "height_attempted": run["coefficient_bound_inclusive"],
                          "status": run["status"], "exclusion_proved": False})
            continue
        path = HERE / "lattices" / (run["name"]+".verified.json")
        data = json.loads(path.read_text(encoding="utf-8"))
        if not data["verified"] or not data["no_nonzero_integer_polynomial_in_budget_proved"]:
            raise ValueError("已完成病例缺少精确排除")
        if data["candidates_within_height"]:
            raise ValueError("实际候选与零候选交付不同")
        if run["degree"] == 64 and F(data["gs_squared_norm_lower"]) <= 10**76:
            raise ValueError("报告中的共同下界不成立")
        cases.append({"name": run["name"], "degree_excluded": run["degree"],
                      "height_excluded_inclusive": run["coefficient_bound_inclusive"],
                      "status": "exact_exclusion_verified", "exclusion_proved": True,
                      "proof_file": path.relative_to(HERE).as_posix(),
                      "proof_sha256": sha(path)})
    if not all((accepted["all_completed_cases_verified"], controls["verified"], rejection["verified"], root["verified"])):
        raise ValueError("前置核验缺失")
    sample = json.loads((HERE / "lattices/r-d16.verified.json").read_text(encoding="utf-8"))
    write_new(HERE / "rejected-r-d16.json",
              {"is_target_defining_polynomial": False, "status": "rejected_exactly_nonzero_and_height_exceeded",
               "degree": 16, "coefficient_order": "low_to_high", "budget_height": str(10**24),
               **sample["shortest_basis_row"], "source_proof_sha256": sha(HERE / "lattices/r-d16.verified.json")})
    write_new(HERE / "delivery.json",
              {"schema": "target-reconstruction-delivery-v1", "client_date": "2026-10-11",
               "timezone": "Asia/Shanghai", "phase_search_and_report_complete": True,
               "exact_algebraic_representation_obtained": False,
               "minimal_polynomial_obtained": False,
               "refined_multivariate_root_box_halfwidth": "10^-5400",
               "refined_box_is_same_parent_root_proved": True,
               "cases": cases, "pslq_cases": runs["pslq"],
               "q_candidates_within_budget": [], "accepted_q": None,
               "coordinate_and_multiplier_rational_expressions": None,
               "modulo_q_original_equations": {"status": "not_performed_no_q", "checked": 0, "passed": None},
               "selected_q_real_root_maps_into_parent_box": {"status": "not_performed_no_q", "passed": None},
               "auxiliary_selection_denominators_nonzero_verified": ["tangent6", "tangent0", "mobius_r"],
               "radius_minimal_polynomial": None,
               "source_minimal_polynomial_degree_guessed": False,
               "complete_groebner_basis_computed": False,
               "known_algebraic_controls_passed": [v["name"] for v in controls["reports"]],
               "rejection_tests_passed": len(rejection["damaged_certificates"]),
               "search_seconds": runs["seconds"],
               "raw_artifacts_losslessly_compressed": True,
               "raw_artifacts_bytes": compression["original_bytes"],
               "compressed_artifacts_bytes": compression["archive_bytes"],
               "prior_screenshot_claim_proved_by_original_experiment": False,
               "finite_exclusion_does_not_rule_out_other_heights_degrees_or_generators": True})
    files = []
    for path in sorted(HERE.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or path.name in ("SHA256SUMS.json", "delivery-replayed.json"):
            continue
        files.append({"path": path.relative_to(HERE).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)})
    write_new(HERE / "SHA256SUMS.json", {"schema": "target-reconstruction-manifest-v1", "files": files})
    print(json.dumps({"exact_exclusions": sum(v["exclusion_proved"] for v in cases),
                      "unfinished": [v["name"] for v in cases if not v["exclusion_proved"]],
                      "files": len(files), "bytes": sum(v["bytes"] for v in files)}), flush=True)


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
