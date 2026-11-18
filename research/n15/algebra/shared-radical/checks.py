"""回放活动证明，拒绝错误分支、规范形、因子及伪临界判据。"""

import argparse
from copy import deepcopy
import time

from bridge_support import HERE, PARENT, TARGET, AXIS, read, write, sha, bounded
from prove_shared import verify as verify_selection
from build_shared import construct
from verify_shared import check_fiber, verify_target
from verify_quotient import verify_symbolic, verify_samples
from verify_system import execute
from analyze_factors import FiniteAlgebra, columns, append_span, stable_power, check_factor


def algebra_controls():
    from flint import fq_default_ctx
    field = fq_default_ctx(65521, 1, "a", fq_type="FQ_NMOD")
    zero, one = field(0), field(1)
    relations = []
    for j in range(5):
        # 首层含非零幂零元，后四层取分裂的单根二次式。
        a, b = [zero]*2**j, [zero]*2**j
        if j:
            b[0] = one
        relations.append((tuple(a), tuple(b)))
    algebra = FiniteAlgebra(field, relations)
    t = tuple([zero, one]+[zero]*30)
    raw_rank = append_span([None]*32, columns(algebra, t))
    stable_rank = append_span([None]*32, columns(algebra, stable_power(algebra, t)))
    unit = tuple([one]+[zero]*31)
    unit_rank = append_span([None]*32, columns(algebra, stable_power(algebra, unit)))
    if (raw_rank, stable_rank, unit_rank) != (16, 0, 32):
        raise AssertionError("稳定幂未正确排除幂零部分或保留单位部分")
    a, b = relations[0]
    relations[0] = (a, (one,))
    split = FiniteAlgebra(field, relations)
    basis = [None]*32
    rank_H = append_span(basis, columns(split, t))
    rank_joint = append_span(basis, columns(split, split.subtract(unit, t)))
    if (rank_H, rank_joint) != (16, 32):
        raise AssertionError("未区分不同分支上的独立零点与同分支公共零点")
    return {"verified": True, "nilpotent_raw_rank": raw_rank, "nilpotent_stable_rank": stable_rank,
            "unit_stable_rank": unit_rank, "split_individual_zero_dimension": 16,
            "split_joint_zero_dimension": 0}


def increment(encoded):
    if encoded["n"]:
        encoded["n"][0] = (encoded["n"][0]+1)%65521
    else:
        encoded["n"] = [1]


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    destination = HERE / args.output
    if destination.exists():
        raise FileExistsError("拒绝覆盖已有共享根检查")
    with bounded() as snapshot:
        tick = time.perf_counter()
        source = HERE / "system-runs/20261015-selected"
        circuit = read(source / "circuit.json.gz")
        certificate = read(source / "final-target-certificate.json.gz")
        root, system = read(TARGET / "refined-root.json"), read(AXIS)
        folder = PARENT / "fiber-runs/20261015-p65521"
        metadata, fiber = read(folder / "field.json"), read(folder / "fiber-0.json.gz")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(65521)(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        quotient_folder = HERE / "quotient-runs/20261015-jets-p65521"
        data, norm = read(quotient_folder / "quotient.json.gz"), read(quotient_folder / "norm.json.gz")
        positive = [verify_selection(root)["verified"], verify_target(certificate, circuit, root)["verified"],
                    check_fiber(circuit, fiber, field, system)["verified"], verify_symbolic(data, norm, circuit)["verified"],
                    verify_samples(data, norm, circuit, [7, 11])["verified"]]
        algebra_report = algebra_controls()
        positive.append(algebra_report["verified"])
        if not all(positive):
            raise ValueError("共享根正例回放失败")
        rejected = []

        def reject(name, action):
            try:
                action()
            except (ValueError, ZeroDivisionError) as error:
                rejected.append({"name": name, "reason": str(error)})
                return
            raise AssertionError("错误材料未被拒绝："+name)

        reject("opposite_shared_sign", lambda: verify_target(certificate, construct(-1), root))
        bad_root = deepcopy(root)
        bad_root["xnum"][13] = str(int(bad_root["xnum"][13])+int(bad_root["Qx"])//10000)
        reject("unbound_target_radius", lambda: verify_target(certificate, circuit, bad_root))
        bad_certificate = deepcopy(certificate)
        bad_certificate["inverse_numerators"] = [["0"]*9 for _ in range(9)]
        reject("zero_jacobian_preconditioner", lambda: verify_target(bad_certificate, circuit, root))
        bad_circuit = deepcopy(circuit)
        bad_circuit["nodes"][0] = {"op": "add", "args": [1, 2]}
        reject("forward_circuit_reference", lambda: execute(bad_circuit, [field(1)]*9))
        zero = next(i for i, node in enumerate(circuit["nodes"]) if node == {"op": "const", "value": 0})
        bad_circuit = deepcopy(circuit)
        bad_circuit["critical_locus_nonzero_guards"] = [zero]
        reject("zero_radius_direction_guard", lambda: verify_target(certificate, bad_circuit, root))
        bad_fiber = deepcopy(fiber)
        bad_fiber["left"][1] = deepcopy(bad_fiber["left"][0])
        reject("duplicate_carrier_root", lambda: check_fiber(circuit, bad_fiber, field, system))
        bad = deepcopy(data)
        bad["gradient_u_v_w_r"].pop()
        reject("missing_total_derivative", lambda: verify_symbolic(bad, norm, circuit))
        bad = deepcopy(data)
        increment(bad["H"][0])
        reject("changed_closure_coefficient", lambda: verify_symbolic(bad, norm, circuit))
        bad = deepcopy(data)
        increment(bad["gradient_u_v_w_r"][0][0])
        reject("changed_derivative_coefficient", lambda: verify_symbolic(bad, norm, circuit))
        bad = deepcopy(data)
        increment(bad["relations"][0]["a"][0])
        reject("changed_tower_relation", lambda: verify_symbolic(bad, norm, circuit))
        bad = deepcopy(data)
        increment(bad["geometry"]["c4"][0][0])
        reject("changed_coordinate_coefficient", lambda: verify_symbolic(bad, norm, circuit))
        bad_norm = deepcopy(norm)
        increment(bad_norm["value"])
        reject("changed_norm_coefficient", lambda: verify_symbolic(data, bad_norm, circuit))
        bad_norm = deepcopy(norm)
        bad_norm["numerator_factors"][0]["multiplicity"] += 1
        reject("wrong_factor_multiplicity", lambda: verify_symbolic(data, bad_norm, circuit))
        bad_norm = deepcopy(norm)
        bad_norm["numerator_factors"][0].update(coefficients=[0, 0, 1], degree=2)
        reject("composite_factor_label", lambda: verify_symbolic(data, bad_norm, circuit))
        reject("factor_with_coefficient_pole", lambda: check_factor(data, norm["numerator_factors"][7]))
        report = {"verified": True, "positive_replays": len(positive), "negative_controls_rejected": rejected,
                  "finite_algebra_controls": algebra_report, "driver_sha256": sha(__file__),
                  "circuit_sha256": sha(source / "circuit.json.gz"), "quotient_sha256": sha(quotient_folder / "quotient.json.gz"),
                  "norm_sha256": sha(quotient_folder / "norm.json.gz"), "seconds": time.perf_counter()-tick, "limits": snapshot()}
        write(destination, report)
        print({"verified": True, "positive_replays": len(positive), "negative_controls_rejected": len(rejected),
               "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="checks-verified.json")
    run(parser.parse_args())
