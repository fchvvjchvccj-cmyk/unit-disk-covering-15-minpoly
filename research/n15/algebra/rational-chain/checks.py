"""回放活动证书，并确认错误分支、导数、图表和预条件器会被拒绝。"""

import argparse
from copy import deepcopy
import time
from unittest.mock import patch

from support import HERE, AXIS, TARGET, FACTOR, read, write, sha, bounded
import prove_branch
import verify_fiber
import verify_system
import target_chart
import quadratic_system


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    destination = HERE / args.output
    if destination.exists():
        raise FileExistsError("拒绝覆盖已有接触链检查")
    with bounded() as snapshot:
        tick = time.perf_counter()
        system, root = read(AXIS), read(TARGET / "refined-root.json")
        fiber_folder = HERE / "fiber-runs/20261015-p65521"
        metadata, fiber = read(fiber_folder / "field.json"), read(fiber_folder / "fiber-0.json.gz")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        sources = [HERE / "system-runs" / name for name in ("20261015-critical9-final", "20261015-critical10-final")]
        circuits = [read(path / "circuit.json.gz") for path in sources]
        certificates = [read(path / "20261015-target-bound/certificate.json.gz") for path in sources]
        positive = []
        positive.append(prove_branch.verify(root, system)["verified"])
        positive.append(quadratic_system.prove_left_identity()["verified"])
        positive.append(verify_fiber.verify(fiber, field, system)["verified"])
        for data, certificate in zip(circuits, certificates):
            positive.append(verify_system.verify(data, fiber, field)["verified"])
            positive.append(target_chart.verify(certificate, data, root)["verified"])
        if not all(positive):
            raise ValueError("正例证书回放失败")
        rejected = []

        def reject(name, action):
            try:
                action()
            except (ValueError, ZeroDivisionError) as error:
                rejected.append({"name": name, "reason": str(error)})
                return
            raise AssertionError("错误材料未被拒绝："+name)

        bad = deepcopy(system)
        bad["equations"][1][0]["coefficient"] += 1
        reject("changed_integer_geometry_identity", lambda: prove_branch.source_identities(bad))
        bad_root = deepcopy(root)
        bad_root["xnum"][13] = str(int(bad_root["xnum"][13])+int(bad_root["Qx"])//10000)
        reject("wrong_radius_in_root_image", lambda: prove_branch.verify(bad_root, system))
        bad_fiber = deepcopy(fiber)
        bad_fiber["left"][1] = deepcopy(bad_fiber["left"][0])
        reject("duplicate_left_carrier_root", lambda: verify_fiber.verify(bad_fiber, field, system))
        bad_fiber = deepcopy(fiber)
        gradient = bad_fiber["bridge_values"][0]["gradient"][0]
        if gradient:
            gradient[0] = (gradient[0]+1)%metadata["prime"]
        else:
            gradient.append(1)
        reject("changed_bridge_tangent", lambda: verify_fiber.verify(bad_fiber, field, system))
        bad_circuit = deepcopy(circuits[0])
        zero = next(i for i, node in enumerate(bad_circuit["nodes"]) if node == {"op": "const", "value": 0})
        bad_circuit["outputs"]["critical_equations"][1] = zero
        reject("missing_total_derivative", lambda: verify_system.verify(bad_circuit, fiber, field))
        reject("changed_circuit_binding", lambda: target_chart.verify(certificates[0], bad_circuit, root))
        bad_circuit = deepcopy(circuits[0])
        bad_circuit["nodes"][0] = {"op": "add", "args": [1, 2]}
        reject("forward_reference_in_circuit", lambda: verify_system.execute(bad_circuit, [field(1)]*9))
        bad_certificate = deepcopy(certificates[0])
        bad_certificate["inverse_numerators"] = [["0"]*9 for _ in range(9)]
        reject("zero_jacobian_preconditioner", lambda: target_chart.verify(bad_certificate, circuits[0], root))
        bad_certificate = deepcopy(certificates[1])
        bad_certificate["dependencies"]["left_factor_selection"]["sha256"] = "0"*64
        reject("unbound_parent_factor_proof", lambda: target_chart.verify(bad_certificate, circuits[1], root))
        bad_factor = deepcopy(read(FACTOR))
        bad_factor["factor_polynomials"][1][0]["coefficient"] += 1
        with patch.object(quadratic_system, "read", return_value=bad_factor):
            reject("wrong_common_height_factor", quadratic_system.prove_left_identity)
        # 非零的偏导证明必须使用所绑定的实际电路，不能换成恒零门。
        bad_circuit = deepcopy(circuits[1])
        zero = next(i for i, node in enumerate(bad_circuit["nodes"]) if node == {"op": "const", "value": 0})
        bad_circuit["critical_locus_nonzero_guards"] = [zero]
        reject("zero_radius_direction_guard", lambda: target_chart.verify(certificates[1], bad_circuit, root))
        report = {"verified": True, "positive_replays": len(positive), "negative_controls_rejected": rejected,
                  "root_sha256": sha(TARGET / "refined-root.json"), "source_system_sha256": sha(AXIS),
                  "circuit_sha256": [sha(path / "circuit.json.gz") for path in sources],
                  "driver_sha256": sha(__file__), "seconds": time.perf_counter()-tick, "limits": snapshot()}
        write(destination, report)
        print({"verified": True, "positive_replays": len(positive), "rejected": len(rejected),
               "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="checks-final-verified.json")
    run(parser.parse_args())
