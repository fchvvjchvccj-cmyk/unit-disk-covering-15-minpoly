"""拒绝错误的投影电路、恢复公式、临界乘子和目标证书。"""

import argparse
from copy import deepcopy
import time

from descent_support import HERE, SHARED, TARGET, FACTOR, AXIS, read, write, sha, bounded
from verify_projection import verify as verify_target
from prove_relative import verify as verify_identities
from verify_geometry import product_check, closed_check
from quadratic_pair import Pair
from chain_fiber import Jet
from fiber import compute
from replay_norm import compare
from verify_quotient import setup, decode_scalar


def pair_controls(field):
    scalar = lambda value: Jet(field(value))
    a = Jet(field(3), list(map(field, [1, 0, 0, 0])))
    b = Jet(field(5), list(map(field, [0, 1, 0, 0])))
    V = Jet(field(4), list(map(field, [2, 3, 5, 7])))
    q = Pair(a, b, V)
    count = 0
    for sign in (1, -1):
        root = Jet(field(sign*2), [d/(field(sign*4)) for d in V.gradient])
        original = a+b*root
        for actual, expected in [(q.inverse().evaluate(root), 1/original),
                                 ((q*q).evaluate(root), original*original)]:
            if (actual.value, actual.gradient) != (expected.value, expected.gradient):
                raise ValueError("系数对运算的值或四个导数不一致")
            count += 1
    V2 = q.cast(4)
    outer = Pair(q, q.cast(1), V2)
    product = outer*outer.inverse()
    if not (product.a-q.cast(1)).is_zero() or not product.b.is_zero():
        raise ValueError("嵌套二次层的求逆没有得到单位")
    return {"verified": True, "single_level_value_and_jet_checks": count, "nested_inverse_checked": True}


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    with bounded() as snapshot:
        tick = time.perf_counter()
        source = HERE / "system-runs/20261015-projected4-final"
        data, certificate = read(source / "circuit.json.gz"), read(source / "target-final-certificate.json.gz")
        root, axis = read(TARGET / "refined-root.json"), read(AXIS)
        folder = SHARED.parent / "rational-chain/fiber-runs/20261015-p65521"
        metadata = read(folder / "field.json")
        field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(65521)(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
        fiber = read(folder / "fiber-0.json.gz")
        closed_field = fq_default_ctx(65521, 8, "a", fq_type="FQ_NMOD")
        closed = compute(closed_field, [2, 3, 1, 1552], read(FACTOR))
        positive = [verify_identities()["verified"], verify_target(certificate, data, root)["verified"],
                    product_check(data, fiber, field, axis)["verified"], closed_check(data, closed, closed_field, axis)["verified"]]
        pairs = pair_controls(field)
        positive.append(pairs["verified"])
        if not all(positive):
            raise ValueError("投影正例回放失败")
        print({"stage": "positive_replays", "count": len(positive), "seconds": time.perf_counter()-tick}, flush=True)
        rejected = []

        def reject(name, action):
            try:
                action()
            except (ValueError, ZeroDivisionError) as error:
                rejected.append({"name": name, "reason": str(error)})
                return
            raise AssertionError("错误材料未被拒绝："+name)

        zero = next(i for i, node in enumerate(data["nodes"]) if node == {"op": "const", "value": 0})
        bad = deepcopy(data)
        bad["outputs"]["critical_equations"][1] = zero
        reject("changed_critical_derivative_binding", lambda: verify_target(certificate, bad, root))
        bad = deepcopy(data)
        bad["closed_chart_nonzero_guards"].pop()
        reject("missing_recovery_guard", lambda: verify_target(certificate, bad, root))
        wrong = deepcopy(certificate)
        wrong["circuit_content_sha256"] = "0"*64
        reject("wrong_circuit_content_binding", lambda: verify_target(wrong, data, root))
        wrong_root = deepcopy(root)
        wrong_root["xnum"][13] = str(int(wrong_root["xnum"][13])+int(wrong_root["Qx"])//10000)
        reject("unbound_radius_root", lambda: verify_target(certificate, data, wrong_root))
        wrong = deepcopy(certificate)
        wrong["inverse_numerators"] = [["0"]*4 for _ in range(4)]
        reject("zero_four_dimensional_preconditioner", lambda: verify_target(wrong, data, root))
        bad = deepcopy(data)
        bad["outputs"]["critical_equations"][0] = zero
        reject("wrong_generic_norm_product", lambda: product_check(bad, fiber, field, axis))
        bad = deepcopy(data)
        bad["outputs"]["recovered_left_and_y"]["s"] = zero
        reject("wrong_closed_generator_recovery", lambda: closed_check(bad, closed, closed_field, axis))
        bad = deepcopy(data)
        bad["outputs"]["geometry"]["c0"][0] = zero
        reject("wrong_closed_geometry_recovery", lambda: closed_check(bad, closed, closed_field, axis))
        bad = deepcopy(data)
        bad["outputs"]["critical_multiplier"] = zero
        reject("wrong_closed_critical_multiplier", lambda: closed_check(bad, closed, closed_field, axis))
        scalar = lambda value: Jet(field(value))
        radicand = scalar(4)
        reject("singular_quadratic_pair_inverse", lambda: Pair(scalar(2), scalar(1), radicand).inverse())
        foreign = Pair(scalar(1), scalar(1), scalar(5))
        reject("mixed_quadratic_domains", lambda: Pair(scalar(1), scalar(1), radicand).cast(foreign))
        setup(65521)
        functions = read(HERE / "norm-runs/20261015-p65521-final/functions.json.gz")
        expected = [decode_scalar(row, 65521) for row in functions["G_and_gradient_u_v_w_r"]]
        actual = list(expected)
        actual[1] += 1
        reject("changed_symbolic_norm_derivative", lambda: compare(actual, expected))
        report = {"verified": True, "positive_replays": len(positive), "negative_controls_rejected": rejected,
                  "pair_arithmetic_controls": pairs, "circuit_sha256": sha(source / "circuit.json.gz"),
                  "certificate_sha256": sha(source / "target-final-certificate.json.gz"),
                  "driver_sha256": sha(__file__), "seconds": time.perf_counter()-tick, "limits": snapshot()}
        write(HERE / args.output, report)
        print({"verified": True, "positive_replays": len(positive), "negative_controls_rejected": len(rejected),
               "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="checks-verified.json.gz")
    run(parser.parse_args())
