"""用可解样例及故意损坏的证书检验接受端是否会误接收。"""

from copy import deepcopy
import argparse
import json
from pathlib import Path
import tempfile
import time

from common import HERE, AXIS, TARGET, bounded, read, write


def must_reject(action, name, reports):
    try:
        action()
    except (ValueError, ZeroDivisionError) as error:
        reports.append({"control": name, "rejected": True, "reason": str(error)})
    else:
        raise AssertionError("接受端错误地接受了损坏证书："+name)


def line_controls(reports):
    from flint import nmod_poly
    from line_profile import module_certificate
    from verify_line import verify
    prime = 65521
    P = lambda a: nmod_poly(a, prime)
    X, N, D = P([0, 1]), P([3, 2, 1]), P([1, 4])
    samples = []
    for t in range(2, 10):
        dn = int(D(t))
        f = int(N(t))*pow(dn, -1, prime)%prime
        g = int((N.derivative()*D-N*D.derivative())(t))*pow(dn*dn, -1, prime)%prime
        samples.append({"r": t, "value": f, "derivative_r": g})
    M, Y = P([1]), P([])
    for sample in samples:
        t = sample["r"]
        Li = P([1])
        for other in samples:
            if other["r"] != t:
                Li *= (X-other["r"])*pow(t-other["r"], -1, prime)
        Y += ((1-2*int(Li.derivative()(t))*(X-t))*sample["value"]+(X-t)*sample["derivative_r"])*Li**2
        M *= (X-t)**2
    for budget, expected in ((2, False), (1, True)):
        certificate = module_certificate(M, Y, prime, budget, len(samples))
        report = verify(certificate, samples)
        if report["no_bounded_rational_function_matches_recorded_jets"] != expected:
            raise AssertionError("已知二次/一次有理函数的次数判断错误")
        reports.append({"control": f"known_rational_degree_budget_{budget}", "passed": True})
    certificate = module_certificate(M, Y, prime, 2, len(samples))
    bad = deepcopy(certificate)
    bad["unimodular_transform"][0][0] = [0]
    must_reject(lambda: verify(bad, samples), "line_corrupted_transform", reports)
    bad = deepcopy(samples)
    bad[0]["derivative_r"] = (bad[0]["derivative_r"]+1)%prime
    must_reject(lambda: verify(certificate, bad), "line_corrupted_derivative", reports)


def lattice_controls(reports):
    from flint import fmpz_mat
    from height_search import verify
    entries = [10**6, 5*10**5, 25*10**4]
    matrix = fmpz_mat([[int(i == j) for j in range(3)]+[entries[i]] for i in range(3)])
    basis, U = matrix.lll(transform=True)
    inverse = U.inv()
    fixture = {"flint_version": "control-only", "basis": [[int(basis[i, j]) for j in range(4)] for i in range(3)],
               "inverse_transform": [[int(inverse[i, j]) for j in range(3)] for i in range(3)]}
    root = {"interval_bits": 128, "xnum": [0]*13+[1], "Qx": 2, "rho": "0"}
    config = {"degree": 2, "seed": "r", "scale_digits": 6, "height_digits": 1}
    with tempfile.TemporaryDirectory(prefix="chain-controls-") as temporary:
        incoming, outgoing = Path(temporary)/"input.txt", Path(temporary)/"basis.json"
        incoming.write_text("3\n"+"\n".join(map(str, entries)), encoding="ascii")
        outgoing.write_text(json.dumps(fixture), encoding="utf-8")
        report = verify(incoming, outgoing, root, config)
        if report["no_nonzero_integer_polynomial_in_budget_proved"]:
            raise AssertionError("错误地排除了已知关系2r-1=0")
        reports.append({"control": "known_lattice_relation_2r_minus_1", "passed": True})
        bad = deepcopy(fixture)
        bad["inverse_transform"][0][0] += 1
        outgoing.write_text(json.dumps(bad), encoding="utf-8")
        must_reject(lambda: verify(incoming, outgoing, root, config), "lattice_corrupted_inverse", reports)
        bad = deepcopy(fixture)
        bad["basis"][0][-1] += 1
        outgoing.write_text(json.dumps(bad), encoding="utf-8")
        must_reject(lambda: verify(incoming, outgoing, root, config), "lattice_row_outside_original_lattice", reports)
        compact = {k: v for k, v in fixture.items() if k != "inverse_transform"}
        outgoing.write_text(json.dumps(compact), encoding="utf-8")
        if verify(incoming, outgoing, root, config)["no_nonzero_integer_polynomial_in_budget_proved"]:
            raise AssertionError("紧凑接受端错误地排除了已知关系")
        reports.append({"control": "known_compact_lattice_relation", "passed": True})
        compact["basis"][0] = [2*c for c in compact["basis"][0]]
        outgoing.write_text(json.dumps(compact), encoding="utf-8")
        must_reject(lambda: verify(incoming, outgoing, root, config), "compact_lattice_nonunimodular_transform", reports)


def fiber_controls(reports):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    from verify_fiber import verify
    folder = HERE / "fiber-runs/20261015-p65521-degree8"
    metadata = read(folder / "field.json")
    field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
    data, system = read(folder / "fiber-0.json.gz"), read(AXIS)
    bad = deepcopy(data)
    bad["left0"][1] = deepcopy(bad["left0"][0])
    must_reject(lambda: verify(bad, field, system), "fiber_duplicate_left_root", reports)
    bad = deepcopy(data)
    bad["left0"][0]["point"][0]["gradient"][0][0] = (bad["left0"][0]["point"][0]["gradient"][0][0]+1)%metadata["prime"]
    must_reject(lambda: verify(bad, field, system), "fiber_corrupted_tangent", reports)
    bad = deepcopy(data)
    bad["norm_h3_mod_p"] = (bad["norm_h3_mod_p"]+1)%metadata["prime"]
    must_reject(lambda: verify(bad, field, system), "fiber_corrupted_norm", reports)


def chart_controls(reports):
    from critical_chart import verify
    certificate = read(HERE / "critical-runs/20261015-target/certificate.json.gz")
    root, system = read(TARGET / "refined-root.json"), read(AXIS)
    bad = deepcopy(certificate)
    bad["inverse_numerators"] = [["0"]*10 for _ in range(10)]
    must_reject(lambda: verify(bad, root, system), "chart_singular_preconditioner", reports)


def run(args):
    with bounded() as snapshot:
        tick = time.perf_counter()
        reports = []
        line_controls(reports)
        lattice_controls(reports)
        fiber_controls(reports)
        chart_controls(reports)
        write(HERE / args.output, {"verified": True, "controls": reports,
                                               "seconds": time.perf_counter()-tick, "limits": snapshot()})
        print({"verified": True, "control_count": len(reports), "seconds": time.perf_counter()-tick}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="controls-verified.json")
    run(parser.parse_args())
