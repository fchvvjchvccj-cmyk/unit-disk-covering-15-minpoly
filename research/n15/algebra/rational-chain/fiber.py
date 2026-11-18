"""在秩六十四的已选有理链载体上计算闭环范数及参数导数。"""

import argparse
import json
import time

from support import HERE, AXIS, FACTOR, read, write, sha, bounded
from geometry import reconstruct, boundary
from chain_fiber import Jet, left_endpoint, right_endpoint, prime_scalar


def encode(root):
    return {key: [v.encode() for v in value] if isinstance(value, list)
            else value.encode() if isinstance(value, Jet) else value for key, value in root.items()}


def compute(field, values, factor, save_points=True):
    Jet.sqrt_cache = {}
    Jet.sqrt_calls = Jet.sqrt_hits = 0
    params = [Jet(field(v), [field(int(i == j)) for j in range(4)]) for i, v in enumerate(values)]
    u, v, w, r = params
    tick = time.perf_counter()
    left0 = left_endpoint(params, factor, field)
    left = [{**root, **reconstruct(root["point"], u, [v, w], r)} for root in left0]
    right = right_endpoint(params)
    if len(left) != 8 or len(right) != 8:
        raise ValueError("有理链有限载体未覆盖预期分支")
    reconstructed = time.perf_counter()
    product = r.cast(1)
    bridge_values, zeros = [], []
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            value = boundary(a["c4"], b["point"], r)
            product *= value
            if value.value.is_zero():
                zeros.append([i, j])
            if save_points:
                bridge_values.append(value.encode())
    result = {"parameters_u_v_w_r": values, "rank_before_closure": 64,
              "norm_h3_mod_p": prime_scalar(product.value),
              "norm_gradient_u_v_w_r_mod_p": [prime_scalar(a) for a in product.gradient],
              "vanishing_bridge_pairs": zeros, "reconstruction_seconds": reconstructed-tick,
              "bridge_product_seconds": time.perf_counter()-reconstructed,
              "native_square_roots": Jet.sqrt_calls}
    if save_points:
        result.update(left=[encode(a) for a in left], right=[encode(a) for a in right], bridge_values=bridge_values)
    return result


def run(args):
    from flint import fq_default_ctx, fmpz_mod_poly_ctx
    source = HERE / "fiber-runs" / args.tag
    if source.exists():
        raise FileExistsError("拒绝覆盖已有有理链纤维")
    source.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        reference = None
        if args.reference:
            reference = HERE.parent / "chain-radius/fiber-runs" / args.reference
            metadata = read(reference / "field.json")
            field = fq_default_ctx(modulus=fmpz_mod_poly_ctx(metadata["prime"])(metadata["modulus_low_to_high"]), fq_type="FQ_NMOD")
            prime, degree = metadata["prime"], metadata["degree"]
        else:
            prime, degree = args.prime, args.degree
            field = fq_default_ctx(prime, degree, "a", fq_type="FQ_NMOD")
        write(source / "field.json", {"prime": prime, "degree": degree,
                                     "modulus_low_to_high": list(map(int, field.modulus().coeffs())),
                                     "source_system_sha256": sha(AXIS), "left_factor_sha256": sha(FACTOR),
                                     "generator_sha256": sha(__file__), "geometry_sha256": sha(HERE / "geometry.py"),
                                     "old_fiber_reference": str(reference.relative_to(HERE.parent)) if reference else None})
        records, failures = [], []
        for i, parameters in enumerate(([2, 3, 1, 2], [3, 5, 2, 3], [2, 5, 3, 1], [7, 11, 3, 2])):
            try:
                data = compute(field, parameters, read(FACTOR))
            except (ValueError, ZeroDivisionError) as error:
                failures.append({"parameters": parameters, "reason": str(error)})
                continue
            write(source / f"fiber-{i}.json.gz", data)
            record = {k: v for k, v in data.items() if k not in ("left", "right", "bridge_values")}
            records.append(record)
            print(json.dumps(record), flush=True)
            if len(records) == 2:
                break
        write(source / "summary.json", {"fibers": records, "rejected": failures,
                                      "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                      "radius_minimal_polynomial_obtained": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--degree", type=int, default=24)
    parser.add_argument("--reference")
    run(parser.parse_args())
