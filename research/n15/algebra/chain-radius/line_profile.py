"""沿参数直线取完整范数的一阶数据，并用多项式模基认证有界有理拟合。"""

import argparse
from math import lcm
import time

from common import HERE, AXIS, bounded, read, write, sha
from chain_fiber import FACTOR, Jet, left_quartic, fiber


def splitting_degree(prime_field_square, values, factor):
    from flint import fq_default_poly_ctx
    params = [Jet(prime_field_square(v)) for v in values]
    u, v, w, r = params
    L = (v-u)**2+w*w
    eta = L*(4*r*r-L)
    Jet.sqrt_cache = {}
    degrees = []
    for t in eta.sqrt_pair():
        a, b, coefficients = left_quartic(params, t, factor["factor_polynomials"][factor["target_selected_factor_index"]])
        fp = fq_default_poly_ctx(prime_field_square)([c.value for c in coefficients])
        if fp.degree() != 4:
            raise ValueError("四次首项在此参数点退化")
        unit, factors = fp.factor()
        if any(m != 1 for f, m in factors):
            raise ValueError("四次载体在此参数点不分离")
        degrees.append([f.degree() for f, m in factors])
    initial = lcm(*(2*d for row in degrees for d in row))
    # 两次几何步各至多再开一个二次扩张；右端在八次扩域中分裂。
    return lcm(4*initial, 8), degrees


def row_degree(row):
    return max(p.degree() for p in row)


def leading_position(row):
    degree = row_degree(row)
    return max(i for i, p in enumerate(row) if p.degree() == degree)


def reduce_module(M, Y, prime):
    from flint import nmod_poly
    zero, one, X = nmod_poly([], prime), nmod_poly([1], prime), nmod_poly([0, 1], prime)
    basis = [[M, zero], [Y, one]]
    transform = [[one, zero], [zero, one]]
    operations = 0
    while leading_position(basis[0]) == leading_position(basis[1]):
        a, b = (0, 1) if row_degree(basis[0]) >= row_degree(basis[1]) else (1, 0)
        j = leading_position(basis[a])
        degree = basis[a][j].degree()-basis[b][j].degree()
        coefficient = int(basis[a][j].leading_coefficient())*pow(int(basis[b][j].leading_coefficient()), -1, prime)%prime
        multiplier = coefficient*X**degree
        basis[a] = [v-multiplier*w for v, w in zip(basis[a], basis[b])]
        transform[a] = [v-multiplier*w for v, w in zip(transform[a], transform[b])]
        operations += 1
    return basis, transform, operations


def encode_matrix(matrix):
    return [[list(map(int, p.coeffs())) for p in row] for row in matrix]


def module_certificate(M, Y, prime, budget, count):
    basis, transform, operations = reduce_module(M, Y, prime)
    degrees = list(map(row_degree, basis))
    return {"prime": prime, "sample_count": count, "degree_budget_inclusive": budget,
            "modulus": list(map(int, M.coeffs())), "hermite_interpolant": list(map(int, Y.coeffs())),
            "basis": encode_matrix(basis), "unimodular_transform": encode_matrix(transform),
            "row_degrees": degrees, "leading_positions": list(map(leading_position, basis)),
            "no_bounded_rational_function_matches_all_recorded_jets": min(degrees) > budget,
            "row_operations": operations,
            "does_not_bound_target_radius_minimal_degree": True}


def run(args):
    from flint import fq_default_ctx, nmod_poly
    out = HERE / "line-runs" / args.tag
    if out.exists():
        raise FileExistsError("拒绝覆盖已有直线实验")
    out.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        factor = read(FACTOR)
        small_field = fq_default_ctx(args.prime, 2, "b", fq_type="FQ_NMOD")
        fields = {}
        write(out / "budget.json", {"prime": args.prime, "fixed_u_v_w": [2, 3, 1],
                                     "maximum_degree_budget": args.degree, "maximum_wall_seconds": args.seconds,
                                     "source_system_sha256": sha(AXIS), "factor_sha256": sha(FACTOR),
                                     "driver_sha256": sha(__file__), "fiber_source_sha256": sha(HERE / "chain_fiber.py")})
        M, Y, X = nmod_poly([1], args.prime), nmod_poly([], args.prime), nmod_poly([0, 1], args.prime)
        samples, rejected, checkpoints = [], [], []
        next_budget = 16
        for radius in range(2, min(args.prime, 10000)):
            if time.perf_counter()-tick > args.seconds:
                break
            values = [2, 3, 1, radius]
            try:
                degree, pattern = splitting_degree(small_field, values, factor)
                if degree not in fields:
                    fields[degree] = fq_default_ctx(args.prime, degree, "a", fq_type="FQ_NMOD")
                    write(out / f"field-{degree}.json", {"degree": degree, "prime": args.prime,
                                                        "modulus_low_to_high": list(map(int, fields[degree].modulus().coeffs()))})
                data = fiber(fields[degree], values, factor, save_points=False)
            except (ValueError, ZeroDivisionError) as error:
                rejected.append({"r": radius, "reason": str(error)})
                continue
            value, derivative = data["norm_h3_mod_p"], data["norm_gradient_u_v_w_r_mod_p"][3]
            at, mv = Y(radius), M(radius)
            c0 = (value-int(at))*pow(int(mv), -1, args.prime)%args.prime
            updated = Y+c0*M
            c1 = (derivative-int(updated.derivative()(radius)))*pow(int(mv), -1, args.prime)%args.prime
            Y = updated+c1*M*(X-radius)
            M *= (X-radius)**2
            samples.append({"r": radius, "value": value, "derivative_r": derivative,
                            "all_derivatives_u_v_w_r": data["norm_gradient_u_v_w_r_mod_p"],
                            "splitting_degree": degree, "quartic_factor_degrees_over_Fp2": pattern,
                            "zero_bridge_pairs": data["vanishing_bridge_pairs"],
                            "evaluation_seconds": data["reconstruction_seconds"]+data["bridge_product_seconds"]})
            if len(samples) == next_budget+1:
                certificate = module_certificate(M, Y, args.prime, next_budget, len(samples))
                path = out / f"module-d{next_budget}.json.gz"
                write(path, certificate)
                checkpoints.append({k: certificate[k] for k in ("sample_count", "degree_budget_inclusive", "row_degrees", "no_bounded_rational_function_matches_all_recorded_jets")})
                write(out / f"samples-d{next_budget}.json.gz", samples)
                print(checkpoints[-1], flush=True)
                if next_budget >= args.degree:
                    break
                next_budget = min(2*next_budget, args.degree)
        write(out / "samples.json.gz", samples)
        write(out / "summary.json", {"sample_count": len(samples), "rejected_parameters": rejected,
                                      "checkpoints": checkpoints, "seconds": time.perf_counter()-tick,
                                      "limits": snapshot(), "radius_polynomial_obtained": False,
                                      "scope": "固定u,v,w的有限域范数函数次数探测；不推断实际半径的代数次数"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--degree", type=int, default=256)
    parser.add_argument("--seconds", type=float, default=600)
    run(parser.parse_args())
