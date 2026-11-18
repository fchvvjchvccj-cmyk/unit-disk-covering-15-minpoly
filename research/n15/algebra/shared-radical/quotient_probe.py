"""在Fp(r)上保留32维二次塔，直接计算闭环规范形、全导数和范数。"""

import argparse
from pathlib import Path
import subprocess
import sys
import time

from bridge_support import HERE, PARENT, AXIS, read, write, sha, bounded, binary_environment
from chain_fiber import Jet
from geometry import dot, boundary
from tower import Algebra, Element, Rational, TowerStop


def model(algebra, values, derivatives, checkpoint):
    params = [Jet(algebra.scalar(value), [algebra.scalar(int(derivatives and i == j)) for j in range(4)])
              for i, value in enumerate(values)]
    u, v, w, r = params

    def extend(a, b, name):
        g = algebra.add_level(a.value, b.value)
        partial = [c+d*g for c, d in zip(a.gradient, b.gradient)]
        if all(value.is_zero() for value in partial):
            root = Jet(g)
        else:
            inverse = (2*g-b.value).inverse()
            root = Jet(g, [value*inverse for value in partial])
        checkpoint("root_"+name)
        return root

    dx, L = v-u, (v-u)**2+w*w
    s = extend(-(L-2*r*dx)/(L+2*r*dx), 4*r*w/(L+2*r*dx), "s")
    aa, bb = 2*r*(1-s*s)/(1+s*s), 4*r*s/(1+s*s)
    px, py = u+aa/2, bb/2
    z = extend(1-py*py, py.cast(0), "z")
    ax = (px+z)/2
    h = extend(r*r-(ax-px)**2-py*py, 2*py, "h")
    a = [ax, h]
    T, K = dot(a, a), 1+dot(a, a)-r*r
    c3 = [v+K*ax/T-z-px, w+K*h/T-2*py]
    TT = dot(c3, c3)
    c4 = [v+(1-r*r)*c3[0]/TT-ax, w+(1-r*r)*c3[1]/TT-h]
    checkpoint("left_c4")
    x = (1+r+v+dx-aa)/2
    y = extend(r*r-(x-1)**2, x.cast(0), "y")
    beta, Sx = r*r-1, 2*x+r*r-1
    ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
    Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
    bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
    ld, le, l0 = bc-ac, be-ec, bc0+beta
    e0, e1 = -l0/le, -ld/le
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    d = extend(-C/A, -B/A, "d")
    e = e0+e1*d
    H = boundary(c4, [d, e], r)
    checkpoint("closure_reduced")
    return H, {"c0": [ax, -h], "c3": c3, "c4": c4, "c5": [d, e], "c6": [x, y]}


def worker(args):
    from flint import nmod_poly
    target = Path(args.worker)
    tick = time.perf_counter()
    Rational.prime = args.prime
    Rational.zero_poly, Rational.one_poly = nmod_poly([], args.prime), nmod_poly([1], args.prime)
    Rational.deadline, Rational.degree_limit = tick+args.seconds, args.max_degree
    Rational.metrics = {"max_coefficient_degree": 0, "rational_additions": 0, "rational_multiplications": 0}
    algebra = Algebra()
    events = []

    def checkpoint(stage):
        event = {"stage": stage, "seconds": time.perf_counter()-tick, "levels": len(algebra.relations),
                 "recursive_products": algebra.products, "recursive_inversions": algebra.inversions, **Rational.metrics}
        events.append(event)
        print(event, flush=True)

    report = {"prime": args.prime, "fixed_u_v_w": args.fixed, "derivatives_computed": args.derivatives,
              "coefficient_field": "Fp(r)", "tower_variables": ["s", "z", "h", "y", "d"], "rank": 32,
              "source_system_sha256": sha(AXIS), "driver_sha256": sha(__file__), "algebra_sha256": sha(HERE / "tower.py")}
    try:
        values = [Rational(value) for value in args.fixed]+[Rational(nmod_poly([0, 1], args.prime))]
        H, geometry = model(algebra, values, args.derivatives, checkpoint)
        data = {**report, "relations": [{"a": [value.encode() for value in a], "b": [value.encode() for value in b]}
                                        for a, b in algebra.relations],
                "H": H.value.encode(5),
                "geometry": {name: [point.value.encode(5) for point in coordinates] for name, coordinates in geometry.items()}}
        if args.derivatives:
            data["gradient_u_v_w_r"] = [g.encode(5) for g in H.gradient]
        write(target / "quotient.json.gz", data)
        norm = algebra.norm(H.value.promote(5).coefficients, 5)
        checkpoint("norm_reduced")
        unit, factors = norm.n.factor()
        denominator_unit, denominator_factors = norm.d.factor()
        norm_data = {"value": norm.encode(), "numerator_degree": int(norm.n.degree()), "denominator_degree": int(norm.d.degree()),
                     "numerator_unit": int(unit), "numerator_factors": [{"coefficients": list(map(int, f.coeffs())), "degree": int(f.degree()), "multiplicity": int(m)} for f, m in factors],
                     "denominator_unit": int(denominator_unit), "denominator_factors": [{"coefficients": list(map(int, f.coeffs())), "degree": int(f.degree()), "multiplicity": int(m)} for f, m in denominator_factors],
                     "radius_minimal_polynomial_obtained": False,
                     "scope": "固定u,v,w的Fp(r)闭环范数，保留r为符号变量；不是实际目标半径的消元多项式。"}
        write(target / "norm.json.gz", norm_data)
        report.update(status="completed", norm_numerator_degree=int(norm.n.degree()), norm_denominator_degree=int(norm.d.degree()))
    except (TowerStop, ZeroDivisionError) as error:
        report.update(status="stopped", reason=str(error))
    checkpoint("finished")
    report.update(events=events, seconds=time.perf_counter()-tick)
    write(target / "worker.json", report)


def run(args):
    target = HERE / "quotient-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有二次塔探测")
    target.mkdir(parents=True)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(target), "--prime", str(args.prime),
               "--seconds", str(args.seconds), "--max-degree", str(args.max_degree), "--fixed", *map(str, args.fixed)]
    if args.derivatives:
        command.append("--derivatives")
    with bounded() as snapshot:
        tick = time.perf_counter()
        with (target / "worker.log").open("x", encoding="utf-8") as stream:
            try:
                child = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=binary_environment(),
                                       timeout=args.seconds+5, check=False)
                status, returncode = "worker_completed", child.returncode
            except subprocess.TimeoutExpired:
                status, returncode = "hard_timeout", None
        report = {"controller_status": status, "returncode": returncode, "seconds": time.perf_counter()-tick,
                  "limits": snapshot(), "log_sha256": sha(target / "worker.log")}
        if (target / "worker.json").exists():
            report["worker"] = read(target / "worker.json")
        write(target / "controller.json", report)
        print({k: v for k, v in report.items() if k not in ("worker", "limits")}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag")
    parser.add_argument("--prime", type=int, default=65521)
    parser.add_argument("--fixed", type=int, nargs=3, default=[2, 3, 1])
    parser.add_argument("--seconds", type=int, default=90)
    parser.add_argument("--max-degree", type=int, default=4000)
    parser.add_argument("--derivatives", action="store_true")
    parser.add_argument("--worker")
    arguments = parser.parse_args()
    worker(arguments) if arguments.worker else run(arguments)
