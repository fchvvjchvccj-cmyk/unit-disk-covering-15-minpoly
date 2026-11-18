"""在小素域逐层枚举六个二次方程，寻找可回代原系统的光滑临界种子。"""

import argparse
from collections import Counter
import time

from support import HERE, AXIS, FACTOR, read, write, sha, bounded
from chain_fiber import Jet
from geometry import dot, boundary
from verify_system import execute
from target_chart import Dual


def square_roots(value):
    if value.value.is_zero() or not value.value.is_square():
        return []
    return value.sqrt_pair()


def scalar(value):
    coefficients = list(map(int, value.to_list()))
    if any(coefficients[1:]):
        raise ValueError("小素域点没有降到基域")
    return coefficients[0] if coefficients else 0


def left_points(params, stats):
    u, v, w, r = params
    dx, L = v-u, (v-u)**2+w*w
    for t in square_roots(L*(4*r*r-L)):
        try:
            aa, bb = dx-w*t/L, w+dx*t/L
            s, px, py = bb/(2*r+aa), u+aa/2, bb/2
        except ZeroDivisionError:
            stats["left_initial_pole"] += 1
            continue
        for z in square_roots(1-py*py):
            ax = (px+z)/2
            for offset in square_roots(r*r-(ax-px)**2):
                h = py+offset
                try:
                    a = [ax, h]
                    T, K = dot(a, a), 1+dot(a, a)-r*r
                    c3 = [v+K*ax/T-z-px, w+K*h/T-2*py]
                    TT = dot(c3, c3)
                    c4 = [v+(1-r*r)*c3[0]/TT-ax, w+(1-r*r)*c3[1]/TT-h]
                    # 保持旧四次图表的逆映射和边界交点公式有效。
                    for guard in (ax-u, h, 2*ax-u):
                        if guard.value.is_zero():
                            raise ZeroDivisionError("旧左端图表主元为零")
                    yield {"ports": [s, z, h], "c0": [ax, -h], "c3": c3, "c4": c4}
                except ZeroDivisionError:
                    stats["left_reconstruction_pole"] += 1


def right_points(params, stats):
    u, v, w, r = params
    E, beta = (u-v)**2+w*w, r*r-1
    try:
        roots = square_roots(4*w*w*r*r/E-w*w)
    except ZeroDivisionError:
        stats["right_initial_pole"] += 1
        return
    for root in roots:
        x = (1+r+v+root)/2
        for y in square_roots(r*r-(x-1)**2):
            try:
                Sx = 2*x+beta
                ac, ec = 2-4*x*x/Sx, -4*x*y/Sx
                Z, Pz = 2*x-1-r, (v-(2*x-1-r))*(v-u)
                bc, be, bc0 = -(Z+2*v-u), -(w*w-Pz)/w, Z*(2*v-u)
                ld, le, l0 = bc-ac, be-ec, bc0+beta
                e0, e1 = -l0/le, -ld/le
                A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
                if A.value.is_zero():
                    raise ZeroDivisionError("右端二次首项为零")
                for disc in square_roots(B*B-4*A*C):
                    d = (-B+disc)/(2*A)
                    e = e0+e1*d
                    yield {"ports": [x, y, d], "c5": [d, e], "c6": [x, y]}
            except ZeroDivisionError:
                stats["right_reconstruction_pole"] += 1


def eliminate(matrix, right, prime):
    rows = [[int(v)%prime for v in row]+[int(b)%prime] for row, b in zip(matrix, right)]
    columns, pivot_row, pivots = len(matrix[0]), 0, []
    for column in range(columns):
        pivot = next((i for i in range(pivot_row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [v*inverse%prime for v in rows[pivot_row]]
        for i in range(len(rows)):
            if i != pivot_row and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [(a-factor*b)%prime for a, b in zip(rows[i], rows[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    if any(all(v == 0 for v in row[:columns]) and row[-1] for row in rows):
        raise ValueError("模素数线性方程不相容")
    if len(pivots) != columns:
        raise ValueError("模素数矩阵没有满列秩")
    solution = [0]*columns
    for i, j in enumerate(pivots):
        solution[j] = rows[i][-1]
    return solution


def row_value(row, point, prime):
    return sum(term["coefficient"]*product(point, term["powers"], prime) for term in row)%prime


def product(point, powers, prime):
    result = 1
    for j, e in powers:
        result = result*pow(point[j], e, prime)%prime
    return result


def row_derivative(row, point, column, prime):
    result = 0
    for term in row:
        exponent = dict(term["powers"]).get(column, 0)
        if exponent:
            powers = [[j, e-int(j == column)] for j, e in term["powers"] if e-int(j == column)]
            result += term["coefficient"]*exponent*product(point, powers, prime)
    return result%prime


def verify_seed(seed, data, system, field):
    prime, ports = seed["prime"], list(map(field, seed["ports"]))
    dimension = len(ports)
    registers = execute(data, [Dual(v, [field(int(i == j)) for j in range(dimension)]) for i, v in enumerate(ports)])
    equations = data["outputs"]["finite_relations"]+data["outputs"]["critical_equations"]
    if any(not registers[n].value.is_zero() for n in equations):
        raise ValueError("种子不是十条临界方程的公共零点")
    for n in data["all_divisor_nodes"]+data["tower_leading_coefficient_nodes"]+data["simple_root_derivative_nodes"]+data["critical_locus_nonzero_guards"]:
        if registers[n].value.is_zero():
            raise ValueError("种子落在局部图表的坏集")
    J = [[scalar(v) for v in registers[n].gradient] for n in equations]
    eliminate(J, [0]*dimension, prime)
    geometry = []
    for name in ("c0", "c3", "c4", "c5", "c6"):
        geometry += [scalar(registers[n].value) for n in data["outputs"]["geometry"][name]]
    geometry += seed["ports"][-4:]
    if any(row_value(row, geometry, prime) for row in system["equations"][:11]):
        raise ValueError("临界种子没有满足原十一条几何方程")
    geometry_jacobian = [[row_derivative(row, geometry, column, prime) for column in range(14)]
                         for row in system["equations"][:11]]
    theta = eliminate([[geometry_jacobian[i][j] for i in range(11)] for j in range(14)],
                      [-int(j == 13) for j in range(14)], prime)
    point = geometry+theta
    if any(row_value(row, point, prime) for row in system["equations"]):
        raise ValueError("回代乘子没有满足原25条整数方程")
    full_jacobian = [[row_derivative(row, point, j, prime) for j in range(25)] for row in system["equations"]]
    eliminate(full_jacobian, [0]*25, prime)
    return {"verified": True, "critical_jacobian_rank": dimension, "original_geometric_jacobian_rank": 11,
            "original_kkt_jacobian_rank": 25, "original_kkt_coordinates": point,
            "hensel_lift_unique_in_residue_class": True,
            "identified_with_certified_real_target": False}


def scan(prime, data, system, deadline):
    from flint import fq_default_ctx
    field = fq_default_ctx(prime, 1, "a", fq_type="FQ_NMOD")
    stats, seeds, failed = Counter(), [], []
    tick = time.perf_counter()
    Jet.sqrt_cache, Jet.sqrt_calls, Jet.sqrt_hits = {}, 0, 0
    for radius in range(2, prime-1):
        for uu in range(prime):
            for vv in range(prime):
                for ww in range(1, prime):
                    if time.perf_counter() > deadline:
                        return {"prime": prime, "complete_parameter_scan": False, "counts": dict(stats),
                                "seeds": seeds, "rejected_candidates": failed, "seconds": time.perf_counter()-tick}
                    stats["parameter_points"] += 1
                    values = list(map(field, (uu, vv, ww, radius)))
                    params = [Jet(v, [field(int(i == j)) for j in range(4)]) for i, v in enumerate(values)]
                    try:
                        left, right = list(left_points(params, stats)), list(right_points(params, stats))
                    except (ZeroDivisionError, ValueError):
                        stats["parameter_degeneracy"] += 1
                        continue
                    stats["left_branches"] += len(left)
                    stats["right_branches"] += len(right)
                    for a in left:
                        for b in right:
                            stats["bridge_pairs"] += 1
                            H = boundary(a["c4"], b["c5"], params[3])
                            if not H.value.is_zero():
                                continue
                            stats["closed_pairs"] += 1
                            if any(not value.is_zero() for value in H.gradient[:3]):
                                continue
                            stats["critical_candidates"] += 1
                            seed = {"prime": prime, "ports": [scalar(v.value) for v in a["ports"]+b["ports"]+params],
                                    "radius": radius}
                            try:
                                certificate = verify_seed(seed, data, system, field)
                                seed.update(certificate)
                                seeds.append(seed)
                            except (ValueError, ZeroDivisionError) as error:
                                failed.append({**seed, "reason": str(error)})
    return {"prime": prime, "complete_parameter_scan": True, "counts": dict(stats),
            "seeds": seeds, "rejected_candidates": failed, "seconds": time.perf_counter()-tick}


def run(args):
    target = HERE / "seed-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有模素数种子搜索")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        source = HERE / "system-runs/20261015-critical10/circuit.json.gz"
        data, system = read(source), read(AXIS)
        reports = []
        for prime in args.primes:
            report = scan(prime, data, system, time.perf_counter()+args.seconds)
            reports.append(report)
            write(target / f"p{prime}.json", report)
            print({"prime": prime, "completed": report["complete_parameter_scan"], "seeds": len(report["seeds"]),
                   "counts": report["counts"], "seconds": report["seconds"]}, flush=True)
        write(target / "summary.json", {"reports": [{"prime": r["prime"], "verified_seeds": len(r["seeds"]),
                                                     "complete_parameter_scan": r["complete_parameter_scan"]} for r in reports],
                                        "source_system_sha256": sha(AXIS), "circuit_sha256": sha(source),
                                        "driver_sha256": sha(__file__), "seconds": time.perf_counter()-tick, "limits": snapshot(),
                                        "scope": "只搜索r不等于0或正负1、w非零且所用二次层和旧左端图表正则的Fp有理点；不声称代数闭包内完整解集或目标半径分量。"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--primes", type=int, nargs="+", default=[7, 11])
    parser.add_argument("--seconds", type=int, default=60)
    run(parser.parse_args())
