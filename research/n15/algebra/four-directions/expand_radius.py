"""尝试在特征零展开消元式，记录实际到达的阶段和资源用量。"""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from chart import HERE
from context import write
from windows_compute_limits import apply_limits


def equations(r, a, e, d, left):
    u = 1+r*(a-1)
    v = a+r*(1-a)
    f = (u+r*r-1)*d*d+r*(1+u)*d+r*r
    g = a*r*r*d*d+r*(a+v)*d+v+a*(r*r-1)
    w = r*d*d+u*d+r
    wbar = a*r*d*d+v*d+a*r
    h0 = -a*d*w*f*g-a*a*r*r*f*f
    h1 = (w*wbar+(1-2*r*r)*d*d*a)*f*g
    h2 = -d*wbar*f*g-r*r*d**4*g*g
    anchor = h0+left*(h1+left*h2)
    tau = r*d
    t = d*(e-1)-r*a*d-r*(d*d+1-d)
    unit = 1-r*r
    pa, pb = -r*r*t*e, left*r*r*t
    pc, pd = -tau*left*r*unit*e, tau*r*unit
    k = left*(tau*unit-r*t*e)-tau*unit*e+r*t
    x, y, z, w = tau*(pb+pc), tau*(pb+pd), tau*(pa-pb), tau*k+pb*t
    c0, c3 = z*z, (x-y)**2
    c1 = t*w*z+tau*(w*w-4*x*y-2*(x+y)*z+z*z)
    c2 = t*t*x*y-tau*t*w*(x+y)+tau*tau*(w*w+c3-4*x*y-2*(x+y)*z)
    pencil = (c1*c1*c2*c2-4*tau**3*c3*c1**3-4*c0*c2**3-
              27*tau**6*c0*c0*c3*c3+18*tau**3*c0*c1*c2*c3)
    divisor = tau**10
    reduced = pencil//divisor
    if reduced*divisor != pencil:
        raise ValueError("铅笔清分母后的精确除法失败")
    return anchor, reduced


def note(stage, polynomial=None):
    record = {"stage": stage}
    if polynomial is not None:
        record.update(terms=len(polynomial), degrees=list(map(int, polynomial.degrees())))
    print(json.dumps(record), flush=True)


def expand(output, checkpoint):
    from flint import fmpz_mpoly_ctx
    ring = fmpz_mpoly_ctx.get(("r", "A", "E", "D", "L"))
    r, a, e, d, left = ring.gens()
    note("build_integer_equations")
    h, branch = equations(r,a,e,d,left)
    if checkpoint:
        from structured_expansion import carrier_components, save_components
        note("reduce_E_before_L", branch)
        pair, exponent = carrier_components(h, branch, note)
        del branch, h
        print(json.dumps(save_components(checkpoint, pair, exponent)), flush=True)
        note("eliminate_E_by_norm", pair.a)
        norm = pair.norm()
    else:
        note("eliminate_L", branch)
        carrier = h.resultant(branch, "L")
        f = carrier//d**4
        if f*d**4 != carrier or f.degrees()[3] != 104:
            raise ValueError("第一结式的阶数或清分母不符")
        note("eliminate_E", f)
        norm = (a*e*e-r*(a*a-1)*e-a).resultant(f, "E")
    note("eliminate_D", norm)
    plane = norm.resultant(norm.derivative("D"), "D")
    if plane.is_zero():
        raise ValueError("投影多项式恒零")
    note("plane_squarefree", plane)
    _, layers = plane.factor_squarefree()
    simple = next(p for p, multiplicity in layers if multiplicity == 1)
    # 对 A 的内容在实际根附近是单位；保留在结式中不妨碍其消失。
    note("eliminate_A", simple)
    relation = simple.resultant(simple.derivative("A"), "A")
    if relation.is_constant() or any(relation.degrees()[1:]):
        raise ValueError("未得到半径单变量关系")
    coefficients = [0]*(int(relation.degrees()[0])+1)
    for monomial, coefficient in relation.terms():
        coefficients[int(monomial[0])] = int(coefficient)
    write(output, {"variable": "r", "coefficients_ascending": coefficients, "irreducibility_proved": False})
    note("radius_relation", relation)


def run(seconds, output, log, checkpoint):
    if output.exists() or log.exists():
        raise FileExistsError("目标关系文件或试算记录已存在")
    snapshot, close = apply_limits()
    started = time.perf_counter()
    try:
        cmd = [sys.executable, "-X", "utf8", "-B", __file__, "--worker", "--output", str(output)]
        if checkpoint:
            cmd.extend(["--checkpoint", str(checkpoint)])
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as child:
            try:
                raw, _ = child.communicate(timeout=seconds)
                status = "completed" if child.returncode == 0 else "failed"
            except subprocess.TimeoutExpired:
                child.kill()
                raw, _ = child.communicate()
                status = "time_limit"
        report = {"status": status, "exit_code": child.returncode, "wall_seconds": time.perf_counter()-started,
                  "time_limit_seconds": seconds, "limits": snapshot(), "output_created": output.is_file(),
                  "observations": raw.decode("utf-8", errors="replace").splitlines()}
        write(log, report)
        print(json.dumps(report, ensure_ascii=False))
    finally:
        close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--seconds", type=int, default=180)
    parser.add_argument("--checkpoint", type=Path)
    args = parser.parse_args()
    if args.worker:
        expand(args.output, args.checkpoint)
    elif args.log:
        run(args.seconds, args.output, args.log, args.checkpoint)
    else:
        parser.error("需要指定 --log")
