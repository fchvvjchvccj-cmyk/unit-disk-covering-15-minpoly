"""精确证明左右二次根共享，并在原认证根盒上选择目标的同号分量。"""

import argparse
from fractions import Fraction as F
import time

from bridge_support import HERE, PARENT, TARGET, AXIS, PARENT_SYSTEM, read, write, sha, bounded


def identity():
    import sympy as S
    s, dx, w, r = S.symbols("s dx w r")
    A, L = 1+s*s, dx*dx+w*w
    aa, bb = 2*r*(1-s*s)/A, 4*r*s/A
    g0 = A*L-2*r*dx*(1-s*s)-4*r*w*s
    other = A*L-2*r*dx*(1-s*s)+4*r*w*s
    difference = S.together(A*A*(L*(dx-aa)**2-w*w*(4*r*r-L))-g0*other)
    if not S.Poly(difference.as_numer_denom()[0], s, dx, w, r).is_zero:
        raise ValueError("共享判别式的整数恒等式失败")
    return {"verified": True,
            "identity": "(1+s^2)^2*[L*(dx-aa)^2-w^2*(4*r^2-L)]=G0*[L*(1+s^2)-2*r*dx*(1-s^2)+4*r*w*s]",
            "definitions": {"dx": "v-u", "L": "dx^2+w^2", "aa": "2*r*(1-s^2)/(1+s^2)"}}


def verify(root):
    from target_chart import image_box, dependencies
    data = read(PARENT_SYSTEM)
    parents = dependencies(data)
    arithmetic, ports = image_box(root, 1024, data["variables"])
    s, _, _, x, _, _, u, v, w, r = ports
    aa = 2*r*(1-s*s)/(1+s*s)
    C, R = v-u-aa, 2*x-1-r-v
    other, selected = (R+C).value, (R-C).value
    if other[0] <= 0 <= other[1]:
        raise ValueError("不能排除共享根的另一符号分量")
    if not selected[0] <= 0 <= selected[1]:
        raise ValueError("目标根与同号分量不相容")
    encode = lambda interval: [str(F(a, arithmetic.scale)) for a in interval]
    return {"verified": True, "integer_identity": identity(), "parent_proofs": parents,
            "opposite_factor_R_plus_C_interval": encode(other),
            "selected_factor_R_minus_C_interval_consistency_only": encode(selected),
            "selected_sign": 1, "c6x_formula": "(1+r+v+(v-u-2*r*(1-s^2)/(1+s^2)))/2",
            "generic_selected_rank_before_closure": 32,
            "argument": "原G0为零使右端二次式精确分解为L*(R-C)*(R+C)。L和另一因子在认证根盒非零，故R=C，c6x可有理恢复。精确零来自因子恒等式和原根方程，区间包含零仅作相容检查。",
            "scope": "选择原64维载体中包含目标根的32维共享根分量，不是半径的代数次数或最小多项式。"}


def run(args):
    target = HERE / "proof-runs" / args.tag
    if target.exists():
        raise FileExistsError("拒绝覆盖已有共享根证明")
    target.mkdir(parents=True)
    with bounded() as snapshot:
        tick = time.perf_counter()
        report = verify(read(TARGET / "refined-root.json"))
        report.update(source_system_sha256=sha(AXIS), root_sha256=sha(TARGET / "refined-root.json"),
                      parent_circuit_sha256=sha(PARENT_SYSTEM), verifier_sha256=sha(__file__),
                      seconds=time.perf_counter()-tick, limits=snapshot())
        write(target / "verified.json.gz", report)
        print({"verified": True, "rank": 32, "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    run(parser.parse_args())
