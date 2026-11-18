"""核验显式有理端点公式使用的系数、根箱、Newton 常数和共轭根间距界。"""

import argparse
from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import tempfile

from chart import HERE
from context import ROOT, sha, write
from verify import require

CORE = ROOT / "n15/proofs/global-final/replay/core"
sys.path.insert(0, str(CORE))
import system15 as system
import verify_root15


def certify():
    _, center, inverse, _ = system.load()
    # 前 71 式与后 77 式的系数直接取自同一组约束和驻定条件。
    first, last = [], [int(i == system.TID) for i in range(system.P)]
    for hessian, linear, constant in system.CONS:
        require(all(a % 2 == 0 for a in hessian.values()), "二次项缺少整数系数")
        first.append(abs(constant) + sum(map(abs, linear.values()))
                     + sum(map(abs, hessian.values())) // 2)
        for i, a in linear.items():
            last[i] += abs(a)
        for (i, _), a in hessian.items():
            last[i] += abs(a)
    l1 = max([2, *first, *last])
    y_norm = max(sum(map(abs, row)) for row in inverse)
    require(l1 <= 1024 and y_norm < 152, "实际方程系数或预条件器超出推导界")
    with tempfile.TemporaryDirectory() as scratch:
        previous = verify_root15.R
        try:
            verify_root15.R = Path(scratch)
            shutil.copyfile(CORE/"root15.json", Path(scratch)/"root15.json")
            root = verify_root15.verify()
        finally:
            verify_root15.R = previous
    outer, inner = root["boxes"]
    require(Q(outer["rho"]) == Q(1,10**60) and Q(outer["contraction"]) < Q(1,2)
            and Q(inner["rho"]) == Q(1,10**80), "两层原根箱或导数界不符")
    r0, upper = map(Q, root["r_interval"])
    require(Q(3,10) < r0 < upper < Q(1,3) and upper-r0 == Q(1,10**60),
            "正平方根初值的误差界不符")
    require(len(center) == 148 and system.TID == 76, "根坐标顺序不符")
    size = 3**149
    log_height = size * (size.bit_length() + 1 + 10*299)
    log_sep = (size*size + 2*size)*(log_height+2)
    exponent, steps = 10**223, 768
    require(log_height < 10**76 and log_sep < 10**222 and (1 << steps) > exponent+1,
            "根分离或迭代次数界不符")
    require(2*(1+304*l1) < 10**10, "联合 Newton 误差界不符")
    return {"system_sha256": sha(CORE/"system15.py"), "root_sha256": root["root_sha256"],
            "augmented_equations":149, "coefficient_l1_upper_actual":l1,
            "preconditioner_infinity_norm":str(y_norm), "large_box_contraction":outer["contraction"],
            "positive_radius_initial":str(r0), "squared_radius_index":system.TID,
            "deformation_algebra_rank":str(size),"log2_annihilator_height_upper":str(log_height),
            "negative_log2_separation_upper":str(log_sep), "newton_steps":steps,
            "endpoints":"r_768 - 2^(-10^223), r_768 + 2^(-10^223)",
            "newton_error_constant":10**10, "rational_endpoint_formula_proved":True,
            "endpoints_expanded":False, "minimal_polynomial_coefficients_computed":False}


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _, close=apply_limits()
    try:
        write(args.output,certify())
    finally:
        close()
