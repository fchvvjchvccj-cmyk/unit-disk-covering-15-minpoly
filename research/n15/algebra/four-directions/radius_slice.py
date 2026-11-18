"""编译有限域截面求值器，并用独立的有理函数运算复核范数和重数。"""

import argparse
from fractions import Fraction as Q
import json
import os
from pathlib import Path
import subprocess
import tempfile

from chart import HERE
from context import ROOT, write
from compare import initialize
from eliminate import primitive_resultant
from verify import require

PRIME = 998244353


def verify(coefficients, fixed, vary_radius):
    from flint import nmod_poly
    _, rational = initialize(PRIME, 180)
    d = rational(nmod_poly([0, 1], PRIME))
    fixed_value = fixed.numerator*pow(fixed.denominator, -1, PRIME) % PRIME
    polynomial = nmod_poly(coefficients, PRIME)
    checks = []
    from sympy.ntheory import sqrt_mod
    for point in (3, *range(19, 100)):
        a, r = (fixed_value, point) if vary_radius else (point, fixed_value)
        s = r*(a-pow(a, -1, PRIME)) % PRIME
        root = sqrt_mod((s*s+4) % PRIME, PRIME)
        if not root:
            continue
        e = (s+int(root))*pow(2, -1, PRIME) % PRIME
        result = pow(a, 14904, PRIME)*(pow(r, 9936, PRIME) if vary_radius else 1) % PRIME
        for branch in (e, -pow(e, -1, PRIME) % PRIME):
            scale = rational(d.n**36*(d.n+r)**16)
            f = scale*primitive_resultant(rational(branch), rational(a), d)
            require(f.d.degree() == 0 and f.n.degree() <= 104, "截面清分母或次数错误")
            val = f.n.resultant(f.n.derivative()) if f.n.degree() == 104 else 0
            result = result*int(val) % PRIME
        require(result == polynomial(point), "独立结式值不符")
        checks.append([point, result])
        if len(checks) == 6:
            break
    require(len(checks) == 6 and any(v for _, v in checks), "独立检查点不足")
    content, layers = polynomial.factor_squarefree()
    rebuilt = nmod_poly([int(content)], PRIME)
    for f, m in layers:
        rebuilt *= f**m
    require(rebuilt == polynomial, "有限域截面回乘失败")
    layer = next(f for f, m in layers if m == 1)
    factors = layer.factor()[1]
    check = nmod_poly([1], PRIME)
    for factor, m in factors:
        check *= factor**m
    require(check == layer, "一重层因式回乘失败")
    return {"prime": PRIME, "fixed_parameter": "A" if vary_radius else "r", "fixed_value": str(fixed),
            "normalized_by_A_power": 14904, "normalized_by_r_power": 9936 if vary_radius else 0,
            "interpolation_points": 57961 if vary_radius else 34777, "withheld_points": 8, "independent_resultants": checks,
            "degree": polynomial.degree(), "squarefree_layers": [[f.degree(), m] for f, m in layers],
            "multiplicity_one_factor_degrees": sorted(f.degree() for f, _ in factors),
            "multiplicity_one_resultant": int(layer.resultant(layer.derivative())),
            "coefficients_ascending": list(map(int, polynomial.coeffs())),
            "radius_minimal_polynomial_obtained": False}


def run(fixed, vary_radius):
    chain = ROOT / ".tools/msolve-paged/chain/ucrt64/bin"
    env = {**os.environ, "PATH": str(chain)+os.pathsep+os.environ["PATH"]}
    with tempfile.TemporaryDirectory(dir=ROOT / ".tools", prefix="cover15-slice-") as temp:
        exe = Path(temp) / "slice.exe"
        subprocess.run([str(chain / "g++.exe"), "-O3", "-std=c++17", str(HERE / "radius_slice.cpp"),
                        "-lflint", "-o", str(exe)], env=env, check=True)
        residue = fixed.numerator*pow(fixed.denominator, -1, PRIME) % PRIME
        if residue in (0, 1, PRIME-1) or (not vary_radius and any((k+residue) % PRIME == 0 for k in range(1,106))):
            raise ValueError("截面参数退化或在固定取值节点上产生零分母")
        raw = subprocess.check_output([str(exe), str(residue), "r" if vary_radius else "A"], env=env)
    return verify(json.loads(raw), fixed, vary_radius)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--radius", type=Q, default=Q(9, 16))
    group.add_argument("--fixed-a", type=Q)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        report = run(args.fixed_a if args.fixed_a is not None else args.radius, args.fixed_a is not None)
        write(args.output, report)
        print(report["degree"], report["squarefree_layers"])
    finally:
        close()
