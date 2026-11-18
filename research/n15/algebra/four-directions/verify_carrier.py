"""逐项重建特征零系数表，并以三个有理方向对独立核对全部 105 个方向系数。"""

import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from itertools import zip_longest

from chart import HERE, read
from binary_projection import polynomial_value
from context import write
from verify import require

TABLE = HERE / "carrier.csv.gz"


def check_rebuild():
    from flint import fmpz_mpoly_ctx
    from expand_radius import equations
    from structured_expansion import carrier_components
    ring = fmpz_mpoly_ctx.get(("r","A","E","D","L"))
    h, branch = equations(*ring.gens())
    fitted, exponent = carrier_components(h, branch, lambda *_: None)
    expected = ((j,*map(int,m),int(c)) for j,p in enumerate((fitted.a,fitted.b)) for m,c in p.terms())
    with gzip.open(TABLE, "rt", encoding="ascii") as stream:
        require(json.loads(next(stream))["A_exponent"] == exponent, "坐标幂次不一致")
        actual = (tuple(map(int,line.split(","))) for line in stream)
        require(all(x == y for x,y in zip_longest(expected,actual)), "重建的全参数系数表不一致")


def verify():
    from flint import nmod_poly
    from sympy import QQ
    from sympy.polys.fields import field
    ring, d = field("D", QQ)
    cases = [(Q(2),Q(3)), (Q(3),Q(2)), (Q(5,2),Q(7,5))]
    data = []
    for e,a in cases:
        r = (e-1/e)/(a-1/a)
        def powers(q, maximum):
            return [q.numerator**i*q.denominator**(maximum-i) for i in range(maximum+1)]
        data.append((powers(r,131), powers(a,78), r.denominator**131*a.denominator**78, [[0]*105 for _ in range(2)]))
    counts = [0,0]
    with gzip.open(TABLE,"rt",encoding="ascii") as stream:
        header = json.loads(next(stream))
        require(header["variables"] == ["r","A","D"] and header["basis"] == ["1","A*E"], "表布局错误")
        for line in stream:
            j, ir, ia, id_, c = map(int, line.split(","))
            require(j in (0,1) and 0 <= ir <= 131 and 0 <= ia <= 78 and 0 <= id_ <= 104,
                    "表中指数越界")
            counts[j] += 1
            for pr,pa,_,sums in data:
                sums[j][id_] += c*pr[ir]*pa[ia]
    for (e,a), (_,_,denom,sums) in zip(cases, data):
        r = (e-1/e)/(a-1/a)
        f = polynomial_value(ring(e),ring(a),d)
        require(f.denom.is_ground and f.numer.degree() == 104, "独立特化次数错误")
        normalization = a**header["A_exponent"]*r**header["r_exponent"]
        for i in range(105):
            actual = (sums[0][i]+a*e*sums[1][i])/denom
            expected = normalization*Q(str(f.numer.get((i,),0)))/Q(str(f.denom[(0,)]))
            require(actual == expected, "与独立特征零方向多项式的系数不符")
    with TABLE.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    prime = 1093
    require(all(prime%i for i in range(2,34)), "不可约性检验需要素数")
    coefficients = list(map(int,read(HERE/"binary.json")["specialized_primitive_coefficients_ascending"]))
    f = nmod_poly(coefficients,prime)
    x = nmod_poly([0,1],prime)
    require(f.degree()==104, "模素数后截面丢失次数")
    power = x
    for j in range(1,105):
        power = power.pow_mod(prime,f)
        if j in (8,52):
            require(f.gcd(power-x).degree()==0, "截面未通过 Rabin 互素条件")
    require(power==x, "截面未通过 Frobenius 恒等式")
    return {"terms_per_component": counts, "sha256": digest, "encoding": header,
            "exact_specializations_E_A": [[str(v) for v in c] for c in cases],
            "coefficients_compared_per_specialization": 105,
            "direction_polynomial_irreducibility":{"specialization_E_A":[2,3],"degree":104,
                 "prime":1093,"coprime_exponents":[8,52],"frobenius_exponent":104}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    from windows_compute_limits import apply_limits
    _, close = apply_limits()
    try:
        if args.rebuild:
            check_rebuild()
        result = verify()
        result["all_coefficients_rebuilt"] = args.rebuild
        write(args.output, result)
        print(result["terms_per_component"])
    finally:
        close()
