"""对重建的整数候选独立检查本原性、模因子次数交集以及小区间的实根数。"""

import argparse
from fractions import Fraction as Q
from functools import reduce
from math import gcd
from pathlib import Path

from chart import HERE,read
from context import write
from sample_plane import digest
from verify import require


def degrees_across_primes(data):
    from flint import nmod_poly,fmpz_mod_poly_ctx
    from sympy import nextprime
    degree=len(data)-1
    common=(1<<(degree+1))-1
    records=[]
    prime=65518
    for _ in range(60):
        prime=int(nextprime(prime))
        p=nmod_poly([v%prime for v in data],prime)
        if p.degree()!=degree or p.gcd(p.derivative()).degree()!=0:
            continue
        unit,factors=p.factor()
        reached=1
        product=nmod_poly([int(unit)],prime)
        degrees=[]
        explicit=[]
        check_ring=fmpz_mod_poly_ctx(prime)
        for f,e in factors:
            require(e==1,"选中的模素数不是平方自由特化")
            coeffs=list(map(int,f.coeffs()))
            require(check_ring(coeffs).is_irreducible(),"分解所得因子未通过独立模不可约检验")
            reached|=reached<<f.degree()
            product*=f
            degrees.append(f.degree())
            explicit.append(coeffs)
        require(product==p,"模不可约因子回乘失败")
        common&=reached
        records.append({"prime":prime,"factor_degrees":sorted(degrees),"possible_degrees_count":common.bit_count(),
                        "unit":int(unit),"monic_factors_ascending":explicit})
        if common.bit_count()==2:
            break
    return {"primes":records,"remaining_possible_degrees":[i for i in range(degree+1) if (common>>i)&1],
            "irreducibility_proved":common.bit_count()==2}


def isolation(data):
    from flint import arb,arb_poly,ctx,fmpz_poly
    root=read(HERE.parents[1]/"proofs/global-final/master-20261118.json")["decimal_intervals"]
    bounds=list(map(Q,(root["r_lower"],root["r_upper"])))
    degree=len(data)-1
    norm=sum(map(abs,data))
    sign=lambda y:1 if y>0 else -1 if y<0 else 0
    for precision in (norm.bit_length()+degree+512,2*(norm.bit_length()+degree)+1024):
        with ctx.workprec(precision):
            lo,hi=[arb(q.numerator)/q.denominator for q in bounds]
            poly=arb_poly(fmpz_poly(data))
            signs=sign(poly(lo)),sign(poly(hi))
            if 0 in signs:
                continue
            if signs[0]==signs[1]:
                return {"endpoints":list(map(str,bounds)),"endpoint_signs":list(signs),"unique_root_proved":False}
            mid,width=(lo+hi)/2,(hi-lo)/2
            k=min(degree,(norm.bit_length()+degree)//190+32)
            for terms in sorted({k,min(2*k,degree),min(4*k,degree)}):
                # 仅展开 P' 在中点的前若干项。余项用整数系数的一范数统一控制。
                t=[arb(0) for _ in range(terms)]
                for i in range(degree,0,-1):
                    for j in range(terms-1,0,-1):
                        t[j]=t[j]*mid+t[j-1]
                    t[0]=t[0]*mid+i*data[i]
                error=arb(0)
                for z in reversed(t[1:]):
                    error=(error+abs(z))*width
                if terms<degree:
                    error+=arb(degree*norm)*(1<<degree)*width**terms/(1-width)
                if abs(t[0])>error:
                    return {"endpoints":list(map(str,bounds)),"endpoint_signs":list(signs),
                            "derivative_sign":sign(t[0]),"taylor_terms":terms,"precision_bits":precision,
                            "unique_root_proved":True}
    return {"endpoints":list(map(str,bounds)),"unique_root_proved":False}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        data=read(args.source)["candidate_primitive_coefficients_ascending"]
        require(len(data)>1 and data[-1]>0 and reduce(gcd,data)==1,"整数候选不满足规范化条件")
        write(args.output,{"source_sha256":digest(args.source),"degree":len(data)-1,
                           "irreducibility":degrees_across_primes(data),"root_count":isolation(data),
                           "geometric_identity_proved":False})
    finally:
        close()
