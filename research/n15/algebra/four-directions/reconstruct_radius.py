"""从已提升的首一半径系数恢复公共分母，逐项回模并检查目标实区间。"""

import argparse
from functools import reduce
from math import gcd,lcm
from pathlib import Path

from chart import HERE,read
from context import write
from sample_plane import PRIME,digest
from verify import require


def symmetric(value,modulus):
    value%=modulus
    return value-modulus if value>modulus//2 else value


def separated_denominator(value,modulus):
    # 大部分商很小。出现大商时，当步余数和 Bézout 系数给出短的有理重建。
    u,v=modulus,value%modulus
    b,c=0,1
    best_q,den=0,0
    while v:
        q,rem=divmod(u,v)
        if q>best_q and q.bit_length()>48 and gcd(c,modulus)==1:
            best_q,den=q,abs(c)//gcd(abs(c),v)
        u,v,b,c=v,rem,c,b-q*c
    return den


def candidate_scales(seq,m):
    from flint import fmpz_mat
    yield 1
    d=len(seq)-1
    seen={1}
    den=1
    for i in dict.fromkeys([*range(min(d,12)),*range(max(0,d-12),d),*[j*d//13 for j in range(1,13)]]):
        part=separated_denominator(seq[i],m)
        if part:
            den=lcm(den,part)
            if den not in seen and den<m:
                seen.add(den)
                yield den
    for n in (4,8,12):
        for indices in (list(range(n)),list(range(d-n,d)),
                        [j//2 if j%2==0 else d-j//2-1 for j in range(n)]):
            lattice=[[1,*[symmetric(seq[i],m) for i in indices]]]
            for j in range(n):
                lattice.append([0,*[m*int(j==k) for k in range(n)]])
            basis=fmpz_mat(lattice).lll()
            for row in range(min(3,n+1)):
                scale=abs(int(basis[row,0]))
                if scale and scale not in seen:
                    seen.add(scale)
                    yield scale


def at_target(data):
    from flint import arb,arb_poly,ctx,fmpz_poly
    root=read(HERE.parents[1]/"proofs/global-final/master-20261118.json")["decimal_intervals"]
    ctx.prec=512
    interval=arb(root["r_lower"]).union(arb(root["r_upper"]))
    val=arb_poly(fmpz_poly(data))(interval)
    return val.contains(0)


def reconstruct(path,output):
    src=read(path)
    require(src["prime"]==PRIME and src["original_equations_verified"],"缺少原临界方程的提升验证")
    seq=src["modular_monic_coefficients_ascending"]
    m=PRIME**src["exponent"]
    bits=m.bit_length()
    require(seq[-1]==1 and len(seq)==src["degree"]+1,"首一半径层的格式错误")
    best,profile=None,(0,0)
    for scale in candidate_scales(seq,m):
        trial=[symmetric(scale*v,m) for v in seq]
        heights=sorted(abs(v).bit_length() for v in trial)
        metric=(bits-heights[-1],bits-heights[(len(heights)*3)//4])
        if best is None or metric>profile:
            best,profile=trial,metric
        if metric[0]>=48:
            break
    content=reduce(gcd,best)
    normalized=[c//content for c in best]
    if normalized[-1]<0:
        normalized=[-c for c in normalized]
    valid_mod=all((a-normalized[-1]*b)%m==0 for a,b in zip(normalized,seq))
    survives=profile[0]>=48 and valid_mod and at_target(normalized)
    item={"source_sha256":digest(path),"prime":PRIME,"exponent":src["exponent"],"degree":len(seq)-1,
          "modulus_bits":bits,"minimum_headroom_bits":profile[0],"upper_quartile_headroom_bits":profile[1],
          "passes_common_height_modular_and_real_tests":survives,"characteristic_zero_identity_proved":False}
    if survives:
        item["candidate_primitive_coefficients_ascending"]=normalized
    write(output,item)
    print("precision",bits,"headroom",profile,"candidate",survives,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        reconstruct(args.source,args.output)
    finally:
        close()
