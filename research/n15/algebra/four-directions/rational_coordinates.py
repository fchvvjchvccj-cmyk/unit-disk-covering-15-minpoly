"""由各个提升分支的迹插值方向坐标，使用半径多项式的导数作公共分母。"""

import argparse
from functools import reduce
from math import gcd
from pathlib import Path

from chart import HERE,read
from context import write
from padic_radius import coefficients
from reconstruct_radius import candidate_scales,symmetric
from sample_plane import WORK,PRIME,digest
from verify import require


def reconstruct_coordinate(residues,m,lead):
    seq=[(v*lead)%m for v in residues]
    bits=m.bit_length()
    for scale in candidate_scales(seq,m):
        array=[symmetric(scale*v,m) for v in seq]
        height=max(abs(v).bit_length() for v in array)
        if bits-height<48:
            continue
        common=reduce(gcd,array,scale)
        denominator=scale//common
        nums=[v//common for v in array]
        require(all((a-denominator*v)%m==0 for a,v in zip(nums,seq)),
                "坐标规范化后的同余式不一致")
        return {"integer_denominator":denominator,"numerator_ascending":nums,
                "numerator_height_bits":max(abs(v).bit_length() for v in nums)}
    return None


def trace_coordinates(source,level,candidate):
    from flint import fmpz_mod_poly_ctx
    m=PRIME**level
    ring=fmpz_mod_poly_ctx(m)
    folder=WORK/("hensel-"+digest(source)[:12])
    states=[read(folder/f"d{s['degree']}-k{level}.json.gz") for s in read(source)["regular_plane_fibers"]]
    factors=[ring(s["factor_ascending"]) for s in states]
    total=ring.one()
    for f in factors:
        total*=f
    data=read(candidate)["candidate_primitive_coefficients_ascending"]
    require(ring(data)==total*data[-1],"特征零半径候选与本次提升不符")
    entries={"g":ring.zero(),**{k+str(j):ring.zero() for k in ("A","E","L","B") for j in range(2)}}
    for state,f in zip(states,factors):
        cof,rest=divmod(total,f)
        require(rest.is_zero(),"半径因子未整除首一乘积")
        deriv=f.derivative()
        for name in entries:
            values=state["D_quadratic_linear"] if name=="g" else state[name[0]][int(name[1])]
            entries[name]+=cof*(ring(values)*deriv%f)
    found={}
    degree=total.degree()
    for name,poly in entries.items():
        raw=coefficients(poly)
        raw.extend([0]*(degree-len(raw)))
        found[name]=reconstruct_coordinate(raw,m,data[-1])
        print(name,found[name]["numerator_height_bits"] if found[name] else None,flush=True)
    return {"candidate_sha256":digest(candidate),"source_sha256":digest(source),"degree":degree,"prime":PRIME,
            "exponent":level,"coordinates":found,"exact_characteristic_zero_coordinates_proved":False}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"mod-anchors.json")
    p.add_argument("--level",type=int,required=True)
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        write(args.output,trace_coordinates(args.source,args.level,args.candidate))
    finally:
        close()
