"""严格夹住整数多项式的实根，把其五个有理函数坐标送入原构型的唯一根盒。"""

import argparse
from fractions import Fraction as Q
from pathlib import Path

from chart import read
from context import write
from sample_plane import digest
from scalar_coordinates import NAMES
from verify import require


def arb_q(a):
    from flint import arb
    return arb(a.numerator)/a.denominator


def inner_box(p,coordinates,coarse,root):
    from flint import arb,arb_poly,ctx,fmpz_poly
    lower,upper=map(Q,coarse)
    centers=list(map(Q,root["center"]))
    width=Q(root["width"])
    ip=fmpz_poly(p)
    for digits in (1200,2400,4800):
        with ctx.workprec(4*digits+9000):
            f=arb_poly(ip)
            df=arb_poly(ip.derivative())
            mid=arb_q((lower+upper)/2).mid()
            for _ in range(12):
                deriv=df(mid)
                if deriv.contains(0):
                    break
                mid=(mid-f(mid)/deriv).mid()
            scale=10**digits
            integer=(mid*scale).floor().unique_fmpz()
            if integer is None:
                continue
            a,b=Q(int(integer)-2,scale),Q(int(integer)+2,scale)
            if not lower<a<b<upper:
                continue
            lo,hi=arb_q(a),arb_q(b)
            if not (f(lo)>0 and f(hi)<0):
                continue
            beta=lo.union(hi)
            dvalue=df(beta)
            if dvalue.contains(0):
                continue
            vals=[beta]
            for name in NAMES:
                rec=coordinates[name]
                v=arb_poly(fmpz_poly(rec["numerator_ascending"]))(beta)
                vals.append(v/(rec["integer_denominator"]*dvalue))
            if not all(x>arb_q(c-width) and x<arb_q(c+width) for x,c in zip(vals,centers)):
                continue
            short=[]
            for x,c in zip(vals,centers):
                # 用具体有理端点记录嵌入方向，重新验证小区间的包含。
                k=int((x.mid()*10**55).floor().unique_fmpz())
                aa,bb=Q(k-1,10**55),Q(k+2,10**55)
                require(c-width<aa<bb<c+width and x>arb_q(aa) and x<arb_q(bb),"实嵌入记录过窄")
                short.append([str(aa),str(bb)])
            return {"root_bracket":list(map(str,(a,b))),"endpoint_signs":[1,-1],
                    "coordinate_intervals":short,"precision_bits":ctx.prec}
    raise ValueError("未把该实根的所有坐标严格放进指定源根盒")


def bind(candidate,check,coords_file,target,identity,output):
    rec=read(candidate)
    verified=read(check)
    coords=read(coords_file)
    source=read(target)
    require(verified["source_sha256"]==coords["candidate_sha256"]==digest(candidate)
            and verified["irreducibility"]["irreducibility_proved"] and verified["root_count"]["unique_root_proved"],
            "整系数多项式、坐标、不可约性或根数未绑定")
    require(source["source_image_inside"] and Q(source["contraction"]["contraction_upper"])<1
            and Q(source["contraction"]["inclusion_ratio_upper"])<1, "实构型的局部根盒无效")
    embedded=inner_box(rec["candidate_primitive_coefficients_ascending"],coords["coordinates"],
                       verified["root_count"]["endpoints"],source)
    proved=False
    if identity:
        cert=read(identity)
        require(cert["candidate_sha256"]==digest(candidate) and cert["coordinates_sha256"]==digest(coords_file)
                and cert["irreducibility_sha256"]==digest(check) and cert["all_six_rational_function_identities_proved"]
                and cert["sum_lower_bits"]>cert["upper_bits"],"特征零方程没有证成")
        proved=True
    write(output,{"candidate_sha256":digest(candidate),"irreducibility_sha256":digest(check),
                  "coordinate_sha256":digest(coords_file),"root_box_sha256":digest(target),
                  "identity_sha256":digest(identity) if identity else None,"real_embedding":embedded,
                  "real_embedding_in_unique_source_box":True,"root_equals_global_optimal_radius_proved":proved})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("candidate","check","coordinates","target","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--identity",type=Path)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        bind(args.candidate,args.check,args.coordinates,args.target,args.identity,args.output)
    finally:
        close()
