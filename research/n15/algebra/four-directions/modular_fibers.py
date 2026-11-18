"""逐个不可约模半径因子回代对称平面曲线，精确恢复其重方向。"""

import argparse
import gzip
from pathlib import Path
import struct

from chart import HERE, read
from context import write
from sample_plane import PRIME, digest
from verify import require


def read_poly(name):
    from flint import nmod_poly
    raw=gzip.decompress((HERE/name).read_bytes())
    return nmod_poly([x[0] for x in struct.iter_unpack("<I",raw)],PRIME)


def coefficients(f):
    return list(map(int,f.coeffs()))


def recover():
    from flint import nmod_poly, fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
    info=read(HERE/"plane-symmetric.json")
    matrix=HERE/"plane-symmetric.bin.gz"
    require(info["sha256"]==digest(matrix),"对称方向表摘要不符")
    width=info["degrees"][0]+1
    rows=[nmod_poly(list(a),PRIME) for a in struct.iter_unpack(f"<{width}I",gzip.decompress(matrix.read_bytes()))]
    require(len(rows)==info["degrees"][1]+1,"方向系数个数错误")
    source=read(HERE/"mod-radius.json")["files"]["multiplicity-one"]
    require(digest(HERE/source["file"])==source["sha256"],"半径一重层摘要不符")
    simple=read_poly(source["file"])
    unit,factors=simple.factor()
    product=nmod_poly([int(unit)],PRIME)
    records=[]
    for f,m in sorted(factors,key=lambda item:(item[0].degree(),coefficients(item[0]))):
        modulus=fmpz_mod_poly_ctx(PRIME)(coefficients(f))
        require(m==1 and modulus.is_irreducible(),"输入层的不可约分解不符")
        product*=f
        field=fq_default_ctx(modulus=modulus,fq_type="FQ_NMOD",check_modulus=False)
        ring=fq_default_poly_ctx(field)
        def specialize(p):
            return field(coefficients(p%f))
        curve=ring([specialize(p) for p in rows])
        gcd=curve.gcd(curve.derivative())
        entry={"factor_ascending":coefficients(f),"degree":f.degree(),
               "plane_X_degree":curve.degree(),"direction_gcd_degree":gcd.degree()}
        if gcd.degree()==1:
            x=-gcd[0]/gcd[1]
            slope=ring([specialize(p.derivative()) for p in rows])(x)
            entry.update(X=list(map(int,x.to_list())),X_is_plus_or_minus_two=x in (field(2),field(-2)),
                         ordinary_double=not curve.derivative().derivative()(x).is_zero(),
                         transverse_in_radius=not slope.is_zero())
            require(curve(x).is_zero() and curve.derivative()(x).is_zero(),"方向重根回代失败")
        records.append(entry)
        print(f.degree(),curve.degree(),gcd.degree(),flush=True)
    require(product==simple,"最后一重层的不可约因子回乘失败")
    return {"prime":PRIME,"input_sha256":source["sha256"],"curve_sha256":info["sha256"],
            "irreducible_factors":records,"characteristic_zero_factor_identified":False}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        write(args.output,recover())
    finally:
        close()
