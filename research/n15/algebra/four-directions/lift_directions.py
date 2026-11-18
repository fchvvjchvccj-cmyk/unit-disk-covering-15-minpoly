"""从半径与实方向的模坐标恢复两种镜像，并回代特征零的三个临界方程。"""

import argparse
import gzip
from pathlib import Path

from chart import HERE, read
from context import write
from eliminate import Quadratic, primitive_resultant
from modular_fibers import coefficients
from sample_plane import PRIME, digest
from verify import independent_variables, require
from verify_carrier import TABLE


class Residue:
    """在固定模多项式的商环中运算，所有除法均检查逆元。"""

    def __init__(self, value, modulus):
        self.poly=value%modulus
        self.modulus=modulus

    def cast(self, value):
        return value if isinstance(value,Residue) else Residue(self.poly*0+value,self.modulus)

    def __add__(self, other):
        return Residue(self.poly+self.cast(other).poly,self.modulus)

    __radd__=__add__

    def __neg__(self):
        return Residue(-self.poly,self.modulus)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return (-self)+other

    def __mul__(self, other):
        return Residue(self.poly*self.cast(other).poly,self.modulus)

    __rmul__=__mul__

    def __truediv__(self, other):
        divisor=self.cast(other).poly
        if divisor.gcd(self.modulus).degree()!=0:
            raise ZeroDivisionError("方向环中的除数不是单位")
        return Residue(self.poly*divisor.inverse_mod(self.modulus),self.modulus)

    def __rtruediv__(self, other):
        return self.cast(other)/self


def critical_equations(values):
    e,a,d=independent_variables(values)
    r=(e-1/e)/(a-1/a)
    f=primitive_resultant(e,a,d)
    fe,fa,fd=f.gradient
    re,ra,_=r.gradient
    return [f.value,fd,fe*ra-fa*re]


def read_table():
    from flint import nmod_poly
    rows={}
    with gzip.open(TABLE,"rt",encoding="ascii") as stream:
        next(stream)
        for line in stream:
            part,ir,ia,d,c=map(int,line.split(","))
            key=part,d,ia-12
            if key not in rows:
                rows[key]=[0]*84
            rows[key][ir-48]=c%PRIME
    return {key:nmod_poly(data,PRIME) for key,data in rows.items()}


def recover_one(record,table):
    from flint import nmod_poly, fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
    f=nmod_poly(record["factor_ascending"],PRIME)
    x=nmod_poly(record["X"],PRIME)
    field=fq_default_ctx(modulus=fmpz_mod_poly_ctx(PRIME)(coefficients(f)),fq_type="FQ_NMOD",check_modulus=False)
    ring=fq_default_poly_ctx(field)
    def in_field(p):
        return field(coefficients(p%f))
    rr,xx=ring([field.gen()]),ring([in_field(x)])
    adir=Quadratic(ring.zero(),ring.one(),-ring.one(),xx)
    one=nmod_poly([1],PRIME)
    zero=one*0
    powers=[(one,zero),(zero,one)]
    for _ in range(2,67):
        h,k=powers[-1]
        powers.append((-k,(h+x*k)%f))
    arrays=[[[zero for _ in range(105)] for _ in range(2)] for _ in range(2)]
    for (part,d,a),coefficient in table.items():
        for j in range(2):
            arrays[part][j][d]+=coefficient*powers[a][j]
    p,q=[adir.pair(*(ring([in_field(poly) for poly in cols]) for cols in part)) for part in arrays]
    aa=adir*adir
    n=p*(p+adir.cast(rr)*(aa-1)*q)-aa*q*q
    total=n.norm()
    g=total.gcd(total.derivative())
    out={"factor_ascending":record["factor_ascending"],"degree":f.degree(),"direction_norm_D_degree":total.degree(),
         "direction_gcd_D_degree":g.degree()}
    if g.degree()!=2 or g[0]!=field.one() or g.gcd(g.derivative()).degree()!=0:
        return out
    def entry(poly):
        return Residue(poly,g)
    def pack(v):
        return [list(map(int,v.poly[j].to_list())) for j in range(2)]
    try:
        d=entry(ring.gen())
        a=-entry(n.a)/entry(n.b)
        e=-(entry(p.a)+a*entry(p.b))/(a*(entry(q.a)+a*entry(q.b)))
        radius=(e-1/e)/(a-1/a)
        require((a+1/a-entry(xx)).poly.is_zero() and (radius-entry(rr)).poly.is_zero(),"杆方向或半径回代失败")
        equations=critical_equations(independent_variables([e,a,d]))
        values=[v.value.poly.is_zero() for v in equations]
        jac=[v.gradient for v in equations]
        determinant=sum((jac[0][i]*(jac[1][(i+1)%3]*jac[2][(i+2)%3]-jac[1][(i+2)%3]*jac[2][(i+1)%3])
                         for i in range(3)),d*0)
        out.update(D_quadratic_linear=list(map(int,g[1].to_list())),A=pack(a),E=pack(e),
                   original_equations_zero=values,nonsingular=determinant.poly.gcd(g).degree()==0)
    except ZeroDivisionError:
        out["exceptional_denominator"]=True
    return out


def recover(report,output):
    require(digest(TABLE)==read(HERE/"carrier-verified.json")["sha256"],"整系数表摘要错误")
    results=read(report)
    table=read_table()
    records=[]
    for row in results["irreducible_factors"]:
        if (row.get("direction_gcd_degree")==1 and row["plane_X_degree"]==656
                and not row["X_is_plus_or_minus_two"]):
            record=recover_one(row,table)
            records.append(record)
            print(record["degree"],record["direction_gcd_D_degree"],record.get("original_equations_zero"),flush=True)
    write(output,{"prime":PRIME,"parent_sha256":digest(report),"regular_plane_fibers":records})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"mod-fibers.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        recover(args.source,args.output)
    finally:
        close()
