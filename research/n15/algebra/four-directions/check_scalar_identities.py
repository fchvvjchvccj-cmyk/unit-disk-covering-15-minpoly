"""在独立有理函数域复核实半角闭环的清分母恒等式及弦长公式。"""

import argparse
from pathlib import Path

from chart import HERE
from context import write
from scalar_incidence import constraints,integer_constraints
from verify import require


def check():
    from sympy import QQ
    from sympy.polys.fields import field
    ring,r,e,z,d,l,b=field("r,e,z,d,l,b",QQ)
    v=[r,e,z,d,l,b]
    u,h,q=constraints(v)
    ui,hi,qi=integer_constraints(v)
    den=(1-z)*(1-z*d*d)
    rq=1-z*d*d+r*(1+z*d*d)
    require(u==ui and h*(den**4)*(rq**2)==hi and q*(1-z*b*b)**2*den**4==qi,
            "通用实半角清分母公式不相等")
    _,z,l,v,wx,wy=field("z,l,v,x,y",QQ)
    wx2=wx*wx-z*wy*wy+1
    squared=(wx-(1+z*l*l)/(1-z*l*l))**2-z*(wy-2*l/(1-z*l*l))**2
    expected=wx2-2*wx+4*z*wy*l-z*(wx2+2*wx)*l*l
    require((1-z*l*l)*squared==expected,"边界弦长的通用恒等式错误")
    require((1+z*l*v)**2-z*(l+v)**2==(1-z*l*l)*(1-z*v*v),"弦法向量恒等式错误")
    line=lambda x,y:(1+z*l*v)*x-z*(l+v)*y+z*l*v-1
    for point in (l,v):
        require(line((1+z*point*point)/(1-z*point*point),2*point/(1-z*point*point))==0,
                "弦直线没有同时通过两端")
    return {"field":"QQ(r,e,z,d,l,b)","cleared_equations_terms":[len(f.numer) for f in (ui,hi,qi)],
            "H_multiplier":"((1-z)(1-z*d^2))^4*(1-z*d^2+r*(1+z*d^2))^2",
            "Q_multiplier":"(1-z*b^2)^2*((1-z)(1-z*d^2))^4",
            "both_chord_endpoints":True,"squared_chord_length":True,"squared_normal_length":True}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        write(args.output,check())
    finally:
        close()
