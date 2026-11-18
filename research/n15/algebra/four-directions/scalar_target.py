"""将原最优构型映入实半角临界系统，并核验局部唯一性所用的全部主元。"""

import argparse
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path

from chart import HERE,circle_closure,direction,forward,read
from context import write
from scalar_incidence import integer_constraints,integer_program
from verify import away,box_values,encode,independent_variables,require
from verify_chart import certify_system


@lru_cache(maxsize=1)
def tape():
    return integer_program()


def equations(params):
    ops,outs=tape()
    values=[]
    for op,*args in ops:
        if op=="var":
            v=params[args[0]]
        elif op=="const":
            v=args[0]
        elif op=="neg":
            v=-values[args[0]]
        elif op=="add":
            v=values[args[0]]+values[args[1]]
        else:
            v=values[args[0]]*values[args[1]]
        values.append(v)
    return [values[i] for i in outs]


def image(params):
    a,b,c,d=params
    _,ex,_,_,_,_,_,r=forward(params)
    _,ay=direction(a)
    left=circle_closure(params)[1][0]
    return [r,r*ay/((1+ex)*a),-a*a,d/a,left[1]/((1+left[0])*a),b/a]


def certify():
    old=read(HERE/"target-fiber.json")
    params=list(map(Q,old["refined_center"]))
    source=box_values(params,Q(old["width"]),1600)
    actual=image(source)
    center=[Q(round(x*10**65),10**65) for x in image(params)]
    width=Q(1,10**28)
    require(all(mid-width<Q(encode(x)[0])<=Q(encode(x)[1])<mid+width
                for mid,x in zip(center,actual)),"原构型没有落入实半角根盒")
    source_pivot=forward(independent_variables(source))[-1].gradient[2]
    source_guard=away(source_pivot,"来源方向 C 的半径偏导")
    box=box_values(center,width,1024)
    r,e,z,d,l,b=box
    u,h,q=integer_constraints(independent_variables(box))
    kappa=1-r*r
    dx=(1+z*d*d)/(1-z*d*d)
    ax=(1+z)/(1-z)
    t=2*z*e/(1-z)-ax-(2*dx-1)
    guards={name:away(value,name) for name,value in {
        "U_e":u.gradient[1],"H_l":h.gradient[4],"Q_l":q.gradient[4],
        "K":r*r*t*(z*l*e-1),"r":r,"kappa":kappa,"z":z,
        "1_plus_r_Dx":1+r*dx,"den_A":1-z,"den_D":1-z*d*d,
        "den_B":1-z*b*b,"den_L":1-z*l*l,"den_E":1-z*e*e}.items()}
    return {"center":list(map(str,center)),"width":str(width),"source_image_inside":True,
            "source_radius_partial_C":source_guard,"pivots":guards,
            "contraction":certify_system(equations,center,width,bits=1024,inverse_digits=130),
            "variables":["r","e","z","d","l","b"]}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        write(args.output,certify())
    finally:
        close()
