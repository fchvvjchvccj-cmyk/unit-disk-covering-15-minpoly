"""恢复五个实半角坐标的有理单变量表示，并检查六个整数临界式的模回代。"""

import argparse
from pathlib import Path

from chart import HERE,read
from context import write
from modular_fibers import coefficients
from padic_radius import Algebra
from rational_coordinates import reconstruct_coordinate
from sample_plane import WORK,PRIME,digest
from scalar_incidence import coordinates,integer_program
from verify import require

NAMES=["e","z","d","l","b"]


def run(source,candidate,level):
    folder=WORK/("hensel-"+digest(source)[:12])
    states=[read(folder/f"d{s['degree']}-k{level}.json.gz") for s in read(source)["regular_plane_fibers"]]
    data=read(candidate)["candidate_primitive_coefficients_ascending"]
    tape,out=integer_program()
    m=PRIME**level
    from flint import fmpz_mod_poly_ctx
    ring=fmpz_mod_poly_ctx(m)
    total=ring(data)*pow(data[-1],-1,m)
    combined=[ring.zero() for _ in NAMES]
    for state in states:
        field=Algebra(state,level)
        v=coordinates(field,state)
        nodes=[]
        uses={}
        for op,*args in tape:
            if op not in ("var","const"):
                for j in args:
                    uses[j]=uses.get(j,0)+1
        for j in out:
            uses[j]=uses.get(j,0)+1
        for op,*args in tape:
            if op=="var":
                value=v[args[0]]
            elif op=="const":
                value=field.poly([args[0]])
            elif op=="neg":
                value=-nodes[args[0]]
            elif op=="add":
                value=nodes[args[0]]+nodes[args[1]]
            else:
                require(op=="mul","整系数电路含有除法")
                value=field.mul(nodes[args[0]],nodes[args[1]])
            nodes.append(value)
            if op not in ("var","const"):
                for j in args:
                    uses[j]-=1
                    if not uses[j]:
                        nodes[j]=None
        require(all(nodes[i].is_zero() for i in out),"实半角临界方程没有精确消失")
        cofactor,rem=divmod(total,field.f)
        require(rem.is_zero(),"整数半径多项式被此次提升排除")
        derivative=field.f.derivative()
        for i,value in enumerate(v[1:]):
            combined[i]+=cofactor*field.mul(derivative,value)
    found={}
    for name,poly in zip(NAMES,combined):
        seq=coefficients(poly)
        seq.extend([0]*(total.degree()-len(seq)))
        found[name]=reconstruct_coordinate(seq,m,data[-1])
        print(name,None if not found[name] else found[name]["numerator_height_bits"],flush=True)
    return {"candidate_sha256":digest(candidate),"source_sha256":digest(source),"degree":total.degree(),
            "prime":PRIME,"exponent":level,"coordinates":found,"integer_equations_modularly_verified":True}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"mod-anchors.json")
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--level",type=int,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        write(args.output,run(args.source,args.candidate,args.level))
    finally:
        close()
