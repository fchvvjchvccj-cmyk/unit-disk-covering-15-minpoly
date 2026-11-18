"""从方向纤维恢复左接点和中间杆，使素数幂提升直接使用低次接触方程。"""

import argparse
from pathlib import Path

from chart import HERE, read
from context import write
from critical_circuit import evaluate
from eliminate import Quadratic, anchor_quadratic, endpoint, primitive_branch, rod_quartic
from padic_radius import Algebra, exact_zero, values
from modular_fibers import coefficients
from sample_plane import PRIME,digest
from verify import require


def linear_gcd(polynomial):
    a=list(polynomial)
    b=[i*q for i,q in enumerate(a) if i]
    while b:
        inverse=b[-1].inverse()
        while len(a)>=len(b):
            n=len(a)-len(b)
            c=a[-1]*inverse
            for j in range(len(b)):
                a[n+j]-=c*b[j]
            while a and exact_zero(a[-1]):
                a.pop()
        a,b=b,a
    require(len(a)==2,"中间杆的重因子并非一次")
    return -a[0]/a[1]


def run(source,output):
    rows=[]
    for row in read(source)["regular_plane_fibers"]:
        require(row["nonsingular"] and all(row["original_equations_zero"]),"初始模几何不满足临界方程")
        base=Algebra(row,1)
        e,a,d=values(base,row)
        r,t,v,w,wbar=endpoint(e,a,d)
        h0,h1,h2=anchor_quadratic(v,w,wbar,r)
        l=Quadratic(r*0,r*0+1,-h0/h2,-h1/h2)
        branch=primitive_branch(l,l.cast(e),l.cast(r),l.cast(t))
        left=-branch.a/branch.b
        b=linear_gcd(rod_quartic(left,e,r,t))
        require(exact_zero(h0+left*(h1+left*h2)),"左接点未满足二次方程")
        require(all(exact_zero(f) for f in evaluate([e,a,d,left,b])),"五元接触临界方程回代失败")
        pack=lambda x:[coefficients(x.a),coefficients(x.b)]
        rows.append({**row,"L":pack(left),"B":pack(b)})
        print(row["degree"],flush=True)
    write(output,{"prime":PRIME,"source_sha256":digest(source),"regular_plane_fibers":rows})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"mod-directions.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        run(args.source,args.output)
    finally:
        close()
