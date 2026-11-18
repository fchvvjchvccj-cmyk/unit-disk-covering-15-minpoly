"""用非奇异临界方程倍增素数幂精度，同时将半径保持为本原坐标。"""

import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from pathlib import Path
import tempfile
import time

from chart import HERE, read
from context import write
from critical_circuit import evaluate
from modular_fibers import coefficients
from sample_plane import PRIME, WORK, digest
from verify import require


class Algebra:
    """先模首一半径多项式，再加入满足 D²+gD+1=0 的有向根。"""

    def __init__(self,state,level):
        from flint import fmpz_mod_poly_ctx
        self.level=level
        self.mod=PRIME**level
        self.poly=fmpz_mod_poly_ctx(self.mod)
        self.f=self.poly(state["factor_ascending"])
        self.g=self.poly(state["D_quadratic_linear"])
        self.size=self.f.degree()
        self.zero=self.poly.zero()
        self.one=self.poly.one()
        self.r=self.poly.gen()%self.f
        self.reverse_inv=self.f.reverse().inverse_series_trunc(max(1,self.size-1))

    def reduce(self,p):
        if p.degree()<self.size:
            return p
        k=p.degree()-self.size+1
        q=p.reverse().mul_low(self.reverse_inv,k).reverse(k-1)
        return (p-q*self.f).truncate(self.size)

    def mul(self,x,y):
        return self.reduce(x*y)

    def inv(self,x):
        from flint import nmod_poly
        if x.degree()==0:
            return self.poly([pow(int(x[0]),-1,self.mod)])
        small=nmod_poly(coefficients(x),PRIME)
        f=nmod_poly(coefficients(self.f),PRIME)
        v=self.poly(coefficients(small.xgcd(f)[1]%f))
        power=1
        while power<self.level:
            v=self.mul(v,2-self.mul(x,v))
            power*=2
        require(self.mul(x,v)==self.one,"半径商环中的逆元校验失败")
        return v

    def scalar(self,value):
        return Number(self,self.poly([value]),self.zero)

    def unpack(self,values):
        return Number(self,*(self.poly(col) for col in values))


class Number:
    """素数幂上的镜像二次代数元素，系数已对半径多项式规约。"""

    def __init__(self,ring,a,b):
        self.ring,self.a,self.b=ring,a,b

    def __add__(self,v):
        if isinstance(v,int):
            return Number(self.ring,self.a+v,self.b)
        return Number(self.ring,self.a+v.a,self.b+v.b)

    __radd__=__add__

    def __neg__(self):
        return Number(self.ring,-self.a,-self.b)

    def __sub__(self,v):
        return self+-v

    def __rsub__(self,v):
        return -self+v

    def __mul__(self,v):
        if isinstance(v,int):
            return Number(self.ring,self.a*v,self.b*v)
        m=self.ring.mul
        if self.b.is_zero():
            return Number(self.ring,m(self.a,v.a),m(self.a,v.b))
        if v.b.is_zero():
            return v*self
        ac,bd=m(self.a,v.a),m(self.b,v.b)
        cross=m(self.a+self.b,v.a+v.b)-ac-bd
        return Number(self.ring,ac-bd,cross-m(self.ring.g,bd))

    __rmul__=__mul__

    def inverse(self):
        m=self.ring.mul
        other=self.a-m(self.ring.g,self.b)
        inv=self.ring.inv(m(self.a,other)+m(self.b,self.b))
        return Number(self.ring,m(other,inv),-m(self.b,inv))

    def __truediv__(self,v):
        return self*(self.ring.scalar(v) if isinstance(v,int) else v).inverse()

    def __rtruediv__(self,v):
        return self.inverse()*v


def values(base,state):
    point=[base.unpack(state["E"]),base.unpack(state["A"]),Number(base,base.zero,base.one)]
    if "L" in state:
        point.extend(base.unpack(state[k]) for k in ("L","B"))
    return point


def exact_zero(x):
    return x.a.is_zero() and x.b.is_zero()


def solve(matrix,b):
    n=len(b)
    rows=[[*row,v] for row,v in zip(matrix,b)]
    for j in range(n):
        for i in range(j,n):
            if exact_zero(rows[i][j]):
                continue
            try:
                inverse=rows[i][j].inverse()
                break
            except (ZeroDivisionError,ValueError):
                continue
        else:
            raise ValueError("临界方程的雅可比缺少可逆主元")
        rows[i],rows[j]=rows[j],rows[i]
        rows[j][j+1:]=[x*inverse for x in rows[j][j+1:]]
        for i in range(n):
            if i==j:
                continue
            value=rows[i][j]
            if not exact_zero(value):
                for k in range(j+1,n+1):
                    rows[i][k]-=value*rows[j][k]
    return [row[-1] for row in rows]


def double(state,level):
    small,big=Algebra(state,level),Algebra(state,2*level)
    at_small=values(small,state)
    n=len(at_small)
    vec=evaluate(at_small,True)
    zero=at_small[0]*0
    vec=[zero+q if isinstance(q,int) else q for q in vec]
    require(all(exact_zero(v) for v in vec[:n]),"当前精度未满足原临界方程")
    point=values(big,state)
    residual=evaluate(point)
    def divide(c):
        vals=coefficients(c)
        require(all(v%small.mod==0 for v in vals),"残差不被旧精度整除")
        return small.poly([v//small.mod for v in vals])
    high=[Number(small,-divide(v.a),-divide(v.b)) for v in residual]
    delta=solve([vec[(i+1)*n:(i+2)*n] for i in range(n)],high)
    def lift(c):
        return big.poly(coefficients(c))*small.mod
    updated=[x+Number(big,lift(y.a),lift(y.b)) for x,y in zip(point,delta)]
    e,a,d=updated[:3]
    radius=(e-1/e)/(a-1/a)
    y=d+1/d
    require(radius.b.is_zero() and y.b.is_zero(),"镜像间半径或方向实部不一致")
    shift=radius.a-big.r
    require(all(int(c)%small.mod==0 for c in shift.coeffs()),"新半径与旧本原坐标不相合")
    def transport(c):
        return c-big.mul(c.derivative(),shift)
    inv=big.inv(d.b)
    def change_basis(v):
        b=big.mul(v.b,inv)
        return [coefficients(transport(v.a-big.mul(b,d.a))),coefficients(transport(b))]
    out={"degree":state["degree"],"factor_ascending":coefficients(transport(big.f)),
         "D_quadratic_linear":coefficients(transport(-y.a)),"A":change_basis(a),"E":change_basis(e)}
    if n==5:
        out.update(L=change_basis(updated[3]),B=change_basis(updated[4]))
    require(len(out["factor_ascending"])==big.size+1 and out["factor_ascending"][-1]==1,
            "坐标回换没有保留首一半径多项式")
    return out


def advance(job):
    state,current,work=job
    cache=work/f"d{state['degree']}-k{2*current}.json.gz"
    if cache.is_file():
        data=read(cache)
    else:
        data=double(state,current)
        temp=cache.with_name(cache.stem+".partial.gz")
        write(temp,data)
        temp.replace(cache)
    print("lift",2*current,state["degree"],flush=True)
    return data


def check_endpoint(state_and_level):
    state,current=state_and_level
    alg=Algebra(state,current)
    point=values(alg,state)
    residual=evaluate(point)
    require(all(exact_zero(v) for v in residual),"终点回代原方程失败")
    radius=(point[0]-1/point[0])/(point[1]-1/point[1])
    require(radius.b.is_zero() and radius.a==alg.r,"终点半径坐标回代失败")
    return state["degree"]


def run(source,level,output,workers=1,work_root=WORK):
    initial=read(source)
    rows=initial["regular_plane_fibers"]
    require(all(r.get("original_equations_zero")==[True]*3 and r["nonsingular"] for r in rows),
            "初始纤维未全部通过原方程")
    work=work_root/("hensel-"+digest(source)[:12])
    work.mkdir(parents=True,exist_ok=True)
    degree=sum(r["degree"] for r in rows)
    stages=[]
    tick=time.perf_counter()
    current=1
    with (ProcessPoolExecutor(max_workers=workers) if workers>1 else nullcontext()) as pool:
        dispatch=pool.map if pool else map
        while current<level:
            ordered=sorted(rows,key=lambda row:-row["degree"])
            rows=list(dispatch(advance,((s,current,work) for s in ordered)))
            current*=2
            stages.append(current)
        checked=list(dispatch(check_endpoint,((s,current) for s in rows)))
        require(sum(checked)==degree,"终点核验未覆盖全部模半径")
    from flint import fmpz_mod_poly_ctx
    ring=fmpz_mod_poly_ctx(PRIME**current)
    total=ring.one()
    for state in rows:
        total*=ring(state["factor_ascending"])
    require(total.degree()==degree,"提升后丢失半径")
    write(output,{"prime":PRIME,"exponent":current,"degree":degree,"source_sha256":digest(source),
                  "modular_monic_coefficients_ascending":coefficients(total),"stages":stages,
                  "original_equations_verified":True,"seconds":time.perf_counter()-tick,
                  "characteristic_zero_polynomial_verified":False})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE/"mod-directions.json")
    p.add_argument("--level",type=int,required=True)
    p.add_argument("--workers",type=int,choices=range(1,5),default=1)
    p.add_argument("--fresh",action="store_true",help="在独立目录从最初的模纤维重算")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    if args.level<1 or args.level&(args.level-1) or args.output.exists():
        raise ValueError("精度指数须为二的幂，输出须为新文件")
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        with (tempfile.TemporaryDirectory() if args.fresh else nullcontext(str(WORK))) as scratch:
            run(args.source,args.level,args.output,args.workers,Path(scratch))
    finally:
        close()
