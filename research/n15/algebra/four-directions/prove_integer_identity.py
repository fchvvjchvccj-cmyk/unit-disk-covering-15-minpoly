"""用整数特化的商群阶和 Sylvester 行列式界证明六个整系数余式为零。"""

import argparse
from dataclasses import dataclass
from math import prod
import os
from pathlib import Path
import subprocess
import tempfile
import time

from chart import HERE,read
from context import write
from sample_plane import CHAIN,WORK,digest
from scalar_coordinates import NAMES
from scalar_incidence import integer_constraints,integer_program
from verify import independent_variables,require


@dataclass
class Size:
    degrees:list
    coordinate_degree:int
    norm:int

    def derivative(self,j):
        ds=self.degrees.copy()
        ds[j]=max(0,ds[j]-1)
        return Size(ds,max(0,self.coordinate_degree-int(j>0)),self.norm*self.degrees[j])

    def __add__(self,other):
        return Size(list(map(max,zip(self.degrees,other.degrees))),
                    max(self.coordinate_degree,other.coordinate_degree),self.norm+other.norm)

    __sub__=__add__

    def __mul__(self,other):
        return Size([a+b for a,b in zip(self.degrees,other.degrees)],
                    self.coordinate_degree+other.coordinate_degree,self.norm*other.norm)


def equation_bounds():
    from flint import fmpz_mpoly_ctx
    ring=fmpz_mpoly_ctx.get(("r",*NAMES))
    data=[]
    for poly in integer_constraints(ring.gens()):
        terms=list(poly.terms())
        data.append(Size(list(map(int,poly.degrees())),max(sum(map(int,m[1:])) for m,c in terms),
                         sum(abs(int(c)) for _,c in terms)))
    u,h,q=data
    tangent=[h.derivative(i)*q.derivative(4)-q.derivative(i)*h.derivative(4) for i in (1,2,3)]
    return [u,h,q,q.derivative(5),tangent[2],tangent[1]*u.derivative(1)-tangent[0]*u.derivative(2)]


def log_upper(n):
    return max(0,n-1).bit_length()


def resultant_bounds(p,coords):
    degree=len(p)-1
    lp=log_upper(sum(map(abs,p)))
    norm_deriv=sum(i*abs(v) for i,v in enumerate(p))
    extra=[]
    for name in NAMES:
        item=coords[name]
        den=item["integer_denominator"]
        n=sum(map(abs,item["numerator_ascending"]))
        # 各输入共享 P'，只有坐标的整数分母需要按各自的多次数清去。
        divisor=den*norm_deriv
        extra.append(log_upper(den)+log_upper(max(1,(n+divisor-1)//divisor)))
    reports=[]
    for size in equation_bounds():
        degree_s=size.degrees[0]+size.coordinate_degree*(degree-1)
        log_s=(log_upper(size.norm)+size.coordinate_degree*log_upper(norm_deriv)
               +sum(i*j for i,j in zip(size.degrees[1:],extra)))
        reports.append({"degrees":size.degrees,"coordinate_total_degree":size.coordinate_degree,
                        "source_l1_bound":size.norm,"substituted_degree_bound":degree_s,
                        "substituted_log2_l1_bound":log_s,"resultant_log2_bound":degree_s*lp+degree*log_s})
    return reports


def prepare_input(folder,p,coords):
    path=folder/"proof-input.txt"
    tape,outputs=integer_program()
    codes={"var":0,"const":1,"neg":2,"add":3,"mul":4}
    with tempfile.NamedTemporaryFile(mode="w",dir=folder,encoding="ascii",newline="\n",delete=False) as fp:
        temp=Path(fp.name)
        def poly(data):
            fp.write(str(len(data))+"\n")
            fp.writelines(str(c)+"\n" for c in data)
        poly(p)
        for name in NAMES:
            fp.write(str(coords[name]["integer_denominator"])+"\n")
            poly(coords[name]["numerator_ascending"])
        fp.write(str(len(tape))+"\n")
        for op,*args in tape:
            fp.write(" ".join(map(str,[codes[op],args[0],args[1] if len(args)>1 else 0]))+"\n")
        fp.write(" ".join(map(str,outputs))+"\n")
    temp.replace(path)
    return path


class Local:
    """用于独立正反控制，在一个整数特化模数上直接作前向微分。"""
    def __init__(self,n,m):
        self.n,self.m=n%m,m

    def __add__(self,b):
        return Local(self.n+(b.n if isinstance(b,Local) else b),self.m)

    __radd__=__add__

    def __neg__(self):
        return Local(-self.n,self.m)

    def __sub__(self,b):
        return self+-b

    def __rsub__(self,b):
        return -self+b

    def __mul__(self,b):
        return Local(self.n*(b.n if isinstance(b,Local) else b),self.m)

    __rmul__=__mul__


def independent_checks(p,coords,first,cutoff,early):
    from flint import fmpz_poly,fmpz
    from sympy import primerange
    poly=fmpz_poly(p)
    dp=poly.derivative()
    forbidden=fmpz(prod(primerange(2,cutoff+1)))*abs(int(p[-1]))
    for name in NAMES:
        forbidden*=coords[name]["integer_denominator"]
    evidence=[]
    for i,bits in enumerate(early):
        k=first+i
        m=abs(poly(k))
        der=dp(k)
        g=m.gcd(forbidden*der)
        while g>1:
            m//=g
            g=m.gcd(g*g)
        require(int(m).bit_length()-1==bits,"整数模数的独立位数核验失败")
        vals=[Local(fmpz(k),m)]
        for name in NAMES:
            item=coords[name]
            value=fmpz_poly(item["numerator_ascending"])(k)
            vals.append(Local(value*pow(der*item["integer_denominator"],-1,m),m))
        def residual(x):
            u,h,q=integer_constraints(independent_variables(x))
            v=[h.gradient[j]*q.gradient[4]-q.gradient[j]*h.gradient[4] for j in (1,2,3)]
            return [u.value,h.value,q.value,q.gradient[5],v[2],v[1]*u.gradient[1]-v[0]*u.gradient[2]]
        require(all(v.n==0 for v in residual(vals)),"独立微分得到非零余数")
        wrong=vals.copy()
        wrong[1]=wrong[1]+1
        require(any(v.n!=0 for v in residual(wrong)),"破坏坐标的反控制未被拒绝")
        evidence.append({"k":k,"bits":bits,"altered_coordinate_rejected":True})
    return evidence


def prove(candidate,certificate,coordinate_file,points,threads,output):
    p=read(candidate)["candidate_primitive_coefficients_ascending"]
    checked=read(certificate)
    require(checked["source_sha256"]==digest(candidate) and checked["irreducibility"]["irreducibility_proved"],
            "缺少与当前 P 绑定的有理不可约性证书")
    packet=read(coordinate_file)
    require(packet["candidate_sha256"]==digest(candidate) and all(packet["coordinates"].get(k) for k in NAMES),
            "五个实坐标尚未全部重建")
    coords=packet["coordinates"]
    reports=resultant_bounds(p,coords)
    bound=max(q["resultant_log2_bound"] for q in reports)
    source=HERE/"integer_identity.cpp"
    folder=WORK/("identity-"+digest(coordinate_file)[:12]+"-"+digest(source)[:12])
    folder.mkdir(parents=True,exist_ok=True)
    exe=folder/"identity.exe"
    env={**os.environ,"PATH":str(CHAIN)+os.pathsep+os.environ["PATH"]}
    if not exe.is_file():
        subprocess.run([str(CHAIN/"g++.exe"),"-std=c++17","-O3","-fopenmp",str(source),"-lgmp","-o",str(exe)],
                       env=env,check=True)
    path=prepare_input(folder,p,coords)
    first,cap=65536,100000
    expected=min(points,cap)
    tick=time.perf_counter()
    bits=[]
    with subprocess.Popen([str(exe),str(path),str(first),str(cap),str(expected),str(threads),str(bound)],
                          stdout=subprocess.PIPE,env=env,text=True) as process:
        total=0
        for line in process.stdout:
            v=int(line)
            require(v>=0,"有效商群阶数界为负")
            bits.append(v)
            total+=v
            if len(bits)%1024==0:
                print("points",len(bits),"bits",total,"needed",bound,flush=True)
        require(process.wait()==0,"整数点的原式同余检验失败")
    require(len(bits)<=cap and (len(bits)==expected or total>bound),"整数特化数量不完整")
    controls=independent_checks(p,coords,first,cap,bits[:3])
    write(output,{"degree":len(p)-1,"candidate_sha256":digest(candidate),"irreducibility_sha256":digest(certificate),
                  "coordinates_sha256":digest(coordinate_file),"native_source_sha256":digest(source),
                  "first_integer":first,"stripped_prime_cutoff":cap,"bounds":reports,
                  "integer_point_bits":bits,"independent_controls":controls,"sum_lower_bits":total,"upper_bits":bound,
                  "all_six_rational_function_identities_proved":total>bound,
                  "seconds":time.perf_counter()-tick})
    print(len(bits),total,bound,total>bound,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate",required=True,type=Path)
    p.add_argument("--irreducibility",required=True,type=Path)
    p.add_argument("--coordinates",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--points",type=int,default=100000)
    p.add_argument("--threads",type=int,choices=range(1,7),default=4)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        prove(args.candidate,args.irreducibility,args.coordinates,args.points,args.threads,args.output)
    finally:
        close()
