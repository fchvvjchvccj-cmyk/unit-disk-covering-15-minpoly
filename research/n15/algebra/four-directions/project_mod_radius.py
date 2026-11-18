"""对已插值的模平面多项式完整计算方向结式，归档原系数与平方自由层。"""

import argparse
from array import array
from fractions import Fraction as Q
import gzip
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
import time

from chart import HERE, read
from context import write
from fit_plane import mirrored_layer, one_layer
from sample_plane import WORK, CHAIN, PRIME, digest
from verify import require


def coefficients(path):
    import flint
    vals=array("I")
    with path.open("rb") as stream:
        vals.fromfile(stream,path.stat().st_size//4)
    if sys.byteorder!="little":
        vals.byteswap()
    return flint.nmod_poly(list(vals),PRIME)


def save_layer(path, poly):
    vals=array("I",map(int,poly.coeffs()))
    if sys.byteorder!="little":
        vals.byteswap()
    with path.open("xb") as stream:
        stream.write(gzip.compress(vals.tobytes(),compresslevel=6,mtime=0))


def compute(matrix, report_path, prefix, threads):
    from flint import nmod_poly
    info=read(report_path)
    require(info["variables"]==["r","X"] and digest(matrix)==info["sha256"], "方向曲线文件未绑定")
    n,width=info["degrees"][1],info["degrees"][0]+1
    source=HERE/"radius_discriminant.cpp"
    job=WORK/("radius-"+info["sha256"][:16]+"-"+digest(source)[:16])
    job.mkdir(exist_ok=True)
    env={**os.environ,"PATH":str(CHAIN)+os.pathsep+os.environ["PATH"]}
    exe=job/"project.exe"
    subprocess.run([str(CHAIN/"g++.exe"),"-O3","-fopenmp","-std=c++17",str(source),"-lflint","-o",str(exe)],
                   env=env,check=True)
    raw=job/"matrix.bin"
    with gzip.open(matrix,"rb") as inp,raw.open("wb") as out:
        shutil.copyfileobj(inp,out)
    require(raw.stat().st_size==4*width*(n+1),"方向曲线矩阵维数不符")
    # 小控制例同时经过首项消失的半径，检验齐次结式的零值以及网格次序。
    with tempfile.TemporaryDirectory(dir=job) as tmp:
        temp=Path(tmp); small=temp/"input.bin"; testout=temp/"out.bin"
        small.write_bytes(struct.pack("<12I",1,2,0,1, 0,3,4,0, PRIME-1,1,0,0))
        subprocess.run([str(exe),str(small),"4","2",str(temp/"blocks"),str(testout),"1"],env=env,check=True)
        test=coefficients(testout)
        for r in range(20):
            f=nmod_poly([1+2*r+r**3,3*r+4*r*r,r-1],PRIME)
            require(test(r)==(f.resultant(f.derivative()) if f.degree()==2 else 0),"小截面结式不一致")
    output=job/"result.bin"
    tick=time.perf_counter()
    subprocess.run([str(exe),str(raw),str(width),str(n),str(job/"blocks"),str(output),str(threads)],
                   env=env,check=True)
    result=coefficients(output)
    checks=[]
    for ratio,name in ((Q(9,16),"slice-9-16"),(Q(2,5),"slice-2-5"),(Q(3,7),None)):
        r=ratio.numerator*pow(ratio.denominator,-1,PRIME)%PRIME
        with raw.open("rb") as stream:
            line=[int(nmod_poly(list(struct.unpack(f"<{width}I",stream.read(4*width))),PRIME)(r))
                  for _ in range(n+1)]
        f=nmod_poly(line,PRIME)
        if name:
            require(f==mirrored_layer(one_layer(name))*line[-1] and line[-1]!=0, "额外截面方向层错误")
        direct=int(f.resultant(f.derivative())) if f.degree()==n else 0
        require(direct==result(r), "最后单变量模结式的独立检查失败")
        checks.append([str(ratio),direct])
    repeated=result.gcd(result.derivative())
    radical=result//repeated
    simple=radical//radical.gcd(repeated)
    layers={"result":result,"radical":radical,"multiplicity-one":simple}
    files={}
    for name,poly in layers.items():
        polyfile=prefix.with_name(prefix.name+"-"+name+".bin.gz")
        save_layer(polyfile,poly)
        files[name]={"degree":poly.degree(),"sha256":digest(polyfile),"file":polyfile.name}
    return {"prime":PRIME,"variable":"r","input_matrix_sha256":info["sha256"],"files":files,
            "degree_bound":(2*n-1)*(width-1),"independent_radii":checks,"seconds":time.perf_counter()-tick,
            "exact_for_supplied_modular_matrix":True,"certified_characteristic_zero_relation":False}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix",type=Path,required=True)
    parser.add_argument("--matrix-report",type=Path,required=True)
    parser.add_argument("--prefix",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=4)
    args=parser.parse_args()
    if args.output.exists() or any(args.prefix.parent.glob(args.prefix.name+"-*.bin.gz")):
        raise FileExistsError("末级模结果已存在")
    from windows_compute_limits import apply_limits
    snapshot,close=apply_limits()
    try:
        report=compute(args.matrix,args.matrix_report,args.prefix,args.threads)
        report["limits"]=snapshot()
        write(args.output,report)
    finally:
        close()
