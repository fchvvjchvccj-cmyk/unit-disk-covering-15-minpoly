"""对有限域截面插值，分别记录建模样本、留出样本与额外非连续检验点。"""

import argparse
from fractions import Fraction as Q
import gzip
from pathlib import Path
import struct
import subprocess
import tempfile

from chart import HERE, read
from context import write
from sample_plane import WORK, CHAIN, PRIME, digest, rows, sample_path
from verify import require


def one_layer(name):
    from flint import nmod_poly
    data = read(HERE / (name+".json.gz"))
    poly = nmod_poly(data["coefficients_ascending"], PRIME)
    return next(f for f,m in poly.factor_squarefree()[1] if m == 1)


def mirrored_layer(f):
    from flint import nmod_poly
    reciprocal = f.gcd(f.reverse())
    if reciprocal.degree()%2 or reciprocal != reciprocal.reverse():
        raise ValueError("对称因子不满足偶次互反关系")
    half = reciprocal.degree()//2
    # (A^k+A^-k) 对 X=A+A^-1 满足二阶递推；用 Clenshaw 法合并。
    x=nmod_poly([0,1],PRIME)
    b1=b2=nmod_poly([],PRIME)
    for j in range(half,0,-1):
        b0=x*b1-b2+int(reciprocal[half+j])
        b2,b1=b1,b0
    return x*b1-2*b2+int(reciprocal[half])


def fit(archives,output,report,symmetric,training_override=None):
    from flint import nmod_poly
    import os
    if output.exists() or report.exists():
        raise FileExistsError("插值输出已存在")
    paths, seen = [], set()
    for file in archives:
        meta = read(file)
        path=sample_path(meta)
        for r,row in rows(path):
            if r in seen:
                raise ValueError("重复使用半径样本")
            seen.add(r)
        paths.append(path)
    calibration = one_layer("slice-a3")
    if symmetric:
        calibration=calibration.gcd(one_layer("slice-a-third"))
    fixed = list(map(int,calibration.coeffs()))
    width, anchor = (657, 10*pow(3,-1,PRIME)%PRIME) if symmetric else (1321,3)
    training=training_override or len(fixed)
    env = {**os.environ,"PATH":str(CHAIN)+os.pathsep+os.environ["PATH"]}
    with tempfile.TemporaryDirectory(dir=WORK,prefix="fit-") as temp:
        temp = Path(temp)
        exe,raw,at3 = temp/"fit.exe",temp/"matrix.bin",temp/"fixed.bin"
        at3.write_bytes(struct.pack(f"<{len(fixed)}I",*fixed))
        if symmetric:
            converted=temp/"symmetric.bin"
            with converted.open("xb") as stream:
                for path in paths:
                    for r,coeffs in rows(path):
                        if len(coeffs)!=1321:
                            continue
                        layer=mirrored_layer(nmod_poly(list(coeffs),PRIME))
                        vals=list(map(int,layer.coeffs()))
                        if len(vals)==width:
                            stream.write(struct.pack(f"<{width+2}I",r,width,*vals))
            paths=[converted]
        subprocess.run([str(CHAIN/"g++.exe"),"-O3","-std=c++17",str(HERE/"plane_interpolate.cpp"),
                        "-lflint","-o",str(exe)],env=env,check=True)
        run=subprocess.run([str(exe),str(at3),str(raw),str(training),str(width),str(anchor),*map(str,paths)],env=env,
                           capture_output=True,text=True)
        require(run.returncode==0,"平面插值与留出样本不符："+run.stderr)
        used,degree=map(int,run.stdout.split())
        payload=raw.read_bytes()
        require(len(payload)==width*training*4,"系数矩阵长度错误")
        controls=[(Q(9,16),one_layer("slice-9-16")),(Q(2,5),one_layer("slice-2-5"))]
        results=[[] for _ in controls]
        for row in struct.iter_unpack(f"<{training}I",payload):
            p=nmod_poly(list(row),PRIME)
            for dest,(r,_) in zip(results,controls):
                dest.append(int(p(r.numerator*pow(r.denominator,-1,PRIME)%PRIME)))
        for coeffs,(_,expected) in zip(results,controls):
            if symmetric:
                expected=mirrored_layer(expected)
            require(nmod_poly(coeffs,PRIME)==expected*coeffs[-1] and coeffs[-1]!=0,
                    "额外有理半径截面不匹配")
        if degree+1 < training:
            payload=b"".join(struct.pack(f"<{degree+1}I",*row[:degree+1])
                             for row in struct.iter_unpack(f"<{training}I",payload))
        with output.open("xb") as stream:
            stream.write(gzip.compress(payload,compresslevel=6,mtime=0))
    write(report,{"prime":PRIME,"variables":["r","X" if symmetric else "A"],"degrees":[degree,width-1],
                  "storage":f"按方向变量升幂逐行，每行 {degree+1} 个小端 uint32，按 r 升幂",
                  "normalization":"X=10/3 取两个互反方向截面的首一公因子" if symmetric else "A=3 取已计算重数一层的首一半径截面",
                  "training_radii":training,"withheld_radii":used-training,"extra_rational_radii":["9/16","2/5"],
                  "sha256":digest(output),"universal_identity_proved":False,"characteristic_zero":False})


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("archives",type=Path,nargs="+")
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--report",required=True,type=Path)
    p.add_argument("--symmetric",action="store_true")
    p.add_argument("--training",type=int)
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        fit(args.archives,args.output,args.report,args.symmetric,args.training)
    finally:
        close()
