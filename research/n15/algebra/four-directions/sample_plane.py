"""从整数表直接生成半径样本，核对已认证的截面并保留可续算的数据。"""

import argparse
import gzip
import hashlib
import os
from pathlib import Path
import struct
import subprocess
import time
import zipfile

from chart import HERE, read
from context import ROOT, write
from radius_slice import PRIME
from verify_carrier import TABLE
from verify import require

WORK = ROOT / ".tools/four-directions-plane"
CHAIN = ROOT / ".tools/msolve-paged/chain/ucrt64/bin"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def prepare():
    WORK.mkdir(exist_ok=True)
    source = HERE / "plane_samples.cpp"
    binary = WORK / ("samples-"+digest(source)[:16]+".exe")
    env = {**os.environ, "PATH": str(CHAIN)+os.pathsep+os.environ["PATH"]}
    if not binary.exists():
        subprocess.run([str(CHAIN/"g++.exe"), "-O3", "-fopenmp", "-std=c++17", str(source),
                        "-lflint", "-o", str(binary)], env=env, check=True)
    target = WORK / "coefficients.bin"
    verified = read(HERE / "carrier-verified.json")
    require(digest(TABLE) == verified["sha256"], "特征零表摘要不匹配")
    buffer = bytearray()
    with gzip.open(TABLE,"rt",encoding="ascii") as stream:
        next(stream)
        for line in stream:
            j,r,a,d,c = map(int,line.split(","))
            buffer.extend(struct.pack("<5I", j,r-48,a-12,d,c%PRIME))
    if not target.exists():
        with target.open("xb") as stream:
            stream.write(buffer)
    require(digest(target) == hashlib.sha256(buffer).hexdigest(), "模系数缓存不匹配")
    return binary,target,env


def rows(path):
    with path.open("rb") as stream:
        while header := stream.read(8):
            r,n = struct.unpack("<2I",header)
            yield r,struct.unpack(f"<{n}I",stream.read(4*n))


def sample_path(meta):
    path=ROOT/meta["sample_file"]
    expected=meta["sha256"]
    if path.is_file() and digest(path)==expected:
        return path
    with zipfile.ZipFile(HERE/"plane-samples.zip") as archive:
        raw=archive.read(path.name)
    require(hashlib.sha256(raw).hexdigest()==expected,"打包样本摘要不匹配")
    WORK.mkdir(exist_ok=True)
    local=WORK/(expected+".bin")
    if not local.exists():
        with local.open("xb") as stream:
            stream.write(raw)
    require(digest(local)==expected,"本地样本缓存摘要错误")
    return local


def compute(binary,table,env,start,count,threads,label):
    path = WORK / (label+".bin")
    if path.exists():
        raise FileExistsError(path)
    tick = time.perf_counter()
    subprocess.run([str(binary),str(table),str(path),str(start),str(count),str(threads)],env=env,check=True)
    return path, time.perf_counter()-tick


def crosscheck(binary,table,env,label):
    from fractions import Fraction as Q
    from flint import nmod_poly
    for name in ("slice-9-16","slice-2-5"):
        archived = read(HERE/(name+".json.gz"))
        radius = Q(archived["fixed_value"])
        residue = radius.numerator*pow(radius.denominator,-1,PRIME)%PRIME
        path,_ = compute(binary,table,env,residue,1,1,label+"-"+name)
        [(r, coeffs)] = rows(path)
        poly = nmod_poly(archived["coefficients_ascending"],PRIME)
        layer = next(f for f,m in poly.factor_squarefree()[1] if m == 1)
        require(r == residue and list(coeffs)==list(map(int,layer.coeffs())), "网格算法与归档截面不同")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", required=True,type=int)
    parser.add_argument("--count", required=True,type=int)
    parser.add_argument("--threads", default=4,type=int)
    parser.add_argument("--crosscheck", action="store_true")
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", required=True,type=Path)
    args=parser.parse_args()
    if Path(args.label).name!=args.label or args.output.exists():
        raise ValueError("标签需是文件名，输出需为新文件")
    from windows_compute_limits import apply_limits
    snapshot, close = apply_limits()
    try:
        binary,table,env = prepare()
        if args.crosscheck:
            crosscheck(binary,table,env,args.label)
        path, elapsed = compute(binary,table,env,args.start,args.count,args.threads,args.label)
        sizes = {r:len(c) for r,c in rows(path)}
        require(set(sizes)==set(range(args.start,args.start+args.count)), "样本记录未覆盖指定范围")
        write(args.output, {"sample_file":str(path.relative_to(ROOT)),"sha256":digest(path),"count":args.count,
              "first_radius":args.start,"regular_rows":sum(n==1321 for n in sizes.values()),
              "crosschecked":args.crosscheck,"seconds":elapsed,"limits":snapshot()})
    finally:
        close()
