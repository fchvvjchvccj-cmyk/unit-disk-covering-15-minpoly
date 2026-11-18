"""在 W=A E 的二次商代数中先规约，再用两个根的幂和计算齐次结式。"""

from collections import defaultdict
import gzip
from itertools import zip_longest
from pathlib import Path

from eliminate import Quadratic


def split(polynomial, base):
    rows = defaultdict(dict)
    for powers, coefficient in polynomial.terms():
        r, a, e, d, left = map(int, powers)
        rows[left, e][r, a, d] = coefficient
    return {key: base.from_dict(terms) for key, terms in rows.items()}


def carrier_components(h, branch, note):
    from flint import fmpz_mpoly_ctx
    base = fmpz_mpoly_ctx.get(("r", "A", "D"))
    r, a, d = base.gens()
    terms = split(branch, base)
    e_degree = max(i for _, i in terms)
    w = Quadratic(base.constant(0), base.constant(1), a*a, r*(a*a-1))
    powers = [w.cast(base.constant(1))]
    for _ in range(e_degree):
        powers.append(powers[-1]*w)
    for i in range(e_degree+1):
        powers[i] = powers[i]*a**(e_degree-i)
    coefficients = []
    for j in range(9):
        p = w.cast(base.constant(0))
        for i, power in enumerate(powers):
            if (j,i) in terms:
                p += power*terms[j,i]
        coefficients.append(p)
    hs = split(h, base)
    h0, h1, h2 = [hs[j,0] for j in range(3)]
    traces = [base.constant(2), -h1]
    for _ in range(2,9):
        traces.append(-h1*traces[-1]-h0*h2*traces[-2])
    h0powers = [h0**i for i in range(9)]
    res = w.cast(base.constant(0))
    for j, p in enumerate(coefficients):
        # 齐次结式按较高指数 j 合并；不用生成两套通用大型符号矩阵。
        row = p*h0powers[j]
        for i in range(j):
            row += coefficients[i]*(h0powers[i]*traces[j-i])
        res = res*h2+p*row
        note("homogeneous_sum_"+str(j), res.a)
    divisor = d**4
    p, q = res.a//divisor, res.b//divisor
    if p*divisor != res.a or q*divisor != res.b or max(p.degrees()[2], q.degrees()[2]) != 104:
        raise ValueError("结构化第一结式未通过精确除法或次数核验")
    return w.pair(p,q), 16+2*e_degree


def save_components(path, pair, exponent):
    """按行流式保存整数，避免为大型多项式再建立一份 Python 系数表。"""
    import hashlib
    import json
    path = Path(path)
    temp = path.with_suffix(path.suffix+".partial")
    if path.exists() or temp.exists():
        raise FileExistsError(path)
    with temp.open("xb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=1, mtime=0) as stream:
            stream.write((json.dumps({"A_exponent": exponent, "r_exponent": 24, "basis": ["1","A*E"],
                                      "variables": ["r","A","D"]})+"\n").encode())
            for index, poly in enumerate((pair.a,pair.b)):
                for (r,a,d), coeff in poly.terms():
                    stream.write(f"{index},{r},{a},{d},{coeff}\n".encode("ascii"))
    # 精确回读与内存中的规范项序逐项比较；仅在检查结束后发布文件。
    with gzip.open(temp, "rt", encoding="ascii") as stream:
        header = json.loads(next(stream))
        if header["A_exponent"] != exponent:
            raise ValueError("中间输入的幂次错误")
        expected = ((j,*map(int,m),int(c)) for j,poly in enumerate((pair.a,pair.b)) for m,c in poly.terms())
        actual = (tuple(map(int,line.split(","))) for line in stream)
        if any(left != right for left,right in zip_longest(expected, actual)):
            raise ValueError("中间输入回读失败")
    temp.rename(path)
    with path.open("rb") as stream:
        return {"path": str(path), "bytes": path.stat().st_size, "sha256": hashlib.file_digest(stream,"sha256").hexdigest()}
