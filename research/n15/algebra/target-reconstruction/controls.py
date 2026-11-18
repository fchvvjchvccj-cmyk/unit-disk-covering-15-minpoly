"""用已知不可约多项式的实根检验同一搜索与接受端；这些数不是N=15结果。"""

from fractions import Fraction as F
import json
import os
import subprocess
import sys
import time

from intervals import Dyadic
from refine import HERE, ROOT, sha, write_new
from verify_lattices import gs_lower_bounds, polynomial_interval
from windows_compute_limits import apply_limits


def remainder(coefficients, modulus):
    output = list(map(F, coefficients))
    while output and output[-1] == 0:
        output.pop()
    while len(output) >= len(modulus):
        multiplier = output[-1]/modulus[-1]
        offset = len(output)-len(modulus)
        for i, coefficient in enumerate(modulus):
            output[offset+i] -= multiplier*coefficient
        while output and output[-1] == 0:
            output.pop()
    return output


def main():
    import mpmath as mp
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    started = time.perf_counter()
    try:
        output = HERE / "controls"
        output.mkdir(exist_ok=True)
        chain = ROOT / ".tools/msolve-paged/chain/ucrt64"
        env = os.environ.copy()
        env["PATH"] = str(chain / "bin")+os.pathsep+env["PATH"]
        binary = ROOT / ".tools/target-relation-lll.exe"
        reports = []
        for name, leading, constant, true_degree, search_degree, pslq_dps in (
                ("sqrt2", 1, 2, 2, 16, 360),
                ("root17_of_2", 1, 2, 17, 32, 1440),
                ("height20_degree7", 10**20, 3, 7, 16, 720)):
            begin = time.perf_counter()
            known = [-constant]+[0]*(true_degree-1)+[leading]
            digits = (search_degree+1)*38+50
            arithmetic = Dyadic(19000)
            with mp.workdps(5760):
                value = mp.root(mp.mpf(constant)/leading, true_degree)
                center = F(int(mp.nint(value*10**5500)), 10**5500)
                radius = F(1, 10**5400)
                left, right = arithmetic.rational(center-radius), arithmetic.rational(center+radius)
                low = polynomial_interval(known, left, arithmetic)
                high = polynomial_interval(known, right, arithmetic)
                if not low[1] < 0 < high[0] or not center-radius > 0:
                    raise ValueError("已知根的严格正实隔离失败")
                entries = [int(mp.nint(10**digits*value**i)) for i in range(search_degree+1)]
            input_path = output / (name+".input.txt")
            input_path.open("x", encoding="ascii").write(str(search_degree+1)+"\n"+"\n".join(map(str, entries))+"\n")
            path = output / (name+".lattice.json")
            with input_path.open("rb") as incoming, path.open("xb") as outgoing:
                subprocess.run([str(binary)], stdin=incoming, stdout=outgoing, stderr=subprocess.PIPE,
                               env=env, check=True, timeout=60)
            proposal = json.loads(path.read_text(encoding="utf-8"))
            basis = [list(map(int, row)) for row in proposal["basis"]]
            inverse = [list(map(int, row)) for row in proposal["inverse_transform"]]
            for i, row in enumerate(basis):
                if row[-1] != sum(x*y for x, y in zip(row[:-1], entries)):
                    raise ValueError("对照格输入回乘错误")
                for j, column in enumerate(zip(*inverse)):
                    if sum(x*y for x, y in zip(row[:-1], column)) != (i == j):
                        raise ValueError("对照格不幺模")
            scalar = (arithmetic.rational(center-radius)[0], arithmetic.rational(center+radius)[1])
            power = arithmetic.one
            for entry in entries:
                error = arithmetic.sub(arithmetic.integer(power, 10**digits), arithmetic.rational(entry))
                if arithmetic.absolute(error) > arithmetic.scale:
                    raise ValueError("对照幂输入误差超预算")
                power = arithmetic.mul(power, scalar)
            relations = [row[:-1] for row in basis if max(map(abs, row[:-1])) <= 10**24
                         and not remainder(row[:-1], known)]
            if not relations:
                raise ValueError("LLL未找回已知关系")
            # 低维对照实际检验排除门，不能误报不存在已知关系。
            diagonal, scale = gs_lower_bounds(basis, 1024)
            threshold = ((search_degree+1)+(search_degree+1)**2)*10**48
            excluded = min(v[0] for v in diagonal) > threshold*scale
            if excluded:
                raise ValueError("排除门误报已知关系不存在")
            with mp.workdps(pslq_dps):
                scalar_mp = mp.root(mp.mpf(constant)/leading, true_degree)
                pslq = mp.pslq([scalar_mp**i for i in range(true_degree+1)],
                              tol=mp.mpf(10)**(-pslq_dps+60), maxcoeff=10**24+1, maxsteps=4000)
            if pslq is None or remainder(pslq, known):
                raise ValueError("PSLQ未找回已知关系")
            report = {"name": name, "is_n15_target": False, "verified": True,
                      "known_minimal_polynomial_low_to_high": list(map(str, known)),
                      "irreducibility_reason": "Eisenstein at 2" if constant == 2 else "Eisenstein at 3",
                      "true_degree": true_degree, "search_degree": search_degree,
                      "coefficient_budget_inclusive": str(10**24), "scale_digits": digits,
                      "root_interval": [str(center-radius), str(center+radius)],
                      "strict_root_isolation_verified": True,
                      "lll_exact_relations_count": len(relations),
                      "recovered_lll_relation": list(map(str, min(relations, key=lambda row: sum(v*v for v in row)))),
                      "all_recovered_relations_mod_known_polynomial_zero": True,
                      "false_exclusion": excluded, "pslq_candidate": pslq,
                      "pslq_dps": pslq_dps, "pslq_maxsteps": 4000,
                      "input_sha256": sha(input_path), "lattice_sha256": sha(path),
                      "seconds": time.perf_counter()-begin}
            reports.append(report)
            print(json.dumps({k: report[k] for k in ("name", "verified", "seconds", "pslq_candidate")}), flush=True)
        write_new(HERE / "controls-verified.json",
                  {"verified": True, "reports": reports, "is_n15_target": False,
                   "binary_sha256": sha(binary), "seconds": time.perf_counter()-started,
                   "windows_job_limits": limits()})
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
