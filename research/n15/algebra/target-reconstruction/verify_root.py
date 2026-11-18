"""用纯整数向外舍入区间重验新盒的自映、压缩及父盒包含。"""

import json
from fractions import Fraction as F
import sys
import time

from intervals import Dyadic, derivatives, evaluate_rows
from refine import HERE, ROOT, sha, write_new
from windows_compute_limits import apply_limits


def verify(root):
    started = time.perf_counter()
    parent_path = HERE.parent / "minpoly-final/axis25-root.json"
    system_path = HERE.parent / "minpoly-final/axis25-system.json"
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    system = json.loads(system_path.read_text(encoding="utf-8"))
    if root["system_sha256"] != sha(system_path) or root["parent_root_sha256"] != sha(parent_path):
        raise ValueError("父系统或根证书摘要改变")
    if root["variables"] != system["variables"] or len(root["xnum"]) != 25:
        raise ValueError("变量布局错误")
    center = [F(int(v), int(root["Qx"])) for v in root["xnum"]]
    radius = F(root["rho"])
    if radius <= 0 or root["interval_bits"] < 19000:
        raise ValueError("新盒半宽或区间精度错误")
    parent_center = [F(int(v), int(parent["Qx"])) for v in parent["xnum"]]
    distance = max(abs(v-w)+radius for v, w in zip(center, parent_center))
    if distance >= F(parent["rho"]):
        raise ValueError("新盒不严格包含于父唯一根盒")
    arithmetic = Dyadic(root["interval_bits"])
    exact = [arithmetic.rational(v) for v in center]
    box = [(arithmetic.rational(v-radius)[0], arithmetic.rational(v+radius)[1])
           for v in center]
    inverse = [[arithmetic.rational(F(int(v), int(root["Qy"]))) for v in row]
               for row in root["Ynum"]]
    if len(inverse) != 25 or any(len(row) != 25 for row in inverse):
        raise ValueError("预条件矩阵维数错误")
    terms = [[(t["coefficient"], tuple(map(tuple, t["powers"]))) for t in row]
             for row in system["equations"]]
    values = evaluate_rows(arithmetic, terms, exact)
    flat = [column for row in derivatives(terms, 25) for column in row]
    jacobian_flat = evaluate_rows(arithmetic, flat, box)
    jacobian = [jacobian_flat[i*25:(i+1)*25] for i in range(25)]
    contraction = []
    residual = []
    for i, row in enumerate(inverse):
        error = 0
        for j in range(25):
            product = arithmetic.zero
            for k in range(25):
                product = arithmetic.add(product, arithmetic.mul(row[k], jacobian[k][j]))
            error += arithmetic.absolute(arithmetic.sub(arithmetic.one if i == j else arithmetic.zero, product))
        contraction.append(error)
        value = arithmetic.zero
        for coefficient, f in zip(row, values):
            value = arithmetic.add(value, arithmetic.mul(coefficient, f))
        residual.append(arithmetic.absolute(value))
    q = F(max(contraction), arithmetic.scale)
    inclusion = max(F(e, arithmetic.scale)/radius+F(c, arithmetic.scale)
                    for e, c in zip(residual, contraction))
    if not q < 1 or not inclusion < 1:
        raise ValueError("未通过严格自映及压缩")
    return {"verified": True, "equations_checked": 25,
            "parent_box_strict_containment": True,
            "parent_distance_upper": str(distance),
            "rho": str(radius), "interval_bits": arithmetic.bits,
            "contraction_upper": str(q), "inclusion_ratio_upper": str(inclusion),
            "contraction_less_than_1e_minus_140": q < F(1, 10**140),
            "inclusion_ratio_less_than_1e_minus_90": inclusion < F(1, 10**90),
            "system_sha256": sha(system_path), "parent_root_sha256": sha(parent_path),
            "refined_root_sha256": sha(HERE / "refined-root.json"),
            "seconds": time.perf_counter()-started,
            "argument": "T(z)=z-Yf(z)严格自映并压缩新盒；新盒包含于已有唯一根盒，故为同一目标根。"}


def main():
    limits, close = apply_limits(memory_bytes=8_000_000_000)
    try:
        root = json.loads((HERE / "refined-root.json").read_text(encoding="utf-8"))
        report = verify(root)
        report["windows_job_limits"] = limits()
        write_new(HERE / "refined-root-verified.json", report)
        print(json.dumps({k: report[k] for k in ("verified", "seconds", "contraction_less_than_1e_minus_140", "inclusion_ratio_less_than_1e_minus_90")}), flush=True)
    finally:
        close()


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    main()
