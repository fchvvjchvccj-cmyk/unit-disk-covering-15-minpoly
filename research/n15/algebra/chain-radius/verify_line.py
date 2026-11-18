"""重验Hermite数据和多项式模基的幺模证书，不信任产生端的次数结论。"""

import argparse
import time

from common import HERE, bounded, read, write, sha


def verify(certificate, samples):
    from flint import nmod_poly
    prime = certificate["prime"]
    P = lambda a: nmod_poly(a, prime)
    M, Y = P(certificate["modulus"]), P(certificate["hermite_interpolant"])
    B = [[P(a) for a in row] for row in certificate["basis"]]
    U = [[P(a) for a in row] for row in certificate["unimodular_transform"]]
    X = P([0, 1])
    expected = P([1])
    selected = samples[:certificate["sample_count"]]
    if len(selected) != certificate["sample_count"] or len({s["r"] for s in selected}) != len(selected):
        raise ValueError("采样数量不足或采样点重复")
    for sample in selected:
        t = sample["r"]
        if int(Y(t)) != sample["value"] or int(Y.derivative()(t)) != sample["derivative_r"]:
            raise ValueError("Hermite插值没有匹配记录的精确取值与导数")
        expected *= (X-t)**2
    if M != expected or Y.degree() >= M.degree():
        raise ValueError("模数不是全部采样点的二重乘积，或插值次数错误")
    if len(B) != 2 or len(U) != 2 or any(len(row) != 2 for row in B+U):
        raise ValueError("模基证书维数不符")
    determinant = U[0][0]*U[1][1]-U[0][1]*U[1][0]
    if determinant.degree() != 0 or determinant.is_zero():
        raise ValueError("多项式变换没有非零常数行列式")
    for i in range(2):
        if B[i] != [U[i][0]*M+U[i][1]*Y, U[i][1]]:
            raise ValueError("返回模基与原始同余模不等价")
    degrees = [max(p.degree() for p in row) for row in B]
    positions = [max(i for i, p in enumerate(row) if p.degree() == d) for row, d in zip(B, degrees)]
    if positions[0] == positions[1]:
        raise ValueError("两个模基行没有不同的最高项位置")
    excluded = min(degrees) > certificate["degree_budget_inclusive"]
    if excluded != certificate["no_bounded_rational_function_matches_all_recorded_jets"]:
        raise ValueError("产生端的范围结论与接受端不一致")
    return {"verified": True, "sample_count": len(selected), "basis_row_degrees": degrees,
            "leading_positions": positions, "degree_budget": certificate["degree_budget_inclusive"],
            "no_bounded_rational_function_matches_recorded_jets": excluded,
            "argument": "N-DY被M整除的全体向量由幺模变换后的两行生成；不同最高位置给出可预测次数，任何非零组合的最大分量次数不低于两行次数的最小值。",
            "scope": "只针对记录的有限域直线数据，不是最优半径次数界，也不是特征零消元恒等式"}


def run(args):
    source = HERE / "line-runs" / args.tag
    with bounded() as snapshot:
        tick = time.perf_counter()
        samples = read(source / "samples.json.gz")
        reports = []
        for path in sorted(source.glob("module-d*.json.gz")):
            report = verify(read(path), samples)
            report.update(certificate_sha256=sha(path), certificate=path.name)
            reports.append(report)
            print(report, flush=True)
        if not reports:
            raise ValueError("没有可重放的模基证书")
        write(source / args.output, {"verified": True, "reports": reports,
                                          "samples_sha256": sha(source / "samples.json.gz"),
                                          "seconds": time.perf_counter()-tick, "limits": snapshot()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--output", default="verified.json")
    run(parser.parse_args())
