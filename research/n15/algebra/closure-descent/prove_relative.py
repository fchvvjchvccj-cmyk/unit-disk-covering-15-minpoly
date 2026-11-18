"""以通用整数多项式恒等式核验闭环线性化、相对范数和临界乘子。"""

import argparse
import time

from descent_support import HERE, SOURCE, read, write, sha, bounded


def verify():
    import sympy as S
    ax, ay, d, e0, e1, ac, ec, beta = S.symbols("ax ay d e0 e1 ac ec beta")
    e = e0+e1*d
    A, B, C = 1+e1*e1, 2*e0*e1+ac+ec*e1, e0*e0+ec*e0-beta
    g = A*d*d+B*d+C
    ka, kb = ax*ax+ay*ay-beta, d*d+e*e-beta
    area = ax*e-ay*d
    H = (ka*e-kb*ay)**2+(ax*kb-ka*d)**2-4*area*area
    linear_kb = -ac*d-ec*e
    nx, ny = ka*e-linear_kb*ay, ax*linear_kb-ka*d
    reduced = S.Poly(S.expand(nx*nx+ny*ny-4*area*area), d)
    U0, U1, U2 = [reduced.nth(j) for j in range(3)]
    J0, J1 = U0-U2*C/A, U1-U2*B/A
    multiplier = U2+A*(-2*ay*nx+2*ax*ny+g*(ax*ax+ay*ay))
    residue = S.Poly(S.cancel(A*(H-J0-J1*d)-g*multiplier), ax, ay, d, e0, e1, ac, ec, beta)
    if not residue.is_zero:
        raise ValueError("通用闭环线性化恒等式失败")
    j0, j1, qa, qb, qc = S.symbols("j0 j1 qa qb qc")
    resultant = qa*j0*j0-qb*j0*j1+qc*j1*j1
    identity = resultant-qa*(j0+j1*d)*(j0+j1*(-qb/qa-d))-j1*j1*(qa*d*d+qb*d+qc)
    if S.cancel(identity) != 0:
        raise ValueError("二次关系的相对范数恒等式失败")
    a, b, t, V, da, db, dV = S.symbols("a b t V da db dV")
    norm = a*a-V*b*b
    if S.expand(norm-(a+b*t)*(a-b*t)-b*b*(t*t-V)) != 0:
        raise ValueError("系数对范数恒等式失败")
    norm_differential = 2*a*da-dV*b*b-2*V*b*db
    selected_differential = da+db*t+b*dV/(2*t)
    if S.cancel(norm_differential.subs({a: -b*t, V: t*t})+2*b*t*selected_differential) != 0:
        raise ValueError("相对范数的临界微分乘子失败")
    dA, dB, dC, dj0, dj1 = S.symbols("dA dB dC dj0 dj1")
    derivative = (dA*j0*j0+2*qa*j0*dj0-dB*j0*j1-qb*dj0*j1-qb*j0*dj1
                  +dC*j1*j1+2*qc*j1*dj1)
    dd = -(dA*d*d+dB*d+dC)/(2*qa*d+qb)
    source_derivative = dj0+dj1*d+j1*dd
    critical_identity = derivative+j1*(2*qa*d+qb)*source_derivative
    if S.cancel(critical_identity.subs({j0: -j1*d, qc: -qa*d*d-qb*d})) != 0:
        raise ValueError("首级闭环临界微分乘子失败")
    return {"verified": True, "characteristic_zero_universal_identities": 5,
            "closure_linearization": "A*(H-J0-J1*d)=g*[U2+A*(-2*ay*nx+2*ax*ny+g*(ax^2+ay^2))]",
            "quadratic_resultant": "F=A*J0^2-B*J0*J1+C*J1^2=A*(J0+J1*d)*(J0+J1*(-B/A-d))+J1^2*g",
            "symmetric_norm": "N(a+b*t)=(a+b*t)*(a-b*t)+b^2*(t^2-V)=a^2-V*b^2",
            "critical_multiplier": "在a+b*t=0,t^2=V上，dN=(-2*b*t)*d(a+b*t)",
            "first_critical_multiplier": "在H=0,g=0上，dF=-J1*(2*A*d+B)*dH",
            "local_chart_argument": "每级一次项系数及二次根简单导数为单位时，可由零式有理恢复被消去的根，范数另一个共轭因子也是单位；临界方程与闭环式同乘非零因子，逐级保留目标局部临界环。",
            "scope": "通用特征零恒等式和单位条件下的局部等价；未证明全局零维、目标半径多项式或不可约性。"}


def run(args):
    with bounded() as snapshot:
        tick = time.perf_counter()
        report = verify()
        report.update(source_circuit_sha256=sha(SOURCE), verifier_sha256=sha(__file__),
                      seconds=time.perf_counter()-tick, limits=snapshot())
        write(HERE / args.output, report)
        print({"verified": True, "universal_identities": 5, "seconds": report["seconds"]}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="relative-identities.json.gz")
    run(parser.parse_args())
