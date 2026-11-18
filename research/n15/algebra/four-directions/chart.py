"""用四条有向单位杆在有理数域上参数化八变量接触载体。"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
TRANSFER = HERE.parent / "contact-transfer-20261024"
sys.path.insert(0, str(TRANSFER))
from context import RESULTS, execute, prune, read


def direction(s):
    square = s*s
    return (1-square)/(1+square), 2*s/(1+square)


def forward(parameters):
    (ax, ay), (bx, by), (cx, cy), (dx, dy) = map(direction, parameters)
    k = ax + 2*cx + 2*(dx-bx) - 1
    r = -2*k/(k*k+ay*ay)
    u = 1+r*(2*(dx-bx)-1)
    px, py = u+r*ax, r*ay
    v, w = px+r*bx, py+r*by
    z, h, y = 1+r*k, py+r*cy, r*dy
    t = 2*r*((v-u)*ay-w*ax)
    return [t, z, h, y, u, v, w, r]


def inverse(point):
    t, z, h, y, u, v, w, r = point
    dx = v-u
    length = dx*dx+w*w
    px, py = u+(dx-w*t/length)/2, (w+dx*t/length)/2
    x6 = (1+r+u+2*(v-px))/2
    return [py/(r+px-u), (w-py)/(r+v-px),
            (h-py)/(r+(z-px)/2), y/(r+x6-1)]


def old_critical():
    source = read(RESULTS / "primitive16.json.gz")
    return prune(source, source["outputs"]["critical_equations"])


def critical(parameters, circuit):
    data, outputs = circuit
    values = execute(data, forward(parameters))
    return [values[i] for i in outputs]


def circle_closure(parameters):
    """通过两个确定的边界接点与右端内接点的外接半径闭合接触链。"""
    (ax, ay), (bx, by), (cx, cy), (dx, dy) = map(direction, parameters)
    _, z, h, y, u, v, w, r = forward(parameters)
    px, py = u+r*ax, r*ay
    a = (px+r*cx, h)
    one = 1+0*r
    c6 = (1+r*dx, y)
    q0 = (z, py)

    def dot(a, b):
        return a[0]*b[0]+a[1]*b[1]

    def other_anchor(c, q):
        factor = (1-r*r)/dot(c, c)+1
        return (factor*c[0]-q[0], factor*c[1]-q[1])

    q1 = other_anchor(a, q0)
    c3 = (q1[0]+v-px, q1[1]+w-py)
    q2 = other_anchor(c3, q1)
    factor = (1-r*r)/dot(c3, c3)
    c4 = (v-a[0]+factor*c3[0], w-h+factor*c3[1])
    left, right = other_anchor(c4, q2), other_anchor(c6, (one, 0*r))
    inner = (px+2*r*bx, py)
    d = (left[0]-inner[0], left[1]-inner[1])
    e = (right[0]-inner[0], right[1]-inner[1])
    gap = (left[0]-right[0], left[1]-right[1])
    area = d[0]*e[1]-d[1]*e[0]
    return dot(d, d)*dot(e, e)*dot(gap, gap)-4*r*r*area*area, (left, right, inner, area, q2)


def geometric_critical(parameters):
    from verify import independent_variables
    variables = independent_variables(parameters)
    closure, _ = circle_closure(variables)
    radius = forward(variables)[-1]
    # 第二方向作为非零主元，保留固定半径的三个独立切向条件。
    pivot = 1
    return [closure.value] + [
        radius.gradient[pivot]*closure.gradient[j]-closure.gradient[pivot]*radius.gradient[j]
        for j in (0, 2, 3)]
