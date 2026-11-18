"""以有理互补交点公式重建目标接触链，不在中间两步开平方。"""


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]


def add(a, b):
    return [a[i]+b[i] for i in range(2)]


def subtract(a, b):
    return [a[i]-b[i] for i in range(2)]


def scale(a, k):
    return [v*k for v in a]


def boundary(a, b, r):
    ka, kb = 1+dot(a, a)-r*r, 1+dot(b, b)-r*r
    nx, ny = ka*b[1]-kb*a[1], a[0]*kb-b[0]*ka
    area = a[0]*b[1]-a[1]*b[0]
    return nx*nx+ny*ny-4*area*area


def circumcenter(a, b, c):
    d, e = subtract(b, a), subtract(c, a)
    determinant = 2*(d[0]*e[1]-d[1]*e[0])
    db, dc = dot(b, b)-dot(a, a), dot(c, c)-dot(a, a)
    return [(db*e[1]-dc*d[1])/determinant, (d[0]*dc-e[0]*db)/determinant]


def common_boundary(a, b, r):
    ka, kb = 1+dot(a, a)-r*r, 1+dot(b, b)-r*r
    determinant = 2*(a[0]*b[1]-a[1]*b[0])
    return [(ka*b[1]-kb*a[1])/determinant, (a[0]*kb-b[0]*ka)/determinant]


def previous_boundary(a, u, r):
    t = 2*a[0]-u
    x = (t*t+1-r*r)/(2*t)
    return [x, (1+dot(a, a)-r*r-2*a[0]*x)/(2*a[1])]


def reconstruct(c0, u, b, r):
    a = [c0[0], -c0[1]]
    q0 = previous_boundary(a, u, r)
    p0 = circumcenter(a, [u, 0*u], b)
    q1 = subtract(scale(a, (1+dot(a, a)-r*r)/dot(a, a)), q0)
    p1 = subtract(add(a, b), p0)
    c3 = subtract(add(q1, p1), a)
    c4 = subtract(add(b, scale(c3, (1-r*r)/dot(c3, c3))), a)
    return {"c3": c3, "c4": c4, "q0": q0, "p0": p0, "q1": q1, "p1": p1}
