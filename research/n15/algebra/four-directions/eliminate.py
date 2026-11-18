"""在单位圆的复坐标中消去中间两杆和左接点，系数运算仅使用有理数。"""


def endpoint(e, a, d):
    r = (e-1/e)/(a-1/a)
    shift = d+1/d-1
    transfer = (e-1)/r-a-shift
    right = (1+r*d)/(1+r/d)
    inner = 1+r*(shift+a)
    inner_bar = 1+r*(shift+1/a)
    return r, transfer, right, inner, inner_bar


def anchor_quadratic(right, inner, inner_bar, r):
    k = inner-right
    h = right*inner_bar-1
    product = k*h
    return [-inner*product-r*r*k*k,
            (inner*inner_bar+1)*product-2*r*r*product,
            -inner_bar*product-r*r*h*h]


def rod_quartic(left, e, r, transfer):
    unit = 1-r*r
    a = -r*r*transfer*e
    b = left*r*r*transfer
    c = -left*r*unit*e
    d = r*unit
    k = left*(unit-r*transfer*e)-unit*e+r*transfer
    p0, p1, p2 = b+d, k+b*transfer, b+c
    q0, q1, q2 = a+d, k+a*transfer, a+c
    return [p0*q0, p0*q1+p1*q0, p0*q2+p1*q1+p2*q0+(a-b)*(a-b),
            p1*q2+p2*q1, p2*q2]


def discriminant27(coefficients):
    e, d, c, b, a = coefficients
    i = 12*a*e-3*b*d+c*c
    j = 72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c
    return 4*i*i*i-j*j


def pencil_discriminant(x, y, z, w, t):
    """两条二次曲线的铅笔判别式，剥去解线性方程引入的四次幂。"""
    c0 = z*z
    c1 = t*w*z+w*w-4*x*y-2*(x+y)*z+z*z
    c2 = t*t*x*y-t*w*(x+y)+w*w+(x-y)*(x-y)-4*x*y-2*(x+y)*z
    c3 = (x-y)*(x-y)
    return c1*c1*c2*c2-4*c3*c1*c1*c1-4*c0*c2*c2*c2-27*c0*c0*c3*c3+18*c0*c1*c2*c3


def primitive_branch(left, e, r, transfer):
    unit = 1-r*r
    a, b = -r*r*transfer*e, left*r*r*transfer
    c, d = -left*r*unit*e, r*unit
    k = left*(unit-r*transfer*e)-unit*e+r*transfer
    return pencil_discriminant(b+c, b+d, a-b, k+b*transfer, transfer)


class Quadratic:
    """在 X²=linear*X+constant 的二维代数中保留常数项与一次项。"""

    def __init__(self, a, b, constant, linear):
        self.a, self.b, self.constant, self.linear = a, b, constant, linear

    def pair(self, a, b):
        return Quadratic(a, b, self.constant, self.linear)

    def cast(self, other):
        if isinstance(other, Quadratic):
            return other
        return self.pair(other, self.a*0)

    def __add__(self, other):
        other = self.cast(other)
        return self.pair(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return self.pair(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.cast(other))

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        other = self.cast(other)
        product = self.b*other.b
        return self.pair(self.a*other.a+product*self.constant,
                         self.a*other.b+self.b*other.a+product*self.linear)

    __rmul__ = __mul__

    def norm(self):
        return self.a*(self.a+self.linear*self.b)-self.constant*self.b*self.b


def eliminant(e, a, d):
    r, transfer, right, inner, inner_bar = endpoint(e, a, d)
    h0, h1, h2 = anchor_quadratic(right, inner, inner_bar, r)
    left = Quadratic(0*r, 1+0*r, -h0/h2, -h1/h2)
    quartic = rod_quartic(left, left.cast(e), left.cast(r), left.cast(transfer))
    return discriminant27(quartic).norm()


def primitive_resultant(e, a, d):
    r, transfer, right, inner, inner_bar = endpoint(e, a, d)
    h0, h1, h2 = anchor_quadratic(right, inner, inner_bar, r)
    left = Quadratic(0*r, 1+0*r, -h0/h2, -h1/h2)
    branch = primitive_branch(left, left.cast(e), left.cast(r), left.cast(transfer))
    power = h2*h2
    power = power*power
    return (power*power)*branch.norm()
