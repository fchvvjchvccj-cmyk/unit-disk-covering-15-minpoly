"""素域有理函数上的迭代二次代数；每次运算立即按单首项关系约分。"""

import time


class TowerStop(Exception):
    pass


class Rational:
    prime = None
    zero_poly = None
    one_poly = None
    deadline = None
    degree_limit = None
    metrics = None

    def __init__(self, numerator, denominator=None, cancel=True):
        from flint import nmod_poly
        n = nmod_poly([numerator], self.prime) if isinstance(numerator, int) else numerator
        d = self.one_poly if denominator is None else denominator
        if d.is_zero():
            raise ZeroDivisionError("有理函数分母为零")
        if n.is_zero():
            d = self.one_poly
        elif cancel and not d.is_one():
            common = n.gcd(d)
            if not common.is_one():
                n, d = n//common, d//common
        leading = int(d.leading_coefficient())
        if leading != 1:
            inverse = pow(leading, -1, self.prime)
            n, d = n*inverse, d*inverse
        self.n, self.d = n, d
        degree = int(max(n.degree(), d.degree()))
        self.metrics["max_coefficient_degree"] = max(self.metrics["max_coefficient_degree"], degree)
        if degree > self.degree_limit:
            raise TowerStop("单个有理系数次数超过上限")
        if time.perf_counter() > self.deadline:
            raise TowerStop("二次塔运算达到墙钟上限")

    def cast(self, other):
        return other if isinstance(other, Rational) else Rational(other)

    def is_zero(self):
        return self.n.is_zero()

    def is_one(self):
        return self.n == self.d

    def __neg__(self):
        return Rational(-self.n, self.d, cancel=False)

    def __add__(self, other):
        other = self.cast(other)
        self.metrics["rational_additions"] += 1
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        if self.d == other.d:
            return Rational(self.n+other.n, self.d)
        common = self.d.gcd(other.d)
        a, b = self.d//common, other.d//common
        return Rational(self.n*b+other.n*a, a*other.d)

    __radd__ = __add__

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        self.metrics["rational_multiplications"] += 1
        if self.is_zero() or other.is_zero():
            return Rational(0)
        if self.is_one():
            return other
        if other.is_one():
            return self
        g, h = self.n.gcd(other.d), other.n.gcd(self.d)
        return Rational((self.n//g)*(other.n//h), (self.d//h)*(other.d//g), cancel=False)

    __rmul__ = __mul__

    def inverse(self):
        if self.is_zero():
            raise ZeroDivisionError("有理函数为零，不能求逆")
        return Rational(self.d, self.n, cancel=False)

    def __truediv__(self, other):
        return self*self.cast(other).inverse()

    def __rtruediv__(self, other):
        return self.cast(other)*self.inverse()

    def encode(self):
        return {"n": list(map(int, self.n.coeffs())), "d": list(map(int, self.d.coeffs()))}


class Algebra:
    def __init__(self):
        self.relations = []
        self.zero, self.one = Rational(0), Rational(1)
        self.products = self.inversions = 0

    def scalar(self, value):
        return Element(self, 0, (value if isinstance(value, Rational) else Rational(value),))

    def add_level(self, a, b):
        level = len(self.relations)
        a, b = self.scalar(0).cast(a).promote(level), self.scalar(0).cast(b).promote(level)
        self.relations.append((a.coefficients, b.coefficients))
        coefficients = [self.zero]*(2**(level+1))
        coefficients[2**level] = self.one
        return Element(self, level+1, tuple(coefficients))

    def add(self, a, b):
        return tuple(x+y for x, y in zip(a, b))

    def subtract(self, a, b):
        return tuple(x-y for x, y in zip(a, b))

    def multiply(self, a, b, level):
        self.products += 1
        if all(x.is_zero() for x in a) or all(x.is_zero() for x in b):
            return (self.zero,)*len(a)
        if level == 0:
            return (a[0]*b[0],)
        if all(x.is_zero() for x in a[1:]):
            return tuple(a[0]*x for x in b)
        if all(x.is_zero() for x in b[1:]):
            return tuple(b[0]*x for x in a)
        half = len(a)//2
        A, B, C, D = a[:half], a[half:], b[:half], b[half:]
        relation_a, relation_b = self.relations[level-1]
        AC, BD = self.multiply(A, C, level-1), self.multiply(B, D, level-1)
        middle = self.subtract(self.subtract(self.multiply(self.add(A, B), self.add(C, D), level-1), AC), BD)
        return self.add(AC, self.multiply(BD, relation_a, level-1))+self.add(middle, self.multiply(BD, relation_b, level-1))

    def inverse(self, coefficients, level):
        self.inversions += 1
        if level == 0:
            return (coefficients[0].inverse(),)
        half = len(coefficients)//2
        A, B = coefficients[:half], coefficients[half:]
        if all(x.is_zero() for x in B):
            return self.inverse(A, level-1)+(self.zero,)*half
        a, b = self.relations[level-1]
        conjugate_A = self.add(A, self.multiply(B, b, level-1))
        norm = self.subtract(self.multiply(A, conjugate_A, level-1), self.multiply(self.multiply(B, B, level-1), a, level-1))
        invnorm = self.inverse(norm, level-1)
        return self.multiply(conjugate_A, invnorm, level-1)+tuple(-x for x in self.multiply(B, invnorm, level-1))

    def norm(self, coefficients, level):
        if level == 0:
            return coefficients[0]
        half = len(coefficients)//2
        A, B = coefficients[:half], coefficients[half:]
        a, b = self.relations[level-1]
        lower = self.subtract(self.multiply(A, self.add(A, self.multiply(B, b, level-1)), level-1),
                              self.multiply(self.multiply(B, B, level-1), a, level-1))
        return self.norm(lower, level-1)


class Element:
    def __init__(self, algebra, level, coefficients):
        self.algebra, self.level, self.coefficients = algebra, level, coefficients
        self.cached_inverse = None

    def cast(self, other):
        if isinstance(other, Element):
            if other.algebra is not self.algebra:
                raise ValueError("混用了两个二次塔")
            return other
        return self.algebra.scalar(other)

    def promote(self, level):
        if level < self.level:
            raise ValueError("不能直接降去二次生成元")
        return self if level == self.level else Element(self.algebra, level, self.coefficients+(self.algebra.zero,)*(2**level-len(self.coefficients)))

    def pair(self, other):
        other = self.cast(other)
        level = max(self.level, other.level)
        return self.promote(level), other.promote(level), level

    def is_zero(self):
        return all(x.is_zero() for x in self.coefficients)

    def __add__(self, other):
        a, b, level = self.pair(other)
        return Element(self.algebra, level, self.algebra.add(a.coefficients, b.coefficients))

    __radd__ = __add__

    def __neg__(self):
        return Element(self.algebra, self.level, tuple(-x for x in self.coefficients))

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        a, b, level = self.pair(other)
        return Element(self.algebra, level, self.algebra.multiply(a.coefficients, b.coefficients, level))

    __rmul__ = __mul__

    def inverse(self):
        if self.cached_inverse is None:
            self.cached_inverse = Element(self.algebra, self.level, self.algebra.inverse(self.coefficients, self.level))
        return self.cached_inverse

    def __truediv__(self, other):
        return self*self.cast(other).inverse()

    def __rtruediv__(self, other):
        return self.cast(other)*self.inverse()

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inverse()**(-exponent)
        result, value = self.cast(1), self
        while exponent:
            if exponent & 1:
                result *= value
            exponent >>= 1
            if exponent:
                value *= value
        return result

    def encode(self, level=None):
        element = self.promote(level) if level is not None else self
        return [x.encode() for x in element.coefficients]
