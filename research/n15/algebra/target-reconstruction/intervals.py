"""固定二进制分母的向外舍入区间；全部接受比较只使用整数。"""

from fractions import Fraction


class Dyadic:
    def __init__(self, bits):
        self.bits = bits
        self.scale = 1 << bits
        self.one = (self.scale, self.scale)
        self.zero = (0, 0)

    def rational(self, value):
        value = Fraction(value)
        numerator = value.numerator * self.scale
        return (numerator // value.denominator,
                -((-numerator) // value.denominator))

    @staticmethod
    def add(a, b):
        return a[0] + b[0], a[1] + b[1]

    @staticmethod
    def neg(a):
        return -a[1], -a[0]

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    @staticmethod
    def integer(a, k):
        return (a[0] * k, a[1] * k) if k >= 0 else (a[1] * k, a[0] * k)

    def mul(self, a, b):
        products = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
        return min(products) >> self.bits, -((-max(products)) >> self.bits)

    def div(self, a, b):
        if b[0] <= 0 <= b[1]:
            raise ValueError("区间分母包含零")
        if b[1] < 0:
            return self.div(self.neg(a), self.neg(b))
        products = [(x*self.scale, y) for x in a for y in b]
        return (min(x//y for x, y in products),
                max(-((-x)//y) for x, y in products))

    def power(self, a, exponent):
        value = self.one
        for _ in range(exponent):
            value = self.mul(value, a)
        return value

    @staticmethod
    def absolute(a):
        return max(abs(a[0]), abs(a[1]))

    def polynomial(self, terms, values):
        powers = [[self.power(v, k) for k in range(7)] for v in values]
        result = self.zero
        for coefficient, exponent in terms:
            term = self.rational(coefficient)
            for index, degree in exponent:
                term = self.mul(term, powers[index][degree])
            result = self.add(result, term)
        return result


def derivatives(terms, dimension):
    output = []
    for polynomial in terms:
        columns = [[] for _ in range(dimension)]
        for coefficient, powers in polynomial:
            for column, degree in powers:
                exponent = tuple((i, e-(i == column)) for i, e in powers
                                 if e-(i == column))
                columns[column].append((coefficient*degree, exponent))
        output.append(columns)
    return output


def evaluate_rows(arithmetic, rows, values):
    powers = [[arithmetic.power(v, k) for k in range(7)] for v in values]
    output = []
    for row in rows:
        result = arithmetic.zero
        for coefficient, exponent in row:
            term = arithmetic.rational(coefficient)
            for index, degree in exponent:
                term = arithmetic.mul(term, powers[index][degree])
            result = arithmetic.add(result, term)
        output.append(result)
    return output
