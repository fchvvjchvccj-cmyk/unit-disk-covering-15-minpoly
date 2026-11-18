"""以系数对保留一次二次扩张，避免展开下层参数多项式。"""


def zero(value):
    if hasattr(value, "gradient"):
        return value.value.is_zero() and all(c.is_zero() for c in value.gradient)
    return value.is_zero()


class Pair:
    def __init__(self, a, b, radicand):
        self.a, self.b, self.radicand = a, b, radicand

    def cast(self, value):
        if isinstance(value, Pair):
            if value.radicand is not self.radicand:
                if isinstance(self.a, Pair):
                    return Pair(self.a.cast(value), self.a*0, self.radicand)
                raise ValueError("混用了不同二次生成元")
            return value
        return Pair(self.a.cast(value), self.a*0, self.radicand)

    def __add__(self, value):
        value = self.cast(value)
        return Pair(self.a+value.a, self.b+value.b, self.radicand)

    __radd__ = __add__

    def __neg__(self):
        return Pair(-self.a, -self.b, self.radicand)

    def __sub__(self, value):
        return self+-self.cast(value)

    def __rsub__(self, value):
        return self.cast(value)+-self

    def __mul__(self, value):
        value = self.cast(value)
        if zero(self.b):
            return Pair(self.a*value.a, self.a*value.b, self.radicand)
        if zero(value.b):
            return Pair(self.a*value.a, self.b*value.a, self.radicand)
        return Pair(self.a*value.a+self.radicand*self.b*value.b,
                    self.a*value.b+self.b*value.a, self.radicand)

    __rmul__ = __mul__

    def inverse(self):
        if zero(self.b):
            return Pair(1/self.a, self.b, self.radicand)
        denominator = self.a*self.a-self.radicand*self.b*self.b
        return Pair(self.a/denominator, -self.b/denominator, self.radicand)

    def __truediv__(self, value):
        return self*self.cast(value).inverse()

    def __rtruediv__(self, value):
        return self.cast(value)*self.inverse()

    def evaluate(self, root):
        return self.a+self.b*root

    def norm(self):
        return self.a*self.a-self.radicand*self.b*self.b

    def is_zero(self):
        return zero(self.a) and zero(self.b)

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
