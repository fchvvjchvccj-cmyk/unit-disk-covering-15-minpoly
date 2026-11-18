"""整数常数的共享有理电路；显式记录所有除法主元。"""


class Circuit:
    def __init__(self, names):
        self.names, self.nodes, self.memo, self.divisors = list(names), [], {}, set()

    def node(self, op, args=(), value=None):
        args = tuple(args)
        if op in ("add", "mul"):
            args = tuple(sorted(args))
        key = (op, args, value)
        if key not in self.memo:
            self.memo[key] = len(self.nodes)
            record = {"op": op}
            if args:
                record["args"] = list(args)
            if value is not None:
                record["value"] = value
            self.nodes.append(record)
        return Expression(self, self.memo[key])

    def constant(self, value):
        if not isinstance(value, int):
            raise TypeError("电路常数只接受整数")
        return self.node("const", value=value)

    def variable(self, index):
        return self.node("var", value=index)


class Expression:
    def __init__(self, circuit, index):
        self.circuit, self.index = circuit, index

    def cast(self, other):
        if isinstance(other, Expression):
            if other.circuit is not self.circuit:
                raise ValueError("混用了不同电路")
            return other
        return self.circuit.constant(other)

    def constant_value(self):
        node = self.circuit.nodes[self.index]
        return node["value"] if node["op"] == "const" else None

    def is_zero(self):
        return self.constant_value() == 0

    def __add__(self, other):
        other = self.cast(other)
        a, b = self.constant_value(), other.constant_value()
        if a is not None and b is not None:
            return self.cast(a+b)
        if a == 0:
            return other
        if b == 0:
            return self
        return self.circuit.node("add", (self.index, other.index))

    __radd__ = __add__

    def __neg__(self):
        value = self.constant_value()
        if value is not None:
            return self.cast(-value)
        node = self.circuit.nodes[self.index]
        return Expression(self.circuit, node["args"][0]) if node["op"] == "neg" else self.circuit.node("neg", (self.index,))

    def __sub__(self, other):
        other = self.cast(other)
        return self.cast(0) if self.index == other.index else self+-other

    def __rsub__(self, other):
        return self.cast(other)-self

    def __mul__(self, other):
        other = self.cast(other)
        a, b = self.constant_value(), other.constant_value()
        if a is not None and b is not None:
            return self.cast(a*b)
        if a == 0 or b == 0:
            return self.cast(0)
        if a == 1:
            return other
        if b == 1:
            return self
        return self.circuit.node("mul", (self.index, other.index))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        if other.is_zero():
            raise ZeroDivisionError("电路构造出现恒零分母")
        if other.constant_value() == 1:
            return self
        self.circuit.divisors.add(other.index)
        return self.circuit.node("div", (self.index, other.index))

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def inverse(self):
        return 1/self

    def __pow__(self, exponent):
        if exponent < 0:
            return 1/self**(-exponent)
        result, value = self.cast(1), self
        while exponent:
            if exponent & 1:
                result *= value
            exponent >>= 1
            if exponent:
                value *= value
        return result


def evaluate(data, values):
    registers = []
    for node in data["nodes"]:
        op = node["op"]
        args = [registers[j] for j in node.get("args", [])]
        if op == "const":
            value = values[0]*0+node["value"]
        elif op == "var":
            value = values[node["value"]]
        elif op == "add":
            value = args[0]+args[1]
        elif op == "mul":
            value = args[0]*args[1]
        elif op == "neg":
            value = -args[0]
        elif op == "div":
            value = args[0]/args[1]
        else:
            raise ValueError("未支持的有理电路指令")
        registers.append(value)
    return registers
