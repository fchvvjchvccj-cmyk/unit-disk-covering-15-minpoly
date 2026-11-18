"""从原特征零表达式编译临界方程与雅可比，共享相同的乘积、逆和高阶偏导。"""

from functools import lru_cache
from collections import Counter

from eliminate import endpoint, primitive_resultant


class Graph:
    def __init__(self):
        self.rows=[]
        self.nodes={}

    def node(self,op,*args):
        key=(op,*args)
        if key not in self.nodes:
            self.nodes[key]=len(self.rows)
            self.rows.append(key)
        return Symbol(self,self.nodes[key])

    def constant(self,c):
        return self.node("const",c)

    @lru_cache(maxsize=None)
    def derivative(self,i,j):
        op,*args=self.rows[i]
        wrap=lambda k:Symbol(self,k)
        der=lambda k:self.derivative(k,j)
        if op=="const":
            return self.constant(0)
        if op=="var":
            return self.constant(int(args[0]==j))
        a=args[0]
        if op=="neg":
            return -der(a)
        if op=="inv":
            return -der(a)*wrap(i)*wrap(i)
        b=args[1]
        if op=="add":
            return der(a)+der(b)
        return der(a)*wrap(b)+wrap(a)*der(b)

    def program(self,outputs):
        active=set()
        def visit(i):
            if i in active:
                return
            active.add(i)
            op,*args=self.rows[i]
            if op not in ("var","const"):
                for j in args:
                    visit(j)
        for v in outputs:
            visit(v.i)
        numbers={j:i for i,j in enumerate(sorted(active))}
        rows=[]
        for j in sorted(active):
            op,*args=self.rows[j]
            rows.append((op,*([numbers[k] for k in args] if op not in ("var","const") else args)))
        return rows,[numbers[v.i] for v in outputs]


class Symbol:
    def __init__(self,graph,i):
        self.graph,self.i=graph,i

    def coerce(self,v):
        return v if isinstance(v,Symbol) else self.graph.constant(v)

    def __add__(self,v):
        v=self.coerce(v)
        sa,sb=self.graph.rows[self.i],self.graph.rows[v.i]
        if sa==("const",0):
            return v
        if sb==("const",0):
            return self
        if sa[0]==sb[0]=="const":
            return self.graph.constant(sa[1]+sb[1])
        return self.graph.node("add",*sorted((self.i,v.i)))

    __radd__=__add__

    def __neg__(self):
        row=self.graph.rows[self.i]
        if row[0]=="const":
            return self.graph.constant(-row[1])
        if row[0]=="neg":
            return Symbol(self.graph,row[1])
        return self.graph.node("neg",self.i)

    def __sub__(self,v):
        return self+-self.coerce(v)

    def __rsub__(self,v):
        return (-self)+v

    def __mul__(self,v):
        v=self.coerce(v)
        sa,sb=self.graph.rows[self.i],self.graph.rows[v.i]
        if sa==("const",0) or sb==("const",1):
            return self
        if sb==("const",0) or sa==("const",1):
            return v
        if sa[0]==sb[0]=="const":
            return self.graph.constant(sa[1]*sb[1])
        if sa==("const",-1):
            return -v
        if sb==("const",-1):
            return -self
        return self.graph.node("mul",*sorted((self.i,v.i)))

    __rmul__=__mul__

    def __truediv__(self,v):
        return self*self.graph.node("inv",self.coerce(v).i)

    def __rtruediv__(self,v):
        return self.coerce(v)/self


def incidence(graph,variables):
    e,a,d,left,b=variables
    r,t,v,w,wbar=endpoint(e,a,d)
    k,j=w-v,v*wbar-1
    h=(w-left)*(left*wbar-1)*k*j-r*r*(k+left*j)*(k+left*j)
    a0,b0=-r*r*t*e,left*r*r*t
    c0,d0=-left*r*(1-r*r)*e,r*(1-r*r)
    k0=left*(1-r*r-r*t*e)-(1-r*r)*e+r*t
    z=a0-b0
    u=(b0+c0)*b*b+(k0+b0*t)*b+b0+d0
    q=u*(u+z*(b*b+t*b+1))+z*z*b*b
    dh=[graph.derivative(h.i,i) for i in range(4)]
    dq=[graph.derivative(q.i,i) for i in range(5)]
    cross=[dh[i]*dq[3]-dq[i]*dh[3] for i in range(3)]
    return [h,q,dq[4],cross[2],cross[0]*graph.derivative(r.i,1)-cross[1]*graph.derivative(r.i,0)]


@lru_cache(maxsize=4)
def compile_system(jacobian,with_anchors=False):
    graph=Graph()
    variables=[graph.node("var",j) for j in range(5 if with_anchors else 3)]
    if with_anchors:
        outputs=incidence(graph,variables)
    else:
        e,a,d=variables
        r=(e-1/e)/(a-1/a)
        f=primitive_resultant(e,a,d)
        fe,fa,fd=(graph.derivative(f.i,j) for j in range(3))
        cross=fe*graph.derivative(r.i,1)-fa*graph.derivative(r.i,0)
        outputs=[f,fd,cross]
    if jacobian:
        outputs += [graph.derivative(q.i,j) for q in outputs for j in range(len(variables))]
    return graph.program(outputs)


@lru_cache(maxsize=4)
def kernel(jacobian,with_anchors):
    tape,indices=compile_system(jacobian,with_anchors)
    uses=Counter(j for op,*args in tape if op not in ("var","const") for j in args)
    uses.update(indices)
    slots,free={},[]
    size=0
    code=["def execute(values):"]
    for i,(op,*args) in enumerate(tape):
        names=[f"v{slots[j]}" for j in args] if op not in ("var","const") else []
        if op=="var":
            expr=f"values[{args[0]}]"
        elif op=="const":
            expr=str(args[0])
        elif op=="neg":
            expr="-"+names[0]
        elif op=="inv":
            expr="1/"+names[0]
        else:
            expr=names[0]+("+" if op=="add" else "*")+names[1]
        if names:
            for j in args:
                uses[j]-=1
                if uses[j]==0:
                    free.append(slots[j])
        if free:
            slots[i]=free.pop()
        else:
            slots[i]=size
            size+=1
        code.append(f" v{slots[i]}={expr}")
    code.append(" return ["+",".join(f"v{slots[j]}" for j in indices)+"]")
    namespace={}
    exec("\n".join(code),namespace)
    return namespace["execute"]


def evaluate(values,jacobian=False):
    return kernel(jacobian,len(values)==5)(values)
