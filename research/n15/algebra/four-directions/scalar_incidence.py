"""在以 A 为参照的实半角坐标中写出两个有理闭环，保留参数 r。"""

from critical_circuit import Graph


def constraints(v):
    r,e,z,d,l,b=v
    def unit(h):
        square=z*h*h
        return (1+square)/(1-square),2*h/(1-square)
    ax,ay=unit(1)
    dx,dy=unit(d)
    bx,by=unit(b)
    u=r*(1-z*e*e)-e*(1-z)
    transfer=z*e*ay-ax-(2*dx-1)
    twice_cx=2*bx+transfer
    kappa=1-r*r
    mx=2*kappa+r*r*transfer*twice_cx
    nx=-2*r*(transfer+kappa*bx)
    ny=2*r*kappa*by
    coefficient=r*r*transfer*(z*l*e-1)
    constant=l*(mx+nx-z*ny*e)-ny-(mx-nx)*e
    q=(twice_cx*twice_cx-4)*coefficient*coefficient-z*constant*constant
    right=r*dy/(1+r*dx)
    wx=1+r*(2*dx-1+ax)
    wy=r*ay
    base=wx*wx-z*wy*wy+1
    def chord(h):
        return base-2*wx+4*z*wy*h-z*(base+2*wx)*h*h
    line=(1+z*l*right)*wx-z*(l+right)*wy+z*l*right-1
    h=chord(l)*chord(right)-4*r*r*line*line
    return u,h,q


def integer_constraints(v):
    r,e,z,d,l,b=v
    den_a,den_d,den_b=1-z,1-z*d*d,1-z*b*b
    dx,bx=1+z*d*d,1+z*b*b
    tden=den_a*den_d
    shift=2*dx-den_d
    tn=(2*z*e-(1+z))*den_d-shift*den_a
    cn=2*bx*tden+tn*den_b
    cd=den_b*tden
    kappa=1-r*r
    mx=2*kappa*den_b*tden*tden+r*r*tn*cn
    nx=-2*r*(tn*den_b+kappa*bx*tden)*tden
    ny=4*r*kappa*b*tden*tden
    kc=r*r*tn*(z*l*e-1)
    constant=l*(mx+nx-z*ny*e)-ny-(mx-nx)*e
    q=(cn*cn-4*cd*cd)*kc*kc-z*constant*constant
    wx=tden+r*(shift*den_a+(1+z)*den_d)
    wy=2*r*den_d
    base=wx*wx-z*wy*wy+tden*tden
    c0,c1,c2=base-2*wx*tden,4*z*wy*tden,-z*(base+2*wx*tden)
    rp,rq=2*r*d,den_d+r*dx
    left=c0+l*(c1+l*c2)
    right=c0*rq*rq+rp*(c1*rq+rp*c2)
    line=(rq+z*l*rp)*wx-z*(l*rq+rp)*wy+tden*(z*l*rp-rq)
    h=left*right-4*r*r*line*line*tden*tden
    return r*(1-z*e*e)-e*(1-z),h,q


def integer_program():
    graph=Graph()
    v=[graph.node("var",i) for i in range(6)]
    u,h,q=integer_constraints(v)
    de=lambda f,i:graph.derivative(f.i,i)
    hl,ql=de(h,4),de(q,4)
    tang=[de(h,i)*ql-de(q,i)*hl for i in (1,2,3)]
    return graph.program([u,h,q,de(q,5),tang[2],tang[1]*de(u,1)-tang[0]*de(u,2)])


def coordinates(algebra,state):
    from padic_radius import values,Number
    e,a,d,l,b=values(algebra,state)
    ref=(a-1)/(a+1)
    def half(v):
        return (v-1)/((v+1)*ref)
    points=[Number(algebra,algebra.r,algebra.zero),half(e),ref*ref,half(d),half(l),half(b)]
    if any(not v.b.is_zero() for v in points):
        raise ValueError("半角比值未落回实方向子环")
    return [v.a for v in points]
