"""Native-source-first independent shell and velocity reference. No candidate imports."""
from itertools import product
import sympy as S

xyz = S.Matrix(S.symbols("x y z"))
normal = S.Matrix(S.symbols("nx ny nz"))
source = S.Matrix((3,-2,4))
target = S.Matrix((1,0,-1))
inner, outer = S.Rational(1,2), S.Rational(3,2)
v1, v2 = S.Matrix((2,-1,3)), S.Matrix((-1,4,2))
radial_rates = S.Rational(-2,5), S.Rational(3,7)

def mean(expr, coords):
    result = 0
    for powers, coefficient in S.Poly(expr,*coords).terms():
        if all(n%2==0 for n in powers):
            result += coefficient*S.prod(S.factorial2(n-1) for n in powers)/S.factorial2(sum(powers)+1)
    return result

def signed_surface(evaluand, radius, orientation):
    result = []
    for value in evaluand:
        pull = value.subs(dict(zip(xyz,target+radius*normal)), simultaneous=True)
        result.append(4*radius**2*orientation*mean(pull,list(normal)))
    return S.Matrix(result)     # units of pi

def shells():
    # Independent realization of radial profiles plus affine remainder, nonlinear p, variable viscosity.
    W = S.Matrix([[1,2,-1],[3,1,2],[-2,4,1]])*xyz + S.Matrix((2,3,5))
    V = 2*(xyz-source)/3-(xyz-target)/5+W
    rho = 2+(xyz-source).dot(xyz-source)/7 + 3+(xyz-target).dot(xyz-target)/11 + 5+xyz[0]-2*xyz[2]
    eta,zeta = 2+xyz[1], 3+xyz[0]
    grad = V.jacobian(xyz)
    visc = eta*(grad+grad.T-2*S.eye(3)*grad.trace()/3)+zeta*S.eye(3)*grad.trace()
    def tensor(offset):
        return S.Matrix(3,3,lambda i,j: offset + 3*i+j+1+
              sum(((1+offset+3*i+j)**(l+1))*x for l,x in enumerate(xyz)))
    other = tensor(0)
    delta = tensor(10)
    current = rho*V*V.T+rho**2*S.eye(3)-visc+other
    both = [S.zeros(3,1) for _ in range(3)]
    for b,ori,rate in zip((inner,outer),(-1,1),radial_rates):
        rest=signed_surface(-current*normal,b,ori)
        moving=signed_surface(rho*V*(v2.dot(normal)+rate),b,ori)
        response=signed_surface(-delta*normal,b,ori)
        for i,vec in enumerate([rest,moving,response]): both[i]+=vec
        print('radius orientation',b,ori,'rest/pi',list(rest),'transport/pi',list(moving),
              'crossing/pi',list(rest+moving),'other response/pi',list(response))
    print('net rest transport crossing other-response in pi',*[list(v) for v in (both[0],both[1],both[0]+both[1],both[2])])
    trace_rows = S.Matrix([sum(delta[i,j].diff(xyz[j]) for j in range(3)) for i in range(3)])
    print('independent divergence test',list(both[2]+S.Rational(4,3)*(outer**3-inner**3)*trace_rows))
    assert both[2] == -S.Rational(4,3)*(outer**3-inner**3)*trace_rows

def velocity():
    # rho1=r1^2, V1=x-X1, direct integral over translated shell, both centre jets retained.
    delta=target-source
    v=4*(outer**3-inner**3)/3
    m=4*(outer**5-inner**5)/15
    jac=(v*delta.dot(delta)+5*m)*S.eye(3)+2*v*delta*delta.T
    for mode,w1,w2,rates in (('moving',v1,v2,radial_rates),('target stationary',v1,S.zeros(3,1),radial_rates),
                            ('radii stationary',v1,v2,(0,0)),('target and radii stationary',v1,S.zeros(3,1),(0,0)),
                            ('source stationary',S.zeros(3,1),v2,radial_rates)):
        boundary=sum((4*r**2*o*b*(delta.dot(delta)+S.Rational(5,3)*r**2) for r,o,b in zip((inner,outer),(-1,1),rates)),S.Integer(0))*delta
        print(mode,'rate/pi', list(jac*(w2-w1)+boundary),'X1 part/pi',list(-jac*w1),'X2 part/pi',list(jac*w2))
    t=S.symbols('t')
    print('stationary centre/radius reference jets, value derivative acceleration',[[S.diff(3+i+(i+2)*t**2,t,k).subs(t,0) for k in (0,1,2)] for i in range(8)])
    print('force dimensions density,length,time',(1,4,-2),'momentum',(1,4,-1),'vector zero',[(0,0,0)],'empty set',[])
    d,b=S.symbols('d b',positive=True)
    # Integral_{-1}^1 (d^2+b^2-2dbz)^(-1) dz; 2pi azimuth.
    print('exponent -2 angular integral',2*S.pi*S.log((d+b)/(d-b))/(d*b),'samples d=5 b=2',
          (S.pi*S.log(S.Rational(7,3))/5).evalf(36))
    print('both signed log intervals orientation +,-',*[S.log(1+q*S.Rational(3,5))-S.log(1+q*S.Rational(1,5)) for q in (1,-1)])
    print('full azimuth sine sine, sine cosine',S.integrate(S.sin(t)**2,(t,0,2*S.pi)),S.integrate(S.sin(t)*S.cos(t),(t,0,2*S.pi)))
if __name__=="__main__":
    shells()
    velocity()
