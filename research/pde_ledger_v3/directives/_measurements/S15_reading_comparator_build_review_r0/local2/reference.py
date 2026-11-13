"""Independent contractions and spherical polynomial integration for moving S15 surfaces."""
import sympy as s

n = s.symbols('n1:4')
N = s.Matrix(n)
eps = s.symbols('eps')
X = s.Matrix([2, -1, 3])
C = s.Matrix([2, 3, -1])
x = s.Matrix(s.symbols('x1:4'))
v = s.Matrix([1 + x[1], 2 + x[2], 3 + x[0]])
rho = 5 + x[0] - 2*x[2]
grad = v.jacobian(x)
I = s.eye(3)
viscous = (2+x[0])*(grad+grad.T-s.Rational(2,3)*grad.trace()*I)+(3+x[1])*grad.trace()*I
other = s.Matrix(3, 3, lambda i,j: (i+2)*x[j] + (2*j+1)*x[i] + (i-j)*x[(j+1)%3])
pieces = {
    'carried': rho*v*v.T,
    'pressure': rho**2*I,
    'viscous': -viscous,
    'other': other
}
G = rho*v

def mean(expr):
    total = 0
    for powers, coeff in s.Poly(expr, *n).terms():
        if any(p % 2 for p in powers):
            continue
        total += coeff * s.prod(s.factorial2(p-1) for p in powers)/s.factorial2(sum(powers)+1)
    return total

def integrate(vec, b):
    return tuple(mean(e.subs(dict(zip(x, X+b*N))))*b**2 for e in vec)

def flux(b, rate, sigma):
    # Values are divided by 4*pi. The chart contributes b^2, outward n=sigma*N.
    out = {k:integrate(-M*(sigma*N), b) for k, M in pieces.items()}
    out['centre'] = integrate(G*(C.dot(sigma*N)), b)
    out['radius'] = integrate(G*rate*sigma*N.dot(N), b)
    out['rest'] = tuple(map(sum, zip(*(out[k] for k in pieces))))
    out['transport'] = tuple(map(sum,zip(out['centre'],out['radius'])))
    out['crossing'] = tuple(map(sum,zip(out['rest'],out['transport'])))
    return out

shells = [flux(1,s.Rational(1,2),-1), flux(2,s.Rational(-2,3),1)]
for k in shells[0]:
    print('shell / (4*pi)', k, 'in',shells[0][k],'out',shells[1][k],'both',tuple(map(sum,zip(*(f[k] for f in shells)))))
print('one-sided other[0,1]+=3*x2; sign and measure', [
    integrate(-s.Matrix([3*x[1]*ori*n[1],0,0]),b) for b,ori in [(1,-1),(2,1)]])
# The same admissible positive radial interval controls all powers, including log.
d,b=s.symbols('d b',positive=True)
a=s.symbols('a',real=True)
z=s.symbols('z',positive=True)
# r^2=d^2+b^2-2*d*b*cos(theta), d cos(theta)=-r dr/(d*b).
print('family scalar shell integral/(2*pi*b^2)', 'generic', s.integrate(z**(a+1),z,conds='none'), '/(d*b) from d-b to d+b', 'a=-2', s.integrate(z**(-1),z))
# Independent separation variation for an analytic field containing the moving radius and
# a free response: X_1=(5,0,0), X_2=(0,0,0), x=(0,2,0). d>R=3>r_in=1.
r=s.sqrt((5*(1+eps))**2+4)
Q=(2+3*eps)*r**2 + (7-5*eps)**2
Qgeo=2*r**2+7**2
geo=Qgeo.diff(eps).subs(eps,0)+2*Q.subs(eps,0)
resp=(Q-Qgeo).diff(eps).subs(eps,0)
print('test Q,dQ,geo,resp,partition',Q.subs(eps,0),Q.diff(eps).subs(eps,0),geo,resp,(Q.diff(eps)+2*Q).subs(eps,0)-geo-resp)
T1,T2,T11,T22=s.symbols('v1 v2 acc1 acc2')
moment=3*T1*T2+5*T1*T11+7*T2*T22+11*T11*T22
print('independent velocity cuts', [s.expand(moment-moment.subs(k,0)) for k in (T1,T2)])
print('density,length,time exp', 'surface rest/transport/geo/resp/velocity-part', [1,4,-2], 'integrated momentum',[1,4,-1])
