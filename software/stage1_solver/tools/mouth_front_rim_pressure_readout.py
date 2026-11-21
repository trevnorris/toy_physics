"""Recover a pressure-jump dipole from a tension-supported rim with mixed harmonics."""
import sympy as sp


r, radius, sigma = sp.symbols('r R sigma', positive=True)
mu, alpha, beta = sp.symbols('mu alpha beta', real=True)
pressure = 0
shape = 0
for ell, amplitude in ((1, alpha), (3, beta)):
    pattern = sp.legendre(ell, mu)
    load = amplitude*(r/radius)**ell*pattern
    # Solve for the coefficient of the regular polynomial satisfying h(R)=0.
    k = sp.Symbol('k')
    trial = k*(r/radius)**ell*(radius**2-r**2)*pattern
    lap = (sp.diff(r**2*trial.diff(r),r)+sp.diff((1-mu**2)*trial.diff(mu),mu))/r**2
    coefficient, = sp.solve(sp.factor((-sigma*lap-load)/((r/radius)**ell*pattern)),k)
    profile = trial.subs(k,coefficient)
    assert profile.subs(r,radius)==0 and sp.simplify(-sigma*lap.subs(k,coefficient)-load)==0
    pressure += load
    shape += profile

inward_slope = -shape.diff(r).subs(r,radius)
s1 = sp.factor(3*sp.integrate(mu*inward_slope,(mu,-1,1))/2)
endpoints = sp.factor((inward_slope.subs(mu,1)-inward_slope.subs(mu,-1))/2)
moment = sp.factor(2*sp.pi*sp.integrate(r**3*mu*pressure,(mu,-1,1),(r,0,radius)))
recovered = sp.factor(sigma*4*sp.pi*radius**3*s1/3)
print('angular dipole of inward slope =',s1)
print('near-minus-far inward slope / 2 =',endpoints)
print('pressure first moment, reconstructed from rim =',moment,recovered)
assert sp.simplify(moment-recovered)==0 and sp.simplify(endpoints-s1-beta*radius/(9*sigma))==0
assert moment.subs(alpha,0)==0 and endpoints.subs(alpha,0)!=0
