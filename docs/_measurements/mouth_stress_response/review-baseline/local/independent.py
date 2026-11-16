"""Cartesian polynomial witnesses on the unit three-ball. Set hbar^2/(4m)=1."""
import sympy as sp

x, s, r, u, n = sp.symbols("x s r u n", real=True)
d = sp.Rational(1, 10)


def dx(f):
    # s=x^2+y^2+z^2.
    return sp.diff(f, x) + 2*x*sp.diff(f, s)


def volume(f):
    # Exact monomial integrals divided by 4 pi. Average (x/r)^i on S^2.
    return sum(c / ((i+1)*(i+2*j+3)) for (i,j),c in sp.Poly(f,x,s).terms()
               if i % 2 == 0)


def linear(f):
    return sp.diff(f, u).subs(u, 0)


def cart_force(n_field, b_field):
    return volume(linear(n_field*dx(b_field)))


def radial(n_field, b_field):
    # Radial dipole amplitudes read on +x, where mu=1.
    def axis(f):
        return f.subs({s:r**2, x:r})
    n0 = axis(n_field.subs(u, 0))
    b0 = axis(b_field.subs(u, 0))
    n1 = axis(linear(n_field))
    b1 = axis(linear(b_field))
    rim = (r**2*n0*b1).subs(r, 1)/3
    interior = sp.integrate(r**2*(n1*sp.diff(b0,r)-b1*sp.diff(n0,r)),(r,0,1))/3
    return rim, interior


def value(label, actual):
    print(label, sp.factor(actual), flush=True)


den = 1+s+u*x
general = 3+s**2+u*x*(2-s)
rim, interior = radial(den, general)
direct = cart_force(den, general)
open_endpoint = rim
value("open ball Cartesian dF/(4*pi)", direct)
value("explicit endpoint", rim)
value("radial interior", interior)
value("difference", direct-rim-interior)

# Conditionally local laws. The second has indistinguishable rest calibrations.
A = 1+s
B = A+d
weight = (1-s)**3
base = weight*(s+n**2)
altered = base+weight*(n-A)*(n-B)
results = []
checks = [direct-rim-interior]
for name, law in [("base",base), ("altered",altered)]:
    measured = law.subs(n,den)
    rim, interior = radial(den, measured)
    Cartesian = cart_force(den, measured)
    # e_n = b in chosen units. If e=e(s,n), then at fixed field
    # d_x e = 2x e_s and its linear response is 2x e_sn d_U n.
    energy_slope = volume(2*x**2*sp.diff(law,s).subs(n,A))
    value(name+" Cartesian", Cartesian)
    value(name+" fixed-field translation of local energy", energy_slope)
    value(name+" radial endpoint", rim)
    checks.extend([Cartesian-energy_slope, rim, Cartesian-rim-interior])
    results.append(Cartesian)
for state in (A,B):
    equal = (altered-base).subs(n,state)
    value("law difference at rest n="+str(state),equal)
    checks.append(equal)
secant_error = ((base.subs(n,B)-base.subs(n,A))/(B-A)-sp.diff(base,n).subs(n,A))
value("base two-state secant minus its derivative at A",secant_error)
value("altered minus base susceptibility at A", sp.diff(altered-base,n).subs(n,A))
value("altered minus base force slope despite identical two rest traces", results[1]-results[0])
assert all(sp.simplify(a)==0 for a in checks)
assert open_endpoint != 0 and results[1]!=results[0] and sp.simplify(secant_error)!=0
