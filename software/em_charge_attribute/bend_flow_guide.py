"""A two-cut ideal-flow boundary example on the sharp and annular bends."""

import mpmath as mp
import sympy as sp


def verify(label, expr, expected):
    value = sp.factor(expr)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def guide_work():
    area, ratio, rho, rate, ell = sp.symbols("A m rho D ell", positive=True)
    vi = rate/(rho*area)
    vo = vi/ratio
    pin = rho*(vo**2-vi**2)/2  # outlet and surrounds at zero gauge pressure
    ti = rate*vi + pin*area
    to = rate*vo
    f0 = rate**2/(rho*area)
    verify("Bernoulli plus pressure, inlet thrust per reference thrust",ti/f0,(1+ratio**-2)/2)
    verify("outlet thrust per reference thrust",to/f0,1/ratio)
    verify("normal force per reference thrust",(ti-to)/f0,(1-1/ratio)**2/2)
    # Z=ell h at a sharp mouth. For local inlet and axial outlet the same
    # throughflow supplies both the bend source and its translation work.
    g, dx, dh = sp.symbols("g dx dh", real=True)
    work = (-ti*ell*g)*dx + (ti-to)*ell*(dh+g*dx)
    verify("sharp emitted load, per ell times reference thrust",work.diff(dh)/(ell*f0),
           (1-1/ratio)**2/2)
    verify("sharp translation response, same units",work.diff(dx)/(ell*g*f0),-1/ratio)
    verify("sharp signed response/emission",work.diff(dx)/(g*work.diff(dh)),
           -2*ratio/(ratio-1)**2)
    verify("sharp both cuts follow bend",((ti-to)*ell*g + ell*g*(to-ti)),0)
    verify("head and kinetic energy flux difference at rest",
           rate*(vi**2/2+pin/rho-vo**2/2),0)
    # Moving cut, outward directions referred to each stream's own positive
    # axis. rx=U dot n_in, sx=U dot n_out; pA U.n is part of the power.
    rx,sx,U2 = sp.symbols("U_i U_o U2",real=True)
    power = rate*((U2+2*vi*rx+vi**2)-(U2+2*vo*sx+vo**2))/2 + pin*area*(vi+rx)
    verify("moving-port power minus complete thrust work",power-(ti*rx-to*sx),0)


def ring(a):
    """Exact radial basis at high precision for l=0,1 with shell matching."""
    a = mp.mpf(a)
    b = a + mp.mpf(1)/20
    span = b**3-a**3
    kappa, km = mp.mpf(7)/8, mp.mpf(1)
    beta = mp.sqrt(3*km/(4*mp.pi*span*kappa))
    # Unit-source monopole and its actual values at both edges, needed
    # to cross-check the localized restoring-force contribution.
    c = 1/(mp.cosh(beta*(b-a))+a*beta*mp.sinh(beta*(b-a)))
    def h0(r):
        u = c*(a*mp.cosh(beta*(r-a))+mp.sinh(beta*(r-a))/beta)
        return (1-u/r)/km
    z = 4*mp.pi*kappa*b*h0(b)
    # Dipole basis obtained by differentiating exp(+-beta*r)/r.
    def basis(r,t):
        return mp.exp(t*r)*(t*r-1)/r**2
    def derivative(r,t):
        return mp.exp(t*r)*(r*r*t*t-2*r*t+2)/r**3
    def fundamental(r):
        return mp.matrix([[basis(r,beta),basis(r,-beta)],[derivative(r,beta),derivative(r,-beta)]])
    final = fundamental(b)*mp.lu_solve(fundamental(a),mp.matrix([a,1]))
    norm = 3/(final[1]+2*final[0]/b)  # outside f=r+dip/r^2
    fa, fb = norm*a, norm*final[0]
    gain = (b**2*fb-a**2*fa)/span
    restoring = -km*(b**2*h0(b)*fb-a**2*h0(a)*fa)/span
    print(f"radius {int(a)}, emission z={float(z):.12f}, slope gain u={float(gain):.12f},"
          f" passive contact correction={float(restoring):.12f}, u/z={float(gain/z):.12f}")
    if abs(gain+restoring-z) > mp.mpf("1e-44"):
        raise ArithmeticError("direct annulus force disagrees with exterior stress")
    return z,gain


if __name__ == "__main__":
    guide_work()
    mp.mp.dps = 55
    data = [ring(a) for a in (1,2)]
    # Both physical mouths use the same widening ratio. Compare all four
    # ordered forces even after choosing the intake ratio to equalize L.
    comparison = [1-5*u/z for z,u in data]  # m=2; Ti/(Ti-To)=5
    print("area ratio 2, B/L at radii 1 and 2 =", [float(v) for v in comparison])
    d_ratio = mp.sqrt(8*data[0][0]/data[1][0])
    print("mass-rate doubling ratio required to equalize L =", float(d_ratio))
    print("equal L, ordered mixed-force magnitude ratio target small/target large =",
          float(comparison[0]/comparison[1]))
    if not all(v<0 for v in comparison) or comparison[0]==comparison[1]:
        raise ArithmeticError("signed numerical control")
