"""Steady S12 material/order balances and an axial redistribution in a 4D receiver."""
import sympy as sp


def check(label, actual, target):
    actual = sp.factor(sp.simplify(actual))
    residual = sp.simplify(actual-target)
    print(f"{label}: {actual} | expected {target} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def main():
    x = sp.symbols("x y z w", real=True)
    r2 = sum(t*t for t in x)
    b = sp.Symbol("B", real=True)
    nb, mass, radius = sp.symbols("n_b m R", positive=True)
    chi = sp.Function("chi")(*x)
    potential = b*x[0]/r2**2
    dj = [nb*sp.diff(potential, t) for t in x]
    div = lambda v: sum(sp.diff(component,t) for component,t in zip(v,x))
    check("change in material divergence away from mouth", div(dj), 0)
    compensation = [-chi*component for component in dj]
    # Compare separately differentiated advected and nonadvected order currents.
    advection = div([chi*component for component in dj])
    check("change in n Gamma", advection+div(compensation), 0)
    print("order advection change at (1,0,0,0) for chi=x:",
          sp.factor(advection.subs(chi,x[0]).doit().subs(dict(zip(x,[1,0,0,0])))), flush=True)
    # Contract the Cartesian gradient with the radial normal before projecting.
    radial = sp.factor(sum(t*component for t,component in zip(x,dj))/sp.sqrt(r2))
    mu = sp.Symbol("mu", real=True)
    radial_mass = sp.simplify(mass*radial.subs(
        {x[0]:radius*mu,x[1]:radius*sp.sqrt(1-mu**2),x[2]:0,x[3]:0},simultaneous=True))
    weight = 4*sp.pi*radius**3*sp.sqrt(1-mu**2)
    front = sp.integrate(radial_mass*weight,(mu,0,1))
    back = sp.integrate(radial_mass*weight,(mu,-1,0))
    check("facing hemisphere mass-rate change", front, -4*sp.pi*mass*nb*b/radius)
    check("away hemisphere mass-rate change", back, 4*sp.pi*mass*nb*b/radius)
    check("total mass-rate change", front+back, 0)
    # Distinct angular distributions can even have identical half totals.
    profile = mu+sp.Rational(5,4)*sp.chebyshevu(3,mu)
    for name,lo,hi in (("facing",0,1),("away",-1,0)):
        check(f"mixed profile {name} integral", sp.integrate(profile*weight,(mu,lo,hi)),0)
    area = 2*sp.pi**2*radius**3
    check("mixed profile dipole projection", 4*sp.integrate(mu*profile*weight,(mu,-1,1))/area,1)


if __name__=="__main__":
    main()
