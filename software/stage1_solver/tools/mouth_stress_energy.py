"""Energy current and boundary-work limits of the supplied neutral GNLS sector.

No complete live brane constitutive law is assumed. Bulk identities are exact
for the neutral action with general U(n), including an explicit V_conf(x,t).
The open interface is not assigned a boundary condition or energy source.
"""
import sympy as sp


def emit(label, value):
    print(f"{label}: {sp.factor(value)}", flush=True)


def check(label, lhs, rhs):
    residual = sp.simplify(lhs-rhs)
    print(f"{label}: {sp.factor(lhs)} | {sp.factor(rhs)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def carrier():
    t=sp.Symbol("t",real=True)
    x=sp.symbols("x0:4",real=True)
    hbar,m=sp.symbols("hbar m",positive=True)
    c=hbar**2/(2*m)
    a,b,vc=(sp.Function(name)(t,*x) for name in ("a","b","V"))
    n=a*a+b*b
    u=sp.Function("U")
    nvar=sp.Symbol("nvar",positive=True)
    chemical=sp.diff(u(nvar),nvar).subs(nvar,n)
    energy=c*sum(sp.diff(a,z)**2+sp.diff(b,z)**2 for z in x)+vc*n+u(n)
    current=[-2*c*(sp.diff(a,t)*sp.diff(a,z)+sp.diff(b,t)*sp.diff(b,z)) for z in x]
    lap_a=sum(sp.diff(a,z,2) for z in x)
    lap_b=sum(sp.diff(b,z,2) for z in x)
    rates={sp.diff(a,t):(-c*lap_b+(vc+chemical)*b)/hbar,
           sp.diff(b,t):(c*lap_a-(vc+chemical)*a)/hbar}
    balance=sp.expand(sp.diff(energy,t)+sum(sp.diff(current[i],z) for i,z in enumerate(x)))
    on_shell=sp.simplify(balance.subs(rates))
    emit("neutral carrier energy density",energy)
    emit("neutral carrier energy current, component 0",current[0])
    check("neutral carrier local energy source",on_shell,n*sp.diff(vc,t))

    # Independently put the complex current in amplitude/phase variables.
    number,nt,nx,theta,tt,tx=sp.symbols("n n_t n_x theta theta_t theta_x",real=True)
    field=sp.sqrt(number)*sp.exp(sp.I*theta)
    field_bar=sp.sqrt(number)*sp.exp(-sp.I*theta)
    psi_t=sp.diff(field,number)*nt+sp.diff(field,theta)*tt
    bar_t=sp.diff(field_bar,number)*nt+sp.diff(field_bar,theta)*tt
    psi_x=sp.diff(field,number)*nx+sp.diff(field,theta)*tx
    bar_x=sp.diff(field_bar,number)*nx+sp.diff(field_bar,theta)*tx
    flux=sp.simplify(-c*(bar_t*psi_x+psi_t*bar_x))
    emit("polar energy current, one component",flux)
    mu,v,wn=sp.symbols("mu v W_N",real=True)
    stationary=sp.factor(flux.subs({nt:0,tt:-mu/hbar,tx:m*v/hbar}))
    check("stationary energy per number current",stationary,mu*number*v)
    # The face speed subtracts energy transport; it does not replace the
    # material velocity in the pressure/stress work part of the current.
    en=sp.Symbol("E_density",real=True)
    emit("moving normal-cut energy current",flux-en*wn)
    number_rate,mu0,speed=sp.symbols("J_N mu0 U",real=True)
    controlled_power=number_rate*(mu0+speed**2/2)
    emit("implemented held-feed power increment",controlled_power-controlled_power.subs(speed,0))


def stress_and_profile():
    # Stress power requires the actual application velocity, not just the
    # integrated force. The two-patch witness is algebraic, not a mouth law.
    t1,t2,u1,u2,vm=sp.symbols("t1 t2 u1 u2 V",real=True)
    work=t1*u1+t2*u2
    force=t1+t2
    relative=t1*(u1-vm)+t2*(u2-vm)
    check("stress power split",work,vm*force+relative)
    traction,common,slip=sp.symbols("T u s",real=True)
    witness={t1:traction,t2:-traction,u1:common+slip,u2:common-slip}
    emit("zero-resultant witness, force",force.subs(witness))
    emit("zero-resultant witness, power",work.subs(witness))

    # For a prescribed translating confinement potential, the local energy
    # source is paired with the opposite reaction on its geometric coordinate.
    t,x=sp.symbols("t x",real=True)
    center=sp.Function("X")(t)
    potential=sp.Function("Vconf")(x-center)
    number=sp.Symbol("n",positive=True)
    local_source=number*sp.diff(potential,t)
    force_on_profile=number*sp.diff(potential,x)
    check("translated-profile energy transfer",local_source,-sp.diff(center,t)*force_on_profile)


if __name__=="__main__":
    print("SymPy version =",sp.__version__,flush=True)
    carrier()
    stress_and_profile()
