"""Scoped energy consequences of the existing real-fraction order model.

Homogeneous parcel diagnostics are not throat solutions. The pair's required
power is an input from the separate ideal monopole benchmark. No conversion
kinetics, transmission efficiency, or physical energy source is selected.
"""
import sympy as sp


def emit(label, value):
    print(f"{label}: {sp.factor(value)}", flush=True)


def check(label, lhs, rhs):
    residual = sp.factor(lhs-rhs)
    print(f"{label}: {sp.factor(lhs)} | {sp.factor(rhs)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def order_witness():
    chi = sp.Symbol("chi", real=True)
    a, loss = sp.symbols("a gamma", positive=True)
    well = a*chi**2*(1-chi)**2
    mu = sp.diff(well, chi)
    # A homogeneous instantaneous balance with J_chi=0 and Gamma_bal=-gamma.
    # Same imposed signed order-loss rate at two states, with no kinetic law.
    power = -loss*mu
    emit("bare order derivative", mu)
    emit("same negative balance rate, power at chi=3/4", power.subs(chi,sp.Rational(3,4)))
    emit("same negative balance rate, power at chi=1/4", power.subs(chi,sp.Rational(1,4)))
    emit("bare well barrier above either endpoint", well.subs(chi,sp.Rational(1,2))-well.subs(chi,1))
    integrated = sp.integrate(mu,(chi,1,0))
    check("integrated bare order change",integrated,well.subs(chi,0)-well.subs(chi,1))

    # Shear is a separate live state: partial order work is only one term in
    # the total energy derivative along a path where the shear energy changes.
    shear = sp.Function("S")(chi)
    extra = sp.Function("G")(chi)
    functional = well+chi*shear+extra
    fixed_state_order_work = mu+shear+sp.diff(extra,chi)
    check("additional live shear-path work",sp.diff(functional,chi)-fixed_state_order_work,
          chi*sp.diff(shear,chi))
    emit("complete homogeneous endpoint drop with live shear and extra energy",
         functional.subs(chi,1)-functional.subs(chi,0))

    # The conserved-number energy reference stays symbolic, on the n=5 EOS.
    number = sp.Symbol("n", positive=True)
    nf, nb, eos = sp.symbols("n_feed n_bulk K", positive=True)
    reference = sp.Symbol("C_ref", real=True)
    density_energy = eos*number**5/4+reference*number
    pressure = sp.expand(number*sp.diff(density_energy,number)-density_energy)
    enthalpy = sp.factor((density_energy+pressure)/number)
    head_difference = sp.factor(enthalpy.subs(number,nf)-enthalpy.subs(number,nb))
    emit("n=5 constituent enthalpy",enthalpy)
    emit("n=5 constituent enthalpy difference",head_difference)
    check("common conserved-number reference cancellation",sp.diff(head_difference,reference),0)
    check("equal-density bare constituent head difference",head_difference.subs(nb,nf),0)


def pair_requirement():
    d1,d2,sep,rho3,rho4=sp.symbols("D1 D2 d rho3 rho4",positive=True)
    v1,v2=sp.symbols("V1 V2",real=True)
    # Supplied ideal common-domain load, not a result for physical mouths.
    attraction=d1*d2*(1/(4*sp.pi*rho3*sep**2)+1/(2*sp.pi**2*rho4*sep**3))
    required_per_mouth=sp.factor(attraction*(v1-v2))
    emit("required pair-minus-alone net energy per mass, mouth 1",required_per_mouth/d1)
    emit("required pair-minus-alone net energy per mass, mouth 2",required_per_mouth/d2)
    eps,di,ddi,ell,dell=sp.symbols("eps D delta_D ell delta_ell",real=True)
    supply=(di+eps*ddi)*(ell+eps*dell)
    variation=sp.diff(supply,eps).subs(eps,0)
    emit("variation of conversion power",variation)
    check("fixed rate and unchanged net energy yield",variation.subs({ddi:0,dell:0}),0)
    check("fixed rate, live net energy yield",variation.subs(ddi,0),di*dell)


if __name__=="__main__":
    print("SymPy version =",sp.__version__,flush=True)
    order_witness()
    pair_requirement()
