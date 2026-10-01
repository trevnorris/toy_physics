#!/usr/bin/env python3
"""Does a live normal drain open a reflection-ODD (twist-parity) bulk channel that H's
parity block-diagonality cannot forbid?

Small symbolic model, built only from (i) the 4D inviscid Euler + continuity equations
(the supplied rest-frame bulk is their U=0 limit), (ii) the integral momentum balance of
a face crossed by mass flux m (Rankine-Hugoniot), with the supplied face traction NORMAL
only (t_s = -(dp + Lambda_X A) n_s, S11c-a :353-354), and (iii) an open-control-volume
momentum balance for the slab material that carries the trapped twist mode.

Class P: fields independent of x3; mirror x3 -> -x3 flips v3 / u3 only.
Bulk background: uniform density rho_m, uniform NORMAL velocity U (drain leaving the
upper face, U>0 outward).  U = 0 is the supplied rest frame.
"""
import sympy as sp

x1, x2, w, t = sp.symbols("x1 x2 w t", real=True)
rho_m, c, U, k2, omega = sp.symbols("rho_m c_s0 U k2 omega", real=True)
A = sp.symbols("A")
# perturbation fields (class P: no x3 dependence)
dr, v1, v2, v3, vw = [sp.Function(n)(x1, x2, w, t) for n in ("drho", "v1", "v2", "v3", "vw")]
dp = c**2 * dr  # barotropic
Ub = sp.Matrix([0, 0, 0, U])          # background 4-velocity (x1,x2,x3,w)
dv = sp.Matrix([v1, v2, v3, vw])
coords = (x1, x2, None, w)            # d/dx3 == 0 on class P

def d(f, i):
    return 0 if coords[i] is None else sp.diff(f, coords[i])

conv = lambda f: sp.diff(f, t) + U * sp.diff(f, w)
cont = sp.expand(conv(dr) + rho_m * sum(d(dv[i], i) for i in range(4)))
mom = [sp.expand(rho_m * conv(dv[i]) + d(dp, i)) for i in range(4)]
print("BULK_CONTINUITY", cont)
for lab, eq in zip(("x1", "x2", "x3", "w"), mom):
    print("BULK_MOMENTUM_" + lab, eq)

# mirror parity of each linearised bulk equation's dependence
mirror = {v3: -v3}
odd_in_even = [sp.simplify(eq.subs(mirror) - eq) for eq in (cont, mom[0], mom[1], mom[3])]
even_in_odd = sp.simplify(mom[2].subs(mirror) + mom[2])
print("EVEN_EQUATIONS_CHANGE_UNDER_MIRROR", odd_in_even)
print("ODD_EQUATION_PLUS_ITS_MIRROR", even_in_odd)

# time-harmonic odd bulk channel: v3 = a(w) exp(i(k2 x2 - omega t))
aw = sp.Function("a")(w)
ansatz = aw * sp.exp(sp.I * (k2 * x2 - omega * t))
ode = sp.simplify(mom[2].subs(v3, ansatz).doit() / sp.exp(sp.I * (k2 * x2 - omega * t)))
print("ODD_CHANNEL_ODE", ode)
sol_U = sp.dsolve(sp.Eq(ode, 0), aw)
print("ODD_CHANNEL_SOLUTION_U", sol_U)
ode0 = sp.simplify(ode.subs(U, 0))
print("ODD_CHANNEL_ODE_AT_U0", ode0, " -> nonzero omega forces a(w)=0:",
      sp.solve(sp.Eq(ode0, 0), aw))
# energy flux of the odd channel through a w=const plane (kinetic energy advected):
amp = sp.symbols("amp", positive=True)
flux = sp.simplify(sp.Rational(1, 2) * rho_m * U * amp**2 / 2)   # time-average of (1/2) rho v3^2 U
print("ODD_CHANNEL_TIME_AVERAGED_ENERGY_FLUX", flux, " at U=0:", flux.subs(U, 0))

# Face jump (control volume of zero thickness straddling the upper face, flat face, normal w):
#  tangential (x3) momentum:  m*(v3_bulk_face - v3_slab) = (tangential stress jump) = 0,
#  because the supplied face traction is purely normal and the bulk is inviscid.
m, v3_bulk_face, u3_t = sp.symbols("m v3_bulk_face u3_t")
tangential_stress_jump = 0
jump = sp.Eq(m * (v3_bulk_face - u3_t), tangential_stress_jump)
print("FACE_TANGENTIAL_JUMP", jump)
print("FACE_TRACE_SOLUTION_m_nonzero", sp.solve(jump.subs(m, sp.Symbol("m0", nonzero=True)), v3_bulk_face))
print("FACE_TRACE_SOLUTION_m_zero", sp.solve(jump.subs(m, 0), v3_bulk_face), "(unconstrained / no coupling)")
D = sp.Matrix([[sp.diff(m * (v3_bulk_face - u3_t), u3_t)]])
print("ODD_TO_ODD_FACE_BLOCK d(jump)/d(u3_t)", D[0], " mirror parity of the block entry:",
      sp.simplify(D[0].subs({u3_t: -u3_t}) - D[0]) == 0)

# Open control volume (the slab material in the trapped-mode region): mass in at rate m with
# zero twist velocity (replenishment from outside the mode), mass out through the face at rate m
# carrying the slab twist velocity V (the jump above).  Restoring force -K X.
M, K, X = sp.symbols("M K X", positive=True)
m0 = sp.symbols("m0", nonnegative=True)
s = sp.symbols("s")
# d(M V)/dt = (m*0) - m*V - K X, dM/dt = 0
char = sp.expand(M * s**2 + m0 * s + K)
roots = sp.solve(char, s)
print("CHARACTERISTIC", char)
print("ROOTS", roots)
print("DECAY_RATE (minus real part)", sp.simplify(-sp.re(roots[0].subs({M: 1, K: 4}))), "at M=1,K=4; general:",
      sp.simplify(-(roots[0] + roots[1]) / 2))
Vt = sp.symbols("V")
E_rate = sp.expand(Vt * (-m0 * Vt - K * X) + K * X * Vt)   # d/dt (1/2 M V^2 + 1/2 K X^2)
print("MODE_ENERGY_RATE", E_rate, "  bulk receives (advected KE flux)", sp.Rational(1, 2) * m0 * Vt**2)
# FORM control: same model with the drain switched off
print("CONTROL_m0_ZERO_ROOTS", sp.solve(char.subs(m0, 0), s))
