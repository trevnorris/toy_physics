# g0_closure_card_v0_checks.py — verification for g0_closure_card_v0.md
# Recovered from the Codex card-build (session 019f80ac-765b-7213-b455-dce77e79e5a2) heredocs;
# re-run reproduces the card's admissibility/instantiability checks. This is the HEADLINE block
# (stability margin, generalized wave speeds, transverse factorization + zero-mode norm, dimensions).
# TODO before full trust: (a) tautology/faithfulness review vs the card's ACTUAL operators;
#   (b) Mathematica dual-engine (.wl); (c) fold in remaining build-log blocks (wall factorization, bulk Hessian).

import sympy as sp

# Numeric/symbolic scalar block.
Aeff = sp.Rational(1)
Mh = sp.Rational(1)
Beff = sp.Rational(2)
Kh = sp.Rational(1)
Chu = sp.Rational(1, 2)
M = sp.diag(Aeff, Mh)
K = sp.Matrix([[Beff, Chu], [Chu, Kh]])
z = sp.symbols('z', real=True)
char = sp.factor((K-z*M).det())
roots = sp.solve(char, z)
margin = sp.factor(Beff*Kh-Chu**2)

# Parent transverse factorization and zero-mode equation.
w, ell = sp.symbols('w ell', positive=True, real=True)
f0 = 1/(ell*sp.cosh(w/ell)**2)
V = (4-6/sp.cosh(w/ell)**2)/ell**2
Of = sp.simplify(-sp.diff(f0,w,2)+V*f0)
q = 2*sp.tanh(w/ell)/ell
factor_potential = sp.simplify(q**2-sp.diff(q,w)-V)
N0 = sp.integrate(2*f0**2, (w,-sp.oo,sp.oo))

# (L,T,M) dimensional algebra.
D = {
    'E': sp.Matrix([2,-2,1]), 'L': sp.Matrix([1,0,0]),
    'T': sp.Matrix([0,1,0]), 'M': sp.Matrix([0,0,1]),
}
E,L,T,Mass = D['E'],D['L'],D['T'],D['M']
r = sp.zeros(3,1); H=-L; h=r; u=L
Zchi=Mass-2*L; kchi=Mass-2*T; lam=Mass-2*L-2*T
M4=Mass; K4=E
Aeff_d=Mass-3*L; Mh_d=Mass-L; Beff_d=Mass-L-2*T
Kh_d=Mass+L-2*T; Chu_d=Mass-2*T
Km=E+2*L; Jm=E+L; eta=-3*L
rho=-4*L; S=-4*L-T; vel=L-T; m=Mass
energy_density=E-4*L; momentum=m+rho+vel
checks = {
 'wall_kin': Zchi-2*T,
 'wall_grad': kchi-2*L,
 'wall_pot': lam,
 'H_kin': M4+2*H-2*T,
 'H_grad': K4+2*(H-L),
 'H_pot': K4+2*H-2*L,
 'u_kin': Aeff_d+2*(u-T),
 'u_stiff': Beff_d+2*(u-L),
 'h_kin': Mh_d-2*T,
 'h_stiff': Kh_d-2*L,
 'mix': Chu_d+(u-L)+(h-L),
 'mouth_robin_integrand': eta+Km+2*H,
 'mouth_source_integrand': eta+Jm+H,
 'mass_source': S,
 'momentum_source': m+S+vel,
 'energy_source': E+S,
}
targets = {
 **{k: E-4*L for k in ('wall_kin','wall_grad','wall_pot','H_kin','H_grad','H_pot')},
 **{k: E-3*L for k in ('u_kin','u_stiff','h_kin','h_stiff','mix')},
 'mouth_robin_integrand': E-3*L,
 'mouth_source_integrand': E-3*L,
 'mass_source': -4*L-T,
 'momentum_source': momentum-T,
 'energy_source': E-4*L-T,
}
residuals = {k: tuple(checks[k]-targets[k]) for k in checks}
assert Mh > 0 and margin > 0
assert all(rt > 0 for rt in roots)
assert Of == 0 and factor_potential == 0
assert sp.simplify(N0-8/(3*ell)) == 0
assert all(v == (0,0,0) for v in residuals.values())
print('PASS_STABILITY margin =', margin)
print('PASS_GENERALIZED_SPEEDS det(K-zM) =', char)
print('c_squared =', roots, 'numeric =', [sp.N(x,12) for x in roots])
print('PASS_TRANSVERSE Of0 =', Of, 'factor_residual =', factor_potential, 'N0 =', N0)
print('PASS_DIMENSIONS residuals =', residuals)
