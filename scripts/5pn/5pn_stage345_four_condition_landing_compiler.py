#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner('STAGE 345 — FOUR-CONDITION LANDING COMPILER')

chi0, deltaU = sp.symbols('chi_0 delta_U', positive=True, real=True)
epsEta, eps = sp.symbols('epsilon_eta epsilon', positive=True, real=True)
ZW, Lambda = sp.symbols('Z_W Lambda', positive=True, real=True)
chiQ = sp.symbols('chi_Q', positive=True, real=True)
Xi_slip, dPi_tan = sp.symbols('Xi_slip deltaPi_tan', real=True)
sigma_star = sp.symbols('sigma_star', positive=True, real=True)

# drift coordinates
u_chi, u_delta, u_ZW = sp.symbols('dlnchi_0 dlndelta_U dlnZ_W', real=True)
u_epsEta, u_eps, u_Lambda = sp.symbols('dlnepsilon_eta dlnepsilon dlnLambda', real=True)

Rtr = (1 + chi0/(1 + deltaU)) / (1 + chi0)
Rtarget = Lambda * (1 - epsEta) * (1 - eps)**2 / (ZW * (1 + chi0)**2)

# first three reduced finish-line conditionals
C1 = sp.simplify(((chi0*u_chi + deltaU*u_delta)/(1 + chi0 + deltaU)
                  - (chi0*u_chi)/(1 + chi0)
                  - (deltaU*u_delta)/(1 + deltaU)))
C2 = sp.simplify(u_Lambda - u_ZW - (epsEta/(1 - epsEta))*u_epsEta - 2*chi0/(1 + chi0)*u_chi - 2*eps/(1 - eps)*u_eps)
C3 = u_epsEta

u_delta_req = sp.simplify(-(1 + deltaU)/(1 + chi0) * u_chi)
u_Lambda_req = sp.simplify(u_ZW + 2*chi0/(1 + chi0)*u_chi + 2*eps/(1 - eps)*u_eps)

# fourth condition on natural source-map branch
DeltaQ = sp.simplify(-sigma_star/(1 - sigma_star) * Xi_slip * dPi_tan)
NQ_minus_1 = sp.simplify(1/(1 + DeltaQ) - 1)

print('Four reduced finish-line conditionals:')
print('  C1 = d ln R_tr')
sp.pprint(C1)
print('  C2 = d ln R_target')
sp.pprint(C2)
print('  C3 = d ln epsilon_eta')
sp.pprint(C3)
print('  C4 = N_Q - 1')
sp.pprint(NQ_minus_1)

# landing substitutions
landing_subs = {
    u_delta: u_delta_req,
    u_epsEta: 0,
    u_Lambda: u_Lambda_req,
    Xi_slip: 0,
}

expect_zero('C1 on orbit-lock surface', C1.subs(landing_subs))
expect_zero('C2 on orbit-lock surface', C2.subs(landing_subs))
expect_zero('C3 on orbit-lock surface', C3.subs(landing_subs))
expect_zero('C4 on exact lower parent compensation family', NQ_minus_1.subs(landing_subs))

# canonical outgoing alternative
expect_zero('C4 on canonical outgoing branch chi_Q = 1', (1/chiQ - 1).subs(chiQ, 1))

print('\nCombined landing theorem (first-order / within closure):')
print('  If the branch lies on the exact orbit-lock surface')
print('      dln delta_U =')
sp.pprint(u_delta_req)
print('      dln Lambda =')
sp.pprint(u_Lambda_req)
print('      dln epsilon_eta = 0')
print('  and the outgoing branch is either')
print('      (a) canonical with chi_Q = 1, or')
print('      (b) tangent to the exact lower parent compensation family with Xi_slip = 0,')
print('  then all four reduced finish-line conditionals vanish simultaneously.')

print('\nInterpretation:')
print('  The four-condition finish line is algebraically reachable inside the current')
print('  reduced hierarchy. The remaining open question is branch realization by the')
print('  completed moving-throat PDE, not another reduced-sector contradiction.')
