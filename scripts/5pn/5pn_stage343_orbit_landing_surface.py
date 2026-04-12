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


banner('STAGE 343 — EXACT ORBIT-LOCK LANDING SURFACE')

chi0, deltaU = sp.symbols('chi_0 delta_U', positive=True, real=True)
epsEta, eps = sp.symbols('epsilon_eta epsilon', positive=True, real=True)
ZW, Lambda, zeta = sp.symbols('Z_W Lambda zeta', positive=True, real=True)

# logarithmic branch drifts
u_chi, u_delta, u_ZW = sp.symbols('dlnchi_0 dlndelta_U dlnZ_W', real=True)
u_epsEta, u_eps = sp.symbols('dlnepsilon_eta dlnepsilon', real=True)
u_Lambda = sp.symbols('dlnLambda', real=True)

Rtr = (1 + chi0/(1 + deltaU)) / (1 + chi0)
Rtarget = Lambda * (1 - epsEta) * (1 - eps)**2 / (ZW * (1 + chi0)**2)

print('R_tr =')
sp.pprint(sp.simplify(Rtr))
print('R_target =')
sp.pprint(sp.simplify(Rtarget))
print('epsilon_eta =')
sp.pprint(epsEta)

# exact first-order observable drifts
# use logarithmic differentials directly
# d ln R_tr
term_Rtr = ((chi0*u_chi + deltaU*u_delta)/(1 + chi0 + deltaU)
            - (chi0*u_chi)/(1 + chi0)
            - (deltaU*u_delta)/(1 + deltaU))
term_Rtr = sp.simplify(sp.factor(term_Rtr))

term_Rtarget = sp.simplify(
    u_Lambda - u_ZW - (epsEta/(1 - epsEta))*u_epsEta - 2*chi0/(1 + chi0)*u_chi - 2*eps/(1 - eps)*u_eps
)

print('\nd ln R_tr =')
sp.pprint(term_Rtr)
print('d ln R_target =')
sp.pprint(term_Rtarget)
print('d ln epsilon_eta =')
sp.pprint(u_epsEta)

# support-blindness of orbit packet
expect_zero('d_zeta ln R_tr', sp.diff(sp.log(Rtr), zeta))
expect_zero('d_zeta ln R_target', sp.diff(sp.log(Rtarget), zeta))
expect_zero('d_zeta ln epsilon_eta', sp.diff(sp.log(epsEta), zeta))

# exact orbit-lock surface relations
u_delta_req = sp.simplify(-(1 + deltaU)/(1 + chi0) * u_chi)
term_Rtr_sub = sp.simplify(term_Rtr.subs(u_delta, u_delta_req))
expect_zero('d ln R_tr after tracking law', term_Rtr_sub)

u_Lambda_req = sp.simplify(u_ZW + 2*chi0/(1 + chi0)*u_chi + 2*eps/(1 - eps)*u_eps)
term_Rtarget_sub = sp.simplify(term_Rtarget.subs({u_epsEta: 0, u_Lambda: u_Lambda_req}))
expect_zero('d ln R_target after target law + dressing lock', term_Rtarget_sub)

print('\nExact orbit-lock landing surface:')
print('  1. (1 + delta_U) dln chi_0 + (1 + chi_0) dln delta_U = 0')
print('  2. dln Lambda - dln Z_W - 2 chi_0 dln chi_0/(1+chi_0) - 2 epsilon dln epsilon/(1-epsilon) = 0')
print('  3. dln epsilon_eta = 0')
print('\nEquivalent solved form:')
print('  dln delta_U =')
sp.pprint(u_delta_req)
print('  dln Lambda =')
sp.pprint(u_Lambda_req)
print('  dln epsilon_eta = 0')

print('\nInterpretation:')
print('  These three equations are the exact first-order landing surface for')
print('  delta ln R_tr = delta ln R_target = delta ln epsilon_eta = 0,')
print('  and they are support-blind with respect to zeta.')
