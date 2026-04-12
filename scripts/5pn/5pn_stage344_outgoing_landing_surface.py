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


banner('STAGE 344 — EXACT OUTGOING LANDING SURFACE')

chiQ = sp.symbols('chi_Q', positive=True, real=True)
DeltaQ = sp.symbols('Delta_Q', real=True)
sigma_star = sp.symbols('sigma_star', positive=True, real=True)
Xi_slip, dPi_tan = sp.symbols('Xi_slip deltaPi_tan', real=True)

NQ_from_chi = sp.simplify(1 / chiQ)
print('Natural source-map relation:')
print('  N_Q = 1/chi_Q =')
sp.pprint(NQ_from_chi)

expect_zero('N_Q - 1 at chi_Q = 1', sp.simplify(NQ_from_chi.subs(chiQ, 1) - 1))

print('\nOutgoing-defect parameterization:')
print('  Delta_Q = chi_Q - 1')
NQ_from_Delta = sp.simplify(1 / (1 + DeltaQ))
print('  N_Q =')
sp.pprint(NQ_from_Delta)
print('  N_Q - 1 =')
sp.pprint(sp.simplify(NQ_from_Delta - 1))

# first-order Family-1 fixed-point bridge
DeltaQ_from_slip = sp.simplify(-sigma_star/(1 - sigma_star) * Xi_slip * dPi_tan)
print('\nFirst-order outgoing bridge on the compensated Family-1 branch:')
print('  Delta_Q =')
sp.pprint(DeltaQ_from_slip)
expect_zero('Delta_Q on-family if Xi_slip = 0', DeltaQ_from_slip.subs(Xi_slip, 0))
expect_zero('Delta_Q on-family if deltaPi_tan = 0', DeltaQ_from_slip.subs(dPi_tan, 0))

NQ_minus_1_from_slip = sp.simplify(1/(1 + DeltaQ_from_slip) - 1)
print('  N_Q - 1 =')
sp.pprint(NQ_minus_1_from_slip)
expect_zero('N_Q - 1 on-family if Xi_slip = 0', NQ_minus_1_from_slip.subs(Xi_slip, 0))

print('\nInterpretation:')
print('  On the natural source-map branch, the last reduced outgoing condition is')
print('      N_Q = 1  <=>  chi_Q = 1.')
print('  On the exact lower parent compensation family, the first-order outgoing defect')
print('  vanishes automatically whenever Xi_slip = 0, and then N_Q - 1 = 0 follows.')
