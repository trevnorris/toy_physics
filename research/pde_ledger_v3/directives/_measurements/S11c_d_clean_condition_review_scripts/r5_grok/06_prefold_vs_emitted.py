#!/usr/bin/env python3
"""Mechanical lookup: the four named operator paths at the cited creation
site versus the emitted overwrite site. No CAS derivation.
"""
from pathlib import Path

engine = Path(
    "/var/projects/toy_physics/research/pde_ledger_v3/scripts/"
    "S11c_b_brane_operator_sympy_audit.py"
).read_text()
lines = engine.splitlines()

ranges = {
    "cited_creation_2416_2433": (2416, 2433),
    "fold_overwrite_U_2967_2999": (2967, 2999),
    "emitted_THETA_3095_3101": (3095, 3101),
    "emit_site_4156_4183": (4156, 4183),
    "s11ca_76_100": None,
    "s11ca_101_105": None,
}

print("--- engine cited creation site ---")
for i in range(2415, 2434):
    print(f"{i}:{lines[i-1]}")

print("--- emitted THETA_BALANCE overwrite ---")
for i in range(3094, 3103):
    print(f"{i}:{lines[i-1]}")

print("--- operator_from_density THETA assignment ---")
for i, line in enumerate(lines, start=1):
    if 'operator = {' in line and 2400 < i < 2440:
        print(f"{i}:{line}")
    if i >= 2420 and i <= 2433 and "THETA" in line:
        print(f"{i}:{line}")

s11ca = Path(
    "/var/projects/toy_physics/research/pde_ledger_v3/scripts/"
    "S11c_a_interface_geometry_sympy_audit.py"
).read_text().splitlines()
print("--- S11c-a 76-105 ---")
for i in range(76, 106):
    print(f"{i}:{s11ca[i-1]}")

print("CITED_CREATION_HAS_MU_THETA_AS_THETA_EXPANDED", "mu_theta_amplitude" in "\n".join(lines[2415:2433]))
print("EMITTED_THETA_USES_MASS_BALANCE", "mass_balance" in "\n".join(lines[3094:3102]))
print("S11CA_76_100_HAS_X1", any("x1 =" in s11ca[i] for i in range(75, 100)))
print("S11CA_101_HAS_X1", "x1 =" in s11ca[100])
