#!/usr/bin/env python3
"""Mechanical lookups on the S11c-b SymPy engine: index names, Kronecker-only
basis, zeta_c row, constraint fold. No CAS of the full operator."""
from pathlib import Path
import re

py = Path("/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py").read_text()
wl = Path("/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl").read_text()

print("=== PY index convention ===")
for pat in [
    r"^DIRECTIONS = .*$",
    r"u = tuple\(",
    r'inherited_symbol\(f"u_\{a\}"',
    r"for a in range\(1, 4\)",
    r"str\(i \+ 1\)",
]:
    m = re.search(pat, py, re.M)
    print(f"PAT {pat!r} -> {m.group(0) if m else None}")

print("PY_u_symbol_names_are_1_based_u_1_u_2_u_3 =", 'f"u_{a}"' in py and "range(1, 4)" in py)
print("PY_DIRECTIONS_loop_is_0_based =", "DIRECTIONS = range(3)" in py)

print("\n=== WL index convention ===")
print("WL_spatialCoordinates =", re.search(r"spatialCoordinates = .*", wl).group(0))
print("WL_directions =", re.search(r"directions = .*", wl).group(0))
print("WL_braneDimension =", re.search(r"braneDimension = .*", wl).group(0))

print("\n=== Kronecker vs Levi-Civita in PY basis constructors ===")
print("delta_contractions_defined =", "def delta_contractions(" in py)
print("enumerate_uniform_uses_delta =", "delta_contractions((left, right))" in py)
print("enumerate_new_uses_delta =", "delta_contractions((spurion, left, right))" in py)
print("LeviCivita_in_py =", "LeviCivita" in py)
print("eps_ijk_in_py =", "eps_ijk" in py or "epsilon_ijk" in py)
print("LeviCivita_in_wl =", "LeviCivita" in wl or "Signature[" in wl)

print("\n=== zeta_c / center row / constraint ===")
print("CENTER_FACE_GENERALIZED_ROW =", "CENTER_FACE_GENERALIZED_ROW" in py)
print("constraint_fold_from_source =", "def constraint_fold_from_source(" in py)
print("pin B comment in py =", "constraint" in py.lower())

print("\n=== S11c-b memory history pointers ===")
rec = Path("/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md").read_text()
for phrase in ["16 GB", "15.6 GB", "2 GiB", "30 GB", "64 GB"]:
    print(f"record_has_{phrase!r} =", phrase in rec)
