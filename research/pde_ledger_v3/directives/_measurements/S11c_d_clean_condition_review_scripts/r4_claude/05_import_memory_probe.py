#!/usr/bin/env python3
"""Leg script 05: peak RSS of importing a COPY of S11c_a_exports.py, which the existing
S11c-b SymPy engine imports at module load (S11c_b_brane_operator_sympy_audit.py:37,
`from S11c_a_exports import LEDGER as INCOMING_LEDGER`).  Directive B's whitelist lets the
SymPy builder import symbol definitions from 'the existing SymPy engine'; the guard caps the
whole job at 2 GiB.  Prints measured numbers only."""
import resource, sys, time
sys.path.insert(0, "/tmp/s11cd_clean_review_r4_claude/import_probe")
t0 = time.time()
import S11c_a_exports
print("IMPORT_SECONDS", round(time.time() - t0, 1))
print("LEDGER_ENTRIES", len(S11c_a_exports.LEDGER))
print("PEAK_RSS_MIB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1))
print("GUARD_CAP_MIB", 2 * 1024)
