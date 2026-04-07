# 4D 2.5PN SymPy Archive

The paper-facing SymPy verification layer has two levels.

- `2_5pn_master_sympy_audit.py` covers the core theorem chain used in the main text:
  - Burke--Thorne benchmark
  - low-frequency odd-channel bookkeeping
  - scalar-sector audit and demotion
  - dipole/vector-sector audit and demotion
  - surviving quadrupole representation/source-map checks

- `2_5pn_master_session_sympy_audit.py` is the session-completion superset:
  - everything in the master audit
  - plus the quadrupole-normalization package
  - canonical/source-port normalization cleanup
  - convention-invariant products
  - grouped real `P_2` formulas
  - the minimal isotropic quadrupole module
  - the single-pole sufficiency/extraction ledger

Both files are standalone and do not import local project modules.
