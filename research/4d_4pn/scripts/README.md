# 4d_4pn SymPy Archive

This directory contains the referee-facing SymPy verification archive for the
`4d_4pn` paper.

## Layout

- `4pn_*_audit.py`: one file per theorem step in the local and hereditary 4PN chain

## Intended Use

Run the paper-supportive suite script by script, or replay the full stage list in order:

```bash
python3 research/4d_4pn/scripts/4pn_onebody_audit.py
python3 research/4d_4pn/scripts/4pn_quartic_legendre_audit.py
python3 research/4d_4pn/scripts/4pn_local_scaffold_audit.py
python3 research/4d_4pn/scripts/4pn_local_target_import_audit.py
python3 research/4d_4pn/scripts/4pn_local_hamiltonian_to_ordinary_audit.py
python3 research/4d_4pn/scripts/4pn_hamiltonian_chart_generic_frame_lift_audit.py
python3 research/4d_4pn/scripts/4pn_hamiltonian_chart_canonical_slice_audit.py
python3 research/4d_4pn/scripts/4pn_generic_frame_ordinary_translation_audit.py
python3 research/4d_4pn/scripts/4pn_generic_frame_aligned_seed_lift_audit.py
python3 research/4d_4pn/scripts/4pn_local_referee_master_sympy_audit.py
python3 research/4d_4pn/scripts/4pn_tail_bridge_audit.py
python3 research/4d_4pn/scripts/4pn_tail_hereditary_bridge_audit.py
python3 research/4d_4pn/scripts/4pn_conditional_referee_master_audit.py
```

## Scope

This archive verifies the symbolic theorem chain used by the paper within the declared
closure hierarchy. It does not claim a stronger result than the manuscript states.
