# 4d_4pn Mathematica Archive

This directory contains the referee-facing Mathematica verification layer for the
`4d_4pn` paper.

Design rules:

- every paper-supportive `.wl` script should run with `math -script`
- every paper-supportive `.wl` script should end with a trailing `Output:` block
- the embedded transcript should match live stdout once the environment-specific OpenMP
  shared-memory warning is ignored

Current manifest roles:

- `4pn_onebody_mathematica_audit.wl`: exact one-body 4PN gate
- `4pn_quartic_legendre_mathematica_audit.wl`: quartic Legendre compiler
- `4pn_local_scaffold_target_mathematica_audit.wl`: local scaffold and fixed-target import
- `4pn_hamiltonian_chain_mathematica_audit.wl`: Hamiltonian local chain through the canonical slice
- `4pn_ordinary_chain_mathematica_audit.wl`: ordinary-chart translation and aligned-seed lift chain
- `4pn_tail_bridge_mathematica_audit.wl`: GR tail shadow and hereditary bridge
- `4pn_conditional_referee_master_mathematica_audit.wl`: conditional full-4PN closure replay

The current manifest is kept in `wl_notes.txt`, and `update_output_blocks.py` refreshes the
embedded transcripts sequentially.
