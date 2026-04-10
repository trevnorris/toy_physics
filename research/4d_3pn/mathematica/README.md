# 4d_3pn Mathematica Archive

This directory will contain the referee-facing Mathematica verification layer for the
`4d_3pn` paper.

Current contents:

- `3pn_onebody_mathematica_audit.wl`
- `3pn_grouped_p2_mathematica_audit.wl`
- `3pn_comparable_mass_mathematica_audit.wl`
- `3pn_com_linear_map_mathematica_audit.wl`
- `3pn_com_gr_target_mathematica_audit.wl`
- `3pn_generic_frame_projection_repair_contact_mathematica_audit.wl`
- `3pn_generic_frame_lift_compiler_mathematica_audit.wl`
- `3pn_grouped_p2_closure_mathematica_audit.wl`
- `3pn_scalar_kinetic_final_mathematica_audit.wl`
- `wl_notes.txt`
- `update_output_blocks.py`

The current Mathematica archive now spans the full referee-facing theorem chain:

- one-body gate and grouped kickoff
- comparable-mass COM compiler foundations
- imported GR COM target
- generic-frame projection / seed repair / contact orbit
- Hamiltonian-first lift and fixed-chart compiler
- grouped real `P_2` no-go and richer exact closure
- scalar geometry completion, pure-kinetic collapse, and the final three-lane theorem

The design target matches the archive rules already used in the other `research/*/mathematica`
paper bundles: every referee-facing script should run with `math -script`, end with a
trailing `Output:` block, and support direct comparison between the embedded transcript and
live stdout.
