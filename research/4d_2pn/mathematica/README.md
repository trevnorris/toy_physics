# 4D 2PN Mathematica Archive

The paper-supportive Wolfram verification layer has two levels.

- `4d_gravity_2pn_master_harness.wl` is the compact referee harness cited in the manuscript.
- `4d_2pn_full_verification_harness.wl` is the manifest-driven suite runner for the broader symbolic derivation ledger listed in `wl_notes.txt`.

The suite manifest intentionally excludes utility files that are useful for local exploration but are not part of the paper's referee-facing evidence:

- `4d_2pn_quicklook.wl`
- `black_hole_scaling.wl`
- `lensing_from_flow.wl`
- `nonlinear_profile_solver.wl`

Those utilities still run, but they should not be cited as support for the conservative two-body 2PN paper.
