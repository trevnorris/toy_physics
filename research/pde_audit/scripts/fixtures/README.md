# Fixture Bundle

This ZIP contains the scripted fixture/sample JSON artifacts created during the V2-21 through V2-22C stages.

## Included files

- `stage_v2_21_sample_branch_manifest.json` — sample branch manifest fixture for the extraction fixture.
- `stage_v2_22a_profile_input_manifest.json` — sample profile-input manifest for the profile-to-coefficient adapter.
- `stage_v2_22b_solver_output_schema.json` — schema for expected solver output.
- `stage_v2_22b_sample_solver_output_valid.json` — valid sample solver-output fixture.
- `stage_v2_22b_sample_solver_output_invalid_hardcap.json` — intentionally invalid hard-cap negative-control fixture.
- `stage_v2_22c_valid_solver_packet.json` — valid open-throat solver packet fixture used in the end-to-end smoke pipeline.
- `stage_v2_22c_invalid_hardcap_solver_packet.json` — intentionally invalid hard-cap solver packet.

These are script-generated fixtures / controls, not claimed physical branch-discovery outputs.
