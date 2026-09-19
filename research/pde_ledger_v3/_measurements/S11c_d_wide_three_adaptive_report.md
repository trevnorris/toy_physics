# S11c-d independent wider-box triple outer quadrature

Production is accepted and annex-verified. All four workers and the supervisor
exited zero with empty stderr after **22h45m**. Independent adaptive GK21
replaced only the outer Gauss216 rule at cutoffs3/4 on both approved fields.
Inner32/32, source/profile256/512, position bounds48/14 and regulator0.2 stayed
fixed; all40 single/all30 paired rows were held in every80-row action.

Each adaptive solve completed after315 points,1260 in total. The largest raw
integral difference from Gauss was2.18e-13 and the largest complete-action
difference1.39e-17. Reported adaptive unit-frame error estimates reach3.24e-11
against requested1e-10; these are numerical estimates, not rigorous bounds.
The sampled cutoff3-to4 action change remains1.13467310948e-7.

Saved-operand acceptance verifies all92 current/frozen sources,46266 worker
artifacts including44990 partials, eight exact caches and1680 zero held-term
scalars. All56185 metadata paths,9258 tags,4627 write keys, native row/source/
profile/field/ordered-limit/current-setting joins and pre/post packet hashes
pass. Actual measure mutations, adaptive interval partitions and saved sums
pass. Checks/stdout are identical; acceptance repeats no numerical integration.
Peak worker/coordinator RSS was133008/252144 KiB.

The2.81 MB transcript is `scripts/out/S11c_d_wide_three_adaptive.out`, published
at `4ba95bcf`. Its MD5E annex key, symlink and independent SHA256 are verified:
`5045909379429de3ee29153db4af1c4cdc05544d679df6f1e4b000d0040e3ce3`.
The full accepted record is `S11c_d_wide_three_adaptive_checkpoint.json`.
Preflight acceptance remains `bdb86303`; its coarse smoke is instrument evidence.
The native engine, approved inputs, prior operands and retained contract suffix
remain unchanged. New-stage provenance includes v10 plus nonlinearPoleV2;
this numerical-action calculation invokes no pole constructor.

The user-approved `exploratoryAcceptanceV1` at `30440cac` governs future work.
This result closes the planned fixed-box quadrature sequence. Next follow
`S11c_d_focused_completion_plan.md` to build the first complete finite scattering
solve, then check its important observables at practical stated precision.
Finite Gaussian action checks alone establish no scattering error bound,
infinite tail, Abel limit or pole result. Rigorous full-operator certification
is optional follow-up; no additional broad refinement sweep is queued.
