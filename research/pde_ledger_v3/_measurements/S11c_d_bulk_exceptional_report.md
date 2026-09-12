# S11c-d independent bulk exceptional-frequency construction

Completed the accepted finite-slice construction after checkpoint `3332fe0a`.
That commit saved nineteen selected files in Git and the reference `.out`
through DataLad/git-annex. The subsequent build and regenerated outputs are
uncommitted. S10/Lean work was not changed, and nothing was pushed.

[Bulk geometry](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3484)
now derives branch and entry-denominator conditions directly from the reduced
physical matrix and radical relation. It retains square-free multiplicities,
leading coefficients, discriminants, intersections, and exact real-root
isolation. [Real-normal projection](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3455)
separates shared real/imaginary curves from residual intersections, so a
real-axis denominator crossing is enumerated without requiring coalescence.
The end-mode determinant is a separate intersection operand, never a gate.
The dimensionless [synthetic control](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_projection_checks.json)
computes frequency one and two real normal lifts with zero equation and
projection-reconstruction residuals; it makes no physical claim.

[Target construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3547)
evaluates every computed branch/denominator lift at each nonnegative critical
frequency and retains singular targets individually. Exact rational witnesses
sample each enumerated positive-frequency interval and each real-normal cell;
both algebraic radical lifts are evaluated. [Region and bank construction](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3704)
uses the inherited Fourier seed and records separate paths, matrix jumps,
refinements and domain failures. [Matrix arithmetic](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_d_mixing_scattering_sympy_audit.py:3669)
uses 40/60 digits at stored points. Continuation coordinates remain 53-bit;
arithmetic refinement does not improve their spectral accuracy.

The [full inventory](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_full_exceptional_inventory_summary.json)
contains twelve bulk and twelve end-mode slices. Across the bulk packets it
records 84 critical-point evaluations: twelve finite branch matrices and 72
singular-denominator records. Thirty-six frequency intervals contain sixty
real-normal cells, with 120 radical lifts, 408 matrix/inverse evaluations,
288 transported bank paths and 144 bank pairs. All 96 interval witnesses lie
inside their bounds; all 48 target root censuses have zero degree residual.
The inherited end slice supplies 312 exact generic minor identities and 24
threshold records. Each threshold has rank three and full left/right bases
of rank/nullity two; its normal-derivative pairing has rank zero. Generalized
normal modes remain unresolved. Twelve zero-frequency intersections retain
separate unresolved records.

In the declared unit frame, the reference bulk branch collision occurs at
frequency coefficient `sqrt(5)`, with `k=q=0` and physical-matrix rank five.
It is not an end-mode root. The separate reference normal threshold is
`sqrt(3/40)`. Targeted radical and full-subspace residuals are zero. The largest
recorded new left inverse residual is `3.652e-59` at restored dimension
`[1,0,0]`; the inverse arithmetic refinement at `[3,2,-1]` is `2.426e-38`.
These are observations, not a global accuracy certificate or acceptance bound.

Four focused runs and the single four-case production run use engine SHA
`bc3257ab3c846896e5af7f6b8221044bcb4113d21e4b7c8050d95ab9a93bcdb4`.
Production took 9,632.055 seconds with 1,730,232 KiB peak child RSS, exit zero,
empty stderr, and all 24 source/input pins unchanged. Its 124,804 tags have no
duplicates; the inventories found no metadata gaps, unmatched objects or
nonfinite exceptional fingerprints. Solved reduced/pencil dimensional
residuals are empty; the raw input-dimensional equations remain emitted.

[Preservation](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_full_payload_preservation.json) accounts for all
109,612 prior tags: 108,352 payloads are text-identical; 936 changes add native
fields or reorder associations; 74 reorder mappings/closure edges with matched
metadata; 246 change one injectively renamed SymPy dummy identity; four update
run provenance/checkpoint data. No prior native value changed, no prior tag is
missing, and no difference is unclassified. All 432 native modes, 432 normal
residues, 1,728 contours and 504 prior sheet paths remain. The 24 unresolved
branch-path intersections and 48 unresolved native sheet labels remain explicit.
All 96 Fourier carriers and their reconstruction residuals are preserved.

The [codec check](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_full_codec_inventory.json) restores every payload
and all 124,802 source-index assignments exactly. The expanded transcript is
217,079,452 bytes; the published main `.out` is 58,095,644 bytes (SHA
`459c047acee94c5f307bd914809747954d4fb2b2d46bc36c78ace8e1b500e6c2`).
Atomic replacement preserved the previous 189,492,142-byte annex payload.
Four focused transcripts are also published under `scripts/out/`.
[Run provenance](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_bulk_exceptional_runs.json) retains commands, hashes,
resources and development failures. Large exact coefficients now use the
arbitrary-precision carrier fingerprint path, avoiding binary64 overflow;
the earlier failed serialization attempts are retained in that history.

These are bound-carrier frequency slices, real-axis samples and selected
complex paths. They do not prove constant matrix properties on entire cells,
a global parameter/sheet atlas, exceptional-point physical sheet assignment,
or section 3b's profile-dependent frequency bound poles. Constant-end
normal-momentum poles remain separate. All ten broad TODOs and the export
remain open. The next construction is generalized threshold modes and their
exceptional-point sheet continuation, before completing current/flux
normalization and variable-profile scattering. No upstream operator discrepancy
was established here. Section 1 premises and the carried c2 operand/sign and
shear-normalization debts remain supplied. No review, comparator, Wolfram or
downstream run was performed; the solver/export contract is byte-identical.
