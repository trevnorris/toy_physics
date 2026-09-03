# Independent physics review — #90 PY §3c coupling content build

You are one of two independent build-review legs. The artifact is a Codex-written change to a SymPy engine (S11c-b
coupling kernel: extract the face + response coupling). Derive independently; ablate a COPY; report literal stdout.
A prose re-derivation is discarded: write your own script and save both script and literal stdout to named absolute
paths, and report them. Report a finding only if it catches a way the physics could be WRONG. Report observed
payloads as findings — do NOT use payload membership or a residual value as a pass/fail exit condition.

## Artifact
The working-tree change (vs `git HEAD`) to:
- `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (the SymPy engine)
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (§0 clarity pin)

Read the diff: `git -C /var/projects/toy_physics diff HEAD -- research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`

## What the build was supposed to do (settled physics — implement, not re-adjudicate)
§3c mandates INCLUDE/INCLUDE: the coupling kernel carries the reversible tilted-face geometry (from the S11c-a
T-a..T-i shape-derivative substrate) AND the irreversible face response (§1c `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)`, Λ symbolic).
PY was bulk-only. The correct fix (NOT a raw-bundle extraction — that is the §3c-forbidden parallel route): compute
the face GENERALIZED-FORCE rows from the consumed virtual work (coefficients of the independent virtual displacements
`δ_vu`/`δ_ve_W`; live `μ_θ` bound via `MU_THETA_FACE_BINDING`, not the reserved `mu_theta_L/M` placeholder), ADD them
to the constraint-reduced operator rows (origin FACE, NOT through the θ-constraint fold), then apply the SAME weak
restriction to the full operator. Directive: `directives/S11c_b_90_coupling_content_directive.md`. #84 records:
`directives/_measurements/S11c_b_coupling_84_{diagnosis,consult}.md`. You are NOT asked whether to include
face/response — that is settled.

## What you are handed
The engine; the S11c-a substrate engine `scripts/S11c_a_interface_geometry_sympy_audit.py` (source of the T-a..T-i
substrate, virtual work, closure/traction, and the reserved `mu_theta_L/M`); the spec `S11c_b_SHARED_PHYSICS.md`
§0/§1c/§3b/§3c and `S11c_a_SHARED_PHYSICS.md` T-i (`:448`). Copy the engine to `/tmp` and ablate the COPY; never
modify the working tree.

## Required probes — each an ABLATION with literal stdout (a FORM ablation is MANDATORY)
1. **Provenance: the face content is COMPUTED INTO the operator, not extracted from the raw bundle.** Confirm the
   emitted `S11CB_COUPLING_KERNEL` face terms trace to face generalized-force rows added to the operator rows, and
   that the weak restriction is applied to those operator rows (not to `FACE_FLUX_BOUNDARY_OPERANDS` directly). If
   the raw bundle is weak-restricted, that is the §3c parallel-route defect.
2. **FORM-ablate the face SOURCE.** On the COPY, corrupt/drop a virtual-work / traction / `Λ` source feeding the
   face rows, re-run, and REPORT the literal diff of the kernel. Separately ablate a NON-face substrate key (e.g. a
   projection object) and report its diff. (Interpretation is the orchestrator's; you report the diffs.)
3. **Delete the face-row feed** on the COPY and report whether the face-origin kernel terms vanish while the bulk
   `γ·profile-jet` terms remain.
4. **Λ symbolic.** Confirm no bulk-response solve, DtN, impedance, or `Z` appears — `Λ_I(ω)` stays a symbolic factor.
5. **Exact-once.** Confirm `evolution_mass_balance` (already the θ-row) is not re-contracted/double-counted, and the
   already-weak virtual-work density is not multiplied by a second test field.
6. **`μ_θ` bound.** Confirm the reserved `mu_theta_L`/`mu_theta_M` placeholder does NOT survive in the coupling
   blocks (the live `mu_theta_amplitude`/`MU_THETA_FACE_BINDING` is bound into the face affinity).
7. **No `ζ_c` over-reach + no channel filter + adjointness over enlarged blocks.** Confirm no center-face channel is
   added; the extraction drops no content by construction (not a presence quota); the adjointness object is formed
   from the complete (bulk+face) blocks.
8. **Controls + #88 + no leak.** Confirm `operator_from_density`/`committed_strong_rows` are unchanged (raw #88
   refs); the §5 controls re-run and PRINT residuals (no assert); the §0 amendment matches §1c/T-i; every emitted
   payload is a CAS object; report any `assert` preceding a guarded value.

## Operational
Kernel-free SymPy ablations run directly (no Mathematica seat). Smoke on a SINGLE `(branch, representative)` case
(e.g. `task_coupling_kernel` on one case, or `S11CB_PRIMARIES_ONLY=1`); ⛔ do NOT run the full package loop (OOMs
~3.5h). Save every ablation script + literal stdout to named absolute paths and report them.

## Output
Per probe: CLEAR / FINDING (the specific physics-wrong path it catches), backed by your script path + literal stdout
+ the computing line. If all CLEAR, say what each ablation would have caught had it been wrong. Do not edit the tree.
