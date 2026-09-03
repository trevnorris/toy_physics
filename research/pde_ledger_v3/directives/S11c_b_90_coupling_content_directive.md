# Builder directive — #90: PY §3c coupling content (compute the face rows into the operator, then restrict) + §0 pin

## Status (read first)

Orchestrator-written builder decision list. **v2** — v1 was REJECTED by both decision legs (Codex + Grok,
independently, computation-backed): v1 told the builder to weak-restrict the raw `FACE_FLUX_BOUNDARY_OPERANDS`
bundle (the §3c-forbidden parallel route), overloaded the token `A_T`, left `μ_θ` unbound, risked double-counting
the θ face flux and re-pairing the already-weak virtual work, over-reached into non-coupling substrate, overstated
the §0 pin, and leaked expected kernel content. v2 folds both legs' convergent replacement (leg logs
`~/.s11_build/S11c_b_90_coupling/decision_{codex,grok}.log`, `.log` gitignored; key claims orchestrator-verified,
rule 13). It goes to Codex build + 2 build legs; NOT re-legged (rule 7: one pass, fold, go).

**SETTLED (do not re-litigate — adversarially confirmed Codex+Grok ×2, user-endorsed):** §3c mandates
**INCLUDE / INCLUDE** — the coupling kernel carries the reversible tilted-face geometry (from the S11c-a T-a..T-i
shape-derivative substrate) AND the irreversible face response (the §1c `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)` flat-face closure).
**WL is spec-correct; PY UNDER-EXTRACTS — its coupling is bulk-stored-energy only.** Records
`directives/_measurements/S11c_b_coupling_84_{diagnosis,consult,basis_verification}.md`. Controlling argument
(verified verbatim): §3c's "weak variational restriction under the stored-energy/kinetic pairing" is the EXTRACTION
INSTRUMENT, ⛔ not a content filter; §3c requires a block "of the §3b operator itself", ⛔ "not from a parallel
direct-variation route", ⛔ "do not filter the kernel to a single channel"; §3b consumes T-a..T-i "for every
boundary/face contribution TO THE OPERATOR"; TERM_ORIGINS classify bulk-energy/face-flux/advective. T-i
(`S11c_a_SHARED_PHYSICS.md:448-452`) fixes that the supplied `Λ` is NOT the §0/S11c-c bulk solve (`δp=Z·v_bulk`)
⇒ keep `Λ(ω)` SYMBOLIC (no bulk elimination / DtN).

## The root (measured on the folded engine)

`task_coupling_kernel` (`:3874`) → `build_kernel` (`:3254`) → `weak_operator_blocks` (`:3404`) contracts the
trial/test fields ONLY against the operator's EXPANDED bulk-EL rows (`U/THETA/E_W_BALANCE.EXPANDED`; the θ-row is
now the imported `evolution_mass_balance`). The T-substrate is attached AFTER the blocks as a raw side payload
(`:3528-3585`, selecting only `traction`/`virtual_work_shape_deriv`/`evolution_mass_balance` — it SKIPS
`closure_shape_deriv`, the T-i Λ_A/Λ_V response). So the face generalized forces are never turned into operator
rows and never weakly restricted — that is why the emitted `S11CB_COUPLING_KERNEL` is bulk-only. WL does the
conversion PY skips: `faceGeneralizedRows` (`mathematica…:1137-1176`) takes the virtual work, forms the
coefficients of the independent virtual displacements, and ADDS `U_FACE`/`EW_FACE` to the operator rows before the
weak pairing. `bulk_kernel_from_density` (`:3177`) / `paired_kernel_from_density` (`:3227`) are dead code — the
§3c-forbidden parallel density-variation route; the emit must not call them.

## What is SUPPLIED (unfalsifiable within this build — state so; do NOT verify against WL)

- The §3c CONTENT verdict (INCLUDE/INCLUDE). You implement it; you do not re-derive whether to include.
- The S11c-a T-a..T-i shape-derivative substrate (the tilted-face geometry, `traction`, `virtual_work_shape_deriv`,
  `closure_shape_deriv`, the flux/mass content) and the §1c `Λ_I(ω)` kernels — supplied objects PY imports; consume
  them, do not recompute.
- The §3c EXTRACTION INSTRUMENT (the weak trial/test restriction already applied to the bulk rows) and its
  prohibitions.

## The changes to make

### (1) §0 clarity pin — spec amendment to `S11c_b_SHARED_PHYSICS.md` §0

Add at the permeability/memory-kernel exclusion (apply verbatim — folded from both legs):

> The §0 exclusion of permeability/memory kernels applies only to kernels derived by the S11c-c curved-bulk
> response solve (`δp = Z·v_bulk`); it does not remove the supplied flat-face kernels `Λ_I(ω)=Λ_I⁰/(1−iωτ_I)` of
> §1c. Those kernels stay symbolic in each §3b face/flux contribution that consumes the S11c-a closure, flux,
> traction, or T-i shape-derivative (no bulk elimination / DtN), and therefore in any §3c block extracted from
> those contributions. The §0 exclusion of a face-parity impedance object is unchanged.

### (2) PY engine — compute the face rows INTO the §3b operator, then restrict the full operator

The operator's face/flux contribution is the generalized FACE FORCE, ⛔ NOT the raw `FACE_FLUX_BOUNDARY_OPERANDS`
bundle. Implement, as an object (not a recipe):

1. **Form the face generalized-force rows** from the consumed T-a..T-i virtual work: the coefficients of the
   independent virtual displacements `δ_vu` and `δ_ve_W` in the supplied virtual work (the S11c-a
   `virtual_work_shape_deriv` / `traction` / `closure_shape_deriv` content), together with the face-flux content
   already carried by the sourced mass-evolution row. **Bind the live `μ_θ` operand** (`MU_THETA_FACE_BINDING`,
   `mu_theta_amplitude`) into the face affinity BEFORE forming those forces; ⛔ do not leave the S11c-a reserved
   `mu_theta_L`/`mu_theta_M` placeholder in the extracted block (T-i keeps `μ_θ` as S11c-b's named operand).
2. **ADD** those computed face rows to the constraint-reduced bulk operator rows (`U_BODY_BALANCE`/`E_W_BALANCE`,
   and the mass-evolution row's face flux) with provenance origin FACE — ⛔ ADDED to, not substituted for, the
   bulk rows, and ⛔ NOT passed through the θ-constraint fold (external virtual work is not the internal
   constraint reaction). Emit them as part of `S11CB_SLAB_OPERATOR` and derive their `S11CB_SLAB_OPERATOR_TERM_ORIGINS`
   entries. Face rows join at the SAME `background_depth` cascade as the operator rows they are added to.
3. **Apply the SAME weak restriction** already used on the bulk `EXPANDED` rows to the full (bulk + face) operator
   rows; that single restriction produces the enlarged `S11CB_COUPLING_KERNEL`. Keep `Λ_I(ω)` SYMBOLIC (no bulk
   elimination / DtN / impedance). The bulk `γ·profile-jet` content already extracted STAYS.
4. ⛔ **Do NOT** multiply the virtual-work density by a second test field (it is already a weak-form density —
   `measure × traction · virtual_displacement`; the instrument acts on the OPERATOR ROWS, not on the density).
   ⛔ **Do NOT** re-contract `evolution_mass_balance` (it is already the `THETA_BALANCE` `EXPANDED` row and is
   already weakly restricted). ⛔ **Do NOT** feed `virtual_constraint`, projection, or `face_shift` into the kernel
   (already consumed / not coupling content). ⛔ **Do NOT** add a `ζ_c` (center-face) channel — §3c's block is
   transverse ↔ `{θ, e_W, u_L}` only (WL's `CENTER_FACE_GENERALIZED_ROW` is in the operator but NOT in the coupling
   extraction).
5. **Rebuild the pairing-based adjointness object over the ENLARGED blocks** (construct both cross-sector blocks
   first, then form the adjointness operand from the complete blocks); ⛔ do not reuse the pre-face adjointness
   object. If the two blocks are adjoint by construction, emit them and state there is no independent second route.
6. **`S11CB_COUPLING_KERNEL_TERM_ORIGINS`**: classify each weakly-restricted term by its CONSTRUCTION origin
   (bulk-energy / face-virtual-work / face-flux / advective), ⛔ NOT by whether a term contains `Λ` or a
   trial-potential symbol (that is a content filter, not provenance).

## Ripples the builder must handle (named as obligations, from both legs)

- **Cache keys.** `KERNEL_BLOCK_CACHE`/`KERNEL_ORIGIN_CACHE` `core_key` currently drops `representative` for the
  Eulerian route; the face bundle is representative-filtered, so include `representative` and every face-fold input
  that affects a result once face rows enter the cache.
- **Depth cascade.** Strong rows build at `min(background_depth, STRONG_ROW_JET_DEPTH)`; the outer weak curl/div is
  performed at the requested coupling depth. Keep the two-stage contract; face rows join at the operator cascade.
- **Tower-depth control.** `task_tower_depth_control` runs `coupling_outer_rows` on `strong_rows` and never
  `build_kernel`; extend it so the face-origin content is in the compared object. PRINT the residual; do not assert.
- **§5c uniform-limit + §5a/§5b controls.** Re-run through the enlarged operator; PRINT operand A, operand B, and
  `A−B`. For the uniform comparison compare only the semantic cross-sector block operands (not labels/source/
  provenance metadata) — a container-shape change vs S11b must not read as a physics residual, and the residual
  value is NOT an acceptance gate.
- **Dead functions.** The emit must not call `bulk_kernel_from_density`/`paired_kernel_from_density`; remove them or
  leave them unused (they encode the forbidden parallel route).

## Script obligations (the three clauses)

1. PRINT computed objects; never state conclusions. 2. PRINT residuals; do not `assert` — emit a difference only
between INDEPENDENT routes; a difference zero by construction is operand theatre (emit the objects and say there is
no second route). 3. Interpretation belongs to the step record.

The face+response coupling terms must be REACHED BY COMPUTATION from the consumed T-a..T-i virtual work and the weak
restriction — ⛔ never typed. The only hand-combination is constructing the energy/ansatz and consuming the SUPPLIED
substrate/kernels.

## Acceptance — construction coverage only (rule-5-clean: NO expected value anywhere)

⛔ **No acceptance criterion references a target kernel term, a term count, a specific structure, a channel, a sign,
a zero/nonzero residual, or agreement with the Wolfram engine.** The builder's job ENDS at compute-and-print; the
cross-engine diff happens on the orchestrator's side, where a mismatch is a FINDING. Do NOT iterate toward any
recorded value. A build is accepted when it:
1. PRINTS the computed folded face operator rows, both complete §3c weak blocks, their computed `TERM_ORIGINS`, and
   the pairing-based adjointness operand (or an honest "no independent second route"), per anchoring and density
   representative, multigraded and dimensioned;
2. shows in the code/dataflow that the weak blocks are obtained from the emitted §3b operator rows and that the
   pre-existing bulk-row builder remains an input (not replaced);
3. keeps every supplied `Λ_I(ω)` symbolic and introduces no bulk-response solve, DtN/impedance, parallel
   density-variation kernel, global projector, or channel filter;
4. recomputes and PRINTS every affected control's operand A, operand B, and `A−B` (no asserted residual);
5. states which inputs are SUPPLIED and unfalsifiable in this build.

## ⚠ Operational (build hygiene — not physics)

⛔ Do NOT run the full package loop (it OOMs ~3.5h). Verify with a targeted smoke: build the coupling kernel for
one `(branch, representative)` case (e.g. via `task_coupling_kernel` on a single case, or `S11CB_PRIMARIES_ONLY=1`)
and PRINT the objects. The in-band full run + heavy controls are DEFERRED (`DEFERRED_HEAVY_RUNS.md`); the
DELIVERABLE is the edited engine (+ the §0 amendment), not a fresh `.out`.

## For the BUILD legs (fresh Claude agent + Grok — Codex-written)

Derive independently; write a script + literal stdout; report payloads as findings WITHOUT using payload membership
or a residual value as a build exit condition. Checks: (a) confirm the source-substrate → folded-operator-face-row
→ weak-restriction provenance (the face rows are computed into the operator, not the raw bundle weak-restricted);
(b) FORM-ablate the virtual-work/traction/`Λ` SOURCE feeding the face rows and REPORT the literal diff (ablating a
non-face key such as projection is not required to move it; bulk-origin terms stay under the face ablation);
(c) delete the face-row feed and REPORT whether the face-origin terms vanish; (d) `Λ` stays SYMBOLIC (no `Z`/DtN/
impedance appears); (e) exact-once consumption — `evolution_mass_balance` and the virtual work are not double-paired
or double-counted; (f) `μ_θ` is bound (no reserved `mu_theta_L/M` placeholder survives in the block); (g) no
acceptance-value leak. Do not edit the working tree.

## Method notes for the DECISION legs
(Superseded — this is v2, already legged. Do not re-leg.)
