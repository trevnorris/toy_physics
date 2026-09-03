# Builder directive — PY constraint-fold of the S11c-b slab U-momentum row (spec-pin B)

## Status (read first)

Orchestrator-written builder decision list for the follow-on of the SETTLED spec-pin (B). **v2** — v1 was
REJECTED by both decision legs (Codex + Grok, independently): v1 folded §1c:143's `(uniform linearisation)` as the
constraint (a rule-17 freeze of `W_bg`), double-counted the mass row, changed `operator_from_density` in place
(breaking #88), truncated the coupling/tower depth cascade, and leaked a prior measurement. This v2 folds both legs'
convergent, computation-backed replacement text (leg logs `~/.s11_build/S11c_b_constraint_fold/decision_{codex,grok}.log`,
`.log` gitignored). The four load-bearing claims were orchestrator-verified against the code (rule 13). v2 goes to
Codex build and its own two BUILD legs; it is NOT re-legged (rule 7: one pass, fold, go).

**What is SETTLED (do not re-litigate — it has its own 2 decision legs + a Codex verification):** the S11c-b slab
U-momentum row is the CONSTRAINT-REDUCED in-plane momentum balance, **not** the raw held-fixed `δU/δu`. Pinned (B).
Records: `directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`;
`directives/S11c_b_jet_depth_spec_pin_decision.md` (FOLDED VERDICT). Decisive governing text:
`S11b_SHARED_PHYSICS.md:337,341-345,356-357,426` (the material virtual constraint `δ_vΣ_mat=0`; "Do NOT vary U with
θ held fixed"; "the **same** multiplier supplies the in-plane restoring force and the thickness term"; the mass
balance with `J_±` is a SEPARATE evolution equation restored AFTER the variation; the `−∇(δU/δθ)` convention check).
⇒ **WL is already correct; PY is the engine that must change.** WL is NOT touched by this build.

## The binding constraint is the MATERIAL, non-uniform object — NOT the uniform three-term formula

`S11c_b_SHARED_PHYSICS.md:143` displays two objects on one line: the binding constraint `δ_vΣ_mat = 0`, and then,
explicitly labelled, the **`(uniform linearisation)`** `δ_vθ + δ_ve_W + ∇_x·δ_vu = 0`. The three-term formula is the
UNIFORM special case; it is **not** the constraint on the non-uniform background and must not be substituted for it
(doing so freezes `W_bg` inside the constraint — rule 17). The actual object is already built and imported: the
S11c-a substrate `virtual_constraint` (T-g), constructed per branch / density representative / route from the
material mass (`S11c_a_interface_geometry_sympy_audit.py` `virtual_constraint_route`, ~:1207-1259; PY imports it via
`SHAPE_SUBSTRATE_KEYS` `:335`). Its live-`W_bg` form carries non-uniform coefficients that the uniform formula
drops. Construct the fold from that imported operand (read its coefficients from the operand itself); do not
re-derive the constraint and do not substitute the uniform formula.

## What is SUPPLIED (unfalsifiable within this build — state so in the report; do NOT verify against WL)

- The material virtual constraint `δ_vΣ_mat = 0` as the imported S11c-a `virtual_constraint` operand (per branch,
  representative, route), and the sourced mass balance `∂_tΣ + ∇_x·(Σv) = −(J₊+J₋)` as the imported S11c-a
  `evolution_mass_balance` operand.
- S11b's binding virtual-displacement rule and the convention check (`S11b_SHARED_PHYSICS.md:337,341-345,356-357,426`).
- The pin itself (object (B)). You implement it; you do not re-adjudicate it.
- `μ_θ ≡ (δU/δθ)|_{u,e_W,… held fixed}` (§1c:128,131) — the CONSTITUTIVE derivative, computed held-fixed BEFORE any
  constraint. PY already computes it (`operator_from_density` returns `mu_theta_amplitude`) and emits it separately
  (`MU_THETA_FACE_BINDING`, `MU_THETA_OPERATOR`). It stays on the held-fixed path, unchanged.

## The target artifact

`scripts/S11c_b_brane_operator_sympy_audit.py` (the SymPy engine) and a one-paragraph amendment to
`directives/S11c_b_SHARED_PHYSICS.md` §3b. No other file changes.

## (1) §3b spec amendment

Add exactly this paragraph to `S11c_b_SHARED_PHYSICS.md` immediately after the §3b prose that ends "…bulk-energy vs
face/flux vs advective, per anchoring and density representative, multigraded and dimensioned." and before the
`⇒ S11CB_SLAB_OPERATOR …` emit block (apply verbatim — it is the pinned text, distilled from both decision legs):

> **The slab momentum and thickness rows are the CONSTRAINT-REDUCED equations under S11b's binding material
> virtual-displacement rule, not the held-fixed variational derivatives.** First compute the constitutive operand
> `μ_θ ≡ (δU/δθ)|_{u,e_W,… held fixed}`, before applying any constraint. Then, for each anchoring, density
> representative, and construction route, obtain the applicable material virtual-constraint operand `δ_vΣ_mat = 0`
> from the supplied S11c-a substrate (§1c; the three-term relation on `S11c_b_SHARED_PHYSICS.md:143` is only the
> uniform linearisation and does not replace it on the non-uniform background), solve it for `δ_vθ`, substitute that
> into the internal virtual variation, and extract the coefficients of the independent `δ_vu` and `δ_ve_W` after
> compact-support in-plane integration by parts (§3c interior convention). The same computed `μ_θ` generates both
> the in-plane and the thickness reaction (`S11b_SHARED_PHYSICS.md:341-345`). The `u`-row is therefore not
> `δU/δu|_{θ,e_W fixed}`. In the convention check the physical in-plane restoring-force contribution is `−∇(δU/δθ)`
> (`S11b_SHARED_PHYSICS.md:426`) — a check the variation must be able to fail, not a term to type, and not to be
> confused with the sign of the internal Euler–Lagrange residual an engine stores. The `θ` slot of the slab operator
> is the supplied sourced mass-evolution equation `∂_tΣ + ∇_x·(Σv) = −(J₊+J₋)` (§1c), restored after the variation
> (`S11b:356-357`), not `μ_θ = 0`. `S11CB_MU_THETA_OPERATOR` remains the separate held-fixed constitutive operand
> used by the face affinity; it is neither the mass-evolution row nor the vector momentum reaction.

## (2) PY engine — keep the raw constructor; add a computed constraint-fold layer

- **Keep `operator_from_density` (`:1968-2062`) as the held-fixed RAW energy-EL constructor**, including its return
  of `mu_theta_amplitude`. This preserves the frozen `committed_strong_rows` (`:2123-2175`) and the #88 raw-EL
  reference (`S11c_b_88_blast_radius.py` imports `operator_from_density` directly). ⛔ Do NOT fold the constraint
  inside `operator_from_density`, and ⛔ do NOT fold it inside `committed_strong_rows`.
- **Add a separate constraint-fold constructor, called explicitly by `build_operator`**, receiving the raw internal
  `U`/`e_W` EL amplitudes, the same `mu_theta_amplitude`, and the selected `virtual_constraint` operand for
  `(branch, DELTA_W, representative, route)`. Apply `material_pullback` for the MATERIAL route BEFORE computing that
  route's raw EL, `μ_θ`, and reduced rows.
- **Construct the reduced rows GENERICALLY**: form the internal virtual variation, solve the supplied constraint for
  `δ_vθ`, substitute, and extract the coefficients of the independent virtual fields with a GENERIC interior-IBP
  coefficient extractor. ⛔ Do NOT encode a `U`-specific gradient formula or an `e_W`-specific multiplier formula —
  the reaction (including its non-uniform coefficients) must be REACHED BY COMPUTATION from the selected constraint
  operand and the same `mu_theta_amplitude`, never typed.
- **Second construction (transcription check only):** also build the reduced variation via a Lagrange-multiplier
  route and PRINT both operands and their symbolic difference; do NOT assert it. ⚠ **CORRECTED (build-leg finding,
  Claude agent + orchestrator CAS proof):** for a constraint LINEAR in `δ_vθ` these two routes are ALGEBRAICALLY
  IDENTICAL (`λ = −ε·μ_θ/(∂C/∂δ_vθ)`, and back-substitution cancels the `δ_vθ` term), so their residual is `0` by
  construction for ANY constraint — including a wrong one. It is therefore a TRANSCRIPTION-CONSISTENCY check between
  two coded implementations, ⛔ **NOT an independent physics route, and must not be cited as one.** The genuine
  independent check of the reduced rows is the blind cross-engine (WL) comparator (deferred `row_residual`); within
  PY the row is protected by the constraint-SOURCE ablation (a wrong constraint MOVES the rows), not by this
  residual. Emit it under a name that says so (`CONSTRAINT_FOLD_TRANSCRIPTION_RESIDUAL`).
- **Depth cascade (rule 17 — the reaction raises the jet order at every stage; nothing may be capped below its
  input's order).** `μ_θ` is already a one-divergence object; its gradient in the `U`-row raises it again, and §3c
  differentiates the §3b `U`-row once more in `build_kernel` (`operator_divergence`, ~:2616). Set
  `STRONG_ROW_JET_DEPTH = 3`, `COUPLING_JET_DEPTH = 4`, `DEPTH_CONTROL_JET_DEPTH = 5` (`:59-61`); extend the
  background-profile jet tables `BACKGROUND_PROFILE_JETS` (`W_PROFILE_JETS`, `M_PROFILE_JETS`) through the control
  depth and the wave-field `DERIVATIVE_MAP` (`:648-673`) through the needed spatial jet; keep the chosen / shallower
  / deeper calls in `task_tower_depth_control` (`:3740+`) structurally distinct after the shift (strong 3, coupling
  4, control 5).

## (3) PY engine — the `θ`-row is the imported mass-evolution operand; `μ_θ` stays separate

- Select the existing S11c-a `evolution_mass_balance` substrate case `(branch, DELTA_W, representative)` as the
  emitted `θ` strong row. That operand is ALREADY the complete sourced balance
  (`DENSITY_TIME + VELOCITY_DILATATION + BACKGROUND_ADVECTION + TRUE_AREA_FACE_FLUX`, S11c-a:1313-1318). ⛔ Do NOT add
  `ADVECTIVE_MASS_OPERAND` to it — that double-counts `BACKGROUND_ADVECTION`. Retain `ADVECTIVE_MASS_OPERAND` only as
  the named background-advection origin/control operand, and select the matching `evolution_term_origins` case for
  the mass-row provenance.
- Select the `DELTA_W` DOF case explicitly (the substrate carries both `DELTA_W` and `ZETA_C`); the slab DOFs are
  `{u,θ,e_W}` (§3b:274), so the thickness DOF is `DELTA_W`.
- Keep `mu_theta_amplitude`, `MU_THETA_FACE_BINDING`, and `MU_THETA_OPERATOR` on the original held-fixed constitutive
  path, unchanged.

## (4) Provenance, controls, diagnostics, and every consumer of the changed rows

- **Provenance stays `BULK_ENERGY`.** The constraint reaction is generated from a bulk-energy constitutive
  derivative; keep the constrained internal `U`/`e_W` contributions under the existing `BULK_ENERGY` origin (matching
  WL, where the reaction is entirely in `BULK_ENERGY`). ⛔ Do NOT add a new top-level `CONSTRAINT_REACTION` origin —
  that is over-reach and inconsistent with WL. In the per-term `energy_origins` loop (`:2289-2308`), compute each
  term's raw EL and its reduced `U`/`e_W` contribution through the same fold constructor; ⛔ do NOT install the
  global mass row once per energy term. Record accumulation / advection / face-flux once from the selected
  `evolution_term_origins` / face substrate. Update `COUPLING_KERNEL_TERM_ORIGINS` to consume the reduced
  bulk-energy rows and the mass origins; do not silently drop a new origin.
- **Diagnostics.** Print `RAW_BULK_U_EL` (the inertia-free held-fixed bulk EL) and `CONSTRAINT_REACTION_U` (produced
  by the virtual-constraint elimination constructor) as engine-local diagnostics; if their difference is printed,
  subtract like-for-like internal bulk operands so kinetic/face content cancels by construction. These are NOT new
  normative §3b provenance categories and are NOT exported as cross-engine primaries without a separately reviewed
  schema change.
- **Corruption/ablation propagate at the SOURCE.** `corrupt_material_constraint` and the §5b form ablations must
  perturb the selected `virtual_constraint` and `evolution_term_origins` SOURCES BEFORE the reduced rows / mass row
  are constructed (§5a:412-416 requires propagation through every construction). ⛔ Do NOT keep the current path that
  only zeroes the already-built detached `ADVECTIVE_MASS_OPERAND` (`:2267-2268`).
- **Every consumer of the changed rows must read the canonical constrained `build_operator` path**: the homogeneity
  inventory (`:3913`), the representation/form controls, the tower-depth control, `uniform_transverse_dispersion`
  (`:3872`, which today calls the raw constructor directly and reads `U` `EXPANDED[1]` — a computed transverse
  elimination must not inject a longitudinal restoring force there), and the §5c uniform-limit adapter (`:3891`),
  which needs a semantic row adapter now that the `θ` slot is mass-evolution: the load-bearing (B) regression is
  `U` and thickness vs the S11b `inplane_eom` / `thickness_eom`; the S11b third slot is the VIRTUAL constraint, not
  mass-evolution, so that residual component is not like-with-like and must not be treated as a build target.
  ⛔ `committed_strong_rows` and the #88 raw constructor stay UNCHANGED.

## Script obligations (the three non-negotiable clauses)

1. The script PRINTS computed objects; it never states conclusions. Any `emit`/`Print` payload is a CAS object.
2. PRINT residuals; do not `assert` them. The two-route (elimination vs Lagrange-multiplier) reaction residual is a
   genuine independent-route check — emit its symbolic difference; do not assert it zero.
3. Interpretation belongs to the step record, not the script.

The ONLY place physical symbols may be combined by hand is in constructing the energy/ansatz. The constraint, the
mass balance, and `μ_θ` are SUPPLIED objects; the reduced rows and their coefficients (including the reaction's
non-uniform coefficients and any sign) must be REACHED BY COMPUTATION from those supplied objects and S11b's rule.
Every control re-enters the chain at the energy/constraint source, never at a result.

## Acceptance — construction coverage, not row value (rule-5-clean: no expected value anywhere)

⛔ **There is NO acceptance criterion referencing a target row content, a jet-atom count, a specific background-jet
atom, a sign, or "agreement with the Wolfram engine".** The builder's job ENDS at compute-and-print. The
cross-engine diff (`scripts/S11c_b_row_residual.py`, deferred to the ≥64 GB box) happens on the orchestrator's side,
where a mismatch is a FINDING, not a build failure. Do NOT iterate the engine toward any recorded value.

A passing build PRINTS, per anchoring and density representative, multigraded and dimensioned: the canonical
constrained `S11CB_SLAB_OPERATOR` (constraint-reduced `U`/`e_W` rows; the selected sourced mass-evolution `θ`-row);
`S11CB_SLAB_OPERATOR_TERM_ORIGINS` with the reduced rows under `BULK_ENERGY` and the mass origins recorded once;
the separate held-fixed `S11CB_MU_THETA_OPERATOR` (unchanged); the `RAW_BULK_U_EL` / `CONSTRAINT_REACTION_U`
diagnostics and both independently-constructed fold operands and their residual; `S11CB_COUPLING_KERNEL` /
`…_TERM_ORIGINS` rebuilt from the reduced operator at the raised depth; and every §5 control re-run through the
changed constructors and PRINTED with no asserted residual, sign, jet membership, atom count, or cross-engine
equality. Source inspection must show both reactions originate from the selected `virtual_constraint` operand and
the same computed `mu_theta_amplitude`, not from a typed `U`/`e_W` reaction formula. State in the report which
objects are SUPPLIED (the constraint, the mass balance, `μ_θ`, the pin) and therefore not verified by this build.

## ⚠ Operational (build hygiene — not physics)

⛔ Do NOT run the full package loop — the full engine run OOMs (~3.5h, memory-bound). Verify your edit with a
targeted smoke test: import the engine and build one representative case (e.g. `S11CB_PRIMARIES_ONLY=1`, one
`(branch, representative)`), or exercise only the changed constructors on a single case, and PRINT their objects.
The in-band full run, the heavy equivalence controls, and the `.out` regeneration are DEFERRED to the ≥64 GB box
(`DEFERRED_HEAVY_RUNS.md`); the DELIVERABLE of this build is the EDITED ENGINE, not a fresh `.out`.

## Follow-on (NOT part of this build — record it, do not execute it here)

- **#88 re-adjudication.** #88's ENERGY-BASIS-COMPLETION disturbance measurement STANDS (a statement about the
  unconstrained EL, which this build preserves as `operator_from_density` / `committed_strong_rows`); the
  KINETIC-family / full-row adjudication is REDONE after PY carries the constraint reaction; #88's `θ` result is a
  disturbance of `MU_THETA_OPERATOR`, not of an independent `θ` equation. Directive
  `directives/S11c_b_88_blast_radius_build_directive.md`, record `_measurements/S11c_b_88_blast_radius_result.md`.
  This build does not touch #88; it only must not contradict it (hence `operator_from_density`/`committed_strong_rows`
  stay the raw reference).
- **Cross-engine confirmation.** After this build, both `.out` are regenerated (WL needs the ≥64 GB box,
  `DEFERRED_HEAVY_RUNS.md`) and `scripts/S11c_b_row_residual.py` compares the two constraint-reduced `U`-rows.

## For the BUILD legs (fresh Claude agent + Grok — this artifact is Codex-written)

Derive independently. The load-bearing checks: (a) the reduced `U`/`e_W` rows are reached from the IMPORTED
`virtual_constraint` operand (form-ablate the constraint SOURCE — a wrong constraint must move the reaction; a
uniform-frozen constraint must be detectably different); (b) `operator_from_density` and `committed_strong_rows`
still return the raw held-fixed EL (so #88's reference is intact); (c) the elimination and Lagrange-multiplier routes
are genuinely independent (one-sided corruption moves only one route); (d) the depth cascade retains, not truncates,
the reaction's jets (ablate a depth cap and show a jet is lost); (e) the `θ`-row is the imported
`evolution_mass_balance` with no double-counted advection; (f) no acceptance-value leak. A prose re-derivation is
discarded — write a script and save its literal stdout.
