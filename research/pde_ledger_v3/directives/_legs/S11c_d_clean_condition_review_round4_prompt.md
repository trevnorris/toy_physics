# Independent physics review

## Artifact

A two-file packet, **version 4** (round 4 of review-until-clear). Authorship is mixed: v1–v3 were orchestrator-written; v4 was revised by a Codex author. v1, v2 and v3 and their review rounds are preserved at commits `7b38e9dc`, `e5f0dfda` and `c1e96e76`.

- **A.** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_clean_condition.md`
  - **Role:** a physics-bearing proposal/record. It states a symmetry-based "clean condition" under which light
    in this model does not leak from the brane slab into the bulk at linear order. It proposes two requirements
    (R-LEAK-1, R-LEAK-2) and re-scopes S11c-d.
- **B.** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md`
  - **Role:** a physics-bearing builder directive that a Codex builder will execute.
  - **Part 1** computes the complete S11c-b linear operator on perturbations independent of one interface
    coordinate, with FORM controls.
  - **Part 2** is a scope note for a round background with the drain flow live.

## What to check

Report anything that changes **what is computed** or **what may be claimed**. In particular:

1. **Is hypothesis H correct for THIS model's operator?** H is in A §1 (H-planar and H-round). The operator's
   constituents are:
   - `u`: three in-plane components, no `w`-component;
   - `θ`;
   - the independent face variables `ζ_+`, `ζ_-`, including the centre shift `ζ_c`;
   - the tilted-face geometry;
   - the face responses `Λ_A𝒜_s`, `Λ_V V_s`, `Λ_X𝒜_s`;
   - the two anchorings `LAB_HELD` and `MATERIAL_ADVECTED`;
   - the constraint fold (pin B);
   - the bulk fields `φ`/`δp`;
   - the background flow.

   Is any premise missing from A's table P1–P4 that could break H? Is any stated premise wrong? Is the parity
   bookkeeping for toroidal vs spheroidal vs scalar fields at each `(ℓ, m)` correct? Is the claim right that
   rotation about the interface normal brings every oblique perturbation to the `z`-independent form? Check this
   against what the background actually depends on.
2. **Coverage and requirements.**
   - Does A's coverage table (§2) overclaim or underclaim?
   - Are R-LEAK-1, its linear falsifiers (F1, F2, F2b, F5, F6) and its nonlinear gates (N1, N2) correctly
     scoped? Consider spin carriers, chiral trapped shear, and second order.
   - Is "transfers to the calibrated model" justified as stated?
   - Are the rejected options in §4 fairly characterised?
   - Does A contradict any governing statement it cites?
3. **Citations.** Every repo citation in A and B must match the verbatim excerpts in the two grounding files
   `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md` (for A) and
   `/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md` (for B),
   and the files themselves. Report any mismatch, with both quotes.
5. **The round-3 fold.** The disposition (handed below) has a round-3 section listing ten accepted findings
   (R3-1 … R3-10) and the resolution owed for each. The v4 change log says what was changed. For each finding,
   check whether v4 actually resolves it. Also check:
   - whether any v4 change introduced a new defect;
   - whether any round-1 or round-2 resolution regressed.
4. **Directive B, as a builder packet.**
   - Does B name the **object**, not a recipe?
   - Does B carry the governing operator accurately? Check it against the S11c-b spec and step record.
   - Does B leak an expected value or the shape of the answer to the builder — through its wording, tag-naming
     rules, file names, or control descriptions?
   - Would its FORM controls K1/K2 expose a block that is zero **by construction** (i.e. unimplemented), and are
     they genuinely FORM, not COEFFICIENT, changes?
   - Does "construct for all fields at once" prevent a tautological result?
   - Are the model point and every freeze declared?
   - Is the blind Wolfram requirement sound?
   - Is it feasible on a 30 GB machine, given the S11c-b Wolfram build's recorded memory history?
   - Does it keep the builder in its lane (build → run → report → stop)?
   - Is Part 2 scoped so that the missing drain-flow specification is surfaced rather than invented by the
     builder?

## What you are handed

- The two artifacts above.
- The grounding files `research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition.md` and
  `research/pde_ledger_v3/directives/_measurements/S11c_d_zinvariant_operator_blocks_directive.md`.
- The disposition, rounds 1–3: `research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_review_disposition.md`.
- The v4 change log `research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v4_changelog.md` and the
  v4 author's scripts and stdout in `research/pde_ledger_v3/directives/_measurements/S11c_d_clean_condition_v4_author_scripts/`.
- `AGENTS.md` and `scripts/s11c_guarded_run.py` (the execution guard B requires).
- The governing sources, all under `/var/projects/toy_physics/`:
  - `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md`
  - `research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md`
  - `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
  - `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
  - `research/pde_ledger_v3/scripts/S11c_b_exports.py`
  - `research/pde_ledger_v3/V3_STEP_PLAN.md`
  - `docs/toy_model_ontology_summary.md`
  - `docs/native_light_em_and_vortex_throat_interpretation.md`
  - `research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md`
  - `research/pde_ledger_v3/directives/S11c_d_SCATTERING_FORM_AMENDMENT.md`
- The rest of the repository is readable.

## Required method

**The artifact is a DOCUMENT packet.** Read the governing sources first: the S11c-b spec field content and
energy-basis construction, the S11c-b step record, and the S11c-b SymPy engine's operator construction. Form your
own view of what symmetry the operator, as actually built, has. Only then read A and B. Quote both sides for
every finding.

**Derive it yourself, and show the computation.** A prose derivation is worth nothing. For any physics claim you
make, write your own script and save both the script and its literal stdout to named absolute paths under
`/tmp/s11cd_clean_review_r4_<your-leg-name>/`, and report those paths. Without them, the claim is discarded. Physics
claims here include:
- that a reflection does or does not map the field sectors as A says;
- that a face-response or anchoring term is or is not reflection-even;
- that a proposed control term is parity-odd.

Where possible, build a small symbolic model of the relevant terms. Better still, read the actual term
structure off the S11c-b SymPy engine or its exports, and test it.

Run any test of the existing engines on a **copy** under `/tmp`; ⛔ never modify the working tree. For any CAS
kernel:
- ⛔ wrap every kernel run in `timeout 600`; a 600 s hit is a FAILED ablation, so report it and move on;
- ⛔ never raise the timeout, and ⛔ never run more than one Wolfram kernel at a time (the licence has two seats);
- ⭐ save every script and its literal stdout to named absolute paths, and report those paths.

⚠ A suspended S11c-d job (PID 4097233) holds memory on this machine. ⛔ Do not signal, resume or touch it. Keep
any computation small.

## Physics filter

Report a finding only if it catches a way the physics could be wrong, or a way the packet would make the builder
compute or the record claim the wrong thing. Do not report "the script would be wrong on a different input."

## Ablation sandbox

Copy anything you execute or modify to `/tmp/s11cd_clean_review_r4_<your-leg-name>/` and work on the copy.
⛔ Never modify the working tree.

## Bounds

Write your report and exit. ⛔ Do not spawn agents or build `run_all` / `watch` / supervisor orchestration.
Iterating to clearance is the orchestrator's job, not the leg's.
