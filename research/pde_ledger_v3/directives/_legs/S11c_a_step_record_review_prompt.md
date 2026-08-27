# Independent review — the S11c-a step record (physics-bearing prose)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_a_interface_shape_derivatives.md`

The orchestrator-written step record for S11c-a (the tilted-face interface shape derivatives). It is the
INTERPRETATION layer: it states conclusions the two CAS engines are forbidden to state, so every conclusion
must be faithful to what the engines actually COMPUTED. Its central claim is strong — **"the two
independently-built engines AGREE on T7; every nonzero cross-engine residual is a representational identity or
the known CONORMAL §3c form; no genuine physics disagreement survives."** Your job is to try to BREAK that
claim: find any residual the record calls representational that is actually a dropped-dependence physics
finding, any result stated that no engine emitted, any overstated/understated conclusion, or a mis-stated
count, limit, or provenance.

## Required method — DOCUMENT review, sources first
Read the SOURCES OF TRUTH and form your own view BEFORE reading the record; then check each claim against the
emitted objects. Where a claim is checkable against the run, VERIFY it yourself (grep the tag / reduce the
residual) and quote both the record and the emitted object. Where you derive or reduce anything, **save the
script AND its literal stdout to named absolute paths under `~/.s11_build/` (prefix `srr_`) and report them —
a physics claim with no computation behind it is discarded.** ⛔ Do not modify the working tree.

## Sources of truth (the authority is the emitted objects, not the record)
- Spec `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`.
- The committed run `~/.s11_build/comparator_run.out` (6714 per-case `operand_A`=PY / `operand_B`=WL /
  `A_minus_B` triples + per-family `ACCOUNTING`; `RUN_ACCOUNTING` at the end).
- The two engines: `scripts/S11c_a_interface_geometry_sympy_audit.py` and
  `mathematica/S11c_a_interface_geometry_mathematica_audit.wl`, and their transcripts under `scripts/out/`
  and `mathematica/out/`.
- The comparator `scripts/S11c_a_cross_engine_comparator.py` (reuse its `combine_bound_integrals` /
  collapse helpers if you reduce residuals; you may also write your own).
- `CLAUDE.md`. You have the full repo; there is no do-not-read list.
The record CITES measurement scripts in `directives/_measurements/S11c_a_*`; you may inspect them, but the
authority is the emitted objects — reproduce the record's claims independently, do not take the measurement
files on trust.
⚠ You are deliberately NOT handed the comparator build directive or the orchestrator's working notes: a record
can be internally consistent and still misrepresent what the engines computed — that is what this review
catches.

## Claims to verify (each against the run; try to refute)
1. **μ_θ representational identity.** The record says TRACTION/CLOSURE_SHAPE_DERIV/VIRTUAL_WORK_SHAPE_DERIV
   residuals are exactly `Σ coeff·(mu_theta − mu_theta(x1,x2,x3,time))` and collapse to 0 for all 48 cases,
   and that WL NEVER differentiates μ_θ. Verify: grep the run for any jet-suffixed `mu_theta` or
   `Derivative(...mu_theta...)` (the record claims zero of each); reduce at least one TRACTION residual
   yourself and confirm the collapse. Does the spec (§1c/§2a/§3a/§3c/§4 T-i) actually make μ_θ's args inert,
   or is there a supplied in-plane/time dependence PY dropped?
2. **δρ_4D representational identity — the load-bearing one.** The record says the PROJECTION density-time
   residual is `delta_rho_4D_bulk_t` (PY, a single named symbol) minus WL's explicit §3a expansion, and
   vanishes for all four (anchoring × representative) combinations under the spec relation
   `δρ_4D = ρ_4D,bg^{0,α}·θ`, INCLUDING the MATERIAL_ADVECTED × RHOBR_CONSTANT advection. **Adversarial
   check:** the whole benign verdict for MAT×RHOBR flips to a FINDING if PY's `delta_rho_4D_bulk` means the
   densification about the *local advected* background (`ρ_bg^{0,M}·θ`, which drops the advection) rather than
   the perturbation of the full field. The record closes this by pointing at PY's EVOLUTION family
   (EVOLUTION_MASS_BALANCE + EVOLUTION_TERM_ORIGINS agree with WL 8/8 to 0, PY computing an explicit
   `BACKGROUND_ADVECTION` term). **Verify that closure independently:** confirm the EVOLUTION A_minus_B are 0,
   confirm PY's engine computes the advection concretely there, and decide from the spec whether that makes
   PY's projection symbol the full δρ. Is the record right, or is MAT×RHOBR a PY-side dropped-advection
   finding?
3. **The three §3c-class fixes** (shifted-trace PY `c36beac4`, background-current WL `6fae82b8`,
   current-freezing PY `49b5c525`) — confirm each is described faithfully (which engine, what defect, spec
   basis) and not overstated.
4. **The control battery.** HOLD controls (REP_INVARIANCE_RESIDUAL, UNIFORM_LIMIT_RESIDUAL): the record says
   0 genuine nonzero cross-engine (with the 240 UNIFORM unjoined being a pairing/coverage gap, not physics).
   BITE controls (CONTROL_FORM_RESIDUAL, CONTROL_INDEPENDENCE_RESIDUAL): the record says every object bites in
   ≥1 case (no dead control). Spot-check both claims against the run; is the "0 genuine nonzero" for the HOLD
   controls right, and is any object's form-ablation control actually dead?
5. **Counts / limits / provenance.** 39 families, 6714 triples, 47 PY tags; the standing `O(v₀|q_n|/ω)`
   limit; the provenance note (the committed PY transcript `afdc8158` is a regen equal up to `srepr` ordering,
   not the scratch run). Flag any number or limit that the run/`git` does not support.

## Physics filter
Report a finding only if it catches the record misrepresenting the physics or the computed result — an
overstated agreement, a residual mis-adjudicated as representational, a result no engine emitted, a dead
control called live, a dropped standing limit — not prose style.

## Output
Per item: verdict + the emitted object (tag + value) or your script path + literal stdout, with exact quotes
from the record and the source. End with a one-line bottom line: is the step record faithful to the computed
physics as written, or the specific corrections it needs.
