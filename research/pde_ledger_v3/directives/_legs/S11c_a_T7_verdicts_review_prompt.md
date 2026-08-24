# Independent review — S11c-a T7 step-1 ADJUDICATION VERDICTS

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_T7_adjudication_verdicts.md`
(+ twin `…/_measurements/S11c_a_T7_adjudication_verdicts.md`).

## What this is
Two independently-built engines (SymPy + Wolfram) compute the S11c-a interface shape-derivatives and
DISAGREE on case structure. A step-0 matrix classified the divergences; this step-1 document ADJUDICATES
the four deep families against the CLOSED spec, deciding **which engine is correct** for each — and thus
which engine gets patched. Your job: are these verdicts **correct**, or is any engine being told to
"fix" something where it is actually right (or the spec is silent)? A wrong verdict here patches a correct
engine and corrupts the reconciliation for the whole S11c family.

## Sources of truth (form your OWN view BEFORE trusting the verdicts)
- PY stream `/home/trevnorris/.s11_build/S11c_a_sympy_engine.out` (sha `6386471…`); WL stream
  `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out` (sha `82062bd…`).
- Engines: PY `scripts/S11c_a_interface_geometry_sympy_audit.py` (`9b6438fa`); WL
  `mathematica/S11c_a_interface_geometry_mathematica_audit.wl` (`ddecdbc2`).
- Spec (CLOSED) `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`) — read §2b, §2d,
  §3b, §4 (esp. T-d line ~419), §5b, §5c yourself.
- The step-0 matrix `S11c_a_T7_adjudication_matrix.md` for context (already 2-legged; do not re-review it).

## Required method — COMPUTE + READ THE SPEC, do not accept prose
⛔ A prose re-derivation is worth nothing. Write your OWN scripts (save script + literal stdout to named
paths under `/home/trevnorris/.s11_build/`; ⛔ do NOT write in the working tree or modify any engine).
Do NOT run the orchestrator's scripts as your check — inspect them only after your own number is in hand.

Verify each verdict INDEPENDENTLY:
1. **VERDICT B (density → PY correct).** Re-derive: within the engine emitting both reps, diff the two
   rep expressions per (branch,face,dof). Confirm WL KINEMATIC_BALANCE / RELATIVE_FLUX are rep-IDENTICAL
   and PY PROJECTION dynamic/shape/residual DIFFER (static identical). Then read §3b:351 — do flux/kinematic
   use `ρ_m` (a rep-independent bound constant) rather than the density representatives ρ_4D/ρ_br? Read
   §2b:229 and §4:398. Is "PY correct, WL fix both ways" right — or is WL's rep-independence a WL BUG
   (it forgot to weight by the rep-dependent density) rather than correct physics? Decide from the spec's
   flux/projection definitions, not from the engines agreeing.
2. **VERDICT C (virtual work → WL correct, PY missing).** Confirm WL's 8 off-diagonal (phys≠virt)
   virtual-work cases are genuinely nonzero AND that the FULL object (traction·δ_vx·measure, not only
   `SHAPE_DERIVATIVE.EXPRESSION`) carries physical-DOF dependence — else WL's physical×virtual matrix is
   partly redundant and "PY missing real physics" overstates it. Read §4 T-d (~line 417-419): does "which
   pairings occur is part of the computation" require emitting the full matrix, or is the diagonal a valid
   computed result? Decide whether PY's hardcoded diagonal (`:919-924`) is a defect.
3. **VERDICT H.1 (coverage → PY correct, WL missing 5).** Confirm WL omits form-ablation + uniform-limit
   for {EVOLUTION_TERM_ORIGINS, PROJECTION_STATIC/DYNAMIC/RESIDUAL/TERM_ORIGINS}. Read §5b:482-493 and
   §5c:495-499: do "every T-object / each S11c-a object" REQUIRE covering the projection sub-operands
   separately, or is ablating PROJECTION_SHAPE_DERIV sufficient (⇒ WL is fine and the spec is
   under-specified)? Say which, with the spec line.
4. **VERDICT BG (background state → WL correct, PY under-emits).** Confirm §2d:251 requires the boundary
   loads / zero-conditions IN `S11CA_BACKGROUND_STATE`, that WL emits them, and that PY does not (and does
   not carry them under another field name). Is this a real spec-required field gap, or a permissible
   placement choice (PY in ADMISSIBILITY_PREMISE)?
5. **Mis-assignment hunt.** Is any "PY correct" actually WL correct (or vice versa)? Is any item I called
   a SHALLOW reconciliation actually engine-physics (or the reverse)? Is BACKGROUND_DENSITY_MAP's PY branch
   axis real (does the density map depend on the anchoring branch per §2b/§2c)?

## Physics filter
Report a finding only if it changes a VERDICT (which engine is correct), the FIX direction, or the
shallow/deep classification — i.e. something that would patch the wrong engine or miss a required change.

## Rule-2 / rule-5
- rule 2: any verdict claim without a command behind it in the twin? Flag it.
- rule 5: does the document leak an expected OUTPUT VALUE (not just "which engine is correct" — that is
  the adjudication's job — but a specific target expression the patched engine must reproduce)? It must
  not; the patch directives will specify what to compute, not the answer. Flag any leak.

## Output
A short verdict per family: CONFIRM or the exact correction (which engine is actually correct, with your
independent computation + spec line). Name your script + stdout paths. Focus on catching a verdict that
would patch a correct engine — that is the failure that corrupts the reconciliation.
