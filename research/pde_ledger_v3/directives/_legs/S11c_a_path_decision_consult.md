# Strategic decision consultation — S11c-a T7 reconciliation: which path forward?

You are one of two independent advisors (the other is the other CAS assistant) asked by the orchestrator for
a JUDGMENT on how to proceed. This is NOT a code/directive review — it is a decision consultation. Read the
situation, verify whatever you wish against the cited artifacts, and give a direct recommendation. If none of
the options is right, say so and propose a better one. If you need clarification, or think we should NOT do
something, say it plainly. Disagreement with the other advisor is useful — do not hedge toward an imagined
consensus.

## Context
S11c-a step T7 (interface shape-derivative) is reconciled across two independent CAS engines — SymPy "PY" and
blind Wolfram "WL" — under a strict two-engine method (the disagreement IS the measurement; neither engine
sees the other; a fold that turns a real disagreement into false agreement is the worst outcome). Three
§3c-class engine implementation defects were already found and fixed (PY shifted-trace; WL free-premise
background current; PY current-freezing = "finding #3", where PY had frozen the perturbation current as a
w-constant so the ∂_w δj_w term vanished). Remaining goal: the full cross-engine T7 result (which of 39
tag-families agree/disagree) via a REVIEWED comparator, then a step record that pins the reconciliation schema
for the sibling steps S11c-b..e.

## The situation (all measured this session)
1. **The reviewed comparator's build directive has failed TWO independent-leg reviews.** rev-1 (~14 defects)
   and rev-2 (~8 fundamental defects — corroborated by BOTH a Codex leg and a Grok leg, all orchestrator
   rule-13-verified against the payloads): the case-keying regressed to 0 joins (it discarded the working
   axis-typing that strips WL's `DOF_`/`VIRTUAL_DOF_` prefixes and maps PY's positional DOF/VDOF); several
   per-family extractors named the wrong payload leaf (PY projection is `VALUE` not `EXPRESSION`;
   KINEMATIC/VIRTUAL_CONSTRAINT/FACE_SHIFT have bespoke schemas); the perturbation-current name-fold could
   still manufacture false agreement if implemented like the existing FIELD map (which does
   `AppliedUndef → bare Symbol`, i.e. strips arguments); and one projection instrument
   (`~/.s11_build/S11c_a_run_projection.py`, the window→`Dwin` + integral bridge) was missed entirely. The
   legs effectively supplied the full corrected spec. This is 2 failed directive designs.
2. **Provenance settled:** the committed PY engine is mathematically deterministic (a fresh regen from HEAD
   equals the prior transcript up to `srepr` term-ordering only); the regen will be committed as the PY input.
3. **The one open PHYSICS question — the projection does NOT trivially reconcile post-fix.** Scouting the
   FIXED projection (the window→`Dwin` bridge + an arg-preserving current name-fold + PROPER integral
   linearity):
   - the in-plane current divergence AGREES — the naive 8/8-nonzero there was an integral-linearity
     false-nonzero (∂_i δj_i pulled out of the w-integral vs kept inside), which cancels once w-constant
     factors are pulled out of the integral;
   - a structured residual SURVIVES:
     `(Dwin(0,1) − Dwin(1,0)) · (delta_j_bulk_4(w) − delta_j_bulk_4(x1,x2,x3,{w,t}))` + a
     `delta_rho_4D_bulk_t` density-time term. i.e. PY's perturbation NORMAL current is a function of **w
     only**; WL's is a full **x,w,t** field.
   - rule-13 characterization: WL NEVER differentiates the normal current w.r.t. x1/x2/x3 anywhere in the
     projection operands (0 occurrences), so its x,t arguments appear to be INERT SPECTATORS ⇒ the difference
     is very likely REPRESENTATIONAL (benign), not a real §1b restriction. BUT this residual sits on the
     orchestrator's UNREVIEWED scouting instruments (an ad-hoc current fold + an integral-linearity canon),
     AND "the x,t args are inert, so identify WL's field with PY's w-only current" is exactly the name-map
     trap that hid finding #3 — legitimate only if the args are truly inert per the SPEC, which is a
     from-spec ADJUDICATION question, not a solo call.

## The three options
1. **Adjudicate the projection residual FIRST** — a focused from-spec resolution (a runnable-CAS consult + 2
   independent adjudication legs, exactly how finding #3 was resolved): does spec §1b require the perturbation
   normal current δj_w to carry x,t-dependence, or are WL's x,t genuinely inert spectators? Settle the one
   physics question cheaply, THEN build the durable reviewed comparator to confirm + cover the rest. Risk: the
   residual rests on unreviewed tools, so adjudicating presupposes it is not a fold/canon artifact.
2. **Build the reviewed comparator NOW** — fold the legs' full corrected spec into a rev-3 directive → 2
   directive legs → Codex build → 2 build legs → run. This surfaces the projection residual (and all 39
   families) UNDER REVIEW, removing the "is my scouting an artifact" doubt, and is the durable instrument that
   pins the schema for S11c-b..e. Cost: 3rd directive attempt + the full build pipeline (though the 2 leg
   reviews already did the hard grounding, so rev-3 is largely transcription of a verified spec).
3. **Orchestrator extends the working reconciliation to full 39-family coverage himself** (fix the axis-typed
   keying + per-family extractors + the arg-preserving current fold; verify joins>0), producing the complete
   measured T7 table NOW — but orchestrator-written and therefore UNREVIEWED (preliminary), to be
   adjudicated/productionized after. Fastest to a complete picture; weakest independence.

## What I am asking
- Which option do you recommend, and why — or a better path?
- Is it sound to adjudicate the projection residual (Option 1) while it rests on unreviewed scouting tools, or
  must a reviewed instrument confirm the residual first (Option 2)?
- Is "inert spectator ⇒ identify WL's x,w,t field with PY's w-only current" a legitimate reconciliation or the
  name-map trap in disguise — and how should the current fold treat it so it can neither smuggle agreement nor
  manufacture a false disagreement (e.g. surface the raw residual AND the residual-after-dropping-inert-args,
  and let the spec decide)?
- Anything material I am missing.

## Artifacts you may verify (read-only; do not modify anything)
- SCOUT record: `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §23 (this session), §§17–22 (prior rounds).
- rev-2 directive: `research/pde_ledger_v3/directives/S11c_a_comparator_build_directive.md`; its 2 leg reviews:
  `~/.s11_build/comp_dir_rev2_codex.log`, `~/.s11_build/comp_dir_rev2_grok.txt`.
- Projection bridge: `~/.s11_build/S11c_a_run_projection.py`; reconciliation base:
  `~/.s11_build/S11c_a_reconcile_fixed.py`, `S11c_a_cov_all.py`, `S11c_a_scratch_loader.py`.
- Engine transcripts: PY `~/.s11_build/S11c_a_py_regen.out` (fixed, committed-engine output); WL
  `research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`.
- Spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (§1b projection; §3c trace; §5–§8).

Give a direct recommendation with reasons. Two to four paragraphs is plenty.
