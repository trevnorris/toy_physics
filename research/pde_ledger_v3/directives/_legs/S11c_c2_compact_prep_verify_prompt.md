# Compact-prep verify — are this session's state docs + memory correct and in the clear?

You are Codex, verifying that the session's state records are accurate before a compaction. ⛔ Document-only; ⛔ do
NOT modify the tree, run no CAS. Read the sources + `git log`, confirm each claim, flag any over-claim/error with
`file:line`. Working dir `/var/projects/toy_projects`? no — `/var/projects/toy_physics`; paths are under
`research/pde_ledger_v3/` unless absolute. Relevant commits: `git log --oneline -8` (expect
`a6737a27` STATUS+N6-design, `cfb2494c` F/G-deferred, `30d4b72d` §5c-corrected, and earlier `6f8dbd34` CLAUDE.md).

## The claimed arc (verify each against the commits + records)
1. **CLAUDE.md CAS-authorship rule** committed `6f8dbd34`: the orchestrator NEVER authors a CAS instrument — Codex
   writes it, reviewed under G1; the biased `verify_F`/`verify_EG` were the E over-clear. Verify the CLAUDE.md
   E1/G4/L-CAS wording + gate record `_measurements/CLAUDE_md_cas_authorship_gate_record.md`.
2. **E/N6 = a §5c SPEC DEFECT** (⛔ not a `representation_pullback` fix): §5c mis-specified N6 as **cross-anchoring**;
   the real N6 (parent S11c-a §5a + sibling c1 §5a) is **Eulerian-vs-material-coordinate within a FIXED anchoring**;
   corrected `30d4b72d`. The current `REP_INVARIANCE_RESIDUAL` is the wrong object; **the control was NEVER RUN**.
   Verify §5c (`directives/S11c_c2_SHARED_PHYSICS.md` ≈303-370) + `_measurements/S11c_c2_N6_spec_adjudication_record.md`.
3. **F/G DEFERRED** committed `cfb2494c`: retired `verify_F`/`verify_EG` F/G conclusions **WITHDRAWN**; numeric
   re-grounding **owed**; increment VALUES unaffected; full-symbolic **tractability wall**. Verify
   `_measurements/S11c_c2_FG_regrounding_deferred.md` (and that the committed `scripts/S11c_c2_FG_diagnostic_sympy.py`
   is marked reference-only / not a cleared deliverable).
4. **N6 build design** `_measurements/S11c_c2_N6_build_design.md`: carrier-first + exact finite-field PIT + native
   `SLAB_M` + full map-back (trial + **test-covector** + Jacobian) + retire the audit's cross-anchoring loop
   (script:1085-1107). Verify it is faithful to §5c and does not over-claim (it is a DESIGN, unbuilt).
5. **STATUS.md** top clause + the state memory
   `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md` (frontmatter
   `description:` + the body "NEXT") accurately reflect 1–4, with **no stale STEP-0 language** and no over-claim.

## Flag specifically
- Any **over-claim**: F/G stated as "resolved" (should be deferred/withdrawn); the N6 control stated as "passed/held"
  (should be "never run"); the CAS-authorship rule misstated; the built F/G script implied usable/reviewed.
- Any **factual error** vs the commits/records (SHAs, what was corrected, what is owed).
- Any **owed item omitted**: numeric F/G re-grounding; the N6 build; MEMORY.md compaction; the 499 MB `.out` durability.

## Output
End with **CLEAR TO COMPACT** (docs accurate + open items correctly recorded) or the exact fix list. Keep it brief.
