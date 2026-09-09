# QUESTION-vet: S11c-c2 N6 cross-engine reconcile framing

## Your role
You are one of two INDEPENDENT legs vetting a **question**, not answering it. The artifact is an
orchestrator-written, physics-bearing **reconcile-framing** document. Your job: decide whether it asks the
RIGHT question, at the RIGHT (retained) order, without pre-adjudicating the answer or smuggling in a fold. You
are NOT asked to run the collapse test (that instrument does not exist yet and will be Codex-written + reviewed).

## Artifact under vet
`/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_question.md`

## Sources of truth (read these FIRST; form your own view before judging the framing)
- The two engines whose cross-engine residuals are the subject:
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_covariance_sympy.py` (SOURCE ACTUAL/PREDICTED/BASELINE, Φ, R_cov)
  - `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py` (CARRIER Eulerian/Material, bridges, R_N6)
  - `research/pde_ledger_v3/mathematica/S11c_c2_N6_mathematica_audit.wl` (the blind WL side of every residual)
- The comparator that produced the residuals: `research/pde_ledger_v3/scripts/S11c_c2_N6_cross_engine_comparator.py`
- The run evidence (committed): `research/pde_ledger_v3/_measurements/S11c_c2_N6_comparator_run_data.md` and
  `…run_tally.txt`.
- The per-engine N6 resolution the disposition must stay consistent with:
  `research/pde_ledger_v3/_measurements/S11c_c2_N6_RESOLVED.md`.
- The comparator directive review (shows what was and was NOT pre-registered — rule 5):
  `research/pde_ledger_v3/_measurements/S11c_c2_N6_comparator_directive_review.md`.
- The c1 reconcile precedent (the staged representational bridge, and its UNDECIDED-surfacing pattern):
  `research/pde_ledger_v3/_measurements/S11c_c1_comparator_reconcile.md`.

## What to check — report a finding only if it changes what is asked or what may be claimed
1. **Right question, not a proxy.** Does the framing target the ACTUAL open cross-engine gap? In particular §2's
   claim: within each engine `C_E−C_M=0` and `R_cov=ACTUAL−PREDICTED=0`, so those differences are 0 cross-engine
   TRIVIALLY, and that says nothing about whether the OPERANDS (carrier, source, Φ) agree cross-engine. Verify
   this structural argument against the actual engine code (cite lines). If it is wrong, that is a MUST-level
   finding — the whole disposition rests on it.
2. **The two-errors-cancel risk in channel (b).** §3(b) asks for a SINGLE identity that simultaneously collapses
   ACTUAL/PREDICTED/BASELINE, else `R_cov=0` may be a cancellation of two genuinely-different sources. Is that
   the correct discriminator? Is there a weaker or stronger test that is the RIGHT one?
3. **Retained order (§5).** Confirm the collapse test must run at the retained η/σ_W rectangle and that a
   `σ_W→0` (or any other) reduction would be a proxy that repeats the measured L-CAS over-clear. Check the engine
   actually retains that order (cite the normalization / a retained-order object).
4. **Candidate identities (§4).** Are these the right physics identities to TEST for collapse? Any missing (an
   identity the residuals plainly need), or any that would make a collapse VACUOUS / trivially true, or any that
   is actually a rule-5 fold that must NOT be pre-registered? Ground each in the engine constructions.
5. **No pre-adjudication / no leak (§6).** Does the document anywhere pre-decide that a channel reconciles or
   genuinely differs, or leak an expected outcome value? (It must not.)
6. **Scope fork (§8).** Given the objects' measured tractability (~330 MB peak, ~80 MB largest — from
   `run_data.md`), is a FAITHFUL retained-order collapse test actually available on a normal box (path A), or is
   there a reason the retained-order test is itself heavy/must be deferred and the channels SURFACED as UNDECIDED
   (path B)? Give your reasoned recommendation with the evidence.
7. **Census/minor (§3d).** Are `ACTUAL_CONTROL_PARAMETERS` (4) and `PHI_DOMAIN_CENSUS` (4) correctly treated as
   provisional bookkeeping to confirm, not assumed non-load-bearing? Check what those tags actually carry
   (`cov:105-124,146-153`).
8. **Architecture (§7).** Is the disposition pipeline rule-compliant (Codex-written+G1-reviewed instrument; the
   disposition adjudication gets legs; mechanical-fact-lookup vs instrument boundary correct)?

## Method / evidence discipline
- Read the source of truth first, form your own view, THEN judge the framing. Quote both sides for every finding.
- ⛔ A prose assertion is worth nothing here: ground every claim in a cited file+line or the committed tally.
  Where a claim is about what the code computes, name the line that computes it.
- ⛔ You are NOT asked to determine whether any residual is representational vs genuine — you cannot without the
  instrument. Assess only whether the framing correctly LEAVES that open and asks the right question.
- If you believe the question is right as framed, say so explicitly and name the one change (if any) that would
  most strengthen it. If you would re-frame it, give the exact re-framing.
- End with a one-line verdict: **QUESTION SOUND** or **QUESTION NEEDS-WORK**, and for the scope fork, **A** or
  **B** with one sentence why.
