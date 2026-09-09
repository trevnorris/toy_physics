# Independent review — S11c-c2 N6 reconcile DISPOSITION record

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_N6_reconcile_disposition.md`

This is an orchestrator-written **governing disposition** of a dual-engine (SymPy + blind Wolfram) cross-engine
comparison. It records the decision to STOP a "collapse-instrument" series and carry an unresolved cross-engine
operand agreement as a DEBT (called "Path B"). Your job is to check that the disposition is **HONEST** — faithful
to its sources, free of over-claim and under-claim, and structurally correct — ⛔ NOT to re-run the physics build.

## What to check (report a finding only if it changes what the disposition may claim)
1. **Numerical fidelity.** Every count in the artifact (§1, §2, §8) must match the committed tally
   `_measurements/S11c_c2_N6_comparator_run_tally.txt` LITERALLY. Read the tally yourself and verify: the matched
   vanishings (`R_COV*` 160/0, `CARRIER_BRIDGE_RESIDUAL` 320/0, `SOURCE_CONTROL_DELTA` 160/0, support 400/0), and
   the surfaced-nonzero operand gap (`CARRIER_EULERIAN`/`CARRIER_MATERIAL` 40, `SOURCE_ACTUAL`/`BASELINE`/
   `PREDICTED` 76, `FROZEN_PHI` 18, census `ACTUAL_CONTROL_PARAMETERS` 4 + `PHI_DOMAIN_CENSUS` 4). Flag any
   mismatch.
2. **The `(0)−(0)` claim (§1).** The artifact says the matched cross-engine zeros are `(0)−(0)` — a dual-engine
   confirmation of a VANISHING statement, ⛔ NOT operand agreement — because each matched object (`R_cov`, carrier
   bridge `C_E−C_M`) is already 0 within each engine. Verify this is correct from the two engine sources
   (`scripts/S11c_c2_N6_{covariance,reconcile}_sympy.py`, `mathematica/S11c_c2_N6_mathematica_audit.wl`) and the
   question doc §2. Is there any matched key where the cross-engine 0 is genuinely operand agreement (not
   `(0)−(0)`) that the artifact wrongly lumps in? Is there any surfaced-nonzero operand the artifact wrongly
   claims agreed?
3. **The structural obstruction (§3) — the load-bearing claim.** The artifact says the collapse bridge is the
   WRONG OBJECT for TWO independent reasons: (a) the emitted source operands are `μ = EL(energy density)` and
   `EL(T·L) ≠ T(EL·L)` for a position-dependent `T = R_W = W_0/W_bg`, so a density-level coefficient table applied
   AFTER Euler–Lagrange differentiation misses product-rule terms; (b) `R_W·W_bg = W_0` is a product of two
   η-dependent series and both engines expand `W_bg→W_0(1+η·w1)` BEFORE grade extraction, so an ungraded identity
   cannot be reproduced per-grade without cross-grade convolution (forbidden by the no-grade-mixing contract).
   ⭐ **Verify reason (a) YOURSELF by direct reasoning** — EL non-commutation with multiplication by a varying
   field is a standard variational fact; check the astra E04/E14 counterexample in
   `_measurements/S11c_c2_N6_reconcile_disposition_PATH_B.md:14-22` reproduces (does a σ_W^1 term survive?). Check
   reason (b) against the cited engine lines (`.wl:127-141`, `S11c_c2_N6_reconcile_sympy.py:245-248`) and the v3
   review `_legs/S11c_c2_N6_reconcile_collapse_v3_review_grok.md`. Is either reason overstated, or is the
   conclusion ("upstream of EL; a replay would not reconcile the already-emitted `.out`") stronger than the
   reasons support?
4. **No over-claim.** Flag any place the artifact upgrades: a `(0)−(0)` vanishing → "operand AGREE"; the correct
   ungraded density table (§4) → "the sources agree"; "representational-difference-UNADJUDICATED" → "known to be
   just thickness" (§5 forbids this — is the guard actually honored in the artifact's own wording?); the debt
   dismissed as "c2 already has everything it needs" (§6 forbids it, citing S11c-d gradient-driven mixing —
   `S11c_decisions.md:83`). Also flag the reverse: any place it UNDER-claims (e.g. calling per-engine covariance
   "weak N6" when two engines confirmed it).
5. **Completeness + carry-open.** Are all 3 N6 premise caveats (Φ physical-correctness; V transform `V_E≡V_M`;
   block leakage) carried (§7), plus the debt, the un-inspected leftover SHAPE, and the earlier carries (2 S11c-b
   signs / 6 §3d / c1 ENERGY)? Is anything the disposition needs to say missing?
6. **Scope discipline.** Does the artifact anywhere re-litigate the N6 covariance resolution (Reading B,
   `R_cov=0`, `d21c8ff5`) or c1 (both must STAND), or pre-adjudicate a carried-open item that should only be
   surfaced?

## What you are handed (read these as your source of truth, BEFORE forming a verdict on the artifact)
- The tally + run record: `_measurements/S11c_c2_N6_comparator_run_{tally.txt,data.md}`.
- The two engines: `scripts/S11c_c2_N6_covariance_sympy.py`, `scripts/S11c_c2_N6_reconcile_sympy.py`,
  `mathematica/S11c_c2_N6_mathematica_audit.wl`; the comparator `scripts/S11c_c2_N6_cross_engine_comparator.py`.
- The framing + planning docs: `_measurements/S11c_c2_N6_reconcile_question.md`,
  `_measurements/S11c_c2_N6_reconcile_disposition_PATH_B.md`, `_measurements/S11c_c2_N6_RESOLVED.md`.
- The closure evidence: `_measurements/S11c_c2_N6_reconcile_collapse_directive_gate.md`;
  `_legs/S11c_c2_N6_reconcile_collapse_v3_review_{opus,grok}.md`;
  `_legs/S11c_c2_N6_reconcile_strategy_consult_{astra,grok}.md`.
- Everything in this repo is readable; ground every claim in a file+line, ⛔ never in paraphrase.

## Required method (DOCUMENT review)
Read the source-of-truth files FIRST (the tally, the two engines, the question doc, the consult reports), form
your OWN view of what the comparison established and what it did not, and ONLY THEN read the disposition artifact.
Quote BOTH sides for every finding (the source line + the artifact line). ⚠ This reading order is a method
request, ⛔ not a blindness control — you receive everything at once.

⭐ For the ONE load-bearing physics question (reason 3a, EL non-commutation), do the small symbolic calc yourself
and show it, ⛔ do not assert a conclusion in prose. A prose "I checked and it's right" is discarded.

## Physics filter
Report a finding only if it catches a way the DISPOSITION could mislead a future reader about what the N6
cross-engine comparison established — an over-claim, an under-claim, a wrong number, a misstated obstruction, or a
missing carry-open. ⛔ Do not report "the disposition would be wrong if the run had produced different data" or
style preferences.

## Output
End with a one-line verdict: **SOUND** (nothing outstanding changes what the disposition may claim) or
**NOT-SOUND** (list the must-fix findings). For each finding: the artifact line, the source line, and what must
change.
