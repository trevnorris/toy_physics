# Independent review — S11c-c2 self-energy fold STEP RECORD

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_c2_self_energy_fold.md`

This is an orchestrator-written **step record** — the interpretation layer for S11c-c2 (the self-energy fold that
closes the S11c-b slab operator with c1's curved-bulk response). Your job is to check that it faithfully represents
its sources and does NOT over- or under-claim what S11c-c2 established. This is a DOCUMENT review, ⛔ NOT a re-run
of any build.

## What to check (report a finding only if it changes what the record may CLAIM)
1. **Source fidelity.** Every claim about what was built/reviewed/resolved must match its cited source. Read the
   sources yourself and verify: the self-energy fold physics review (what is SOUND vs OPEN — the adjudication was
   itself CORRECTED, so the record must NOT propagate the `8f3a017f` commit subject's "0 defects"); the N6 per-engine
   resolution (Reading B; `R_N6=18/288` nonzero AND `R_cov` no-nonzero, BOTH preserved); the N6 cross-engine
   disposition (Path B, the operand DEBT); the §5c correction; F/G paused. Flag any claim not grounded in its source.
2. **Per-engine vs cross-engine honesty — the load-bearing distinction.** The record claims the self-energy fold is
   SymPy-only (no WL self-energy engine, no self-energy comparator), and that the ONLY cross-engine work on this box
   is the N6 representation-invariance thread. Verify this against the actual artifacts
   (`ls scripts/S11c_c2_*`, `ls mathematica/S11c_c2_*`). Is there any cross-engine self-energy check the record
   missed, or any cross-engine claim it makes that no engine/comparator supports?
3. **The `(0)−(0)` N6 point.** The record says the matched covariance-channel cross-engine zeros are `(0)−(0)` (a
   dual-engine confirmation of a VANISHING, ⛔ NOT operand agreement) while the surfaced operand residuals (carrier
   40 / source 76 / Φ 18) are UNADJUDICATED and carried as a DEBT. Verify against the committed tally
   `_measurements/S11c_c2_N6_comparator_run_tally.txt` and the disposition. Is the split stated correctly?
4. **The `I_{M→E}` fix.** The record says `I_{M→E} = extract(close(SLAB_M)−SLAB_M)` is the MATERIAL-anchoring
   increment (native material sources differenced directly, no `T`/pullback), NOT a "mapped-to-Eulerian operand", and
   that the frame-change faithfulness is checked separately by `R_cov`. Verify against `S11c_c2_N6_RESOLVED.md:54`,
   the SHARED_PHYSICS §5c definition, and `S11c_c2_N6_diagnostic_sympy.py`. Is the terminology fix correct?
5. **The forbidden over/under-claims.** Flag any occurrence of: "weak N6" (the record must NOT — two engines
   confirmed covariance); "c2 already has everything it needs" (forbidden — S11c-d uses gradient-driven mixing, so
   the debt is material); "known to be just thickness" (forbidden — the leftover SHAPE was not inspected); "0
   defects" (the corrected adjudication forbids it); any surfaced carry-open item that is PRE-ADJUDICATED rather than
   surfaced (the 2 S11c-b signs, the 6 §3d re-adjudications, c1 ENERGY must be surfaced, ⛔ not decided).
6. **Completeness + carries.** Are all carries present and correctly attributed: the cross-engine DEBT + un-inspected
   SHAPE; the 3 N6 premise caveats; the 2 S11c-b sign conventions (which do NOT cancel from c2's residual — verify
   the kinetic-vs-face-force distinction); the 6 §3d re-adjudications; c1 ENERGY; the F-wording OWED (§5e still says
   "must vanish"); F/G PAUSED; the census 6-row crosswalk? Is anything the record needs to say missing, or any SHA
   wrong (`git show --stat` the ones you doubt)?

## What you are handed (read these as your source of truth, BEFORE forming a verdict)
- The N6 governing disposition `_measurements/S11c_c2_N6_reconcile_disposition.md`; the tally
  `_measurements/S11c_c2_N6_comparator_run_tally.txt` + run record `..._comparator_run_data.md`.
- The N6 per-engine resolution `_measurements/S11c_c2_N6_RESOLVED.md`; the self-energy physics review adjudication
  `_measurements/S11c_c2_physics_review_adjudication.md`.
- The physics authority `directives/S11c_c2_SHARED_PHYSICS.md` (esp. §0–1, §3d, §5).
- The predecessor step records `steps/S11c_c1_curved_bulk_closure.md`,
  `steps/S11c_b_variable_coefficient_operator.md` (for the 2 S11c-b signs + the c1 UNDECIDED imports).
- The engines/comparator: `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`,
  `scripts/S11c_c2_N6_{covariance,reconcile,diagnostic}_sympy.py`,
  `mathematica/S11c_c2_N6_mathematica_audit.wl`, `scripts/S11c_c2_N6_cross_engine_comparator.py`.
- Everything in this repo is readable; ground every finding in a file+line, ⛔ never in paraphrase.
- ⛔ You are NOT handed the build directives (a record can satisfy a directive and still misrepresent its source);
  form your own view from the sources above.

## Required method (DOCUMENT review)
Read the source-of-truth records FIRST, form your OWN view of what S11c-c2 established and what it left open, and
ONLY THEN read the step record. Quote BOTH sides for every finding (the source line + the record line). ⚠ This
reading order is a method request, ⛔ not a blindness control — you receive everything at once. A measured claim that
cites a run must carry a real command/record; a prose "I checked" is discarded.

## Physics filter
Report a finding only if it catches a way the record could mislead a future reader (an over-claim, an under-claim, a
wrong SHA/number, a pre-adjudicated carry, a missing open item), ⛔ not style preferences or "would be wrong on
different data."

## Output
End with a one-line verdict: **SOUND** (nothing outstanding changes what the record may claim) or **NOT-SOUND**
(list the must-fix findings). For each finding: the record line, the source line, and what must change.
