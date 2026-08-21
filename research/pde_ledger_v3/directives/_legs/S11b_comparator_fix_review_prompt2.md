# Independent review — S11b comparator AFTER fix round 1 (Codex-written instrument)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py`
(just repaired in place; ~1357 lines). It classifies each shared `PY_S11B_<QUANTITY>` / `WL_S11B_<QUANTITY>`
object AGREE / DISAGREE / UNDECIDED / UNCOMPARED. Its ONE cardinal sin: reporting a genuine
content/sign/dimension/structure difference as anything other than DISAGREE.

## What to check
The fix closed three holes and added a precedence rule + per-leaf budget + streaming. Verify (a) the three
holes are truly closed, (b) NO NEW hole was opened by the ~274 new lines, and — the adversarial core — (c)
can ANY code path still make a genuine disagreement come out as AGREE / UNCOMPARED / UNDECIDED?

## What you are handed
- The artifact above.
- The fix directive it must satisfy: `research/pde_ledger_v3/directives/S11b_comparator_fix_round1_directive.md`
  (the unifying rule: classify each LEAF, record class = most severe leaf under
  `DISAGREE > UNCOMPARED > UNDECIDED > AGREE`; the three fixes; per-leaf budget; streaming). Read it as the
  source of truth, then read the code.
- The two real transcripts, to check no crash / no hidden disagreement:
  `scripts/out/S11b_interface_coupling_law_sympy_audit.out`,
  `mathematica/out/S11b_interface_coupling_law_mathematica_audit.out`. ⚠ Correctness verdict comes from your
  ablations/probes, ⛔ not from whether real output "looks right."

## Required method — SCRIPT; ablation mandatory
Construct your own PY/WL fixtures (placeholder symbols/values); run the comparator; save each fixture + its
literal stdout to named absolute paths and report them. A prose claim is discarded.

⛔⛔ Probe the precedence rule hardest — it is the new load-bearing core:
- A record with a DISAGREEing leaf AND (a status token / a boolean leaf / a budget-exceeded heavy leaf /
  a collision) — must be DISAGREE every time. Try each combination and nesting depth.
- `factor`/`cancel` on a differing leaf that returns a nonzero residual the code might mis-read as zero.
- A per-leaf budget: does a timed-out leaf ever land AGREE or bury a cheap DISAGREEing sibling? Force a heavy
  leaf beside a cheap differing sibling.
- Collision: function↔function AND function↔symbol targets → must flag TRANSLITERATION_COLLISION (UNCOMPARED),
  never AGREE. Also verify a NON-colliding rename still residuals (no false collision).
- Tuple promotion: Symbol-keyed and mixed-key tuples → STRUCTURE DISAGREE; all-`Str` → AGREE (not broken).
- UNDECIDED only when the severest leaf is a lone status/coverage token; never above a real difference.
Also FORM-ablate the precedence ordering (swap two levels) and confirm a fixture catches it. And check the
build's extended acceptance is decisive (name any test a bad fix survives). Confirm NO regression: the 7 core
rules and the repoint ablation still behave.

## Physics filter
Report a finding only if it is a way the comparator could be WRONG — above all, report agreement (or a
non-DISAGREE bucket) on a genuine disagreement, or a crash/hang on real input. A null result needs shown
ablations.

## Ablation sandbox
Copy to /tmp and ablate the COPY. ⛔ Never modify the working tree.

## Report — under ~30 lines
Numbered findings (code line, fixture/stdout path, wrong classification, fix); end with a one-line verdict:
sound to FREEZE and run, or repair again?
