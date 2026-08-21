# Independent review — S11b comparator AFTER fix round 2 (Codex-written instrument)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_cross_engine_comparator.py` (~1413 lines).
It classifies each shared `PY_S11B_<QUANTITY>` / `WL_S11B_<QUANTITY>` object AGREE / DISAGREE / UNDECIDED /
UNCOMPARED. ONE cardinal sin: a genuine content/sign/dimension/structure difference reported as anything
other than DISAGREE (including buried as UNCOMPARED or UNDECIDED).

## Fix round 2 closed two defects — verify them AND that no new hole opened
- **Defect 1** (was: over-broad collision). The collision guard must flag a collision ONLY on an actual
  same-kind post-transliteration merge: group by (kind, target); a collision needs 2+ distinct SOURCE NAMES
  in one (kind,target) group. So `f_a`/`fA` SAME KIND → TRANSLITERATION_COLLISION/UNCOMPARED; a function
  `fA` + a symbol `fA` (same name, different kind) → NOT a collision; a cross-kind distinct pair
  (`f_a` function + `fA` symbol) → NOT a collision (residual normally). ⛔ Probe: does any genuinely-different
  record get buried as UNCOMPARED by a FALSE collision? Does any REAL same-kind collision escape?
- **Defect 2** (was: render hang). Classification (NAME + CLASSIFICATION + reason) must be emitted and
  FLUSHED before the operands are rendered, and huge operands are rendered in bounded form. ⛔ Probe: can a
  large operand ever delay or alter the classification? Can the bound ever drop/alter a classification (it
  may only truncate the DISPLAYED operand)?

## Sources of truth / method
- Fix directive: `research/pde_ledger_v3/directives/S11b_comparator_fix_round2_directive.md`; the round-1
  contract: `S11b_comparator_fix_round1_directive.md` (the precedence rule
  `DISAGREE > UNCOMPARED > UNDECIDED > AGREE`).
- SCRIPT review: construct your own PY/WL fixtures (⚠ WL uses `f[x]` for application, `f(x)` is `f*x`), run
  the comparator, FORM-ablate on /tmp copies, save fixture+stdout to named paths and report them.
- ⛔⛔ Adversarial core, probe hardest: any combination where a genuine disagreement comes out as AGREE /
  UNCOMPARED / UNDECIDED. Re-confirm the three round-1 holes stay closed and the 7 core rules + repoint
  ablation + precedence do not regress. Check the strengthened acceptance actually RE-RUNS the tool and
  FAILS when a fix is reverted (name any test a bad fix survives).
- ⚠ KNOWN, not to be reported as the cardinal sin: the real pair is large (MB-scale objects); parsing is
  slow (a full run is many minutes; background it) and some WL forms (`Unequal[...]`, chained bracketed
  `Inequality[...]`, one un-parseable AST) land UNCOMPARED. UNCOMPARED-for-a-parse-gap is SAFE (never a
  misclassification). Report a parse gap ONLY if it produces a WRONG classification (a false AGREE), not for
  slowness or an honest UNCOMPARED.

## Report — under ~25 lines
Numbered findings (code line, fixture/stdout path, wrong classification, fix); end with a one-line verdict:
sound to FREEZE and run, or repair again?
