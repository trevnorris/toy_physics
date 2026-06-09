---
unit_id: 218
batch: VI.1
verifier_model: Opus 4.8 (1M context) [claude-opus-4-8[1m]]
verify_date: 2026-06-09T17:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 218 (CHECKPOINT)

The directive carries exactly one finding (F1 = `stale_output`) and it explicitly requires NO Codex
edit — the orchestrator's standard re-run regenerates the `.txt`. I confirmed (A) the output refresh
landed clean, (B) the checkpoint higher bar holds on the refreshed artifacts (both engines
substantive, independent, paper-exact), (C) the budget tally is internally consistent and matches the
appendix, and (D) `material_change: false`. Scripts-only, no execution; I read the scripts, the
committed outputs, the exec logs, and the part-06 appendix budget rows.

## Per-finding outcomes

### F1 — stale_output (no Codex action)

**Classification:** resolved

**What changed:**
No script edit (correct — the directive forbids one; there is no `stage_218_diff.patch`, confirming
Codex touched nothing). The fix is the orchestrator's re-run. The committed SymPy output
(`scripts/output/…sympy_audit.txt`) line 3 now reads
`STAGE 218 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE`; the F1-stale
`STAGE 201` banner is GONE. The `.py` banner string (py:136) and the refreshed `.txt` now agree.

**Assessment:**
Correct and complete. Every numeric/logical line in the refreshed `.txt` is unchanged from the stale
capture (the auditor predicted this: the numbering-reconciliation commit only flipped the banner
`201`→`218`): strata `30`, all `= 0` residuals, all `True`/`contradicted=True` flags, regime counts
`192/192/64-64`, budget `1140/324/1500/1464/2640`. Pure label refresh, no content move.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `stage_218_sympy.log`:
- L8: `STAGE 218 — FULL SUPPORT-<=5 COMPLETION AND LOCAL MIXED-RAY SEARCH CLOSURE` (banner corrected)
- L17: `#proper nonempty support strata = 30`; L18-21: `axes - 5 = 0`, `quadruple packet count - 5 = 0`
- L76-103: all 14 splice branches `counterexample = False` / `contradicted = True`
- L150: `All Stage 218 imported ledger, splice, classification, and budget checks verified.`; exit 0

**Mathematica:** exit=0. Notable lines from `stage_218_mathematica.log`:
- L16: `support-cardinality tally = {{1, 5}, {2, 10}, {3, 10}, {4, 5}}`; L24: `M1 2^5 - 2 - 30 = 0`
- L33: `PASS: M2 tau<=5 best is the two-ledger support ceiling`; L35/L37 M3 lower/upper branches close
- L48-50: independent witnesses `improvement 256 / no-improvement 192 / overlap 65-63` (differ from
  `.py`'s `192/192/64-64` by design — independent window generators)
- L86-99: `M5 lifted per-envelope bound - 162 = 0` … `M5 fallback total - 2640 = 0`, all PASS
- L101: `All Stage 218 Mathematica claims M1-M5 verified.`; exit 0

**Output freshness:** confirmed. Both `.txt` mtimes are 2026-06-09 16:51:54, newer than the `.py`
(2026-06-03 15:59:11) and the `.wl` (2026-06-02 12:30:19). The committed SymPy `.txt` banner is
`STAGE 218`; F1 cleared.

## Checkpoint higher-bar assessment

**(A) Output refresh landed clean:** YES. SymPy committed output now reads STAGE 218 (was STAGE 201);
Mathematica was already fresh; both `.txt` mtimes postdate their scripts; arbiter/exec logs both exit 0.

**(B) Both engines substantive, independent, paper-exact on the refreshed artifacts:** YES.
- M1 boundary: `.wl` `Subsets[axes,{4}]` + `Boole[ContainsAll]` + `Sort[Tally[...]]` + `2^5-2-30`
  (wl:42-62) vs `.py` `itertools.combinations` + `set.issubset` (py:159-185). Different generators
  and predicates → independent.
- M2-M3 splice (load-bearing heart): `.wl` proves the Min ceiling and the `lo<=best<=hi` bracket by
  `Resolve[ForAll[{...}, Implies[Element[...,Reals] && ..., ...]], Reals]` real quantifier elimination
  (wl:83-126) vs `.py` `sp.simplify_logic(sp.And(lo>best, lo<=best)) is sp.false` finite propositional
  contradiction (py:111-124). Real-QE vs finite propositional → the named independence signature.
- M4 regimes: `.wl` `makeBoundaryWindows[start,width,gap]` generated witnesses (wl:144-176) give
  `256/192/65-63`; `.py` hand-listed windows give `192/192/64-64`. Same theorem, genuinely different
  witness data → independent (the differing counts prove the witnesses were not shared).
- M5 budget: shared product arithmetic in both engines (benign, posited bookkeeping anchored to the
  appendix; not the stage's theorem).
This matches the pass-2 orchestrator's ground-truth read; the re-author is genuinely independent
(unlike the V.3-200 surviving-transliteration case).

**(C) Budget internally consistent + appendix-exact:** YES. Both engines COMPUTE (not hardcode) the
products and sums and assert `prod - value == 0`:
`prod(3,3,3,3,2)=162` → `2*162=324`; `prod(5,5,5,6)=750` → `2*750=1500`; `1140+324=1464`;
`1140+1500=2640` (py:336-366, wl:217-239). Appendix part-06 matches exactly: `3^4·2=162` (L1200),
`2×162=324` (L1205), `750`/envelope → fallback `1500` (L1207), `600+540=1140` (L1079),
`1140+324=1464` (L1261), `1140+1500=2640` (L1266); the 217/218 summary rows (L65/L67) carry `162`
and `1464`. Internally consistent and paper-exact.

## Material-change assessment

`material_change`: false. F1 is a label-only banner refresh with NO numeric/logical change; no
script edit was applied (no diff). Every emitted deliverable value is byte-identical to the prior
capture and every value reconciles with the appendix. No downstream unit's carried value moves.

## Side observations (non-blocking)

- The dead numbering labels (`source_stage: 200` at py:330, the "Stage 249" comment py:347-348/the
  upstream `600+5*2*54` note, and `source_stage: 198/200` metadata) are never read or asserted on —
  they are the known incomplete-renumber artifact and are correctly deferred to the dedicated
  numbering pass. They do not touch any verified quantity. Non-blocking.
- The two engines' M4 regime counts differ by design (`256/192/65-63` vs `192/192/64-64`); this is
  the intended independent-witness shape, not an `engine_disagreement` — the regime CONCLUSIONS
  (5.1 all support5 / 5.2 all boundary / 5.3 both nonzero, no ties) agree in both.

## Verdict justification

The single finding (F1 `stale_output`) is resolved by the orchestrator re-run: the committed SymPy
output banner is now `STAGE 218` (stale `STAGE 201` gone), both `.txt` mtimes postdate their scripts,
and both engines exit 0 with all PASS lines and no FAIL. No Codex edit was applied (correct — none was
required; no diff patch exists), so there is no regression surface. The checkpoint higher bar holds on
the refreshed artifacts: both engines are substantive and exercise each of the three load-bearing
objects by genuinely different operations (`Subsets`/`ContainsAll`/`Tally` vs `itertools`/`issubset`;
`Resolve[ForAll[...],Reals]` real-QE vs `simplify_logic` finite contradiction; generator-windows
`256/192/65-63` vs hand-windows `192/192/64-64`), the lone shared object (the budget product) is a
posited bookkeeping tally that is internally consistent and appendix-exact (`162/324/750/1500/1140/
1464/2640`), and `material_change` is false. Unlike the V.3-200 case, this re-author is genuinely
independent. Verdict: verified, checkpoint bar CLEARED.
