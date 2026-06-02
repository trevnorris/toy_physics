---
unit_id: 218
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 5
findings_total: 5
material_change: false
---

# Verification — unit 218

Checkpoint stage. Higher bar applied; no rubber-stamp. All five findings assessed against the current
script files, the diff, and the captured exec logs (scripts not re-run by the verifier).

## Per-finding outcomes

### F1 — mathematica_transliteration (CRITICAL)

**Classification:** resolved

**What changed:**
The `.wl` (`mathematica/...mathematica_audit.wl`) was fully re-authored. The diff deletes every
transliterated helper (`packetInterval`, `boundaryBest`, `classifyFamily`, the `Association`-based
`quadrupleFaces` with `"omit_"` packets, the `Min @@ Append[...]` flattening block, and the
byte-identical `improvementFamily`/`noImprovementFamily`/`overlapFamily` numbers) and replaces them
with native Mathematica derivations.

**Assessment — a reviewer CANNOT line-pair `.wl` to `.py`. Three corresponding sections:**

1. **M1 incidence.**
   - `.py:181-190`: imperative Python loop `for subset in proper_faces: covering_faces =
     [name for name,packet in quadruple_faces.items() if set(subset).issubset(packet["support"])]
     ... if len(covering_faces) != expected: raise`.
   - `.wl:42-64`: native set combinatorics — `quadrupleFaces = Subsets[axes,{4}]`,
     `properStrata = Subsets[axes,{1,4}]`, covering count via
     `Total[Boole[ContainsAll[#,support]] & /@ quadrupleFaces]`, strata count via
     `Tally[Length /@ properStrata]`, and `AllTrue[incidenceRows, #[[3]]==#[[4]]&]`.
     These are exactly the directive's named primitives (`Subsets`/`ContainsAll`/`Boole`+`Total`/`Tally`);
     there is no Python-loop counterpart. `grep` confirms 15 uses of native primitives and zero
     `classifyFamily`/`packetInterval`/`boundaryBest`.

2. **M2/M3 ceiling + splice.**
   - `.py:87-133` (`certify_splice_bracket`): per-packet **boolean contradiction** proof —
     `sp.simplify_logic(sp.And(lo > best, Le(lo,best))) is sp.false` (lower) and the `best > hi` analog
     (upper).
   - `.wl:83-129`: **quantifier elimination** — `Resolve[ForAll[{...},Implies[Element[...,Reals] &&
     lo<=best, !(lo>best)]], Reals]` per packet, plus the two ceiling branches
     `Resolve[ForAll[..., Implies[..., Min[a,b]==a]]]`. Different engine route (real-domain `Resolve`
     vs. SymPy `simplify_logic`). The forbidden `Min[Min[...]]==Min[...]` flattening identity is gone
     from both files (`grep` finds no `Min[Min` and no `tau_le5_*_flat` / `expect_equal(...flat...)`).

3. **M4 regimes — non-byte-identical witnesses.**
   - `.py:247-273`: still the original hardcoded dicts `improvement_family={"support_le3":(10,11),
     "omit_lambda":(12,13),...}`, `support5_int:(2,3,4)` — improvement total 192.
   - `.wl:144-182`: witnesses **regenerated from the regime hypotheses** via
     `makeBoundaryWindows[start,width,gap]` — improvement `makeBoundaryWindows[20,2,3]` +
     `support5_int->Range[1,4]`; no-improvement `makeBoundaryWindows[2,2,3]` + `Range[40,42]`; overlap
     `ConstantArray[{4,8},5]` + `{3,7}`. None of these are the `.py` numbers, and the resulting outcome
     counts differ (improvement 256 / overlap 65+63 in `.wl` vs 192 / 64+64 in `.py`), proving an
     independent witness set rather than a transcription.

   M4 still asserts the F5 exhaustive equalities (`support5 wins equal total`, `boundary wins equal
   total`, `support5+boundary == total`, `ties==0`) and gates the regime hypotheses (5.1
   `Max[support5] < Min[boundary lo]`, etc.). Per-claim PASS/FAIL with `Exit[1]` on failure preserved.

   `wl_now_independent: yes.`

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Both the `.py:183-185` `Min`-flatten `expect_equal` and the in-loop construction guard (`.py` old
`:70-71` `full_lo <= full_value <= full_hi` inside `classify_family`) were removed. `classify_family`
now only tallies wins (`.py:56-76`). The splice is proved by the new `certify_splice_bracket`
(`.py:87-133`) and its `.wl` `Resolve` mirror (`.wl:104-140`).

**Assessment — falsifiable in BOTH engines under the named mutations:**
- The per-packet proofs show `And(lo>best, lo<=best) ≡ false` and `And(best>hi, best<=hi) ≡ false` —
  genuine contradictions, not auto-flattening; the log prints each branch (`lower/upper splice
  counterexample branch <name> = False`, then `... contradicted`).
- The concrete endpoint probes make the `max(hi)` mutation observable:
  `upper_probe = (3,6,7,8,9,10,11)`; the check asserts `certified_min_endpoint(upper_probe)==3` AND
  `certified_min_endpoint(upper_probe) < max(...)=11` (log: "upper splice probe rejects max endpoint =
  True"). Mutating the upper splice to `max_p(hi_p)` would yield 11, failing
  `... probe equals least packet hi = 0` and the "rejects max endpoint" gate. A lo/hi swap makes the
  lower probe equal 3 (not its required 1), failing `lower splice probe equals least packet lo`. The
  `.wl` mirror (`M3 upper endpoint audit probe`, `M3 upper endpoint is sharper than max-hi mutation`)
  fails identically. `f2_falsifiable: yes.`

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
`.py:347-366`: `support_le4_budget` is now the literal paper value `1140` (not `600 + 5*2*54`); the
`600`/`54` decomposition survives only as an inline comment. Gating assertions compare directly to the
paper constants: `support<=4 == 1140`, `support-five lifted == 324`, `fallback == 1500`,
`preferred == 1464`, `fallback total == 2640`, with `lifted_per_envelope == 162` and
`projected == 750` retained as factorization checks. The `.wl:217-239` mirror does the same
(`supportLe4Budget = 1140` literal, then `M5 ... - 1140/324/1500/162/750/1464/2640` all `== 0`).

**Assessment:** A `600/54` transcription error that still summed to 1140 can no longer slip through —
`1140` is the gate, not a constructed sum. No value changed (1140/324/1464/2640/162/750 all preserved,
log confirms). Correct.

### F4 — paper_misalignment (notes 230→162; USER-RESOLVED direction (a))

**Classification:** resolved

**What changed (script side, which is all the verifier checks):** nothing. The SymPy `.py:344`
`expect_zero("lifted compiler bound - 162", lifted_per_envelope - 162)` is unchanged — diff line 676
shows it as a context line (leading space), confirming Codex did not touch the script's `162`. The
notes `230→162` edit is out of the verifier's scripts-only scope and was already reviewed by the
orchestrator. `py_162_unchanged: yes.`

### F5 — insufficient_verification

**Classification:** resolved

**What changed:**
`.py:278-326` and `.wl:184-213`: each regime now asserts the **exhaustive** count plus a regime-
hypothesis gate, replacing the old `> 0` "wins exist" assertions.
- Improvement (5.1): hypothesis `max(support5_int) < min over boundary lows`; then `boundary==0`,
  `tie==0`, `sum==total`, `support5==total` (log: 192/0/0, all four pass).
- No-improvement (5.2): hypothesis `min(support5_int) > max over boundary highs`; then `support5==0`,
  `tie==0`, `boundary==total` (log: 0/192/0).
- Overlap (5.3): interleave hypothesis; then `support5>0`, `boundary>0`, `tie==0`,
  `support5+boundary==total` (log: 64/64/0).

**Assessment — counts are exhaustive and non-trivial.** The regime-hypothesis gate is asserted
*separately* and passes (`... regime hypothesis = True`), so the equalities are not vacuous: the
windows genuinely sit in their regime, and a misclassifier or a wrong window would break the exact
equality (e.g. one stray boundary win would make `support5 != total`). The `.wl` carries the same
exhaustive equalities on its independent witness set (256/0/0, 0/192/0, 65/63/0).
`f5_exhaustive_counts: yes.`

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `#proper nonempty support strata = 30`; per-subset incidence printed for all 30 strata.
- `lower/upper splice branch <name> contradicted = True` for all 7 packets (×2); `support<=5 upper
  splice probe rejects max endpoint = True`.
- `improvement {support5:192,boundary:0,tie:0}`, `no-improvement {0,192,0}`, `overlap {64,64,0}` with
  the exhaustive equalities at `= 0`.
- Budget: `support<=4 paper budget - 1140 = 0`, `... - 324/1500/1464/2640 = 0`.

**Mathematica:** exit=0, 31 PASS, 0 FAIL. Notable lines:
- `M1 every proper support has 5-k covering quadruple faces`, `M1 proper strata tally is 5+10+10+5`,
  `M1 ... - 30`, `M1 2^5 - 2 - 30` all PASS.
- `M2 tau<=5 best is the two-ledger support ceiling`, `M3 lower/upper splice counterexample branches
  close`, `M3 upper endpoint is sharper than max-hi mutation` PASS.
- `M4 improvement outcomes = {support5->256, boundary->0, tie->0}` (independent witness set, distinct
  from `.py`'s 192), regime 5.1/5.2/5.3 hypothesis + exhaustive-count PASS.
- `M5 ... - 162/324/750/1500/1140/1464/2640` all PASS.

PASS lines: 31 (Mathematica). They cover M1–M5 plus the F2 contradiction/endpoint substance and the
F5 exhaustive counts. Nothing was weakened to go green — assertions were strengthened (equalities and
contradiction proofs replacing `>0` and `Min`-flattening).

**Output freshness:** confirmed. `.py` mtime 12:27:42, `.wl` mtime 12:30:19; both saved `.txt`
outputs mtime 12:33:13 — newer than the scripts, regenerated post-fix.

## Material-change assessment

`material_change`: false. The edits are a de-transliteration of the `.wl`, stronger/falsifiable splice
and regime checks, a direct paper-anchored budget block, and a notes typo fix. No derived value a
downstream unit consumes changed (162/324/750/1500/1140/1464/2640 and the `min`/`min` splice form all
preserved; logs reproduce identical SymPy numerics).

## Side observations (non-blocking)

- The `.py` banner still reads "STAGE 201 …" (cosmetic, pre-existing; the `.wl` banner correctly says
  "Stage 218"). Not part of any finding; not blocking.
- `.py` `expect_equal` (`:31`) and `packet_interval`/`full_best` (`:41,:52`) are now dead helpers
  (defined, unused). Harmless leftover; not blocking.
- The budget comment at `.py:347` references "Stage 249's upstream decomposition" while elsewhere the
  card cites Stage 215 for the 600/54 split — a comment-only provenance label, no effect on any gate.

## Verdict justification

All five findings are `resolved`. The critical F1 is genuinely closed: the re-authored `.wl`
establishes M1 by `Subsets`/`ContainsAll`/`Boole`/`Tally`, M2/M3 by `Resolve`/`ForAll` quantifier
elimination, and M4 by hypothesis-generated witnesses whose counts differ from the `.py` — a reviewer
cannot line-pair the two engines. F2 replaces the `Min`-flatten and construction guard with per-packet
contradiction proofs plus endpoint probes that fail under the `max(hi)` mutation or a lo/hi swap, in
both engines. F3 gates directly on the paper constants. F4's script `162` is untouched (notes-only,
already reviewed). F5 pins each regime to an exhaustive count with a separately-asserted, satisfied
regime hypothesis. Both engines exit 0 (31 Mathematica PASS, 0 FAIL either side), outputs are fresh,
and nothing was weakened to pass. Checkpoint bar met; `material_change: false`.
