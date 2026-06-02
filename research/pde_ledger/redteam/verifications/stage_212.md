---
unit_id: 212
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T12:15:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 212

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the new second-engine audit
`mathematica/moving_throat_pde_stage212_full_primitive_triple_ranking_theorem_mathematica_audit.wl`
(206 lines) covering claim manifest M1–M4, each with a fail-guard that calls
`Exit[1]` on mismatch and prints a `PASS:` line on success.

**Assessment:**
The `.wl` is a GENUINELY INDEPENDENT derivation versus the SymPy
`.py` (`scripts/.../sympy_audit.py`), not a transliteration. Corresponding sections:

- **M1 combinatorial ledger.** `.py` (lines 43–44, 55–65) uses
  `itertools.combinations`, Python `set(pair).issubset(tri)`, and generator
  `sum(...)` incidence loops. The `.wl` (lines 42–64) uses entirely different
  native primitives: `Subsets[axes,{2}]`/`Subsets[axes,{3}]`,
  `Count[triples, t_ /; SubsetQ[t, pair]]`,
  `Count[triples, t_ /; MemberQ[t, axis]]`, with `Length` compared against
  `Binomial[5,2]`/`Binomial[5,3]`. Non-tautological (target counts are
  independently recomputed via `Binomial`, not hardcoded).

- **M3 interval/order theorems (the load-bearing independence test).** `.py`
  (lines 98–185) proves these by FINITE INTEGER-RANGE ENUMERATION — nested
  `for ... in range(...)` loops over 3136/924/462/896/3136/924/462 integer
  samples. The `.wl` (lines 83–186) instead proves each as a universally
  quantified real inequality via SYMBOLIC QUANTIFIER ELIMINATION:
  `Resolve[ForAll[{vars}, hypotheses \[Implies] conclusion], Reals]` for M3a
  (`Min[betaLo,iotaLo] <= Min[bBest,iBest] <= Min[betaHi,iotaHi]`), M3b
  interior (`iotaHi<betaLo \[Implies] iBest<bBest`) and boundary
  (`iotaLo>betaHi \[Implies] bBest<iBest`), M3c ranking
  (`U1<L2 \[Implies] x<y`), and M3d sandwich/improvement/no-improvement —
  each guarded `If[res =!= True, fail; Exit[1]]`. This is the different,
  stronger decomposition the directive mandated; NOT the `.py` integer loop.

- **M4 finite evaluation budget.** `.wl` lines 190–200 use
  `pairwiseTotal = 10*12` (=120), `tripleInteriorTotal = 10*48` (=480),
  `fullTotal = 120+480` (=600), each checked with `expectExact` against literals
  120/480/600. The pairwise total is **120, not the stale notes value 188** —
  confirmed in both the script and the exec log (`M4 pairwise total = 120`,
  `PASS: M4 10*12`).

- **M2** is a parallel native-Min-flattening check in both engines
  (`Min[iota,Min[a,b,c]] === Min[iota,a,b,c]`, log lines 31/33 both `True`),
  which is the appropriate independent confirmation of the splice flattening.

All seven M3 `Resolve` calls return `True`; none is tautological (each carries
genuine interval-membership + ordering hypotheses driving a strict/non-strict
conclusion that `Reals` quantifier elimination actually decides).

### F2 — paper_misalignment (script-side banner label fix)

**Classification:** resolved

**What changed:**
`scripts/.../sympy_audit.py:35` banner changed
`STAGE 195 — FULL PRIMITIVE-TRIPLE RANKING ...` →
`STAGE 212 — FULL PRIMITIVE-TRIPLE RANKING ...`.

**Assessment:**
The captured diff (`stage_212_diff.patch`) shows exactly ONE changed line — the
banner — with no other `.py` line touched. The SymPy script is otherwise
unmodified, as the directive required. SymPy output header now reads
`STAGE 212 ...` (exec log line 8; saved `.txt` confirmed). No `STAGE 195`
residue remains in either script.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `#pairs - binomial(5,2) = 0`, `#triples - binomial(5,3) = 0`
- `verified local full-simplex interval theorem on 3136 ordered integer samples`
- `pairwise budget - 120 = 0`, `full budget - 600 = 0`
- `All Stage 212 identities and interval theorems verified.`

**Mathematica:** exit=0. Notable lines:
- `PASS: M1 #pairs`, `PASS: M1 every primitive pair lies in three triples`,
  `PASS: M1 every primitive axis lies in six triples`
- `M3a local full-simplex sandwich = True` … `M3d no-improvement order = True`
  (all seven M3 `Resolve` results `True`, each followed by `PASS:`)
- `M4 pairwise total = 120`, `PASS: M4 10*12`, `PASS: M4 120+480`
- `Stage 212 Mathematica audit passed.`

Total PASS lines in the Mathematica log: 16 (M1×4, M2×2, M3×7, M4×3). All four
manifest groups M1–M4 are covered; none tautological.

**Output freshness:** confirmed. Both saved `.txt` outputs
(`mathematica/output/...txt` and `scripts/output/...txt`) carry mtime
2026-06-02 11:59, newer than both scripts' mtime 2026-06-02 11:38 — regenerated
post-fix.

## Material-change assessment

`material_change`: false.

F1 adds a brand-new independent verifier that reproduces the already-published
combinatorial/interval/budget results; it changes no derived value. F2 is a
display-only banner string. R1/R2 are notes-prose only (orchestrator already
reviewed; not in scripts scope). No downstream-relied result is altered.

## Side observations (non-blocking)

- The `.wl` defines an unused helper `expectTrue` (FullSimplify-based) and a
  one-arg `pass[]`; M2's two checks use an inline `If[... =!= True ...]` guard
  rather than `expectTrue`. Cosmetic only; behaviour is correct and guarded.
- The `.py` M3 enumeration uses small integer ranges (0–5 / 0–7), so it is a
  sampling argument, whereas the `.wl` `Resolve` is a true all-reals proof —
  the `.wl` is strictly the stronger of the two engines here, which is the
  intended dual-engine benefit. Not a defect.

## Verdict justification

Both findings are resolved: F2 is a clean one-line banner fix with no collateral
edits, and F1 delivers a genuinely independent Mathematica audit whose M3
theorems are proved by `Resolve[ForAll[...], Reals]` quantifier elimination
(distinct from the SymPy integer enumeration), whose M1 uses native
`Subsets`/`Binomial`/`Count`/`SubsetQ`/`MemberQ`, and whose M4 budget correctly
uses 120 (not the stale 188). Both engines exit 0 with all M1–M4 checks passing
and non-tautological, and the saved outputs were regenerated post-fix. No
regressions in the diff; no material change.
