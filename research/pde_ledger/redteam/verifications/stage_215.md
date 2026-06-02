---
unit_id: 215
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 215

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created the new file `mathematica/moving_throat_pde_stage215_full_primitive_quadruple_ranking_theorem_mathematica_audit.wl` (257 lines) at the registered Target path. It implements the full M1–M7 claim manifest with guarded `expectZero`/`expectTrue` helpers that call `Exit[1]` on failure (lines 36–46), and exits 0 only after `Stage 215 Mathematica audit passed.` The orchestrator's re-run shows exit 0 with 19 PASS lines.

**Assessment:**
Correct and complete. Manifest coverage maps cleanly to the 19 PASS lines: M1 (5 PASS, `mathematica:78-82`), M2 (2, `:98-105`), M3 (1, `:109-120`), M4 (2, `:124-147`), M5 (2 — pairwise `:151-162` and the new 5-interval unique-winner `:169-191`), M6 (3, `:195-230`), M7 (4, `:234-251`). Every manifest item is an explicit guarded check; none is tautological — the `Resolve[ForAll[...], Reals]` propositions are falsifiable (a false proposition would return a non-`True` residual and trip `fail` → `Exit[1]`), and M7 reconstructs the constants from factors rather than asserting bare literals.

**Genuine independence (load-bearing) — confirmed.** The `.wl` derives each result by a different route than the SymPy `.py`:
- Interval/order theorems M3–M6: `.wl` uses continuous quantifier elimination, e.g. `Resolve[ForAll[{betaLo,betaHi,iotaLo,iotaHi,bBest,iBest}, Implies[... , Min[betaLo,iotaLo] <= Min[bBest,iBest] <= Min[betaHi,iotaHi]]], Reals]` (`mathematica:109-120`). The `.py` instead enumerates nested integer `range()` lattices, e.g. `for beta_lo_v in range(0,6): ... if not (lo_full <= simplex_best <= hi_full): raise` (`scripts:97-109`). Continuous QE vs. finite integer sampling = genuinely different proofs (and the continuous route also closes the `.py` sampling-coverage gap, as the directive intended).
- Combinatorics M1: `.wl` uses `Subsets[axes,{3}]`/`Subsets[axes,{4}]` with `Count[triples, tri_ /; SubsetQ[quad,tri]]` and `Boole[MemberQ[...]]` (`mathematica:55-69`). The `.py` uses `itertools.combinations(axes,3)` with `set(tri).issubset(quad)` (`scripts:39-67`). Different primitives, not a verbatim port of the Python loop.
- Unique-winner M5(ii): `.wl` states it directly over five real intervals as a 15-variable `Resolve[ForAll[..., Implies[premise, conclusion]], Reals]` per star with `minOtherLower = Min @@ Delete[lVars, star]` (`mathematica:169-190`); the `.py` does integer-lattice enumeration over `itertools.product(interval_samples, repeat=5)` (`scripts:145-159`).

This is independent derivation, not transliteration. No collateral edits in the `.wl` (it is a new file).

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy diff (`stage_215_diff.patch`) shows exactly the three directive-named edits, no collateral:
1. Banner `STAGE 198` → `STAGE 215` (`scripts:35`). Confirmed in the exec log header (`STAGE 215 — FULL PRIMITIVE-QUADRUPLE RANKING ...`).
2. New explicit min-over-others unique-winner check (`scripts:145-159`): `other_lows = [...]; if intervals[star][1] < min(other_lows): ... assert x_star < x_p`. Exercised on 542760 samples per the log.
3. Budget literals reconstructed: `support_le3_budget = 10*12 + 10*48` and `quad_eval_per_envelope = (3*3*3)*2` (`scripts:198-199`), so the surviving `expect_zero` product checks now confirm reconstructed (not self-asserted) values.

**Assessment:**
Both previously-partial sub-claims are now substantively exercised and falsifiable. (a) The unique-winner check is the explicit "minimum over the other four quadruples" statement, distinct from the pre-existing pairwise loop, and is mirrored as the load-bearing M5(ii) `Resolve` proposition in the `.wl`. (b) The budget constants are reconstructed from their paper-stated factor decompositions in both engines (`.wl` M7: `Product[degreePattern]` = 54 and `Dot[{10,10},{12,48}]` = 600). Neither new check is vacuous: the unique-winner assertion would raise on any counterexample where a non-winner value undercut the certified winner, and the budget reconstructions would fail if any factor were mistyped. The optional items 1–2 were applied rather than blocked, so the SymPy engine independently closes the same gaps the M5/M7 manifest carries.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `verified unique certified winner theorem on 542760 min-over-others samples` (the new M5(ii) check)
- `quadruple interior budget - 540 = 0` and `full support<=4 budget - 1140 = 0` (now from reconstructed factors)
- `All Stage 215 identities and interval theorems verified.`
- Header banner reads `STAGE 215 — ...` (stale 198 corrected).

**Mathematica:** exit=0. 19 PASS lines, 0 FAIL. Notable:
- `PASS: M5 unique certified winner over five quadruple intervals` with `M5 unique-winner quantified results = {True, True, True, True, True}`
- `PASS: M7 per-envelope 3*3*3*2 factor product`, `PASS: M7 support<=3 10*12+10*48 factor sum`, `PASS: M7 full support<=4 budget - 1140`
- `Stage 215 Mathematica audit passed.`

**Output freshness:** confirmed re-generated post-fix. `.py` mtime 11:05:30 and `.wl` mtime 11:07:23; both saved `.txt` outputs (sympy + mathematica) mtime 11:12:53 — newer than both scripts.

## Material-change assessment

`material_change`: false. No derived result changed — the combinatorial counts (10/5), the certified-interval/min structure, and the budget constants (54, 600, 540, 1140) are confirmed, not altered. The reconstructions reproduce the same numbers the `.py` previously hardcoded. Adding the second engine and the new checks introduces no downstream-relied value change.

## Side observations (non-blocking)

- The SymPy script's M7 prints/asserts `quadruple interior budget - 540` and `full support<=4 budget - 1140` but does not separately echo a `54 = 3*3*3*2` / `600 = 10*12+10*48` equality line (the reconstruction is implicit in the variable assignment). This is fine — the `.wl` M7 makes both factor equalities explicit guarded checks. No action needed.
- The directive front-matter shows `findings_applied: 2`, `findings_blocked: 0`; consistent with both `## Applied:` blocks. No `## Blocked:` block present.

## Verdict justification

Both findings are `resolved`. F1's new `.wl` exists at the registered path, exits 0 with 19 PASS lines covering the full M1–M7 manifest, uses guarded `Exit[1]`-on-failure checks, and is a genuinely independent derivation (continuous `Resolve[ForAll,Reals]` quantifier elimination + `Subsets`/`Count` combinatorics + `Product`/`Dot` budget reconstruction) versus the SymPy integer-lattice/`itertools` route — not a transliteration. F2's two previously-partial sub-claims are now substantively and falsifiably exercised in both engines (the explicit min-over-others unique-winner statement and the factor-reconstructed budget constants), and the stale STAGE 198 banner is corrected. Both engines exit 0, outputs are fresh, the diff shows only directive-named edits with no collateral or regressions, and no downstream-relied result changed. Verdict: verified.
