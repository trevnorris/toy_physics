---
unit_id: 206
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 206

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/..._sympy_audit.py:130-147`. The old Section V (which bound `H0_stage205 = sp.log(Phi0)` and then re-typed the identical denominator, so `TauStage189 - TauStage188` was identically 0 before `simplify`) is replaced by an equality between two independently-built operands:
- left: `TauBracketLog = sp.simplify(Tau.subs({H0: sp.log(Phi0), k: -L0, c: L1}))` — the Section-I generic root map `Tau = 2*H0/(k + sqrt(k**2 - 2*c*H0))` (positive-`k` form, line 36) specialized to the log branch via `k -> -L0`.
- right: `TauLog2 = -2*sp.log(Phi0)/(L0 + sp.sign(L0)*sp.sqrt(L0**2 - 2*L1*sp.log(Phi0)))` — a fresh transcription of the Stage-239 appendix predictor in the oriented `K0 + sgn(K0)*sqrt(...)` form, with `sign(L0)` kept general.
The check is `expect_zero("Stage 206/239 log-predictor collapse", sp.simplify(TauBracketLog - sp.refine(TauLog2, sp.Q.negative(L0))))`.

**Assessment:**
The two operands are genuinely non-aliased: they are assembled from two different closed forms (`2*H0/(k+sqrt)` vs `-2*H0/(K0+sgn(K0)*sqrt)`), and the exec output prints them in visibly distinct shapes — left as `-2*log(Φ0)/(L0 - sqrt(L0² - 2*L1*log(Φ0)))`, right as `-2*log(Φ0)/(L0 + sqrt(...)*sign(L0))`. Reconciliation requires the `sgn` convention to be resolved via `refine(..., Q.negative(L0))`, so the check exercises the sign convention rather than assuming `sgn(L0) = -1`. The check is falsifiable (deleting a distinctive subterm of either side breaks the simplification to 0). Output line "Stage 206/239 log-predictor collapse = 0" prints residual 0. Resolved.

### F2 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New untracked file `mathematica/..._mathematica_audit.wl` (233 lines), covering claim manifest M1–M7 with hard-failing PASS/FAIL checks (`fail` → `Exit[1]`).

**Assessment — independence (load-bearing):**
The `.wl` is a genuinely independent derivation, not a port of the `.py`:
- **M1 (root map):** `solveRoots = stripConditional[tau /. Solve[q[c, tau] == 0, tau]]` then physical branch chosen by `SelectFirst[solveRoots, TrueQ[FullSimplify[linearLimit[#] == H0/k, ...]] &]` — i.e. roots are obtained by `Solve` and selected by a `Limit[..., c->0] == H0/k` test. The SymPy script instead *seeds* with the typed closed form `Tau = 2*H0/(k+sqrt(...))`. Different decomposition.
- **M7 (turning root):** `turningRoots = stripConditional[tau /. Solve[turningEquation, tau]]` then `SelectFirst[..., TrueQ[Refine[# > 0, ...]] &]`, vs the SymPy script directly typing `sqrt(2*H0/a)`.
- **M3/M6 (sign + width):** `dRootDc = D[rootMap, c]` with `expectTrue["M3 curvature derivative sign", dRootDc > 0, branchAssumptions[c]]` via `Refine` (a sign assertion the SymPy file does not make); M6 uses `Series` + separate `Coefficient[widthSeries, eta, 1]` and `Coefficient[widthSeries, eta, 2]` extraction (eta² cancellation checked explicitly), a different mechanism than SymPy's `series(...).removeO()` minus-leading comparison.
Roots via `Solve`, derivatives via `D[]`, width via `Series`/`Coefficient`, limit via `Limit`, signs via `Refine` — exactly the native-primitive decomposition the directive's anti-transliteration guard required, distinct from the SymPy choreography. Not a transliteration. Mathematica exec exits 0 with 21 PASS, 0 FAIL covering M1–M7. Resolved.

### F3 — paper_misalignment (USER-RESOLVED direction (a): additive script checks, both engines, no prose edit)

**Classification:** resolved

**What changed:**
Both engines gained M(F3a) pairwise ray-ordering and M(F3b) admissibility checks. SymPy: `..._sympy_audit.py:149-226` (Section VI). Mathematica: `..._mathematica_audit.wl:176-230` (F3a/F3b blocks).

**Assessment — M(F3a) non-vacuity (load-bearing):**
- SymPy encodes the ordering by *constructive* slack: `tau_hi_a = tau_star_a + slack_a_hi`, `tau_lo_b = tau_hi_a + strict_sep`, `tau_star_b = tau_lo_b + slack_b_lo` with `slack_* >= 0` and `strict_sep > 0`, giving `ordering_gap = slack_a_hi + slack_b_lo + strict_sep`; `ask(Q.positive(gap))` is `True` and `ask(Q.nonpositive(gap))` is `False`. Non-vacuity is a concrete relaxed counterexample (separation dropped): hypotheses `0<=1<=2` and `0<=1/2<=2` hold while the conclusion `1 < 1/2` is `False`. Output: "without separation counterexample hypotheses = True", "... conclusion = False" — separation is load-bearing.
- Mathematica encodes it by a *different mechanism*: `Resolve[ForAll[vars, Implies[premise, conclusion]], Reals]`. With `tauHiA < tauLoB` → `True` ("F3a pairwise certified interval ordering = True"); with that hypothesis dropped → `False` ("F3a ordering without separation is not tautological = False"). The implication is True with separation and NOT-True without it, in both engines, via two distinct mechanisms (Reduce/Resolve-over-Reals vs SymPy assumption discharge). Non-vacuous.

**Assessment — M(F3b) discrimination:**
Both engines assert three cases: monotone-admissible bracket → True, turning-admissible bracket → True, single-clause violation (`tau_hi = 2 > T_valid = 1`, all other clauses satisfied) → False. The predicate flips on exactly the one violated clause (`tau_hi <= T_valid`), so omitting/flipping a clause changes the result. Both True and False cases asserted in `.py` (`if monotone_bad is not sp.S.false: raise`) and `.wl` (`expectFalse[...]`). Discriminating.

**No prose edit:** `git status` confirms only `scripts/...sympy_audit.py` (modified, tracked) and `mathematica/...audit.wl` (new, untracked) changed; `git diff --name-only HEAD` over `paper/` and `notes/` returns nothing for 206. The captured diff patch touches only the `.py` (the `.wl` is a new file, so it does not appear in the working-tree diff against HEAD, consistent with untracked status). No paper.tex / notes edit was made. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Stage 206/239 log-predictor collapse = 0` (F1 residual zero, two independently-built operands).
- `pairwise ordering theorem = True`, `pairwise ordering negation satisfiable = False` (F3a discharged).
- `pairwise theorem without separation counterexample conclusion = False` (F3a non-vacuity).
- `local sieve monotone admissible case = True` / `turning admissible case = True` / `single-clause violation case = False` (F3b discriminates).
- `STAGE 189 SYMPY AUDIT PASSED` (banner text still reads "STAGE 189"; cosmetic only — see side observations).

**Mathematica:** exit=0, 21 PASS / 0 FAIL. Notable lines:
- `M1 selected root = (k - Sqrt[-2*c*H0 + k^2])/c` (root obtained by Solve + branch-select, not typed) → `PASS: M1 closed root - Solve-selected branch`.
- `PASS: M3 curvature derivative sign`, `PASS: M4 strict endpoint descent on simple-root branch` (Refine sign claims).
- `F3a pairwise certified interval ordering = True` → PASS; `F3a ordering without separation is not tautological = False` → PASS (non-vacuity).
- `F3b ... = True/True/False` → all PASS (discrimination).
- `STAGE 206 MATHEMATICA AUDIT PASSED`.

**Output freshness:** confirmed. `scripts/output/...sympy_audit.txt` and `mathematica/output/...audit.txt` both mtime 2026-06-02 11:32:11, newer than the SymPy script (11:24:37) and the `.wl` (11:26:34). Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

No edit changed a derived result. F1 replaced a vacuous self-comparison with a substantive identity that still evaluates to residual 0 — the carried Stage-239 log predictor and the Section-I root map are confirmed equal, no constant changed. F2/F3 are purely additive verification (new `.wl`, new Section VI). The load-bearing root map itself was already correctly verified pre-fix and is unchanged. No downstream unit relies on a value that moved.

## Side observations (non-blocking)

- The SymPy script's top/bottom banners still read "STAGE 189 — CERTIFIED RAY RANKING" / "STAGE 189 SYMPY AUDIT PASSED" (lines 28, 228) and the exec-log header echoes them, even though the unit is 206 and the Section V/VI prose was correctly updated to "Stage 239". Pre-existing mislabel, cosmetic, does not affect any assertion. The Mathematica banner correctly says "STAGE 206". Not blocking.
- `slack_sep` is declared (line 154) but unused; `strict_sep` carries the separation. Harmless dead symbol.

## Verdict justification

All three findings are resolved. F1's Section-V collapse now compares two independently-constructed expressions (Section-I root map specialized to the log branch vs a fresh Stage-239 predictor transcription) with `sign(L0)` kept general and reconciled by `refine` under `L0 < 0`; the operands are non-aliased, the check is falsifiable, and the residual prints 0. F2's new `.wl` is a genuine second engine that derives M1–M7 via Solve/Limit/D/Series/Coefficient/Refine in a different decomposition than the SymPy seed-the-closed-form approach — not a transliteration. F3 adds the pairwise-ordering theorem (proved True with the separation hypothesis and demonstrated NOT-True without it, in both engines via different mechanisms) and a discriminating admissibility predicate (True on admissible, False on a single-clause violation) to both engines, with no paper.tex or notes edit. Both engines exit 0 (Mathematica 21 PASS / 0 FAIL), and outputs are fresh. Verdict: verified.
