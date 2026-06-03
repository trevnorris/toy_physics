---
unit_id: 225
batch: VII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T17:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 225

## Per-finding outcomes

### F1 — missing_verification_script (subtype: missing_mathematica)

**Classification:** resolved

**What changed:**
New independent Mathematica audit created at
`mathematica/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.wl`
(465 lines), covering the full claim manifest M1–M8. Exec log shows 61 PASS / 0 FAIL, exit 0. Output `.txt` regenerated at 17:05:50, newer than the `.wl` (16:51:09).

**Assessment:**
The script is a genuine independent re-derivation, not a SymPy transliteration, through a different decomposition on every load-bearing item:
- **M1 slopes** use `Coefficient[Normal[Series[expr,{eps,0,1}]], eps]` (`epsSlope`, wl:92,134-137) — Series/Coefficient, NOT SymPy's `diff(...).subs(eps,0)` (py:36-38). Residuals against independently-typed closed forms = 0.
- **M2** independently solves the compensation surface via `Solve` (wl:155-156), checks both `(D4/D0)D01` and `(u2^2-u4)D01` forms, AND verifies the one-pole reduction by substituting `D4 -> -3 D0 u2^2` into the *solved* `d41Comp` (wl:162-164) — plus an explicit negative control `expectNonZero["M2 wrong -2 coefficient is not an identity"]` (wl:165-168, residual `-(D01*D2^2)/D0^2` ≠ 0). So the `-3` is load-bearing in the second engine too.
- **M3** re-derives every drift via `epsSlope` of dressed primitives and asserts against closed forms typed in a *different algebraic grouping* than SymPy (e.g. z2Closed/z4Closed at wl:251-264). Non-tautological (residual zero is the test).
- **M6/M7** use native `MatrixRank`, `NullSpace`, `Det`. The sieve checks verify basis-invariant facts: rank==2, nullity==3 (wl:360-361,397-399), the null-residual `mixedMatrix . v == {0,0}` (wl:410-416, residual {0,0}), the 2×2 BdG determinant formula + nonzero sample value, and the scalar Ξ1(v1)=1.36026… (a basis-invariant quantity). The notes basis is reconstructed via `LinearSolve` on the free block but is checked by *lying in the nullspace*, not echoed.
- **M4/M8** carried targets verified against exact decimals (`expectClose`); kappa = `2 Sqrt[2]/Pi` (wl:293); d0Compat/kCompat/nullspace computed in-script from sample rules, not echoed.
No collateral edits; only the M1–M8 audit. resolved.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/...stage225..._sympy_audit.py:68-72`. The old `0==0` self-check
(`one_pole_D41 = (-3*u2**2)*D01; assert one_pole_D41 + 3*u2**2*D01 == 0`)
was replaced with:
```python
D4_one_pole = -3 * D0 * u2**2
one_pole_D41 = sp.simplify((D41_comp).subs(D4, D4_one_pole))
assert sp.simplify(one_pole_D41 - (-3 * u2**2) * D01) == 0
```

**Assessment:**
Genuinely de-tautologized. `D41_comp` (py:62) is the result of `sp.solve` on the compensation constraint — it equals `D4*D01/D0`, derived not hardcoded. `u2 = -D2/D0` (py:32) is independent of `D4`, so `D4_one_pole = -3*D0*u2^2 = -3*D2^2/D0` is a real substitution. The assert checks `D41_comp(D4 -> -3 D0 u2^2) == -3 u2^2 D01`, which holds iff the one-pole constraint `D4/D0 = u2^2 - u4 = -3 u2^2` truly reduces the surface to `-3 u2^2 D01`. Perturbing the coefficient (e.g. to `-2`) breaks the assert — the `-3` is load-bearing. The Mathematica M2 block independently confirms this with a positive check (residual 0) and a negative control (`-2` residual ≠ 0). Not a definitional round-trip. resolved.

### F3 — paper_misalignment (subtype: notes_contradicts_script), user-resolved

**Classification:** resolved

**What changed:**
Notes-only renumber in
`notes/stages/moving_throat_pde_stage225_..._sympy_audit.md` (42-line diff per diffstat).

**Assessment:**
Confirmed (notes scope is normally outside the verifier, but F3 is a user-authorized notes-only edit per the directive `## RESOLVED` block and I verified only the presence/absence of the stale tokens, not prose content). `grep -nE "Stage[ -]?24[012]|stage24[012]"` returns **NO_STALE_REFERENCES_FOUND**. Canonical references now present: self-reference is `Stage 225` (lines 1, 50, 542, 546); cited inputs are `Stage 223` (compatibility branch/sample) and `Stage 224` (transported ceilings); the §8 supporting-file citation at line 590 is the existing `...stage225..._sympy_audit.py`. The SymPy script's correct comments (223/224, py:195/240/392) are unchanged, and the paper card `stage_225.tex` + appendix `stage_appendix_part07.tex` are untouched (git status clean). Values unchanged. resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `On a one-pole base branch: D41 = -3*D01*D2**2/D0**2` (the reduced form, now actually derived via subs, not printed-back ansatz)
- `Verified concrete Stage 223 compatibility point:` with `D0 = 24.23730998862225`, `P0_target,compat = 0.00206979231806288270`
- `All Stage 225 symbolic and numerical checks passed.`

**Mathematica:** exit=0, 61 PASS, 0 FAIL. Notable lines:
- `PASS: M2 one-pole D41 reduction` and `M2 wrong -2 coefficient is not an identity residual = -((D01*D2^2)/D0^2)` → `PASS` (negative control fires correctly)
- `PASS: M4 D0` (actual 24.2373099886223…, residual ~5e-14), `PASS: M4 exact one-pole identity` (residual 0)
- `PASS: M7 mixed matrix rank two`, `PASS: M7 mixed nullity three`, `PASS: M7 Xi1(v1)` (1.36026097…)
- `All Stage 225 Mathematica checks passed.`

**Output freshness:** confirmed. SymPy `.txt` (17:05:50) > `.py` (16:46:15); Mathematica `.txt` (17:05:50) > `.wl` (16:51:09). Both regenerated post-fix.

## Material-change assessment

`material_change`: **false**.

No derived result changed. F1 added a corroborating second engine (same values, independently re-derived). F2 changed only *how* an existing identity is asserted — the printed reduced form `D41 = -3*D01*D2^2/D0^2` is unchanged. F3 is notes-only stage-number text; the carried literals (D0, K_compat, P0 from 223; the four budgets from 224) and all Ξ1 / nullspace / window values are unchanged. Downstream units depending on stage 225 outputs see identical numbers. No `upstream_stale` propagation is warranted on numerical grounds; the orchestrator's blanket `upstream_stale` flag (if any) can be treated as conservative here.

## Side observations (non-blocking)

- The new `.wl` uses `Rationalize`/`exactDecimal` to lift the carried decimal targets to exact rationals for `expectClose`/`expectVectorClose`; comparisons are at 16-digit precision with tol 1e-12 (M6 sample det relaxed to 1e-15, appropriate for its ~5e-5 magnitude). Consistent with the kappa = `2 Sqrt[2]/Pi` exact-arithmetic mandate. Not a finding.
- M7 reconstructs the notes-normalized basis via `LinearSolve[Transpose[freeBlock], UnitVector[3,j]]`; this is a presentation convenience and is guarded by the basis-invariant nullspace-residual check, so it does not smuggle in echoed SymPy basis vectors. Not a finding.

## Verdict justification

All three findings are resolved. F1's new Mathematica audit independently re-derives M1–M8 through Series/Coefficient, NullSpace/MatrixRank/Det, and exact `2 Sqrt[2]/Pi` arithmetic — a genuinely different decomposition than the SymPy diff/subs route — with 61 PASS / 0 FAIL and verifying basis-invariant facts (ranks, nullities, residuals, scalar Ξ1 values), not echoed vectors. F2 is genuinely de-tautologized: the one-pole `D4 = -3 D0 u2^2` constraint is substituted into the `solve`-derived `D41_comp` and the `-3` is load-bearing (negative control confirms on both engines). F3's notes renumber removed every stale `Stage 240/241/242` / `stage242` token (grep-confirmed) without touching scripts, paper card, or appendix; values unchanged. Exec logs pass, outputs are fresh, no regressions. `material_change: false`.
