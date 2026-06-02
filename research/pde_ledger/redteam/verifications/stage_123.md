---
unit_id: 123
batch: IV.3
verifier_model: claude-opus-4-8
verify_date: 2026-06-01
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 123

> Retro-sweep dual-engine retrofit. Stage 123 was already audited+verified in batch IV.3
> (SymPy-only; the `228 → 160` notes typo fix and banner relabel landed there). The single
> finding here is the directive's **F1 = missing_mathematica**: Codex added a brand-new
> independent-route Mathematica `.wl`. Verified against the DIRECTIVE F1 + claim manifest
> M1–M4 (no fresh auditor report exists; stale IV.3 report content ignored).

## Per-finding outcomes

### F1 — missing_mathematica (dual-engine retrofit)

**Classification:** resolved

**What changed:**
Codex created the new second engine
`mathematica/moving_throat_pde_stage123_parent_normalized_branch_values_mathematica_audit.wl`
(untracked new file, 126 lines). No other file was modified — the reference SymPy
`scripts/moving_throat_pde_stage123_parent_normalized_branch_values_sympy_audit.py` is
unchanged vs HEAD, and `paper/` and `notes/` are untouched (confirmed via `git diff --stat HEAD`
and `git status --porcelain`). The captured `stage_123_diff.patch` is empty precisely because
the `.wl` is untracked (git diff does not show untracked files); the file itself was read
directly and the saved output `.txt` was re-generated post-creation.

**Assessment:** Correct and complete. Detailed checks below.

**1. Exec + PASS count.** Mathematica exec log exits 0 with exactly 6 substantive PASS lines:
`Xi_v law`, `Xi_T law`, `Xi_v(F1)`, `Xi_T(nat)`, `Xi_T(-)`, `Xi_T(+)` — meets the ≥6 bar
(M1, M2, M3, M4×3). Closes with "Stage 123 Mathematica audit passed." Each PASS is gated by
`expectZero` which `Exit[1]`s on a non-zero residual, so the passes are real, not cosmetic.

**2. Genuine independent route (not a transliteration).** The two laws EMERGE from inversions:
- M1: `rParent = lambda/Sqrt[Ks*Kq]` (`.wl:68`) is inverted for `vw` via
  `Reduce[rr == rParent && Element[vw, Reals], vw, Reals]` + a `LogicalExpand`/`Cases`
  branch extractor (`uniqueBranch`, `.wl:34-45,73-74`) — a different decomposition than the
  `.py`'s `sp.solve(sp.Eq(r, r_from_parent), vw)[0]`. The derived `v_w0(r)` and
  `Xi_v(r) derived = (-3*Sqrt[3/10]*Pi^(3/2)*rr)/16` are printed as intermediate results,
  then checked against the anchored `xiVExpected = -(3 Sqrt30 Pi^(3/2)/160) rr` via
  `expectZero` (residual = 0). `Sqrt[3/10] = Sqrt[30]/10`, so the two forms are identical.
- M2: `gParentLocked` (`.wl:83-84`) is inverted for `Tm` via
  `Reduce[gg == gParentLocked && Tm > 0 && gg > 0, Tm, Reals]` + `uniqueBranch`, then
  substituted into `xiTDef`. The derived `Xi_T(g) = 3*Sqrt[3/(10 Pi)]/gg` equals the anchored
  `xiTExpected = (3 Sqrt30/(10 Sqrt Pi))/gg` (`Sqrt[3/(10 Pi)] = Sqrt[30]/(10 Sqrt[Pi])`).
The numerics (M3/M4) are COMPUTED from the derived `xiVFromParent`/`xiTFromParent` at the
upstream-anchored `rF1`, `gMinus`, `gPlus` (`.wl:97-120`) via `printNumeric` =
`N[RootReduce[...], 20]` — not re-typed targets. The closed forms are NOT self-checked against
re-typed strings; they are derived then matched to the anchored building blocks. Anti-transliteration
requirement satisfied.

**3. THE CRITICAL GUARD — un-squared λ sign: PASS.**
- `lambda = -(8*Sqrt[2]/3)*q*vw*a^2*ell*Sqrt[L]` (`.wl:66`) — NEGATIVE λ convention (stage 118),
  LINEAR in `vw`.
- `rParent = lambda/Sqrt[Ks*Kq]` (`.wl:68`) — λ is un-squared (linear), so the minus propagates.
  λ is NOT squared anywhere (no `r_c = λ²/(K_s K_q)` confusion).
- `Element[vw, Reals]` only (`.wl:58`); `vw` is NOT assumed positive — every other symbol is
  `> 0`. Matches the `.py` declaration `v_w0` real-only. The reality-only assumption is what lets
  the inversion keep the correct sign.
- Result: derived `Xi_v(r) = (-3 Sqrt[3/10] Pi^(3/2) rr)/16` carries the LEADING MINUS, and
  `Xi_v(F1) numeric = -1.01675633282525946...` is NEGATIVE — exactly the load-bearing sign.
  Had λ been squared or the minus dropped, this would have come out +1.0168 and the
  `Xi_v law` residual check would have FAILED. It did not.
- Healing lock `c_s -> hbar/(2 m ell)` is applied ONLY inside the Ξ_T inversion (`.wl:84`,
  on `gParentUnlocked`), not on the Ξ_v branch — giving Ξ_T its `Sqrt[Pi]`. Correct per stage 118.

**4. Cross-engine agreement: PASS.** All four numerics match the SymPy reference to ≥18 digits,
including the Xi_v sign:

| Quantity   | SymPy ref                | Mathematica (`.wl`)                          |
|------------|--------------------------|----------------------------------------------|
| Xi_v(F1)   | -1.0167563328252594644   | -1.01675633282525946441973678979837042557    |
| Xi_T(nat)  |  0.92705808485565499282  |  0.92705808485565499282126256474757339116    |
| Xi_T(-)    |  1.2229751770146391627   |  1.2229751770146391627143360114746526363     |
| Xi_T(+)    |  0.33133452164460908424  |  0.33133452164460908424123901424332098173    |

The negative Xi_v matches the `.py`'s residual form (`Xi_v + 3 sqrt(30) pi^(3/2) r/160`, i.e.
`Xi_v = -...`). Both engines independently produce the same closed forms and decimals.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy run was part of this retrofit (the `.py` is the unchanged
reference engine; `stage_123_sympy.log` is empty/absent — expected, only the new `.wl` was run).
The committed SymPy output `.txt` was used to confirm cross-engine agreement.

**Mathematica:** exit=0. Notable lines:
- `Xi_v(r) derived = (-3*Sqrt[3/10]*Pi^(3/2)*rr)/16` then `Xi_v law = 0` / `PASS: Xi_v law`
- `Xi_T(g) derived = (3*Sqrt[3/(10*Pi)])/gg` then `Xi_T law = 0` / `PASS: Xi_T law`
- `Xi_v(F1) numeric = -1.0167563328252594644...` / `PASS: Xi_v(F1)` (NEGATIVE — guard holds)
- `Xi_T(-) numeric = 1.2229751770146391627...` / `PASS: Xi_T(-)`
- `Stage 123 Mathematica audit passed.`

**Output freshness:** confirmed. The `.wl` mtime is `2026-06-01 15:53:17`; the saved
`mathematica/output/...stage123...txt` mtime is `2026-06-01 16:06:06` (newer), and its contents
match the exec log line-for-line. The output was re-generated post-creation.

## Material-change assessment

`material_change`: **false**. This adds a second (Mathematica) verification engine only; no
derived result changed. The `.wl` reproduces the existing SymPy laws and the four branch
numerics exactly — no value that any downstream unit depends on was altered. No downstream
units need re-audit on account of this change.

## Side observations (non-blocking)

- The `Xi_v(F1)` / `Xi_T(*)` `expectZero` checks compare the derived law evaluated at the branch
  point against the *expected* law evaluated at the same branch point (`.wl:102-119`). On its own
  that would be slightly self-referential, BUT the substantive independence comes from (a) the
  preceding `Xi_v law` / `Xi_T law` checks that already prove derived ≡ expected, and (b) the
  `printNumeric` lines that compute the actual ≥15-digit decimals from the DERIVED expression
  and match the SymPy reference. The chain is therefore a genuine check, not tautological. Noted
  only for transparency — not blocking.
- `gMinus`/`gPlus` are imported in the stage-122 surd form `(2 Sqrt[4107-100 Pi^2] ± 37 Sqrt[3])/(20 Pi)`
  (`.wl:99-100`), whereas the `.py` writes them as `rF ∓ sqrt(1+rF^2)/2`; both forms numerically
  agree (verified via the matching decimals), and importing upstream anchors as givens is expected.

## Verdict justification

The single retro-sweep finding (F1 missing_mathematica) is fully resolved. Codex added a genuine
independent-route Mathematica `.wl` that derives both laws via `Reduce`-based inversions (distinct
decomposition from the SymPy `sp.solve[0]` route), preserves the load-bearing un-squared negative-λ
sign so `Xi_v(F1)` comes out NEGATIVE (-1.01675633...), applies the healing lock only inside the
Ξ_T inversion, exits 0 with all six PASS lines, and agrees with the SymPy reference to ≥18 digits
on all four branch numerics. The reference `.py` and all paper/notes are untouched, and the saved
output was re-generated post-fix. No regressions. Verdict: verified; material_change false.
