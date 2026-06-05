---
unit_id: 043
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 043

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
SymPy `.py` only (per directive — the `.wl` self-labels were already canonical):
- line 3 docstring `Moving-throat PDE — Stage 26 SymPy audit.` → `... Stage 43 SymPy audit.`
- subbanners `26.1`–`26.5` → `43.1`–`43.5` (script lines 57/83/117/150/164, now `43.k`).

**Assessment:**
The captured diff shows exactly these six numeric-token swaps and nothing else for F1; surrounding strings are preserved verbatim (e.g. the already-canonical `STAGE 43` banner at line 55 is untouched, and there are no cross-refs in this `.py`). `grep -E "Stage 26|26\.[0-9]|STAGE 26|STAGE 026"` over both scripts returns no matches (exit 1), confirming no residual stale self-labels. The refreshed transcripts now print the canonical labels: sympy `.txt` line 8 `STAGE 43 — ...`, subbanners `43.1`–`43.5` (txt lines 12/34/60/94/104); mathematica `.txt` line 8 `STAGE 043 — ...`. All PASS/`= 0` lines remain. Label-only, no result change.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
The redundant baseline block (which substituted the same literal `8/pi^2` into both sides of the already-proven structural identity) was replaced in both engines with a check that DERIVES the baseline from the stage's frozen overlap constant:
- SymPy (lines 145-148): `sigma_value = 88/(9 pi^2)`; `B_value = simplify((9/11)*sigma_value)`; `expect_zero("baseline B = 8/pi^2 from frozen sigma", B_value - 8/pi^2)`.
- Mathematica (lines 127-130): `sigmaValue = 88/(9 Pi^2)`; `bValue = FullSimplify[(9/11) sigmaValue]`; `expectZero["baseline B = 8/Pi^2 from frozen sigma", bValue - 8/Pi^2]`.
The structural-form check (sympy 141 / wl 123, `M_supp structural form (free baseline)`) is retained verbatim. The trailing `M_supp at baseline B = 8/pi^2` evaluation is kept as a sanity print but now uses the derived `B_value`.

**Assessment:**
Non-tautological and genuinely exercises the value. `B_value` is COMPUTED from `sigma_value` via `(9/11)*88/(9 pi^2) = 88/(11 pi^2) = 8/pi^2`; the assertion then subtracts the independent literal `8/pi^2`. If `kappa0^2` were mis-stated (wrong ratio `kappa1^2/sigma`, or wrong `sigma`), `B_value - 8/pi^2` would be a concrete nonzero rational-times-`1/pi^2` and the check would fail — so it can fail and thus tests something. The value `8/pi^2 = (9/11)*88/(9 pi^2)` matches the notes' frozen constants (`sigma = 88/(9 pi^2)`, `kappa0^2/sigma = 9/11`) exactly, so no new paper_misalignment is introduced. Both engines emit the derived value then `= 0` PASS (sympy txt 77-78; wl txt 58-60). The structural-form PASS still prints (sympy txt 76; wl txt 56-57). The Mathematica directive labeled the old block `F3`; the diff corrects the comment to `F2` — cosmetic, no behavior change.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `baseline B = kappa0^2 = (9/11) sigma = 8/pi**2` (txt 77) — derived value printed.
- `baseline B = 8/pi^2 from frozen sigma = 0` (txt 78) — new check passes.
- `M_supp structural form (free baseline) = 0` (txt 76) — prior structural check intact.
- `mismatch formula = 0` (txt 110) — deliverable 7 unaffected (the legitimate "mismatch" derivation, not a failure token).

**Mathematica:** exit=0. Notable lines:
- `baseline B = kappa0^2 = (9/11) sigma = 8/Pi^2` (txt 58) and `PASS: baseline B = 8/Pi^2 from frozen sigma` (txt 59-60).
- `PASS: M_supp structural form (free baseline)` (txt 57) — prior check intact.
- `Stage 043 Mathematica audit passed.` (txt 90).
- The two `Limit::alimv` warnings (txt 28, 40) are the pre-existing endpoint-limit warnings, not errors; both endpoint checks still PASS (txt 31-34).

**Output freshness:** confirmed fresh. Scripts: `.py` 2026-06-05 07:57:22, `.wl` 2026-06-05 07:57:31. Outputs: both `.txt` 2026-06-05 08:08:02 — newer than their scripts. Exec-log timestamps (08:07:5x) corroborate the post-fix re-run.

## Material-change assessment

`material_change`: false. F1 is label/banner only. F2 only re-points an existing assertion from a guaranteed-true substitution to a genuine derivation of the same already-correct value `8/pi^2`; no derived result that downstream units consume is altered (the baseline `kappa0^2 = 8/pi^2` is unchanged — it is now derived rather than assumed). De-tautologizing an already-correct value is not a material change.

## Side observations (non-blocking)

- The in-source `# F2:` comment in the diff replaces the old `# F3:` comment; the directive's old block was mislabeled `F3`. Purely cosmetic, correctly aligns with this directive's finding numbering.
- The trailing `M_supp at baseline B = 8/pi^2` evaluation remains structurally redundant after the structural-form check (as the original report noted), but it is now a labeled sanity print downstream of the load-bearing new check and is harmless; the directive explicitly permitted keeping it. Not a blocker.

## Verdict justification

Both findings are resolved. F1 is a clean label-only swap (six `26→43` tokens in the `.py` docstring/subbanners, no `.wl` change needed), verified by the diff, a clean grep for residual stale labels, and the refreshed transcript banners. F2's new check genuinely derives `B = kappa0^2 = (9/11)*sigma = 8/pi^2` from the frozen `sigma = 88/(9 pi^2)` and asserts equality with the independent literal — non-vacuous (it can fail if the overlap constants were mis-stated), matches the notes exactly, with the structural-form check retained. Both engines exit 0, all prior PASS lines remain, and outputs are freshly regenerated. No regressions in the diff. Verdict: verified, material_change false.
