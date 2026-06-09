---
unit_id: 173
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T17:20:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 173

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex re-authored the derivation portion of
`mathematica/moving_throat_pde_stage173_axisymmetric_loading_mismatch_mathematica_audit.wl`
(only file touched per the diff). Four distinct structural changes, all matching the directive:

1. **D[]-route slope derivation (wl L32-43).** The `d0A/d2A/d4A/n0A` named-ratio +
   `Series[...]`→`Coefficient[...,eps,1]` pipeline is deleted. `u21/u41/p1` are now derived by
   direct analytic differentiation of the inline ratios at `eps=0`, e.g.
   `u21 = FullSimplify[(D[-(d2 + eps*lam*d21)/(d0 + eps*lam*d01), eps] /. eps -> 0)/lam, ...]`.
   This is a genuinely distinct mechanism from the `.py`'s `series(...).diff(eps).subs(eps,0)`.
2. **SymPy-mirroring intermediate names gone.** No `d0A`, `d2A`, `d4A`, `n0A` anywhere in the
   file (confirmed by reading the full 102-line file). Canonical/branch helpers renamed to
   `canonicalBranchRules`, `branchU2Slope`, `branchU4Slope`, `staticPressureRatio`.
3. **Even-preserving collapse via linear coefficient extraction (wl L69-78).** `Solve[...]/First`
   is gone; `u21ZeroD21` and `d41Even` are now derived as
   `-(Numerator[Together[...]] /. var->0)/Coefficient[Numerator[Together[...]], var]` —
   manual linear root extraction, not `Solve`.
4. **Carry-forward Print block (wl L95-96).** The byte-identical 7-line `.py` carry-forward block
   is replaced with a single Mathematica-native summary line.

**Assessment:**
Correct and complete. The directive's four required changes are all present, and the diff shows no
collateral edits beyond the named derivation region (L29-96). All six `expectZero` targets are
unchanged (u2 slope identity, u4 canonical formula, P1/P0 formula, hidden-even residual,
`D21 + D01/9`, `D41 + D01/27`) — the deliverables and their RHS target formulas are byte-identical
to pre-fix. The new D[]-route assertions remain non-tautological: each slope LHS is independently
produced by analytic differentiation and compared to the paper's boxed operator-variable formula;
the general-form outputs (`(d01*d2 - d0*d21)/d0^2`, etc.) are nonzero, confirming genuine `eps`
dependence. The `.py` was not touched (empty py diff-stat). No new constants introduced.

## Exec log assessment

**SymPy:** exit=n/a. No sympy exec log exists for this unit — expected, since F1 only re-authored
the `.wl` and the `.py` is unchanged (git diff-stat on the `.py` is empty).

**Mathematica:** exit=0. Notable lines:
- `u2^(1) general = (d01*d2 - d0*d21)/d0^2` (nonzero general form → non-trivial slope)
- `PASS: u2 slope identity`, `PASS: u4 canonical formula`, `PASS: P1/P0 formula`
- `PASS: hidden-even residual`, `PASS: D21 + D01/9`, `PASS: D41 + D01/27`
- `D21 from u2^(1)=0 = -1/9*d01`, `D41 on even-preserving branch = -1/27*d01`
- `Mathematica summary: verified slope identities, even-preserving collapse, and Xi_load lanes above.`
All six `expectZero` residuals = 0; `Xi_load = -(d01/d0) + n01/n0` with lanes (1, 1/2, -1) intact.

**Output freshness:** confirmed. `.wl` source mtime 1780958927 < `.txt` output mtime 1780960402,
so the saved Mathematica output was regenerated post-fix.

## Material-change assessment

`material_change`: false. The edit is method-only: the derived slope laws, canonical formulas,
even-preserving collapse values, and `Xi_load` lanes are all preserved bit-for-bit (same residuals,
same printed forms, same target RHS in every `expectZero`). No derived result that downstream units
depend on has changed. No downstream re-audit needed on math grounds.

## Side observations (non-blocking)

- The general-form outputs printed by the new D[]-route (L10-12 of the log) match the auditor's
  recorded engine cross-check table exactly, an independent confirmation that the re-derivation
  reaches the same algebra by a different path.

## Verdict justification

The single finding (F1, mathematica_transliteration, medium) is fully resolved. The `.wl` now
derives the first-order slopes by direct analytic `D[ratio, eps]/.eps->0` differentiation — a
mechanism genuinely distinct from the `.py`'s truncated-series coefficient extraction — the
`Series`/`Coefficient` pipeline and the SymPy-mirroring `d0A/d2A/d4A/n0A` intermediates are gone,
the even-preserving collapse uses linear coefficient extraction instead of `Solve[...]/First`, and
the byte-identical carry-forward block is replaced by a native one-liner. All six `expectZero`
checks pass with unchanged targets, the exit code is 0, the output is fresh, and the `.py` is
untouched. The deliverables are preserved, so `material_change` is false. Verdict: `verified`.
