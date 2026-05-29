---
unit_id: 122
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T13:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 122

Scope note: this is the second-pass remediation of the single carried finding F1
(`tautological_check`) on the traction-ratio assertions. The original auditor report
listed three findings (F1 paper_misalignment routed to user, F2 insufficient_verification,
F3 missing_mathematica informational); the directive for this iteration carries exactly one
finding — F1 = the tautological traction-ratio check — and that is what is verified here.
The F2-added assertions (compensation quadratic, defect closed form, natural-off-compensation
nonzero) are present and re-checked as regression guards.

## Per-finding outcomes

### F1 — tautological_check (traction-ratio identities)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:43-57`
and `:73-74`. The reciprocal definitions
`T_ratio_minus = sp.simplify(1/gminus)` / `T_ratio_plus = sp.simplify(1/gplus)` were
replaced with an independent construction from the stage-119 proportionality law `g = C/T_m`:
a fresh symbol `C = sp.symbols("C", positive=True)` (line 48), `Tm_nat = C/g_nat`,
`Tm_minus = C/gminus`, `Tm_plus = C/gplus` (lines 49-51), and
`T_ratio_minus = sp.simplify(Tm_minus/Tm_nat)`, `T_ratio_plus = sp.simplify(Tm_plus/Tm_nat)`
(lines 52-53). The two assertions (lines 73-74) now read
`expect_zero("traction ratio (-) identity", T_ratio_minus - 1/gminus)` and
`expect_zero("traction ratio (+) identity", T_ratio_plus - 1/gplus)`, with labels unchanged.
The diff (`stage_122_diff.patch`) is exactly these two hunks — no stray edits, no collateral
changes elsewhere in the script.

**Assessment:**

1. **No X−X re-tautologization — confirmed.** The left-hand side `T_ratio_minus` is built
   from the primitives `C`, `gminus`, and `g_nat` (the line-34 natural-branch ansatz
   `g_nat = sp.Integer(1)`). The comparison target `1/gminus` appears only as the right-hand
   side of the assertion residual, never as the primitive the LHS was defined from. The
   former reciprocal-by-construction loop is gone.

2. **C genuinely cancels — confirmed.** `(C/gminus)/(C/g_nat)` reduces to `g_nat/gminus`
   only because the SAME symbolic `C` appears in numerator and denominator. The exec log
   print line `T_m(-)/T_m(nat) = -20*pi/(-2*sqrt(4107 - 100*pi**2) + 37*sqrt(3))` shows no
   surviving `C` — the constant cancelled. This is the substantive content: if the traction
   normalization were branch-DEPENDENT (`C_minus/gminus` vs `C_nat/g_nat`), the ratio would
   carry `C_minus/C_nat` and the residual against `1/gminus` would be nonzero, so the check
   now exercises branch-independence of the background-fixed traction constant.

3. **Can it now fail? — yes.** The residual simplifies to
   `g_nat/gminus - 1/gminus = (g_nat - 1)/gminus`, which is `0` ONLY because `g_nat = 1`.
   If `g_nat` were mis-stated to any normalization ≠ 1, the residual would be nonzero and the
   check would FAIL. Honest strength assessment: this is a consistency check of two stated
   assumptions — the natural-branch normalization ansatz `g_nat = 1` (line 34) and the
   branch-independence of the traction constant `C` — leveraging the stage-119 proportionality
   `g ∝ 1/T_m` taken as given. It is NOT a fully independent 4D traction re-derivation; the
   `g = C/T_m` law itself is imported, not re-derived here. It is nonetheless a genuine
   de-tautologization: the old form could never fail under any mis-statement of either input,
   whereas the new form fails if `g_nat ≠ 1` or if `C` were made branch-dependent. Moderate
   strength, correctly scoped to what this stage owns.

**Side note (informational, non-blocking):** because `expect_zero` runs
`sp.simplify(sp.expand(...))`, a `g_nat`-misstatement test would in practice catch the failure;
the residual `(g_nat-1)/gminus` does not vanish for any `g_nat != 1`. The check is therefore
live, not merely structurally non-tautological.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `T_m(-)/T_m(nat) = -20*pi/(-2*sqrt(4107 - 100*pi**2) + 37*sqrt(3))` — `C` cancelled,
  ratio is `g_nat/gminus`.
- `numeric T ratio (-) = 1.3192001633911203307`, `numeric T ratio (+) = 0.35740427386078899977`
  — unchanged from the directive's expected transcript (as predicted, since `g_nat = 1`).
- `traction ratio (-) identity = 0`, `traction ratio (+) identity = 0` — both green.
- Regression guards all green: `gminus exact form = 0`, `gplus exact form = 0`,
  `compensation quadratic at gminus = 0`, `compensation quadratic at gplus = 0`,
  `defect closed form = 0`, and `natural off compensation =
  (-12321 + 80*pi*sqrt(4107 - 100*pi**2))/(100*pi**2)` (nonzero, correctly `100*pi**2`,
  NOT `168*pi**2`).
- Total `=0` assertion lines: 7. One `expect_nonzero` line. No `168` anywhere in script or log.

**Mathematica:** exit=n/a. Stage 122 is SymPy-only; no `.wl` exists (F3 of the original
report is informational, acknowledged by the paper card, not in this directive).

**Output freshness:** confirmed. Script mtime `2026-05-29 12:54:56`, output `.txt` mtime
`2026-05-29 13:14:39` (newer), matching the exec-log header timestamp `13:14:36`. Output
was re-generated post-fix.

## Material-change assessment

`material_change`: false.

The numeric traction-ratio values are unchanged (`1.3192...` and `0.357404...`) because
`g_nat = 1` makes the new `g_nat/g±` symbolically equal the prior `1/g±`. No derived result
that downstream units consume has changed — only the derivation route of the assertion was
made independent. The compensated couplings `gminus`/`gplus`, the defect, and the closed
forms are all byte-for-byte identical to the pre-fix transcript. No downstream re-audit is
warranted on numeric grounds.

## Side observations (non-blocking)

- Line 22 banner reads "STAGE 122" correctly (the original report noted a stale "STAGE 105"
  banner; it is now correct in the current script — not load-bearing either way).
- The `100π²` vs `168π²` paper/notes discrepancy (original F1 paper_misalignment) is a
  prose/notes matter routed to the user and outside scripts-only scope; the script's
  `100π²` form remains self-consistent and reproduces all numerics. Not a script issue.

## Verdict justification

The single carried finding (tautological traction-ratio check) is genuinely resolved. The
assertion now derives `T_ratio_±` from the independent law `g = C/T_m` using a fresh
symbolic `C` (which provably cancels in the exec output) and the line-34 ansatz `g_nat = 1`,
and compares to `1/g±` only as the target — no X−X loop. The residual `(g_nat-1)/g±` is live:
it fails if `g_nat ≠ 1` or if `C` were branch-dependent, so the check exercises the
natural-branch normalization and traction-constant branch-independence rather than nothing.
This is a moderate-strength consistency check (not a full 4D traction re-derivation), which
is the correct scope for this stage. The diff is exactly the C/g derivation with no stray
edits; all regression guards remain green; script exits 0; outputs are fresh; numerics
unchanged. Verdict: verified, material_change false.
