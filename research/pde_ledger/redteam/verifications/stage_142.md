---
unit_id: 142
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T18:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 142 (REMEDIATION)

This supersedes the stale 2026-05-27 verdict that PRE-DATED the authoritative
recovery review. The authoritative findings are in
`redteam/codex_reviews/stage_142.md` (verdict FINDINGS: R1 tautological,
R2 transliteration). The directive `redteam/directives/stage_142.md` was rewritten
against that review; only TWO findings were open (R1, R2). F2-kept (five
external-decimal anchors) and F5 (banners) are pre-resolved per the directive's
`## RESOLVED:` blocks. The crux: confirm the still-tautological R_q check was
relabeled and a genuinely independent anchor added (R1), and the self-series
transliteration was removed and replaced by an independent integral (R2).

## Per-finding outcomes

### F1 (R1) — tautological_check

**Classification:** resolved

**What changed:**
- SymPy `scripts/...stage142...sympy_audit.py:75-95`: the original
  `R_q(Pi_*)=1/4` check is KEPT but explicitly relabeled "solver-consistency"
  (lines 75-83), tol `1e-15`. A NEW non-tautological anchor (lines 85-95)
  evaluates `R_q` at `Pi_ext = sp.Float("1.50882951349315558300555075595", 30)`
  annotated `# Stage 131 Pi_* (independent)`, asserting `|R_q(Pi_ext) - 1/4| <= 1e-12`.
- Mathematica `mathematica/...stage142...mathematica_audit.wl:112-119`: the
  solver-consistency `expectApprox` is relabeled (line 114) at tol `10^-15`; a new
  anchor (lines 115-119) sets
  `piExt = N[Rationalize[1.50882951349315558300555075595, 0], 30]` and asserts
  `R_q(Pi_ext) = 1/4 (independent anchor)` at tol `10^-12`.

**Assessment:** Correct and genuinely non-tautological.
- (a) The SymPy `Pi_ext` literal is `1.50882951349315558300555075595`. I
  cross-checked it against Stage 131's own output
  (`scripts/output/...stage131...sympy_audit.txt:2`: `Pi_* = 1.50882951349315558300555075595`)
  — exact match; the comment cites Stage 131.
- (b) The Mathematica `piExt` uses the identical Stage-131 literal.
- (c) Both PASS: sympy log L33 `R_q(Pi_ext) - 1/4 (independent anchor) =
  1.479...E-31`; mathematica log L50-51 `R_q(Pi_ext) = 1/4 (independent anchor)
  residual = 0; PASS`.
- The independence is real because `Pi_ext` is Stage 131's value, found by a
  structurally different route (cleared-denominator FindRoot, batch-4-verified),
  and is DISTINCT from 142's own nsolve output. Sympy log L26 shows 142's own
  nsolve `Pi_*` = `1.50882951349315552747043511772`, diverging from 131's
  `...558300555075595` at digit ~16 — exactly as the anti-tautology guard warned.
  The anchor does NOT use 142's own line-81 literal `...5274704351177`, and does
  NOT re-solve `gPi=g_-`. `R_q(Pi_ext)` lands on 1/4 only because the hardcoded
  `gPi`/`r` genuinely cross `g_-` at this externally fixed bias; a `gPi` sign typo
  shifts `R_q(Pi_ext)` by O(1), not O(1e-12).
- Tolerance: solver-consistency uses `1e-15` (sympy L82) / `10^-15` (wl L114),
  NOT the over-tight `1e-20`/`10^-20`. Sympy residual is `1.945e-18` (log L32),
  inside `1e-15`. Not reverted.
- No collateral edits; matches `## Applied: F1` ("deviation: none").

### F2 (R2) — transliteration (self-series)

**Classification:** resolved

**What changed:**
- Mathematica `mathematica/...stage142...mathematica_audit.wl:60-76`: the
  `gPiSeries = Normal[Series[gPi,...]]` self-comparison (former wl:64-73) is GONE.
  Replaced by an independent projection integral:
  `sigmaPi = piM*Exp[-piM*zVar]/(1 - Exp[-piM])` (L69, Stage 129 §2 source law),
  `gPiFromSource = FullSimplify[Integrate[sigmaPi*Cos[Pi*zVar/2], {zVar,0,1},
  Assumptions -> piM > 0], ...]` (L70-73), then
  `expectZero["g_Pi closed form = integral of mouth-source law (Stage 130 §1)",
  gPiFromSource - gPi]` (L74-76).
- SymPy script untouched for R2 (correct — it never had a `gPiSeries` check).

**Assessment:** Correct and genuinely independent.
- (a) The integrand is `sigmaPi * Cos[Pi*zVar/2]`. Here `Pi` is the geometric
  constant (the bias variable is `piM`, used throughout), so the projection shape
  is `cos(pi z / 2)` — the correct first D/N derivative shape, NOT the bias
  variable. The integrand never references the hardcoded `gPi`, so a typo in
  `gPi` is not shared with the integral. `grep` confirms zero `Series[` /
  `gPiSeries` occurrences remain.
- (b) PASS with symbolic-zero residual: mathematica log L18 the integral
  evaluates to `(2*piM*(Pi + 2*E^piM*piM))/((-1 + E^piM)*(Pi^2 + 4*piM^2))`,
  exactly matching the hardcoded `gPi` (log L28); L19 reports the `expectZero`
  residual `= 0`; L20 `PASS`. This is symbolic closure, not a small-`piM`
  numeric coincidence.
- Matches `## Applied: F2` ("deviation: none").

### F2-kept (RESOLVED) — five external-decimal anchors

**Classification:** resolved

**What changed:** Intact and passing. SymPy L97-101: `g_-^{F1}` (tol 1e-25),
`Pi_*`, `S_q(Pi_*)`, `Sigma_0(Pi_*)`, `That(Pi_*)` (tol 1e-12 each); log L34-38
residuals ~1e-29. Mathematica L120-124 mirror; log L52-61 all PASS.

**Assessment:** Untouched, un-re-toleranced, all passing — the genuine numeric
anchors of the stage. Confirmed kept per directive.

### F5 (RESOLVED) — banner mismatch

**Classification:** resolved

**What changed:** No edit required; already correct. `grep` for `STAGE 125` /
`Stage 125` across both files returns zero hits. Banners read `STAGE 142`
throughout (sympy L28/L103; wl L39/L126; wl L138 "Stage 142 Mathematica audit
passed.").

**Assessment:** Correct.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L20: `R_q(g_minus)-1/4 = 0` (symbolic identity).
- L32: `R_q(Pi_*) - 1/4 (solver-consistency) = 1.945...E-18` (within 1e-15).
- L33: `R_q(Pi_ext) - 1/4 (independent anchor) = 1.479...E-31` (anchor passes).
- L34-38: five external-decimal residuals ~1e-29, all passing.

**Mathematica:** exit=0, 11 PASS / 0 FAIL. Notable lines:
- L18-20: integral `= (2*piM*(Pi + 2*E^piM*piM))/((-1 + E^piM)*(Pi^2 + 4*piM^2))`;
  `g_Pi closed form = integral of mouth-source law (Stage 130 §1) = 0`; `PASS`.
- L48-49: `R_q(Pi_*) numeric = 1/4 (solver-consistency)` residual 0 → PASS (tol 10^-15).
- L50-51: `R_q(Pi_ext) = 1/4 (independent anchor)` residual 0 → PASS.
- L52-61: five external-decimal anchors all PASS.

**Output freshness:** Confirmed regenerated post-fix. Scripts last modified
2026-05-29 16:51 (both); saved outputs `scripts/output/...sympy_audit.txt` and
`mathematica/output/...mathematica_audit.txt` both modified 2026-05-29 18:02 —
newer than the scripts and consistent with the exec-log timestamps (17:58).

## Material-change assessment

`material_change`: false.

No derived value changed. F1 relabeled a tautological solver-residual check (kept,
not removed) and added an independent anchor; F2 removed a transliterated
self-series check and replaced it with an independent integral. Canonical derived
results (`Pi_*`, `Sigma0_*`, `That_*`, the five anchors) are unchanged — the edits
only strengthened verification, altering no formula or numeric output. Downstream
units depending on 142's ledger are unaffected.

## Side observations (non-blocking)

- The SymPy script carries no independent `gPi` re-derivation (R2 was
  Mathematica-only per directive); its independence rests on F1's external anchor
  and the F2-kept decimal targets. By design, not a gap.
- The Mathematica `Pi_*` external-anchor residual is ~5.5e-17 (log L54), inside
  the 10^-12 tolerance; the recorded decimal target (`...5274704351177`, 142's
  own nsolve flavor) differs from FindRoot's higher-precision `...5830055...` at
  digit ~16 — the same benign divergence the F1 guard documents. Affects no verdict.

## Verdict justification

All four required closures hold. F1's tautological `R_q(Pi_*)=1/4` was kept but
relabeled solver-consistency, and a genuinely independent anchor evaluating `R_q`
at Stage 131's structurally-different-route `Pi_*`
(`1.50882951349315558300555075595`, confirmed against 131's own output, distinct
from 142's own nsolve `...552747...`) passes in both engines. F2's self-series
transliteration is gone, replaced by a symbolic projection integral of the
Stage-129 mouth-source law against `cos(pi z/2)` that closes to the hardcoded
`gPi` with symbolic-zero residual and shares no literal with `gPi`. The F2-kept
decimal anchors are intact and passing; F5 banners read STAGE 142; tolerances are
`1e-15`/`10^-15` (not reverted). SymPy exits 0, Mathematica exits 0 with 11 PASS /
0 FAIL, outputs freshly regenerated, and no derived value changed. Verdict:
verified.
