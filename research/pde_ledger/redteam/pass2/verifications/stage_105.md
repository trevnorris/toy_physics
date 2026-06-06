---
unit_id: 105
batch: IV.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 105 (CHECKPOINT)

## Per-finding outcomes

### F1 — insufficient_verification (checkpoint tautology)

**Classification:** resolved

**What changed:**

SymPy (`scripts/...stage105...sympy_audit.py:49-68`): a new block constructs the
outgoing l=2 spherical-Hankel DtN operator from an EXPLICIT closed form —
`j2_out = (3/z³ − 1/z) sin z − 3 cos z / z²`, `y2_out = −((3/z³ − 1/z) cos z + 3 sin z / z²)`,
`h2_out = j2_out + i·y2_out` (these are the standard `j_2`, `y_2`; verified correct by hand),
then `Lambda_out = z·d/dz ln h2_out`, series-expands it, and asserts (can-fail)
`Lambda_2^out series fingerprint == -3 + z²/3 + z⁴/9 + i z⁵/9` (line 60). It then normalizes
`Y_out = -3/Lambda_out`, reads the imaginary z⁵ coefficient as a DERIVED quantity and asserts
it equals `1/27` (line 65), maps `z = aω/c_s` to obtain `outgoing_omega5_coeff` (line 66), and
solves `[retarded iω⁵ coeff] = outgoing_omega5_coeff` for `chi_Q`, asserting `chi_Q − 1 == 0`
(lines 68-71).

Mathematica (`mathematica/...stage105...mathematica_audit.wl:52-73`): same logical chain but
built from the NATIVE primitive `SphericalHankelH1[2, zOut]` → `lambdaOut = zOut·D[hankelOut]/hankelOut`,
series → fingerprint `expectZero` against `-3 + zOut²/3 + zOut⁴/9 + i zOut⁵/9` (line 58),
`yOut = -3/lambdaOut`, derived z⁵ coeff `== 1/27` (line 63), z→ω map (line 64), then
`Reduce[imagPart5 - outgoingOmega5Coeff == 0, chiQ]` → `chiQ == 1` (lines 66-73).

**Assessment:**

Correct and complete. The load-bearing canonical match (py:68 / wl:67) now references the
DERIVED `outgoing_omega5_coeff` / `outgoingOmega5Coeff`, not a typed literal. I grepped the full
sources: the only remaining `27*c_s**5` / `27*cSound^5` occurrences are (a) py:37/wl:37 — the
`sigma_Q^can = (9/8)/Omega_Q⁵` D1 definition check (a genuinely derived-vs-literal compare, NOT
the chi_Q pin), and (b) py:47/wl:50 — the retarded `imag ω⁵ coeff = chi_Q·a⁵/(27c_s⁵)` D2 series
check (the auditor's A4/B2, a real series exercise). The typed literal is GONE from the canonical
chi_Q match RHS as required. The fingerprint asserts are genuinely can-fail: a wrong Hankel form
or sign would break `Lambda series == -3 + z²/3 + z⁴/9 + i z⁵/9` and the derived-z⁵-coeff `== 1/27`.
`chi_Q` stays a free positive symbol up to the solve; the solve has a unique root 1 because the LHS
`chi_Q·K` is now matched against an INDEPENDENTLY-derived `K = a⁵/(27c_s⁵)` rather than a copy — a
transcription error in the true DtN coefficient would now surface as `chi_Q ≠ 1` or as a fingerprint
failure. No collateral edits: pre-existing sigma/ω²/ω⁴/imag-ω⁵ checks and the deformed branch are
untouched.

### F2 — mathematica_transliteration (subsumed by F1)

**Classification:** resolved

**What changed:**

The two engines now reach `Lambda_2^out` by VISIBLY DIFFERENT routes: SymPy via an explicit
`j_2 + i y_2` spherical-Bessel closed form built from sin/cos rationals; Mathematica via the native
`SphericalHankelH1[2, z]` special-function primitive. These are structurally independent (hand-typed
rational-trig form vs. built-in Hankel), not a line-by-line port.

**Assessment:**

Resolved. The retarded/canonical half no longer mirrors the `.py` step-for-step on a shared typed
RHS — each engine independently constructs the fingerprint and the match target. Combined with the
already-independent deformed half (SymPy `series(-3/Lam_def)` vs. Mathematica polynomial inversion of
`Lambda·Y = -3`, left UNCHANGED per directive criterion 5), the `.wl` is now independent across both
the canonical and deformed sections. The F2 requirement (subsumed by F1) is satisfied.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
`Lambda_2^out series fingerprint = 0`; `Y_2^out imag z^5 coefficient - 1/27 = 0`;
`chi_Q from exact outgoing match = 1`; all pre-existing `omega^2/omega^4/imag omega^5 coefficient = 0`
and `sigma_Q^can ... = 0` still present; deformed z²/z⁴/z⁵ coeff = 0.

**Mathematica:** exit=0. Notable lines:
`Lambda_2^out(z) from SphericalHankelH1 = -3 + zOut^2/3 + zOut^4/9 + (I/9)*zOut^5`;
`PASS: Lambda_2^out series fingerprint`; `PASS: Y_2^out imag z^5 coefficient - 1/27`;
`chi_Q condition from outgoing match (Reduce) = chiQ == 1`; `PASS: chi_Q - 1`; all four deformed
`PASS` lines and the sigma/ω²/ω⁴/imag-ω⁵ `PASS` lines retained.

**Output freshness:** confirmed. Committed `.txt` outputs (both at 2026-06-06 01:57:20) are newer
than their scripts (py 01:46:30, wl 01:46:42). The committed outputs contain the new fingerprint
lines and match the orchestrator's exec logs byte-for-byte on the relevant lines.

## Material-change assessment

`material_change`: false. Every emitted deliverable value is unchanged: `Omega_Q = 3c_s/(2a)`,
`sigma_Q^can = 4a⁵/(27c_s⁵)`, retarded series coeffs, `chi_Q = 1`, deformed `Y_2^def = 1 + z²/9 +
4z⁴/81 + iξ_Q z⁵/27`, `chi_Q = ξ_Q`. The edit adds the fingerprint derivation lines and re-grounds
the match METHOD only. Downstream consumers (097/100/106, 107-113) depend on the VALUE `chi_Q = 1`,
which is preserved. No downstream re-audit needed on value grounds.

## Side observations (non-blocking)

- py:47/wl:50 (the retarded `imag ω⁵ coeff = chi_Q·a⁵/(27c_s⁵)` D2 check) still carries the literal
  on its RHS, but this is by design: it asserts the STRUCTURE of the retarded series (coefficient is
  `chi_Q × a⁵/(27c_s⁵)`), which is a separate deliverable from the canonical pin. The auditor agreed
  (A4/B2 "anchored: yes"). Not a defect.
- The SymPy `h2_out = j2_out + i·y2_out` differs from the textbook overall phase of `h_2^(1)`
  (which carries an extra `i^{-(l+1)}` / `(-i)` factor in some conventions), but the DtN ratio
  `z·h'/h` is phase-invariant, so the fingerprint and `Y_2^out` are convention-independent — the
  series came out exactly `-3 + z²/3 + z⁴/9 + i z⁵/9` as the can-fail assert requires. Correct.

## Checkpoint higher-bar verdict

**CLEARED.** With the fix, `chi_Q = 1` is now FORCED by the actual outgoing l=2 DtN fingerprint:
the match RHS is the Hankel-derived `a⁵/(27c_s⁵)` (via a can-fail `Lambda_2^out` fingerprint assert
and a derived `Y_2^out` z⁵ coefficient = 1/27), NOT a typed copy of the answer. The pin is exercised,
not self-matched. The card's Check #3 ("outgoing l=2 DtN fingerprint against the normalized z-expansion")
is genuinely performed in both engines by visibly distinct constructions.

## Verdict justification

Both findings are `resolved`. The canonical `chi_Q = 1` pin is now derived from the exact spherical-Hankel
DtN operator in both engines via genuinely independent constructions (explicit `j_2 + i y_2` in SymPy,
native `SphericalHankelH1` in Mathematica), the fingerprint and derived-coefficient asserts are can-fail,
the typed `a⁵/(27c_s⁵)` literal is gone from the load-bearing match RHS, all pre-existing checks and emitted
values are unchanged, the deformed branch is untouched, and both engines exit 0 with fresh committed outputs.
The checkpoint higher bar is cleared. Verdict: **verified**, `material_change: false`.
