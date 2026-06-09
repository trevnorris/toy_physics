---
unit_id: 165
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 165

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/...stage165...wl:35-47`. The two literal coefficient sums (`eqR = dZ + 2*dcs + ...`; `eqG = ...`) were removed and replaced with seven primitive log-drift scalings (`dEll = -dcsw`; `dLam = 2*da + dEll + dLW/2 + dv`; `dKq = dZ + 2*dcs - 2*dLW`; `dKs = 2*da - drho - dEll`; `dJs = 2*da + dEll`; `dGs = dT + dJs`; `dGq = dZ - 3*dLW/2`). The two channels are then assembled and FullSimplify'd: `chanR = dKs + dKq - 2*dLam`, `chanG = dGq + dKs - dGs - dLam`, with `eqR = (chanR == 0); eqG = (chanG == 0)`. Two new `expectZero` checks (wl:44-45) compare the assembled channels against the expected coefficient vectors.

**Assessment:**
Genuinely independent, NOT a relabel. The target coefficient vectors are no longer typed anywhere as the equation definitions — they EMERGE from linear combinations of the parent-coupling (K_s, K_q, λ, J_s, g_s, g_q) primitive scalings via FullSimplify. I re-derived the assembly independently: `chanR - targetR = 0` and `chanG - targetG = 0` both hold exactly, so the §3/§4 structural route reproduces the Stage-164 fixed-r/fixed-g conditions. The two new "matches channel" `expectZero` checks are non-tautological (they would fail if any primitive scaling coefficient were mis-assembled) and both PASS in the exec log. The downstream Solve/DN-substitution/ratio/product/n=5 blocks are unchanged, and all four original `expectZero` checks still PASS.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
Mirrored numeric prefactor block added to both engines (py:99-120 via diff; wl:96-114). Each pins `gstar = 0.758035078944663` (legitimate decimal — no closed form) and derives `rstar = Sqrt[4107 - 100*Pi^2]/(10*Pi)` symbolically (SymPy: `sp.sqrt(4107 - 100*sp.pi**2)/(10*sp.pi)`), NOT a pasted decimal. The four prefactors (`Tm_pref`, `v_pref`, `ratio_pref`, `prod_pref`) are computed from the §5/§6 closed forms and compared to the notes decimals (1.2715890393387603, 1.1428896163056477, 0.8987885086678338, 1.4532859092683434) at tol 1e-12.

**Assessment:**
Non-tautological and present in BOTH engines. Each check subtracts a separately-transcribed notes decimal from an independent symbolic prefactor evaluated at exact `rstar`; a slip in either the prefactor algebra or the quoted decimal would trip the assertion. Exec logs show all four PASS in both engines with residuals ~1e-16 to ~5e-16 (well under 1e-12). `rstar` is derived from the canonical Family-1 closed form as required, so the prefactors evaluate to full precision rather than to a truncated paste.

## Exec log assessment

**SymPy:** exit=0. `Tm_pref diff = 5.8e-16`, `v_pref diff = 2.8e-16`, `ratio_pref diff = 2.8e-16`, `prod_pref diff = 5.4e-16`; all original drift identities print `= 0`.

**Mathematica:** exit=0. `PASS: eqR matches Stage-164 fixed-r channel`, `PASS: eqG matches Stage-164 fixed-g channel`, four `PASS: *_pref numeric check`, plus all original PASS lines; closes "Stage 165 Mathematica audit passed."

**Output freshness:** both outputs are newer than their scripts — wl 16:08:18 vs txt 17:12:06; py 16:09:19 vs txt 17:11:51. Fresh.

## Material-change assessment

`material_change`: false. F1 only reroutes the `.wl` derivation of eqR/eqG (same result, confirmed identical to the prior literal form). F2 adds numeric confirmations of notes-only §5-6 constants without moving any derived value. No carry-forward formula changed; no downstream unit depends on a new or altered result.

## Side observations (non-blocking)

The SymPy `expect_close` truncates `value_num` display (e.g. `0.898788508667833515...`) vs Mathematica's full-precision print; cosmetic only, both pass the 1e-12 gate. The seven primitive scalings encode notes-§3/§4 proportionalities; only their linear log-drift coefficients are load-bearing and those are exactly determined, so the lack of normalization constants is harmless.

## Verdict justification

Both findings are `resolved`. F1's re-author is genuinely independent — the Stage-164 conditions now emerge by FullSimplify from primitive parent-coupling scalings (independently re-derived to match the targets exactly), not re-typed literals, and the assembled-vs-target `expectZero` guards PASS. F2 adds non-tautological four-prefactor numeric checks to both engines with `rstar` derived from `Sqrt[4107-100*Pi^2]/(10*Pi)`, all passing at 1e-12. Both engines exit 0, outputs are fresh, no regressions in the diff, and the change is method/coverage-only (material_change false).
