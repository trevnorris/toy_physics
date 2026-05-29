---
unit_id: 118
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 118

## Per-finding outcomes

### F1 — paper_misalignment (λ sign, resolved to direction (a))

**Classification:** resolved

**What changed:**
Both engines now carry the MINUS sign on λ consistently.

- SymPy `scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:88`:
  `lam_uniform = sp.simplify(-qstar * v0 * I_sq_uniform)` (minus present).
- SymPy:102-103 closure target:
  `expect_zero("lambda uniform closure", lam_uniform + 8*sp.sqrt(2)*qstar*v0*a**2*ell*sp.sqrt(L_W)/3)`
  — i.e. asserts `lam_uniform = -(8√2/3) q_* v0 a² ℓ √L_W`.
- Mathematica `mathematica/...mathematica_audit.wl:96`:
  `lamUniform = FullSimplify[-qStar*v0*iSqUniform, Assumptions -> $Assumptions]` (minus present).
- Mathematica:110 closure target:
  `expectZero["lambda uniform closure", lamUniform + (8*Sqrt[2]*qStar*v0*a^2*ell*Sqrt[lW])/3]`.

**Assessment:**
Correct and internally consistent. The whole point of direction (a) was self-consistency
with the script's own section IV, which *independently derives* the bilinear coefficient
(not by restating a target): SymPy:74-77 expands
`(1/2)m(rho0 + s·varrho_s)(v0 - (qstar/m)·q·A_q)²`, extracts the `s·q` coefficient via
`.coeff(s,1).coeff(q,1)`, and asserts `sq_coeff + qstar*varrho_s*v0*A_q == 0` — output
confirms `sq coefficient = -A_q*qstar*v0*varrho_s`. Mathematica:78-81 mirrors this via
`Coefficient[Coefficient[expr, sAmp, 1], qAmp, 1]`. Section V's λ now agrees with that
independently-derived minus, and the saved outputs print the resolved
`lambda (uniform-core closure) = -8*sqrt(2)*sqrt(L_W)*a**2*ell*qstar*v0/3`
(SymPy) and `(-8*Sqrt[2]*a^2*ell*Sqrt[lW]*qStar*v0)/3` (Mathematica). The MINUS is now
genuinely carried through to the printed closure value, not merely flipped inside an
assertion target. Directive records `## Applied: F1` with deviation: none. The directive
and original report both record that the notes box the minus in three places
(notes 169-171/176-180/196-197) and that downstream stage 123 already consumes UN-squared
λ with the minus — consistency is therefore restored rather than newly introduced.
No collateral edits beyond the F1/F2 surface.

### F2 — insufficient_verification (3 new asserts)

**Classification:** resolved

**What changed:**
Three new `expect_zero` / `expectZero` lines present and passing in both engines:
- `K_q closed form` (SymPy:97, Math:105)
- `g_s closed form` (SymPy:100, Math:108)
- `lambda from bilinear` (SymPy:104, Math:111)

**Assessment (per-assertion non-triviality, adversarial):**

- **`g_s closed form`** — GENUINE. `J_s` is built from the section-I tanh moment:
  SymPy:84 `J_s = sp.simplify(4*sp.pi*a**2*ell*I_f)`, and `I_f` is the independently
  computed integral asserted = 1/3 in section I (SymPy:30/34). `g_s = Tm * J_s` (line 85).
  The assertion `g_s - Tm*4*pi*a²*ell/3` (line 100) is therefore NOT X−X: it verifies that
  the independently-derived `I_f = 1/3` flows through into `g_s`'s boxed closed form. A wrong
  upstream `I_f` would make this fail. Same structure in Mathematica (jS from iF, line 92).
  Closes the original "g_s = Tm·J_s by construction" gap.

- **`K_q closed form`** — present and passing, but BY CONSTRUCTION (X−X). `K_q` is *defined*
  at SymPy:82 as exactly `(Zq/mu0)*(pi²·c_s²/(4·L_W²))`, and the assertion target (line 97) is
  the literal same expression, so the residual is identically zero regardless of any upstream
  value. The original auditor had suggested the stronger anchor `K_q - (Zq/mu0)*chi_grad`
  (which would chain through the independently-asserted D/N stiffness moment); the directive
  instead specified the literal closed form, and Codex applied the directive as written. This
  is a faithful application of the directive, but the assertion adds no independent
  verification power beyond restating the definition. I flag this as a non-blocking weakness,
  not a rework item: the directive (post Claude+Codex consult) explicitly specified this exact
  target, and the underlying D/N stiffness `chi_grad = pi²/(4·L_W²)` IS independently asserted
  upstream (`D/N stiffness check`, line 50), so the K_q boxed form is still indirectly anchored
  in the script. No deviation from directive.

- **`lambda from bilinear`** (`lam_uniform + qstar*v0*J_s*I_q`, SymPy:104 / Math:111) —
  MODEST non-triviality. `I_sq_uniform` is *defined* as `J_s*I_q` (SymPy:87), and
  `lam_uniform = -qstar*v0*I_sq_uniform`, so the residual reduces to
  `-qstar*v0*(J_s*I_q) + qstar*v0*J_s*I_q = 0` purely from the definition of `I_sq_uniform`.
  It does NOT independently verify that `I_sq_uniform` equals the 4D overlap integral
  ∫d⁴X varrho_s A_q — that factorization is asserted by construction, not re-derived. What it
  DOES verify is the SIGN coupling between the bilinear form and the uniform closure (that the
  minus is carried as `λ = -q_* v0 J_s I_q`), which is exactly the F1 concern, and it confirms
  `J_s·I_q` simplifies to the same `I_sq_uniform` SymPy/Mathematica each compute. Judgment: it
  is a sign/factorization-consistency check, weaker than a fully independent overlap derivation
  but not vacuous — it pins the F1-resolved minus to the bilinear chain in both engines.

Net: all three deliverables that were "by construction / missing" in the original report now
have explicit assertions; `g_s` is genuinely anchored, `lambda from bilinear` is a meaningful
sign-consistency check, `K_q` is faithful-to-directive but by construction. F2's stated goal
(make the verification surface match the claim surface) is met. Directive records
`## Applied: F2`, deviation: none.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `sq coefficient = -A_q*qstar*v0*varrho_s` then `bilinear sq coefficient = 0`
- `lambda (uniform-core closure) = -8*sqrt(2)*sqrt(L_W)*a**2*ell*qstar*v0/3`
- `lambda uniform closure = 0`, `lambda from bilinear = 0`, `K_q closed form = 0`,
  `g_s closed form = 0`, plus `All Stage 118 symbolic checks passed.`

PASS-line count (SymPy uses `<name> = 0` convention; counting load-bearing `expect_zero`
checks that print `= 0`): 12 zero-residual checks —
D/N norm, D/N stiffness, K_s closed form, healing-lock K_s, bilinear sq coefficient,
K_q closed form, g_q closed form, J_s closed form, g_s closed form, I_q closed form,
lambda uniform closure, lambda from bilinear (plus the section-I `I_f/I_g` guard which
raises on mismatch). All print `= 0`.

**Mathematica:** exit=0. Notable lines:
- `sq coefficient = -(aQ*qStar*v0*varrhoS)`, `PASS: bilinear sq coefficient`
- `lambda (uniform-core closure) = (-8*Sqrt[2]*a^2*ell*Sqrt[lW]*qStar*v0)/3`
- `PASS: lambda uniform closure`, `PASS: lambda from bilinear`, `PASS: K_q closed form`,
  `PASS: g_s closed form`, and `Stage 118 Mathematica audit passed.`

PASS-line count (Mathematica `PASS:` lines): 14 —
I_f-1/3, I_g-4/15, D/N norm check, D/N stiffness check, K_s closed form, healing-lock K_s,
bilinear sq coefficient, K_q closed form, g_q closed form, J_s closed form, g_s closed form,
I_q closed form, lambda uniform closure, lambda from bilinear. No FAIL lines; `Exit[0]`.

Engines use genuinely different mechanics for the load-bearing steps (SymPy `.coeff().coeff()`
vs Mathematica `Coefficient[Coefficient[...]]`; SymPy `sp.integrate` vs Mathematica
`Integrate`; independent symbolic kernels). Both reach identical final values; the previously
shared sign error is now corrected in both. Adversarial check "could a wrong upstream value
still pass?": for `K_s closed form`, `healing-lock K_s`, `bilinear sq coefficient`, `g_s closed
form`, and the D/N/tanh moments — no, these chain through independently computed integrals/
coefficients. For `K_q closed form` and `lambda from bilinear` — these are by-construction /
factorization-restatement (assessed above), but neither is a regression and both are faithful
to the directive.

**Output freshness:** confirmed. mtimes —
SymPy script 1779921289 < output 1780080483 (output newer);
Mathematica script 1779921298 < output 1780080483 (output newer). Both saved `.txt` outputs
were re-generated post-fix and contain the resolved minus-sign λ value and all new assertions.
Note: the committed `stage_118_diff.patch` is 0 bytes — EXPECTED per the remediation context
(the edits were committed earlier in the tainted IV.3 pass; this remediation reconciles +
records them). The edits are verified present in the current script bodies and in the saved
outputs, which is the substantive check.

## Material-change assessment

`material_change`: false.

This is a sign-consistency repair confined to stage 118's section V λ closure. Per the
directive and original report, downstream stage 123 already consumed the UN-squared λ with
the MINUS sign (stage123 sympy:28,32), so making 118 use the minus brings 118 *into* agreement
with the already-existing downstream convention — it does not change any value a downstream
unit currently depends on. Elsewhere λ enters only squared (Schur-complement `K_s K_q + λ²`
denominators are sign-invariant); the lone signed appearance is the numerator cross term, which
the downstream stages already treat consistently with the minus. The F2 additions are new
assertions only (no derived value changed). Hence no carried value changes downstream. (The
orchestrator may still mark units > 118 `upstream_stale: true` as routine bookkeeping, but I
flag no specific re-audit need.)

## Side observations (non-blocking)

1. `K_q closed form` (SymPy:97 / Math:105) is an X−X by-construction assertion: the target
   restates the definition at SymPy:82 literally. It adds no independent verification power
   (the auditor's alternative `K_q - (Zq/mu0)*chi_grad` would have chained through the
   independently-asserted D/N stiffness). Faithful to the directive as written; flagged for
   awareness only, not a rework trigger. The D/N stiffness moment underlying K_q is
   independently asserted upstream (line 50), so K_q's boxed form is still indirectly anchored.

2. `lambda from bilinear` verifies the sign/factorization relationship `λ = -q_* v0 J_s I_q`
   but treats `I_sq_uniform = J_s·I_q` as a definition rather than independently verifying the
   4D overlap split. Meaningful as a sign check (the F1 concern) but not a full overlap
   re-derivation.

3. Both engines remain a near-line-by-line transliteration of the same algebraic path (as the
   original auditor noted under "Independent-derivation check"); the symbolic work is still done
   by each kernel independently, so the cross-check retains value. Not a finding; carried over
   from the original report's judgment.

## Verdict justification

F1 is fully resolved to direction (a): both SymPy and Mathematica now define λ with the MINUS
sign and assert the minus closure, matching the script's own independently-derived section-IV
bilinear coefficient `-q_* v0 varrho_s A_q` and the boxed notes value; the saved outputs print
the corrected `-(8√2/3)…` closure. F2 is resolved: the three previously-missing/by-construction
deliverables now have explicit assertions — `g_s closed form` chains genuinely through the
independent `I_f`, `lambda from bilinear` is a meaningful minus-sign/factorization-consistency
check, and `K_q closed form` is faithful to the directive though by construction (flagged
non-blocking). Both engines exit 0 with all checks passing (SymPy 12 zero-residual checks,
Mathematica 14 PASS lines, no FAIL), outputs are fresh, and the 0-byte diff is the expected
tainted-pass artifact. No regressions. This is a sign-consistency repair that aligns 118 with
the already-minus downstream stage 123, so `material_change: false`. Verdict: verified.
