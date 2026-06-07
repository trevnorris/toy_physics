---
unit_id: 129
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T22:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 129

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage129_mouth_boundary_layer_mathematica_audit.wl:37-44`
adds an independent derivation block placed after the `jSigma` definition (l.34)
and before the existing `Print["sigma_Pi(z) = "...]` (now l.46), exactly as the
directive specified. The block is:

- l.38 `sol = DSolve[thetaSigma*sigmaFn'[z] + v1*sigmaFn[z] == 0, sigmaFn, z];`
- l.39 `sigmaGen = sigmaFn[z] /. First[sol];`
- l.40 `cVal = C[1] /. First[Solve[Integrate[sigmaGen, {z, 0, lM}] == 1, C[1]]];`
- l.41-42 `sigmaDerived = FullSimplify[(sigmaGen /. C[1] -> cVal) /. v1 -> piM*thetaSigma/lM, ...];`
- l.43 `Print["Independently derived sigma = ", fmt[sigmaDerived]];`
- l.44 `expectZero["derived profile matches boxed sigma_Pi", sigmaDerived - sigma];`

The hard-coded `sigma` (l.33) and the three original checks (normalization l.48,
zero-flux l.52, ODE residual l.56) are untouched and now serve as confirmations.

**Assessment:**
Genuinely independent, not cosmetic. The source of `sigmaDerived` is a real
`DSolve` solve of the zero-flux ODE `Θσ' + V₁σ = 0` over an UNDEFINED function
symbol `sigmaFn` (distinct from `sigma`; `sigmaFn` appears nowhere else and is
not pre-assigned, so `sigmaFn'[z]` is a genuine non-trivial derivative). The
integration constant is fixed by a real `Solve` of the normalization integral,
not retyped. The exec log (l.10) prints the independently-derived form
`(E^(piM - (piM*z)/lM)*piM)/((-1 + E^piM)*lM)` — the boxed form rescaled by
`e^Π/e^Π` — and the comparison `sigmaDerived - sigma` reduces to exactly 0
(log l.11-12, `PASS`). This is non-tautological: a wrong ODE coefficient or
wrong normalization would leave a nonzero residual. F1 is truly resolved — the
`.wl` now derives D1 by an independent route rather than copying the `.py`'s
closed form. No collateral edits beyond the directive (the `mu` line at l.35 is
the F2 edit).

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy `scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py:25-29`
(after the residual print at l.24, per directive):
```
J_from_mu = -M*sigma*sp.diff(mu, z)
mu_link_res = sp.simplify(J_from_mu - J)
print("Onsager current from mu identity residual =", mu_link_res)
if mu_link_res != 0:
    raise AssertionError("Onsager current does not match -M*sigma*d(mu)/dz.")
```
This uses the previously dead `mu` (l.22) and the existing `J` (l.14) — no new
symbols.

Mathematica `..._mathematica_audit.wl:35` adds
`mu = thetaSigma*Log[sigma/sigmaStar] + v1*z;`, and l.57 adds
`expectZero["Onsager current from mu identity", (-mobility*sigma*D[mu, z]) - jSigma];`
placed after the existing ODE-residual check (l.56), using `jSigma` (the
un-substituted reduced current), NOT `jSub` — exactly as the directive required.

**Assessment:**
Correct and non-tautological in both engines. The check compares two structurally
different constructions: `−Mσ∂_zμ` built from the chemical-potential definition
`μ = Θ log(σ/σ_*) + V₁z` against the independently-asserted reduced current
`J = −M(Θσ' + V₁σ)`. They coincide only because `∂_zμ = Θσ'/σ + V₁`; a wrong
coefficient on the entropic (`Θ`) or potential (`V₁`) term of `μ` would make the
residual nonzero. So this genuinely exercises deliverable D4 (the μ→Onsager-current
link that motivates the whole stage) and is not a tautology. SymPy exec log l.10
shows residual `= 0` and the `raise` guard makes it a hard assertion; Mathematica
log l.23-24 shows `Onsager current from mu identity = 0` / `PASS`. The identity is
asserted BEFORE any `v1` substitution in both engines (`J`/`jSigma`), which is
correct since it holds for all `v1`. `mu` is no longer dead code. F2 resolved in
both engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- l.6 `sigma_Pi(z) = Pi*exp(Pi*(L - z)/L)/(L*(exp(Pi) - 1))` (unchanged boxed form)
- l.9 `Stationary zero-flux ODE residual = 0`
- l.10 `Onsager current from mu identity residual = 0` (new F2 check, passes)
No exceptions raised; clean exit 0.

**Mathematica:** exit=0. Notable lines:
- l.10 `Independently derived sigma = (E^(piM - (piM*z)/lM)*piM)/((-1 + E^piM)*lM)` (new F1 derivation source)
- l.12 `PASS: derived profile matches boxed sigma_Pi` (F1 independence check)
- l.24 `PASS: Onsager current from mu identity` (new F2 check)
- l.16/l.19/l.22 the three original checks still `PASS` with residual 0
- l.29 `Stage 129 Mathematica audit passed.`
All five `expectZero` checks pass; clean exit 0.

**Output freshness:** confirmed re-generated post-fix. Script mtimes are both
`2026-06-06 17:02:05`; committed `.txt` outputs are both `2026-06-06 21:47:39`
(newer). The committed `scripts/output/...sympy_audit.txt` and
`mathematica/output/...mathematica_audit.txt` contents match the captured exec
logs line-for-line (including the two new check lines). Fresh.

## Material-change assessment

`material_change`: false.

No deliverable VALUE changed. The boxed profile `σ_Π(z) = Π e^{−Πz/L}/(L(1−e^{−Π}))`,
the unit normalization `= 1`, the zero-flux current `J_σ = 0`, the stationary ODE
residual `= 0`, and the bias parameter `Π_m = V₁L/Θ_σ` are all byte-identical to
the pre-fix outputs. The edits only ADD an independent re-derivation of the
already-boxed profile (F1, which confirms the existing value) and ADD a μ→J
identity check (F2). The original `sigma` definition and the three original checks
are untouched. Downstream units that depend on stage 129's deliverables see no
changed result, so no narrow or broad re-audit is warranted on value grounds.

## Side observations (non-blocking)

- The new Mathematica derivation introduces three new top-level symbols
  (`sol`, `sigmaGen`, `cVal`, `sigmaDerived`) and `sigmaFn`/`C[1]`. `sigmaFn` is
  clean (used only in the DSolve block) and `C[1]` is Mathematica's standard
  DSolve constant; none collide with the `Clear[]` list (l.28) or downstream
  usage. No issue, noted only for completeness.
- A2/A3 (sympy) and A5/A6 (math) remain mutually redundant (the ODE-residual check
  is the zero-flux check with the `−M` factor stripped), as the original auditor
  already noted. Not a regression and out of scope for these two findings.

## Verdict justification

Both findings are fully resolved. F1's independence is genuine: the `.wl` now
sources the profile from a real `DSolve` of the zero-flux ODE plus a `Solve`-fixed
normalization constant, then `expectZero`-confirms it equals the boxed form — a
distinct derivation route, not a re-typed closed form. F2 adds a non-tautological
μ→Onsager-current identity check (`−Mσ∂_zμ − J`) in BOTH engines, exercising the
previously-unverified deliverable D4 and resurrecting the formerly-dead `mu`. Both
engines exit 0 with all checks PASS, outputs are fresh, the diff shows no collateral
edits, and no deliverable value changed (`material_change: false`). Verdict: verified.
