---
unit_id: 161
batch: IV.6
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 161

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The `.wl` (`mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl`) was
re-authored along a derivation route independent of the SymPy script. Concretely (per the diff patch
and the current file):

- **Linearization (wl:40-46):** the `eps -> 0` derivative trick (old wl:40-41
  `bWPert = ...; dBW = D[bWPert, eps] /. eps -> 0`) is replaced by a first-order series:
  `bWPertSeries = Normal[Series[bW /. {epsKappa->eps*depsK, epsGamma->eps*depsG, rc->rcStar+eps*drc}, {eps,0,1}]]`
  then `dBW = FullSimplify[Coefficient[bWPertSeries, eps]]`. The asserted target at wl:48 is unchanged:
  `expectZero["linearized slippage law", dBW - (1 + rcStar)*(depsG - depsK)/9]`.
- **d eps_kappa (wl:55-64):** the old direct-differentiation of `epsKExact` plus the
  `12 lW^2 -> Pi^2 a^2 (1+rc)` substitution and `PolynomialRemainder` choreography (old wl:50-61) is
  replaced by an explicit logarithmic derivative:
  `depsKBranch = FullSimplify[D[Log[kappaRatio], lW]*dLW + D[Log[kappaRatio], a]*da + D[Log[kappaRatio], rc]*drc]`
  with `kappaRatio = 12 lW^2/(Pi^2 a^2 (1+rc))`. The residual is now `FullSimplify[depsKBranch - depsKTarget]`
  (no `PolynomialRemainder`). Asserted target `depsKTarget = 2 dLW/lW - 2 da/a - drc/(1+rc)` is unchanged (wl:62-64).
- **d eps_gamma (wl:66-79):** old `epsGExact` direct differentiation + `gamma0Sym -> (1+rc)/9` (old wl:48-51)
  is replaced by a logarithmic-variation route:
  `depsGLog = FullSimplify[D[Log[gamma0Sym/(1+rc)], gamma0Sym]*dgamma0 + D[Log[gamma0Sym/(1+rc)], rc]*drc]`,
  then `depsGBranch = FullSimplify[depsGLog /. gamma0Sym -> (1+rc)/9]`. Asserted target at wl:76-79 unchanged.
- **difference identity (wl:80-83):** old `PolynomialRemainder`-based residual is replaced by plain
  `FullSimplify[(depsGBranch - depsKBranch) - (9 dgamma0/(1+rc) - 2 (dLW/lW - da/a))]`. Asserted at wl:83, unchanged.
- **Collapse (wl:90-105):** `stage160Prefactor = FullSimplify[9 sigmaStar/((1-sigmaStar)(1+rcStar))]` is built
  as a separate named quantity, `upsilonPi = (1+rcStar)(xiGamma-2 xiL)/9` separately, and
  `deltaQ = -stage160Prefactor*upsilonPi*dPiTan` (old inline `-9 sigmaStar upsilonPi dPiTan/((1-sigmaStar)(1+rcStar))`).
  Both `collapsed Delta_Q law` and `collapsed N_Q-1 law` asserted targets unchanged (wl:98-105).
- **Ported comment removed:** old wl:72 `(* On the branch d gamma_0 = (1+r_c)/9 * d ln gamma_0. *)` is
  deleted in the diff and is absent from the current file.

**Assessment:**
The edit correctly and fully addresses F1. The `.wl` no longer mirrors the SymPy `eps -> 0` derivative
trick (now `Series`/`Coefficient`), no longer re-types the SymPy `epsk_exact`/`epsg_exact`
direct-differentiation + branch-substitution + `PolynomialRemainder` choreography (now logarithmic
derivatives `D[Log[...], ...]` with plain `FullSimplify` residuals), and the verbatim-ported comment is
gone. This is a genuinely different algebraic route to the same identities — `Log`-derivative is an
independent linearization mechanism from SymPy's direct `diff` + closed-form substitution, and the
`Series` expansion is independent of the `diff(...).subs(eps,0)` trick. No `PolynomialRemainder` survives.

The asserted set is unchanged: 9 `expectZero` checks (exact similarity-defect decomposition; linearized
slippage law; d eps_kappa identity; d eps_gamma = d ln gamma0 - d ln(1+r_c); difference identity;
collapsed Delta_Q law; collapsed N_Q-1 law; Xi_gamma=2Xi_L => Delta_Q=0; Xi_gamma=2Xi_L => N_Q-1=0),
matching the SymPy script's nine and the auditor's A10-A18.

No new tautology introduced. The collapse construction (A15/A16) was already flagged
"construction-guaranteed" by the auditor and remains exactly that — the directive explicitly noted
making the prefactor an explicit factor "slightly hardens the cancellation" but keeps the same identity.
It does not become *more* tautological; the dependence structure is identical to before. All other
`expectZero`s remain substantive (derivatives are wrt symbols the expressions genuinely depend on).

Collateral edits: none beyond the named route changes. The carry-forward `Print` block, the similarity-
preservation substitution checks, and the printed `(1+r_F1^2)/9` and `Delta_Q in (dThat,dS)` results are
untouched in value. No deviation reported by Codex; diff confirms.

## Exec log assessment

**SymPy:** exit=0. The SymPy script was NOT edited (directive forbade it); its log confirms the
unchanged reference transcript:
- `exact similarity-defect decomposition = 0`, `linearized slippage law = 0`, `d eps_kappa identity = 0`
- `(1+r_F1^2)/9 = 0.462362334687868748`
- `Delta_Q in (dThat,dS) = ... 5.35223887169623*Xi_gamma*dThat... - 10.7044777433925*Xi_L*dThat...`

**Mathematica:** exit=0. All 9 `expectZero` print `PASS:` and the run ends `Stage 161 Mathematica audit passed.`:
- `PASS: exact similarity-defect decomposition` ... `PASS: Xi_gamma = 2 Xi_L => N_Q - 1 = 0` (9 PASS lines).
- `(1+r_F1^2)/9 = 0.46236233468786880105466350900749063734` (matches directive's required `0.46236233468786880…`).
- `Delta_Q in (dThat,dS) = ... 5.352238871696225*...xiGamma ... - 10.70447774339245*...xiL ...` and the
  `-1.16275838754222*dS*...xiGamma` term — all three required coefficients match exactly.
- The `d eps_kappa` intermediate PRINT changed FORM only: old direct-diff form
  `(24*lW*dLW - pi^2 a^2 drc - 2 pi^2 a da (1+rc))/(pi^2 a^2 (1+rc))` (still the SymPy form, log L19) is now
  the reduced log-derivative form `(-2*da)/a + (2*dLW)/lW - drc/(1 + rc)` (Mathematica log L21). These are
  algebraically equal (multiply through), and the `d eps_kappa identity = 0` check passes in both engines —
  confirmed a form-only improvement, NOT a value change. Engine agreement preserved.

**Output freshness:** confirmed. `stat` mtimes: `.wl` = 1780940984; committed `.txt` outputs (both
mathematica and sympy) = 1780941392, i.e. newer than the script. Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

Only the `.wl` derivation route changed. Every asserted identity, every expected-zero target, and every
printed numeric deliverable (`(1+r_F1^2)/9 = 0.46236233468786880…`; Delta_Q coefficients
`5.352238871696225` / `10.70447774339245` / `-1.16275838754222`) is byte-for-byte unchanged from the prior
transcript. The SymPy script and its output are untouched. No derived result that downstream units depend
on has changed, so no downstream unit is affected.

## Side observations (non-blocking)

- The `d eps_gamma` print in the new `.wl` reads `(9*dgamma0 - drc)/(1 + rc)` (log L25), identical in value
  to the SymPy `(9*d_gamma0 - dr_c)/(r_c + 1)` — consistent across engines.
- The `dPiTanExpr` literals (`0.832409471081635`, `1.16275838754222`, `6.42981496203006`) are still re-typed
  Stage-159/160 transport coefficients in both engines; the directive deliberately left these as printed
  carry-forward values (item 5: "unchanged in VALUE"), and the auditor classified them as printed/INTERNAL,
  not an independence defect. Not a regression; noted only for completeness.

## Verdict justification

The single finding F1 (mathematica_transliteration) is fully resolved: the `.wl` reaches the same nine
asserted identities by genuinely independent routes (`Series`/`Coefficient` linearization and `D[Log[...]]`
logarithmic defect variations replacing the `eps->0` trick, the re-typed closed-form differentiation, and
the `PolynomialRemainder`/`.subs(poly,0)` branch choreography), the verbatim-ported comment is gone, no new
tautology is introduced, the script exits 0 with all 9 `expectZero` printing PASS, and every printed numeric
deliverable is unchanged so engine agreement with SymPy is preserved. The `d eps_kappa` intermediate print
changed form only (reduced log-derivative form, same value, identity still passes). Outputs are fresh.
Verdict: verified, material_change: false.
