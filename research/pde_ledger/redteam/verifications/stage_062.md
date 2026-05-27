---
unit_id: 062
batch: III.3
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T19:30:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 062

## Per-finding outcomes

### F1 — paper_misalignment / script_missing_paper_claim (second equality + Cauchy-Schwarz bound)

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:88-100`: added the second-equality assertion and the Cauchy-Schwarz parameterization check. Specifically:
  - line 89: `C_sp_sq = Osp**2 / (Nss * Npp)`
  - line 90: `G_micro_factored = (rho_star * g_phi**2 * Npp / (m * cs_star_sq * KX)) * C_sp_sq`
  - line 91: `expect_zero("Second equality of boxed G_micro: closed vs factored form", G_micro_closed - G_micro_factored)`
  - line 94-97: introduce `theta` symbol, substitute `Osp -> cos(theta) * sqrt(Nss * Npp)` into `C_sp_sq`, then `expect_zero("C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta)", C_sp_sq_cos_simplified - sp.cos(theta)**2)`.
  - line 99 is a documentation tautology (`assert ... if False else True`) which is harmless but worth noting as inert.
- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:93-103`: mirrors the SymPy structure:
  - line 94: `cSpSq = oSP^2 / (nSS * nPP);`
  - line 95: `gMicroFactored = (rhoStar * gPhi^2 * nPP / (m * csStarSq * kX)) * cSpSq;`
  - line 96: `expectZero["Second equality of boxed G_micro: closed vs factored form", gClosed - gMicroFactored]`
  - line 99-102: `theta = Symbol["theta"]; cSpSqCos = cSpSq /. oSP -> Cos[theta] * Sqrt[nSS * nPP]; expectZero["C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta)", FullSimplify[cSpSqCos - Cos[theta]^2, ...]]`

**Assessment:**

The two new assertions per engine genuinely exercise the second equality and the Cauchy bound:

1. **Second equality non-tautology.** `G_micro_closed = rho_* g_phi^2 Osp^2/(m cs_star^2 KX Nss)` has no `Npp` and `G_micro_factored = (rho_* g_phi^2 Npp/(m cs_star^2 KX)) * (Osp^2/(Nss Npp))` has an explicit `Npp` in both numerator and denominator. The residual is zero only after `Npp` cancellation against the `C_sp_sq` denominator. A typo such as `Npp -> Npp**2` in either the prefactor or the coherence factor would break the assertion. The check is built from the same fundamental quantities (`Osp`, `Nss`, `Npp`, plus the per-stage constants) but follows two different algebraic groupings — exactly what the paper's second equality is meant to claim.

2. **Cauchy parameterization non-tautology.** Substituting `Osp -> cos(theta) sqrt(Nss Npp)` makes `C_sp_sq` reduce to `cos(theta)^2 (Nss Npp) / (Nss Npp) = cos(theta)^2` only after the `Nss Npp` factor cancels — a non-trivial simplification step. Both engines confirm via FullSimplify/simplify with `theta` real that the residual is exactly zero. This is the standard symbolic encoding of "0 <= C_sp_sq <= 1" since `cos(theta)^2` is manifestly in `[0, 1]`. The auditor's resolution direction (a) named this exact parameterization.

3. **Side note (inert line).** Sympy line 99 `assert sp.simplify(sp.cos(theta)**2 - 0) >= 0 if False else True` is effectively `assert True` — the auditor's directive called for this to be "documentation pinned" and Codex implemented it as such. It contributes no verification, but it does not falsely claim to. It's a harmless documentation marker. Not a regression.

4. **Engine-cross check on the F1 additions.** SymPy log line 24-25 shows `Second equality of boxed G_micro: closed vs factored form = 0` and `C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta) = 0`. Mathematica log lines 35-38 confirm both `PASS:` for the same two assertions. Both engines independently arrive at the same conclusion.

### F2 — paper_misalignment / notes_contradicts_script (sigma-phi coupling sign)

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:71-75`: `S_parent` cross-coupling sign flipped from `+ Lambda_phi * sigma * phi` to `- Lambda_phi * sigma * phi` (diff line 58-59).
- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:72`: `sParent` cross-coupling sign flipped from `+ lambdaPhi*sigma*phi` to `- lambdaPhi*sigma*phi` (diff line 9 -> 17).

**Assessment:**

The on-shell `sigma_star(phi)` flips sign as predicted by the auditor:

- SymPy log line 20: `sigma_star = O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)` — now **positive**, matching the notes' `F_red` with the `-Lambda sigma phi` convention.
- Mathematica log line 24: `sigmaStar = (gPhi*oSP*phi*rhoStar)/(csStarSq*m*nSS)` — same, **positive**.

The gain identity is invariant under this sign flip (it enters as `Lambda^2`), so `G_micro from action = O_sp**2*g_phi**2*rho_star/(K_X*N_ss*cs_star_sq*m)` is unchanged from the previous iteration. Both `G_micro from parent action vs closed form = 0` and the Mathematica equivalent still pass. The `S_eff_phi` form (`K_X*phi**2/2 - O_sp**2*g_phi**2*phi**2*rho_star/(2*N_ss*cs_star_sq*m)`) is also unchanged (also a squared dependence on Lambda).

The note in the directive was `direction: F3 (susceptibility route)` in the apply log, but the diff clearly shows the sigma-phi sign also flipped in both engines — F2 was applied alongside F3. This matches the user's stated context that F2 (the sigma-phi coupling flip to `-Lambda_phi sigma phi`) was applied. Resolution direction (a) was implemented: notes' sign is canonical, scripts changed.

### F3 — mathematica_transliteration (susceptibility route)

**Classification:** resolved

**What changed:**

- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:66-69`: inserted the susceptibility-route block before the action-coefficient block.
  - line 68: `chiSigmaEff = 1/thetaSigma;`
  - line 69: `gainViaSusceptibility = FullSimplify[chiSigmaEff * lambdaPhi^2 / kX, Assumptions -> $Assumptions];`
- lines 83-84: added `Print` statements for `chi_sigma^(eff)` and `G_micro via susceptibility route`.
- lines 88-89: added two new `expectZero` calls:
  - `expectZero["G_micro via susceptibility route vs closed form", gainViaSusceptibility - gClosed]`
  - `expectZero["G_micro: action route equals susceptibility route", gainFromAction - gainViaSusceptibility]`

**Assessment:**

The `gainViaSusceptibility` route is genuinely algebraically distinct from `gainFromAction`:

- `gainViaSusceptibility = (1/thetaSigma) * lambdaPhi^2 / kX`. This uses only arithmetic on the predefined `thetaSigma = (m csStarSq/rhoStar) nSS` and `lambdaPhi = gPhi oSP`. **No `Solve`, no `Coefficient`, no `Series`.**
- `gainFromAction = (kX - 2 * Coefficient[sEff, phi, 2]) / kX`, where `sEff` itself comes from `sParent /. sigma -> sigmaStar` with `sigmaStar = sigma /. First[Solve[D[sParent, sigma] == 0, sigma]]`. This chain goes through `Solve` + on-shell substitution + `Coefficient` extraction.

So `gainViaSusceptibility - gainFromAction = 0` is a non-trivial cross-check: the two computations only share the four atomic symbols (`thetaSigma`, `lambdaPhi`, `kX`, `gClosed`), not any intermediate step. A sign flip on `lambdaPhi`'s contribution in either route would be caught (in the susceptibility route, the gain is `lambdaPhi^2` so sign-insensitive on that; in the action route, the sign enters via `sigmaStar`'s sign, and the resulting `Coefficient[sEff, phi, 2]` still uses `lambdaPhi^2` — so this particular sign vulnerability is not caught by this cross-check, but a wrong factor like `lambdaPhi -> lambdaPhi^2` would be caught immediately).

Both new assertions pass:
- Mathematica log line 22: `chi_sigma^(eff) = rhoStar/(csStarSq*m*nSS)` — clean reciprocal.
- log line 23: `G_micro via susceptibility route = (gPhi^2*oSP^2*rhoStar)/(csStarSq*kX*m*nSS)` — matches `gClosed` exactly.
- log line 27-28: `G_micro via susceptibility route vs closed form = 0` / `PASS`
- log line 29-30: `G_micro: action route equals susceptibility route = 0` / `PASS`

This means the Mathematica engine now performs a genuinely independent re-derivation, satisfying the framework's second-engine policy.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- log line 23: `G_micro from parent action vs closed form = 0`
- log line 24: `Second equality of boxed G_micro: closed vs factored form = 0`
- log line 25: `C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta) = 0`
- log line 20: `sigma_star = O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)` — positive sign (matches notes' F_red after F2 fix)
- log line 29: `kappa solved from Xi_micro = Xi_target: K_X*L**2/T_X`

**Mathematica:** exit=0. Notable lines:
- log line 27-28: `G_micro via susceptibility route vs closed form = 0` / `PASS`
- log line 29-30: `G_micro: action route equals susceptibility route = 0` / `PASS`
- log line 31-32: `Mathematica two-route consistency = 0` / `PASS` (the prior internal-consistency check, preserved)
- log line 33-34: `gMicro from parent action vs closed form = 0` / `PASS`
- log line 35-36: `Second equality of boxed G_micro: closed vs factored form = 0` / `PASS`
- log line 37-38: `C_{sigma phi}^2 via Cauchy parameterization equals cos^2(theta) = 0` / `PASS`
- log line 24: `sigmaStar = (gPhi*oSP*phi*rhoStar)/(csStarSq*m*nSS)` — positive sign (matches the notes after F2 fix)

**Output freshness:** confirmed. Script mtimes are 1779820790 (sympy), 1779820805 (wl); output mtimes are 1779820820 (sympy), 1779820826 (wl). Outputs were re-generated after the edits.

## Material-change assessment

`material_change`: false.

Rationale: the boxed gain magnitude `G_micro = rho_* g_phi^2 Osp^2 / (m cs_star^2 KX Nss)` is unchanged by all three fixes. F1 adds verification of an algebraically equivalent rewriting (no new value computed). F2 changes a sign convention internally and flips `sigma_star`'s sign, but `G_micro` is quadratic in Lambda and so invariant. F3 adds a redundant route that reaches the same value. Downstream consumers (063 amplitude/coherence threshold, 064 equilibrium source profile) that import `G_micro` or `Xi_micro = kappa G_micro` see no change. The only downstream-visible artifact would be a consumer that imports `sigma_star(phi)` with its sign (rather than `|sigma_star|`); the auditor flagged this as the F2 risk, and the fix now aligns the script with the notes — if any downstream stage had silently used the old (incorrect) sign, that stage would now disagree. Looking at the audit unit ordering (063, 064 are downstream), neither imports `sigma_star(phi)` directly per the auditor's own statement that F2 was "invariant for the gain magnitude" and the prefix "invisible to the script"; downstream stages key off `G_micro`/`Xi_micro`. So material_change remains false.

## Side observations (non-blocking)

- The SymPy line 99 `assert sp.simplify(sp.cos(theta)**2 - 0) >= 0 if False else True` is a no-op (effectively `assert True`). It's labeled "tautology pinned to documentation" so the intent is clear. Not a finding, just noting it for completeness.
- The directive's `## Apply log` block (lines 165-171 of the directive) lists `direction: F3 (susceptibility route)` but the diff clearly applied F2 (the sigma-phi sign flip) and F1 (second equality + Cauchy bound) in addition. This is a metadata bookkeeping inconsistency in the apply log, not a verification problem — the actual edits all match the auditor's resolution directions. The apply log should arguably be augmented to list `direction: F1 (a), F2 (a), F3` to fully reflect what was applied, but this is cosmetic.

## Verdict justification

All three findings are `resolved`. F1 and F2 were paper_misalignment findings the orchestrator was holding for user resolution; the user chose direction (a) for both (script-side fix), and Codex applied them correctly. F3 was a script-side rewrite, and Codex added the susceptibility-route block that is structurally and algebraically distinct from the action-coefficient route. Both engines exit 0; all new assertions (`Second equality of boxed G_micro`, `C_{sigma phi}^2 via Cauchy parameterization`, `G_micro via susceptibility route vs closed form`, `G_micro: action route equals susceptibility route`) are non-tautological and pass. The on-shell `sigma_star(phi)` sign matches the notes' `F_red` after the F2 flip. No regressions in the diff: the only edits are exactly the lines named in the directive and the F1/F2 additions, with no collateral changes to the EOS block, the inconsistency probe, or the `Xi_micro` block. `material_change: false` because all numeric/symbolic outputs reaching downstream stages are unchanged.
