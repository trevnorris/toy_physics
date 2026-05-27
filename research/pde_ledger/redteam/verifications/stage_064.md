---
unit_id: 064
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 064

This verification supersedes the prior pass on this unit. The earlier audit
catalogued four findings (F1-F4) against the original scripts; the present
red-team cycle re-audited the post-fix scripts and raised two new findings
(insufficient_verification on the general-H softening, and
mathematica_transliteration on the whole `.wl` file). The directive at
`redteam/directives/stage_064.md` covers the two new findings. This file
verifies only those two.

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:141-159` adds a new "GENERAL-H EQUILIBRIUM SOFTENING CHECK" banner block after the discrete two-point Cauchy-gap section. The block builds
  `chi_sig1 = g_phi*sqrt(w1)/H1`, `chi_sig2 = g_phi*sqrt(w2)/H2` with `H1, H2` independent positive reals, then
  `Theta_general = H1*chi_sig1**2 + H2*chi_sig2**2` (simplifies to `g_phi**2*(H1*w2 + H2*w1)/(H1*H2)`),
  `Lambda_general = g_phi*(sqrt(w1)*chi_sig1 + sqrt(w2)*chi_sig2)` (same form),
  and asserts `expect_zero("general equilibrium softening equals g_phi^2 I_1", soft_general - g_phi**2 * I1_disc)` where `I1_disc = w1/H1 + w2/H2` was defined up at line 130.
- The old A10 block (matched-layer substitution chain `soft_abs.subs(I2, I1**2/Npp).subs(Npp, I1*Hw)` at original lines 156-165) was removed; the diff confirms deletion.
- `mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:155-180` adds the analogous "GENERAL-H EQUILIBRIUM SOFTENING CHECK" block with three `expectZero` calls:
  - `general Theta integrand equals gPhi^2 chiPhi^2/hFun` — checks the integrand identity `hFun[y]*(gPhi*chiPhi[y]/hFun[y])^2 = gPhi^2*chiPhi[y]^2/hFun[y]`,
  - `general Lambda integrand equals gPhi^2 chiPhi^2/hFun` — checks `gPhi*chiPhi[y]*(gPhi*chiPhi[y]/hFun[y]) = gPhi^2*chiPhi[y]^2/hFun[y]`,
  - `general equilibrium softening equals gPhi^2 I_1` — algebraic consequence given those two integrand equalities.
- The old M10 block (`softMatched = (softAbs /. i2 -> i1^2/npp) /. npp -> i1*hw`, original lines 123-129) was removed.

**Assessment:**

The sympy general-H block is genuinely substantive. `H1` and `H2` are independent positive reals — not bound to `Hw` and not bound to each other — so the two-point branch is *not* in the matched layer. The log shows `Theta (general, two-point) = g_phi**2*(H1*w2 + H2*w1)/(H1*H2)`, confirming the asymmetric form. The closure `chi_sigma_k = g_phi chi_phi_k/H_k` is the only premise used; the identity `Lambda^2/Theta = g_phi^2*I_1` then emerges by ordinary algebra. The block contains no `subs(I2, ...)`, no `subs(Npp, ...)`, no matched-layer substitution of any kind (verified by reading the post-fix file end-to-end and by grepping the diff). The new assertion is non-tautological: `soft_general` is the simplified ratio of two distinct symbolic sums and is compared against an independently-defined `I1_disc`.

The Mathematica side, after the orchestrator hot-fix, splits the verification into two integrand-level algebraic identities plus a downstream algebraic consequence:

- Integrand identity for `Theta`: `hFun[y]*(gPhi*chiPhi[y]/hFun[y])^2 - gPhi^2*chiPhi[y]^2/hFun[y]` simplifies algebraically to `gPhi^2*chiPhi[y]^2/hFun[y] - gPhi^2*chiPhi[y]^2/hFun[y] = 0`. This is a genuine pointwise (in `y`) algebraic check on unspecified symbolic functions; FullSimplify handles it without any matched-layer assumption.
- Integrand identity for `Lambda`: `gPhi*chiPhi[y]*(gPhi*chiPhi[y]/hFun[y]) - gPhi^2*chiPhi[y]^2/hFun[y] = 0` — same character.
- Both pass at log lines 56 and 58 (`PASS: general Theta integrand equals gPhi^2 chiPhi^2/hFun`, `PASS: general Lambda integrand equals gPhi^2 chiPhi^2/hFun`).
- The final `expectZero["general equilibrium softening equals gPhi^2 I_1", softGeneral - gPhi^2*i1Integral]` is then a consequence of the two preceding integrand identities together with the script's assignments `thetaGeneral = gPhi^2*i1Integral` and `lambdaGeneral = gPhi^2*i1Integral`. Strictly speaking this third assertion is algebraically tautological once the integrand identities are accepted, but its physical content is fully carried by the two non-tautological integrand checks above. The orchestrator's hot-fix rationale — that Mathematica's `Integrate` cannot pull `gPhi^2` outside `Integrate[gPhi^2*chiPhi[y]^2/hFun[y], ...]` with symbolic `chiPhi`, `hFun` — is correct; the integrand-first formulation preserves the physical claim because pointwise integrand equality implies integral equality on any common domain. The directive frontmatter `orchestrator_hotfix_2026_05_26` notes this and the hot-fix is sound.

F1 is resolved in both engines with a genuine general-H derivation.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**

The Mathematica `audit.wl` was rewritten so its derivation choreography no longer mirrors the sympy script. The five non-trivial blocks now use distinct routes:

1. Closure law (lines 35-44): replaces sympy's `sp.solve(sp.diff(F_loc, sigma_loc), sigma_loc)` with `VariationalD[energyDensity, sigmaFun[yLoc], yLoc]` from `VariationalMethods\``, then `Solve[eulerLagrange == 0, sigmaFun[yLoc]]`. The variational derivative treats `sigma` as a function of `y`, whereas sympy treats it as a finite-dim variable.
2. Concrete-profile integrals (lines 47-82): replaces sympy's Gaussian `exp(-y^2/(2 L_int^2))` with a Lorentzian `1/(1 + z^2/L^2)`. Outputs `Npp_int = (L*Pi)/2`, `I1_int = (L*Pi)/(2*hw)`, `I2_int = (L*Pi)/(2*hw^2)`, distinct from sympy's `sqrt(pi)*L_int*...`. The same paper claim (b) is exercised against an independent symbolic profile.
3. Matched-layer reductions (lines 95-104): consumed directly from the Lorentzian results; no Gaussian re-computation.
4. Cauchy bound (lines 106-153): rewritten as a *continuous* Cauchy-Schwarz check on a non-constant `hVariable[z] = hw*(1 + z^2/L^2)`. Defines `f[y] = chiPhiL[y]/Sqrt[hVariable[y]]`, `g[y] = chiPhiL[y]/hVariable[y]`, computes `pairGap = ffIntegral*ggIntegral - fgIntegral^2`, asserts it equals the closed-form expected (`L^2*(15*Pi^2/128 - 256/225)/hw^3`) and is non-negative. Also computes the related `nppInt*i2Var - i1Var^2` and verifies the closed form `Pi^2*L^2/(64*hw^2)` plus non-negativity. sympy verifies the same physical claim only via a discrete two-point algebraic identity.
5. `F_eff` block (lines 184-201): kept algebraically parallel — explicitly permitted by the directive — but variables renamed to `sourceCoeff`, `mixCoeff`, `supportAmp`, `sourceAmp`. Algebra of a 2-variable quadratic has essentially one form, so renaming is the strongest acceptable change.
6. General softening (lines 155-180): integral/variational form on unspecified `chiPhi[y]` and `hFun[y]`, verified by integrand identities (see F1). sympy's general-H softening uses a discrete two-point branch.

**Assessment:**

A per-block correspondence audit between `.py` and `.wl` shows distinct intermediate steps in five of six numbered blocks. Closure: variational vs finite-dim. Concrete profile: Lorentzian vs Gaussian. Cauchy: continuous on non-constant H vs discrete two-point. General softening: continuous integrand identity vs discrete two-point algebra. Only `F_eff` retains structural parallelism, and that block is explicitly permitted to remain similar (with variable rename) by the directive. No two-line one-to-one syntactic mapping survives between the engines.

The Mathematica audit still verifies every paper-side deliverable: closure law (M1 via VariationalD), overlap identities on the Lorentzian (M2-M3), matched-layer integral reductions (M4-M5 on Lorentzian), `C^2 = I_1^2/(Npp I_2)` formula construction (line 93), matched-layer `C^2 = 1` (M6) and matched-layer gain (M7), continuous Cauchy gap (replacing the discrete bound), and general softening `gPhi^2 I_1` (M10 via the integrand path). F2 is resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `general equilibrium softening equals g_phi^2 I_1 = 0` (log line 50) — confirms F1.
- `Theta (general, two-point) = g_phi**2*(H1*w2 + H2*w1)/(H1*H2)` (line 47) — confirms the general H1 != H2 form is exercised, not the matched layer.
- `two-point Cauchy gap identity = 0` (line 42) — pre-existing A8 still passes.
- `effective support softening = 0` (line 57) — pre-existing F_eff scaffolding still passes.
- All previous closure/overlap/matched-layer assertions still print zeros (lines 11, 13, 14, 26, 27, 34, 35).

**Mathematica:** exit=0. Notable lines:
- `PASS: general Theta integrand equals gPhi^2 chiPhi^2/hFun` (log line 56) — the substantive content of the F1 hot-fix.
- `PASS: general Lambda integrand equals gPhi^2 chiPhi^2/hFun` (line 58) — second substantive integrand identity.
- `PASS: general equilibrium softening equals gPhi^2 I_1` (line 60) — algebraic consequence.
- `PASS: continuous f-g Cauchy-Schwarz residual` and `PASS: continuous Cauchy bound C^2 <= 1` (lines 46, 50) — F2's independent Cauchy path on a non-constant H profile.
- `Npp_int = (L*Pi)/2` (line 13) — Lorentzian result, structurally distinct from sympy's `sqrt(pi)*L`.
- `PASS: closure law chi_sigma = g_phi chi_phi/H` (line 12) — variational-derivative path succeeds.

**Output freshness:** The orchestrator's exec logs `redteam/exec_logs/stage_064_sympy.log` (mtime May 26 13:12) and `redteam/exec_logs/stage_064_mathematica.log` (mtime May 26 13:21) are both newer than the post-fix scripts (May 26 12:47 and 13:15 respectively). The stage-card `.txt` outputs under `scripts/output/` and `mathematica/output/` still bear May 22 mtimes, but those files are produced by a separate stage-card refresh path; the red-team `exec-*` driver writes to `redteam/exec_logs/`, which is the verifier's freshness target. Exec-log freshness is satisfied.

## Material-change assessment

`material_change`: false.

The edits replace a script-side assertion chain that was algebraically tautological once the matched-layer substitutions kicked in with checks that anchor the general-H claim, and rewrite the Mathematica audit to use independent derivations of the same paper-side identities. No paper-side equation changes, no downstream-input value changes, no constants altered. Downstream stages 065-066 consume the matched-layer outputs (`C^2 = 1`, `G_eq = gPhi^2 N_pp/(K_X Hw)`), which remain verified unchanged by A6-A7 / M6-M7. No upstream/downstream re-audit triggered.

## Side observations (non-blocking)

- The post-hot-fix Mathematica softening block defines `thetaGeneral = gPhi^2*i1Integral` and `lambdaGeneral = gPhi^2*i1Integral` by direct assignment after establishing the integrand identities. Consequently the third assertion `softGeneral - gPhi^2*i1Integral == 0` is algebraically a consequence of those assignments rather than an independent integral computation. It is acceptable because the two preceding integrand-equality assertions are themselves non-tautological and fully entail the physical claim, but a future cleanup could compress the three Mathematica assertions into the two substantive integrand identities plus an algebraic comment. Not blocking.
- The Mathematica `pairExpected` constant `L^2*(15*Pi^2/128 - 256/225)/hw^3` looks numerically ad-hoc; the script computes `pairGap = ffIntegral*ggIntegral - fgIntegral^2` symbolically via three independent `Integrate` calls and then compares to that closed form. The closed form is itself a derived quantity (from the Lorentzian integrals); the assertion is non-tautological. Noted only because the constant is unusual-looking.
- The `Clear[chiPhi, hFun, ySoft]` in the general-softening block narrows `$Assumptions` to `{ySoft, gPhi}`. The earlier wider `$Assumptions` (including `npp > 0`, `kX > 0`, `i1 > 0`, etc.) is out of scope for this block. The block only references `gPhi` and `ySoft`, so no leakage problem. The subsequent `ELIMINATED-SOURCE SOFTENING CHECK` block redeclares `$Assumptions` to its own scope (`sourceCoeff > 0` etc.), so the narrowing does not propagate undesirably.
- The previously-existing verification file under this same path (older audit pass) catalogued four findings (F1-F4) against the pre-fix scripts. That file has been superseded by this verification, which addresses the two findings raised in the current red-team cycle.

## Verdict justification

Both findings are resolved. The sympy F1 patch exercises the general-H softening identity on a two-point branch with independent `H1, H2`, with no matched-layer substitution anywhere in the new block; the assertion is non-tautological. The Mathematica F1 patch, after the orchestrator hot-fix from integral-comparison to integrand-comparison, anchors the same physical claim to two integrand-level algebraic identities that are genuinely non-trivial (they require canceling `hFun[y]` from the `Theta` integrand and from the `Lambda` integrand pointwise), and the orchestrator hot-fix correctly preserves the physical claim because pointwise integrand equality implies integral equality. F2's rewrite produces a Mathematica audit whose closure, integral profile, Cauchy bound, and softening blocks all follow paths structurally distinct from the sympy choreography; only the `F_eff` block remains algebraically parallel and is permitted by the directive (variables renamed). Both exec logs exit 0, all `expect_zero`/`expectZero` assertions evaluate to 0, no regressions in the diff. Verdict: `verified`.
