---
unit_id: 045
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 045

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

SymPy (`scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:46-68`): added `coupling_density = sp.expand((lam_W*W_sym + lam_phi*phi_sym) * (eta_sym - gamma*U_sym))` and extracted `c_W_eta`, `c_W_U`, `c_phi_eta`, `c_phi_U` via `coupling_density.coeff(...).coeff(...)`. Defined `g_W_ext`, `g_R_ext`, `g_B_ext`, `g_S_ext` from those coefficients (with the documented `-c_W_U`/`-c_phi_U` sign for the `U` legs). Kept the original reference amplitudes and added four `expect_zero("g_X extracted - reference", ...)` assertions plus the cross-product identity on the `_ext` quantities. Lines 70-71 rewired `rho_0`/`sigma_0` through the extracted amplitudes. Lines 114-129 replaced the hand-written `M_tr_expected` with an enumerated `channels = [("W", Z_W, eps_W_split), ("phi", Z_phi, eps_phi_split)]` list and `M_tr_channel_sum = sp.simplify(sum(prefactor * Z_i / (1 - eps_i) for ...))`; the new assertion is `expect_zero("M_tr - channel_sum", M_tr - M_tr_channel_sum)`.

Mathematica (`mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-60` and `:72-91`): mirrored both refactors — `couplingDensity = Expand[...]`, `Coefficient[Coefficient[...]]` extraction, four `g_X extracted - reference` checks, then the cross-product identity on `gBext*gRext - gWext*gSext`; and an analogous `channels = {{ZW, epsWSplit}, {ZPhi, epsPhiSplit}}` / `Total[prefactor*#[[1]]/(1 - #[[2]]) & /@ channels]` form with the renamed assertion `M_tr - channel_sum`.

**Assessment:**

The directive required restructuring the three assertions so the left and right are constructed via different routes. Both scripts comply:

- The `g_B g_R - g_W g_S` assertion now reads `g_B_ext * g_R_ext - g_W_ext * g_S_ext`, with the `_ext` quantities sourced from `coefficient(...)` calls on the polynomial `coupling_density`. The four `expect_zero("g_X extracted - reference", ...)` lines (and WL mirrors) are the structural firewall: any modification to the kernel `(lam_W W + lam_phi phi)(eta - gamma U)` that changes a coupling sign or scaling will fail at least one of those checks before the cross-product line is reached. The cross-product itself remains algebraically tautological under the *current* kernel structure (the auditor anticipated this and accepted the polynomial-extraction route precisely because it shifts the failure mode upstream), and the extracted-vs-reference assertions provide the new adversarial signal.
- `rho_0`/`sigma_0` are recomputed through `_ext` amplitudes, so a kernel-structure error would now propagate into `rho_0 - sigma_0`.
- `M_tr_channel_sum` is genuinely an enumerated sum across a list of channel tuples; if the list were missing a channel (or had a channel with a wrong coupling), `M_tr - channel_sum` would no longer be zero. The pair-up against `M_tr = M_mix + M_supp` (each derived directly, not from the enumerated sum) makes the residual a meaningful check on the additive enumeration.

Output transcripts confirm all five new assertions pass on both engines (`g_W extracted - reference = 0`, ..., `g_B g_R - g_W g_S = 0`, `rho_0 - sigma_0 = 0`, `M_tr - channel_sum = 0`). No collateral edits beyond what the directive prescribed.

### F2 — tautological_check

**Classification:** resolved

**What changed:**

Mathematica (`mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:115-118`):
```
mTrReqSolutions = Solve[collapsedNum == 0, mTrSym];
mTrReq = FullSimplify[mTrSym /. First[mTrReqSolutions], Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```

**Assessment:**

The directive's exact replacement is in place. `collapsedNum` (defined at line 109 as `Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym]`) is linear in `mTrSym`, so `Solve` returns a single root and `First[mTrReqSolutions]` is unambiguous. The output shows the solved expression `M_tr required on tracking branch = (xi*(delta + xi))/(delta + xi + lambda0*rU^2*xi)`, matching (modulo ordering) the SymPy side `xi*(delta + xi)/(R_U**2*lambda_0*xi + delta + xi)`. The residual `G_tr generic formula = 0` is no longer a self-comparison: a wrong coefficient in `collapsedNum` would change the `Solve` output and the residual would be non-zero.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**

Mathematica (`mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:98-113`): replaced direct `branchTrack = Together[FullSimplify[branchEq /. rPhi -> rU, ...]]` route with a series-expansion route:
```
branchNumRaw = Together[branchEq] // Numerator // Expand;
seriesAtTrack = Series[branchNumRaw, {rPhi, rU, 0}] // Normal // Expand;
numTrack = Expand[seriesAtTrack];
branchDenRaw = Together[branchEq] // Denominator;
denTrack = Factor[branchDenRaw /. rPhi -> rU];
```

**Assessment:**

The directive's required `Series[branchNumRaw, {rPhi, rU, 0}] // Normal` construct is present and is structurally distinct from the SymPy `branch_eq.subs(R_phi, R_U)` route. The denominator is still extracted by direct substitution, which is acceptable since the auditor required only "at least one" route change in section 4 and the directive explicitly only refactors the numerator route.

The output `tracking numerator = delta*mMixSym + delta*mSuppSym - delta*xi + mMixSym*xi + mSuppSym*xi + lambda0*mMixSym*rU^2*xi + lambda0*mSuppSym*rU^2*xi - xi^2` matches the SymPy `tracking numerator = Mmix*R_U**2*lambda_0*xi + Mmix*delta + Mmix*xi + Msupp*R_U**2*lambda_0*xi + Msupp*delta + Msupp*xi - delta*xi - xi**2` modulo term ordering and the engine's symbol-name conventions (`mMixSym` vs `Mmix`). The `tracking quadratic collapse = 0` residual survives, confirming both routes produce the same numerator and the series expansion did not introduce extraneous higher-order terms (the `Series[..., 0]` request and the subsequent `Normal` strip retain only the zeroth-order term, which is the desired `rPhi = rU` value for a polynomial expression).

## Exec log assessment

**SymPy:** exit=0. Notable lines (from `scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt`):
- `g_W extracted - reference = 0` (and the three sibling checks)
- `g_B g_R - g_W g_S = 0`
- `rho_0 - sigma_0 = 0`
- `M_tr - channel_sum = 0`
- `tracking quadratic collapse = 0`, `G_tr formula = 0`, `G_tr D/N specialization = 0`, `F_tr normalization law = 0`
- Final line: `All Stage-28 symbolic checks passed.`

No `Traceback` or `AssertionError` in the transcript.

**Mathematica:** exit=0. Notable lines (from `mathematica/output/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.txt`):
- `PASS: g_W extracted - reference` (and the three sibling checks)
- `PASS: g_B g_R - g_W g_S`
- `PASS: rho_0 - sigma_0`
- `PASS: M_tr - channel_sum`
- `M_tr required on tracking branch = (xi*(delta + xi))/(delta + xi + lambda0*rU^2*xi)`
- `PASS: G_tr generic formula`, `PASS: G_tr D/N specialization`, `PASS: F_tr normalization law`
- Final line: `Stage 045 Mathematica audit passed.`

No `$Failed` or `FAIL:` line in the transcript.

**Output freshness:**
- `scripts/.../sympy_audit.py` mtime 2026-05-22 12:47:33; `.txt` 2026-05-22 12:50:29 — output newer than script, fresh.
- `mathematica/.../mathematica_audit.wl` mtime 2026-05-22 12:48:23; `.txt` 2026-05-22 12:50:39 — output newer than script, fresh.

## Material-change assessment

`material_change`: false.

The edits restructure how three identities are *verified*; they do not change any derived numerical or symbolic result that downstream units consume. The displayed forms of `rho_0`, `R_tr`, `M_tr`, `M_tr required on tracking branch`, `G_tr D/N specialization`, and the `F_tr` normalization law are identical to the pre-fix outputs (the pre-fix forms were already algebraically correct; only their verification route was tautological). No downstream unit re-uses an intermediate symbol that has been redefined.

## Side observations (non-blocking)

- The cross-product identity `g_B g_R - g_W g_S` remains algebraically tautological under the current kernel (`g_W_ext = lam_W/sqrt(mu_eta*mu_W)` etc.; cross-product still factors to `gamma*lam_W*lam_phi/sqrt(mu_eta*mu_U*mu_W*mu_phi)` on both sides). The directive accepted this because the four `g_X extracted - reference` assertions move the failure mode upstream into the polynomial extraction. This is the spec the auditor signed off on, so it does not block verification, but a future auditor may want to assert the cross-product on a *deformed* kernel (e.g., with a perturbed `gamma`) to make even the cross-product line non-trivial.
- `coherent normalization residual` is a non-zero printed residual on both engines that involves a free `R_target`/`rTarget` symbol; it is a diagnostic print, not an assertion, and was present before the fix.

## Verdict justification

All three findings are resolved: F1 swaps in a polynomial-extraction route plus channel-sum enumeration on both engines, F2 introduces a `Solve` step on the Mathematica side that replaces the prior self-comparison, and F3 substitutes a `Series[..., {rPhi, rU, 0}]` route for the direct `rPhi -> rU` substitution in the Mathematica branch-numerator extraction. Both exec transcripts exit 0 with every `expectZero`/`expect_zero` reporting residual 0. Outputs are fresh relative to the edited sources. No regressions in the diff, no collateral edits beyond what the directive requested. Verdict: `verified`.
