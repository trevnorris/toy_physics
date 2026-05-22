---
unit_id: 045
batch: III.1
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T18:47:11Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 045

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:46-108`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:34-75`

**Issue:**

Three of the four advertised checks in the SymPy script (and their mirrors in the WL script) are tautological algebraic identities by construction. (1) `g_B g_R - g_W g_S` reduces to `a*b*c*d - a*b*c*d` because both products multiply the same four scalars. (2) `rho_0 - sigma_0` factors through the same cancellation. (3) `M_tr - M_tr_expected` is the distributive law `A*B + A*C - A*(B+C) = 0` since `M_mix` and `M_supp` share the common factor `8*(1+chi_0)^2/(pi^2*(1-eps_eta))`. None of these assertions can fail under any modification of the physical input that preserves the kernel's algebraic structure.

**Required change:**

Restructure the three assertions so that the left-hand expression is derived through a *different* route than the right-hand expression. Concretely:

Step 1 (SymPy, section 1, replace lines 46-51 of the .py):
Before:
```python
g_W = lam_W / sp.sqrt(mu_eta * mu_W)
g_R = gamma * lam_W / sp.sqrt(mu_U * mu_W)
g_B = lam_phi / sp.sqrt(mu_eta * mu_phi)
g_S = gamma * lam_phi / sp.sqrt(mu_U * mu_phi)

expect_zero("g_B g_R - g_W g_S", g_B * g_R - g_W * g_S)
```
After:
```python
W_sym, phi_sym, eta_sym, U_sym = sp.symbols("W_sym phi_sym eta_sym U_sym")
coupling_density = sp.expand((lam_W*W_sym + lam_phi*phi_sym) * (eta_sym - gamma*U_sym))
# Extract bilinear coefficients from the kernel directly.
c_W_eta   = coupling_density.coeff(W_sym).coeff(eta_sym)
c_W_U     = coupling_density.coeff(W_sym).coeff(U_sym)
c_phi_eta = coupling_density.coeff(phi_sym).coeff(eta_sym)
c_phi_U   = coupling_density.coeff(phi_sym).coeff(U_sym)
g_W_ext = c_W_eta   / sp.sqrt(mu_eta * mu_W)
g_R_ext = -c_W_U    / sp.sqrt(mu_U   * mu_W)
g_B_ext = c_phi_eta / sp.sqrt(mu_eta * mu_phi)
g_S_ext = -c_phi_U  / sp.sqrt(mu_U   * mu_phi)
# Reference (the form the script historically used).
g_W = lam_W / sp.sqrt(mu_eta * mu_W)
g_R = gamma * lam_W / sp.sqrt(mu_U * mu_W)
g_B = lam_phi / sp.sqrt(mu_eta * mu_phi)
g_S = gamma * lam_phi / sp.sqrt(mu_U * mu_phi)
# Cross-check: extracted vs reference (catches sign / coefficient errors in the kernel).
expect_zero("g_W extracted - reference", g_W_ext - g_W)
expect_zero("g_R extracted - reference", g_R_ext - g_R)
expect_zero("g_B extracted - reference", g_B_ext - g_B)
expect_zero("g_S extracted - reference", g_S_ext - g_S)
# Now the kernel-derived identity becomes a non-trivial check.
expect_zero("g_B g_R - g_W g_S", g_B_ext * g_R_ext - g_W_ext * g_S_ext)
```
The `-c_W_U` / `-c_phi_U` signs are because the kernel is `(lam_W W + lam_phi phi)(eta - gamma U)`, so the `W*U` and `phi*U` coefficients are negative; we want the magnitudes (the `g_R, g_S` amplitudes), hence the negation.

Step 2 (SymPy, section 1, lines 53-57):
After step 1, replace
```python
rho_0 = sp.simplify(g_R * g_U / (K_U * g_W))
sigma_0 = sp.simplify(g_U * g_S / (K_U * g_B))
```
with
```python
rho_0 = sp.simplify(g_R_ext * g_U / (K_U * g_W_ext))
sigma_0 = sp.simplify(g_U * g_S_ext / (K_U * g_B_ext))
```
Keep the existing `expect_zero("rho_0 - sigma_0", rho_0 - sigma_0)` assertion at line 57 unchanged; it now exercises the extracted amplitudes' identity, which is non-trivial conditional on the kernel structure.

Step 3 (SymPy, section 3, lines 97-108):
Before:
```python
M_mix = sp.simplify(8 * Z_W * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - eps_W_split)))
M_supp = sp.simplify(8 * Z_phi * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta) * (1 - eps_phi_split)))
M_tr = sp.simplify(M_mix + M_supp)
...
M_tr_expected = sp.simplify(
    8 * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta))
    * (Z_W / (1 - eps_W_split) + Z_phi / (1 - eps_phi_split))
)
expect_zero("M_tr - expected", M_tr - M_tr_expected)
```
After: derive M_tr via an enumerated channel sum so that the additive structure is genuinely tested, not just hand-written twice. Replace lines 97-108 with:
```python
channels = [
    ("W", Z_W, eps_W_split),
    ("phi", Z_phi, eps_phi_split),
]
prefactor = 8 * (1 + chi_0) ** 2 / (sp.pi ** 2 * (1 - eps_eta))
M_tr_channel_sum = sp.simplify(
    sum(prefactor * Z_i / (1 - eps_i) for (_, Z_i, eps_i) in channels)
)
M_mix = sp.simplify(prefactor * Z_W / (1 - eps_W_split))
M_supp = sp.simplify(prefactor * Z_phi / (1 - eps_phi_split))
M_tr = sp.simplify(M_mix + M_supp)
print("M_mix  =", M_mix)
print("M_supp =", M_supp)
print("M_tr   =", M_tr)
print("M_tr_channel_sum =", M_tr_channel_sum)
expect_zero("M_tr - channel_sum", M_tr - M_tr_channel_sum)
```
The check now confirms that the explicit two-term decomposition and the channel-enumerated sum agree. Either side can be wrong without the other (e.g., dropping a channel from `channels`).

Step 4 (Mathematica, section 1, lines 34-46): apply the analogous polynomial-extraction route in WL:
```mathematica
ClearAll[Wsym, phisym, etasym, Usym];
couplingDensity = Expand[(lamW*Wsym + lamPhi*phisym)*(etasym - gamma*Usym)];
cWeta = Coefficient[Coefficient[couplingDensity, Wsym], etasym];
cWU = Coefficient[Coefficient[couplingDensity, Wsym], Usym];
cPhiEta = Coefficient[Coefficient[couplingDensity, phisym], etasym];
cPhiU = Coefficient[Coefficient[couplingDensity, phisym], Usym];
gWext = cWeta/Sqrt[muEta*muW];
gRext = -cWU/Sqrt[muU*muW];
gBext = cPhiEta/Sqrt[muEta*muPhi];
gSext = -cPhiU/Sqrt[muU*muPhi];
gW = lamW/Sqrt[muEta*muW];
gR = gamma*lamW/Sqrt[muU*muW];
gB = lamPhi/Sqrt[muEta*muPhi];
gS = gamma*lamPhi/Sqrt[muU*muPhi];
expectZero["g_W extracted - reference", gWext - gW];
expectZero["g_R extracted - reference", gRext - gR];
expectZero["g_B extracted - reference", gBext - gB];
expectZero["g_S extracted - reference", gSext - gS];
rho0 = FullSimplify[gRext*gU/(KU*gWext), Assumptions -> $Assumptions];
sigma0 = FullSimplify[gU*gSext/(KU*gBext), Assumptions -> $Assumptions];
(* rest of section 1 unchanged *)
expectZero["g_B g_R - g_W g_S", gBext*gRext - gWext*gSext];
Print["rho_0 = ", fmt[rho0]];
Print["sigma_0 = ", fmt[sigma0]];
expectZero["rho_0 - sigma_0", rho0 - sigma0];
```
Replace the original lines 34-46 with this block (keep the `R_tr` and range identities at lines 48-56 unchanged).

Step 5 (Mathematica, section M_tr, lines 64-75): analogous channel-sum refactor in WL:
```mathematica
channels = {{ZW, epsWSplit}, {ZPhi, epsPhiSplit}};
prefactor = 8*(1 + chi0)^2/(Pi^2*(1 - epsEta));
mTrChannelSum = FullSimplify[
  Total[prefactor*#[[1]]/(1 - #[[2]]) & /@ channels],
  Assumptions -> $Assumptions
];
mMix = FullSimplify[prefactor*ZW/(1 - epsWSplit), Assumptions -> $Assumptions];
mSupp = FullSimplify[prefactor*ZPhi/(1 - epsPhiSplit), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
Print["M_mix = ", fmt[mMix]];
Print["M_supp = ", fmt[mSupp]];
Print["M_tr = ", fmt[mTr]];
Print["M_tr_channel_sum = ", fmt[mTrChannelSum]];
expectZero["M_tr - channel_sum", mTr - mTrChannelSum];
```
Drop the old `mTrExpected` line entirely; the `M_tr - channel_sum` assertion replaces it.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 045` and `redteam exec-mathematica 045` and confirm:
- The new SymPy script contains the literal strings `coupling_density = sp.expand(`, `c_W_eta`, `c_phi_U`, `M_tr_channel_sum`, and the assertion `g_W extracted - reference`.
- The new WL script contains the literal strings `couplingDensity = Expand[`, `cWeta`, `cPhiU`, `mTrChannelSum`, and the assertion `g_W extracted - reference`.
- Both scripts still exit 0 with residuals `= 0` on every `expectZero` line.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- summary: Replaced tautological coupling and total-baseline checks with kernel coefficient extraction and channel-sum derivations.
- deviation: Preserved the existing `rTr` assignment so the unchanged R_tr range checks retain their input.

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:96-98`

**Issue:**

`mTrReq` is defined as the literal target expression and then checked against the same target expression on the next assertion line, making `mTrReq - <target> == 0` a syntactic self-comparison. The corresponding SymPy section solves the quadratic for `M_tr_sym` and checks the solved expression against the closed form (a genuine test); the Mathematica side has lost that derivation step.

**Required change:**

Replace lines 96-98 of the `.wl` script:

Before:
```mathematica
mTrReq = FullSimplify[xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi), Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```

After:
```mathematica
mTrReqSolutions = Solve[collapsedNum == 0, mTrSym];
mTrReq = FullSimplify[mTrSym /. First[mTrReqSolutions], Assumptions -> $Assumptions];
Print["M_tr required on tracking branch = ", fmt[mTrReq]];
expectZero["G_tr generic formula", mTrReq - xi*(delta + xi)/(delta + (1 + lambda0*rU^2)*xi)];
```

`collapsedNum` is already defined at line 90 (`Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym]`) and is linear in `mTrSym`, so `Solve` returns a single root.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 045` and confirm:
- The new `.wl` contains the literal string `mTrReqSolutions = Solve[collapsedNum == 0, mTrSym]`.
- The `.txt` output still shows `M_tr required on tracking branch = (xi*(delta + xi))/(delta + xi + lambda0*rU^2*xi)` (modulo ordering) and `G_tr generic formula = 0`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- summary: Replaced the literal tracking-branch target with a Solve-derived expression before checking the closed form.
- deviation: none

## F3 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl:82-94`

**Issue:**

Beyond F1 (which already injects an independent polynomial-extraction route into section 1 of the WL script) and F2 (which adds a `Solve` step), section 4 of the WL script remains a token-for-token rewrite of section 4 of the PY script: identical `branchEq` rational, identical `rPhi -> rU` substitution, identical `Numerator` / `Denominator` choreography. The independence policy requires *at least one* substantively different algebraic route in the Mathematica script for the section-4 claim.

**Required change:**

Replace lines 82-94 of the `.wl` script. The new code should arrive at the same `numTrack` via a series-expansion route around `rPhi = rU` rather than direct substitution.

Before:
```mathematica
branchEq = FullSimplify[
  mSuppSym - (xi*(delta + xi) - mMixSym*(delta + (1 + lambda0*rU^2)*xi))/
    (delta + (1 + lambda0*rPhi^2)*xi - mMixSym*lambda0*(rU - rPhi)^2),
  Assumptions -> $Assumptions
];
branchTrack = Together[FullSimplify[branchEq /. rPhi -> rU, Assumptions -> $Assumptions]];
numTrack = Expand[Numerator[branchTrack]];
denTrack = Factor[Denominator[branchTrack]];
collapsedNum = Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym];

Print["tracking numerator = ", fmt[numTrack]];
Print["tracking denominator = ", fmt[denTrack]];
expectZero["tracking quadratic collapse", numTrack + (collapsedNum /. mTrSym -> mMixSym + mSuppSym)];
```

After:
```mathematica
branchEq = FullSimplify[
  mSuppSym - (xi*(delta + xi) - mMixSym*(delta + (1 + lambda0*rU^2)*xi))/
    (delta + (1 + lambda0*rPhi^2)*xi - mMixSym*lambda0*(rU - rPhi)^2),
  Assumptions -> $Assumptions
];
(* Independent route: expand the numerator of branchEq as a polynomial in (rPhi - rU) and read off the zeroth-order term, instead of substituting rPhi -> rU directly. *)
branchNumRaw = Together[branchEq] // Numerator // Expand;
seriesAtTrack = Series[branchNumRaw, {rPhi, rU, 0}] // Normal // Expand;
numTrack = Expand[seriesAtTrack];
branchDenRaw = Together[branchEq] // Denominator;
denTrack = Factor[branchDenRaw /. rPhi -> rU];
collapsedNum = Expand[xi^2 + (delta - mTrSym*(1 + lambda0*rU^2))*xi - delta*mTrSym];

Print["tracking numerator = ", fmt[numTrack]];
Print["tracking denominator = ", fmt[denTrack]];
expectZero["tracking quadratic collapse", numTrack + (collapsedNum /. mTrSym -> mMixSym + mSuppSym)];
```

The `Series[..., {rPhi, rU, 0}] // Normal` step extracts the rPhi -> rU value of the numerator via Taylor expansion rather than direct substitution; algebraically equivalent, structurally distinct from the PY route.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 045` and confirm:
- The new `.wl` contains the literal string `Series[branchNumRaw, {rPhi, rU, 0}]`.
- The `.txt` output still shows `tracking quadratic collapse = 0`.
- The displayed `tracking numerator` form matches the SymPy `tracking numerator` line modulo Mathematica's term ordering.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- summary: Replaced the direct tracking substitution route with a numerator series expansion around `rPhi = rU`.
- deviation: none
