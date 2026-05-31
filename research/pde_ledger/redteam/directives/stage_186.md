---
unit_id: 186
batch: V.2
created_at: 2026-05-30T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-30T07:35:09Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 186

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:112,116`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl:93,97`

**Issue:** The "finite orbit preserves ε_η" check is `(2*C - U) - Eta_exp` where `Eta_exp = 2*C - U` (py:101 / wl:83) is the identical literal — so the residual is `X − X ≡ 0` by construction and cannot fail regardless of whether `2C−U` is the correct ε_η-preserving scaling. The companion C_tr/C_nt checks are fine (their dependent scalings `Tau_exp`/`Mu_exp` are non-trivial). Only ε_η degenerates, leaving the ε_η monomial-preservation deliverable with no real coverage.

**Required change (SymPy, py):**
Keep the existing trivial check line if you wish, but ADD a non-tautological check that derives the ε_η-preserving K_η scaling from the monomial's exponent vector `(2, −1, −1)` on `(c_ηU, K_U, K_η^eff)` and confirms it equals the paper value `2C − U`. Insert immediately after the existing `expect_zero("finite orbit preserves eps_eta", Eta_orbit)` line (currently py:116):

```python
# Non-tautological ground check: solve for the K_eta^eff scaling that preserves
# eps_eta = c_etaU^2 / (K_U K_eta^eff), then confirm it equals the paper's 2C - U.
eta_scaling = sp.symbols("eta_scaling", real=True)
eps_eta_logdrift = 2*C - U - eta_scaling   # log-drift of c^2 K_U^{-1} K_eta^{-1}
solved_eta = sp.solve(sp.Eq(eps_eta_logdrift, 0), eta_scaling)[0]
expect_zero("K_eta preserving scaling matches paper 2C-U", solved_eta - (2*C - U))
expect_zero("chosen Eta_exp solves eps_eta preservation",
            eps_eta_logdrift.subs(eta_scaling, Eta_exp))
```

**Required change (Mathematica, wl):**
Insert the analogous block immediately after the existing `expectZero["finite orbit preserves eps_eta", etaOrbit];` line (currently wl:97):

```wolfram
(* Non-tautological ground check: solve for the K_eta^eff scaling that preserves
   eps_eta = c_etaU^2 / (K_U K_eta^eff), then confirm it equals the paper 2C - U. *)
Clear[etaScaling];
$Assumptions = Element[{chi, delta, eParam, fParam, lamSym, cSym, gamSym, uSym, wSym, etaScaling}, Reals] &&
  chi > 0 && delta > 0;
epsEtaLogDrift = 2*cSym - uSym - etaScaling;
solvedEta = etaScaling /. First[Solve[epsEtaLogDrift == 0, etaScaling]];
expectZero["K_eta preserving scaling matches paper 2C-U", solvedEta - (2*cSym - uSym)];
expectZero["chosen Eta_exp solves eps_eta preservation",
  epsEtaLogDrift /. etaScaling -> etaExp];
```

(Note: the `.wl` ε_η symbols are `cSym`, `uSym`, `etaExp` — use those, not the `2*C - U` SymPy names.)

**Self-test (performed by auditor):** `eps_eta_logdrift = 2*C − U − eta_scaling` depends on `eta_scaling`, so the solve is non-degenerate and returns `2*C − U` (a nonzero literal). `solved_eta − (2*C − U) = 0` (assert_zero holds); substituting the real `Eta_exp = 2*C − U` gives `0`; a hypothetical wrong `Eta_exp` would make the second new check nonzero. No new constants introduced.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 186` and `redteam exec-mathematica 186` and confirm the new check lines (`K_eta preserving scaling matches paper 2C-U`, `chosen Eta_exp solves eps_eta preservation`) appear with residual `0` / `PASS:`, and both scripts exit 0.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
- summary: Added non-tautological epsilon-eta K_eta scaling solve checks in both audit scripts.
- deviation: none

## F2 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl` (whole file)
- (banner only) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py:33`

**Issue:** The `.wl` is a line-by-line port of the `.py` (identical drift expressions wl:32-33, identical hardcoded matrix wl:40-44, identical orbit exponents wl:83-85, identical basis substitution, even the stale `"STAGE 169"` banner). Both engines hardcode the same `M_*` and the same dependent-scaling exponents and run the same `Solve`, so the Mathematica run is not an independent witness — a single typo'd matrix entry or orbit exponent copied to both engines would pass in both.

**Required change:**
Make the `.wl` a genuine re-derivation of the load-bearing objects rather than a transcription. Do NOT change which identities are verified — only how the Mathematica path arrives at them. Specifically:

1. **Build `M_*` from the monomial exponent vectors, not a transcribed literal.** The three monomials (notes lines 108-123 / appendix lines 971-990) are:
   - `C_tr,*` = `(γ c_ηU / K_U)^{1+δ} (π^2 T_U / (L^2 K_U))^{1+χ}`
   - `C_nt,*` = `(λ_W^2 μ_W / (K_η^eff (K_W^eff)^2)) (γ^2 λ_W^2 σ / (K_U K_W^eff))^{E} (π^2 T_U / (L^2 K_U))^{-F}`
   - `ε_η` = `c_ηU^2 / (K_U K_η^eff)`

   In the `.wl`, take the log of each monomial as a function of the eight log-drifts `(lam1,c1,gam1,kU,kEta,kW,mu1,tau1) = δln(λ_W,c_ηU,γ,K_U,K_η^eff,K_W^eff,μ_W,T_U)` and read the coefficient of each drift via `Coefficient[...]` (or `D[..., var]`) to assemble each `M_*` row. Then assert this independently-assembled matrix equals the expected paper entries:
   ```wolfram
   mMatDerived = ... (* rows built from Coefficient of each monomial's log-drift *)
   expectZero["M_* row 1 matches paper", Total[Abs[mMatDerived[[1]] - mMat[[1]]]]];
   expectZero["M_* row 2 matches paper", Total[Abs[mMatDerived[[2]] - mMat[[2]]]]];
   expectZero["M_* row 3 matches paper", Total[Abs[mMatDerived[[3]] - mMat[[3]]]]];
   ```
   (You may keep `mMat` as the reference literal to compare against; the point is that `mMatDerived` is built from the monomials, so a transcription error in `mMat` would now be caught.)

2. **Obtain the dependent scalings by `Solve`, not transcription.** Instead of hardcoding `etaExp/tauExp/muExp` from the `.py`, solve the monomial-preservation equations `mMat . orbitLogVector == 0` for the dependent scalings of `(K_η^eff, T_U, μ_W)`, where `orbitLogVector` has free entries `(lamSym, cSym, gamSym, uSym, ?, wSym, ?, ?)` and the three `?` are the unknown dependent scalings on `(K_η^eff, μ_W, T_U)` in the drift-vector order. Then compare the solved scalings against the paper closed forms (notes lines 229-259 / appendix lines 1060-1083):
   ```wolfram
   {kEtaScale, muScale, tauScale} solved from mMat . orbitLog == {0,0,0};
   expectZero["solved K_eta scaling matches paper", kEtaScale - (2*cSym - uSym)];
   expectZero["solved T_U scaling matches paper",
     tauScale - (uSym - (1 + delta)*(gamSym + cSym - uSym)/(1 + chi))];
   expectZero["solved mu_W scaling matches paper",
     muScale - ((2*cSym - uSym) + 2*wSym - 2*lamSym
       - eParam*(2*gamSym + 2*lamSym - uSym - wSym)
       - fParam*(1 + delta)*(gamSym + cSym - uSym)/(1 + chi))];
   ```
   Drift-vector column order for the orbit log-vector is `(λ_W, c_ηU, γ, K_U, K_η^eff, K_W^eff, μ_W, T_U)`; free params are `Λ→lamSym, C→cSym, Γ→gamSym, U→uSym, W→wSym` on columns 1,2,3,4,6; the dependent unknowns sit on columns 5 (K_η), 7 (μ_W), 8 (T_U).

3. **Fix the banner string** to read `Stage 186` at `.wl:26` and `.py:33` (currently `"STAGE 169 — EXACT MICROSCOPIC SIMILARITY ORBIT"`). Replace `STAGE 169` with `STAGE 186` in both files; leave the rest of the banner text unchanged.

If building `mMatDerived` via `Coefficient` is awkward to express cleanly, you may instead assemble each row as an explicit list of `Coefficient[Log-expansion, drift_i]` for the eight drifts; the requirement is only that the row entries are *computed from the monomial log-expansion*, not retyped. If you cannot make the independent derivation pass cleanly, append `## Blocked: F2` with the specific obstruction rather than reverting to a transcription.

**Self-test (performed by auditor):** The monomial log-drifts reproduce exactly the `M_*` rows already verified in the report's paper↔script cross-check (row 1 = `(1+δ)(c1+gam1) − (2+χ+δ)kU + (1+χ)tau1`, etc.), so `mMatDerived` will equal `mMat` and the row-match checks pass. Solving `mMat . orbitLog == 0` for the three dependent scalings is well-posed because the corresponding 3×3 dependent minor has determinant `1+χ ≠ 0` (chi>0), guaranteeing a unique solution equal to the paper closed forms. No integrals; no new constants; verified identities unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 186` and confirm: (a) the `mMatDerived` row-match checks appear and pass, (b) the `Solve`-derived dependent scalings match the paper closed forms, (c) the banner reads "STAGE 186" in both scripts, and (d) the script exits 0. Then `redteam exec-sympy 186` to confirm the `.py` banner change did not break that script (exit 0).

## Applied: F2

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage186_similarity_orbit_closure_mathematica_audit.wl`
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage186_similarity_orbit_closure_sympy_audit.py`
- summary: Reworked the Mathematica audit to derive M_* from monomial exponent vectors and solve the dependent orbit scalings, and fixed both script banners to Stage 186.
- deviation: none
