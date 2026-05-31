---
unit_id: 182
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-30T01:26:48-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

## Consult resolution (Claude+Codex read-only, 2026-05-30, codex_session 019e77af)

**Q3 — F2 independence (DISPUTE of `## Blocked`; a genuine lighter route exists):** do NOT
default to `## Blocked: F2`. The full 8-log→5-Σ inverse is under-determined, but a concrete
in-scope linear route works: **gauge-fix / `Solve` (or `SolveAlways`/`CoefficientArrays`) five
of the eight microscopic logs from the five Σ definitions, substitute into `xi1Direct`,
`Collect` on the Σ symbols, and `expectZero` each Σ-coefficient against the paper law** (real
fail mode: no solution, leftover free logs, or coefficient mismatch). Derive the split via
`sigmaChi -> (sigmaTr - (1 + chi0) sigmaDel)/(1 + deltaU)`. Stay in scope — do NOT reconstruct
the upstream `zeta` cancellation. Append `## Blocked: F2` ONLY if that concrete linear reduction
genuinely fails — never for a cosmetic rename.

# Codex directive — unit 182

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification (support-blindness checks are vacuous)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py:139-144`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl:107-112`

**Issue:**
The support symbols `lamphi1`, `kphi` are declared (sympy L40; math L28-29) but never wired into any expression. The six support-blindness checks per engine differentiate the abstract-Sigma forms (`Xi1_slip`, `R1_slip`, `Theta1_fact`) with respect to symbols those expressions do not contain, so every check is identically `0` by construction and passes regardless of the physics. The microscopic-log "direct" forms (`Xi1_direct`, `R1_direct`, `Theta1_direct`) likewise do not contain `lamphi1`/`kphi`. The honest in-scope check is a structural free-symbol-absence assertion on the microscopic-log direct forms, which is non-vacuous (it fails if a support drift is ever introduced into the defect construction). Do NOT attempt to reconstruct the full `\zeta`-cancellation here — the notes (§2.2) delegate that mechanism to Stage 249, so building it here would be an out-of-scope scope extension.

**Required change (SymPy):**
Replace lines 139-144:

Before:
```python
expect_zero("dXi_1/dlamphi1", sp.diff(Xi1_slip, lamphi1))
expect_zero("dXi_1/dkphi", sp.diff(Xi1_slip, kphi))
expect_zero("dR_1/dlamphi1", sp.diff(R1_slip, lamphi1))
expect_zero("dR_1/dkphi", sp.diff(R1_slip, kphi))
expect_zero("dTheta_1/dlamphi1", sp.diff(Theta1_fact, lamphi1))
expect_zero("dTheta_1/dkphi", sp.diff(Theta1_fact, kphi))
```

After (assert support symbols are structurally absent from the microscopic-log direct forms):
```python
# Support-blindness: the support-lane drifts never enter the microscopic-log
# defect construction. (The zeta-cancellation mechanism itself lives upstream in
# Stage 249; here we verify only that no support log was wired into Xi_1/R_1/Theta_1.)
support_syms = {lamphi1, kphi}
for label, form in [("Xi_1 direct", Xi1_direct), ("R_1 direct", R1_direct), ("Theta_1 direct", Theta1_direct)]:
    leaked = support_syms & form.free_symbols
    print(f"{label} support-symbol leakage = {leaked}")
    if leaked:
        raise AssertionError(f"{label} unexpectedly depends on support drifts {leaked}")
```

Note: `Theta1_direct` already exists at sympy lines 105-108, so no new construction is needed.

**Required change (Mathematica):**
Replace lines 107-112:

Before:
```
expectZero["dXi_1/dlamphi1", D[xi1Slip, lamphi1]];
expectZero["dXi_1/dkphi", D[xi1Slip, kphi]];
expectZero["dR_1/dlamphi1", D[r1Slip, lamphi1]];
expectZero["dR_1/dkphi", D[r1Slip, kphi]];
expectZero["dTheta_1/dlamphi1", D[theta1Fact, lamphi1]];
expectZero["dTheta_1/dkphi", D[theta1Fact, kphi]];
```

After (assert support symbols are structurally absent from the microscopic-log direct forms via `FreeQ`):
```
(* Support-blindness: no support-lane drift enters the microscopic-log defect
   construction. The zeta-cancellation mechanism lives upstream in Stage 249. *)
Module[{forms = {{"Xi_1 direct", xi1Direct}, {"R_1 direct", r1Direct}, {"Theta_1 direct", theta1Direct}}},
  Do[
    With[{label = forms[[i, 1]], form = forms[[i, 2]]},
      If[FreeQ[form, lamphi1] && FreeQ[form, kphi],
        pass[label <> " support-blind"],
        fail[label <> " support-blind", form]
      ]
    ],
    {i, Length[forms]}
  ]
];
```
`theta1Direct` already exists at math lines 80-84, so no new construction is needed.

**Banner label fix (fold into this finding, same files):**
- sympy line 34: change `banner("STAGE 165 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION")` to `banner("STAGE 182 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION")`.
- math line 26: change `banner["STAGE 165 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION"];` to `banner["STAGE 182 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION"];`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 182` and `redteam exec-mathematica 182` and confirm: (a) the support-blindness section now reports a structural-absence result (empty leakage set / `PASS … support-blind`) rather than six `= 0` derivative lines; (b) both scripts exit 0; (c) the banner reads STAGE 182.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl`
- summary: Replaced vacuous support-derivative checks with structural support-symbol absence checks on the direct forms and fixed both banners to Stage 182.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl:61-103`

**Issue:**
The `.wl` re-types the SymPy script's hand-typed target forms (`xi1Slip` L63-68, `theta1Fact` L85-88, `xi1Split` L95-101) verbatim, with the same coefficients and the same `slipSubs` map (L45-51) as the `.py`. This makes the second engine an echo of the first rather than an independent re-derivation: an algebra error in the SymPy literals would be reproduced, not caught. The defect/`R_1`/`Theta_1`/split constructions and assertion ordering mirror the `.py` line-for-line.

**Required change:**
Make the Mathematica verification of the four-slippage law and the tracking/nontracking split derive its target forms by Mathematica algebra rather than by re-typing the SymPy coefficients. Concretely, for the four-slippage law (currently L63-69):

- Keep `xi1Direct` (L62) built from the Stage-249 microscopic-log expressions.
- Substitute the inverse slip relations into `xi1Direct` and let Mathematica reduce it, then `Collect` on the Sigma symbols and verify it matches the claimed law by coefficient comparison rather than by subtracting a hand-typed `xi1Slip`. For example:
  ```
  invSlip = {
    2*lam1 + mu1 - kEta - 2*kW -> sigmaZ,  (* not directly substitutable; instead express logs via Sigmas *)
  };
  ```
  Since the microscopic logs cannot be inverted one-to-one to Sigmas, the robust independent route is: form `xi1DirectInSigmas = FullSimplify[xi1Direct] /. (solve the five Sigma definitions for five of the eight logs)`, then `Collect[xi1DirectInSigmas, {sigmaZ, sigmaChi, sigmaEps, sigmaDel}]` and assert each coefficient equals the paper's coefficient via `expectZero` on `Coefficient[xi1DirectInSigmas, sigmaChi] - 2*chi0/(1 + chi0)` etc.

  This makes Mathematica *compute* the four coefficients from the direct form instead of being handed `xi1Slip`. Do the analogous coefficient-extraction for `xi1Split` (verify `Coefficient` on each Sigma matches the §6 split coefficients).

If the inverse substitution proves under-determined (eight logs, five Sigmas), an acceptable lighter-weight independence is: keep the existing `expectZero` structure but build the comparison target `xi1Slip` by `Collect`/`Simplify` of `xi1Direct` re-expressed through the Sigma combinations (using `Eliminate`/`Solve` to introduce the Sigma symbols), so the `.wl`'s `xi1Slip` is *derived*, not copied. If you cannot achieve genuine independence without scope creep, append `## Blocked: F2` describing the obstruction rather than making a cosmetic change that does not improve independence.

**Verification command:**
The verifier confirms the `.wl` constructs at least the four-slippage-law and split target coefficients by Mathematica algebra (`Collect`/`Coefficient`/`Solve`) rather than re-typing the SymPy literals, and that `math -script` exits 0 with all checks passing.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl`
- summary: Re-derived the Mathematica four-slippage and tracking-split forms by solving Sigma definitions, collecting coefficients, and checking those coefficients directly.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py:74-79`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl:54-59`

**Issue:**
The six "Physical branch drifts" checks subtract a definition from itself (e.g. sympy L74 `zetaZ - (2*lam1 - keta - kW)` where `zetaZ` was defined as exactly that at L46). They cannot fail and verify nothing.

**Required change:**
Replace the six tautological checks with the notes §3 boxed equalities, which are non-tautological because they relate independently-defined quantities. The notes state `\Sigma_\chi = \chi_1/\chi_0`, `\Sigma_\eta = \eta_1/\epsilon_\eta`, `\Sigma_\delta = \delta_{U,1}/\delta_U`. Since `chi1`, `eta1`, `deltaU1` are defined independently of the Sigma symbols (sympy L48-51), these are real checks.

SymPy — replace lines 74-79:
```python
expect_zero("Sigma_chi = chi_1/chi_0", SigmaChi.subs(slip_subs) - chi1 / chi0)
expect_zero("Sigma_eta = eta_1/eps_eta", SigmaEta.subs(slip_subs) - eta1 / eps_eta)
expect_zero("Sigma_del = delta_U,1/delta_U", SigmaDel.subs(slip_subs) - deltaU1 / deltaU)
expect_zero("Sigma_eps = varepsilon_W/eps_W", SigmaEps.subs(slip_subs) - varepsW / epsW)
expect_zero("Sigma_Z = zeta_Z - omega_W", SigmaZ.subs(slip_subs) - (zetaZ - omegaW))
```
(Note `\Sigma_Z = \zeta_Z - \omega_W` is notes §3 line 241-242; `\Sigma_\epsilon = \varepsilon_W/\epsilon_W` is notes §3 line 256.) `SigmaChi.subs(slip_subs)` expands `SigmaChi` to `gam1 + c1 - kU`; `chi1/chi0` = `(gam1 + c1 - kU)`, so the residual is 0 non-trivially (it would be nonzero if `chi1`'s definition were wrong).

Mathematica — replace lines 54-59 analogously:
```
expectZero["Sigma_chi = chi_1/chi_0", (sigmaChi /. slipSubs) - chi1/chi0];
expectZero["Sigma_eta = eta_1/eps_eta", (sigmaEta /. slipSubs) - eta1/epsEta];
expectZero["Sigma_del = delta_U,1/delta_U", (sigmaDel /. slipSubs) - deltaU1/deltaU];
expectZero["Sigma_eps = varepsilon_W/eps_W", (sigmaEps /. slipSubs) - varepsW/epsW];
expectZero["Sigma_Z = zeta_Z - omega_W", (sigmaZ /. slipSubs) - (zetaZ - omegaW)];
```

If you prefer the conservative path, simply delete lines 74-79 (sympy) / 54-59 (math) entirely; both are acceptable, but the replacement above is preferred because it anchors deliverable (5) to the notes §3 equalities.

**Verification command:**
The verifier confirms the six tautological self-subtractions are gone, replaced by the five notes-§3 equality checks (or deleted), and that both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage182_microscopic_coherent_slippage_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage182_microscopic_coherent_slippage_mathematica_audit.wl`
- summary: Replaced the tautological microscopic drift formula checks with the notes §3 Sigma equality checks.
- deviation: none
