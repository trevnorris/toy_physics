---
unit_id: 185
batch: V.2
created_at: 2026-05-30T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-30T09:04:00-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 185

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond the named edits. Edit only the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py:196-209`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl:188-202`

**Issue:**
The "Theta_1 monomial law" and "Xi_1 monomial law" checks form residuals `C_tr_star*(Sigma_tr - Sigma_tr_direct)` and `A_tr_star*(Sigma_tr - Sigma_tr_direct) + (Sigma_nt - Sigma_nt_direct)`. Because `Sigma_tr == Sigma_tr_direct` and `Sigma_nt == Sigma_nt_direct` are already proven earlier in the script, these residuals are `coeff·0 ≡ 0` for *any* value of the coefficients `C_tr_star`, `A_tr_star`. The observable coefficients are therefore never exercised. Fix by making the two "monomial law" residuals run the coefficient *against the primitive-compiler drift* (`Sigma_tr_compiled`, line 164; `Sigma_nt_compiled`, line 181) instead of `Sigma_tr_direct` / `Sigma_nt_direct`, so the check fails if the primitive-exponent compiler is wrong AND keeps the coefficient load-bearing; and additionally add an explicit coefficient-definition anchor so a wrong `C_tr_star`/`A_tr_star` literal is caught.

**Required change (SymPy):**

In the `banner("Observable triangular law in microscopic monomials")` block (lines 196-209), keep `C_tr_star`, `A_tr_star`, `Theta1`, `Xi1`, `Rcombo` as-is. Change the two existing law checks (lines 202-203) so the comparison form is built from the primitive compiler rather than `_direct`:

Before:
```
expect_zero("Theta_1 monomial law", Theta1 - (-C_tr_star * Sigma_tr_direct))
expect_zero("Xi_1 monomial law", Xi1 - (A_tr_star * Sigma_tr_direct + Sigma_nt_direct))
```
After:
```
expect_zero("Theta_1 monomial law", Theta1 - (-C_tr_star * Sigma_tr_compiled))
expect_zero("Xi_1 monomial law", Xi1 - (A_tr_star * Sigma_tr_compiled + Sigma_nt_compiled))
```
Then, immediately after the `C_tr_star`/`A_tr_star` definitions (after line 198), add two explicit coefficient anchors tying the coefficients to their paper-stated closed forms (appendix `eq:app-part05-Ctr-Atr-defs`):
```
expect_zero(
    "C_tr,* coefficient form",
    C_tr_star - chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)),
)
expect_zero(
    "A_tr,* coefficient form",
    A_tr_star - 2 * chi0s / ((1 + chi0s) * (1 + deltaUs)),
)
```

**Required change (Mathematica):**

In the `banner["Observable triangular law in microscopic monomials"]` block (lines 188-202), change lines 197-198 from `sigmaTrDirect`/`sigmaNtDirect` to `sigmaTrCompiled`/`sigmaNtCompiled`:

Before:
```
expectZero["Theta_1 monomial law", theta1 - (-cTrStar*sigmaTrDirect)];
expectZero["Xi_1 monomial law", xi1 - (aTrStar*sigmaTrDirect + sigmaNtDirect)];
```
After:
```
expectZero["Theta_1 monomial law", theta1 - (-cTrStar*sigmaTrCompiled)];
expectZero["Xi_1 monomial law", xi1 - (aTrStar*sigmaTrCompiled + sigmaNtCompiled)];
```
Then immediately after the `cTrStar`/`aTrStar` definitions (after line 190) add:
```
expectZero[
  "C_tr,* coefficient form",
  cTrStar - chi0s*deltaUs/((1 + chi0s)*(1 + deltaUs)*(1 + chi0s + deltaUs))
];
expectZero[
  "A_tr,* coefficient form",
  aTrStar - 2*chi0s/((1 + chi0s)*(1 + deltaUs))
];
```

Note: `sigmaTrCompiled` / `sigmaNtCompiled` (sympy lines 164/181; wl lines 156/173) already exist; you are only re-pointing the existing law checks at them and adding the two coefficient anchors. Do not introduce new variables beyond the two anchor checks.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- summary: Repointed the Theta_1/Xi_1 monomial laws to the primitive-compiled drifts and added explicit C_tr,* / A_tr,* coefficient-form anchors.
- deviation: none

---

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py:211-220`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl:204-213`

**Issue:**
The appendix attributes to Stage 185 the dependent-block determinant `det M_*^(τ,κη,μ) = 1+χ0* > 0` (the rank-3 / nonsingularity fact Stage 186 consumes). The scripts demonstrate nonsingularity only implicitly via the unique compatibility solve; the determinant value is never asserted. Add an explicit determinant check that BUILDS the 3×3 dependent minor from the already-defined drift expressions `Sigma_tr`, `Sigma_nt`, `Sigma_eta` by differentiating w.r.t. the dependent triple `(tau1, keta, mu1)` — do NOT hardcode the matrix entries — and verify the determinant equals `1+chi0s`.

**Required change (SymPy):**

At the start of the `banner("Exact zero-defect compatibility solve")` block (insert after line 211, before `tau_sol = ...`), add:
```
Mstar_minor = sp.Matrix(
    [
        [sp.diff(Sigma_tr, tau1), sp.diff(Sigma_tr, keta), sp.diff(Sigma_tr, mu1)],
        [sp.diff(Sigma_nt, tau1), sp.diff(Sigma_nt, keta), sp.diff(Sigma_nt, mu1)],
        [sp.diff(Sigma_eta, tau1), sp.diff(Sigma_eta, keta), sp.diff(Sigma_eta, mu1)],
    ]
)
expect_zero("det M_*^(tau,keta,mu) - (1+chi0s)", Mstar_minor.det() - (1 + chi0s))
```

**Required change (Mathematica):**

At the start of the `banner["Exact zero-defect compatibility solve"]` block (insert after line 204, before `tauSol = ...`), add:
```
mStarMinor = {
  {D[sigmaTr, tau1], D[sigmaTr, keta], D[sigmaTr, mu1]},
  {D[sigmaNt, tau1], D[sigmaNt, keta], D[sigmaNt, mu1]},
  {D[sigmaEta, tau1], D[sigmaEta, keta], D[sigmaEta, mu1]}
};
expectZero["det M_*^(tau,keta,mu) - (1+chi0s)", Det[mStarMinor] - (1 + chi0s)];
```

Self-test (already verified by auditor): the differentiated minor is
`[[1+chi0s, 0, 0], [-F_star, -1, 1], [0, -1, 0]]`, whose determinant (cofactor along row 1) is `(1+chi0s)·det([[-1,1],[-1,0]]) = (1+chi0s)·1 = 1+chi0s`. The residual `det - (1+chi0s)` is therefore identically 0 (check passes), and the determinant is a nonzero literal `1+chi0s` (the load-bearing positivity), so the check is non-tautological — a wrong drift coefficient anywhere in the dependent block changes the determinant and fails the check.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- summary: Added the differentiated dependent-block minor determinant check asserting det M_*^(tau,keta,mu) = 1 + chi0s.
- deviation: none

---

## F3 — stale banner label (non-blocking hygiene)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py:41`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl:26`

**Issue:**
Both scripts print `STAGE 168 — DIRECT MICROSCOPIC MONOMIALS`; this is Stage 185 (Stage 168 is "Off-bundle slippage decomposition"). Print label only; no assertion depends on it.

**Required change:**

SymPy line 41:
```
banner("STAGE 185 — DIRECT MICROSCOPIC MONOMIALS")
```
Mathematica line 26:
```
banner["STAGE 185 — DIRECT MICROSCOPIC MONOMIALS"];
```

**Verification command:**
After Codex applies F1-F3, the verifier runs `redteam exec-sympy 185` and `redteam exec-mathematica 185` and confirms: the new `C_tr,* coefficient form`, `A_tr,* coefficient form`, and `det M_*^(tau,keta,mu) - (1+chi0s)` checks appear and read `= 0`; the re-pointed "Theta_1/Xi_1 monomial law" checks still read `= 0`; the banner reads "STAGE 185"; both scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- summary: Updated both printed stage banners from Stage 168 to Stage 185.
- deviation: none

---

## Iteration 2 (delta) — F1 coefficients C_tr,* / A_tr,* are still not load-bearing (CHECKPOINT)

Verifier (`redteam/verifications/stage_185.md`) confirmed **F2 (det minor) and F3 (banner)
RESOLVED**. F1 is still `partial`: the iter-1 remedy did not meet its own goal. The two
"coefficient form" anchors (sympy ~L199–206 / wl ~L191–198) subtract the SAME paper formula from
itself (literal X−X), and the re-pointed law `Theta1 - (-C_tr_star*Sigma_tr_compiled)` multiplies
the coefficient onto `(Sigma_tr - Sigma_tr_compiled)`, which is **proven zero at sympy:166 / wl:158**
— so a wrong `C_tr_star`/`A_tr_star` still passes silently. A normalization like
`C_tr_star*(D/(chi0s*deltaUs)) - 1` does NOT help: since `C_tr_star` IS that ratio, it is also X−X.

**Resolution (Claude+Codex consult, 2026-05-30, codex_session 019e77e6 — route (b)):** reconstruct
the observable defects `Theta_1` and `Xi_1` INDEPENDENTLY from the microscopic slippage drifts
(`Sigma_chi`, `Sigma_delta`, `Sigma_Z`, `Sigma_eps` — already defined/verified in this script),
NOT via `C_tr_star`/`A_tr_star`. Then comparing the independent `Theta_1`/`Xi_1` against
`-C_tr_star*Sigma_tr` / `(A_tr_star*Sigma_tr + Sigma_nt)` makes the residual
`(C_typed − C_true)*Sigma_tr`, which is symbolically NONZERO for a wrong typed coefficient
(`Sigma_tr` is not zero in this block).

**Target:** the "Observable triangular law" block — sympy ~L197–211 / wl mirror ~L189–206.

**Required change (SymPy):** DELETE the two tautological `expect_zero("C_tr,* coefficient form", …)`
and `expect_zero("A_tr,* coefficient form", …)` anchors. Keep `C_tr_star`/`A_tr_star` typed
(L197–198). REPLACE the `Theta1`/`Xi1` definitions (L207–208) with independent slippage-law
constructions and add the independent-law checks:
```python
# chi_1, deltaU_1 are the first-order microscopic drifts (NOT functions of C_tr,*/A_tr,*)
chi1_indep    = sp.simplify(chi0s   * Sigma_chi)
deltaU1_indep = sp.simplify(deltaUs * Sigma_delta)
# Independent observable defects, built from the slippage drifts (appendix Stage 182/183 laws):
Theta1 = sp.simplify(
    -(chi0s * (1 + chi0s) * deltaU1_indep + deltaUs * (1 + deltaUs) * chi1_indep)
    / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs))
)
Xi1 = sp.simplify(
    Sigma_Z + 2 * chi0s / (1 + chi0s) * Sigma_chi + E_star * Sigma_eps
    - 4 * epsWs * deltaUs / (11 * (1 - epss) * (1 + deltaUs)**2) * Sigma_delta
)
expect_zero("Theta_1 independent slippage law", Theta1 - (-C_tr_star * Sigma_tr))
expect_zero("Xi_1 independent slippage law",   Xi1 - (A_tr_star * Sigma_tr + Sigma_nt))
expect_zero("Theta_1 monomial law", Theta1 - (-C_tr_star * Sigma_tr_compiled))
expect_zero("Xi_1 monomial law",   Xi1 - (A_tr_star * Sigma_tr_compiled + Sigma_nt_compiled))
```
**Required change (Mathematica):** mirror exactly — delete the two `coefficient form` anchors,
keep `cTrStar`/`aTrStar` typed, redefine `theta1`/`xi1` as the independent slippage-law forms
(`chi1Indep = chi0s*sigmaChi`, `deltaU1Indep = deltaUs*sigmaDelta`, then the same closed forms with
`FullSimplify[..., Assumptions -> $Assumptions]`), and add the four `expectZero` checks
(`Theta_1 independent slippage law`, `Xi_1 independent slippage law`, and the existing two monomial
laws against `sigmaTrCompiled`/`sigmaNtCompiled`).

**Use the script's ACTUAL in-scope symbol names** for the slippage drifts and parameters
(`Sigma_chi`/`sigmaChi`, `Sigma_delta`/`sigmaDelta`, `Sigma_Z`/`sigmaZ`, `Sigma_eps`/`sigmaEps`,
`E_star`/`eStar`, `epsWs`/`epss` or whatever this script names the loaded-eps and its zeroth value).
Do NOT introduce new free symbols. **Keep F2 (det) and F3 (banner) exactly as-is.**

**Self-test (REQUIRED before finalizing):** confirm (1) the independent `Theta1`/`Xi1` do NOT
reference `C_tr_star`/`A_tr_star` (else still circular); (2) the four new checks reduce to exact 0
(they are linear rational simplifications — should close fast; do NOT raise the 600s cap); (3) reason
that a DELIBERATELY wrong `C_tr_star` (e.g. wrong numerator power) leaves `(C_bad − C_true)*Sigma_tr
≠ 0` and FAILS — do not commit a broken variant, just confirm the fail mode. If the independent
`Xi1` closed form does not reduce against `(A_tr_star*Sigma_tr + Sigma_nt)`, the typed `Xi1` form is
the suspect — re-derive from the slippage laws rather than tweaking until it passes.

**Verification:** both engines show `Theta_1 independent slippage law = 0`, `Xi_1 independent
slippage law = 0`, and the two monomial-law checks `= 0`, all PASS; scripts exit 0; F2 det check and
the banner unchanged.

## Applied: Iteration 2 F1

- files_changed:
  - `scripts/moving_throat_pde_stage185_microscopic_monomials_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage185_microscopic_monomials_mathematica_audit.wl`
- summary: Replaced the circular coefficient anchors with independent slippage-law Theta_1/Xi_1 constructions and load-bearing independent-law checks.
- deviation: none
