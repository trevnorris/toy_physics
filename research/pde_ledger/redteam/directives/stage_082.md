---
unit_id: 082
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 082

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — paper_misalignment (resolved → Codex apply direction (a))

**Subtype:** script_missing_paper_claim

**Paper side:**
- `paper/stages/stage_082.tex:21-25` quote:
  > `zeta_phys(Pe, eta; kappa) = Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y(eta)^2)`

**User direction (approved 2026-05-27):** (a) extend both stage 082 scripts to pin the closed form and add a Family-1 numerical leg matching stage 084 .wl's existing instantiation.

**Required change:**

Stage 084 .wl already implements the closed-form pin (`mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:43`) and the Family-1 Pe→∞ limit check (lines 73-76). The directive for stage 082 is to add a structurally similar block at the end of both stage 082 scripts, AFTER the F3 edits below have been applied.

In `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`, after the F3 block but before the final ledger prints, insert:

```python
# F1 closed-form pin (paper eq. app-stage082-zeta-phys, user direction (a)):
# zeta_phys(Pe, eta; kappa) = Omega_Pe(Pe)^2 * (kappa + pi^2/4) / (kappa + y(eta)^2)
# with Omega_Pe(Pe) = pi*Pe*(2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2)*(exp(Pe) - 1))
# and y(eta) the smallest positive root of y*tan(y) = eta.
# Verify by reproducing the Pe->oo limit at the Family-1 point (eta, kappa) = (37, 12321/5),
# which equals stage 084's zeta_max^(F1) constant.
Pe, kappa_sym, eta_sym, y_sym = sp.symbols("Pe kappa eta y", positive=True, real=True)
Omega_Pe = sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi) / (
    (4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1)
)
zeta_phys_closed = Omega_Pe**2 * (kappa_sym + sp.pi**2 / 4) / (kappa_sym + y_sym**2)
# Pe->oo: Omega_Pe -> pi/2, so zeta_phys_closed -> (pi^2/4) * (kappa + pi^2/4) / (kappa + y^2)
Omega_Pe_limit = sp.limit(Omega_Pe, Pe, sp.oo)
print(f"\nOmega_Pe -> {Omega_Pe_limit} as Pe -> oo")
expect_zero("Omega_Pe -> pi/2 as Pe -> oo", Omega_Pe_limit - sp.pi / 2)

# Family-1 instantiation: solve y*tan(y) = 37 for the smallest positive root in (0, pi/2).
y_F1 = sp.nsolve(y_sym * sp.tan(y_sym) - 37, y_sym, sp.Float("1.527"), prec=30)
print(f"y_F1 (root of y tan y = 37, smallest positive) = {y_F1}")

kappa_F1 = sp.Rational(12321, 5)
zeta_phys_F1_limit = sp.simplify(
    (sp.pi**2 / 4) * (kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2)
)
print(f"zeta_phys(Pe->oo, kappa_F1, y_F1) = {sp.N(zeta_phys_F1_limit, 20)}")

# Compare against the Family-1 limit known from upstream (stages 080/084).
# zeta_max^(F1) is the Pe->oo support-demand limit at Family-1 (per paper appendix part 03).
# Cross-check value from mathematica/stage084 line 75-76: matches upstream constant to 10^-10.
# We carry it forward here with provenance.
# Source: stages 080 and 084 audits, both verified at v2.
zeta_max_F1_reference = sp.Rational(
    "246752922945583500", "100000000000000000"
)  # ≈ 2.467529229455835, 16-digit anchor from stage 084 .wl
# Verify the new closed-form pin matches the upstream constant to 10 digits.
diff_to_reference = abs(sp.N(zeta_phys_F1_limit - zeta_max_F1_reference, 30))
print(f"|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = {diff_to_reference}")
assert diff_to_reference < sp.Float("1e-10"), (
    f"Family-1 closed-form pin disagrees with upstream zeta_max^(F1) by {diff_to_reference}"
)
print("PASS: Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10.")
```

In `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, after the F3 block but before the final ledger prints, insert the parallel block:

```mathematica
(* F1 closed-form pin (paper eq. app-stage082-zeta-phys, user direction (a)):
   zeta_phys(Pe, eta; kappa) = Omega_Pe(Pe)^2 * (kappa + Pi^2/4) / (kappa + y(eta)^2)
   with Omega_Pe(Pe) = Pi*Pe*(2*Pe*Exp[Pe] + Pi) / ((4*Pe^2 + Pi^2)*(Exp[Pe] - 1))
   and y(eta) the smallest positive root of y*Tan[y] = eta.
   Verify by reproducing the Pe->Infinity limit at the Family-1 point (eta, kappa) = (37, 12321/5),
   matching stage 084's zeta_max^(F1) constant. *)
Module[{pe, kappaSym, etaSym, ySym, OmegaPe, zetaPhysClosed, OmegaPeLimit,
         yF1, kappaF1, zetaPhysF1Limit, zetaMaxF1Reference, diffToReference},
  ClearAll[pe, kappaSym, etaSym, ySym];
  OmegaPe = Pi*pe*(2*pe*Exp[pe] + Pi) / ((4*pe^2 + Pi^2)*(Exp[pe] - 1));
  zetaPhysClosed = OmegaPe^2 * (kappaSym + Pi^2/4) / (kappaSym + ySym^2);
  OmegaPeLimit = Limit[OmegaPe, pe -> Infinity];
  Print["Omega_Pe -> ", fmt[OmegaPeLimit], " as Pe -> oo"];
  expectZero["Omega_Pe -> Pi/2 as Pe -> oo", OmegaPeLimit - Pi/2];
  yF1 = ySym /. FindRoot[ySym*Tan[ySym] - 37 == 0, {ySym, 1.527}, WorkingPrecision -> 30];
  Print["y_F1 (root of y tan y = 37, smallest positive) = ", fmt[yF1]];
  kappaF1 = 12321/5;
  zetaPhysF1Limit = FullSimplify[(Pi^2/4) * (kappaF1 + Pi^2/4) / (kappaF1 + yF1^2)];
  Print["zeta_phys(Pe->oo, kappa_F1, y_F1) = ", fmt[N[zetaPhysF1Limit, 20]]];
  (* Reference constant from stage 084 .wl (verified at v2): zeta_max^(F1) ≈ 2.4675292294558... *)
  zetaMaxF1Reference = ToExpression["2.467529229455835`30"];
  diffToReference = Abs[N[zetaPhysF1Limit - zetaMaxF1Reference, 30]];
  Print["|zeta_phys(F1, Pe->oo) - zeta_max^(F1)| = ", fmt[diffToReference]];
  If[TrueQ[diffToReference < 10^-10],
    pass["Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10"],
    fail["Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10", diffToReference]];
];
```

**Verification:**
After Codex applies, both transcripts should print:
- `Omega_Pe -> pi/2 as Pe -> oo` with PASS.
- A printed `y_F1 ~ 1.5269...` (smallest positive root of `y tan y = 37` in `(0, pi/2)`).
- `zeta_phys(Pe->oo, kappa_F1, y_F1) ~ 2.4675292294558...`
- PASS for `Family-1 closed-form pin matches upstream zeta_max^(F1) to 10^-10`.
Both scripts must still exit 0.

## F2 — paper_misalignment (resolved → subsumed by F1's Family-1 pin)

**Subtype:** script_missing_paper_claim

**Paper side:**
- `paper/stages/stage_082.tex:43-47` quote:
  > "The Family--1 specialization is obtained by setting `(eta, kappa) = (37, 12321/5)` and `Xi_F1 = W_wall = 1369*Upsilon_w = 136900*Theta_w`."

**User direction (approved 2026-05-27):**

Sub-question (i) — Family-1 numerical pair `(eta, kappa) = (37, 12321/5)` instantiation: **subsumed by F1**. The new F1 closed-form pin block instantiates `(eta_F1, kappa_F1) = (37, 12321/5)` at the Pe->oo limit and verifies against upstream `zeta_max^(F1)`. No additional script-side change is needed beyond the F1 block above.

Sub-question (ii) — Upsilon_w convention (`117` vs `100`): resolved by the Q2 (stage 075) user direction (a) on 2026-05-27. Paper `stage_075.tex:7,24` and notes `stage075:108,116,124-128` were updated by the orchestrator to `100 Theta_w`; the script's `Upsilon_w = 100 Theta_w` is now consistent with all paper-side documents. No script-side change needed; F2(ii) is closed.

**Required change:** none beyond F1. Codex must append the standard `## Applied: F2` block noting that F2 is subsumed by F1 (no separate edit).

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py:77-87`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl:71-81`

**Issue:**

The derivative checks `dR_quad/dzeta_phys + 1 == 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr == 0` are identically true by construction of `R_quad = zeta_req - zeta_phys` (where `zeta_phys` is a free symbol independent of `Pi_tr`). They verify only that the symbolic differentiation engine works, not any physical content. Replace them with a check that actually exercises the strict-positivity content notes section 4 relies on, namely `dzeta_req/dPi_tr > 0` on the physical branch — equivalently, that `dzeta_req/dPi_tr` simplifies to a positive ratio under the existing positivity assumptions.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py`, replace lines 77-87 (the block from the comment `# Directional content of R_quad: ...` through the second `expect_zero` call) with:

```python
# Directional content of zeta_req: the inverse-map theorem (notes section 4)
# relies on zeta_req being strictly increasing in Pi_tr on the physical branch
# (where Pi_tr, C_mix, eps_blk are positive). Verify by checking that
# d zeta_req / d Pi_tr simplifies to a positive ratio under those assumptions.
dzeta_req_dPi = sp.simplify(sp.diff(zeta_req, Pi_tr))
print(f"\ndzeta_req/dPi_tr = {dzeta_req_dPi}")
# Numerator and denominator should be sign-controlled: numerator equals
# C_mix*(1 - eps_blk) (positive for eps_blk in [0,1)), denominator is a square
# (always positive on the branch). Confirm the numerator and the squared form
# of the denominator are exposed by the simplification.
num, den = sp.fraction(sp.together(dzeta_req_dPi))
expect_zero(
    "numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk)",
    sp.simplify(num - Cmix * (1 - eps_blk)),
)
# Denominator is (C_mix - eps_blk*(2*C_mix - Pi_tr))^2, identically positive
# on the branch domain because it is a square of the same denominator that
# appears in zeta_req itself.
expect_zero(
    "denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))**2",
    sp.simplify(den - (Cmix - eps_blk * (2 * Cmix - Pi_tr)) ** 2),
)
```

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl`, replace lines 71-81 (the block from the comment `(* Directional content of R_quad: ...` through the second `expectZero[...]` call) with:

```mathematica
(* Directional content of zeta_req: the inverse-map theorem (notes section 4)
   relies on zeta_req being strictly increasing in Pi_tr on the physical branch
   (where PiTr, cMix, epsBlk are positive). Verify by checking that
   d zeta_req / d Pi_tr simplifies to a positive ratio under those assumptions. *)
dZetaReqDPi = FullSimplify[D[zetaReq, PiTr], Assumptions -> $Assumptions];
Print["dzeta_req/dPi_tr = ", fmt[dZetaReqDPi]];
(* Numerator and denominator should be sign-controlled. *)
{numD, denD} = {Numerator[Together[dZetaReqDPi]], Denominator[Together[dZetaReqDPi]]};
expectZero[
  "numerator(d zeta_req/d Pi_tr) - C_mix*(1 - eps_blk)",
  FullSimplify[numD - cMix*(1 - epsBlk), Assumptions -> $Assumptions]
];
expectZero[
  "denominator(d zeta_req/d Pi_tr) - (C_mix - eps_blk*(2*C_mix - Pi_tr))^2",
  FullSimplify[denD - (cMix - epsBlk*(2*cMix - PiTr))^2, Assumptions -> $Assumptions]
];
```

**Self-test (auditor performed before writing this directive):**

1. **Variable independence:** `zeta_req` is a function of `Pi_tr`, `C_mix`, `eps_blk` (sympy line 38; mathematica derived via Solve at line 38). The new `sp.diff(zeta_req, Pi_tr)` / `D[zetaReq, PiTr]` differentiates over a real dependency, not over a free symbol. Verified non-trivial.
2. **Numerator factorization:** By quotient rule on `f/g` with `f = Pi_tr - C_mix`, `g = C_mix - eps_blk*(2*C_mix - Pi_tr)`: `f' = 1`, `g' = eps_blk` (w.r.t. Pi_tr). So `(f'g - fg')/g^2 = (g - eps_blk*(Pi_tr - C_mix))/g^2`. Expanding the numerator: `g - eps_blk*f = C_mix - eps_blk*(2*C_mix - Pi_tr) - eps_blk*(Pi_tr - C_mix) = C_mix - 2*eps_blk*C_mix + eps_blk*Pi_tr - eps_blk*Pi_tr + eps_blk*C_mix = C_mix - eps_blk*C_mix = C_mix*(1 - eps_blk)`. Numerator simplifies to `C_mix*(1 - eps_blk)`. Verified.
3. **Denominator:** Squared form of `g` is `(C_mix - eps_blk*(2*C_mix - Pi_tr))^2`. Verified.
4. **Trivial-case pre-check:** At eps_blk = 0, `dzeta_req/dPi_tr = (1*C_mix - (Pi_tr - C_mix)*0)/C_mix^2 = 1/C_mix > 0`. The numerator simplifies to `C_mix*(1-0) = C_mix`, the denominator to `C_mix^2`. Both checks pass at eps_blk = 0.
5. **Path specifications:** SymPy edit targets `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py` (correct directory). Mathematica edit targets `mathematica/moving_throat_pde_stage082_master_quadrupole_residual_mathematica_audit.wl` (correct directory). Both files exist.
6. **Paper round-trip:** The new checks introduce no new paper-side constants; they only factorize `dzeta_req/dPi_tr` into its numerator (`C_mix*(1-eps_blk)`) and denominator (`(C_mix - eps_blk*(2*C_mix - Pi_tr))^2`), both of which follow algebraically from eq. 082-zeta-req. No new paper-misalignment introduced.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 082` and `redteam exec-mathematica 082` and confirm:
- The new prints `"dzeta_req/dPi_tr = ..."` appear in both transcripts.
- The two new `expect_zero` / `expectZero` checks (numerator factorization and denominator-squared form) both produce `0` and both engines report PASS.
- Both scripts exit 0.
