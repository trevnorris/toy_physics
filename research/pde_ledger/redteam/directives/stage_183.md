---
unit_id: 183
batch: V.2
created_at: 2026-05-30T00:43:16-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-30T01:21:36-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 183

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN both affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py:106-114`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl:79-85`

**Issue:** The "Triple-rigidity theorem" block claims to verify the paper's headline equivalence `Θ₁=Ξ₁=R₁=0 ⟺ Σ_tr=Σ_nt=Σ_η=0` (notes Sec. 4.4; appendix eq. `app-part05-zero-normal-form`), but its three assertions only substitute the slippages to zero into linear forms and confirm the output is zero. That exercises only the trivial `⟸` direction, which holds for any linear map regardless of the prefactors, so the checks cannot fail. The non-trivial `⟹` direction (vanishing observables force vanishing slippages) requires the map's prefactors `C_tr`, `A_tr`, and `ε_η/(1−ε_η)` to be nonzero on the branch `χ₀>0, δ_U>0, 0<ε_η<1`; that is never tested. (The block also uses the bare symbol `SigmaNT`/`sigmaNT`, detached from the constructed `SigmaNT_def`/`sigmaNTDef`.) Replace it with a check of the genuine non-trivial content: the three prefactors are nonzero on the branch, which (together with the already-present inverse round-trips A3/A5/A6) certifies the equivalence.

**Required change (SymPy):**

1. Add an `expect_nonzero` helper next to `expect_zero` (after the `expect_zero` definition, around sympy:33). Insert:

```python
def expect_nonzero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(expr)
    print(f"{name} = {expr}")
    if expr == 0:
        raise AssertionError(f"{name} is unexpectedly zero")
```

2. Replace the entire `banner("Triple-rigidity theorem")` block, currently sympy:106–114:

```python
banner("Triple-rigidity theorem")
# The map is triangular, so vanishing observables imply vanishing adapted slippages.
# We verify the forward zero map directly.
Theta_zero = sp.simplify(Theta1.subs({SigmaTr: 0}))
Xi_zero = sp.simplify((A_tr * SigmaTr + SigmaNT).subs({SigmaTr: 0, SigmaNT: 0}))
Rsum_zero = sp.simplify((-eps_eta / (1 - eps_eta) * SigmaEta).subs({SigmaEta: 0}))
expect_zero("Theta_1|(Sigma_tr=0)", Theta_zero)
expect_zero("Xi_1|(Sigma_tr=Sigma_nt=0)", Xi_zero)
expect_zero("(R_1+Xi_1)|(Sigma_eta=0)", Rsum_zero)
```

with:

```python
banner("Triple-rigidity theorem")
# Rigidity (Theta_1=Xi_1=R_1=0 <=> Sigma_tr=Sigma_nt=Sigma_eta=0) holds iff the
# triangular map is invertible on the branch chi0>0, deltaU>0, 0<eps_eta<1, i.e.
# iff each diagonal prefactor is nonzero there. We test that non-trivial content;
# the trivial forward direction is already implied, and the inverse round-trips
# above (Sigma_tr/Sigma_nt/Sigma_eta inverse) confirm full invertibility.
dressing_pref = sp.simplify(eps_eta / (1 - eps_eta))
expect_nonzero("C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch", C_tr)
expect_nonzero("A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch", A_tr)
expect_nonzero("eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch", dressing_pref)
```

(`C_tr` and `A_tr` are already defined at sympy:78–79; `eps_eta`, `chi0`, `deltaU` carry `positive=True`, so simplification yields strictly nonzero literals — `C_tr = chi0*deltaU/((chi0+1)(deltaU+1)(chi0+deltaU+1))`, `A_tr = 2*chi0/((chi0+1)(deltaU+1))`, `eps_eta/(1-eps_eta)`.)

**Required change (Mathematica):**

1. Add an `expectNonzero` helper next to `expectZero` (after the `expectZero` definition, around wl:24). Insert:

```wolfram
expectNonzero[name_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], fail[name, res], pass[name]];
];
```

2. Replace the `banner["Triple-rigidity theorem"]` block, currently wl:79–85:

```wolfram
banner["Triple-rigidity theorem"];
thetaZero = FullSimplify[theta1 /. sigmaTr -> 0, Assumptions -> $Assumptions];
xiZero = FullSimplify[(aTr*sigmaTr + sigmaNT) /. {sigmaTr -> 0, sigmaNT -> 0}, Assumptions -> $Assumptions];
rsumZero = FullSimplify[(-epsEta*sigmaEta/(1 - epsEta)) /. sigmaEta -> 0, Assumptions -> $Assumptions];
expectZero["Theta_1|(Sigma_tr=0)", thetaZero];
expectZero["Xi_1|(Sigma_tr=Sigma_nt=0)", xiZero];
expectZero["(R_1+Xi_1)|(Sigma_eta=0)", rsumZero];
```

with:

```wolfram
banner["Triple-rigidity theorem"];
(* Rigidity holds iff the triangular map is invertible on the branch, i.e. iff
   each diagonal prefactor is nonzero there. We test that non-trivial content;
   the inverse round-trips above confirm full invertibility. *)
dressingPref = FullSimplify[epsEta/(1 - epsEta), Assumptions -> $Assumptions];
expectNonzero["C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch", cTr];
expectNonzero["A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch", aTr];
expectNonzero["eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch", dressingPref];
```

(`cTr` and `aTr` are already defined at wl:55–56; `$Assumptions` (wl:29–30) carries `chi0>0, epsW>0, epsEta>0, deltaU>0`, so the three prefactors simplify to strictly nonzero forms.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 183` and `redteam exec-mathematica 183` and confirm: (a) the "Triple-rigidity theorem" section now prints three nonzero-prefactor lines (`C_tr ... nonzero`, `A_tr ... nonzero`, `eps_eta/(1-eps_eta) ... nonzero`) instead of three lines reading `= 0`; (b) both scripts exit 0 with all checks passing.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage183_triangular_normal_form_mathematica_audit.wl`
- summary: Replaced the tautological zero-substitution rigidity checks with nonzero diagonal-prefactor checks and added matching nonzero helpers in both audit scripts.
- deviation: none
