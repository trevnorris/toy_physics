---
unit_id: 043
batch: III.1
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T18:39:33Z
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 043

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:43-125`

**Issue:** The Mathematica script is a line-by-line port of the SymPy script with snake_case → camelCase identifier renaming. Every intermediate definition (`dU`, `v`, `y`, `sigma0`, `rPhi`, `rPhiExpected`, `dPhi`, `dPhiExpected`, ...), every substitution choice (`kappa1^2 -> (2/11) sigma`, `kappa0^2 -> 8/Pi^2`), and every assertion appears in the same order with the same algebraic form. There is no independent derivation path.

**Required change:**

Edit the Mathematica script so each of the five claim sections (1–5) arrives at its result via a path that is structurally distinct from the SymPy script's algebra. Apply ALL of the following:

1. In section 1 (lines 43–59): replace the manual definition of `dPhi`
   - Before: `dPhi = FullSimplify[kappa0 y1 - kappa1 y0, Assumptions -> $Assumptions];`
   - After: `dPhi = FullSimplify[Det[{{y0, y1}, {kappa0, kappa1}}], Assumptions -> $Assumptions];`

   And replace the explicit `rPhi = FullSimplify[(y1/y0)/(kappa1/kappa0), ...]` with a residue/ratio formulation:
   - Before: `rPhi = FullSimplify[(y1/y0)/(kappa1/kappa0), Assumptions -> $Assumptions];`
   - After: `rPhi = FullSimplify[(y[[2]]/kappa1) / (y[[1]]/kappa0), Assumptions -> $Assumptions];`

2. In section 2 (lines 63–80): replace the direct substitution `sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> sigma - (2/11) sigma}` with a limit-and-coefficient check. Before line 79 (`expectZero["support overlap contraction", ...]`), insert:

   ```mathematica
   sUEndpointZero = FullSimplify[(sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> (9/11) sigma}) /. deltaU -> 0, Assumptions -> $Assumptions];
   sUEndpointZeroExpected = sigma/kU;
   sUEndpointInf = FullSimplify[Limit[sU /. {kappa1^2 -> (2/11) sigma, kappa0^2 -> (9/11) sigma}, deltaU -> Infinity], Assumptions -> $Assumptions];
   sUEndpointInfExpected = (9/11) sigma/kU;
   Print["v.D_U.v at deltaU=0 = ", fmt[sUEndpointZero]];
   Print["v.D_U.v as deltaU->Infinity = ", fmt[sUEndpointInf]];
   expectZero["overlap endpoint deltaU=0", sUEndpointZero - sUEndpointZeroExpected];
   expectZero["overlap endpoint deltaU->Infinity", sUEndpointInf - sUEndpointInfExpected];
   ```

3. In section 3 (lines 84–98): build `mSuppExpected` from a free baseline symbol `bBaseline` rather than copying the SymPy quotient form. Apply the F3 change there (see below).

4. In section 4 (lines 102–114): replace the explicit `dPhiZ = FullSimplify[y0 z1 - y1 z0, ...]` with a determinant form:
   - Before: `dPhiZ = FullSimplify[y0 z1 - y1 z0, Assumptions -> $Assumptions];`
   - After: `dPhiZ = FullSimplify[Det[{{y0, y1}, {z0, z1}}], Assumptions -> $Assumptions];`

5. In section 5 (lines 118–125): add a Series-expansion sanity check on the mismatch BEFORE the final closed-form assertion. Insert immediately after line 122 (after the `mismatchExpected` definition, before the `Print["R_phi - R_U..."]` line):

   ```mathematica
   mismatchLeading = FullSimplify[Series[mismatch, {deltaU, 0, 1}] // Normal, Assumptions -> $Assumptions];
   mismatchLeadingExpected = FullSimplify[deltaU (rho0 - sigma0)/((1 + rho0) (1 + sigma0)), Assumptions -> $Assumptions];
   Print["mismatch leading in deltaU = ", fmt[mismatchLeading]];
   expectZero["mismatch leading-in-deltaU coefficient", mismatchLeading - mismatchLeadingExpected];
   ```

   This adds a genuinely new algebraic test (series expansion vs. coefficient extraction) that the SymPy script does not perform.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 043` and confirm the new output `.txt` contains all of:
- `PASS: overlap endpoint deltaU=0`
- `PASS: overlap endpoint deltaU->Infinity`
- `PASS: mismatch leading-in-deltaU coefficient`

AND the script must exit 0. The script must NOT contain the literal expressions `kappa0 y1 - kappa1 y0` or `y0 z1 - y1 z0`; they must be expressed via `Det[...]`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Reworked the Mathematica support-direction checks to use determinant, endpoint, free-baseline, and series-expansion witnesses.
- deviation: Adjusted `dPhiExpected` to the sign implied by the required determinant orientation.

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:96-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:69-80`

**Issue:** The `A_phi^(eff) - expected` check is tautological because both sides are constructed from the same `SU_expected`: `Aphi_eff = Kphi_eff - cUphi^2 * SU_expected` and `Aphi_eff_expected = Kphi_eff * (1 - eps_phi_split)` with `eps_phi -> cUphi^2 sigma/(KU Kphi_eff)` reduces to the same `Kphi_eff - cUphi^2 * SU_expected`. The residual is identically zero by construction.

**Required change:**

(a) In the SymPy script at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`, insert the following block AFTER the existing line 101 (`expect_zero("A_phi^(eff) - expected", Aphi_eff - Aphi_eff_expected)`), BEFORE `subbanner("26.3 — Exact physical support baseline")` on line 103:

```python
# F2: minimal-overlap baseline (delta_U = 0) reduces to the trivial pole shift.
Aphi_eff_min = sp.simplify(Aphi_eff.subs(delta_U, 0))
Aphi_eff_min_expected = sp.simplify(Kphi_eff - cUphi**2 * sigma / KU)
print("A_phi^(eff) at delta_U=0 =")
sp.pprint(sp.factor(Aphi_eff_min))
expect_zero("A_phi^(eff) at delta_U=0 (minimal)", Aphi_eff_min - Aphi_eff_min_expected)

# F2: the split-vs-minimal overlap ratio is exactly (1 - (2/11) delta_U/(1+delta_U)).
overlap_ratio = sp.simplify((Kphi_eff - Aphi_eff) / (Kphi_eff - Aphi_eff_min))
overlap_ratio_expected = sp.simplify(1 - sp.Rational(2, 11) * delta_U / (1 + delta_U))
print("split-vs-minimal overlap ratio =")
sp.pprint(sp.factor(overlap_ratio))
expect_zero("split-vs-minimal overlap ratio", overlap_ratio - overlap_ratio_expected)
```

(b) In the Mathematica script at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`, insert the following AFTER existing line 80 (`expectZero["A_phi^(eff) - expected", aPhiEff - aPhiEffExpected];`), BEFORE `subbanner["3. Exact physical support baseline"];` on line 82:

```mathematica
aPhiEffMin = FullSimplify[Limit[aPhiEff, deltaU -> 0], Assumptions -> $Assumptions];
aPhiEffMinExpected = FullSimplify[kPhiEff - cUphi^2 sigma/kU, Assumptions -> $Assumptions];
Print["A_phi^(eff) at deltaU=0 = ", fmt[aPhiEffMin]];
expectZero["A_phi^(eff) at deltaU=0 (minimal)", aPhiEffMin - aPhiEffMinExpected];

overlapRatio = FullSimplify[(kPhiEff - aPhiEff)/(kPhiEff - aPhiEffMin), Assumptions -> $Assumptions];
overlapRatioExpected = FullSimplify[1 - (2/11) deltaU/(1 + deltaU), Assumptions -> $Assumptions];
Print["split-vs-minimal overlap ratio = ", fmt[overlapRatio]];
expectZero["split-vs-minimal overlap ratio", overlapRatio - overlapRatioExpected];
```

Note: the Mathematica side uses `Limit[..., deltaU -> 0]` rather than `/. deltaU -> 0` to be structurally distinct from the SymPy `subs(delta_U, 0)` (independent-derivation principle).

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 043` and `redteam exec-mathematica 043` and confirm:
- SymPy output contains: `A_phi^(eff) at delta_U=0 (minimal) = 0` AND `split-vs-minimal overlap ratio = 0`.
- Mathematica output contains: `PASS: A_phi^(eff) at deltaU=0 (minimal)` AND `PASS: split-vs-minimal overlap ratio`.
- Both scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Added minimal-overlap and split-vs-minimal ratio checks for the effective phi pole shift in both scripts.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:107-120`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:84-98`

**Issue:** The `M_supp - expected` check is structurally tautological. `Msupp_cont` is built as `(numer)/(denom)` where `mu_eta * mu_phi` appears in both, cancelling; `Msupp_expected` is built as the post-cancellation form; the hand-substitution `kappa0^2 -> 8/pi^2` is applied to both sides identically. The check verifies `(a/b)*(b/a)*X == X`, not the physics.

**Required change:**

(a) In the SymPy script at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`, REPLACE lines 107-120 (the entire `Msupp_cont` / `Msupp_expected` / `Msupp_cont_eval` block) with the following structure (still ends at `subbanner("26.4 — Exact tracking condition relative to the mixed vector")`):

Before (lines 107–120):
```python
Msupp_cont = sp.simplify(
    (kappa0**2 * cB**2 * (1 + ceU * cUphi / (KU * cB))**2 / (mu_eta * mu_phi))
    / ((Keta_eff * (1 - eps_eta) / mu_eta) * (Kphi_eff * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))) / mu_phi))
)
Msupp_expected = sp.simplify(
    (sp.Rational(8, 1) / sp.pi**2)
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
Msupp_cont_eval = sp.simplify(Msupp_cont.subs(kappa0**2, sp.Rational(8, 1) / sp.pi**2))
print("M_supp =")
sp.pprint(sp.factor(Msupp_cont_eval))
expect_zero("M_supp - expected", Msupp_cont_eval - Msupp_expected)
```

After:
```python
Msupp_cont = sp.simplify(
    (kappa0**2 * cB**2 * (1 + ceU * cUphi / (KU * cB))**2 / (mu_eta * mu_phi))
    / ((Keta_eff * (1 - eps_eta) / mu_eta) * (Kphi_eff * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))) / mu_phi))
)

# F3: M_supp must not depend on the bare-mode masses mu_eta, mu_phi (they must cancel).
expect_zero("M_supp independent of mu_eta", sp.simplify(sp.diff(Msupp_cont, mu_eta)))
expect_zero("M_supp independent of mu_phi", sp.simplify(sp.diff(Msupp_cont, mu_phi)))

# F3: structural form check with a FREE baseline B (no hand-substitution of kappa0^2).
B = sp.symbols("B_baseline", positive=True, real=True)
Msupp_cont_in_B = sp.simplify(Msupp_cont.subs(kappa0**2, B))
Msupp_struct_expected = sp.simplify(
    B
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
print("M_supp structural form (free baseline B) =")
sp.pprint(sp.factor(Msupp_cont_in_B))
expect_zero("M_supp structural form (free baseline)", Msupp_cont_in_B - Msupp_struct_expected)

# F3: baseline value identification, isolated from the structural check.
Msupp_cont_eval = sp.simplify(Msupp_cont_in_B.subs(B, sp.Rational(8, 1) / sp.pi**2))
Msupp_expected = sp.simplify(Msupp_struct_expected.subs(B, sp.Rational(8, 1) / sp.pi**2))
print("M_supp at baseline B = 8/pi^2 =")
sp.pprint(sp.factor(Msupp_cont_eval))
expect_zero("M_supp at baseline B = 8/pi^2", Msupp_cont_eval - Msupp_expected)
```

(b) In the Mathematica script at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`, REPLACE lines 84-98 (the entire `epsPhiSplitPhys` / `mSuppCont` / `mSuppExpected` / `mSuppContEval` block, up to and including `expectZero["M_supp - expected", ...]`) with:

Before (lines 84–98):
```mathematica
epsPhiSplitPhys = FullSimplify[(epsPhiSplit /. epsPhi -> cUphi^2 sigma/(kU kPhiEff)), Assumptions -> $Assumptions];
mSuppCont = FullSimplify[
  (kappa0^2 cB^2 (1 + cEtaU cUphi/(kU cB))^2/(muEta muPhi))/
    ((kEtaEff (1 - epsEta)/muEta) (kPhiEff (1 - epsPhiSplitPhys)/muPhi)),
  Assumptions -> $Assumptions
];
mSuppExpected = FullSimplify[
  (8/Pi^2) (cB^2/(kEtaEff kPhiEff)) (1 + cEtaU cUphi/(kU cB))^2/
    ((1 - epsEta) (1 - epsPhiSplitPhys)),
  Assumptions -> $Assumptions
];
mSuppContEval = FullSimplify[mSuppCont /. kappa0^2 -> 8/Pi^2, Assumptions -> $Assumptions];

Print["M_supp = ", fmt[mSuppContEval]];
expectZero["M_supp - expected", mSuppContEval - mSuppExpected];
```

After:
```mathematica
epsPhiSplitPhys = FullSimplify[(epsPhiSplit /. epsPhi -> cUphi^2 sigma/(kU kPhiEff)), Assumptions -> $Assumptions];
mSuppCont = FullSimplify[
  (kappa0^2 cB^2 (1 + cEtaU cUphi/(kU cB))^2/(muEta muPhi))/
    ((kEtaEff (1 - epsEta)/muEta) (kPhiEff (1 - epsPhiSplitPhys)/muPhi)),
  Assumptions -> $Assumptions
];

(* F3: M_supp must not depend on muEta or muPhi (they must cancel). *)
expectZero["M_supp independent of muEta", FullSimplify[D[mSuppCont, muEta], Assumptions -> $Assumptions]];
expectZero["M_supp independent of muPhi", FullSimplify[D[mSuppCont, muPhi], Assumptions -> $Assumptions]];

(* F3: structural-form check with a free baseline symbol bBaseline. *)
Clear[bBaseline];
$Assumptions = $Assumptions && bBaseline > 0;
mSuppContInB = FullSimplify[mSuppCont /. kappa0^2 -> bBaseline, Assumptions -> $Assumptions];
mSuppStructExpected = FullSimplify[
  bBaseline (cB^2/(kEtaEff kPhiEff)) (1 + cEtaU cUphi/(kU cB))^2/
    ((1 - epsEta) (1 - epsPhiSplitPhys)),
  Assumptions -> $Assumptions
];
Print["M_supp structural form (free baseline) = ", fmt[mSuppContInB]];
expectZero["M_supp structural form (free baseline)", mSuppContInB - mSuppStructExpected];

(* F3: baseline value identification, isolated from the structural check. *)
mSuppContEval = FullSimplify[mSuppContInB /. bBaseline -> 8/Pi^2, Assumptions -> $Assumptions];
mSuppExpected = FullSimplify[mSuppStructExpected /. bBaseline -> 8/Pi^2, Assumptions -> $Assumptions];
Print["M_supp at baseline B = 8/Pi^2 = ", fmt[mSuppContEval]];
expectZero["M_supp at baseline B = 8/Pi^2", mSuppContEval - mSuppExpected];
```

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 043` and `redteam exec-mathematica 043` and confirm:
- SymPy output contains all of: `M_supp independent of mu_eta = 0`, `M_supp independent of mu_phi = 0`, `M_supp structural form (free baseline) = 0`, `M_supp at baseline B = 8/pi^2 = 0`.
- Mathematica output contains all of: `PASS: M_supp independent of muEta`, `PASS: M_supp independent of muPhi`, `PASS: M_supp structural form (free baseline)`, `PASS: M_supp at baseline B = 8/Pi^2`.
- Both scripts exit 0.
- The original `M_supp - expected` line is no longer present (replaced by the four new checks).

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- summary: Replaced the support-baseline comparison with mass-independence, free-baseline structural, and isolated baseline-value checks.
- deviation: none
