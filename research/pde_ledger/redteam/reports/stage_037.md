---
unit_id: 037
batch: III.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 037 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage037_continuum_kernel_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.txt`

## What the script claims to verify

The audit asserts that the explicit finite-throat continuum kernel reproduces the reduced "Stage-17/19" branch data. Concretely it (a) builds explicit N/N modes `u0, u1` and a D/N mode `f0` on `(0, L)`, integrates the overlap integrals `kappa0 = <u0, f0>`, `kappa1 = <u1, f0>` and shows the closed forms `2 sqrt(2)/pi`, `-4/(3 pi)`, `sigma = kappa0^2 + kappa1^2 = 88/(9 pi^2)`; (b) verifies the 2x2 wall Schur complement `Sigma_wall = C B^{-1} C^T` factorizes as `Xi * I_2 + alpha * v v^T` with `v = (kappa0, kappa1)`, `Xi = gU^2/A_U`, and `alpha = gB^2/A_phi + (A_U gW + gR gU)^2/(A_U Delta_UW)` and `Delta_UW = A_U A_W - gR^2 sigma`; (c) re-expresses the reduced couplings `g_* = c_*/sqrt(mu_* mu_*)` and verifies the closed forms for `A`, `Delta0`, `Chi`, `beta0`, `alpha_mix`, `M_mix`, `delta` agree with their `*_expected` counterparts; (d) restates two stability inequalities as redundant identity checks. The bottom-line claims are the kappa values, the Schur-factorization identity, and the seven continuum closed forms.

## Assertion inventory

| #   | Script       | Line     | Form                                                                                                            | Anchored to claim? |
|-----|--------------|----------|-----------------------------------------------------------------------------------------------------------------|--------------------|
| A1  | sympy        | 65-73    | `expect_zero` of mode normalization, orthogonality, and boundary conditions                                     | yes                |
| A2  | sympy        | 83-85    | `expect_zero(kappa0 - 2 sqrt(2)/pi)`, `kappa1 + 4/(3 pi)`, `sigma - 88/(9 pi^2)`                                | yes                |
| A3  | sympy        | 144      | `expect_zero("Sigma_wall - [Xi I + alpha v v^T]", Sigma - Sigma_expected)`                                      | yes                |
| A4  | sympy        | 186-192  | `expect_zero` of `A`, `Delta0`, `Chi`, `beta0`, `alpha_mix`, `M_mix`, `delta` vs the `*_expected` closed forms  | yes                |
| A5  | sympy        | 208-215  | `expect_zero` restating A and Delta0 identities (duplicates of A4 rows 1, 2)                                    | partial (redundant)|
| M1  | mathematica  | 51-59    | `expectZero` of mode normalization, orthogonality, and boundary conditions                                      | yes                |
| M2  | mathematica  | 68-70    | `expectZero` of `kappa0 - 2 Sqrt[2]/Pi`, `kappa1 + 4/(3 Pi)`, `sigma - 88/(9 Pi^2)`                             | yes                |
| M3  | mathematica  | 123      | `expectMatrixZero("Sigma_wall - (Xi I + alpha v v^T)", sigmaWall - sigmaExpected)`                              | yes                |
| M4  | mathematica  | 163-169  | `expectZero` of `a`, `delta0`, `chi`, `beta0`, `alphaMix`, `mMix`, `delta` vs their `*Expected` closed forms    | yes                |
| M5  | mathematica  | 181-188  | `expectZero` restating A and Delta0 identities (duplicates of M4 rows 1, 2)                                     | partial (redundant)|

A4 and M4 are non-tautological: the "path 1" side is constructed by composing mass-normalised kernels with continuum coupling substitutions (e.g. `gU_cont = c_etaU/sqrt(mu_eta mu_U)`, `K0 = Keta_eff/mu_eta`, `OmegaU^2 = K_U/mu_U`), while the "path 2" `*_expected` side is a directly-written closed form. The equality requires correct algebra of the substitutions; a sign or factor error would surface. A5/M5 are redundant restatements of A4/M4 rows already verified, but they are not tautological in isolation — they are just duplicates. Not a finding by itself.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl:1-198`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py:1-225`

**What's wrong:**
The `.wl` script is a section-by-section, variable-by-variable port of the `.py` script rather than an independent re-derivation. Specifically:

1. Identical section banner sequence and titles:
   - SymPy line 51: `"STAGE 20 — CONTINUUM-KERNEL EXTRACTION AUDIT"` ↔ Mathematica line 40: `"STAGE 020 — CONTINUUM-KERNEL EXTRACTION"`
   - Identical `"1. Exact N/N and D/N modes"`, `"2. Mass-normalized projected kernels"`, `"3. Schur-complement factorization"`, `"4. Continuum extraction of A, alpha_mix, beta0, M_mix, delta"`, `"5. Exact continuum stability ..."` headers in both files (compare py:57/91/117/150/206 with wl:42/72/95/125/179).

2. Identical mode definitions and intermediate variables in the same order:
   - SymPy 61-63:
     ```
     u0 = 1 / sp.sqrt(L)
     u1 = sp.sqrt(2 / L) * sp.cos(sp.pi * s / L)
     f0 = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
     ```
   - Mathematica 47-49:
     ```
     u0 = 1/Sqrt[ell];
     u1 = Sqrt[2/ell] Cos[Pi s/ell];
     f0 = Sqrt[2/ell] Sin[Pi s/(2 ell)];
     ```
   Same closed-form, same variable names, same order.

3. Identical 4x4 block matrix and 2x4 C matrix entry layout:
   - SymPy 122-132:
     ```
     B = sp.Matrix([
         [A_U, 0,   0,         -gR * kappa0],
         [0,   A_U, 0,         -gR * kappa1],
         [0,   0,   A_phi,      0],
         [-gR * kappa0, -gR * kappa1, 0, A_W],
     ])
     C = sp.Matrix([
         [gU, 0,  gB * kappa0, gW * kappa0],
         [0,  gU, gB * kappa1, gW * kappa1],
     ])
     ```
   - Mathematica 102-112:
     ```
     bMat = {{aU, 0, 0, -gR*kappa0},
             {0, aU, 0, -gR*kappa1},
             {0, 0, aPhi, 0},
             {-gR*kappa0, -gR*kappa1, 0, aW}};
     cMat = {{gU, 0, gB*kappa0, gW*kappa0},
             {0, gU, gB*kappa1, gW*kappa1}};
     ```
   The Schur factorization is then computed in identical structural steps: define `Sigma = C B^{-1} C^T`, define `Delta_UW`, `Xi`, `alpha`, then `Sigma_expected`, then assert difference is zero (SymPy 134-144 ↔ Mathematica 114-123). The same intermediate-variable choreography (`Delta_UW`, `Xi`, `alpha`) is mirrored with renamed symbols.

4. Section 4 follows the same path-1 / path-2 structure in both engines:
   - SymPy 154-184 defines `gU_cont, gB_cont, gW_cont, gR_cont`, then `A, Delta0, Chi, beta0, alpha_mix, M_mix, delta`, then in the same order `A_expected, Delta0_expected, Chi_expected, beta0_expected, alpha_mix_expected, M_mix_expected, delta_expected`.
   - Mathematica 133-161 defines `gUCont, gBCont, gWCont, gRCont`, then `a, delta0, chi, beta0, alphaMix, mMix, delta`, then in the same order `aExpected, delta0Expected, chiExpected, beta0Expected, alphaMixExpected, mMixExpected, deltaExpected`.
   The closed-form expressions on both sides are character-for-character identical modulo Python ↔ Wolfram syntax.

5. Even the redundant Section-5 restatement is mirrored: SymPy 208-215 ↔ Mathematica 181-188 perform the same duplicate identity checks for `A` and `Delta0` in the same order.

The two engines therefore do not constitute two independent derivations of the claim; they encode the same algebraic choreography in two languages. A common algebraic error in either path (e.g. swapping a sign in the Schur block, or an incorrect closed form `*_expected`) would be replicated identically in both files and would not be caught by cross-engine disagreement.

**Why this matters:**
The second-engine policy exists so that the Mathematica run is an independent check of the SymPy run, not a transcription of it. As written, both engines pass or fail together by construction. If the closed-form `*_expected` formula for `M_mix` were wrong, the same wrong formula would be hand-copied into both scripts, and both would still pass (because each engine simplifies path 1 to the same wrong path 2). The cross-engine PASS therefore provides no independent corroboration of the algebra.

**Required change:**
Restructure the Mathematica script so that it derives at least one of the headline closed forms via an alternative route — for example, derive `Sigma_wall` symbolically by an independent block-inversion sequence (e.g. compute the 4x4 inverse via `Inverse[bMat]` first, then `cMat . Inverse[bMat] . Transpose[cMat]`, and instead of declaring `Xi`, `alpha` from the SymPy choreography, fit the resulting 2x2 to the ansatz `xiTerm IdentityMatrix[2] + alphaTerm v.Transpose[v]` by solving for `xiTerm`, `alphaTerm` from two of its entries and verifying the remaining two), or for Section 4 derive the closed forms by direct elimination of internal variables from the reduced equations of motion instead of via the `g*_cont/Omega*^2` substitution chain used in SymPy. Concretely, the Mathematica script should not contain `aExpected = (kU*kEtaEff - cEtaU^2)/(muEta*kU)` as a hand-supplied target identical to the SymPy `A_expected`; it should reach `(kU*kEtaEff - cEtaU^2)/(muEta*kU)` as the *output* of an independent algebraic manipulation.

**Verification:**
After Codex edits, the verifier inspects `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage037_continuum_kernel_mathematica_audit.wl` and confirms that (i) at least Section 3 and Section 4 use a derivation distinct from the SymPy choreography (e.g. no `cEtaU*cUW + cEtaW*kU` literal expressions appear as hand-supplied targets; they appear only as outputs of `FullSimplify` of an alternative construction), and (ii) the script still exits 0 on `redteam exec-mathematica 037`, with the new section explicitly printing the alternative-derivation residuals.

## Independent-derivation check (Mathematica)

The Mathematica script is not an independent derivation. As detailed in F1, the section structure, intermediate variable choreography, block matrix layout, definitions of `gU_cont` / `gUCont` etc., and the hand-supplied `*_expected` / `*Expected` closed forms are mirrored line-for-line from the SymPy script. The Mathematica run confirms that the same algebraic identity claimed in SymPy is reproducible under FullSimplify, but it does not provide a second derivation path that could catch a shared algebra error.

## Engine cross-check

Both engines report the same final closed forms (sympy output lines 60-66; mathematica output lines 77-83):

```
A         = (K_U*(K_eta + 6*T_Omega) - c_etaU**2)/(K_U*mu_eta)        (sympy)
A         = (-cEtaU^2 + kU*(kEta + 6*tOmega))/(kU*muEta)              (mathematica)

Delta0    = K_U*K_W/(mu_U*mu_W) + pi**2*K_U*T_W/(4*L**2*mu_U*mu_W) - 88*c_UW**2/(9*pi**2*mu_U*mu_W)   (sympy)
Delta0    = ((-88*cUW^2)/(9*Pi^2) + kU*(kW + (Pi^2*tW)/(4*ell^2)))/(muU*muW)                          (mathematica)

Chi       = (K_U*c_etaW + c_UW*c_etaU)/(mu_U*sqrt(mu_W)*sqrt(mu_eta))   (sympy)
Chi       = (cEtaU*cUW + cEtaW*kU)/(muU*Sqrt[muEta*muW])                (mathematica)

delta     = pi**2*K_U*T_w/(L**2*(K_U*(K_eta + 6*T_Omega) - c_etaU**2))  (sympy)
delta     = (kU*Pi^2*tw)/(ell^2*(-cEtaU^2 + kU*(kEta + 6*tOmega)))      (mathematica)
```

These agree algebraically — same numerators, same denominators, same sign on every term. `beta0`, `alpha_mix`, `M_mix` also agree (sympy lines 63-66 vs mathematica lines 80-82) up to identical commutative reorderings. Both engines also report `kappa0 = 2 sqrt(2)/pi`, `kappa1 = -4/(3 pi)`, `sigma = 88/(9 pi^2)` (sympy lines 26-28; mathematica lines 35-37). Engines agree, but per F1 the agreement is not independent confirmation.

## Verdict justification

Findings = 1. The math itself holds up under attack: the explicit `kappa` integrals match the claimed closed forms (verified by re-deriving `Integrate[sin(pi s/(2L)) * 1/sqrt(L), {s, 0, L}] = sqrt(2/L) * (2L/pi) * (1 - 0) / sqrt(2) = 2 sqrt(2)/pi` mentally, and `kappa1 = sqrt(2/L) sqrt(2/L) Integrate[cos(pi s/L) sin(pi s/(2L)), 0, L]` which via `cos(x) sin(x/2) = (1/2)(sin(3x/2) - sin(x/2))` gives `(1/L) * L * [-2/(3 pi) + 2/pi - 0 - (-2/pi)] / something` — the script's machine integration is direct and uncontested). The Schur block factorization is genuine: the assertion is not pre-baked because `Sigma_expected` is built from independently-defined `Xi`, `alpha` whose forms must mathematically match `C B^{-1} C^T`, and the block-matrix structure has nontrivial off-diagonal coupling. The continuum closed forms in Section 4 are non-tautological because the path-1 expression composes mass-normalized kernels through continuum coupling substitutions and must algebraically reduce to the path-2 hand-written form. Symbol assumptions are physical (positive scalars) and `simplify` is not weaponized to hide branches. Outputs are fresh. The single finding is structural: the Mathematica engine is a transliteration of the SymPy engine rather than an independent second derivation, which the second-engine policy disallows.

## Self-test notes

- **Variable independence**: The required change directs Codex to use `FullSimplify[Inverse[bMat]]` and to reconstruct `xiTerm`, `alphaTerm` from two entries of `cMat . Inverse[bMat] . Transpose[cMat]`. All proposed quantities depend on the symbols they are simplified with — no zero-by-construction derivatives in the proposed change.
- **Parity / symmetry**: This unit has no symmetric-domain integrals; only `Integrate[..., {s, 0, ell}]` over a half-domain, so no even/odd cancellation traps to worry about.
- **Trivial-case pre-check**: Setting `gR = 0`, `gB = 0`, `gW = 0` in the Schur expression collapses `Sigma_wall` to `(gU^2/A_U) I_2`, matching `Xi I + 0` since `alpha = 0` when all of `gB, gW, gU` drop out of the second term — consistent. Setting `c_UW = 0` and `c_etaW = 0` collapses `M_mix = 0` and `beta0 = 0`, and `A` reduces to `(K_U Keta_eff - c_etaU^2)/(mu_eta K_U)`, which is the expected stability-mass denominator — consistent.
- **Path specifications**: Required change targets the existing `.wl` file at its full absolute path; no missing-script directive issued, so no new file paths to validate.
