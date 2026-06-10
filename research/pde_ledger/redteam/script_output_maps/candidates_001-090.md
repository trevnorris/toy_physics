# Phase-A classifier input — BAND 001-090 (read-only)

Adjudicate each ACTIONABLE row by CONTENT (which stage owns the cited deliverable?); NEVER offset-sweep. Ground truth = filename `stageNNN_`. Report A1/A2/A3/A4 + the per-ref owning-stage map (or LEAVE/cannot-confirm).

## SELF  (49)

| canon | kind | line | token | cited | Δ | flags | context |
| ---: | --- | ---: | --- | ---: | ---: | --- | --- |
| 002 | py | 101 | `STAGE 001` | 1 | -1 | - | banner("SECTION II — AXISYMMETRIC TWO-MODE REDUCTION FROM THE STAGE 001 ACTION") |
| 003 | py | 26 | `Stage-002` | 2 | -1 | - | Stage-002 breathing-reduction variables carried forward without renaming. |
| 022 | py | 205 | `STAGE-021` | 21 | -1 | - | banner("SECTION IV — STAGE-021 MIXED-SECTOR PROTOTYPE: EXACT N0/N2/N4") |
| 022 | py | 231 | `Stage-021` | 21 | -1 | - | subbanner("IV.2 — Dictionary back to the Stage-021 Maxwell/mixed symbols") |
| 022 | py | 274 | `Stage-021` | 21 | -1 | - | subbanner("V.1 — Stage-021 outgoing l=2 fingerprint anchor") |
| 022 | py | 333 | `stage4` | 4 | -18 | stem | stage4_prototype_transfer_coefficients() |
| 022 | py-out | 152 | `Stage-021` | 21 | -1 | OUT | Stage-021 A coefficient = 0 |
| 022 | py-out | 153 | `Stage-021` | 21 | -1 | OUT | Stage-021 B coefficient = 0 |
| 022 | py-out | 154 | `Stage-021` | 21 | -1 | OUT | Stage-021 G5 coefficient = 0 |
| 022 | wl-out | 63 | `Stage-021` | 21 | -1 | OUT | Stage-021 A coefficient = 0 |
| 022 | wl-out | 65 | `Stage-021` | 21 | -1 | OUT | Stage-021 B coefficient = 0 |
| 022 | wl-out | 67 | `Stage-021` | 21 | -1 | OUT | Stage-021 G5 coefficient = 0 |
| 023 | py-out | 177 | `Stage-022` | 22 | -1 | OUT | Stage-022 Gamma5_port anchor = 0 |
| 023 | wl | 131 | `Stage-003` | 3 | -20 | - | (* Stage-003 carry-forward: B_{A0}, B_{A2}, B_{A4} are the stable-BdG Schur |
| 023 | wl-out | 117 | `Stage-022` | 22 | -1 | OUT | Stage-022 Gamma5_port anchor = 0 |
| 024 | py | 10 | `Stage-023` | 23 | -1 | - | Stage-023 grouped bundle: |
| 024 | py | 29 | `Stage 022` | 22 | -2 | - | Stage 022 grouped prefactor law carried through the Stage 023 bundle |
| 024 | py | 29 | `Stage 023` | 23 | -1 | - | Stage 022 grouped prefactor law carried through the Stage 023 bundle |
| 026 | py | 10 | `Stage-025` | 25 | -1 | - | Stage-025 minimal isotropic normalization closure: |
| 026 | py | 120 | `STAGE-025` | 25 | -1 | - | banner("SECTION III — CONCRETE BRANCH SUBSTITUTION INTO STAGE-025 QUANTITIES") |
| 027 | py | 183 | `STAGE-025` | 25 | -2 | - | banner("SECTION IV — PROFILE-DRESSED STAGE-025/026 BRANCH") |
| 027 | py | 226 | `Stage 026` | 26 | -1 | - | subbanner("IV.2 — Exact recovery of Stage 026 at theta = 0") |
| 037 | py-out | 74 | `Stage-034` | 34 | -3 | OUT | Stage-034/036 branch data exactly. |
| 038 | py | 42 | `Stage-037` | 37 | -1 | - | subbanner("1. Stage-037 continuum formulas") |
| 040 | py | 10 | `Stage-035` | 35 | -5 | - | Stage-035 normalization function deforms to an exact two-vector shape |
| 040 | py | 101 | `Stage 22` | 22 | -18 | - | subbanner("40.3 — Specialization to the split-U continuum of Stage 22") |
| 042 | py | 53 | `Stage-041` | 41 | -1 | - | subbanner("42.1 — Exact selected-mode eigenvector ratio after inserting the Stage-041 support loading") |
| 042 | py | 103 | `Stage 040` | 40 | -2 | - | subbanner("42.3 — Tracking-support collapse back to Stage 040") |
| 045 | py | 128 | `Stage-27` | 27 | -18 | - | banner("4. Exact collapse of the Stage-27 branch equation") |
| 045 | py-out | 43 | `Stage-044` | 44 | -1 | OUT | Stage-044 tracking F collapse = 0 |
| 045 | wl-out | 36 | `Stage-044` | 44 | -1 | OUT | Stage-044 tracking F collapse = 0 |
| 051 | wl | 87 | `Stage 047` | 47 | -4 | - | (* Stage 047/030 coherent forward map: Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2]. *) |
| 069 | py-out | 37 | `Stage-066` | 66 | -3 | OUT | Stage-066 matched branch: |
| 069 | py-out | 40 | `Stage-068` | 68 | -1 | OUT | Stage-068 resonance family: |
| 069 | wl | 99 | `Stage-066` | 66 | -3 | - | (* Stage-066 matched-window derivation via parameterized effective gap. *) |
| 069 | wl | 116 | `Stage-068` | 68 | -1 | - | (* Stage-068 / Stage-068 resonance penalty via band-edge ratio extraction. *) |
| 069 | wl-out | 61 | `Stage-066` | 66 | -3 | OUT | Stage-066 matched branch: |
| 069 | wl-out | 64 | `Stage-068` | 68 | -1 | OUT | Stage-068 resonance family: |
| 070 | wl-out | 27 | `Stage-065` | 65 | -5 | OUT | Stage-065 normalization: J_1 := I_f/H_w (shell measure 4 pi a^2 ell absorbed into J_1) |
| 070 | wl-out | 28 | `Stage-064` | 64 | -6 | OUT | Stage-064 normalization: I_1 := N_phiphi/H_w = (4 pi a^2 ell I_f)/H_w |
| 087 | py | 12 | `stage 65` | 65 | -22 | - | - scripts/moving_throat_pde_stage081_*_sympy_audit.py / .wl  (a.k.a. former stage 65) |
| 087 | py | 13 | `stage 69` | 69 | -18 | - | - scripts/moving_throat_pde_stage082_*_sympy_audit.py / .wl  (former stage 69 closure) |
| 088 | wl | 86 | `Stage-085` | 85 | -3 | - | (* Stage-085 product identity: Pi_tr = rho_alpha * C_mix (verified upstream |
| 089 | py-out | 7 | `Stage-080` | 80 | -9 | OUT | Stage-080 zeta_max = A_F1 pi^2/4 = -3.2381793033374248823225545803690717748724259704392E-101 |
| 089 | py-out | 8 | `Stage-082` | 82 | -7 | OUT | Stage-082 Q(zeta;0)=1+zeta reduction = 0 |
| 089 | wl | 43 | `Stage-079` | 79 | -10 | - | (* Stage-079 Family-1 demand map data. *) |
| 089 | wl | 62 | `Stage-082` | 82 | -7 | - | (* Stage-082 (post-renumber) thresholds at lambda_mu = 1. Independent |
| 089 | wl-out | 11 | `Stage-080` | 80 | -9 | OUT | Stage-080 zeta_max = A_F1 pi^2/4 diff = 0``49.306707698469 |
| 090 | wl | 61 | `Stage 079` | 79 | -11 | - | (* Stage 079 transport map: zeta_req < A_F1 ==> Pe_req = 0. The inequality |

## FORWARD  (13)

| canon | kind | line | token | cited | Δ | flags | context |
| ---: | --- | ---: | --- | ---: | ---: | --- | --- |
| 012 | py | 5 | `Stage-021` | 21 | +9 | - | This script starts one step closer to the Stage-021 / Stage-022 one-port moving-throat |
| 012 | py | 5 | `Stage-022` | 22 | +10 | - | This script starts one step closer to the Stage-021 / Stage-022 one-port moving-throat |
| 012 | py | 55 | `Stage-021` | 21 | +9 | - | # Stage-021 / Stage-022 primitive one-port formulas |
| 012 | py | 55 | `Stage-022` | 22 | +10 | - | # Stage-021 / Stage-022 primitive one-port formulas |
| 012 | wl | 50 | `Stage 021` | 21 | +9 | - | Print["M1 primitive one-port forms (carried from Stage 021 / Stage 022):"]; |
| 012 | wl | 50 | `Stage 022` | 22 | +10 | - | Print["M1 primitive one-port forms (carried from Stage 021 / Stage 022):"]; |
| 012 | wl-out | 2 | `Stage 021` | 21 | +9 | OUT | M1 primitive one-port forms (carried from Stage 021 / Stage 022): |
| 012 | wl-out | 2 | `Stage 022` | 22 | +10 | OUT | M1 primitive one-port forms (carried from Stage 021 / Stage 022): |
| 017 | py | 62 | `Stage-022` | 22 | +5 | - | # Wall-only contributions to the Stage-022 weak-axisymmetric gates. |
| 082 | py | 100 | `stage 084` | 84 | +2 | - | # which equals stage 084's zeta_max^(F1) constant. |
| 082 | py | 129 | `stage 084` | 84 | +2 | - | # Reference constant from stage 084 .wl (verified at v2): zeta_max^(F1) ~ 2.4675292294558... |
| 082 | py | 130 | `Stage 084` | 84 | +2 | - | # Stage 084 lines 73-76 compute the same Pe->oo limit and assert it matches the |
| 082 | wl | 94 | `stage 084` | 84 | +2 | - | matching stage 084's zeta_max^(F1) constant. *) |

## OLD-BACK  (21)

| canon | kind | line | token | cited | Δ | flags | context |
| ---: | --- | ---: | --- | ---: | ---: | --- | --- |
| 021 | py | 9 | `Stage-3` | 3 | -18 | - | This script backs the next step after the matter-coupled Stage-3 reduction: |
| 022 | py | 204 | `stage4` | 4 | -18 | stem | def stage4_prototype_transfer_coefficients() -> dict[str, sp.Expr]: |
| 040 | py | 121 | `Stage18` | 18 | -22 | - | expect_zero("F_U(R_U=1) - Stage18 F", sp.simplify(F_U.subs(R_U, 1) - F_stage18)) |
| 040 | py | 122 | `Stage19` | 19 | -21 | - | expect_zero("G_U(R_U=1) - Stage19 G", sp.simplify(G_U.subs(R_U, 1) - G_stage19)) |
| 040 | py-out | 26 | `Stage 22` | 22 | -18 | OUT | 40.3 — Specialization to the split-U continuum of Stage 22 |
| 040 | py-out | 30 | `Stage18` | 18 | -22 | OUT | F_U(R_U=1) - Stage18 F = 0 |
| 040 | py-out | 31 | `Stage19` | 19 | -21 | OUT | G_U(R_U=1) - Stage19 G = 0 |
| 040 | wl | 128 | `Stage18` | 18 | -22 | - | expectZero["F_U(R_U=1) - Stage18 F", (fU /. rU -> 1) - fStage18]; |
| 040 | wl | 129 | `Stage19` | 19 | -21 | - | expectZero["G_U(R_U=1) - Stage19 G", (gU /. rU -> 1) - gStage19]; |
| 040 | wl-out | 35 | `Stage18` | 18 | -22 | OUT | F_U(R_U=1) - Stage18 F = 0 |
| 040 | wl-out | 36 | `Stage18` | 18 | -22 | OUT | PASS: F_U(R_U=1) - Stage18 F |
| 040 | wl-out | 37 | `Stage19` | 19 | -21 | OUT | G_U(R_U=1) - Stage19 G = 0 |
| 040 | wl-out | 38 | `Stage19` | 19 | -21 | OUT | PASS: G_U(R_U=1) - Stage19 G |
| 045 | py-out | 35 | `Stage-27` | 27 | -18 | OUT | 4. Exact collapse of the Stage-27 branch equation |
| 047 | py | 121 | `Stage-30` | 30 | -17 | - | # Negative control: if the support packet is perturbed off the exact Stage-30 |
| 078 | py | 26 | `Stage-77` | 77 | -1 | - | # Stage-77 family-1 Theta extraction: |
| 078 | py | 32 | `Stage-75` | 75 | -3 | - | # Stage-75 family-1 threshold window: |
| 078 | wl | 38 | `Stage-75` | 75 | -3 | - | their Stage-75 symbolic closed forms (no SymPy input).             *) |
| 078 | wl | 44 | `Stage-75` | 75 | -3 | - | (* Independent closed form for Theta_suff from Stage-75 sympy output line 21: |
| 078 | wl | 48 | `Stage-77` | 77 | -1 | - | (* The chi^2 and Jensen-floor Theta values are adopted from Stage-77 at |
| 078 | wl | 74 | `Stage-77` | 77 | -1 | - | closed form (thetaFailSym) and the chi/J Stage-77 numerics.   *) |

## VARIABLE  (99)

| canon | kind | line | token | cited | Δ | flags | context |
| ---: | --- | ---: | --- | ---: | ---: | --- | --- |
| 004 | py | 24 | `stage005` | 5 | +1 | stem | "moving_throat_pde_stage005_projected_maxwell_covariant_sympy_audit.py", |
| 004 | py | 25 | `stage006` | 6 | +2 | stem | "moving_throat_pde_stage006_projected_maxwell_vector_sympy_audit.py", |
| 004 | py | 26 | `stage007` | 7 | +3 | stem | "moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.py", |
| 021 | py | 3 | `stage4` | 4 | -17 | stem | moving_throat_pde_stage4_maxwell_mixed_sympy_audit.py |
| 022 | py | 270 | `stage4` | 4 | -18 | - | A_stage4 = sp.simplify(Y2_hat.coeff(omega, 2)) |
| 022 | py | 271 | `stage4` | 4 | -18 | - | B_stage4 = sp.simplify(Y2_hat.coeff(omega, 4)) |
| 022 | py | 272 | `stage4` | 4 | -18 | - | G5_stage4 = sp.simplify(Y2_hat.coeff(omega, 5) / I) |
| 022 | py | 275 | `stage4` | 4 | -18 | - | expect_zero("Stage-021 A coefficient", A_stage4 - a**2 / (9 * c_s**2)) |
| 022 | py | 276 | `stage4` | 4 | -18 | - | expect_zero("Stage-021 B coefficient", B_stage4 - 4 * a**4 / (81 * c_s**4)) |
| 022 | py | 277 | `stage4` | 4 | -18 | - | expect_zero("Stage-021 G5 coefficient", G5_stage4 - a**5 / (27 * c_s**5)) |
| 022 | py | 278 | `stage4` | 4 | -18 | - | print("A_stage4 =") |
| 022 | py | 279 | `stage4` | 4 | -18 | - | sp.pprint(A_stage4) |
| 022 | py | 280 | `stage4` | 4 | -18 | - | print("B_stage4 =") |
| 022 | py | 281 | `stage4` | 4 | -18 | - | sp.pprint(B_stage4) |
| 022 | py | 282 | `stage4` | 4 | -18 | - | print("G5_stage4 =") |
| 022 | py | 283 | `stage4` | 4 | -18 | - | sp.pprint(G5_stage4) |
| 022 | py | 285 | `stage4` | 4 | -18 | - | Gamma5_port = G5_stage4 |
| 022 | py | 303 | `stage4` | 4 | -18 | - | K2_target = sp.simplify(NQ_target * A_stage4) |
| 022 | py | 304 | `stage4` | 4 | -18 | - | K4_target = sp.simplify(NQ_target * B_stage4) |
| 022 | py-out | 155 | `stage4` | 4 | -18 | OUT | A_stage4 = |
| 022 | py-out | 161 | `stage4` | 4 | -18 | OUT | B_stage4 = |
| 022 | py-out | 167 | `stage4` | 4 | -18 | OUT | G5_stage4 = |
| 022 | wl | 135 | `Stage4` | 4 | -18 | - | Clear[z, j2, y2, h2, lambda2, lambda2Series, y2, y2Static, y2Hat, y2HatOmega, aStage4, bStage4, g5Stage4]; |
| 022 | wl | 149 | `Stage4` | 4 | -18 | - | aStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 2], Assumptions -> $Assumptions]; |
| 022 | wl | 150 | `Stage4` | 4 | -18 | - | bStage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 4], Assumptions -> $Assumptions]; |
| 022 | wl | 151 | `Stage4` | 4 | -18 | - | g5Stage4 = FullSimplify[Coefficient[Expand[y2Hat], omega, 5]/I, Assumptions -> $Assumptions]; |
| 022 | wl | 154 | `Stage4` | 4 | -18 | - | expectZero["Stage-021 A coefficient", aStage4 - a^2/(9*cS^2)]; |
| 022 | wl | 155 | `Stage4` | 4 | -18 | - | expectZero["Stage-021 B coefficient", bStage4 - 4*a^4/(81*cS^4)]; |
| 022 | wl | 156 | `Stage4` | 4 | -18 | - | expectZero["Stage-021 G5 coefficient", g5Stage4 - a^5/(27*cS^5)]; |
| 022 | wl | 158 | `Stage4` | 4 | -18 | - | gamma5Port = g5Stage4; |
| 022 | wl | 161 | `Stage4` | 4 | -18 | - | k2Target = FullSimplify[ratioTarget*aStage4, Assumptions -> $Assumptions]; |
| 022 | wl | 162 | `Stage4` | 4 | -18 | - | k4Target = FullSimplify[ratioTarget*bStage4, Assumptions -> $Assumptions]; |
| 023 | wl | 215 | `Stage4` | 4 | -19 | - | Clear[z, j2, y2, h2, lambda2, lambda2Series, y2Resp, y2Static, y2Hat, g5Stage4]; |
| 023 | wl | 228 | `Stage4` | 4 | -19 | - | g5Stage4 = FullSimplify[Coefficient[y2Hat, omega, 5]/I, Assumptions -> $Assumptions]; |
| 023 | wl | 229 | `Stage4` | 4 | -19 | - | gamma5Port = g5Stage4; |
| 025 | py | 120 | `stage023` | 23 | -2 | stem | # scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342, |
| 025 | wl | 80 | `stage023` | 23 | -2 | stem | scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342, |
| 040 | py | 113 | `stage18` | 18 | -22 | - | # F_stage18 reproduces the Stage-035 closed-form normalization F(xi, delta) |
| 040 | py | 114 | `stage035` | 35 | -5 | stem | # verified in scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py lines 46-58. |
| 040 | py | 115 | `stage19` | 19 | -21 | - | # G_stage19 reproduces the Stage-036 closed-form loading G(xi, delta) |
| 040 | py | 116 | `stage036` | 36 | -4 | stem | # verified in scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py lines 53-70. |
| 040 | py | 118 | `stage18` | 18 | -22 | - | F_stage18 = sp.simplify((9 * delta + 11 * xi)**4 / (81 * (1 - xi) * (9 * delta**2 + 18 * delta * xi + 11 * xi**2)**2)) |
| 040 | py | 119 | `stage19` | 19 | -21 | - | G_stage19 = sp.simplify(9 * xi * (delta + xi) / (9 * delta + 11 * xi)) |
| 040 | py | 121 | `stage18` | 18 | -22 | - | expect_zero("F_U(R_U=1) - Stage18 F", sp.simplify(F_U.subs(R_U, 1) - F_stage18)) |
| 040 | py | 122 | `stage19` | 19 | -21 | - | expect_zero("G_U(R_U=1) - Stage19 G", sp.simplify(G_U.subs(R_U, 1) - G_stage19)) |
| 040 | py | 132 | `stage18` | 18 | -22 | - | HF = sp.simplify(sp.diff(F_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / F_stage18) |
| 040 | py | 133 | `stage19` | 19 | -21 | - | HG = sp.simplify(sp.diff(G_U.subs(R_U, 1 + eps), eps).subs(eps, 0) / G_stage19) |
| 040 | py | 140 | `stage18` | 18 | -22 | - | HF_direct = sp.simplify(sp.diff(F_general_eps, eps).subs(eps, 0) / F_stage18) |
| 040 | py | 141 | `stage19` | 19 | -21 | - | HG_direct = sp.simplify(sp.diff(G_general_eps, eps).subs(eps, 0) / G_stage19) |
| 040 | wl | 115 | `Stage18` | 18 | -22 | - | (* fStage18 reproduces the Stage-035 closed-form F(xi, delta) verified  *) |
| 040 | wl | 116 | `stage035` | 35 | -5 | stem | (* in mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl lines 50-61. *) |
| 040 | wl | 117 | `Stage19` | 19 | -21 | - | (* gStage19 reproduces the Stage-036 closed-form G(xi, delta) verified  *) |
| 040 | wl | 118 | `stage036` | 36 | -4 | stem | (* in mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl lines 41-58. *) |
| 040 | wl | 120 | `Stage18` | 18 | -22 | - | fStage18 = FullSimplify[ |
| 040 | wl | 124 | `Stage19` | 19 | -21 | - | gStage19 = FullSimplify[9 xi (delta + xi)/(9 delta + 11 xi), Assumptions -> $Assumptions]; |
| 040 | wl | 136 | `Stage18` | 18 | -22 | - | hF = FullSimplify[(D[fU /. rU -> 1 + eps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions]; |
| 040 | wl | 137 | `Stage19` | 19 | -21 | - | hG = FullSimplify[(D[gU /. rU -> 1 + eps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions]; |
| 040 | wl | 141 | `Stage18` | 18 | -22 | - | hFDirect = FullSimplify[(D[fGeneralEps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions]; |
| 040 | wl | 142 | `Stage19` | 19 | -21 | - | hGDirect = FullSimplify[(D[gGeneralEps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions]; |
| 042 | py | 106 | `stage23` | 23 | -19 | - | F_stage23 = sp.simplify( |
| 042 | py | 112 | `stage23` | 23 | -19 | - | expect_zero("tracking collapse", F_track - F_stage23) |
| 042 | wl | 87 | `Stage23` | 23 | -19 | - | fStage23 = FullSimplify[ |
| 042 | wl | 94 | `Stage23` | 23 | -19 | - | expectZero["tracking collapse", fTrack - fStage23]; |
| 045 | py | 173 | `stage044` | 44 | -1 | stem | # See: scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146 |
| 045 | py | 174 | `stage044` | 44 | -1 | - | D_cont_stage044 = sp.simplify( |
| 045 | py | 178 | `stage044` | 44 | -1 | - | F_cont_stage044 = sp.simplify( |
| 045 | py | 181 | `stage044` | 44 | -1 | - | / ((1 - xi) * D_cont_stage044 ** 2) |
| 045 | py | 183 | `stage044` | 44 | -1 | - | F_track_stage044 = sp.simplify(F_cont_stage044.subs(R_phi, R_U)) |
| 045 | py | 189 | `stage044` | 44 | -1 | - | expect_zero("Stage-044 tracking F collapse", F_track_stage044 - F_track_expected) |
| 045 | py | 191 | `stage044` | 44 | -1 | - | F_tr_from_stage044 = sp.simplify( |
| 045 | py | 192 | `stage044` | 44 | -1 | - | F_cont_stage044.subs([(R_phi, R_U), (lam0, lam0_dn)]) |
| 045 | py | 199 | `stage044` | 44 | -1 | - | expect_zero("F_tr collapse from Stage-044 residual", F_tr_from_stage044 - F_tr_expected) |
| 045 | py | 200 | `stage044` | 44 | -1 | - | print("coherent normalization residual =", sp.simplify(R_target - F_tr_from_stage044)) |
| 045 | wl | 121 | `stage044` | 44 | -1 | stem | See: mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:81-89,128-138 *) |
| 045 | wl | 122 | `Stage044` | 44 | -1 | - | dContStage044 = FullSimplify[ |
| 045 | wl | 126 | `Stage044` | 44 | -1 | - | fContStage044 = FullSimplify[ |
| 045 | wl | 129 | `Stage044` | 44 | -1 | - | ((1 - xi) dContStage044^2), |
| 045 | wl | 132 | `Stage044` | 44 | -1 | - | fTrackStage044 = FullSimplify[fContStage044 /. rPhi -> rU, Assumptions -> $Assumptions]; |
| 045 | wl | 138 | `Stage044` | 44 | -1 | - | expectZero["Stage-044 tracking F collapse", fTrackStage044 - fTrackExpected]; |
| 045 | wl | 140 | `Stage044` | 44 | -1 | - | fTrFromStage044 = FullSimplify[ |
| 045 | wl | 141 | `Stage044` | 44 | -1 | - | fContStage044 /. {rPhi -> rU, lambda0 -> lambda0DN}, |
| 045 | wl | 149 | `Stage044` | 44 | -1 | - | expectZero["F_tr collapse from Stage-044 residual", fTrFromStage044 - fTrExpected]; |
| 045 | wl | 150 | `Stage044` | 44 | -1 | - | Print["coherent normalization residual = ", fmt[FullSimplify[rTarget - fTrFromStage044, Assumptions -> $Assumptions]]]; |
| 050 | py | 17 | `stage049` | 49 | -1 | stem | from moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit import twin_support_ratio |
| 069 | py | 25 | `stage066` | 66 | -3 | - | `scripts/moving_throat_pde_stage066_*.py` produces this form from the |
| 069 | py | 30 | `stage068` | 68 | -1 | - | `scripts/moving_throat_pde_stage068_*.py` produces this form from the |
| 070 | py | 63 | `stage48` | 48 | -22 | - | J1_stage48 = sp.simplify(If_sym / Hw) |
| 070 | py | 64 | `stage48` | 48 | -22 | - | expect_zero("I1 / J1 - 4 pi a^2 ell (sech profile)", sp.simplify(I1_constH/J1_stage48 - 4*pi*a**2*ell)) |
| 070 | wl | 78 | `Stage48` | 48 | -22 | - | kappaCmp, WwallNum, WwallCmp, JfromProfile, J1Stage48, XiNum, XiCmp, tol}, |
| 078 | py | 29 | `stage077` | 77 | -1 | stem | # Source: scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt:22-23 |
| 078 | py | 37 | `stage075` | 75 | -3 | stem | # Source: scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt:30,34 |
| 087 | py | 12 | `stage081` | 81 | -6 | - | - scripts/moving_throat_pde_stage081_*_sympy_audit.py / .wl  (a.k.a. former stage 65) |
| 087 | py | 13 | `stage082` | 82 | -5 | - | - scripts/moving_throat_pde_stage082_*_sympy_audit.py / .wl  (former stage 69 closure) |
| 087 | py | 14 | `stage085` | 85 | -2 | stem | - scripts/moving_throat_pde_stage085_quadrupole_demand_cancellation_*  (Pi_tr/C_mix cancellation) |
| 087 | py | 15 | `stage086` | 86 | -1 | stem | - scripts/moving_throat_pde_stage086_family1_loading_ratio_window_*  (Family-1 window) |
| 087 | wl | 51 | `stage085` | 85 | -2 | stem | - mathematica/moving_throat_pde_stage085_quadrupole_demand_cancellation_* |
| 087 | wl | 52 | `stage086` | 86 | -1 | stem | - mathematica/moving_throat_pde_stage086_family1_loading_ratio_window_* |
| 088 | py | 112 | `stage085` | 85 | -3 | - | # identity (verified upstream in scripts/moving_throat_pde_stage085_*). Substitute |
| 089 | py | 61 | `stage082` | 82 | -7 | - | # Source: scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt and the |

## CANON-BACK (267) — 3-digit padded backward refs; default LEAVE, spot-check only

| canon | stem | count |
| ---: | --- | ---: |
| 002 | stage002_breathing_reduction | 8 |
| 003 | stage003_bdg | 1 |
| 014 | stage014_projected_maxwell_mouth_taylor_gate_bridge | 1 |
| 022 | stage022_grouped_p2 | 16 |
| 022 | stage022_grouped_p2_normalization_bridge | 10 |
| 023 | stage023_full_grouped_bundle | 8 |
| 024 | stage024_overlap_isotropy | 6 |
| 025 | stage025_minimal_isotropic_normalization | 1 |
| 026 | stage026_concrete_axial_overlaps | 5 |
| 027 | stage027_nonconstant_axial_family | 10 |
| 028 | stage028_loaded_profile_selection | 5 |
| 029 | stage029_dynamic_loading | 4 |
| 031 | stage031_selected_branch_reachability | 1 |
| 033 | stage033_microscopic_normalization_equation | 1 |
| 036 | stage036_support_feasibility_frontier | 4 |
| 037 | stage037_continuum_kernel | 2 |
| 038 | stage038_dimensionless_continuum_placement | 8 |
| 039 | stage039_split_u_sector | 4 |
| 040 | stage040_generalized_selected_branch | 10 |
| 041 | stage041_rank2_support | 3 |
| 042 | stage042_rank2_selected_mode | 10 |
| 044 | stage044_continuum_selected_rank2 | 5 |
| 045 | stage045_coherent_local_tracking | 15 |
| 047 | stage047_coherent_kernel_map | 2 |
| 050 | stage050_zeta_threshold_comparison | 3 |
| 051 | stage051_lowest_twin_criterion | 4 |
| 055 | stage055_explicit_reachability | 3 |
| 056 | stage056_transport_source_asymmetry | 1 |
| 057 | stage057_physical_parameter_map | 4 |
| 059 | stage059_operator_branch_residual_bounds | 3 |
| 060 | stage060_entropic_microclosure | 2 |
| 061 | stage061_microscopic_gain_thresholds | 1 |
| 063 | stage063_parent_thresholds | 1 |
| 064 | stage064_equilibrium_alignment | 8 |
| 065 | stage065_thin_wall_confinement | 1 |
| 066 | stage066_wall_figure_of_merit | 2 |
| 068 | stage068_resonance_thresholds | 3 |
| 069 | stage069_final_reduced_verdict | 18 |
| 070 | stage070_gnls_wall_shell | 4 |
| 078 | stage078_family1_branch_verdict | 3 |
| 080 | stage080_family1_zeta_thresholds | 5 |
| 081 | stage081_family1_pi_thresholds | 2 |
| 082 | stage082_master_quadrupole_residual | 9 |
| 083 | stage083_family1_direct_operator_window | 4 |
| 086 | stage086_family1_loading_ratio_window | 1 |
| 087 | stage087_outgoing_branch_loading_ratio_finish | 5 |
| 088 | stage088_loading_ratio_from_minimal_module | 2 |
| 089 | stage089_family1_minimal_isotropic_verdict | 26 |
| 090 | stage090_updated_reduced_status | 12 |
