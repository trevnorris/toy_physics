---
unit_id: 065
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage065_thin_wall_confinement.md]
  paper_appendix: present
---

# Audit unit 065 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_065.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage065_thin_wall_confinement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 108; `\input{stages/stage_065}` at line 248)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage065_thin_wall_confinement_mathematica_audit.txt`

## What the paper claims

The stage card's `\stagefield{Output}` reads verbatim: "A physical wall-amplitude interpretation of \(g_\phi\)." The card's body equation states `g_phi = V0/ell` (eq:app-stage065-gphi) under the input confinement `V_conf(r;a) = V0 f((r-a)/ell)`, and the prose declares that "inserting this into the matched-layer gain yields thresholds directly in terms of V0, ell, support stiffness, and wall geometry." The appendix row (line 108 of `stage_appendix_part03.tex`) summarises the deliverable as "Wall amplitude thresholds for V_conf = V0 f((r-a)/ell)." The source notes enumerate eight discrete results: (1) g_phi = V0/ell; (2) I_1 = 4 pi ell [a^2 J_1 + 2 a ell J_2 + ell^2 J_3]; (3) J_2 = 0 for a centred symmetric wall layer; (4) the exact G_eq = 4 pi V0^2 [a^2 J_1/ell + 2 a J_2 + ell J_3]/K_X form; (5) the thin-wall leading gain G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell); (6) V0_fail^2 / V0_suff^2 in terms of K_X; (7) K_X cancellation after inserting kappa = K_X L^2 / T_X; (8) the constant-H reduction V0_fail^2 = H_w T_X ell Pe_req / (4 pi a^2 L^2 I_f Delta_inf) with J_1 = I_f / H_w. The script must therefore exercise these eight identities, and the paper card's bottom line is anchored on (1) plus the downstream chain (5)-(8).

## What the script claims to verify

The SymPy module docstring enumerates six bottom-line claims that map to notes-items (1)-(6) and (8). The script's actual machinery does two things: (i) a symbolic block (lines 51-84, 135-176) that constructs g_phi, I_1, G_eq, G_eq_tw and the V0_fail^2 / V0_suff^2 thresholds and verifies their algebraic decomposition plus the K_X cancellation under kappa = K_X L^2/T_X and the constant-H closed form; (ii) a concrete-profile block (lines 86-133, 178-183) that picks f(u) = exp(-u^2) and constant h' = 1, computes J_1, J_2, J_3, I_f as definite integrals via `sp.integrate`, verifies J_2 = 0 by parity, the (1, 2, 1) polynomial coefficient pattern by direct shell integration, the 1/ell scaling of g_phi via an independent r-derivative of V0 exp(-((r-a)/ell)^2), the relative-correction identity (ell/a)^2 J_3/J_1, and J_1 = I_f / H_w. The Mathematica script (lines 26-56, then 58-128) opens with its own independent shell-integral derivation block using `Integrate[...]` over a Gaussian profile to anchor J_2 = 0, the (1, 2, 1) polynomial expansion, and J_1 = I_f / H_w, and then runs a symbolic block mirroring the SymPy threshold derivation via a different algebraic route (`gFail/gEqCoeff` rather than `sp.solve`). Both scripts exit 0.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| g_phi = V0/ell (paper eq:app-stage065-gphi; notes-item 1) | sympy line 61 + independent r-derivative anchor at lines 109-123 verifying -2 V0 exp(-1)/ell at r=a+ell | match |
| I_1 = 4 pi ell [a^2 J_1 + 2 a ell J_2 + ell^2 J_3] (notes-item 2) | sympy lines 66-67 (symbolic) + lines 125-133 (direct Gaussian shell integral vs polynomial); math lines 46-50 (independent `Integrate`) | match |
| J_2 = 0 symmetric (notes-item 3) | sympy line 100 `expect_zero(..., J2_num)` with J2_num = ∫ xi (exp(-xi^2))'^2 dxi; math line 43 same | match |
| Exact G_eq form (notes-item 4) | sympy line 74 builds G_eq from g_phi^2*I_1/K_X; thin-wall remainder check at lines 81-84 anchors the algebraic decomposition | match |
| G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell) (notes-item 5) | sympy line 76 + thin-wall remainder check + concrete-profile relative-correction check (lines 102-106); math lines 76, 81-84 | match |
| V0_fail^2 / V0_suff^2 in K_X form (notes-item 6) | sympy lines 137-144 via sp.solve; math lines 88-96 via gFail/gEqCoeff (independent path) | match |
| K_X cancellation under kappa = K_X L^2/T_X (notes-item 7) | sympy lines 147-160; math lines 98-112 | match |
| Constant-H thresholds (notes-item 8) | sympy lines 164-183; math lines 116-127, plus independent J_1 = I_f/H_w check via direct integration in both engines | match |

`paper_alignment: aligned` — every deliverable in the paper card body equation and the notes has a corresponding script-side check, on both engines, with concrete-profile anchors backing the symbolic ones.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 81-84 | `expect_zero("thin-wall remainder ...", G_eq_sym - G_eq_tw - 4 pi V0^2 ell J3/KX)` | notes-items 4+5 | partial (algebraic decomposition; backed by A3 concrete check) |
| A2 | sympy | 100 | `expect_zero("J2 vanishes by parity", J2_num)` for J2_num computed via `sp.integrate(xi * (exp(-xi^2))'^2, ...)` | notes-item 3 | yes |
| A3 | sympy | 102-106 | `expect_zero("relative correction (ell/a)^2 J3/J1", rel)` with numeric J1_num, J3_num | notes-items 2+5 | yes |
| A4 | sympy | 122-123 | `expect_zero("g_phi 1/ell scaling", dV/dr at r=a+ell - (-2 V0 exp(-1)/ell))` | notes-item 1 | yes |
| A5 | sympy | 128-133 | `expect_zero("I1 polynomial coefficients (1,2,1)", I1_full - I1_poly)` via direct shell integral | notes-item 2 | yes |
| A6 | sympy | 153-156 | `expect_zero("K_X cancellation in V0_fail^2", V0_fail_sq_k - TX ell Pe_req/(4 pi a^2 L^2 J1 Delta_inf))` | notes-items 6+7 | yes |
| A7 | sympy | 158-160 | analogous K_X cancellation in V0_suff^2 | notes-items 6+7 | yes |
| A8 | sympy | 169-172 | `expect_zero("constant-H fail threshold", V0_fail_const - Hw TX ell Pe_req/(4 pi a^2 L^2 If Delta_inf))` | notes-item 8 | yes (symbolic; subs would propagate consistently but the closed-form anchor catches misplaced H_w) |
| A9 | sympy | 173-176 | constant-H suff analog | notes-item 8 | yes |
| A10 | sympy | 182-183 | `expect_zero("J1 = I_f / H_w under constant compressibility", J1_num - If_num/Hw_num)` | notes-item 8 | partial (H_w=1 numerically, so both sides reduce to ∫ f'^2; verifies dimensional placement but weak in isolation) |
| M1 | math | 43 | `expectZero["independent: J2 vanishes for symmetric profile", j2Num]` via `Integrate` | notes-item 3 | yes |
| M2 | math | 49-50 | `expectZero["I1 polynomial expansion matches direct integral", i1Direct - i1Poly]` | notes-item 2 | yes |
| M3 | math | 55-56 | `expectZero["J1 = I_f/H_w under constant compressibility", j1Num - ifMomDirect/hwSym]` | notes-item 8 | partial (same H_w=1 caveat) |
| M4 | math | 81-84 | `expectZero["thin-wall remainder ..."]` | notes-items 4+5 | partial (backed by independent block) |
| M5 | math | 105-108 | `expectZero["K_X cancellation in V0_fail^2"]` | notes-items 6+7 | yes |
| M6 | math | 109-112 | `expectZero["K_X cancellation in V0_suff^2"]` | notes-items 6+7 | yes |
| M7 | math | 121-124 | `expectZero["constant-H fail threshold"]` | notes-item 8 | yes |
| M8 | math | 125-128 | `expectZero["constant-H suff threshold"]` | notes-item 8 | yes |

The two `partial (H_w=1)` rows (A10, M3) collapse to ∫ f'^2 = ∫ f'^2 because both engines fix the constant compressibility numerically to 1. Considered in isolation, that check is weak; but the symbolic closed-form anchors A8/A9/M7/M8 carry H_w as a free symbol and catch any misplaced H_w in the V0_fail/V0_suff formulas. The A1/M4 thin-wall remainder is an algebraic decomposition (G_eq_sym = G_eq_tw + 4 pi V0^2 ell J3/K_X is implicit in the construction), but the concrete-profile relative-correction check (A3) and the direct shell integral (A5/M2) anchor the same scaling from a non-tautological route. Net effect: no standalone gap.

## Findings

(None.)

## Independent-derivation check (Mathematica)

The Mathematica script opens (lines 26-56) with an "INDEPENDENT SHELL-INTEGRAL DERIVATION (concrete Gaussian profile)" block before the symbolic block. It computes the Gaussian moments directly via Mathematica's `Integrate`:
```
j1Num = Integrate[(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j2Num = Integrate[xi*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j3Num = Integrate[xi^2*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
```
which output J1_num = Sqrt[Pi/2], J2_num = 0, J3_num = (3 Sqrt[Pi/2])/4 — matching SymPy's `sqrt(2)*sqrt(pi)/2`, `0`, `3*sqrt(2)*sqrt(pi)/8` up to algebraic rewriting. The (1, 2, 1) polynomial expansion check at lines 46-50 integrates `(fpProf[xi])^2/hConst * (aSym + ellSym*xi)^2` and compares against `aSym^2*j1Num + 2*aSym*ellSym*j2Num + ellSym^2*j3Num` — this is a Mathematica-native `Integrate` call, not echoing SymPy's `sp.integrate`. The threshold derivation in the symbolic block uses `gFail/gEqCoeff` (line 92) to extract V0^2 from G_eq^(tw) = G_fail rather than SymPy's `sp.solve(sp.Eq(G_eq_tw, G_fail), V0**2)` — different algorithmic route to the same identity. Sections 26-56 versus SymPy lines 86-133 differ: SymPy ordered the Gaussian-block AFTER the symbolic shell decomposition, Mathematica orders it BEFORE; SymPy uses `sp.integrate` with explicit infinity bounds, Mathematica uses `Integrate[..., {xi, -Infinity, Infinity}]` and `FullSimplify` with positivity assumptions. Not a transliteration. Both engines independently arrive at the same closed-form residuals (all zero).

## Engine cross-check

Both engines verify identical claims and agree on every assertion:

- J1: sympy `sqrt(2)*sqrt(pi)/2` vs. math `Sqrt[Pi/2]` (algebraically identical).
- J2: 0 (both).
- J3: sympy `3*sqrt(2)*sqrt(pi)/8` vs. math `(3 Sqrt[Pi/2])/4` (algebraically identical).
- V0_fail^2 with kappa inserted: sympy `Pe_req*T_X*ell/(4*pi*Delta_inf*J1*L**2*a**2)` vs. math `(ell*peReq*tx)/(4*a^2*deltaInf*j1*len^2*Pi)` (same expression up to symbol renaming).
- Constant-H thresholds: sympy `H_w*Pe_req*T_X*ell/(4*pi*Delta_inf*I_f*L**2*a**2)` vs. math `(ell*hw*peReq*tx)/(4*a^2*deltaInf*ifMom*len^2*Pi)` (identical).
- All `expect_zero`/`expectZero` residuals print `0`.
- Both transcripts terminate with success banners and exit status 0.

No engine disagreement.

## Verdict justification

The scripts faithfully exercise every deliverable enumerated in the source notes and the paper-card body equation. The non-tautological backbone consists of: independent symbolic derivative of V0 f((r-a)/ell) at r=a+ell to exhibit the 1/ell scaling (A4); direct shell integration for the (1, 2, 1) polynomial pattern (A5, M2); parity vanishing of J_2 (A2, M1); algebraic K_X cancellation under kappa = K_X L^2 / T_X (A6, A7, M5, M6); and the closed-form constant-H thresholds with H_w as a free symbol (A8, A9, M7, M8). Attacks attempted: (i) checking whether g_phi = V0/ell is just assigned at line 61 and never verified — refuted by the independent r-derivative anchor at lines 109-123 verifying the 1/ell scaling from a concrete profile; (ii) checking whether the thin-wall remainder (A1/M4) is algebraically guaranteed by the construction of G_eq_sym from I_1_sym — yes algebraically, but the concrete-profile relative-correction check (A3) verifies the same scaling from a different route, and A5/M2 verify the underlying polynomial coefficients independently; (iii) checking whether the K_X cancellation is built into `sp.solve` — refuted because V0_fail_sq still carries K_X in its kappa form (K_X*Pe_req*ell/(4*pi*Delta_inf*J1*a^2*kappa)) and only the explicit `subs(kappa -> K_X L^2 / T_X)` cancels it (the residual subtraction would fail if the K_X exponent were wrong); (iv) checking whether the constant-H J_1 = I_f/H_w check at A10/M3 is meaningful — answered: weak in isolation (H_w=1 numerically), but the symbolic A8/A9/M7/M8 carry the burden. Outputs are fresher than scripts (mtimes 1779501135/1779501141 > 1779501055 for both pairs). The docstring at the top of the .py and the opening banner label this audit as "Stage 48" (and the .wl interior banner says "STAGE 048"), reflecting an older stage-number scheme; the file paths, paper card filename, appendix row at line 108, and final Mathematica banner ("Stage 065 Mathematica audit passed.") all correctly identify this as Stage 065 — purely cosmetic, no claim drift. Verdict: clean.

## Self-test notes

Checked: (1) variable independence — `f_loc = exp(-((r_sym - a)/ell)**2)` at line 114 genuinely depends on r_sym, so `sp.diff(V0 * f_loc, r_sym)` at line 115 yields V0*(-2(r-a)/ell^2)*exp(...) (not identically zero); at r=a+ell this evaluates to -2 V0 exp(-1)/ell, matching the asserted residual. (2) Symmetry/parity — for the J_2 check, integrand xi * (exp(-xi^2))'^2 = xi * 4 xi^2 exp(-2 xi^2) = 4 xi^3 exp(-2 xi^2) is odd in xi over a symmetric domain → integral is 0; SymPy and Mathematica both return 0. (3) Trivial-case pre-check — for the relative-correction check at A3, substituting J1_num = sqrt(2 pi)/2 and J3_num = 3 sqrt(2 pi)/8 gives ((4 pi V0^2 ell J3_num)/KX) / (4 pi a^2 V0^2 J1_num/(KX ell)) = ell^2 J3_num/(a^2 J1_num); subtracting (ell^2 J3_num)/(a^2 J1_num) gives 0 as asserted. (4) Path specifications — both engines present at the expected paths (`scripts/`, `mathematica/`); no missing-script finding. (5) Paper round-trip — not applicable, no fix prescribed; paper-card body equation `g_phi = V0/ell` aligns with sympy line 61 and math line 66; all eight notes deliverables are anchored.
