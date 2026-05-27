---
unit_id: 063
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage063_parent_thresholds.md
  paper_appendix: present
---

# Audit unit 063 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_063.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage063_parent_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 104; `\input{stages/stage_063}` at line 244)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage063_parent_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 063 turns the parent-action gain formula (Stage 062's `G_micro = rho_* g_phi^2 O_sp^2 / (m c_s*^2 K_X N_ss)`) and the Stage-44 operator phase diagram into explicit fail/succeed surfaces. The card's `\stagefield{Output}` reads verbatim: "Amplitude thresholds \eqref{eq:app-stage063-gphi-thresholds}, coherence thresholds \eqref{eq:app-stage063-coherence-thresholds}, and Cauchy no-go \eqref{eq:app-stage063-cauchy-nogo}." The three boxed deliverables are:

1. Amplitude thresholds: `g_(phi,fail)^2 = m c_s*^2 K_X N_ss G_fail / (rho_* O_sp^2)` and the `G_suff` sibling.
2. Coherence thresholds: `C_fail^2 = m c_s*^2 K_X G_fail / (rho_* g_phi^2 N_pp)` and the `G_suff` sibling.
3. Cauchy no-go: `C_fail^2 > 1 ⟹ no source profile in the support channel can reach the fail threshold`.

The notes additionally develop the equivalence between overlap-form and coherence-form thresholds via `O_sp^2 = N_ss N_pp C_sp^2`, the best-case gain `G_max(g_phi) = rho_* g_phi^2 N_pp / (m c_s*^2 K_X)` at perfect alignment, and the kappa-insertion `kappa = K_X L^2 / T_X` together with `G_fail = Pe_req/(kappa Delta_inf)`, `G_suff = Pe_req/(kappa Delta_0)` (notes items 4–6, exhibiting explicit `K_X` cancellation). The matched-profile special case (notes section 5) is illustrative but not part of `\stagefield{Output}`.

## What the script claims to verify

The SymPy script and its Mathematica counterpart assert: (1) the overlap-form `g_phi^2` thresholds equal the paper's boxed forms (anchored via `sp.solve(G_micro - G_fail, gphi_sq)` and Mathematica `Reduce[..., gphiSq > 0, ...]`); (2) substituting `O_sp^2 = C_sp^2 N_ss N_pp` reduces the overlap-form thresholds to the coherence-form thresholds (with `N_ss` cancelling); (3) the coherence thresholds `C_fail^2, C_suff^2` are written down explicitly and the ratio `C_suff^2/C_fail^2 = G_suff/G_fail` is verified; (4) `G_micro` factors as `G_max * C_sp^2` and saturates `G_max` at `O_sp^2 = N_ss N_pp` — the Cauchy-Schwarz factorization underwriting the no-go theorem; (5) after inserting `kappa = K_X L^2 / T_X` together with Stage-44's `G_fail/G_suff = Pe_req/(kappa Delta_(inf/0))`, the explicit `K_X` cancels in both overlap- and coherence-form thresholds, leaving `T_X`-controlled prefactors.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Amplitude thresholds (eq:app-stage063-gphi-thresholds) | sympy L45-48 declare `g_fail_sq, g_suff_sq`; L102-115 anchor via `sp.solve(G_micro - G_fail/suff, gphi_sq)` against hand-rearranged form. Mathematica L36-40, L87-99 mirror via `Reduce[..., gphiSq > 0, ...]`. | match |
| Coherence thresholds (eq:app-stage063-coherence-thresholds) | sympy L51-67 verify substitution `O_sp^2 -> C2 N_ss N_pp` collapses overlap form to coherence form; checks `C_suff^2 / C_fail^2 = G_suff/G_fail`. Mathematica L41-55 mirror. | match |
| Cauchy no-go (eq:app-stage063-cauchy-nogo) | sympy L70-74, L119-122 verify factorization `G_micro = G_max * C_sp^2` and saturation `G_micro(O_sp^2 = N_ss N_pp) = G_max`. Combined with the explicit `C_fail^2 = G_fail/G_max` (from L63 vs L70), the no-go `C_fail^2 > 1 ⇔ G_max < G_fail` is an immediate algebraic consequence of the verified identities. Mathematica L57-58, L102-105 mirror. | match (algebraic underwriting of the inequality; the no-go is a one-step corollary of the verified saturation + definition) |
| Notes item 6 (kappa-substitution, K_X cancellation) | sympy L77-100; Mathematica L60-83. Four `expect_zero` calls cover overlap- and coherence-form, fail- and suff-branches. | match (notes-side; supports paper's amplitude thresholds in their Pe_req form) |

The `paper_alignment` field is `aligned`: each of the three boxed paper deliverables is covered by load-bearing assertions in both engines, and notes-side material is additionally exercised.

## Assertion inventory

| #   | Script      | Line     | Form                                                                                  | Exercises which paper claim?                       | Anchored to claim? |
|-----|-------------|----------|---------------------------------------------------------------------------------------|----------------------------------------------------|--------------------|
| A1  | sympy       | 51-55    | `simplify(g_fail_sq[O_sp^2→C^2 N_ss N_pp] - m cs^2 K_X G_fail/(rho_* N_pp C^2)) == 0` | claim 1 ↔ claim 2 equivalence                       | yes                |
| A2  | sympy       | 56-60    | same for `g_suff_sq`                                                                  | claim 1 ↔ claim 2 (suff)                            | yes                |
| A3  | sympy       | 67       | `simplify(C_suff_sq / C_fail_sq - G_suff/G_fail) == 0`                                | claim 2 ratio internal consistency                  | yes                |
| A4  | sympy       | 71-74    | `G_micro[O_sp^2 → C^2 N_ss N_pp] - G_max * C^2 == 0`                                  | factorization underwriting claim 3                  | yes                |
| A5  | sympy       | 80-84    | `g_fail_sq[G_fail → Pe_req/(κ Δ_inf), κ → K_X L^2/T_X]` matches T_X form              | notes item 6, K_X cancellation                      | yes                |
| A6  | sympy       | 85-89    | same for `g_suff_sq` with `Delta_0`                                                   | notes item 6                                        | yes                |
| A7  | sympy       | 91-95    | coherence-form fail threshold with kappa inserted                                     | notes item 6 (coherence side)                       | yes                |
| A8  | sympy       | 96-100   | coherence-form suff threshold with kappa inserted                                     | notes item 6                                        | yes                |
| A9  | sympy       | 105-115  | `sp.solve(G_micro - G_(fail/suff), gphi_sq)` matches `g_(fail/suff)_sq`               | independent derivation of claim 1 from `G_micro`     | yes                |
| A10 | sympy       | 119-122  | `G_micro[O_sp^2 → N_ss N_pp] - G_max == 0`                                            | Cauchy saturation; supports claim 3                  | yes                |
| B1  | mathematica | 42-45    | `expectZero[g_fail_sq - m cs^2 K_X G_fail/(rho_* N_pp c2Def)]` with `c2Def=oSP^2/(nSS nPP)` | claim 1 ↔ claim 2                              | yes                |
| B2  | mathematica | 46-49    | same for `g_suff_sq`                                                                  | claim 1 ↔ claim 2                                    | yes                |
| B3  | mathematica | 55       | `cSuffSq/cFailSq - gSuff/gFail`                                                       | claim 2                                              | yes                |
| B4  | mathematica | 58       | `(G_micro /. oSP^2 → c2 nSS nPP) - gMax c2`                                          | claim 3 factorization                                | yes                |
| B5  | mathematica | 64-68    | overlap-form fail with kappa inserted                                                 | notes item 6                                         | yes                |
| B6  | mathematica | 69-73    | overlap-form suff with kappa inserted                                                 | notes item 6                                         | yes                |
| B7  | mathematica | 74-78    | coherence-form fail with kappa inserted                                               | notes item 6                                         | yes                |
| B8  | mathematica | 79-83    | coherence-form suff with kappa inserted                                               | notes item 6                                         | yes                |
| B9  | mathematica | 88-99    | `Reduce[gMicro == gFail/suff && gphiSq > 0]` matches hand form (positive-root branch) | independent derivation of claim 1 via Reduce         | yes                |
| B10 | mathematica | 102-105  | `(G_micro /. oSP^2 → nSS nPP) - gMax`                                                | Cauchy saturation; supports claim 3                  | yes                |

All rows are non-tautological: each `expect_zero`/`expectZero` requires actual algebraic cancellation (e.g., `N_ss` cancellation in coherence substitution; positive-root selection in Solve/Reduce; factorization of `G_micro` through `O_sp^2 → C^2 N_ss N_pp`; explicit `K_X` cancellation after kappa insertion). The Solve-derivation pair (A9) and Reduce-derivation pair (B9) are the load-bearing anchors that tie the hand-written threshold expressions back to `G_micro`; if `G_micro` were rewritten with the wrong factor of `N_ss` or `N_pp`, A9/B9 would fail.

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script is not a transliteration of the SymPy script. Specifically:

1. **Solver path differs**: SymPy uses `sp.solve(G_micro_gphi - G_fail, gphi_sq)` (algebraic Solve returning the unique root); Mathematica uses `Reduce[(gMicro /. gPhi^2 -> gphiSq) == gFail && gphiSq > 0, gphiSq, Reals]` followed by `Cases[..., HoldPattern[gphiSq == rhs_] :> rhs, Infinity][[1]]` to extract the positive-root branch from the Reals-domain reduction. These are genuinely different algebraic-derivation pipelines — Mathematica's Reduce solves under a domain restriction, SymPy's solve does not.

   - sympy L105: `sol_fail = sp.solve(G_micro_gphi - G_fail, gphi_sq)`
   - mathematica L88: `reduceFail = Reduce[(gMicro /. gPhi^2 -> gphiSq) == gFail && gphiSq > 0, gphiSq, Reals];`

2. **Coherence-substitution direction differs**: Mathematica defines `c2Def = oSP^2/(nSS*nPP)` (L41) and substitutes via FullSimplify against this ratio; SymPy uses `.subs(Osp**2, C2 * Nss * Npp)`. The substitution direction is opposite (Mathematica defines `c2` as the ratio; SymPy substitutes `O_sp^2 = C^2 N_ss N_pp`). Both probe the same identity from opposite sides.

   - sympy L53: `g_fail_sq.subs(Osp**2, C2 * Nss * Npp) - sp.simplify(m * cs_star_sq * KX * G_fail / (rho_star * Npp * C2))`
   - mathematica L44: `gFailSq - FullSimplify[m*csStarSq*kX*gFail/(rhoStar*nPP*c2Def), Assumptions -> $Assumptions]`

3. **Simplification stacks differ**: SymPy uses `sp.simplify(sp.expand(expr))` in `expect_zero`; Mathematica uses `FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions]` with a global positive-real `$Assumptions` block (L30-33). These are independent simplification pipelines.

The structural identities being verified are necessarily the same (they are paper-side identities), but the algebraic machinery used to verify them is independent.

## Engine cross-check

Both engines produce residual `0` for every check:

| Identity                                            | SymPy residual | Mathematica residual |
|-----------------------------------------------------|----------------|----------------------|
| coherence substitution in g_fail^2                  | 0              | 0                    |
| coherence substitution in g_suff^2                  | 0              | 0                    |
| C_suff^2/C_fail^2 - G_suff/G_fail                   | 0              | 0                    |
| G_micro - G_max * C^2                               | 0              | 0                    |
| K_X g_fail threshold with kappa inserted            | 0              | 0                    |
| K_X g_suff threshold with kappa inserted            | 0              | 0                    |
| coherence-form fail threshold with kappa inserted   | 0              | 0                    |
| coherence-form suff threshold with kappa inserted   | 0              | 0                    |
| g_(fail/suff)^2 from Solve/Reduce matches hand form | 0              | 0                    |
| G_max = G_micro at Cauchy saturation                | 0              | 0                    |

Printed values for `g_(phi,fail/suff)^2` and `C_(fail/suff)^2` agree symbolically up to ordering: sympy gives `G_fail*K_X*N_ss*cs_star_sq*m/(O_sp**2*rho_star)`, Mathematica gives `(csStarSq*gFail*kX*m*nSS)/(oSP^2*rhoStar)` — identical up to commutativity.

Output freshness: sympy script mtime 2026-05-22 19:41, output mtime 19:43 (fresh); Mathematica script mtime 19:42, output mtime 19:44 (fresh). No `stale_output`.

## Verdict justification

The audit is clean. Attempted attacks:

- **Tautology probe**: Each `expect_zero` candidate was traced to a substitution or solve step that requires non-trivial cancellation (e.g., `N_ss` cancellation in coherence substitution; positive-root selection in Reduce; explicit `K_X` cancellation in the kappa-inserted forms). The Solve/Reduce derivation (A9/B9) is the load-bearing anchor that connects the hand-rearranged thresholds back to `G_micro` — if `G_micro` were defined with a wrong constituent (e.g., `N_pp` instead of `N_ss`), the Solve output would diverge from `g_fail_sq` and the check would fail.
- **Hardcoded-result probe**: No literal numeric constants appear; all expressions are symbolic in the paper's named constants.
- **Symbol-assumption probe**: All symbols declared positive real, consistent with their physical roles (densities, squared sound speeds, squared overlap norms, dimensionless gains). The Cauchy bound `O_sp^2 <= N_ss N_pp` is not contradicted by `O_sp > 0`; saturation is exercised explicitly at the upper end (`O_sp^2 = N_ss N_pp`).
- **Missing-branch probe**: All three boxed paper deliverables have script-side coverage. The Cauchy no-go inequality is covered via its algebraic underpinnings — `G_micro = G_max * C_sp^2` (A4/B4) plus saturation `G_micro|_{O_sp^2 = N_ss N_pp} = G_max` (A10/B10) plus the explicit definitions `C_fail^2 = G_fail/G_max` (from L63/L70 of sympy) jointly entail `C_fail^2 > 1 ⇔ G_max < G_fail`. This is the standard symbolic-CAS form of an inequality theorem.
- **Mathematica transliteration probe**: SymPy uses `sp.solve`; Mathematica uses `Reduce` with positive-root selection. SymPy substitutes `O_sp^2 → C^2 N_ss N_pp`; Mathematica defines `c2Def = O_sp^2/(N_ss N_pp)` and substitutes in the opposite direction. These are independent algebraic paths confirming the same identity.
- **Paper round-trip**: The script's bottom-line forms match the paper's boxed equations (14), (24), and the algebraic content of (34) verbatim up to symbol renaming (e.g., `O_sp` ↔ `O_{\sigma\phi}`, `cs_star_sq` ↔ `c_{s,*}^2`). The kappa-insertion section matches notes item 6 exactly.

Cosmetic observation (not a finding): the script banners say "STAGE 46" / "STAGE 046", reflecting legacy Stage-numbering from the notes (notes title also mentions "Stage 45/46"). The file names and paper card use Stage 063; load-bearing content (assertion identities) matches the paper's Stage 063 deliverables. This is purely informational lineage from an earlier renumbering.

## Self-test notes

- **Variable independence**: no `sp.diff` or `D[…]` calls anywhere; no derivative-trap applies. The Solve/Reduce calls treat all symbols other than `gphi_sq` as parameters; no constant-derivative pitfall.
- **Symmetry/parity**: no integrals over unbounded domains; not applicable.
- **Trivial-case pre-check**: with `m = K_X = rho_* = cs_star_sq = N_ss = N_pp = O_sp = g_phi = 1` and `G_fail = 1`, `g_fail_sq` evaluates to 1; `sp.solve(G_micro - G_fail, gphi_sq)` returns `[1]`; residual 0. Substituting `O_sp^2 = N_ss N_pp` into `G_micro` gives `G_max` exactly (mechanical cancellation of `N_ss`); residual 0. Substituting `O_sp^2 = C^2 N_ss N_pp` into `g_fail_sq = m cs^2 K_X N_ss G_fail/(rho_* O_sp^2)` gives `m cs^2 K_X G_fail/(rho_* N_pp C^2)` after `N_ss` cancellation; matches the RHS on L54.
- **Paper round-trip**: no fix prescribed; nothing to round-trip.
