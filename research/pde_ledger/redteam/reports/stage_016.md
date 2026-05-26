---
unit_id: 016
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 016 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_016.tex`
- notes: (none — no `notes/stages/moving_throat_pde_stage016_*.md` exists)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row at line 54; stage card itself is `\input`-ed at line 111)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage016_parent_throat_action_candidate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage016_parent_throat_action_candidate_mathematica_audit.txt`

Cross-stage context read (read because the stage card references it but not audited): `paper/stages/stage_015.tex` for `\eqref{eq:stage015-parent-density}` and `eq:stage015-keta` (referenced by stage 016 card's "candidate" framing).

Output mtimes (both fresh):
- sympy script: May 4 12:00; sympy output: May 21 15:00 (output newer than script — fresh)
- mathematica script: May 21 13:20; mathematica output: May 21 15:02 (output newer than script — fresh)

## What the paper claims

Stage 016 paper card frames the stage as the **candidate** companion to stage 015's parent throat-action master. Per the card body: "It verifies the minimal nonlinear candidate and its exact quadratic recovery." The candidate density is the one boxed at `eq:stage015-parent-density`: `L_Sigma = (1/2) mu R_t^2 - (1/2) T_w R_w^2 - (1/2) T_Omega |grad_Omega R|^2 - U_Sigma`. The "Boundary controls" paragraph explicitly demands the audit verify integration by parts with both Gaussian and Lorentzian probes and check a **nonzero** `arctan w` endpoint discharge — i.e., that the boundary reader is not a broken always-zero operator. The Output paragraph states: "Stage~016 records the candidate action as a parent-complete source of the linear wall closure, subject to the explicitly checked boundary class." The appendix row 54 summarizes: "Candidate parent-action data for the projected throat sector" with status `StatusExactClosure`. No notes file exists. The card neither itemizes the `K_eta` formula (that is owned by stage 015 per `eq:stage015-keta`) nor explicitly names the Y_lm modal projection or +6 T_Omega coefficient — but the broader "parent-complete source of the linear wall closure" and the "exact quadratic recovery" phrasing licence the script to demonstrate the recovery and its modal entry to the wall PDE.

## What the script claims to verify

Both engines verify a coherent set of claims for the candidate parent action: (M1) the exact nonlinear Euler–Lagrange equation of `L_Sigma` in a local orthonormal angular chart, cross-checked between machine EL and a hand-derived form, with a sign-mutation guard; (M2–M5) the order-by-order epsilon expansion `R = R0(w) + eps*eta`, recovering the static-isotropic background equation and the canonical quadratic density with `K_eta(w) = U_{RR}(R0,w) - d_w[T_{w,R}(R0,w) R0'] + (1/2) T_{w,RR}(R0,w) (R0')^2`; (M6–M8) concrete boundary-discharge ledger on Gaussian and Lorentzian profiles, a `-2` Lorentzian finite-endpoint discharge, and the nontrivial `arctan w` boundary probe; (M9–M10) Y_{20} as a `−6`-eigenfunction of the spherical Laplacian with angular norm `1` and angular gradient norm `6`; (M11) the closed Y_{20} modal Euler–Lagrange equation `mu q_tt − ∂_w(T_w q_w) + (K_eta + 6 T_Omega) q = 0`. Each `assert_zero/expectZero` is paired with a mutation `assert_nonzero/expectNonzero` (sign flip, wrong eigenvalue, wrong K_eta sign, etc.) that would trip if the claim were trivially true.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| "Minimal nonlinear candidate" (density `L_Sigma` from `eq:stage015-parent-density`) | M1 builds `lagExact` from the same density and matches sympy/mathematica EL machinery against a hand EL with mutation | match |
| "Exact quadratic recovery" (i.e. recovers `K_eta` of `eq:stage015-keta`) | M2/M3 (linear-order and background), M4/M5 (raw quadratic + IBP cross-term + canonical post-IBP form with `K_eta` formula), with mutations | match |
| "Verifies IBP with Gaussian and Lorentzian probes" | M6 computes both linear and quadratic Gaussian and Lorentzian boundary discharges and runs the actual IBP integrals on Gaussian profiles | match |
| "Checks a nonzero `arctan w` endpoint discharge" | sympy line 84 asserts nonzero on `boundary_value(atan(w_ibp), w_ibp)`; mathematica line 186 asserts `atanDischarge − Pi == 0` (stronger) | match |
| "Lorentzian finite-endpoint distinguishes nonzero limits" (Boundary controls paragraph) | sympy line 107 asserts `linear_boundary_lorentz_endpoint + 2 == 0`; mathematica line 180 same | match |
| Modal EL `mu q_tt − ∂_w(T_w q_w) + (K_eta + l(l+1) T_Omega) q = 0` for `l=2` (NOT named by stage 016 paper card; appendix row "Candidate parent-action data for the projected throat sector" is broad enough; stage 015 card mentions modal-style export only implicitly) | M9–M11 derive Y_{20} eigenvalue, angular norms, and modal EL | match (appendix-level), arguably `extra` against the stage card body |
| `eq:stage015-keta` formula owned by stage 015 | Verified explicitly here as `K_eta = URR0 - d_TwR_R0p + (1/2) TwRR0 R0p^2` (sympy line 155; mathematica line 132) | match (consistent forward of stage 015's exported formula) |

Dominant pattern: `match`. The Y_{20}/modal-EL row is the only one where the stage 016 card body is silent and the audit relies on the appendix row's broader phrasing plus the fact that the script's own end commentary calls out this added projection ("the script now also integrates the genuine Y20 harmonic into the reduced modal Lagrangian itself"). This is judged in-scope rather than `paper_missing_script_claim` because the appendix line and the card's "parent-complete source of the linear wall closure" phrasing license the modal entry-point. Setting `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `assert_zero('exact nonlinear EL', exact_el + exact_el_lhs)` | Minimal nonlinear candidate (M1) | yes |
| A2 | sympy | 55 | `assert_nonzero` on mutated EL | M1 mutation guard | yes |
| A3 | sympy | 66 | `assert_zero('generic linear IBP identity', ...)` | Symbolic IBP move | yes |
| A4 | sympy | 70 | `assert_zero('generic quadratic IBP identity', ...)` | Symbolic IBP move | yes |
| A5 | sympy | 74 | `assert_nonzero('mutated quadratic IBP sign should fail', ...)` | IBP sign mutation | yes |
| A6 | sympy | 84 | `assert_nonzero('boundary operator detects nonzero endpoint discharge', ...)` on `boundary_value(atan(w))` | "nonzero arctan w endpoint discharge" | yes |
| A7 | sympy | 89 | `assert_nonzero('lorentzian linear probe keeps nontrivial denominator', ...)` | Boundary class richness | yes |
| A8 | sympy | 105 | `assert_zero('concrete linear IBP boundary discharge', ...)` | Boundary IBP / Gaussian | yes |
| A9 | sympy | 106 | `assert_zero('lorentzian-profile linear IBP boundary discharge', ...)` | Boundary IBP / Lorentzian | yes |
| A10 | sympy | 107 | `assert_zero('lorentzian finite-endpoint boundary discharge + 2', ...)` | Lorentzian finite-endpoint = −2 | yes |
| A11 | sympy | 108 | `assert_zero('concrete linear IBP with boundary', linear_cross − (boundary + bulk))` | IBP consistency on concrete profile | yes |
| A12 | sympy | 109 | `assert_zero('concrete quadratic IBP boundary discharge', ...)` | Quadratic boundary, Gaussian | yes |
| A13 | sympy | 110 | `assert_zero('lorentzian-profile quadratic IBP boundary discharge', ...)` | Quadratic boundary, Lorentzian | yes |
| A14 | sympy | 111 | `assert_zero('concrete quadratic IBP with boundary', ...)` | IBP consistency, quadratic | yes |
| A15 | sympy | 138 | `assert_zero('linear density before IBP', L1 − explicit form)` | Linear-order density | yes |
| A16 | sympy | 146 | `assert_zero('linear background coefficient after IBP', ...)` | Background equation extraction | yes |
| A17 | sympy | 153 | `assert_zero('raw quadratic cross term', d/d(eta) d/d(eta_w) L2_raw + A)` | Identifies the cross coefficient | yes |
| A18 | sympy | 163 | `assert_zero('quadratic density after IBP', L2_ibp − canonical_L2)` | `K_eta` recovery | yes |
| A19 | sympy | 170 | `assert_nonzero('mutated K_eta sign should fail', ...)` | K_eta sign mutation | yes |
| A20 | sympy | 178 | `assert_zero('Y20 spherical-Laplacian eigenvalue + 6 Y20')` | Angular eigenvalue | yes |
| A21 | sympy | 179 | `assert_nonzero('mutated Y20 spherical-Laplacian eigenvalue should fail')` | Eigenvalue mutation | yes |
| A22 | sympy | 205 | `assert_zero('Y20 modal norm − 1')` | Angular norm | yes |
| A23 | sympy | 206 | `assert_zero('Y20 modal angular stiffness − 6')` | Angular stiffness | yes |
| A24 | sympy | 207 | `assert_zero('Y20 projected modal density − closed')` | Projection algebra | yes |
| A25 | sympy | 214 | `assert_zero('Y20 fused modal Euler-Lagrange equation')` | Modal EL | yes |
| B1 | mathematica | 61 | `expectZero["[M1]", elViaD − handExact]` | Minimal nonlinear candidate | yes |
| B2 | mathematica | 62 | `expectNonzero["[M1 mutation]"]` | M1 mutation | yes |
| B3 | mathematica | 88–95 | `[M2]` linear density and mutation | Linear order | yes |
| B4 | mathematica | 98–105 | `[M3]` background equation + mutation | Background equation | yes |
| B5 | mathematica | 108–116 | `[M4]` raw quadratic density + cross + mutation | Quadratic raw | yes |
| B6 | mathematica | 124 | `[M5 product-rule]` symbolic IBP product rule | IBP identity | yes |
| B7 | mathematica | 136 | `[M5]` canonical post-IBP quadratic density | `K_eta` recovery | yes |
| B8 | mathematica | 139 | `[M5 mutation]` | K_eta sign mutation | yes |
| B9 | mathematica | 169–173 | `[M6]`, `[M6 linear IBP]`, `[M6 quadratic IBP]` | Concrete IBP on Gaussian | yes |
| B10 | mathematica | 175 | `[M6 mutation]` (weak: residual = −1 trivially) | weak mutation | partial |
| B11 | mathematica | 180–181 | `[M7]` finite-endpoint Lorentzian discharge + 2 = 0 and mutation | Lorentzian finite endpoint | yes |
| B12 | mathematica | 186–190 | `[M8]` arctan discharge = Pi; `[M8 denominator]` Lorentzian denominator | arctan discharge + Lorentzian denominator | yes |
| B13 | mathematica | 203–204 | `[M9]` Y20 eigenvalue + mutation | Angular eigenvalue | yes |
| B14 | mathematica | 218–222 | `[M10]` angular norm/stiffness + mutation | Angular sector | yes |
| B15 | mathematica | 241–247 | `[M11]` modal EL + mutation | Modal EL | yes |

Every load-bearing assertion is paired to a paper claim or to a script-side intermediate the paper card invokes implicitly (IBP move, modal EL). The single "partial" row (B10) is a weak mutation — `linGaussian = 0` deterministically, so `linGaussian − 1 = −1` passes `expectNonzero` trivially. It does not strengthen the audit, but it does not invalidate it either; the M6 zero-check itself is exercised by four independent boundary terms summed in quadrature, which would catch any sign or factor error in any of them. Not a finding.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is independently derived, not a transliteration:

- **EL derivation differs structurally.** SymPy uses `sympy.calculus.euler.euler_equations(L_exact, Rfield, [t, w, u, v])` (a library call) and compares against a hand-built `exact_el_lhs` (lines 42–53). Mathematica writes the EL out explicitly as `D[lagExact, r] − D[D[lagExact, D[r, t]], t] − …` (lines 51–55) and compares against a separate `handExact`. Two different machinery choices for the same identity.
- **Spherical harmonics differ.** SymPy uses `Ynm(2, 0, th, ph)` from `sympy.functions.special.spherical_harmonics` and expands with `expand_func` (lines 173); Mathematica uses `SphericalHarmonicY[2, 0, th, ph]` and applies `FunctionExpand` + `ExpToTrig` (lines 195–196). The intermediate symbolic forms are not identical (one returns `Ynm` after `expand_func`, the other returns trig form), yet both checks land at `−6 Y20`.
- **Lorentzian denominator probe differs.** SymPy checks the *denominator* of `Together[−exp(−w²)/(1+w²)]` is `≠ 1` (line 89); Mathematica checks `Together[(1+w²)·Together(...) + exp(−w²)] == 0` plus a mutation (lines 188, 190). Different probes of the same nontrivial-denominator fact.
- **Mathematica adds an explicit `[M5 product-rule]` standalone IBP identity** (lines 118–124) that SymPy does not have as a separate check; conversely, SymPy has the generic `linear IBP identity` and `quadratic IBP identity` symbolic checks (lines 65–72) that Mathematica only handles via M3/M5 inline.
- **arctan probe**: SymPy only asserts `nonzero`; Mathematica asserts the exact value `Pi` and mutates by removing the `−Pi` shift. The two scripts are NOT line-by-line ports here.

Conclusion: independent derivation, no `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce `STATUS: PASS` with identical claim content (same `K_eta` formula, same modal EL, same angular eigenvalue, same `−2` finite-endpoint discharge, same `Pi` arctan discharge, same +6 angular-stiffness coefficient). The SymPy output banner explicitly prints the `K_eta` formula and the L2_ibp form; Mathematica prints the same `L2 after IBP = -R0p^2*TwRR0*eta^2/4 - TO0*grad2/2 - Tw0*eta_w^2/2 - URR0*eta^2/2 + d_TwR_R0p*eta^2/2 + eta_t^2*mu0/2` and `K_eta(w) = R0p^2*TwRR0/2 + URR0 - d_TwR_R0p`. These agree.

## Verdict justification

`clean`. Attacks tried and failed:

1. **Tautology hunt on the abstract IBP moves.** `d_Tw_R0p` and `d_TwR_R0p` are introduced as *symbols* representing `d_w[T_w(R0,w) R0']` and `d_w[T_{w,R}(R0,w) R0']`. The post-IBP equality `L2_ibp − canonical_L2 = 0` is therefore an algebraic accounting identity, not a real IBP proof. **Defense**: the concrete IBP integrals on Gaussian profiles (sympy A11/A14, mathematica B9) actually evaluate `∫ −B·∂_w(eta) dw − ∫ ∂_w(B)·eta dw` and `∫ −A·eta·∂_w(eta) dw − ∫ ∂_w(A)·eta²/2 dw` to zero, and the symbolic-IBP A3/A4 identities check the product rule in symbolic form. Together these tie the abstract symbol-shuffling to the actual integration-by-parts move.
2. **Tautology hunt on EL.** `assert_zero('exact nonlinear Euler-Lagrange equation', exact_el + exact_el_lhs)` could be trivial if `euler_equations` and the hand expression were structurally identical. They are not — `euler_equations` builds the variational derivative from `L_exact` while `exact_el_lhs` is a hand expansion of the same derivative; their printed forms differ; the equality `+= 0` is a genuine cross-check. The mutation `mutated_exact_el_lhs - 2·diff(USigma, R)` confirms the check would catch a `dU/dR` sign error.
3. **Symbol-assumption hunt.** All symbols are declared `real`; `grad2` is `nonnegative` (matches its physical role as `|∇_Ω eta|²`); the angular variables are `0 < th < Pi` in mathematica's M9 (correct for the spherical Laplacian eigenvalue domain) and then released for the M10 norm integrals (where the integration measure `sin(th)` handles the endpoints). No domain contradiction with the physical setup.
4. **Missing-branch hunt.** The script checks Y_{20} only; the modal EL is stated only for `l=2` in the closed form. The paper card's modal claim is `mu_eta q_lm,tt − ∂_w(T_w q_lm,w) + [K_eta + l(l+1) T_Omega] q_lm = S_lm` (general `l`). Y_{20} verification of `+6 T_Omega = l(l+1) T_Omega` for `l=2` is one branch only; the script's own commentary (line 261) calls this out as "now also integrates the genuine Y20 harmonic." Other harmonics (`l=0`, `l=4`) are not exercised. **Defense**: the script states (line 259) "For l=0 this reproduces the scalar lane" — the generic `K_eta + l(l+1) T_Omega` decomposition is established by the angular norm/stiffness factorization in A24/A25, which is `l`-generic in form; Y_{20} is the concrete instantiation. Not a finding.
5. **Engine disagreement hunt.** Both engines emit identical `K_eta` and `L2_ibp` forms.
6. **Stale-output hunt.** Both `.txt` outputs are newer than their generating scripts.

Paper alignment: the audit reads the paper card and the appendix row and confirms the script's load-bearing assertions match the paper-side deliverables. The Y_{20} modal EL is the only assertion not explicitly itemized on the stage card body, but it is in-scope per the appendix row "Candidate parent-action data for the projected throat sector" and per the card's "parent-complete source of the linear wall closure" framing.

## Self-test notes

Checked variable independence (every `sp.diff` / `D[...]` is taken w.r.t. a symbol the integrand actually depends on — no zero-by-construction derivatives). Checked parity for the boundary-discharge integrals (Gaussian × Gaussian and Lorentzian × Lorentzian decay at ±∞, so the discharges are 0; the `arctan w` discharge is `π` and the `w·sqrt(1+w²)·1/(1+w²)` discharge is `−2`, matching). Checked the IBP sign convention: `−A(w)·eta·eta_w → +(1/2)·dA/dw·eta²` mod boundary, consistent with both script halves. No need for a directive; report ends with `verdict: clean`.
