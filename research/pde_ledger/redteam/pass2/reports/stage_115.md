---
unit_id: 115
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage115_core_balance_compensation.md]
  paper_appendix: present
---

# Audit unit 115 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_115.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage115_core_balance_compensation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: `\input{stages/stage_115}` line 1264; subsec `app-part04-core-balance` lines 511-573; eqs `core-canonical-conditions`, `core-coupling-balance`, `r-g-parent-ratios`, `parent-compensation-family`; anchor MTDC-T8.3 line 1177)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage115_core_balance_compensation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage115_core_balance_compensation_mathematica_audit.txt`

## What the paper claims

Stage 115 ("Exact Core-Balance Compensation Theorem") states the conditions under which the concrete two-channel core model of Stage 114 lands exactly on the nontrivial compensated canonical branch of Stage 112. The card's body quote is the bottom-line claim: "Coupling balance \(g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2\) realizes the compensated branch." The notes expand this into the full set of deliverables: (1) the canonical-branch conditions \(\rho_c=4\sigma_c,\ \kappa_c=1/3,\ \gamma_c=1/9\); (2) the exact coupling-balance law \(g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2\); (3) the explicit two-branch solution \(g_q=\frac{g_s}{2K_s}(2\lambda\pm\sqrt{K_sK_q+\lambda^2})\); (4) the bare-channel normalizations \(\kappa_0=(1+r_c)/3,\ \gamma_0=(1+r_c)/9\); (5) the exact collapse \(\delta\Lambda_{\rm core}(z)=4\sigma_*-\sigma_*/(1-z^2/3-iz^5/9)\) with \(\sigma_*=g_s^2/(4K_s)\); and (6) the preserved normalized outgoing fingerprint \(\widehat Y_2(z)=1+z^2/9+4z^4/81+iz^5/27+O(z^6)\). The appendix subsection `app-part04-core-balance` independently states the same balance law (eq:app-part04-core-coupling-balance) and gives the equivalent parent-overlap reparametrization \(\mathfrak r:=\lambda/\sqrt{K_sK_q},\ \mathfrak g:=g_q\sqrt{K_s}/(g_s\sqrt{K_q})\) with parent family \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,\ \mathfrak g_\pm=\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\) (eqs r-g-parent-ratios, parent-compensation-family).

## What the script claims to verify

The SymPy script defines \(r_c=\lambda^2/(K_sK_q)\), \(\rho_c=g_s^2/K_s\), \(\sigma_c=(K_sg_q-\lambda g_s)^2/(K_s^2K_q(1+r_c))\), forms `balance_eq = rho_c - 4 sigma_c`, solves it for `g_q` (recovering the two-branch law), and then asserts: (i) `sigma_c` on the balance surface equals `sigma_* = g_s^2/(4K_s)`; (ii) the concrete core deformation `rho_c - sigma_c/(1 - (kappa0/(1+r_c))z^2 - i(gamma0/(1+r_c))z^5)` collapses identically to the Stage-112 target branch when `kappa0=(1+r_c)/3, gamma0=(1+r_c)/9, g_q=` the solved branch; and (iii) the normalized fingerprint of `Lambda_out + target_delta` reproduces `Y_target = 1 + z^2/9 + 4z^4/81 + iz^5/27`. It also prints `sigma_*, kappa0, gamma0`. The Mathematica script verifies the same three identities AND adds an independent parent-overlap route: it builds `parentFamilyResidual = 1 + frakR^2 - 4(frakG-frakR)^2`, proves it equals `balanceEq` up to the nonzero factor `(kS kQ + lam^2)/(gS^2 kQ)`, solves the parent family for a fresh variable `gVar`, matches the closed-form root `frakGMinus = frakR - Sqrt[1+frakR^2]/2` against the Solve output, translates back to `g_q`, and re-verifies `sigma_c = sigma_*`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Balance law \(g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2\) | `balance_eq = rho_c - 4 sigma_c`; numerator zeroed (py L28-32; wl L37-42); printed factored form in both outputs | match |
| Two-branch \(g_q=\frac{g_s}{2K_s}(2\lambda\pm\sqrt{K_sK_q+\lambda^2})\) | `solve(balance_eq,g_q)` (py L32); `Solve[...,gQ]` (wl L40); both outputs print the two roots | match |
| \(\rho_c=4\sigma_c\) i.e. \(\sigma_c=\sigma_*=g_s^2/(4K_s)\) on surface | `expect_zero("sigma_c on balance surface", sigma_c.subs(...)-sigma_star)` (py L41-44; wl L46) | match |
| \(\kappa_0=(1+r_c)/3,\ \gamma_0=(1+r_c)/9\) | `kappa0_can=(1+r_c)/3, gamma0_can=(1+r_c)/9` (py L46-47; wl L88-89); used in collapse; printed | match |
| Exact collapse \(\delta\Lambda_{\rm core}=4\sigma_*-\sigma_*/(1-z^2/3-iz^5/9)\) | `expect_zero("exact collapse...", delta_core - target_delta)` (py L49-61; wl L91-97) | match |
| Fingerprint \(\widehat Y_2=1+z^2/9+4z^4/81+iz^5/27\) | `expect_zero("normalized outgoing fingerprint preserved", ...)` (py L63-70; wl L99-103) | match |
| Parent family \(1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2\), \(\mathfrak g_\pm=\mathfrak r\pm\tfrac12\sqrt{1+\mathfrak r^2}\) (appendix) | `parentFamilyResidual` identity + `frakG_-` root match (wl L48-86) | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `expect_zero(sigma_c|surface - sigma_star)` | \(\sigma_c=\sigma_*\) on balance surface (claim 1/3) | yes |
| A2 | sympy | 49-61 | `expect_zero(delta_core - target_delta)` | exact collapse to Stage-112 branch (claim 5) | yes |
| A3 | sympy | 63-70 | `expect_zero(series(Y_eff) - Y_target)` | preserved fingerprint (claim 6) | yes |
| A4 | sympy | 32-36 | `solve(balance_eq,g_q)` (printed, not asserted) | two-branch law (claim 3) | partial (print) |
| B1 | math | 41 | `If[Length[gQSolutions]=!=2, fail]` | exactly two branches (claim 3) | yes |
| B2 | math | 46 | `expectZero(sigma_c|gQBranch - sigmaStar)` | \(\sigma_c=\sigma_*\) (claim 1/3) | yes |
| B3 | math | 52-59 | `expectZero(parentFamilyResidual - balanceEq·factor)` | balance law ≡ parent family (independent route) | yes |
| B4 | math | 64-65 | `If[Length[frakGRoots]=!=2, fail]` | two parent-family roots | yes |
| B5 | math | 69-75 | `expectZero((root1-frakGMinus)(root2-frakGMinus))` | closed-form \(\mathfrak g_-\) matches Solve | yes |
| B6 | math | 80-86 | `expectZero(sigma_c|gQFromFrakMinus - sigmaStar)` | \(\sigma_*\) via independent reparam (claim 1/3, indep) | yes |
| B7 | math | 91-97 | `expectZero(deltaCore - targetDelta)` | exact collapse (claim 5) | yes |
| B8 | math | 99-103 | `expectZero(series(yEff) - yTarget)` | preserved fingerprint (claim 6) | yes |

All `expectZero`/`expect_zero` checks compare two independently-constructed expressions (a derived/substituted form vs. a separately-written target). None resubstitutes its own solved root into the equation it was solved from in a way that guarantees zero by construction: e.g. A1/B2 substitute the *solved* `g_q` branch into `sigma_c` and compare against the *separately defined* `sigma_* = g_s^2/(4K_s)` — that is a genuine consequence of the balance condition, not an identity. A3/B8 compare a Taylor series of a derived rational function against a separately-typed target polynomial. No tautology.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a **genuine independent route**, not a transliteration. The shared three identities (sigma_c-on-surface, collapse, fingerprint) are unavoidably the same closed forms, but the `.wl` reaches the load-bearing balance result by a structurally different algebra than the `.py`:

- SymPy (py L32): `gq_solutions = sp.solve(sp.Eq(balance_eq, 0), g_q)` — directly solves the balance equation in the original `(K_s,K_q,λ,g_s,g_q)` variables.
- Mathematica (wl L48-86): introduces the parent-overlap reparametrization `frakR = lam/Sqrt[kS*kQ]`, `frakG = gQ*Sqrt[kS]/(gS*Sqrt[kQ])`, forms `parentFamilyResidual = 1 + frakR^2 - 4*(frakG - frakR)^2`, **proves** this equals `balanceEq*(kS*kQ + lam^2)/(gS^2*kQ)` (wl L52-59), then solves the *parent family* for a fresh symbol `gVar` (wl L62-63), matches the hand-derived closed form `frakGMinus = frakR - Sqrt[1+frakR^2]/2` against the Solve output (wl L67-75), and finally translates `frakGMinus` back to `g_q` and re-derives `sigma_c = sigma_*` (wl L76-86).

This reparametrized derivation is exactly the appendix's `eq:app-part04-r-g-parent-ratios` / `eq:app-part04-parent-compensation-family` path (the `.wl` even cites these eq labels in its comment, wl L48-49). It is independent algebra reaching the same physical conclusion, which is precisely what the second-engine policy requires. The `.wl` does NOT reuse the `.py`'s `solve(balance_eq, g_q)` choreography. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at every comparable checkpoint. Side-by-side:

- balance numerator: py prints the expanded `4 K_s^2 g_q^2/... + g_s^2/K_s` form; wl prints factored `(gS^2*kQ*kS - 4*gQ^2*kS^2 + 8*gQ*gS*kS*lam - 3*gS^2*lam^2)/(kS*(kQ*kS + lam^2))`. These are the same rational expression (numerator = `g_s^2(K_sK_q+λ^2) - 4(K_sg_q-λg_s)^2` algebraically; verified by hand-expansion).
- two `g_q` roots: identical in both outputs (`(gS*lam)/kS ± Sqrt[gS^2*kQ*kS + gS^2*lam^2]/(2*kS)` = `g_s/(2K_s)(2λ ± √(K_sK_q+λ^2))`).
- `sigma_c on balance surface = 0`: both.
- `exact collapse = 0`: both.
- `normalized outgoing fingerprint preserved = 0`: both.
- `sigma_* = g_s^2/(4K_s)`, `kappa0 = (1+r_c)/3`, `gamma0 = (1+r_c)/9`: both (py prints `(K_q*K_s + lam**2)/(3*K_q*K_s)`, wl prints the equivalent `(1 + lam^2/(kQ*kS))/3`).

No `engine_disagreement`.

## Verdict justification

`clean`. The paper card (terse body quote + Inputs/Verification/Checks scaffolding) and the notes/appendix fully describe six deliverables; the scripts verify all six with non-tautological zero-residual checks, and the Mathematica script additionally supplies a genuinely independent parent-overlap derivation of the load-bearing balance result. Attacks tried and failed: (1) tautology hunt — every `expectZero` compares a derived/substituted form against a separately defined target, none resubstitutes its solved root into its own source equation; (2) transliteration hunt — the `.wl` uses a distinct reparametrized algebra (parent family) plus an independent root match, not a port of the `.py`'s direct `solve`; (3) symbol-domain hunt — `K_s,K_q>0` and `kappa0,gamma0>0` are consistent with the physical setup (positive geometric/coupling constants; `g_s,g_q` correctly left unrestricted-real since the balance solution can take either sign), no assumption masks a branch; (4) freshness — both `.txt` outputs are newer than their scripts (sympy out 1779921486 > py 1779921193; math out 1779921812 > wl 1779921802) and their content matches the current scripts; (5) value reconciliation — all emitted deliverable values appear correctly in the notes/appendix (see below). I read the card, the notes, and the relevant appendix subsection; the scripts' claims match the paper's claims.

## Self-test notes

Checked traps: (1) variable-independence — no `diff`/`D` derivatives in this stage, so the zero-derivative trap does not apply; the substitutions (`g_q -> gq_branch`, `kappa0/gamma0 -> _can`) put the checked expressions into genuine forms that depend on all of `K_s,K_q,λ,g_s,z`. (2) The series-truncation comparisons (A3/B8) compare degree-≤5 expansions of a rational `Y_eff` against a separately-typed degree-5 polynomial — confirmed the target `1+z^2/9+4z^4/81+iz^5/27` is what the notes box. (3) Trivial-case: substituting the solved `g_q` branch reduces `sigma_c` to `g_s^2/(4K_s)` as a real algebraic consequence (verified by hand-expanding the balance numerator), not a built-in zero. No directive written (zero findings).

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation of every RESULT value the scripts emit:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| balance law `g_s^2(K_sK_q+λ^2)=4(K_sg_q-λg_s)^2` | py L28 `balance_eq`, wl L37-38; sympy out L5-11, math out L5 | card L16 (body quote); notes L43-48 (boxed); appendix eq:app-part04-core-coupling-balance L524-529 | MATCH |
| two-branch `g_q = (g_s/2K_s)(2λ ± √(K_sK_q+λ^2))` | py L32-36, wl L40-42; sympy out L12-22, math out L6 | notes L50-58 (boxed); appendix parent-compensation-family roots L538-543 | MATCH |
| `sigma_* = g_s^2/(4K_s)` | py L40,73 `sigma_star`; wl L45,107; sympy out L28, math out L22 | notes L83 (boxed `σ_*=g_s^2/4K_s`) | MATCH |
| `kappa0 = (1+r_c)/3` | py L46,74; wl L88,108; sympy out L29, math out L23 | notes L67 (boxed `κ_0=(1+r_c)/3`); appendix κ_c=1/3 family | MATCH |
| `gamma0 = (1+r_c)/9` | py L47,75; wl L89,109; sympy out L30, math out L24 | notes L69 (boxed `γ_0=(1+r_c)/9`) | MATCH |
| canonical conditions `ρ_c=4σ_c, κ_c=1/3, γ_c=1/9` | py L28 (`rho_c-4 sigma_c`), L46-47; wl L37,88-89 | notes L20-27 (boxed); appendix eq:app-part04-core-canonical-conditions L515-522 | MATCH |
| exact collapse `δΛ_core = 4σ_* - σ_*/(1-z^2/3-iz^5/9)` | py L58-61 `target_delta`; wl L96-97; sympy out L24, math out L16-17 | notes L77-84 (boxed) | MATCH |
| fingerprint `Ŷ_2 = 1 + z^2/9 + 4z^4/81 + iz^5/27` | py L66 `Y_target`; wl L102; sympy out L25, math out L18-19 | notes L93-95 (boxed) | MATCH |
| parent family `1+𝔯^2=4(𝔤-𝔯)^2`, `𝔤_- = 𝔯 - ½√(1+𝔯^2)` | wl L52,62-67 (independent route); math out L9-13 | appendix eq:app-part04-parent-compensation-family L538-543, r-g-parent-ratios L531-536 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `r_c = λ^2/(K_sK_q)` (intermediate ratio, defined in py L24 / wl L33; appears as `(1+r_c)` factor inside the printed κ0/γ0 forms anyway); `rho_c = g_s^2/K_s`, `sigma_c` (intermediate definitions); `balanceEq*(kS*kQ+lam^2)/(gS^2*kQ)` proportionality factor (wl L57, verification scaffolding); `gQBranch`/`gQFromFrakMinus`/`frakGMinus`/`frakGRoots`/`gVar` (intermediate solve handles); `Lambda_out`/`Lambda_eff`/`lambdaEff`/`Y_eff`/`yEff` (intermediate series objects); all PASS/`= 0` residual flags.

reconciliation: complete; 9 deliverable values checked, 0 misaligned.
