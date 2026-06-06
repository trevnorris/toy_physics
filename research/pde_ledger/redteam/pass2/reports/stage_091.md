---
unit_id: 091
batch: IV.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md]
  paper_appendix: present
---

# Audit unit 091 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_091.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows/eqs referencing stage 091: lines 10-25, 90-153, 218-245, 1175)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.txt`

## What the paper claims

Stage 091 is a geometry-lane firewall ledger step (checkpoint: False). Its `\stagefield{Derivation ledger}` and the quoted Output box state: "Static-geometry plus one-pole grouped carrier forces the \(3/4+1/4\) conservative quadrupole module and hence \(\rho_\alpha=4/3\), \(\zeta_{\rm req}=1/3\)." The card names the forced carrier `\widehat Y_Q^{\rm cons}=3/4+(1/4)/(1-\omega^2/\Omega_Q^2)`. The notes make the derivation explicit: take `K_Q^cons(ω)=K_geom + K_pole/(1-ω²/Ω_Q²)` with static (ω-independent) `K_geom`; expand to get `K0=K_geom+K_pole`, `K2=K_pole/Ω_Q²`, `K4=K_pole/Ω_Q⁴`; insert into the minimal isotropic branch identity `K0·K4=4·K2²`; with `K_pole≠0` this forces `K_geom=3K_pole`, i.e. `K_pole=K0/4`, `K_geom=3K0/4`, hence the `3/4+1/4` module and the corollary loading ratios `rho_alpha=alpha_req/alpha_mix=4/3`, `zeta_req=(alpha_req-alpha_mix)/alpha_mix=1/3`. The appendix (eqs app-part04-geometry-static-split, -kbar-cons) corroborates `K_{g,0}=3K_{\rm pole}`, `K_{\rm pole}=K_0/4`, `K_{g,0}=3K_0/4`, and the same module form. Distinct deliverables: (1) `K_geom=3K_pole` / `K_pole=K0/4`; (2) `Yhat=3/4+1/4/(1-ω²/Ω_Q²)`; (3) `rho_alpha=4/3`, `zeta_req=1/3`. Card Checks: static limit ε2=ε4=0 returns c_pole=1/4; l=0/l=2 orthogonality before firewall; support/source success carries the minimal-module hypothesis.

## What the script claims to verify

Both scripts construct `Kcons = Kgeom + Kpole/(1-ω²/Ω_Q²)` with `Kgeom` a pure constant (encoding the static-geometry hypothesis K_{g,2}=K_{g,4}=0), series-expand to O(ω⁴), and read off K0,K2,K4. They confirm those coefficients match the notes' expressions, compute the branch-identity residual `K0·K4-4·K2²` (which is the non-zero `K_pole(K_geom-3K_pole)/Ω_Q⁴`), solve it for `Kgeom`, and assert the solution equals `3*Kpole`. On that branch they assert `K0=4*Kpole`, normalize and assert `Yhat = 3/4 + 1/4/(1-ω²/Ω_Q²)`, and assert `rho_alpha=K0/Kgeom=4/3`, `zeta_req=(K0-Kgeom)/Kgeom=1/3`. The Mathematica script adds an independent partial-fraction recombination path (starting from `K_geom=3K_pole`, recombining via `Together` rather than Series+Solve) for the `Yhat` deliverable. The SymPy docstring annotates that the l=0/l=2 orthogonality check (card Check item 2) is carried forward from stage 094.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Series coeffs K0=K_geom+K_pole, K2=K_pole/Ω², K4=K_pole/Ω⁴ | `expect_zero` on each (py 54-56; wl 43-45) | match |
| Branch identity forces K_geom=3K_pole | solve(K0K4-4K2²=0,Kgeom); assert =3Kpole (py 58-63; wl 47-53) | match |
| K_pole=K0/4 (c_pole=1/4) | `K0_on_branch - 4*Kpole == 0` (py 67; wl 56) | match |
| Yhat = 3/4 + 1/4/(1-ω²/Ω²) | assert Yhat-Yhat_expected==0 (py 69-72; wl 58-61) + independent recomb (wl 74-78) | match |
| rho_alpha = 4/3 | `rho_alpha - 4/3 == 0` (py 78; wl 67) | match |
| zeta_req = 1/3 | `zeta_req - 1/3 == 0` (py 79; wl 68) | match |
| Check: static limit ε2=ε4=0 → c_pole=1/4 | static case baked in (Kgeom const); c_pole=1/4 ≡ K0=4Kpole assertion | match (implicit; dynamic-ε obstruction is a separate appendix analysis, not this card's target) |
| Check: l=0/l=2 orthogonality before firewall | carried forward to stage 094 (docstring 11-15; card Inputs authorise the scalar ansatz) | match (legitimate carry-forward) |

Dominant pattern: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `K0-(Kgeom+Kpole)==0` | series coeff K0 | yes |
| A2 | sympy | 55 | `K2-Kpole/Ω²==0` | series coeff K2 | yes |
| A3 | sympy | 56 | `K4-Kpole/Ω⁴==0` | series coeff K4 | yes |
| A4 | sympy | 63 | `Kgeom_sol-3Kpole==0` | claim 1 (K_geom=3K_pole) | yes |
| A5 | sympy | 67 | `K0_on_branch-4Kpole==0` | claim 1 (K_pole=K0/4, c_pole=1/4) | yes |
| A6 | sympy | 72 | `Yhat-[3/4+1/4/(1-ω²/Ω²)]==0` | claim 2 | yes |
| A7 | sympy | 78 | `rho_alpha-4/3==0` | claim 3 | yes |
| A8 | sympy | 79 | `zeta_req-1/3==0` | claim 3 | yes |
| A9 | mathematica | 43-45 | `expectZero` on K0,K2,K4 | series coeffs | yes |
| A10 | mathematica | 53 | `kGeomForced-3kPole==0` | claim 1 | yes |
| A11 | mathematica | 56 | `k0OnBranch-4kPole==0` | claim 1 | yes |
| A12 | mathematica | 61 | `yHat-yHatExpected==0` | claim 2 | yes |
| A13 | mathematica | 67-68 | `rhoAlpha-4/3`, `zetaReq-1/3` | claim 3 | yes |
| A14 | mathematica | 78 | `yHatRecomb-yHatTargetRecomb==0` (independent recomb path) | claim 2 (independent route) | yes |

## Findings

(none)

## Independent-derivation check (Mathematica)

The `.wl` mirrors the SymPy choreography for the main block (Series → Coefficient → Solve), which is legitimate cross-engine confirmation. Crucially it does NOT stop there: lines 70-78 add a genuinely independent algebraic route to the central `Yhat` deliverable that bypasses both Series and Solve — it starts from the branch-identity output `K_geom=3K_pole`, builds `kConsBranchDirect = 3*kPole + kPole/(1-omega^2/omegaQ^2)`, normalizes by `k0BranchDirect = 4*kPole`, recombines with `Together`, and checks it equals `Together[3/4 + (1/4)/(1-omega^2/omegaQ^2)]`. This is an algebraically distinct path (rational recombination vs. series+solve) confirming the same bottom line, so the second engine is not a bare transliteration. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical bottom-line results (residuals all 0):
- SymPy (txt 13): `Branch identity = K_pole*(K_geom - 3*K_pole)/Omega_Q**4`; Mathematica (txt 16): `((kGeom - 3*kPole)*kPole)/omegaQ^4` — same nonzero residual.
- SymPy (txt 18): `Yhat = (Omega_Q**2 - 3*omega**2/4)/(Omega_Q**2 - omega**2)`; Mathematica (txt 22): `(3 + (1 - omega^2/omegaQ^2)^(-1))/4` — algebraically equal forms.
- Both: `rho_alpha=4/3`, `zeta_req=1/3`, all `expect_zero`/`expectZero` checks pass; Mathematica exits 0 with "Stage 091 Mathematica audit passed."
Engines agree.

## Verdict justification

Verdict: **clean**. I read the card, notes, and the part-04 appendix firewall block before opening the scripts. Attacks tried and failed: (1) Tautology on the branch identity — the residual `K_pole(K_geom-3K_pole)/Ω⁴` is genuinely nonzero (output line 13), so the solve-for-Kgeom is a real constraint, not a no-op; the K_geom=3K_pole assertion can fail. (2) Hardcoded answer — Yhat_expected, 4/3, 1/3 are written as independent targets and the LHS is derived through substitution of the *solved* Kgeom, not pinned. (3) Symbol-domain masking — `omega` real, `Kgeom,Kpole,OmegaQ` positive&real; positivity is exactly the notes' "K_pole≠0 nontrivial branch" assumption (notes line 101) and masks no branch in the rational/series algebra. (4) Static-geometry baked in — `Kgeom` being a pure constant *is* the stated hypothesis K_{g,2}=K_{g,4}=0 (appendix eq -geometry-static-split is derived under exactly that), so it is faithful, not a hidden cheat; the dynamic-ε obstruction is a separate appendix analysis the card does not assign to this stage. (5) rho_alpha/zeta_req identification — `K_geom↔alpha_mix`, `K0↔alpha_req` reproduces the notes' `alpha_mix/alpha_req=3/4` exactly, so 4/3 and 1/3 are the right quantities. The l=0/l=2 orthogonality Check is a documented carry-forward to stage 094 (authorised by the card Inputs/scalar-ansatz). Both engines present, agree, and outputs are fresh. Every assertion traces to a paper deliverable and is non-tautological.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 8 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| K0 = K_geom + K_pole | py:54, wl:43, out py:7 / wl:7 | notes:82 | MATCH |
| K2 = K_pole/Ω_Q² | py:55, wl:44, out py:8 / wl:8 | notes:84 | MATCH |
| K4 = K_pole/Ω_Q⁴ | py:56, wl:45, out py:9 / wl:9 | notes:86 | MATCH |
| K_geom = 3·K_pole | py:63, wl:53, out py:14 / wl:17 | notes:107; appendix:119 (`K_{g,0}=3K_{\rm pole}`) | MATCH |
| K_pole = K0/4 (K0=4K_pole) | py:67, wl:56, out py:16 / wl:20 | notes:111; appendix:121 (`K_{\rm pole}=\frac{K_0}{4}`); card:13 (c_pole=1/4 in Checks) | MATCH |
| Yhat = 3/4 + (1/4)/(1-ω²/Ω_Q²) | py:72, wl:61/78, out py:18 / wl:22 | notes:119; appendix:225; card:13,16 | MATCH |
| rho_alpha = 4/3 | py:78, wl:67, out py:20 / wl:25 | notes:149; card:16 | MATCH |
| zeta_req = 1/3 | py:79, wl:68, out py:21 / wl:26 | notes:153; card:16 | MATCH |

INTERNAL (scaffolding, no finding): branch-identity residual `K_pole(K_geom-3K_pole)/Ω_Q⁴` (intermediate that drives the solve, also exact in both outputs); series intermediate `Series=...`; pass/fail flags; banner/ledger print lines.

All eight stage-level deliverable values emitted by the scripts reconcile against the notes and/or the stage card and part-04 appendix. No MISMATCH, no MISSING-DELIVERABLE. Outputs are fresh (both .txt mtimes newer than their scripts), so the reconciliation rests on both source and committed output.

## Self-test notes

Checked: (1) variable-independence — no `diff`/`D` derivatives present; series expansion is in `omega` and `Kcons` genuinely depends on `omega`, so coefficients K2,K4 are nonzero literals (confirmed in outputs). (2) symmetry/parity — no unbounded integrals in this stage. (3) trivial-case pre-check — branch-identity residual is nonzero (so the solve is load-bearing) and substituting the solved `K_geom=3K_pole` reduces Yhat to the 3/4+1/4 target (confirmed against output lines 18/22). (4) paper round-trip — every assertion target (3K_pole, 4K_pole, 3/4+1/4 module, 4/3, 1/3) is exactly the value stated in notes/card/appendix; no new constant introduced. Conclusion: clean, no directive needed.
