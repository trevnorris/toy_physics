---
unit_id: 177
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md"]
  paper_appendix: present
---

# Audit unit 177 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_177.tex`
- notes: `/var/projects/toy_physics/../notes/stages/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md` (present — `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage.md`)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 85, 644-653)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage177_weak_axisymmetric_outgoing_slippage_mathematica_audit.txt`

## What the paper claims

Stage 177 collapses the Stage-176 outgoing-slippage bundle (the three microscopic port slippages \(\mathcal M_r,\mathcal I_r,\mathcal H_r\)) to a single weak-axisymmetric scalar amplitude. The card's `\stagefield{Output}` reads verbatim: "Shows all weak-axisymmetric outgoing slippages inherit the same grouped signature and collapse to \(\Xi_1=P_1/P_0\)." The notes enumerate the substantive deliverables: (1) each microscopic slippage inherits the lane signature with exact slopes \(\mathfrak m_r=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\tfrac12\kappa_1\), \(\mathfrak i_r=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r}\), \(\mathfrak h_r=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r}\); (2) the portwise defect collapses to \(\Sigma_{A,r}=\epsilon\lambda_A\sigma_r\) with \(\sigma_r=2\mathfrak m_r+\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\); (3) the grouped trace/anomaly law \(\bar\Xi=0,\ a_\Xi=\tfrac{\epsilon}{4}\Xi_1,\ b_\Xi=\tfrac{3\epsilon}{4}\Xi_1\Rightarrow b_\Xi=3a_\Xi\); (4) the physical identification \(P_1/P_0=\Xi_1\) (via Stage-173 \(P_1/P_0=\Xi_{\rm load}\)). The fixed grouped signature is \(\lambda_{20}=1,\lambda_{21}=\tfrac12,\lambda_{22}=-1\). Sections 8-9 are corollaries (rigidity / dominant-port limits) that follow algebraically and carry no new boxed identity.

## What the script claims to verify

The SymPy docstring lists four checks matching the deliverables: weak-axisymmetric grouped slopes of \(M_r,I_r,H_r\); the portwise defect collapse \(\Sigma_{A,r}=\epsilon\lambda_A\sigma_r\); the grouped trace/anomaly law for \(\Xi_A\); and the outgoing-prefactor-slope identification. Concretely the script builds each microscopic invariant from primitive symbols, perturbs each primitive by \(e^{\epsilon\lambda\,(\text{slope})}\), takes the order-\(\epsilon\) log-slope by `series(... eps, 0, 2)` (SymPy) / `D[...,eps]/.eps->0` (Mathematica), and asserts the slope equals \(\lambda\) times the claimed closed-form slope. The defect amplitude is checked by log-slope of \(\Lambda^2/K\) against \(\lambda\sigma_r\). The trace/anomaly block hard-instantiates \(\lambda_{20,21,22}=1,\tfrac12,-1\) and checks the spherical-harmonic projections. The prefactor block builds \(P_A=N_A/D_A\) and checks its slope equals \(P_0\,\Xi_{\rm load}\) with \(\Xi_{\rm load}=n_1-d_1\).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| \(\mathfrak m_r,\mathfrak i_r,\mathfrak h_r\) slopes inherit lane signature | py 67-69 / wl 57-59 (`d ln M/I/H = lam·{m,i,h}_r`) | match |
| Portwise defect \(\Sigma_{A,r}=\epsilon\lambda_A\sigma_r\) | py 89 / wl 70 (`Sigma - lam·sigma_r`) | match |
| Trace/anomaly \(\bar\Xi=0,\ a_\Xi=\tfrac{\epsilon}{4}\Xi_1,\ b_\Xi=\tfrac{3\epsilon}{4}\Xi_1,\ b_\Xi=3a_\Xi\) | py 106-109 / wl 86-89 | match |
| Identification \(P_1/P_0=\Xi_1\) (=\(\Xi_{\rm load}\)) | py 121-122 / wl 100-101 | match |
| \(\Lambda^2/K=M^2(1+I)^2/(1-H)^2\) factorization | wl 62-63 (extra, sympy uses \(\Lambda\) directly) | extra (independent-route scaffold, harmless) |

`paper_alignment: aligned` — every boxed deliverable in the notes maps to a substantive, non-tautological script check, and the wl `extra` check is a derivation-route bridge, not an unmotivated assertion.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 67 | `simplify(dlnM_exact - lam*m_r)==0` | slope \(\mathfrak m_r\) | yes |
| A2 | sympy | 68 | `... dlnI_exact - lam*i_r ==0` | slope \(\mathfrak i_r\) | yes |
| A3 | sympy | 69 | `... dlnH_exact - lam*h_r ==0` | slope \(\mathfrak h_r\) | yes |
| A4 | sympy | 89 | `Sigma_exact - lam*sigma_r ==0` | portwise defect \(\sigma_r\) | yes |
| A5 | sympy | 106 | `Xi_bar ==0` | \(\bar\Xi=0\) | yes |
| A6 | sympy | 107 | `a_Xi - eps Xi1/4 ==0` | \(a_\Xi\) | yes |
| A7 | sympy | 108 | `b_Xi - 3 eps Xi1/4 ==0` | \(b_\Xi\) | yes |
| A8 | sympy | 109 | `b_Xi - 3 a_Xi ==0` | \(b_\Xi=3a_\Xi\) | partial* |
| A9 | sympy | 121 | `P_slope_exact - P0*Xi_load ==0` | \(P_A\) slope | yes |
| A10 | sympy | 122 | `P_slope/P0 - Xi_load ==0` | \(P_1/P_0=\Xi_{\rm load}\) | yes |
| B1-B3 | math | 57-59 | `expectZero[d ln {M,I,H} - lam·{m,i,h}R]` | slopes | yes |
| B4 | math | 62-63 | `expectZero[lambda0^2/k - M^2(1+I)^2/(1-H)^2]` | factorization bridge | yes |
| B5 | math | 70 | `expectZero[sigmaExact - lam*sigmaR]` | portwise defect | yes |
| B6-B9 | math | 86-89 | trace/anomaly projections | \(\bar\Xi,a_\Xi,b_\Xi,b=3a\) | yes |
| B10-B11 | math | 100-101 | prefactor slope identification | \(P_1/P_0=\Xi_{\rm load}\) | yes |

\*A8 is implied by A6+A7 (\(b=3a\) once both are pinned to \(\Xi_1\)); it is a true consequence, not tautological, and is also independently re-derivable from the raw projections — acceptable as a redundant sanity check, not a finding.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent route**, not a transliteration, for the load-bearing defect check. The decisive evidence is the defect-amplitude block:

- SymPy (py 73-82) computes the slope of \(\Lambda^2/K\) where `Lambda = (OU2*GW + R*GU)/(OU2*OW2 - R**2)` (py 42) — i.e. directly from the algebraic closed form of the load factor, via `sp.series(sp.log((Lambda_p**2/Kp)/(Lambda**2/K)), eps, 0, 2)`.
- Mathematica (wl 66-68) computes the slope of `lambdaSqOverKP = mCalP^2*(1 + iCalP)^2/(1 - hCalP)^2` — i.e. from the **factorized** M/I/H representation, never forming `lambdaP` (it is defined at wl 47 but is NOT used in the sigma block).

These are two genuinely different symbolic objects that must be proven equal; the `.wl` makes that equality an explicit *separate* assertion (wl 62-63, `lambda0^2/k - mCal^2 (1+iCal)^2/(1-hCal)^2 == 0`) which has no SymPy counterpart. So Mathematica reaches \(\sigma_r\) by a path SymPy does not take, and additionally certifies the bridge. Elsewhere the two scripts share the same primitive-symbol perturbation scaffold (`Exp[eps*lam*slope]`) and the same closed-form slopes — unavoidable, since both must model the same physical setup — but the engines differ in the slope-extraction primitive (SymPy `series().removeO()/eps` vs. Mathematica `D[...,eps]/.eps->0`) and in the \(\Lambda^2/K\) construction. Not a port.

## Engine cross-check

Both transcripts report identical results for every shared check: `d ln M/I/H = 0`, `Sigma_{A,r}=lambda_A sigma_r = 0`, `grouped trace vanishes = 0`, `a_Xi - eps Xi1/4 = 0`, `b_Xi - 3 eps Xi1/4 = 0`, `b_Xi - 3 a_Xi = 0`, `P_A slope = P0*Xi_load = 0`, `(P1/P0) - Xi_load = 0`. The Mathematica transcript additionally shows `load-factor factorization ... = 0` (its extra bridge). All `PASS`; `Stage 177 Mathematica audit passed.` Engines agree.

## Verdict justification

`clean`. I read the card, the notes (8 sections of boxed identities), and the part-05 appendix rows (85, 644-653) before the scripts. Every boxed deliverable in the notes — the three microscopic slopes, the portwise \(\sigma_r\), the grouped trace/anomaly law, and the \(P_1/P_0=\Xi_1=\Xi_{\rm load}\) identification — has a matching, substantive, non-tautological assertion in BOTH engines, and the constants \(\lambda_{20,21,22}=1,\tfrac12,-1\) used in the trace block match the notes (§5) and appendix exactly. Attacks tried that failed: (a) tautology hunt — the slopes are extracted by independent series/derivative from `Exp`-perturbed primitives and compared to *separately written* closed forms, so a wrong closed form would not cancel; (b) `Lambda` round-trip — the two engines construct \(\Lambda^2/K\) by different algebra (closed form vs. M/I/H factorization) and the wl proves the bridge, so the sigma check is not "X−X"; (c) symbol-domain attack — `K,OU2,OW2,R,GU,GW` are `positive=True` in both engines, which is required for `log` and `sqrt(K)` to be single-valued and matches the physical setup (squared frequencies, magnitudes); no positivity is smuggled onto the perturbation slopes (`real` only), so no branch is hidden; (d) the `P` block uses `n1-d1` as \(\Xi_{\rm load}\) and recovers \(P_1/P_0\), matching the notes §6 Stage-173 carry. Outputs are fresh (txt mtimes 01:38-01:39 > script mtimes 01:09).

## Value Reconciliation (pass-2 augmentation)

The scripts emit only **symbolic closed-form** deliverables (this is an exact-closure structural stage; there are no numeric constants, benchmarks, or figures-of-merit). Enumerating every labeled/boxed result the scripts assert and locating it in the card / notes / appendix:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\mathfrak m_r=\mathfrak g_{W,r}-\mathfrak o_{W,r}-\tfrac12\kappa_1\) | py 58 / wl 49 / out "m_r = gW1 - oW1 - kappa1/2" | md:124 (boxed) | MATCH |
| \(\mathfrak i_r=\mathfrak r_r+\mathfrak g_{U,r}-\mathfrak o_{U,r}-\mathfrak g_{W,r}\) | py 59 / wl 50 / out "i_r = ..." | md:131 (boxed) | MATCH |
| \(\mathfrak h_r=2\mathfrak r_r-\mathfrak o_{U,r}-\mathfrak o_{W,r}\) | py 60 / wl 51 / out "h_r = ..." | md:138 (boxed) | MATCH |
| \(\sigma_r=2\mathfrak m_r+\frac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\frac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\) | py 84-88 / wl 69 / out "sigma_r = ..." | md:165-172 (boxed) | MATCH |
| \(\Xi_{\rm load}^{(A)}=\epsilon\lambda_A\Xi_1\) | py 98-100 / wl 78-80 / out "Xi_A = eps lambda_A Xi1" | md:187-192 / tex:15 / appx eq:app-part05-Xi1-definition (644-648) | MATCH |
| signature \(\lambda_{20}=1,\lambda_{21}=\tfrac12,\lambda_{22}=-1\) | py 94-96 / wl 75-77 | md:36-38, md:211-217 / appx | MATCH |
| \(\bar\Xi=0\) | py 102,106 / wl 82,86 / out "abar=0" | md:236 | MATCH |
| \(a_\Xi=\tfrac{\epsilon}{4}\Xi_1\) | py 103,107 / wl 83,87 / out "a=eps Xi1/4" | md:237 | MATCH |
| \(b_\Xi=\tfrac{3\epsilon}{4}\Xi_1\) | py 104,108 / wl 84,88 / out "b=3 eps Xi1/4" | md:239 | MATCH |
| \(b_\Xi=3a_\Xi\) | py 109 / wl 89 / out (defect law) | md:243 (boxed) | MATCH |
| \(\Xi_1=P_1/P_0\) | py 121-122 / wl 100-101 / out "Xi1 = P1/P0" | md:67-69, 261-263 / tex:15 / appx eq:app-part05-Xi1-P1P0 (650-652) | MATCH |
| \(\Lambda^2/K=M^2(1+I)^2/(1-H)^2\) | wl 62-63 (extra) | md:151-156, §3 (implicit in \(\Sigma_r\) form) | MATCH (route-bridge; the M/I/H form is the notes' \(\Sigma_r\) definition) |

INTERNAL scaffolding (accounted-for, no prose expected): `Mcal/Ical/Hcal/Lambda` and their perturbed `_p` forms (intermediate invariants), `Kp/GWp/.../OW2p` (perturbation builders), `P0_exact/DA/NA/PA/P_slope_series` (prefactor intermediates), `Xi_load=n1-d1`, `eps/lam/epsilon` symbols, residual `= 0` print values and `PASS`/`FAIL` flags.

reconciliation: complete; 12 deliverable values checked, 0 misaligned.

## Self-test notes

Traps checked: (1) Variable independence — each `series`/`D[...,eps]` acts on an expression in which `eps` genuinely appears (every primitive is multiplied by `Exp[eps*lam*·]`), so no slope is identically zero by construction; the assertions are `expr - lam·closedform == 0`, not `nonzero` checks, so a vanishing-derivative trap cannot produce a silent pass. (2) Domain — `K,OU2,OW2,R,GU,GW` positivity is required (for `Sqrt[K]`, `Log`, and `1-H` sign) and is consistent with the physical setup; perturbation slopes are `real` only, so no spurious positivity hides a branch. (3) Trivial-case — substituting all slopes \(=0\) collapses every residual to \(0\) and the M/I/H↔\(\Lambda\) bridge holds identically, consistent with the PASS transcript. No finding; no directive written.
