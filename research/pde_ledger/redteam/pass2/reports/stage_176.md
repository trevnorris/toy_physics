---
unit_id: 176
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
  notes_stage_files: [moving_throat_pde_stage176_outgoing_load_factorization.md]
  paper_appendix: present
---

# Audit unit 176 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_176.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage176_outgoing_load_factorization.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows 83, 611-639)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage176_outgoing_load_factorization_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage176_outgoing_load_factorization_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: *"Factors the outgoing load into mixed-leg, interference, and hybridization slippages; under rigidity, zero defect requires \(G_W/\Omega_W^2\propto\sqrt K\)."* The notes (the authoritative carrier here, since the .tex card is terse) spell out four deliverables: (1) the **exact factorization** \(\Lambda_r^2/K=\mathcal M_r^2(1+\mathcal I_r)^2/(1-\mathcal H_r)^2\) with \(\mathcal M_r=G_{W,r}/(\Omega_{W,r}^2\sqrt K)\), \(\mathcal I_r=R_rG_{U,r}/(\Omega_{U,r}^2 G_{W,r})\), \(\mathcal H_r=R_r^2/(\Omega_{U,r}^2\Omega_{W,r}^2)\) (notes §2, lines 54-69); (2) the **factored first-order defect** \(\Sigma_r^{(N)}=2\delta\ln\mathcal M_r+\frac{2\mathcal I_r}{1+\mathcal I_r}\delta\ln\mathcal I_r+\frac{2\mathcal H_r}{1-\mathcal H_r}\delta\ln\mathcal H_r\) (§4, lines 130-141); (3) the **expanded primitive-variable transport** law (§5, lines 154-171); (4) the **rigidity corollary** that under \(\delta\ln\mathcal I_r=\delta\ln\mathcal H_r=0\), \(\Sigma_r^{(N)}=2\delta\ln\mathcal M_r\), giving the square-root mixed-leg law (§6-7, lines 176-236). \(\Lambda_r=P_r/\Delta_r\), \(P_r=\Omega_{U,r}^2G_{W,r}+R_rG_{U,r}\), \(\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2\) (§1). The appendix (lines 618-639) restates the same factorization, \(\Sigma^{(N)}\) decomposition, and \(G_W/\Omega_W^2\propto\sqrt K\) law.

## What the script claims to verify

The docstring lists exactly four checks: (1) exact factorization of \(\Lambda^2/K\); (2) the factored first-order defect decomposition; (3) the expanded primitive-variable transport; (4) the rigidity corollary reducing the defect to \(2\delta\ln\mathcal M\). The SymPy script builds \(\Lambda\), \(\mathcal M\), \(\mathcal I\), \(\mathcal H\) from the primitive port symbols and asserts the algebraic factorization is exact (`expect_zero`, line 41-44); it then perturbs every primitive by `exp(eps·d·)`, extracts the first-order log-drift `Sigma_exact` via a Taylor series in `eps`, and checks it equals both the factored form (line 82) and the expanded primitive-variable form (line 92); finally it substitutes the rigidity relations and checks `Sigma_exact` reduces to `2·dlnM` (line 112). The Mathematica script mirrors the same four claims with an independent first-order extraction route.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Exact factorization \(\Lambda^2/K=\mathcal M^2(1+\mathcal I)^2/(1-\mathcal H)^2\) | py:41-44 / wl:37-40 `expect_zero(Lambda^2/K - M^2(1+I)^2/(1-H)^2)` | match |
| (2) Factored first-order defect \(\Sigma^{(N)}\) | py:79-82 / wl:72-74 `Sigma_exact - Sigma_factored == 0` | match |
| (3) Expanded primitive-variable transport | py:84-92 / wl:76-85 `Sigma_exact - Sigma_expanded == 0` | match |
| (4) Rigidity corollary → \(2\delta\ln\mathcal M\) / square-root law | py:106-112 / wl:97-100 `Sigma_exact_rigid - 2·dlnM_rigid == 0` | match |
| Micro drift formulas \(\delta\ln\mathcal M,\mathcal I,\mathcal H\) | py:75-77 / wl:68-70 (used inside checks 2 & 4) | match |

`paper_alignment: aligned` — every notes/appendix deliverable has a faithful, non-tautological script counterpart in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-44 | `simplify(Lambda^2/K - M^2(1+I)^2/(1-H)^2) == 0` | claim 1 (factorization) | yes |
| A2 | sympy | 79-82 | `simplify(Sigma_exact - Sigma_factored) == 0` | claim 2 (factored defect) | yes |
| A3 | sympy | 84-92 | `simplify(Sigma_exact - Sigma_expanded) == 0` | claim 3 (expanded transport) | yes |
| A4 | sympy | 106-112 | `simplify(Sigma_exact_rigid - 2·dlnM_rigid) == 0` | claim 4 (rigidity corollary) | yes |
| B1 | mathematica | 37-40 | `expectZero[lambda^2/k - mCal^2(1+iCal)^2/(1-hCal)^2]` | claim 1 | yes |
| B2 | mathematica | 72-74 | `expectZero[sigmaExact - sigmaFactored]` | claim 2 | yes |
| B3 | mathematica | 76-85 | `expectZero[sigmaExact - sigmaExpanded]` | claim 3 | yes |
| B4 | mathematica | 97-100 | `expectZero[sigmaExactRigid - dlnMRigid]` | claim 4 | yes |

All eight assertions trace to a specific paper deliverable. None is tautological: `Sigma_exact` (the load-bearing object in A2/A3/A4 and B2/B3/B4) is derived *from the perturbed \(\Lambda\)*, not from the candidate `Sigma_factored`/`Sigma_expanded`/`2·dlnM` forms — so each check is `independent_derivation − stated_form == 0`, which genuinely tests the stated form.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent route**, not a transliteration. Three points:

1. **Different first-order extraction method.** SymPy reads the \(\varepsilon^1\) term off a truncated Taylor series:
   `sp.series(sp.log((Lambdap**2/Kp)/(Lambda**2/K)), eps, 0, 2).removeO()/eps` (py:62-72).
   Mathematica instead takes a symbolic derivative and evaluates at the origin:
   `D[Log[((lambdaP^2/kP)/(lambda^2/k))], eps] /. eps -> 0` (wl:62-65).
   These are mathematically distinct procedures for obtaining the linear coefficient; the `.wl` itself documents the divergence in its comment (wl:60-61): *"Independent first-order extraction: Mathematica uses symbolic D[Log,eps] at eps->0; the SymPy audit instead Taylor-expands and reads the eps^1 term."*
2. **Objects built from premises, not ported.** `lambda`, `mCal`, `iCal`, `hCal` are each reconstructed in Mathematica from the same primitive symbols (wl:32-35) rather than importing a pre-simplified SymPy form. The simplifier is `FullSimplify[Together[Expand[...]]]` (wl:20-24) vs SymPy's `simplify(expand(...))` — different engines, no shared black box.
3. The candidate forms `Sigma_factored`, `Sigma_expanded`, `2·dlnM` are necessarily the same expressions in both engines, but that is unavoidable — they *are* the paper's claimed results, and the test in both engines is `independent_object − claimed_form`. The independence lives in `Sigma_exact`, which is computed by genuinely different machinery.

No `mathematica_transliteration` finding.

## Engine cross-check

Both transcripts report identical bottom lines. SymPy output (lines 5, 13, 14, 25): `exact factorization = 0`, `factored first-order defect formula = 0`, `expanded primitive-variable transport = 0`, `rigidity reduction = 0`. Mathematica output (lines 6, 15, 17, 29) reports the same four residuals as `0` plus `PASS:` on each and `Stage 176 Mathematica audit passed.` (line 37). Engines agree. Outputs are fresh: `.py` mtime 2026-05-30 01:10:19, `.txt` 01:38:31; `.wl` mtime 01:10:19, `.txt` 01:38:40 — both transcripts post-date their scripts. No `stale_output`.

## Verdict justification

`clean`. I read the paper card, the notes, and the appendix rows before the scripts. The four notes/appendix deliverables (factorization, factored defect, expanded transport, rigidity corollary) map one-to-one onto eight substantive, non-tautological assertions across the two engines. Attacks tried and failed: (a) tautology — the load-bearing `Sigma_exact` is derived independently of the candidate forms it is checked against, so the `==0` checks can fail; (b) hidden round-trip in rigidity — A4/B4 substitute rigidity into the *independent* `Sigma_exact`, not into the factored form, so it is a real cross-check; (c) transliteration — the two engines use genuinely different first-order extraction routines (Taylor `series` vs symbolic `D`), documented in the `.wl`; (d) symbol-domain abuse — all six primitives carry `positive=True`/`>0`, which is exactly the physical setup (squared frequencies, couplings, stiffness K are positive) and is required to make `1/(OU2·OW2-R²)` and `Sqrt[K]` well-defined; no positivity is used to mask a branch error (the factorization is a pure rational identity). (e) Hand-verified the factorization algebraically: `M²(1+I)²/(1-H)²` collapses to `(OU2·GW+R·GU)²/(K(OU2·OW2-R²)²) = Lambda²/K`. Paper and script claims match exactly.

## Value Reconciliation (pass-2 augmentation)

This is a purely symbolic stage; the scripts emit no named numeric constants — every deliverable is a boxed closed-form identity. The outputs are fresh (no staleness caveat). Enumerating the symbolic deliverables and locating each in the docs:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\Lambda_r^2/K = \mathcal M_r^2(1+\mathcal I_r)^2/(1-\mathcal H_r)^2\) | py:43, wl:39; out py:5 / wl:5-6 | notes:64-69 (boxed); appendix:618-624 | MATCH |
| \(\mathcal M_r=G_{W,r}/(\Omega_{W,r}^2\sqrt K)\) | py:37, wl:33 | notes:56 | MATCH |
| \(\mathcal I_r=R_rG_{U,r}/(\Omega_{U,r}^2 G_{W,r})\) | py:38, wl:34 | notes:57 | MATCH |
| \(\mathcal H_r=R_r^2/(\Omega_{U,r}^2\Omega_{W,r}^2)\) | py:39, wl:35 | notes:58-60 | MATCH |
| \(\delta\ln\mathcal M=\delta\ln G_W-\delta\ln\Omega_W^2-\tfrac12\delta_K\) | py:75, wl:68; out py:17 / wl:21 | notes:115-120 | MATCH |
| \(\delta\ln\mathcal I=\delta\ln R+\delta\ln G_U-\delta\ln\Omega_U^2-\delta\ln G_W\) | py:76, wl:69; out py:18 / wl:22 | notes:121-124 | MATCH |
| \(\delta\ln\mathcal H=2\delta\ln R-\delta\ln\Omega_U^2-\delta\ln\Omega_W^2\) | py:77, wl:70; out py:19 / wl:23 | notes:125-129 | MATCH |
| \(\Sigma^{(N)}=2\delta\ln\mathcal M+\tfrac{2\mathcal I}{1+\mathcal I}\delta\ln\mathcal I+\tfrac{2\mathcal H}{1-\mathcal H}\delta\ln\mathcal H\) | py:79-82, wl:72-74; out py:13 / wl:14-15 | notes:130-141; appendix:626-635 | MATCH |
| Expanded primitive transport (the \(-\delta_K + \tfrac{2}{1+\mathcal I}\delta\ln G_W+\dots\) form) | py:84-92, wl:76-85; out py:14 / wl:16-17 | notes:154-171 (two boxes) | MATCH |
| Rigidity corollary \(\Sigma^{(N)}=2\delta\ln\mathcal M\) | py:106-112, wl:97-100; out py:25 / wl:28-29 | notes:191-203; appendix:636 | MATCH |
| Square-root law \(G_W/\Omega_W^2\propto\sqrt K\) | py:118 (printed conclusion), wl:107 | notes:233-235; appendix:638-639; stage_176.tex:15 (Output) | MATCH |

INTERNAL scaffolding (no finding): `Lambda`, `P_r`/`Delta_r` building blocks, the `exp(eps·d·)` perturbed copies, `Sigma_exact`, the `rigid` substitution dict, `banner`/`expect_zero`/`pass`/`fail`/`fmt` helpers, residual `=0` print lines, `PASS:` flags.

reconciliation: complete; 11 values checked, 0 misaligned.
