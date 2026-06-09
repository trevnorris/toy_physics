---
unit_id: 178
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
  notes_stage_files: [moving_throat_pde_stage178_outgoing_port_coloading_theorem.md]
  paper_appendix: present
---

# Audit unit 178 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_178.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage178_outgoing_port_coloading_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 87 and 654-669)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage178_outgoing_port_coloading_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage178_outgoing_port_coloading_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` says the stage "Rewrites \(\Xi_1\) as the weighted mismatch between actual static outgoing-port loading and wall-baseline loading." The notes make this concrete. The central boxed deliverable is \(\Xi_1=\sum_r\rho_r^{(N)}(\nu_r-\kappa_1)=\bar\nu_N-\kappa_1\), built from: (i) the static outgoing-transfer slope identity \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) for \(N_{A,0}^{(r)}=P_{A,r}^2/\Delta_{A,r}^2\); (ii) closed forms for the numerator/detuning slopes \(\mathfrak p_r=\alpha_r(\mathfrak o_{U,r}+\mathfrak g_{W,r})+\beta_r(\mathfrak r_r+\mathfrak g_{U,r})\) and \(\mathfrak d_r=\chi_r(\mathfrak o_{U,r}+\mathfrak o_{W,r})-2\zeta_r\mathfrak r_r\) with weights \(\alpha_r+\beta_r=1\), \(\chi_r-\zeta_r=1\); and (iii) the exact equivalence to the Stage 177 slippage formula \(\nu_r=\kappa_1+2\mathfrak m_r+\tfrac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\tfrac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\), so \(\sigma_r=\nu_r-\kappa_1\). Section 5's zero-defect / dominant-port / "naive rigidity fails" statements are corollaries that follow algebraically once \(\Xi_1=\bar\nu_N-\kappa_1\) is established; they are not separate computed values.

## What the script claims to verify

The SymPy script (docstring lines 6-13) verifies four identities and the Mathematica script mirrors them plus one extra route. (1) The \(\epsilon\lambda\)-linear log-slope of \(N_{A,0}^{(r)}=P_{A,r}^2/\Delta_{A,r}^2\) equals \(2(\mathfrak p_r-\mathfrak d_r)\) (sympy:52, wl:37). (2) The weighted sum \(\sum_r\rho_r(\nu_r-\kappa_1)\) under \(\rho_3=1-\rho_1-\rho_2\) equals \(\bar\nu_N-\kappa_1\) (sympy:69-72, wl:47-50). (3) The series-derived \(\mathfrak p_r,\mathfrak d_r\) from the actual port data \(P=\Omega_U^2G_W+RG_U\), \(\Delta=\Omega_U^2\Omega_W^2-R^2\) equal the closed-form weighted expressions, and \(\alpha+\beta-1=0\), \(\chi-\zeta-1=0\) (sympy:103-106, wl:76-79). (4) The direct slope \(2(\mathfrak p_r-\mathfrak d_r)\) equals the Stage-177 slippage form \(\kappa_1+2\mathfrak m_r+\dots\) (sympy:135, wl:101). The .wl additionally verifies that a from-scratch log-series \(\nu_r=\mathrm{Coefficient}[\mathrm{Series}[\log(P_A^2/\Delta_A^2)],\epsilon\lambda]\) equals the slippage form (wl:105-109).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) for \(N_{A,0}^{(r)}=P_{A,r}^2/\Delta_{A,r}^2\) (notes §1) | sympy:52 / wl:37 series of \(P_A^2/D_A^2\) | match |
| \(\Xi_1=\sum_r\rho_r(\nu_r-\kappa_1)=\bar\nu_N-\kappa_1\) (Output; notes §2) | sympy:69-72 / wl:47-50 (uses \(\sum\rho=1\) via \(\rho_3{=}1-\rho_1-\rho_2\)) | match |
| \(\mathfrak p_r,\mathfrak d_r\) closed forms + \(\alpha{+}\beta{=}1\), \(\chi{-}\zeta{=}1\) (notes §3) | sympy:103-106 / wl:76-79 | match |
| \(\nu_r=\kappa_1+2\mathfrak m_r+\tfrac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\tfrac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\), \(\sigma_r=\nu_r-\kappa_1\) (notes §4) | sympy:135 / wl:101 | match |
| Independent log-data route for \(\nu_r\) | wl:105-109 (no sympy counterpart) | extra (strengthens, not required) |
| §5 zero-defect / dominant / rigidity corollaries | (algebraic consequences of row 2; not separately asserted) | match (trivial corollary) |

`paper_alignment: aligned`. Every load-bearing boxed result in the notes and the Output field maps to a non-tautological script assertion.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 52 | `expect_zero(nu_from_series - 2(p_r-d_r))` | \(\nu_r=2(\mathfrak p-\mathfrak d)\) | yes |
| A2 | sympy | 69-72 | `expect_zero(Xi - (nu_bar - kappa1))` w/ \(\rho_3{=}1{-}\rho_1{-}\rho_2\) | \(\Xi_1=\bar\nu_N-\kappa_1\) | yes |
| A3 | sympy | 103 | `expect_zero(p_from_series - p_expected)` | \(\mathfrak p_r\) closed form | yes |
| A4 | sympy | 104 | `expect_zero(d_from_series - d_expected)` | \(\mathfrak d_r\) closed form | yes |
| A5 | sympy | 105-106 | `alpha+beta-1`, `chi-zeta-1` | weight constraints | yes |
| A6 | sympy | 135 | `expect_zero(nu_direct - nu_expected)` | Stage-177 equivalence | yes |
| A7 | math | 37 | `expectZero[nuFromSeries - 2(pR-dR)]` | \(\nu_r=2(\mathfrak p-\mathfrak d)\) | yes |
| A8 | math | 47-50 | `expectZero[xi - (nuBar - kappa1)]` | \(\Xi_1=\bar\nu_N-\kappa_1\) | yes |
| A9 | math | 76-79 | `p_r`, `d_r`, weight constraints | \(\mathfrak p,\mathfrak d\) + weights | yes |
| A10 | math | 101 | `expectZero[nuDirect - nuExpected]` | Stage-177 equivalence | yes |
| A11 | math | 105-109 | `expectZero[nuFromData - nuExpected]` (log-series route) | independent \(\nu_r\) derivation | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

Sections 1-3 of the `.wl` are close structural parallels of the `.py`: identical symbol setup, the same `pAr/dAr` (resp. `PAr/DAr`) construction, the same `Series[..,{eps,0,1}]` / `series(..,eps,0,2)` truncation, the same `alpha/beta/chi/zeta` definitions, and the same four `expectZero/expect_zero` residual checks. Taken alone those sections would read as a transliteration. However, section 4 contains a genuinely independent route that the SymPy script lacks (wl:103-109):

```
nuFromData = FullSimplify[
  Coefficient[Normal[Series[Log[pA^2/dA^2], {eps, 0, 1}]], eps*lam], ...];
expectZero["nu via log-data vs slippage", nuFromData - nuExpected];
```

Here \(\nu_r\) is obtained as the \(\epsilon\lambda\) coefficient of the log-series of \(P_A^2/\Delta_A^2\) built directly from the component drifts (`pA`, `dA`), bypassing the `pExpected/dExpected` weighted-average intermediates entirely, and is cross-checked against the slippage form `nuExpected`. The SymPy script reaches `nu_direct` only via `2*p_expected - 2*d_expected` (sympy:128), i.e. through the closed-form weights; it never independently log-differentiates the raw \(P^2/\Delta^2\). The output confirms this route runs and resolves to the same expression (mathematica output line 47: `nuFromData = ...`, line 46 `PASS: nu via log-data vs slippage`). This `Log[]`-based coefficient extraction is a distinct algebraic path (logarithmic-derivative vs ratio-series), so the `.wl` is NOT a pure port — it adds independent verification value. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines pass all shared checks and emit consistent closed forms. The `nu_r`, `sigma_r`, `p_r`, `d_r`, `alpha/beta/chi/zeta` printed forms agree up to ordering/sign-grouping:
- sympy `p_r = (GU*R*(gU+rr)+GW*OU2*(gW+oU))/(GU*R+GW*OU2)` ↔ math `p_r = (gw*(gW+oU)*ou2+gu*r*(gU+rr))/(gw*ou2+gu*r)` — identical.
- sympy `d_r = (OU2*OW2*(oU+oW)-2*R**2*rr)/(OU2*OW2-R**2)` ↔ math `d_r = (ou2*(oU+oW)*ow2-2*r^2*rr)/(ou2*ow2-r^2)` — identical.
- Both `Xi_1` and `nubar_N` expansions match (sympy out 22-23, math out 19-20).
- Section-1 linearized \(N_{A,0}^{(r)}\): sympy `P0r²/D0r² + 2 P0r² eps lam (p_r - d_r)/D0r²` ↔ math `p0r^2/d0r^2 - 2 eps lam p0r^2 (dR - pR)/d0r^2` — identical. `engines_agree: true`.

## Verdict justification

Clean. I read the card, the notes, and the appendix rows (lines 87, 654-669) before the scripts; the script's verified identities match the paper's Output and every boxed notes result. Attacks tried and failed: (1) Tautology — each `expect_zero` compares two independently constructed objects (a series-truncated ratio/log vs a closed-form expression), not `x==x`; A2's collapse genuinely exercises \(\sum\rho=1\) via the `rho3` substitution rather than asserting it for free; A5's `alpha+beta-1`/`chi-zeta-1` are real algebraic facts about the defined weights, not constructions. (2) Symbol-assumption — `positive=True` on `P0r,D0r,OU2,OW2,R,GU,GW` is justified (these are squared magnitudes / positive transfer factors / a coupling, matching the notes' "positive numerator weights"); positivity is needed only to let `simplify` cancel denominators \(P\neq0\), \(\Delta\neq0\), which is physical. (3) Series order — both use first-order \(\epsilon\) (the linear weak-axisymmetric slope is exactly what \(\nu_r\) is defined as); higher orders are irrelevant to a log-slope claim. (4) Misalignment — no constant or identity diverges from the card/notes; the §5 corollaries need no separate check. The `.wl` adds an independent log-data route, so the second-engine policy is satisfied. Outputs are fresh (txt mtimes 01:39 > script mtimes 01:15).

## Self-test notes

Checked: variable-independence (the `Series[...,eps]` expansions all depend on `eps` through `pAr/dAr/pA/dA`, so no identically-zero derivatives; the `eps*lam` coefficient extraction is well-posed). No unbounded integrals, so parity trap N/A. Trivial-case: setting `p_r=d_r` gives `nu_r=0` (A1 consistent); setting all drifts to 0 gives `nu_r=0` and `sigma_r=-kappa1`, matching notes §5.4 \(\Xi_1=-\kappa_1\). Paper round-trip: no fix prescribed, so no new misalignment introduced.

## Value Reconciliation (pass-2 augmentation)

All emitted deliverable values are symbolic closed forms (no numeric constants in this stage). Reconciled against the notes (the natural carrier; the `.tex` card is terse by design and the appendix carries the headline forms).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| \(\nu_r=2(\mathfrak p_r-\mathfrak d_r)\) | py:52,56 / wl:37,40; sympy out 16, math out 12 | notes:92; appendix nu-def 661-663 | MATCH |
| \(\Xi_1=\sum\rho_r(\nu_r-\kappa_1)=\bar\nu_N-\kappa_1\) | py:67-74 / wl:45-52; sympy out 22-23 | notes:134-152; Output (tex:15); appendix:667 | MATCH |
| \(\alpha_r=\Omega_U^2G_W/P\) | py:95,108; sympy out 32 | notes:186 | MATCH |
| \(\beta_r=RG_U/P\) | py:96,109; sympy out 33 | notes:188 | MATCH |
| \(\chi_r=\Omega_U^2\Omega_W^2/\Delta\) | py:97,110; sympy out 34 | notes:213 | MATCH |
| \(\zeta_r=R^2/\Delta\) | py:98,111; sympy out 35 | notes:216 | MATCH |
| \(\mathfrak p_r=\alpha_r(o_U+g_W)+\beta_r(r_r+g_U)\) | py:100,112; sympy out 36 | notes:196-200 | MATCH |
| \(\mathfrak d_r=\chi_r(o_U+o_W)-2\zeta_r r_r\) | py:101,113; sympy out 37 | notes:222-227 | MATCH |
| \(\mathcal I_r=RG_U/(\Omega_U^2G_W)\) | py:121,138; sympy out 43 | notes:271 | MATCH |
| \(\mathcal H_r=R^2/(\Omega_U^2\Omega_W^2)\) | py:122,139; sympy out 44 | notes:273 | MATCH |
| \(\nu_r=\kappa_1+2\mathfrak m_r+\tfrac{2\mathcal I_r}{1+\mathcal I_r}\mathfrak i_r+\tfrac{2\mathcal H_r}{1-\mathcal H_r}\mathfrak h_r\) | py:129-134,140; sympy out 45 | notes:288-296 | MATCH |
| \(\sigma_r=\nu_r-\kappa_1\) | py:137,141; sympy out 46 | notes:301 | MATCH |
| \(\nu_r\) via log-data route (wl-only) | wl:105-110; math out 47,50 | (consistent with notes:92,288-296) | MATCH |

INTERNAL (scaffolding, no prose expected): `eps`, `lam` (perturbation bookkeeping), `PAr/DAr/pA/dA/P/Delta` (intermediate constructions), `nu_from_series/p_from_series/d_from_series/nuFromSeries/pFromSeries/dFromSeries` (series intermediates driving assertions), `Xi`/`nu_bar` raw expansions before collapse, `expect_zero`/`expectZero` residuals and `PASS:` flags.

reconciliation: complete; 13 values checked, 0 misaligned.
