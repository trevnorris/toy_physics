---
unit_id: 180
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage180_effective_transfer_shape_collapse.md]
  paper_appendix: present
---

# Audit unit 180 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_180.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage180_effective_transfer_shape_collapse.md`
- part appendix: `/var/projects/toy_physics/.../paper/appendices/stage_appendix_part05.tex` (rows at lines 91, 671-713, 715 reviewed)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.txt`

## What the paper claims

Stage 180 ("Effective transfer-shape collapse") collapses the Stage-179 weighted sum of port transfer-shape slopes to the logarithmic slope of a **single** effective transfer shape. `\stagefield{Output}`: "Collapses many ports to \(\mathcal T_{{\rm eff},A}^2=N_{A,0}/K_A\) and gives the actual one-port continuum formula." The notes enumerate the distinct deliverables: (1) the multi-port collapse identity \(\Xi_1=2\sum_r\rho_r^{(N)}\tau_r=\delta\ln\mathcal T_{\mathrm{eff},A}^2/(\epsilon\lambda_A)\) with \(\mathcal T_{\mathrm{eff},A}^2=N_{A,0}/K_A\); (2) the one-port continuum shape \(\mathcal T_A^2=\beta_{0,A}/K_{0,A}=\mu_W/K_W\cdot Z_W(1+\rho)^2/(1-\epsilon_W)^2=Z_W(1+\rho)^2/[\Omega_W^2(1-\epsilon_W)^2]\); (3) the selected-branch form \(\mathcal T_A^2=(27\pi^2 G c_s^5/20a^5c^5)(1-\epsilon_\eta)/R_{\mathrm{target}}\); (4) the direct slope law \(\Xi_1=\zeta_W-\omega_W+2\rho_1/(1+\rho)+2\varepsilon_W/(1-\epsilon_W)\); and (5) the selected-branch slope law \(\Xi_1=-\eta_1/(1-\epsilon_\eta)-\mathcal R_1\). The stage is an exact-closure card; both engines required.

## What the script claims to verify

The SymPy docstring enumerates four checks; in code they map to six `expect_zero` assertions: (A1) multi-port collapse \(\Xi_{\rm eff}-2\sum\rho_r\tau_r=0\) via a two-port exponential model; (A2) \(T^2=\beta_0/K_0\) reduces to \(\mu_W/K_W\) form; (A3) that same form equals \(Z_W(1+\rho)^2/[\Omega_W^2(1-\epsilon_W)^2]\) under \(\Omega_W^2=K_W/\mu_W\); (A4) the selected-branch identity (with \(R_{\rm target}\) substituted by its definition on both sides); (A5) the direct weak-axisymmetric slope law via a log-derivative of a perturbed \(T^2\); (A6) the selected-branch slope law via the same log-derivative technique. The `.wl` mirrors all six. These are exactly the five notes deliverables (A2 and A3 together cover deliverable 2).

## Paper ↔ script cross-check

| paper deliverable | script check | status |
|---|---|---|
| (1) multi-port collapse \(\Xi_1=2\sum\rho_r\tau_r=\delta\ln\mathcal T_{\rm eff}^2/\epsilon\lambda\) | A1 (py 41-48 / wl 30-36) | match |
| (2) one-port continuum \(\mathcal T_A^2=\beta_0/K_0=Z_W(1+\rho)^2/[\Omega_W^2(1-\epsilon_W)^2]\) | A2+A3 (py 64-72 / wl 45-51) | match |
| (3) selected-branch \(\mathcal T_A^2=(27\pi^2Gc_s^5/20a^5c^5)(1-\epsilon_\eta)/R_{\rm target}\) | A4 (py 82-89 / wl 58-64) | match |
| (4) direct slope law | A5 (py 100-109 / wl 71-76) | match |
| (5) selected-branch slope law | A6 (py 111-119 / wl 78-81) | match |

Dominant pattern: aligned. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero(Xi_eff - Xi_expected)` | claim 1 | yes |
| A2 | sympy | 66 | `expect_zero(T2_direct - T2_expected)` | claim 2 (first form) | yes |
| A3 | sympy | 68-72 | `expect_zero(T2_direct - ZW(1+rho)^2/[OmegaW2(1-epsW)^2] /. OmegaW2->KW/muW)` | claim 2 (Omega form) | yes |
| A4 | sympy | 86-89 | `expect_zero(T2_direct{Rtgt}-T2_selected{Rtgt})` | claim 3 | yes |
| A5 | sympy | 109 | `expect_zero(Xi_direct - Xi_direct_expected)` | claim 4 | yes |
| A6 | sympy | 119 | `expect_zero(Xi_sel - Xi_sel_expected)` | claim 5 | yes |
| M1-M6 | mathematica | 36,47,48,61,76,81 | `expectZero[...]` mirroring A1-A6 | claims 1-5 | yes |

All six SymPy rows are non-tautological: A1 and A5/A6 build `Xi` from a `D[Log[...]]` derivative of an *independently-perturbed* `T^2`/`Teff2` model and compare to a *separately written-out* closed form — the two sides are constructed by different routes (calculus vs. hand-formula), so a wrong slope law would surface. A2/A3/A4 are genuine algebraic-equality reductions of distinct symbolic forms, not `x==x`.

Note on A3 (py 68-72): Python attribute-access precedence binds `.subs(...)` to the immediately-preceding parenthesized primary `(OmegaW2 * (1 - epsW)**2)`, i.e. the whole denominator group, NOT just `(1-epsW)**2`. So the substitution `OmegaW2 -> KW/muW` is correctly applied to the whole denominator and the residual is genuinely zero. (I initially suspected a precedence no-op bug here; on careful re-reading of the parenthesization it is correct, and the saved `=0` output is consistent with that.) No finding.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl:30-81`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py:41-119`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Every object is the same symbol-renamed construct in the same order:
- multi-port: `.py:41` `Teff2 = T1**2*sp.exp(2*eps*lam*tau1)+T2**2*sp.exp(2*eps*lam*tau2)` ↔ `.wl:30` `teff2 = t1^2*Exp[2*eps*lam*tau1]+t2^2*Exp[2*eps*lam*tau2]`; both then `D[Log[teff2],eps]/.eps->0` over `lam` (`.py:42` / `.wl:31`) and compare to `2*(rho1*tau1+rho2*tau2)` with identical `rho1=T1^2/(T1^2+T2^2)` (`.py:44-46` / `.wl:33-35`).
- one-port: `.py:61-65` `K0=Keta/mueta; beta0=(muW/mueta)*(Keta/KW)*ZW*(1+rho)^2/(1-epsW)^2; T2_direct=beta0/K0` ↔ `.wl:43-45` byte-for-byte the same construction.
- slope laws: `.py:100-104` `T2_pert = ZW*sp.exp(e*lam*zetaW)*(1+rho+e*lam*rho1s)^2/((KW/muW)*sp.exp(e*lam*omegaW)*(1-epsW-e*lam*epsW1)^2)` ↔ `.wl:71-73` identical; both differentiate `Log[t2Pert]` at `e->0`.

Same variable choreography, same intermediate `T2_direct`/`t2Direct`, same `.subs`/`/.` substitution targets, same six checks in the same order, identical trailing Print block. This is echo, not a second derivation from the physical premises.

**Why this matters:**
The two-engine policy requires the second engine to re-derive the result independently so that a systematic error in one engine's algebra (or in the author's hand-written closed form) cannot be silently reproduced by the other. A transliteration provides no such independent cross-check; it only confirms that two CASes agree the author copied the same expressions. Memory flags transliteration as historically UNDER-called.

**Required change:**
Re-author the `.wl` so at least the load-bearing identities are reached by a route the `.py` does not use. Concretely, for the one-port/selected-branch identities (A2-A4) the Mathematica side can build `\mathcal T_A^2` from `\beta_{0,A}` and `K_{0,A}` symbolically and then `Solve`/`Eliminate` against `\Omega_W^2 = K_W/\mu_W` and the `R_target` definition rather than `/.`-substituting the same hand-form; for the slope laws (A1,A5,A6) derive the closed-form slope by `Series[Log[...], {e,0,1}]` coefficient extraction independently and compare to the perturbed model, rather than mirroring the SymPy `D[Log[...]]/.e->0` step verbatim. Keep the assertions; change the derivation path on the Mathematica side. (Re-author-vs-accept is a USER-LEVEL call per process memory; the directive routes this to user/Codex as a re-author request.)

**Verification:**
The verifier confirms the `.wl` no longer mirrors the `.py` line-for-line (different construction of `t2Direct`, `xiDirect`, `xiSel`) and still exits 0 with all six `PASS` lines, AND the refreshed `.wl` output banner reads STAGE 180.

### F2 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.txt:3`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.txt:3`

**What's wrong:**
Both committed saved outputs are stale. (a) mtime: both `.txt` files are `May 30 01:39`; both scripts are `Jun 3 15:59` (the Phase-1 banner relabel commit e2a4780 changed only the banner string `STAGE 163` → `STAGE 180`). (b) content: both `.txt` line 3 still print the OLD banner `STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE`, whereas the current scripts print `STAGE 180` (`.py:33`, `.wl:26`). The numeric checks all still read `= 0` / `PASS`, so the math result is unaffected — only the banner label is stale.

**Why this matters:**
The committed transcript no longer matches what the current script emits (wrong stage banner). Low severity because every assertion still passes and no value changed; this is the known SCRIPT/OUTPUT-BAND staleness (per project memory, committed `.txt` outputs predate the banner fix). It is informational and will clear on the orchestrator's independent re-run, which regenerates the `.txt` with the STAGE 180 banner.

**Required change:**
No source edit needed beyond the F1 re-author. The orchestrator's standard re-run of `python3 <py>` and `math -script <wl>` regenerates both `.txt` files; confirm the regenerated line 3 reads `STAGE 180 — EFFECTIVE TRANSFER-SHAPE COLLAPSE` in both.

**Verification:**
After re-run, both `.txt` line 3 reads STAGE 180; mtimes are newer than the scripts.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration — see F1. Three corresponding sections:
1. `.py:41-42` `Teff2 = T1**2*sp.exp(2*eps*lam*tau1)+T2**2*sp.exp(2*eps*lam*tau2)` / `Xi_eff = sp.diff(sp.log(Teff2),eps).subs(eps,0)/lam` ↔ `.wl:30-31` `teff2 = t1^2*Exp[2*eps*lam*tau1]+t2^2*Exp[2*eps*lam*tau2]` / `xiEff = (D[Log[teff2],eps]/.eps->0)/lam`.
2. `.py:61-65` `K0=Keta/mueta; beta0=(muW/mueta)*(Keta/KW)*ZW*(1+rho)^2/(1-epsW)^2; T2_direct=beta0/K0` ↔ `.wl:43-45` `k0=kEta/muEta; beta0=(muW/muEta)*(kEta/kW)*zW*(1+rho)^2/(1-epsW)^2; t2Direct=beta0/k0`.
3. `.py:100-105` perturbed `T2_pert` and `Xi_direct=sp.diff(sp.log(T2_pert),e).subs(e,0)/lam` ↔ `.wl:71-74` perturbed `t2Pert` and `xiDirect=(D[Log[t2Pert],e]/.e->0)/lam`.
Identical construction order, identical perturbation ansatz, identical differentiation step. Verdict: transliteration (F1).

## Engine cross-check

Both engines verify the same six identities and both report all-zero residuals (sympy `.txt` lines 5,10,11,16,21,22 = `0`; mathematica `.txt` lines 5,11,13,19,25,27 = `0` with `PASS`). They agree. (The agreement is expected and not informative as an independent cross-check precisely because of F1.) No `engine_disagreement`.

## Verdict justification

The math is sound and aligned with the paper: all five notes/appendix deliverables map 1:1 to non-tautological SymPy assertions, the slope laws are built by an independent log-derivative route distinct from their hand-written closed forms, the constant `27\pi^2 G c_s^5/20a^5c^5` and every symbolic form match the `.tex` (appendix line 710) and notes verbatim, and the suspected A3 precedence bug does not exist on careful reading. The verdict is `findings` for two reasons: (F1) the Mathematica script is a line-by-line transliteration of the SymPy script and therefore not the independent second engine the policy requires; (F2) both saved outputs are stale (old `STAGE 163` banner, pre-relabel mtime) — informational, clears on re-run. No paper_misalignment, no tautology, no value mismatch.

## Value Reconciliation (pass-2 augmentation)

All deliverables of stage 180 are **symbolic closed forms**; there are no pinned numeric figures-of-merit (no `Pi_star`/`gamma_0`-style constants). The only literal is the dimensionless coefficient `27/(20)` and `\pi^2` in the front factor of the selected-branch shape, which is a structural part of a symbolic deliverable, not a standalone computed number.

| value (symbolic deliverable) | source (py / wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `T_eff^2 = sum_r T_r^2 = N_0/K` | py:48,122 / wl:36,85; txt line 25/31 | tex `\stagefield{Output}` L15; appendix L690-694; notes L37,76-81 | MATCH |
| `Xi_1 = delta ln(T_eff^2)/(eps lambda)` = `2 sum rho_r tau_r` | py:42-48 / wl:31-36; txt L5 | appendix L699-701; notes L104-109 | MATCH |
| `T^2 = beta0/K0 = muW/KW · ZW(1+rho)^2/(1-epsW)^2` | py:64-66 / wl:45-47; txt L11 | notes L154-161 | MATCH |
| `T^2 = ZW(1+rho)^2 / [OmegaW^2 (1-epsW)^2]` | py:68-72 / wl:48-51; txt L13/L13 | appendix L708; notes L169-173,400 | MATCH |
| `T^2 = (27 pi^2 G c_s^5 / 20 a^5 c^5)·(1-eps_eta)/R_target` | py:85-89 / wl:60-64; txt L16/L19 | appendix L710-711; notes L272-274,402-403 | MATCH |
| direct slope `Xi_1 = zeta_W - omega_W + 2 rho_1/(1+rho) + 2 epsW_1/(1-eps_W)` | py:106-109 / wl:75-76; txt L21 | notes L216-221 | MATCH |
| selected-branch slope `Xi_1 = -eta_1/(1-eps_eta) - R_1` | py:118-119 / wl:80-81; txt L22 | notes L287-292 | MATCH |

INTERNAL (scaffolding, no finding, not expected in prose): two-port model symbols `T1,T2,tau1,tau2,eps,lam`; `rho1,rho2` two-port weights; the perturbation tag symbols `e`, `zetaW,omegaW,rho1s,epsW1,eta1,R1`; `Lambda`/`lambdaNorm` and `Rtarget_def` intermediate; pass/fail `= 0` residual flags.

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

## Self-test notes

I checked: (1) variable-independence of the two `D[Log[...]]` slope derivatives — `T2_pert` and `T2_sel_pert` both genuinely depend on `e`, so the derivatives are non-trivially nonzero before the `/lam` and the slope laws are real exercises, not identically-zero traps. (2) The suspected SymPy A3 precedence no-op: re-read the parenthesization, `.subs` binds to the full `(OmegaW2*(1-epsW)**2)` primary, so the substitution lands correctly and the `=0` is genuine — no finding. (3) Constant round-trip: the `27 pi^2/(20)` front factor matches appendix L710 and notes L272/L402 exactly, so the stale-output regeneration introduces no new paper_misalignment. F1 is a re-author (USER-level) routing, F2 is informational; neither is auto-applicable by Codex as a silent script rewrite without user direction on the re-author.
