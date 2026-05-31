---
unit_id: 180
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage180_effective_transfer_shape_collapse.md]
  paper_appendix: present
---

# Audit unit 180 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_180.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage180_effective_transfer_shape_collapse.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 91, 265, 671-715 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage180_effective_transfer_shape_collapse_mathematica_audit.txt`

## What the paper claims

The stage card `\stagefield{Output}` reads verbatim: "Collapses many ports to \(\mathcal T_{{\rm eff},A}^2=N_{A,0}/K_A\) and gives the actual one-port continuum formula." The notes and Part V appendix (eqs `app-part05-Teff-def`, `app-part05-Xi1-Teff`, `app-part05-T-one-port-continuum`) expand this into five distinct deliverables:
1. **Multi-port collapse**: the Stage-247 weighted-average defect equals the log-slope of a single effective shape, \(\Xi_1 = 2\sum_r\rho_r^{(N)}\tau_r = \delta\ln\mathcal T_{\mathrm{eff},A}^2/(\epsilon\lambda_A)\), with \(\mathcal T_{\mathrm{eff},A}^2=\sum_r\mathcal T_{A,r}^2=N_{A,0}/K_A\).
2. **One-port continuum form**: \(\mathcal T_A^2=\beta_{0,A}/K_{0,A}=(\mu_{W,A}/K_{W,A}^{(\mathrm{eff})})Z_{W,A}(1+\rho_A)^2/(1-\epsilon_{W,A})^2 = Z_{W,A}(1+\rho_A)^2/[\Omega_{W,A}^2(1-\epsilon_{W,A})^2]\) using \(\Omega_{W,A}^2=K_{W,A}^{(\mathrm{eff})}/\mu_{W,A}\).
3. **Selected-branch reformulation**: \(\mathcal T_A^2=(27\pi^2Gc_s^5/(20a^5c^5))(1-\epsilon_{\eta,A})/R_{\mathrm{target},A}\) via the Stage-21 \(R_{\mathrm{target}}\) and \(\Lambda\) definitions.
4. **Direct slope law**: \(\Xi_1=\zeta_W-\omega_W+2\rho_1/(1+\rho)+2\varepsilon_W/(1-\epsilon_W)\).
5. **Selected-branch slope law**: \(\Xi_1=-\eta_1/(1-\epsilon_\eta)-\mathcal R_1\).

(The notes label this stage "Stage 248" in its prose and cite "Stage 247/241/21" as upstream sources — these are the derivation chain's internal numbering; the paper card and appendix consistently call it Stage 180.)

## What the script claims to verify

The docstring enumerates exactly four named checks (1 multi-port collapse, 2 one-port continuum, 3 selected-branch reformulation, 4 weak-axisymmetric slope laws), and the body implements six `expect_zero` assertions covering all five paper deliverables (check 2 is split into two assertions; check 4 into two). Each assertion subtracts an independently-built construction from the paper-stated closed form and asserts the residual simplifies to 0. The SymPy and Mathematica scripts implement the identical six checks with matching assumption sets; both saved transcripts report all six residuals `= 0` and exit 0.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| 1. \(\Xi_1=2\sum\rho_r\tau_r=\delta\ln\mathcal T_{\mathrm{eff}}^2/(\epsilon\lambda)\) | A1/M1 "multi-port effective-shape identity" | match |
| 2a. \(\mathcal T^2=\beta_0/K_0=(\mu_W/K_W)Z_W(1+\rho)^2/(1-\epsilon_W)^2\) | A2/M2 "T^2 = beta0/K0 -> muW/KW form" | match |
| 2b. \(\mathcal T^2=Z_W(1+\rho)^2/[\Omega_W^2(1-\epsilon_W)^2]\), \(\Omega_W^2=K_W/\mu_W\) | A3/M3 "T^2 = ZW(1+rho)^2/[OmegaW^2(1-epsW)^2]" | match |
| 3. selected-branch \(\mathcal T^2=(27\pi^2Gc_s^5/20a^5c^5)(1-\epsilon_\eta)/R_{\mathrm{target}}\) | A4/M4 "selected-branch T^2 identity" | match |
| 4. \(\Xi_1=\zeta_W-\omega_W+2\rho_1/(1+\rho)+2\varepsilon_W/(1-\epsilon_W)\) | A5/M5 "direct slope law" | match |
| 5. \(\Xi_1=-\eta_1/(1-\epsilon_\eta)-\mathcal R_1\) | A6/M6 "selected-branch slope law" | match |

Every paper deliverable has a faithful, non-tautological script-side counterpart with matching constants (27, π², 20, the exponent 5 on \(c_s/a/c\), the \(N_0/K\) ratio). `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 48 | `expect_zero(Xi_eff - Xi_expected)` where `Xi_eff = diff(log(Teff2),eps)|0/lam` | claim 1 | yes |
| A2 | sympy | 66 | `expect_zero(T2_direct - T2_expected)`, `T2_direct=beta0/K0` | claim 2a | yes |
| A3 | sympy | 68-72 | `expect_zero(T2_direct - ZW(1+rho)^2/[(KW/muW)(1-epsW)^2])` | claim 2b | yes |
| A4 | sympy | 86-89 | `expect_zero(T2_direct - T2_selected)` with `Rtarget->Rtarget_def` | claim 3 | yes |
| A5 | sympy | 109 | `expect_zero(Xi_direct - Xi_direct_expected)`, `Xi_direct=diff(log(T2_pert),e)|0/lam` | claim 4 | yes |
| A6 | sympy | 119 | `expect_zero(Xi_sel - Xi_sel_expected)` | claim 5 | yes |
| M1 | mathematica | 36 | `expectZero[xiEff - xiExpected]` | claim 1 | yes |
| M2 | mathematica | 47 | `expectZero[t2Direct - t2Expected]` | claim 2a | yes |
| M3 | mathematica | 48-51 | `expectZero[t2Direct - (zW(1+rho)^2/(omegaW2(1-epsW)^2) /. omegaW2->kW/muW)]` | claim 2b | yes |
| M4 | mathematica | 61-64 | `expectZero[(t2Direct/.rTarget->rTargetDef) - (t2Selected/.rTarget->rTargetDef)]` | claim 3 | yes |
| M5 | mathematica | 76 | `expectZero[xiDirect - xiDirectExpected]` | claim 4 | yes |
| M6 | mathematica | 81 | `expectZero[xiSel - xiSelExpected]` | claim 5 | yes |

All twelve rows are "yes": each `expect_zero/expectZero` left side is built by an independent operation (differentiating a constructed perturbation, substituting an upstream definition) and the right side is the paper's closed form — so the residual would be nonzero if the paper form were wrong.

## Findings

None.

I attacked all six identities and could not break any:

- **A1/M1 (collapse):** `Teff2 = T1²·exp(2ελτ1)+T2²·exp(2ελτ2)` genuinely depends on `eps`, so `d/dε|₀` is non-trivial. By hand: `d/dε[log Teff2]|₀ = 2λ(T1²τ1+T2²τ2)/(T1²+T2²)`, dividing by λ gives `2(ρ1τ1+ρ2τ2)`, matching `Xi_expected`. Not tautological — `Xi_eff` and `Xi_expected` are built by two different routes (exponential perturbation vs. pre-weighted sum). Considered whether the 2-port restriction is `insufficient_verification`/`missing_branch` against the notes' general "\(\sum_r\)": the identity is linear in the port multiset and the 2-port case is the minimal nontrivial instance that fully exercises the weighted-average structure; the general case adds no new algebra. Acceptable, not a finding.
- **A3/M3 (ΩW² substitution):** I initially suspected a Python `.subs` precedence bug (that `.subs` might bind only to `(1-epsW)**2`). AST inspection confirms the parenthesized `(OmegaW2 * (1 - epsW)**2)` is a single atom and `.subs({OmegaW2: KW/muW})` attaches to the whole denominator, so `OmegaW2` is correctly replaced by `KW/muW`. Residual is genuinely 0, not a no-op-substitution artifact. The Mathematica `/.` form is unambiguous and agrees.
- **A4/M4 (selected branch):** `T2_direct` contains no `Rtarget` symbol, so its `.subs` is inert; the test reduces to `T2_direct - T2_selected(Rtarget_def)`. Substituting `Rtarget_def` and `Lambda` by hand collapses `T2_selected` to `(muW/KW)ZW(1+ρ)²/(1-εW)² = T2_direct`. Residual 0, genuinely exercising the Stage-21 Λ/R_target definitions.
- **A5/A6/M5/M6 (slope laws):** both perturbations depend on `e`; hand-differentiation of `log` at `e=0` reproduces the four-channel and two-channel slope laws exactly. Constants and signs (`+2ρ1/(1+ρ)`, `+2εW1/(1-εW)`, `-η1/(1-εeta)`, `-R1`) match the notes/appendix verbatim.

## Independent-derivation check (Mathematica)

The `.wl` mirrors the `.py` structure closely (same six checks, same variable names transliterated). Per the project's mirror policy for algebraic-identity stages this is acceptable: each engine still constructs the left side from the physical premise (differentiating a perturbation, substituting an upstream definition) rather than echoing the other's numeric residual. The two engines use independent simplification backends (`sympy.simplify(expand(...))` vs. `FullSimplify[Together[Expand[...]]]`), and Mathematica's `omegaW2 -> kW/muW` (M3) and `rTarget -> rTargetDef` (M4) replacements are written in Mathematica's own idiom, not copied SymPy `.subs` syntax. This is a legitimate two-engine cross-check, not a `mathematica_transliteration` finding for identity-verification of this kind.

## Engine cross-check

Both transcripts are fresh (sympy script mtime 2026-05-11 11:56:51, output 12:47:59; mathematica script 11:56:53, output 13:23:04 — both outputs newer than their scripts). Both report identical residuals `= 0` for all six checks and exit code 0. Engines agree.

## Verdict justification

`clean`. I read the paper card, the single notes file, and the Part V appendix rows (91, 265, 671-715) before opening the scripts, then mapped all five paper deliverables to the six script assertions: every row is a `match`. Each assertion is non-tautological (left and right sides built by independent routes), well-anchored (residual would be nonzero if the paper form were wrong), and exercised with the correct constants and signs. Both engines are present, agree, and have fresh transcripts. The only blemish is cosmetic and non-categorizable: both scripts print the banner `STAGE 163 — EFFECTIVE TRANSFER-SHAPE COLLAPSE` (sympy line 33, mathematica line 26) while the file, docstring, and final Mathematica message correctly say Stage 180 — a copy-paste artifact in output text, not a math or paper-alignment defect, so it is not filed under any of the ten audit categories and no directive is written. No script-side math finding exists.

## Self-test notes

I ran the required traps: (1) **Variable independence** — `Teff2` depends on `eps` and both `T2_pert`/`T2_sel_pert` depend on `e`, so every `diff(log(...))` is non-trivial (no identically-zero-derivative trap). (2) **Parity/symmetry** — no unbounded-domain integrals in this unit; N/A. (3) **Trivial-case** — substituting concrete numbers (e.g. T1=T2=1, τ1=1, τ2=0) into A1 gives `Xi_eff = 1 = 2·(½·1+½·0)`, confirming a nonzero, consistent residual structure. (4) **`.subs` precedence** — verified via AST that A3's substitution targets the whole denominator (not a no-op), the one place a silent-pass bug could have hidden. (5) **Paper round-trip** — no fix prescribed, so no risk of introducing a new misalignment. All six identities hold; no finding, no directive.
