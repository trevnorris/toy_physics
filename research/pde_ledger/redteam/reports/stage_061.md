---
unit_id: 061
batch: III.3
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-26T17:39:56Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage061_microscopic_gain_thresholds.md
  paper_appendix: present
---

# Audit unit 061 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_061.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage061_microscopic_gain_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 100; stage `\input` at line 240)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage061_microscopic_gain_thresholds_mathematica_audit.txt`

## What the paper claims

Paper card `\stagefield{Output}` states: "The microscopic gain phase diagram \eqref{eq:app-stage061-G-thresholds}--\eqref{eq:app-stage061-phase-diagram}." This bundles two boxed equations: (i) the threshold pair `G_fail = Xi_fail/kappa`, `G_suff = Xi_suff/kappa`, built atop the identity `Xi_micro = kappa G_micro`; (ii) the phase-diagram region statement `G_micro < G_fail => fail`, `G_micro > G_suff => succeed`, `G_fail <= G_micro <= G_suff => root-sensitive`. The supporting notes enumerate six deliverables: (1) `Xi_micro = kappa G_micro` with `G_micro := chi_sigma Lambda_phi^2/K_X`; (2) threshold map in terms of Stage-41's `Delta_0, Delta_inf`; (3) microscopic threshold surfaces for `chi_sigma` and `Lambda_phi^2`; (4) soft-support limits `Delta_0 -> 1/2`, `Delta_inf -> 1`, with `G_fail ~ Pe_req/kappa`, `G_suff ~ 2 Pe_req/kappa`; (5) compliant-mouth limit (`eta -> +infty`) with `Delta_0 -> (1-sech(sqrt(kappa)))/kappa`, `Delta_inf -> tanh(sqrt(kappa))/sqrt(kappa)`, hence `G_fail^(inf) = Pe_req/[sqrt(kappa) tanh(sqrt(kappa))]`, `G_suff^(inf) = Pe_req/[1-sech(sqrt(kappa))]`; (6) combined stiff-support/compliant-mouth asymptotics `G_fail ~ Pe_req/sqrt(kappa)`, `G_suff ~ Pe_req`. The appendix row (part03 line 100) summarises this as "Dimensionless gain phase diagram and geometry-controlled success/failure surfaces."

## What the script claims to verify

Starting from Stage-41's endpoint functions `Delta_0(kappa, eta)`, `Delta_inf(kappa, eta)` and `alpha = sqrt(kappa)`, the SymPy script asserts: (A1) `Xi_micro - kappa G_micro = 0` after `T_X -> K_X L^2/kappa`; (A2-A5) `chi_fail`, `chi_suff`, `Lambda^2_fail`, `Lambda^2_suff` agree with `(K_X/Lam^2) G_fail`, `(K_X/Lam^2) G_suff`, `(K_X/chi) G_fail`, `(K_X/chi) G_suff` after the same substitution; (A6-A7) `Delta_0|_{kappa->0+} = 1/2`, `Delta_inf|_{kappa->0+} = 1`; (A8-A9) `lim_{kappa->0+}(kappa G_fail) = Pe_req`, `lim_{kappa->0+}(kappa G_suff) = 2 Pe_req`; (A10-A11) compliant-mouth endpoint formulas; (A12-A13) compliant-mouth gain formulas; (A14-A15) stiff-support compliant-mouth asymptotics `sqrt(kappa) G_fail -> Pe_req`, `G_suff -> Pe_req`. The Mathematica script (B1-B15) mirrors each check in its own algebra. Banner text says "STAGE 44 / 044" - a vestige of the pre-renumbering numbering (this unit was Stage 44 prior to renumbering to 061); cosmetic only, file paths and assertions are correctly anchored to unit 061.

## Paper <-> script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Notes (1): `Xi_micro = kappa G_micro`, `G_micro := chi_sigma Lambda_phi^2/K_X` | A1 / B1 | match |
| Paper eq (G-thresholds) / Notes (2): `G_fail = Pe_req/[kappa Delta_inf]`, `G_suff = Pe_req/[kappa Delta_0]` | Definition at sympy L46-47 / wl L39-40; exercised downstream via A8-A15 | match (definitional, then exercised by all limit checks) |
| Paper eq (phase-diagram): region statements over the thresholds | No separate check (region definition follows directly from the threshold values) | match (no separate assertion needed) |
| Notes (3): `chi_sigma`, `Lambda_phi^2` threshold surfaces | A2-A5 / B2-B5 | match |
| Notes (4): soft-support endpoint values and leading scalings | A6-A9 / B6-B9 | match |
| Notes (5): compliant-mouth endpoint and threshold formulas | A10-A13 / B10-B13 | match |
| Notes (6): combined stiff-support/compliant-mouth asymptotics | A14-A15 / B14-B15 | match |

Every paper-side deliverable maps onto a non-tautological script-side assertion. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40-41 | `expect_zero(Xi_micro.subs(TX, KX*L^2/kappa) - kappa*G_def)` | Notes (1) / paper eq Xi-micro | yes |
| A2 | sympy | 59-60 | `expect_zero(chi_fail.subs(...) - (KX/Lam^2)*G_fail)` | Notes (3) | yes |
| A3 | sympy | 61-62 | `expect_zero(chi_suff.subs(...) - (KX/Lam^2)*G_suff)` | Notes (3) | yes |
| A4 | sympy | 63-64 | `expect_zero(Lam2_fail.subs(...) - (KX/chi)*G_fail)` | Notes (3) | yes |
| A5 | sympy | 65-66 | `expect_zero(Lam2_suff.subs(...) - (KX/chi)*G_suff)` | Notes (3) | yes |
| A6 | sympy | 73 | `expect_zero(Delta0_k0 - 1/2)` | Notes (4) | yes |
| A7 | sympy | 74 | `expect_zero(Delta_inf_k0 - 1)` | Notes (4) | yes |
| A8 | sympy | 76-77 | `expect_zero(limit(kappa*G_fail,kappa,0+) - Pe_req)` | Notes (4) | yes |
| A9 | sympy | 78-79 | `expect_zero(limit(kappa*G_suff,kappa,0+) - 2 Pe_req)` | Notes (4) | yes |
| A10 | sympy | 86-87 | `expect_zero(Delta0_eta_inf - (1-sech alpha)/kappa)` | Notes (5) | yes |
| A11 | sympy | 88-89 | `expect_zero(Delta_inf_eta_inf - tanh(alpha)/alpha)` | Notes (5) | yes |
| A12 | sympy | 95-96 | `expect_zero(G_fail_inf - Pe_req/(alpha tanh alpha))` | Notes (5) | yes |
| A13 | sympy | 97-98 | `expect_zero(G_suff_inf - Pe_req/(1-sech alpha))` | Notes (5) | yes |
| A14 | sympy | 106 | `expect_zero(limit(z*G_fail_inf_z, z->oo) - Pe_req)` | Notes (6) | yes |
| A15 | sympy | 107 | `expect_zero(limit(G_suff_inf_z, z->oo) - Pe_req)` | Notes (6) | yes |
| B1 | math | 48 | `expectZero["Xi_micro - kappa G_micro", ...]` | Notes (1) | yes |
| B2-B5 | math | 54-57 | `chi_fail/chi_suff/Lambda^2_fail/Lambda^2_suff` cross-link | Notes (3) | yes |
| B6-B7 | math | 63-64 | `delta0K0 - 1/2`, `deltaInfK0 - 1` | Notes (4) | yes |
| B8-B9 | math | 65-66 | `kappa gFail`, `kappa gSuff` soft-support limits | Notes (4) | yes |
| B10-B11 | math | 77-78 | `delta0EtaInf`, `deltaInfEtaInf` formulas | Notes (5) | yes |
| B12-B13 | math | 79-80 | `gFailInf`, `gSuffInf` formulas | Notes (5) | yes |
| B14-B15 | math | 85-86 | stiff-support asymptotics | Notes (6) | yes |

No "no" rows; every assertion is non-tautological (a wrong factor or wrong endpoint would produce a non-zero residual).

## Findings

(None.)

## Independent-derivation check (Mathematica)

The Mathematica `.wl` mirrors the SymPy `.py` in flow: same `Delta_0`/`Delta_inf` closed forms (drawn from Stage 41 in both), same intermediate definitions of `Xi_micro`, `G_micro`, `G_fail`, `G_suff`, same threshold substitutions, same sequence of `expectZero` checks. This parallelism is structural to the task — the stage's job is to take Stage-41 endpoint functions as inputs and verify algebraic/limit identities against the operator phase diagram, so both engines start from the same inputs. However, each engine runs its own `Limit`/`FullSimplify` independently, and the canonical forms diverge: SymPy's `G_fail^(inf)` simplifies to `Pe_req/(sqrt(kappa)*tanh(sqrt(kappa)))` while Mathematica returns `peReq*Coth[Sqrt[kappa]]/Sqrt[kappa]`; SymPy's `G_suff^(inf)` simplifies to `Pe_req*cosh(sqrt(kappa))/(cosh(sqrt(kappa)) - 1)` while Mathematica returns `peReq*(1 + 1/(-1+Cosh[Sqrt[kappa]]))`. The stiff-support representation differs too: SymPy uses `subs(alpha, z)`, Mathematica uses `kappa -> z^2`. The simplification engines reach equivalent results through different canonicalisations. This is acceptable parallel verification, not transliteration.

## Engine cross-check

Both engines exit 0 with PASS on every assertion. Side-by-side selected forms:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `Delta_0` | `eta*(cosh(sqrt(kappa))-1)/(kappa*(eta*cosh(sqrt(kappa))+sqrt(kappa)*sinh(sqrt(kappa))))` | `eta*(-1+Cosh[Sqrt[kappa]])/(eta*kappa*Cosh[Sqrt[kappa]]+kappa^(3/2)*Sinh[Sqrt[kappa]])` |
| `Delta_inf` | `(eta*sinh(sqrt(kappa))+sqrt(kappa)*(cosh(sqrt(kappa))-1))/(sqrt(kappa)*(eta*cosh(sqrt(kappa))+sqrt(kappa)*sinh(sqrt(kappa))))` | `(-1+Cosh[Sqrt[kappa]]+(eta*Sinh[Sqrt[kappa]])/Sqrt[kappa])/(eta*Cosh[Sqrt[kappa]]+Sqrt[kappa]*Sinh[Sqrt[kappa]])` |
| `Delta_0(kappa->0+)` | `1/2` | `1/2` |
| `Delta_inf(kappa->0+)` | `1` | `1` |
| `Delta_0(eta->oo)` | `(1 - 1/cosh(sqrt(kappa)))/kappa` | `(1 - Sech[Sqrt[kappa]])/kappa` |
| `Delta_inf(eta->oo)` | `tanh(sqrt(kappa))/sqrt(kappa)` | `Tanh[Sqrt[kappa]]/Sqrt[kappa]` |
| `G_fail^(inf)` | `Pe_req/(sqrt(kappa)*tanh(sqrt(kappa)))` | `peReq*Coth[Sqrt[kappa]]/Sqrt[kappa]` |
| `G_suff^(inf)` | `Pe_req*cosh(sqrt(kappa))/(cosh(sqrt(kappa))-1)` | `peReq*(1 + 1/(-1+Cosh[Sqrt[kappa]]))` |

All residuals against the target formulas are 0 in both engines. Mathematica's `Limit::alimv` warnings are benign (it ignores assumptions involving the limit variable; direction is set explicitly).

## Verdict justification

Verdict: `clean`. Attacks attempted and rebuffed: (a) tautology probe on A1-A5 - these substitute `T_X -> K_X L^2/kappa` (the kappa definition from Stage 40) and check equality; a wrong factor in `G_micro` or `chi_fail`/etc. would leave a non-zero residual. (b) Sign/factor probe on soft-support limit: leading-order expansion `cosh(sqrt(k))-1 ~ k/2`, `sinh(sqrt(k)) ~ sqrt(k)` gives `Delta_0 -> 1/2`, `Delta_inf -> 1`; matches A6-A7. (c) eta->infty limit: ratio of leading terms reproduces `(1-sech alpha)/kappa` and `tanh(alpha)/alpha`; matches A10-A11. (d) `G_fail^(inf)` consistency: `Pe_req/(kappa * tanh(alpha)/alpha) = Pe_req/(alpha tanh alpha)` since `kappa = alpha^2`; matches A12. (e) Stiff-support asymptotics: `z*Pe_req/(z tanh z) -> Pe_req` and `Pe_req/(1-sech z) -> Pe_req` as z -> oo; matches A14-A15. (f) Branch probe: `Delta_inf` and `Delta_0` are strictly positive for `kappa>0, eta>0`, so no division-by-zero risk. (g) Symbol assumption probe: positivity on `kappa, eta, Pe_req, K_X, T_X, L, chi_sigma, Lambda_phi` is consistent with the physical setup in both engines. (h) Output freshness: sympy output mtime 1778525054 > script mtime 1775068798; mathematica output mtime 1778525675 > script mtime 1778522212; both fresh. (i) Paper alignment: every notes-listed deliverable maps onto a script assertion, and the paper card's boxed Output is fully covered. The `\stagefield{Output}` references the phase-diagram region statement, which is a logical consequence of the threshold definitions and requires no separate equality check. No findings.

## Self-test notes

I checked: (1) variable-dependence trap - no `diff` calls; the `subs(TX -> KX L^2/kappa)` substitution involves `TX` which actually appears in `Xi_micro`/`chi_fail`/`chi_suff`/etc., so the substitution is non-trivial; the `subs(alpha, z)` in A14 puts `z` where `alpha = sqrt(kappa)` actually appears in `G_fail_inf`/`G_suff_inf`. (2) Trivial-case pre-check for A1: set `chi=Lam=KX=L=1`, then `Xi_micro = 1/TX`; sub `TX = 1/kappa` gives `Xi_micro = kappa`, and `kappa*G_def = kappa*1*1/1 = kappa`. Residual 0, as expected. For A6: `eta(cosh(sqrt(k))-1)/(k(sqrt(k) sinh(sqrt(k))+eta cosh(sqrt(k))))` at `k=0+` Taylor-expands to `eta*(k/2)/(k*eta) = 1/2`. Correct. (3) Symmetry/parity - no improper integrals; trap N/A. (4) Path specifications - all script paths verified to exist at `scripts/` and `mathematica/`. (5) Paper round-trip - notes file is fully covered by the script's 15 assertions; no extra deliverable is left unverified, no script assertion lacks paper-side anchor. No traps fired.
