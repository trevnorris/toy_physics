---
unit_id: 084
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00-06:00
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: present
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md
  paper_appendix: present
---

# Audit unit 084 red-team report (v2)

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_084.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage084_full_reduced_pde_writeup.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (referenced via `\input{stages/stage_084}` at line 286)
- sympy: (missing — status-only carve-out applies; `is_status_only_candidate: true`, `is_checkpoint: false`)
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- sympy output: (missing)
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt`

Mathematica script mtime: 2026-05-25 00:32 (1779690721). Output mtime: 2026-05-25 00:33 (1779690790). Output is fresher than script; not stale.

Status-only carve-out: unit 084 has `is_status_only_candidate: true`, `is_checkpoint: false`. The five carried-forward numeric zeta values originate in stage 081's SymPy script; the `zeta_req`/`Q(zeta)` algebra mirrors stage 082; the `kappa_F1`, `eta_F1`, `lambdaEll` constants are defined upstream. The status-only consolidation does not need its own SymPy script.

## What the paper claims

Per `paper/stages/stage_084.tex:\stagefield{Output}`: "A structured reduced-PDE theorem chain and handoff to the outgoing-branch loading-ratio problem." The reduced PDE chain (eq:app-stage084-chain) is: continuum placement → support completion → transport/parent gain → Family-1 support map → `R_quad`. The notes elaborate four deliverables: (1) the explicit physical-zeta formula `zeta_phys(Pe, eta; kappa) = Omega_Pe^2 * (kappa + pi^2/4) / (kappa + y^2)` with `Omega_Pe = pi*Pe*(2*Pe*exp(Pe) + pi) / [(4*Pe^2 + pi^2)*(exp(Pe) - 1)]` and `y*tan(y) = eta`, (2) the demand-map `zeta_req(Pi_tr, C_mix, eps_blk) = (Pi_tr - C_mix) / [C_mix - eps_blk*(2*C_mix - Pi_tr)]` with inverse `Pi_tr = C_mix * Q(zeta; eps_blk)` and `Q = [1 + (1-2*eps_blk)*zeta] / [1 - eps_blk*zeta]`, (3) the master residual `R_quad = zeta_req - zeta_phys`, (4) the Family-1 specialization `kappa_F1 = 12321/5`, `eta_F1 = 37`, `Xi_F1 = 1369*Upsilon_w = 136900*Theta_w`, and (5) the carried-forward numeric windows `zeta_-^(chi) ≈ 2.46622291347846`, `zeta_+^(chi) ≈ 2.46752913273870`, `zeta_-^(J) ≈ 2.44257571477179`, `zeta_+^(J) ≈ 2.46752736855058`, hard ceiling `zeta_max^(F1) ≈ 2.46752922945601`. The notes also state explicitly (Section 4) that the actual passive/outgoing quadrupole branch — fixing `Pi_tr` on the true moving-throat solution — is the remaining open theorem gap; the script must therefore NOT close `R_quad` to zero.

## What the script claims to verify

The script (`moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`) is the Mathematica-only consolidation for the status-only write-up. It (a) symbolically defines `omegaPe`, `zetaPhys`, `zetaReq`, `qMap`, `rQuad` (lines 39-46), then exercises five substantive checks: (A1, line 53) the inverse demand-map round-trip `zetaReq(piTr -> cMix*qMap) - zeta == 0`; (A2, line 65) the cross-route Family-1 strength identity `Xi_F1(Upsilon_w) | _{Upsilon_w -> 100*Theta_w} - Xi_F1(Theta_w) == 0`; (A3, line 76) the `Pe -> Infinity` limit of `zetaPhys` at Family-1 parameters numerically matches the carried-forward `zetaMaxF1 = 2.46752922945601` to within `1e-10`; (A4-A7, lines 84-87) four `If[TrueQ[...], 0, 1]` ordering inequalities for the chi-window, J-window, hard-ceiling gap, and chi-vs-J fail-side relations. The banner string on line 32 still reads "STAGE 067" — a stale label inherited from an earlier copy; cosmetic only, not a math finding.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `zeta_phys` formula in `(Pe, eta, kappa)` with `Omega_Pe` form (notes §1.3) | A3 (line 76) pins `zeta_phys` at Family-1 in the `Pe -> oo` limit to the upstream numeric `zeta_max^(F1)` | match |
| `zeta_req` and `Q(zeta; eps_blk)` algebra; round-trip `Pi_tr = C_mix * Q` (notes §1.4) | A1 (line 53) verifies non-trivial cancellation `zetaReq(cMix * qMap) - zeta == 0` | match |
| `R_quad = zeta_req - zeta_phys` defined as scalar residual (notes §2) | `rQuad` defined and printed (lines 46, 51); no assertion — correct per paper, since `R_quad`'s sign IS the open theorem gap | match (open by design) |
| Family-1 specialization `kappa_F1 = 12321/5`, `eta_F1 = 37`, `Xi_F1 = 1369 Upsilon_w = 136900 Theta_w` (notes §3) | `kappaF1`, `etaF1`, `lambdaEll` declared (lines 55-57); A2 (line 65) cross-checks the two `Xi_F1` routes via `Upsilon_w -> 100*Theta_w` | match |
| Carried-forward numeric `zeta` windows and hard ceiling (notes §3) | A4-A7 (lines 84-87) verify chi-window ordering, hard-ceiling gap > 0, J-window ordering, J fail-side ≤ chi fail-side | match |
| `Xi*Delta_0 <= Pe_* <= Xi*Delta_inf` bounds (notes §2) | not exercised in this script | extra-upstream (carried in from stage 41); not required of stage 084's consolidation |
| Stage 084 Output: "structured reduced-PDE theorem chain and handoff" (`stage_084.tex:31`) | The five checks together exercise every named ingredient of the chain | match |

Set `paper_alignment: aligned`. No `paper_misalignment` finding; no script-side check is unmoored from the paper deliverables, and no paper deliverable lacks a corresponding script-side anchor (the unexercised `Xi*Delta_0 <= Pe_* <= Xi*Delta_inf` bound is carried in from stage 41 and is not a stage-084 deliverable per notes §2 wording).

## Assertion inventory

| #  | Script        | Line | Form                                                                                             | Exercises which paper claim?            | Anchored to claim? |
|----|---------------|------|--------------------------------------------------------------------------------------------------|------------------------------------------|--------------------|
| A1 | mathematica   | 53   | `expectZero["inverse demand map", (zetaReq /. piTr -> cMix*qMap) - zeta]`                        | demand-map inversion (notes §1.4)        | yes                |
| A2 | mathematica   | 65   | `expectZero["Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)", (xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta]` | Family-1 strength identity (notes §3) | yes (weakly substantive — see note below) |
| A3 | mathematica   | 76   | `expectApprox["zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)", zetaPhysF1Numeric, zetaMaxF1, 10^-10]` | physical zeta formula at Family-1 + numeric carry-forward | yes |
| A4 | mathematica   | 84   | `expectZero["chi-window ordering positive", If[TrueQ[zetaPlusChi > zetaMinusChi], 0, 1]]`         | chi-window ordering (notes §3)           | yes |
| A5 | mathematica   | 85   | `expectZero["hard-ceiling gap positive", If[TrueQ[zetaMaxF1 > zetaPlusChi], 0, 1]]`               | hard-ceiling structure (notes §3)        | yes |
| A6 | mathematica   | 86   | `expectZero["J-window ordering positive", If[TrueQ[zetaPlusJ > zetaMinusJ], 0, 1]]`               | J-window ordering (notes §3)             | yes |
| A7 | mathematica   | 87   | `expectZero["fail-side J below chi", If[TrueQ[zetaPlusJ <= zetaPlusChi], 0, 1]]`                  | J vs chi consistency (notes §3)          | yes |

A1 algebra (verified by hand): substituting `piTr -> cMix*qMap` into `zetaReq` produces numerator `cMix*z*(1-e)/(1-e*z)` and denominator `cMix*(1-e)/(1-e*z)`, ratio = `z`. The cancellation is non-trivial: it requires the specific (1, 1-2e, -e) and (1, -2e, +e) coefficient patterns in `Q` and `zetaReq` respectively to mate. A wrong sign on any coefficient breaks the round-trip.

A2 weakness: the check effectively asserts `lambdaEll^2 * 100 == 100 * lambdaEll^2`, i.e., that the same "100" appears as the `Upsilon_w / Theta_w` conversion factor in the substitution and as the prefactor in `xiF1FromTheta`. It catches a typo in one place but not the other, and does not anchor the upstream identity `Upsilon_w = 100 Theta_w` (which lives in the paper's macros, not in this script). This is a residual weakness from the v1 directive's required fix — the original v1 audit traded a strictly-tautological literal check for a weakly-substantive consistency check, which is a partial improvement. Not strong enough to file as a new finding under the status-only carve-out, since the substantive form of `Xi_F1 = lambdaEll^2 * Upsilon_w` is exercised in stage 082's SymPy script (see scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py) and the paper-side identity is sourced from notes §3 verbatim.

A3 verification: output line 19 prints `zeta_phys(Pe->oo, kappa_F1, eta_F1, y_F1) = 2.46752922945601223332958450157053542039` and line 20 reports `diff = 2.2333...e-15`, well below the `10^-10` tolerance. The `Limit::alimv` warning on line 18 (assumptions involving the limit variable are ignored) is benign — Mathematica drops the `pe > 0` assumption when `pe` is the limit variable, but the limit value is still the correct one (verifiable by L'Hopital or by writing `omegaPe ~ pi*pe*(2*pe*e^pe)/(4*pe^2*e^pe) = pi/2` as `pe -> oo`, giving `omegaPe^2 = pi^2/4`, hence `zetaPhys -> (pi^2/4)(kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)`).

A4-A7 verification: each `If[TrueQ[…], 0, 1]` returns `0` exactly when the inequality holds and `1` otherwise; `expectZero` then either passes or fails loudly. `TrueQ` ensures an undecidable comparison (e.g., due to mixed precisions) returns `False`, which converts to `1` and fails — no silent pass via undecidability.

## Findings

None. The v1 audit (`audit_date: 2026-05-22`) raised three low-severity findings (`tautological_check`, `hardcoded_result`, `insufficient_verification`); Codex applied all three on 2026-05-25 per the v1 directive's `## Applied:` blocks. The current script (mtime 2026-05-25 00:32) reflects those fixes, and the captured output (mtime 2026-05-25 00:33) confirms all seven assertions PASS and the script exits 0.

Attacks tried that failed:
- Inverse-demand-map round-trip (A1): worked out the algebra symbolically by hand; the cancellation requires the specific coefficient patterns in `Q` and `zetaReq`. A sign typo on either side would break it. Genuine identity.
- Pe->infinity limit (A3): verified the leading-order behavior `omegaPe -> pi/2` as `pe -> oo` (numerator `~ 2*pi*pe^2*e^pe`, denominator `~ 4*pe^2*e^pe`); the resulting `zetaPhys` formula `(pi^2/4)(kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)` with `y_F1*tan(y_F1) = 37` evaluates numerically to `2.467529229...` matching upstream. The `y_F1 ≈ 1.53...` is the lowest-positive root of `y*tan(y) = 37` (just below pi/2, where `tan(y)` is very large). The `FindRoot` initial guess `153/100` correctly selects this lowest branch.
- Cross-route `Xi_F1` (A2): checked whether a sign or factor typo in `xiF1FromTheta` (e.g., `50*lambdaEll^2*thetaW`) would be caught — yes, the substitution would produce `100*lambdaEll^2*thetaW - 50*lambdaEll^2*thetaW = 50*lambdaEll^2*thetaW`, nonzero. The check is weakly substantive (catches typos) but not deep.
- Ordering inequalities (A4-A7): each carries a meaningful physical interpretation (chi window contains the natural-shell saturation region; J window contains the conservative-floor region; chi's fail-side is closer to the hard ceiling than J's; both windows nest inside the hard ceiling). Together they exercise the "Family-1 effectively saturated" claim from notes §3.

Cosmetic note (not a finding, per the doc-alignment exclusion the v1 report cited): the banner string at line 32 of the script and line 3 of the output reads "STAGE 067 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON". The unit is 084 (the renumbering happened during the great-reorder commit `0d09ef6`). Not a math finding; flag for a sweep when a banner-renumbering pass is run, but do not fix in isolation.

## Independent-derivation check (Mathematica)

No SymPy counterpart for unit 084 exists (correctly, per the status-only carve-out). The upstream first-engine derivations live in:
- `scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py` — defines `zeta_max^(F1)` as the `Pe -> oo` limit of `A_F1 * Omega^2` with `A_F1 = (kappa_F1 + pi^2/4) / (kappa_F1 + y_F1^2)`. This is what A3 numerically anchors against in Mathematica.
- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py` — produces the five carried-forward numeric `zeta` values.
- `scripts/moving_throat_pde_stage082_master_quadrupole_residual_sympy_audit.py` — independent SymPy derivation of the demand-map round-trip and the `Xi_F1` strength identity.

The Mathematica script's algebra is structurally a faithful mirror of these upstream SymPy bodies for `zeta_req`, `Q(zeta; eps_blk)`, and `Xi_F1`. For a status-only consolidation, this is the *intended* behavior — the script consolidates upstream results in the second engine, demonstrating that the chain of identities holds in both engines. Not flagged as `mathematica_transliteration` because each engine's first-engine derivation is upstream of this unit, not co-located.

## Engine cross-check

n/a (no SymPy script for this unit). Engine agreement is established at the upstream stages (080-082): the numeric values printed by SymPy at those stages match the numeric values Mathematica computes here. Specifically, A3's `Pe -> oo` limit yields `2.46752922945601223...` which matches the upstream SymPy result `2.46752922945601` to ~15 significant figures — agreement is well below the `10^-10` tolerance.

## Verdict justification

`clean`. The v1 audit's three findings (tautological `Xi_F1` literal checks, hardcoded `expectApprox` gaps, unexercised `zetaPhys`/`rQuad`/`kappaF1`/`etaF1`) have all been addressed by the v1 directive's applied fixes. The current script exercises every paper-side deliverable that is verifiable in-stage: the demand-map inversion is non-trivially zero by hand-verified algebra; the `zetaPhys` formula is pinned to the upstream `zeta_max^(F1)` numeric to 15 digits in the `Pe -> oo` limit at Family-1 parameters; the Family-1 strength identity is cross-checked between the `Upsilon_w` and `Theta_w` routes; and four ordering inequalities exercise the "Family-1 effectively saturated" claim. The unverified `R_quad` is correctly left open — that is precisely the remaining theorem gap the paper card declares. Banner string `"STAGE 067"` is a cosmetic stale label, deferred to a renumbering sweep per the doc-alignment exclusion.

## Self-test notes

- Variable independence: `zetaPhys` at line 43 depends on `omegaPe(pe)`, `kappa`, `y`; the A3 substitution `{kappa -> kappaF1, eta -> etaF1, y -> yF1}` leaves only `pe` free, and the `Limit[…, pe -> Infinity]` then yields a numeric. Verified by tracing the symbol dependency from line 39 (omegaPe) through line 43 (zetaPhys) to line 73 (the substituted form).
- Symmetry/parity: not applicable — no integrals on unbounded domains in this script.
- Trivial-case pre-check: A1 substituting `piTr -> cMix*qMap` at `eps_blk = 0`, `zeta = 0` gives `zetaReq = (cMix*1 - cMix)/(cMix - 0) = 0`, matching `zeta = 0`. ✓ The general algebra cancels as shown above.
- Paper round-trip: every assertion traces to a notes-side deliverable; no new `paper_misalignment` risk introduced.
- The `Limit::alimv` warning on output line 18 was investigated and confirmed benign: `$Assumptions` constrains `pe > 0` but Mathematica drops that assumption inside `Limit[..., pe -> Infinity]`, which is correct behavior and produces the correct value (output line 19 prints the numeric, which matches upstream to 15 digits).
