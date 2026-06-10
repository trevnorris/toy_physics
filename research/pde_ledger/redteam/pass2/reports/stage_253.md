---
unit_id: 253
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: ["moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 253 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_253.tex` (140 lines, abstract/symbolic)
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows 104, 333–341)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.txt`

## What the paper claims

This is THE FINAL STAGE — a physical-calibration / materials companion (explicitly "not part of the core PDE theorem ladder"). It introduces the unit dictionary (`t^phys = t_* t`, `r^phys = (lambda_phys/lambda_ref) r`, `E^phys = E_* E`) and maps the exact Stage-252 safe-edge lattice rate `gamma_lat,safe^eq = f_lat(s_c) mu_eta (s_0^2 - s_c^2)/s_c` into condensed-matter thresholds via one transport projection factor `Upsilon_lat > 0`. Distinct deliverables: (1) the exact lattice-turnover threshold `(lambda_ep omega_D)_min^(253)` and its event-product form; (2) recovery of the legacy Session-V slice via `Upsilon_lat^(sess) = gamma_lat,safe^eq / gamma_lattice^red`; (3) the harmonic geometric trigger `chi_lambda = 2 lambda_ref / r_turn` separated from the force-matched stiffness `k_eff,req = K_turn E_* / lambda_phys^2` (and `4 K_turn E_*/a_int^2`); (4) the Korringa ceiling `T_max = s_c K_corr / t_*`; (5) the four screening ratios `Pi_ep, Pi_chi, Pi_k, Pi_T`, each `>= 1`. The card is abstract; the benchmark decimals live only in the notes §2.2/§3.2/§4. Notes-side benchmark deliverables: `gamma_lat,safe^eq ≈ 65.45193926`, micro product `119.23361317/zeta_ep`, `Upsilon_lat^(sess) ≈ 13.64824695`, legacy product `8.73618521/zeta_ep`, `chi_lambda ≈ 2.19084649`, `K_turn ≈ 2.73855812`, a_int stiffness coeff `10.95423248`.

## What the script claims to verify

Both engines verify all five symbolic compiler families plus the benchmark specialization layer. SymPy asserts: the threshold/product compilers equal their expected closed forms; the legacy-slice recovery; the chi_lambda geometry ratio, K_turn definition, and the two k_eff,req rewrites; the Korringa ceiling and Pi_T; the Pi_ep/Pi_k ratios; and 9 numeric benchmark asserts. Mathematica mirrors the same families through `expectZero`/`expectNear`. The benchmark numerics are declared inputs (`s_c, s_0, f_lat, mu_eta, gamma_lattice_legacy, lambda_ref, r_turn, K_turn`) feeding derived outputs (`gamma_lat,safe^eq`, products, Upsilon, chi, a_int coeff, V'). Inputs are carry-forward (Stage-252/Session-V benchmark) and so are legitimately literal here.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Threshold `(lambda_ep omega_D)_min^(253)` + product form | py L55–70 / wl L96–123 | match |
| Legacy Session-V recovery via `Upsilon_lat^(sess)` | py L76–88 / wl L127–151 | match |
| chi_lambda geometric trigger = 2 lambda_ref/r_turn | py L100,115 / wl L155–165,199–200 | match |
| Force-matched k_eff,req (+K_turn, +a_int forms) | py L101–118 / wl L166–207 | match |
| Korringa ceiling T_max = s_c K_corr/t_* | py L126–133 / wl L211–221 | match |
| Four screening ratios Pi_ep/Pi_chi/Pi_k/Pi_T >= 1 | py L142–153 / wl L226–243 | match |
| Benchmark values (notes §2.2/§3.2/§4) | py L159–198 / wl L247–286 | match |

paper_alignment: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 68–70 | `simplify(thr - expected)==0` (x3) | claim 1 | yes |
| A2 | sympy | 87–88 | legacy recovery ==0 | claim 2 | yes |
| A3 | sympy | 114–118 | chi/k_eff/K_turn/a_int ==0 | claim 3 | yes |
| A4 | sympy | 133–134 | T_max, Pi_T ==0 | claim 4 | yes |
| A5 | sympy | 152–153 | Pi_ep, Pi_k ==0 | claim 5 | yes |
| A6 | sympy | 190–198 | 9 numeric `abs(.-lit)<tol` | benchmark | yes |
| B1 | math | 118–123 | `expectZero` threshold (Solve-derived) | claim 1 | yes |
| B2 | math | 144–151 | `expectZero` legacy | claim 2 | yes |
| B3 | math | 199–207 | `expectZero` chi/k_eff/K_turn/a_int | claim 3 | yes |
| B4 | math | 221–222 | `expectZero` Korringa/Pi_T | claim 4 | yes |
| B5 | math | 241–243 | `expectZero` Pi_ep/Pi_chi/Pi_k | claim 5 | yes |
| B6 | math | 278–286 | 9 `expectNear` | benchmark | yes |

All rows non-tautological: the `expectZero` left side is built by an independent route (Solve/D) and compared to the posited closed form, so a wrong derivation would yield a nonzero residual.

## Findings

None. The audit is clean.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent re-derivation, NOT a transliteration. Across every load-bearing object the two engines use materially different information flows, with the `.wl` DERIVING (via `Solve`/`D`) what the `.py` POSITS by direct algebra:

- **Lattice-turnover threshold (the central deliverable).** `.py` L51,53 *posits* `gamma_lat_turn_phys = gamma_lat_safe_eq/(Upsilon_lat*t_star)` and `threshold_lambda = gamma_lat_turn_phys/zeta_ep` by division. `.wl` L96–97 instead forms the *physical balance equation* `lambdaBalance = gammaLatTurnPhys == zetaEp lambdaEpOmegaD` and `Solve[lambdaBalance, lambdaEpOmegaD]` — it solves the condensed-matter identification for the threshold rather than dividing. Different operation, same closed form.
- **chi_lambda (cleanest independence).** `.py` L100 hardcodes the analytic harmonic log-slope as `2*lambda_phys/r_turn_phys`. `.wl` L157–165 defines `V[rad] = (1/2) kEff rad^2`, computes `radialLogSlope = D[Log[V[r]], r]` (an actual symbolic derivative → `2/r`), evaluates at `rTurnPhys`, multiplies by `lambdaPhys`. The `.wl` derives the `2/r` factor; the `.py` asserts it.
- **Force-matched stiffness.** `.py` L101 *posits* `k_eff_req = E_star*lambda_ref*Vprime/(lambda_phys*r_turn_phys)`. `.wl` L166–170 forms `forceMatch = EStar(lambdaRef/lambdaPhys)VprimeTurnAbs == kEff rTurnPhys` and `Solve[forceMatch, kEff]`.
- **Korringa ceiling.** `.py` L126 divides `Kcorr/t_cross_phys`. `.wl` L211–214 solves `KCorr == TCeiling tCrossPhys` for `TCeiling`.
- **K_turn → V' and a_int substitutions.** `.py` uses `.subs(...)` (L104–105); `.wl` uses `Solve[KTurnSym==KTurn, VprimeTurnAbs]` and `Solve[aInt==2 lambdaPhys, lambdaPhys]` (L175–186).
- **Pi_ep.** `.py` L142 builds `Upsilon*zeta*lambda*t_star/gamma_safe`; `.wl` L226–229 builds `zetaEp lambdaEpOmegaD/gammaLatTurnPhys` — distinct starting forms, both checked equal to `lambda_ep_omega_D/threshold_lambda`.

This pass-1 re-author from a transliteration is SUFFICIENT (like VI.1-218 / VII.1-221, unlike V.3-200): the Solve/D routes are not cosmetic syntax swaps — they re-derive the load-bearing forms from the physical equalities.

## Engine cross-check

Final symbolic forms agree (modulo CAS factoring: `.py` keeps `(s_0**2 - s_c**2)`, `.wl` factors to `(s0-sc)(s0+sc)` and `.wl` Pi_ep prints `f_lat*mu_eta*s0^2 - f_lat*mu_eta*sc^2` expanded — same expression). All `expectZero` residuals are 0; all `expectNear` diffs are 0 or ~1e-14, well under tolerances. Benchmark numerics agree: micro product `119.23361317476524` (py) vs `119.23361317476522` (wl, diff 1.4e-14 < 1e-9); `gamma_lat,safe^eq 65.45193925961132`; `Upsilon 13.64824695299483`; a_int coeff `10.95423248`; `T_max coeff 0.5489386551062235`. Both scripts report all checks passed.

## Verdict justification

Clean. Both engines present and substantive (no status-only carve-out taken; checkpoint bar MET). Attacks tried and failed: (a) tautology — the `expectZero`/assert left sides are built by an independent route (Solve/D) then differenced against the posited closed form, so they can fail; (b) hardcoded result — the only literals are declared Stage-252/Session-V benchmark INPUTS, with every benchmark OUTPUT derived from them; (c) symbol-assumption error — all symbols are `positive, real`, matching the physical setup (rates, widths, energies, temperatures all > 0); positivity is used only for `FullSimplify` of manifestly-positive quotients, not to hide a branch; (d) variable-independence trap — `D[Log[V[r]], r]` with `V` quadratic in `r` genuinely depends on `r` (→ `2/r`, nonzero), so no identically-zero derivative slipping an assert; (e) paper misalignment — card is abstract and matches notes; benchmark deliverables reconcile (below). Paper, notes, and appendix all read; the script claims match the paper claims exactly.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 11 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `gamma_lat,safe^eq ≈ 65.45193926` | py out L46 / wl out L82 | notes:266 `65.45193926` | MATCH |
| micro product `119.23361317(476524)/zeta_ep` | py out L48 / wl out L84 | notes:274 `119.23361317` | MATCH |
| `Upsilon_lat^(sess) ≈ 13.64824695` | py out L47 / wl out L83 | notes:284 `13.64824695` | MATCH |
| legacy product `8.73618521/zeta_ep` | py out L49 / wl out L85 | notes:244 `8.73618521` | MATCH |
| `gamma_lattice^red = 4.79562976` (input) | py L163 / wl L251 | notes:222,229,237,282 | MATCH |
| `r_turn^phys coeff ≈ 0.9128891530` | py out L50 / wl out L86 | notes:391 `0.9128891530` | MATCH |
| `chi_lambda ≈ 2.19084649` | py out L51 / wl out L87 | notes:398 `2.19084649` | MATCH |
| `K_turn ≈ 2.73855812` (input) | py L166 / wl L254 | notes:403 `2.73855812` | MATCH |
| a_int stiffness coeff `10.95423248` | py out L53 / wl out L89 | notes:419 `10.95423248` | MATCH |
| `T_max coeff = s_c ≈ 0.548938655` | py out L54 / wl out L90 | notes:457 `0.548938655` | MATCH |
| `s_c ≈ 0.5489386551` (input) | py L159 / wl L247 | notes:259,451 | MATCH |

STALE-VALUE SCAN: explicit grep for `187.2336`, `136.2336`, `10.95423247` across card, notes, `.py`, `.wl`, and BOTH outputs returned ZERO hits in every file. The pass-1 NOTES-ONLY corrections hold: notes:274 = `119.23361317` (not 187.x), notes:419 = `10.95423248` (not ...247). The scripts emit `119.23361317476524` / `119.23361317476522` and `10.95423248`; notes match. No stale tracker `136.2336` survives in any stage-253 file. Published card is CLEAN (140 lines, abstract — carries no benchmark decimals, so the value is correctly absent there; pass-1's misattribution to `stage_253.tex:274` was wrong and is NOT re-flagged). `.wl` banner line 65 = `STAGE 253 …` — canonical (no stale stage number). Internal scaffolding (no finding): `s_0 = 6.94311167` input, `|V'_red| = 5.837462857946154` (derived intermediate also asserted, present in both notes-implied and scripts), `mu_eta = 1.0`, `f_lat = 0.75` benchmark inputs.

## Self-test notes

Checked: (1) variable-independence — `D[Log[V[r]], r]` depends on `r` (V quadratic), derivative is `2/r` nonzero, no zero-derivative trap; the `.wl` Solve targets (`lambdaEpOmegaD`, `kEff`, `TCeiling`, `VprimeTurnAbs`, `lambdaPhys`) each appear in their respective equations, so Solve returns a genuine non-trivial root. (2) No unbounded-domain integrals here, so parity trap n/a. (3) Trivial-case: all `expectZero` differences reduce to literal 0 in both outputs; the benchmark `expectNear` diffs are 0 or 1e-14. No directive written (zero findings).
