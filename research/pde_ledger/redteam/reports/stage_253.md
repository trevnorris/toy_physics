---
unit_id: 253
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 253 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_253.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows 104, 333-340, 369 reference this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage253_physical_calibration_and_material_threshold_companion_from_the_stage252_export_and_cold_survival_compiler_mathematica_audit.txt`

## What the paper claims

Stage 253 is the physical-calibration / materials companion that maps the exact Stage-252 microscopic safe-edge lattice rate into condensed-matter threshold variables. The deliverables (paper card subsecs and notes §§1-5) are: (1) the lattice-turnover threshold `(λ_ep ω_D)_min^(253) = f_lat(s_c) μ_η (s_0²−s_c²) / (Υ_lat ζ_ep s_c t_*)` and its event-product form `(λ_ep ω_D t_cross^phys)_min^(253) = f_lat μ_η (s_0²−s_c²)/(Υ_lat ζ_ep s_c²)`; (2) recovery of the legacy Session-V slice via `Υ_lat^(sess) = γ_lat,safe^eq / γ_lattice^red`; (3) the harmonic geometric trigger `χ_{λ,lattice}(r_turn) = 2λ_phys/r_turn^phys = 2λ_ref/r_turn` and the force-matched stiffness `k_eff,req = K_turn E_*/λ_phys² = 4 K_turn E_*/a_int²` with `K_turn = (λ_ref²/r_turn)|V'_red(r_turn)|`; (4) the Korringa ceiling `T_max = K_corr/t_cross^phys = s_c K_corr/t_*`; (5) the four screening ratios `Π_ep, Π_χ, Π_k, Π_T` and the screen `Π_• ≥ 1`. The card's final-screen equation (the stage's bottom-line `Output`) is the inequality stack `Π_ep ≥ 1, Π_χ ≥ 1, Π_k ≥ 1, Π_T(T) ≥ 1`, explicitly stated to be a reduced threshold *criterion*, not a claim that any host passes. Benchmark subsections (§2.2, §3.2, §4) carry declared decimal values, including the disputed `(λ_ep ω_D t_cross^phys)_min^(micro) ≈ 187.23361317/ζ_ep`.

## What the script claims to verify

The SymPy docstring (lines 7-13) states five verified items: the exact lattice-turnover threshold, the Session-V slice recovery via Υ_lat, the χ_λ + force-matched k_eff compiler, the Korringa ceiling, and the exact screening-ratio identities, followed by the benchmark readback. The assertions test exactly these as symbolic identities: each compiler form is built by composing the calibration definitions and then `simplify`-compared against the closed form the paper states (e.g. line 68/115/116/133/152). Section 6 (lines 159-198) then re-evaluates the benchmark decimals and asserts them against hardcoded targets. The Mathematica script (`expectZero`/`expectNear`) covers the identical list section-for-section. Neither script tests any `Π ≥ 1` inequality (correctly — those are host-screening criteria with no host data in scope).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (λ_ep ω_D)_min and product threshold (§2) | py 68-70 / wl 113-118 | match |
| Session-V slice recovery via Υ_lat (§2.1) | py 87-88 / wl 139-146 | match |
| χ_λ geometry ratio = 2λ_ref/r_turn (§3) | py 115 / wl 177 | match |
| force-matched k_eff,req, K_turn, a_int form (§3.1) | py 116-118 / wl 180-184 | match |
| Korringa ceiling T_max = s_c K_corr/t_* (§4) | py 133 / wl 195 | match |
| four screening ratio identities Π_ep,Π_χ,Π_k,Π_T (§5) | py 152-153, 143 / wl 215-217 | match |
| final-screen inequalities Π_• ≥ 1 (§5.4, card Output) | none (criteria, no host data) | match (correctly not provable) |
| benchmark γ_safe ≈ 65.45, Υ^(sess) ≈ 13.65, χ ≈ 2.19, K_turn, T_max coeff | py 190-198 / wl 252-260 | match |
| benchmark micro product ≈ 187.23361317 (§2.2) | py 192 asserts 119.23361317; wl 254 asserts 119.23361317 | mismatch |
| benchmark a_int stiffness coeff 10.95423247 (§3.2) | py 196 asserts 10.95423248; wl 258 asserts 10.95423248 | mismatch (1-ulp) |

`paper_alignment: misaligned` — the symbolic theorem path matches exactly, but a load-bearing benchmark decimal in §2.2 disagrees with both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 68 | `simplify(threshold_lambda − expected) == 0` | claim 1 (λω threshold) | yes |
| A2 | sympy | 69 | `simplify(threshold_product − expected) == 0` | claim 1 (product) | yes |
| A3 | sympy | 70 | `simplify(threshold_product − γ_safe/(Υ ζ s_c)) == 0` | claim 1 (reduced) | yes |
| A4 | sympy | 87-88 | legacy slice recovery `== γ_legacy/(ζ t_*)`, `== γ_legacy/(ζ s_c)` | claim 2 | yes |
| A5 | sympy | 114 | `simplify(r_turn_phys − lambda_phys*r_turn/lambda_ref) == 0` | claim 3 (radius map) | **no (tautological)** |
| A6 | sympy | 115 | `simplify(chi − 2*lambda_ref/r_turn) == 0` | claim 3 (χ ratio) | yes |
| A7 | sympy | 116-118 | k_eff_req / K_turn / a_int forms | claim 3 (stiffness) | yes |
| A8 | sympy | 133-134 | `T_max == s_c K_corr/t_*`, `Pi_T == K_corr/(T t_cross)` | claim 4 | yes |
| A9 | sympy | 152-153 | `Pi_ep == λω/threshold_lambda`, `Pi_k == k_eff/k_eff_req` | claim 5 | yes |
| A10 | sympy | 190-198 | 9 numeric benchmark asserts vs hardcoded targets | benchmark | partial (see F1) |
| B1 | mathematica | 113-118 | `expectZero` lattice-turnover trio | claim 1 | yes |
| B2 | mathematica | 139-146 | `expectZero` legacy recovery pair | claim 2 | yes |
| B3 | mathematica | 176 | `expectZero["turning-point radius map", rTurnPhys − lambdaPhys rTurn/lambdaRef]` | claim 3 | **no (tautological)** |
| B4 | mathematica | 177-184 | χ, k_eff, K_turn, a_int | claim 3 | yes |
| B5 | mathematica | 195-196 | Korringa ceiling + Pi_T | claim 4 | yes |
| B6 | mathematica | 215-217 | Pi_ep, Pi_chi, Pi_k | claim 5 | yes |
| B7 | mathematica | 252-260 | 9 `expectNear` benchmark checks | benchmark | partial (see F1) |

23 SymPy asserts and 25 Mathematica checks — both far exceed the 5-item docstring; count is sufficient.

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** value_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_253.tex:274` (and notes line 274)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage253_..._sympy_audit.py:192`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_..._mathematica_audit.wl:254`

**What's wrong:**
The paper card §2.2 states the raw-microscopic (Υ_lat=1) event-product threshold as
`(λ_ep ω_D t_cross^phys)_min^(micro) ≈ 187.23361317/ζ_ep` (stage_253.tex:274, notes:274).
But this quantity is exactly `γ_lat,safe^eq / s_c = 65.45193926 / 0.5489386551`. Both engines compute and assert `119.23361317`:
- py:192 `assert abs(product_micro_num - 119.23361317476524) < 1e-9` (output line 56: `micro product thresh = 119.23361317477 / zeta_ep`)
- wl:254 `expectNear["micro product threshold", productMicroNum, 119.23361317476524, 10^-9]` (output line 92: `119.23361317476522`)
The shared trailing digits `…23361317` between `187.23361317` and `119.23361317` indicate a leading-digit transcription typo in the paper/notes (119 → 187). The script value is the mathematically correct one. (Note: §2.2 also displays an intermediate constant; `gamma_lat,safe^eq ≈ 65.45193926`, which both engines confirm — only the final 187 is wrong.)

A second, much smaller mismatch in the same family: the a_int stiffness coefficient. Paper §3.2 / notes line 419 give `10.95423247`; both scripts assert `10.95423248` (= `4 × K_turn = 4 × 2.73855812`, exact). This is a 1-ulp last-digit rounding difference, not a real disagreement, but the published digit `…247` should be reconciled to `…248` for consistency.

**Why this matters:**
A reader auditing the benchmark would find the paper's headline `187` number cannot be reproduced from the stage's own inputs. Because this is a published-card value, per audit policy the resolution direction (fix paper vs. fix script) is the user's call, not Codex's. Both engines agree on `119.23…`, so the weight of evidence is that the paper is wrong, but I do not prescribe overwriting the card.

**Required change:**
See `## Resolve before fix_loop` in the directive. Do NOT auto-edit either side.

**Verification:**
After user chooses a direction: if paper is wrong, the card line 274 (and notes 274) becomes `119.23361317` and no script change is needed (scripts already assert 119.23…); if (implausibly) the scripts are wrong, py:192 / wl:254 targets change and the scripts are re-run. The a_int coefficient sub-item: reconcile card/notes `10.95423247` to `10.95423248`.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/.../mathematica/moving_throat_pde_stage253_..._mathematica_audit.wl` (whole file)

**What's wrong:**
The `.wl` is a section-for-section, line-for-line port of the `.py`, not an independent re-derivation. Identical decomposition, identical intermediate variable choreography, identical substitution order. Three corresponding excerpts:

- py:49-51
  ```
  gamma_lat_safe_eq = sp.simplify(f_lat * mu_eta * (s0**2 - s_c**2) / s_c)
  t_cross_phys = sp.simplify(t_star / s_c)
  gamma_lat_turn_phys = sp.simplify(gamma_lat_safe_eq / (Upsilon_lat * t_star))
  ```
  wl:90-95
  ```
  gammaLatSafeEq = FullSimplify[fLat muEta (s0^2 - sc^2)/sc, ...];
  tCrossPhys = FullSimplify[tStar/sc, ...];
  gammaLatTurnPhys = FullSimplify[gammaLatSafeEq/(UpsilonLat tStar), ...];
  ```
- py:104-105 (substitution sequence)
  ```
  k_eff_req_Kturn = sp.simplify(k_eff_req.subs(Vprime_turn_abs, K_turn_sym * r_turn / lambda_ref**2))
  k_eff_req_aint = sp.simplify(k_eff_req_Kturn.subs(lambda_phys, a_int / 2))
  ```
  wl:160-167 (same two substitutions, same order)
  ```
  kEffReqKTurn = FullSimplify[kEffReq /. VprimeTurnAbs -> KTurnSym rTurn/lambdaRef^2, ...];
  kEffReqAInt = FullSimplify[kEffReqKTurn /. lambdaPhys -> aInt/2, ...];
  ```
- py:142 vs wl:200-203 — `Pi_ep` written with the identical factor ordering.

Every `assert simplify(a−b)==0` maps to a single `expectZero[name, a−b]`; every benchmark `assert abs(x−t)<tol` maps to a single `expectNear[name, x, t, tol]`. This violates the second-engine independence policy: the second engine should derive the compiler from the physical premises by a different route (e.g., native dimensional substitution, or building r_turn^phys from the energy/length maps and verifying χ and k_eff via independent `Solve`/limit primitives), not echo the SymPy algebra.

**Why this matters:**
A transliteration cannot catch an algebra error that is shared by construction — it only re-confirms SymPy's own steps in a second syntax. For a checkpoint capstone, the dual-engine guarantee is supposed to be a genuine cross-check.

**Required change:**
Re-author `mathematica/moving_throat_pde_stage253_..._mathematica_audit.wl` as an independent route: derive r_turn^phys and χ_λ from the length map `r^phys = (λ_phys/λ_ref) r` and the harmonic log-derivative `∂_r ln V = 2/r_phys` directly (not by transcribing the cancelled `2 λ_ref/r_turn`), build k_eff,req by force-matching `E_* (λ_ref/λ_phys)|V'_red| = k_eff r_turn^phys` symbolically and solving for k_eff, and verify the threshold/Korringa/screening forms by independent algebra. Keep the same final claims (so engine_disagreement stays clean), but reach them by a non-mirrored path.

**Verification:**
The re-authored `.wl` still exits 0 with all `expectZero`/`expectNear` passing, but its intermediate construction no longer mirrors the `.py` line-by-line; verifier confirms a distinct decomposition (e.g., a `Solve[...]` step where the `.py` used `subs`).

### F3 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage253_..._sympy_audit.py:114`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_..._mathematica_audit.wl:176`

**What's wrong:**
`r_turn_phys` is defined as `sp.simplify(lambda_phys * r_turn / lambda_ref)` (py:99) and then asserted equal to the literal same expression at py:114:
```
r_turn_phys = sp.simplify(lambda_phys * r_turn / lambda_ref)
...
assert sp.simplify(r_turn_phys - lambda_phys * r_turn / lambda_ref) == 0
```
The residual is `expr − expr ≡ 0` by construction; the assertion cannot fail regardless of the physics. The Mathematica twin (wl:150 define, wl:176 `expectZero["turning-point radius map", rTurnPhys - lambdaPhys rTurn/lambdaRef]`) is the same tautology. This is the one place where the radius-map deliverable (`r_turn^phys = (r_turn/λ_ref) λ_phys`) is "checked" — and the check is vacuous. (The downstream χ and k_eff checks at py:115-116 DO exercise the map non-tautologically via cancellation, so the map is indirectly covered; but the dedicated radius assertion itself proves nothing.)

**Why this matters:**
A reader trusting the assertion list would believe the radius-map equation was independently verified; it is not. Low severity because χ/k_eff implicitly re-use the map and would fail if it were mis-typed there.

**Required change:**
Make the radius-map check non-tautological by anchoring it to the physical-unit-map premises rather than restating the result. Concretely, define `r_turn_phys` from the chain `r^phys = (λ_phys/λ_ref) r` evaluated at `r = r_turn` (the map), and separately assert it equals the *independent* expectation derived from the inverse map used in §3.1 (py:101/wl:153 use `r = (λ_ref/λ_phys) r_phys`). I.e. verify round-trip consistency: substituting `r_turn_phys` back through the inverse map returns `r_turn`:
- py:114 → `assert sp.simplify((lambda_ref/lambda_phys) * r_turn_phys - r_turn) == 0`
- wl:176 → `expectZero["turning-point radius round-trip", (lambdaRef/lambdaPhys) rTurnPhys - rTurn]`
This is a real check: it fails if the forward and inverse length maps are inconsistent.

**Verification:**
py:114 / wl:176 now reference the inverse-map composition; residual is still 0 (forward∘inverse = identity) but the assertion is no longer `expr−expr`. Scripts exit 0.

### F4 — stale label (informational)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage253_..._mathematica_audit.wl:65`

**What's wrong:**
The Mathematica banner reads `banner["STAGE 236 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"]` (wl:65; reproduced verbatim in the output transcript line 11). This is a stale stage number — the unit is 253. Cosmetic only (a `Print` string, not load-bearing), but it would mislead anyone reading the transcript and is the kind of stale-renumber artifact this project tracks.

**Why this matters:**
A reviewer scanning the Mathematica transcript sees "STAGE 236" and may mis-file the artifact. No mathematical impact.

**Required change:**
wl:65 → `banner["STAGE 253 — PHYSICAL CALIBRATION AND MATERIAL-THRESHOLD COMPANION"];`

**Verification:**
Re-run output banner line reads "STAGE 253".

## Independent-derivation check (Mathematica)

Not independent — see F2. The `.wl` reproduces the `.py` section ordering, variable names (transliterated camelCase), intermediate-expression sequence, and substitution order exactly. It is a transliteration mirror, flagged as `mathematica_transliteration`.

## Engine cross-check

Both engines agree at every check. Symbolic forms match (SymPy `f_lat*mu_eta*(s_0**2 - s_c**2)/s_c` vs Mathematica `(fLat*muEta*(s0 - sc)*(s0 + sc))/sc` — same value, factored differently). Benchmark numerics agree to ≤1.4e-14 (Mathematica output lines 100-117). No `engine_disagreement`. Critically, both engines assert the SAME `119.23361317` micro-product value, so the paper's `187.23361317` disagrees with BOTH engines, strengthening F1.

## Verdict justification

The symbolic theorem path (the five calibration/compiler identities and the four screening-ratio definitions) holds up under attack in both engines and faithfully matches the paper card and notes — I tried to break the χ-ratio cancellation (py:115), the lambda_ref²/lambda_phys² squaring in k_eff (py:116), the factor-of-4 a_int rewrite (py:118), and the Pi_ep↔threshold reciprocity (py:152), and all are genuine, non-tautological, well-anchored checks. The verdict is `findings`, not `clean`, because of four issues: a high-severity `paper_misalignment` where the paper's headline micro-product benchmark `187.23361317` cannot be reproduced (both engines give `119.23361317`; routed to user); a `mathematica_transliteration` (the second engine is a line-by-line port, not an independent route); one genuinely tautological assertion (py:114 / wl:176, the radius map asserted against its own definition); and a stale "STAGE 236" banner. Outputs are fresh (script mtime 11:58 < output mtimes 12:52/13:27). Not `stop_cold`: the paper_misalignment is a benchmark decimal that does not feed any downstream unit's algebra (Stage 253 is an end-of-ladder companion, and the disputed value is a calibration readback, not a forward constant), so it is not CRITICAL_DOWNSTREAM; and nothing is internally contradictory, so not UNFIXABLE.

## Self-test notes

Checked the variable-independence trap: F3's proposed round-trip check `(lambda_ref/lambda_phys)*r_turn_phys - r_turn` genuinely depends on lambda_ref, lambda_phys, r_turn (no spurious derivative), and reduces to 0 only because forward∘inverse = identity — non-vacuous. Checked the inequality-strictness warning: confirmed neither script tests any `Π ≥ 1`, and that is correct (screening criteria, no host data), so no non-strict-threshold defect exists. Checked hardcoded-result trap: the §6 benchmark literals are declared benchmark inputs (s_c, s_0, K_turn, etc.) re-combined arithmetically and asserted against hardcoded outputs — legitimate readback, except the paper-side `187` target (F1). Confirmed paper round-trip: F3's fix introduces no new constant and stays on the published maps; F4 is a string-only change.
