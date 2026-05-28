---
unit_id: 166
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage166_bundle_inversion_four_drifts.md]
  paper_appendix: present
---

# Audit unit 166 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_166.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage166_bundle_inversion_four_drifts.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 63, 334-350)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}`: "Inverts the four observables \((\Theta_w,K_s,K_q,P_0)\) to determine the four remaining microscopic drifts." The notes spell out the deliverables. (D1) Four exact log-linear branch laws (notes §1): `dlnΘ_w = 2 dlnρ_w`, `dlnK_s = 2 dln a + dlnρ_w`, `dlnK_q = dlnZ_q + 2 dln c_s − 2 dln a`, `dlnP_0 = 5(dln c_s − dln a)`, themselves linearizations of `Θ_w = 25 λ²ρ_w²`, `K_s ∝ a²ρ_w`, `K_q ∝ Z_q c_s²/L_W²` (with `dln L_W = dln a`), `P_0 = (54G/5c⁵) c_s⁵/a⁵`. (D2) The exact inversion (notes §2): `dlnρ_w = ½ dlnΘ_w`, `dln a = ½ dlnK_s − ¼ dlnΘ_w`, `dln c_s = ½ dlnK_s − ¼ dlnΘ_w + ⅕ dlnP_0`, `dlnZ_q = dlnK_q − ⅖ dlnP_0`. (D3) The full-bundle rewrite (notes §3) with `dlnP_0 = dlnN_0 − dlnD_0`. (D4) The frozen-wall corollary (notes §4) at `dlnΘ_w = 0`, plus the explicit numeric `Θ_w^(χ) ≈ 4.06863235008162 λ²` and `ρ_w^(χ) = √(Θ_w^(χ)/25) λ⁻¹ ≈ 0.403417022451042 λ⁻¹`. This is not a checkpoint (`is_checkpoint: False`).

## What the script claims to verify

The docstring (sympy lines 6-12) states four checks: (1) solve the log inversion system `Θ_w,K_s,K_q,P_0 → ρ_w,a,c_s,Z_q`; (2) verify the forward transport identities after substitution; (3) verify the equivalent bundle form with `P_0 = N_0/D_0`; (4) report the frozen-wall simplification and `ρ_w^(χ)`. The script encodes the four log-laws as `eq1..eq4` (sympy 40-43), solves with `sp.solve` for `(drho, da, dcs, dZ)`, re-substitutes the solution into each law and asserts the residual is zero (53-56), asserts the bundle-form expressions equal the paper §3 targets (67-74), asserts the `dTheta→0` frozen expressions equal the paper §4 targets (82-85), and computes `ρ_w^(χ) = √(Θ_chi/25)` numerically (88-90). The Mathematica `.wl` mirrors this step-for-step.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Verdict |
|---|---|---|
| D1 — four log-linear branch laws | `eq1..eq4` (sympy 40-43) are *assumed* (the paper's §1 boxed linearized laws), not re-derived from the power-law forms; the constants 25, 54G/5c⁵ never enter | match (inputs carried as §1 boxed forms) |
| D2 — inversion: drho = ½ dTheta | printed (sympy 47) but **no assertion** against `dTheta/2` | partial |
| D2 — inversion: da = ½ dKs − ¼ dTheta | printed (48); only the `dTheta→0` limit `da = dKs/2` is asserted (frozen, 83); the `−¼ dTheta` coefficient is never asserted | partial |
| D2 — inversion: dcs = ½ dKs − ¼ dTheta + ⅕ dP | the `dTheta` and `dKs` coefficients are pinned by the bundle assertion (69: target includes `− dTheta/4`); dP coeff pinned by bundle/frozen | match |
| D2 — inversion: dZ = dKq − ⅖ dP | bundle assertion (73) pins both coefficients | match |
| D3 — bundle rewrite (dP→dN0−dD0) | asserted exactly (sympy 67-74) | match |
| D4 — frozen-wall corollary | asserted exactly (sympy 82-85) | match |
| D4 — ρ_w^(χ) numeric | computed `sqrt(Theta_chi/25)`, prints 0.403417022451042 = paper value | match |

Dominant pattern: aligned. Two D2 coefficients (drho's `½ dTheta`, da's `−¼ dTheta`) are only printed, not asserted — this is the `insufficient_verification` finding below, not a misalignment (the printed values are correct and match the paper).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `expect_zero(eq1.lhs.subs(sol)-eq1.rhs.subs(sol))` | D1/D2 forward consistency | partial (re-subs of solve output; near-tautological) |
| A2 | sympy | 54 | same for Ks law | D1/D2 | partial (near-tautological) |
| A3 | sympy | 55 | same for Kq law | D1/D2 | partial (near-tautological) |
| A4 | sympy | 56 | same for P0 law | D1/D2 | partial (near-tautological) |
| A5 | sympy | 67-70 | `dcs_bundle − (dKs/2 − dTheta/4 + (dN0−dD0)/5) == 0` | D3 + D2 dcs coeffs | yes |
| A6 | sympy | 71-74 | `dZ_bundle − (dKq − 2(dN0−dD0)/5) == 0` | D3 + D2 dZ coeffs | yes |
| A7 | sympy | 82 | `drho|frozen == 0` | D4 (drho at dTheta=0) | yes (but does not pin drho's dTheta slope) |
| A8 | sympy | 83 | `da|frozen − dKs/2 == 0` | D4 + D2 da dKs coeff | yes (dTheta slope NOT pinned) |
| A9 | sympy | 84 | `dcs|frozen − (dKs/2 + dP/5) == 0` | D4 + D2 dcs dKs,dP coeffs | yes |
| A10 | sympy | 85 | `dZ|frozen − (dKq − 2dP/5) == 0` | D4 + D2 dZ coeffs | yes |
| B1-B10 | mathematica | 48-79 | identical structure (`expectZero`) | same deliverables | same verdicts |

No assertion of the form `sol[drho] − dTheta/2 == 0` or `sol[da] − (dKs/2 − dTheta/4) == 0` exists in either script. The `½ dTheta` slope of drho and the `−¼ dTheta` slope of da are printed (sympy 47-48; math 42-43) but never asserted against the paper §2 boxed targets.

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py:47-56,82,83`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:42-51,76,77`

**What's wrong:**
The paper's primary deliverable (notes §2) is the general inversion, including `δlnρ_w = ½ δlnΘ_w` and `δln a = ½ δlnK_s − ¼ δlnΘ_w`. The script's only load-bearing assertions on `drho` and `da` are:
- the forward-verification block (sympy 53-56), which re-substitutes `sol` (the output of `sp.solve`) back into the same `eq1..eq4` it was solved from. For an exactly-solved linear system, sympy guarantees this residual is zero by construction, so these four checks are near-tautological — they confirm sympy's solver is self-consistent, not that `sol` equals the paper's boxed forms.
- the frozen-wall block (sympy 82-83 / math 76-77), which asserts `drho|_{dTheta=0} = 0` and `da|_{dTheta=0} = dKs/2`. Setting `dTheta → 0` deletes precisely the `dTheta`-dependent terms, so these assertions never exercise the `½ dTheta` slope of `drho` nor the `−¼ dTheta` slope of `da`.

Consequently the two coefficients `drho = ½ dTheta` and `da = … − ¼ dTheta` are only `print`ed (sympy 47-48), never asserted against the paper §2 targets. A solver bug or paper-mismatch in those two slopes would pass silently. (The `dcs` `−¼ dTheta` and `dZ` coefficients are, by contrast, pinned by the bundle assertions A5/A6, which retain the full `dTheta` term.)

**Why this matters:**
Two of the four boxed §2 results — the paper's headline deliverable — have no non-tautological script check. The verdict "Stage 166 inversion verified" overstates what the assertions actually pin.

**Required change:**
Add two assertions in each engine that compare the *general* solution to the paper §2 boxed forms (with `dTheta` present):
- sympy: `expect_zero("drho law", sol[drho] - sp.Rational(1,2)*dTheta)` and `expect_zero("da law", sol[da] - (sp.Rational(1,2)*dKs - sp.Rational(1,4)*dTheta))`. (Optionally also `dcs` and `dZ` general forms for completeness, but those are already pinned via A5/A6.)
- mathematica: `expectZero["drho law", drhoSol - dTheta/2]` and `expectZero["da law", daSol - (dKs/2 - dTheta/4)]`.
Place these in the first block, after the `Print` of the solved forms (sympy after line 50; math after line 45).

**Verification:**
New lines `drho law = 0` and `da law = 0` (sympy) / `PASS: drho law`, `PASS: da law` (mathematica) appear in the refreshed transcripts; both scripts still exit 0.

### F2 — mathematica_transliteration

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:26-79`
- (compare) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py:33-85`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation. Corresponding sections:
- Setup: sympy `eq1..eq4` (40-43: `dTheta == 2*drho`, etc.) ↔ mathematica `eq1..eq4` (31-34, identical RHS) with the same ten symbol names (`dTheta, dKs, dKq, dP, drho, da, dcs, dZ, dN0, dD0`).
- Solve: sympy `sp.solve((eq1..eq4),(drho,da,dcs,dZ))` (45) ↔ mathematica `Solve[{eq1..eq4},{drho,da,dcs,dZ},Reals]` (35) — same system, same unknowns, same order.
- Same `dP -> dN0 - dD0` bundle substitution (sympy 62-63 ↔ math 54-55), same `dTheta -> 0` frozen specialization (sympy 77 ↔ math 68-71), same target-form assertions, same numeric `Sqrt[Theta/25]`, same printed carry-forward formulas.
- Shared copy-paste artifact: **both** scripts open with the wrong banner `"STAGE 149 — …"` (sympy line 33; math line 26) for unit 166, strong evidence of duplication rather than independent authorship.

**Why this matters:**
The second-engine policy requires both engines to derive the result independently so a transcription error in one is caught by the other. A line-by-line port cannot catch a setup error common to both (e.g., a mistyped log-law coefficient would be copied verbatim). Mitigating factor: this stage is a 4×4 linear-algebra inversion with essentially one canonical setup, so the independent-derivation bar is intrinsically low here — hence low severity, not medium.

**Required change:**
Minimal, since the algebra is canonical: (a) fix the stale `"STAGE 149"` banner to `"STAGE 166"` in both `.py:33` and `.wl:26` so the two are not byte-for-byte duplicates of a wrong label; (b) in the Mathematica script, have the second engine cross-check the inversion by an independent route rather than echoing `Solve` over the same four equations — e.g., construct the 4×4 coefficient matrix `M` mapping `(drho,da,dcs,dZ) → (dTheta,dKs,dKq,dP)` and verify `Inverse[M]` reproduces the paper §2 coefficients, then confirm `M . (paper-§2 solution vector) == (dTheta,dKs,dKq,dP)`. This gives a genuinely distinct derivation (matrix inverse vs. equation solve) of the same claim.

**Verification:**
`.wl:26` and `.py:33` banners read "STAGE 166"; the Mathematica transcript shows a matrix-inverse cross-check (e.g., `Inverse check = 0` / `PASS`) distinct from the SymPy `solve` path; script still exits 0.

## Independent-derivation check (Mathematica)

Transliteration confirmed. The `.wl` reproduces the `.py` choreography one-for-one: identical symbols, identical `eq1..eq4`, identical `Solve`/`solve` over the same unknowns in the same order, identical bundle/frozen specializations, identical target-form assertions, identical numeric, and even the identical wrong `"STAGE 149"` banner. See F2 for quoted correspondences. Both engines do, however, produce algebraically identical final forms (verified in the cross-check below), so the transliteration has not introduced a discrepancy — the concern is methodological coverage, not a numeric disagreement.

## Engine cross-check

The two engines agree on every solved form and every assertion:

| Quantity | SymPy output | Mathematica output | Agree? |
|---|---|---|---|
| drho | `dTheta/2` | `dTheta/2` | yes |
| da | `dKs/2 - dTheta/4` | `(2*dKs - dTheta)/4` | yes (identical) |
| dcs | `dKs/2 + dP/5 - dTheta/4` | `(10*dKs + 4*dP - 5*dTheta)/20` | yes (identical) |
| dZ | `dKq - 2*dP/5` | `dKq - (2*dP)/5` | yes |
| forward laws | all 0 | all 0 / PASS | yes |
| bundle dcs, dZ | 0, 0 | 0, 0 / PASS | yes |
| frozen drho,da,dcs,dZ | 0,0,0,0 | 0,0,0,0 / PASS | yes |
| ρ_w^(χ) | 0.403417022451042341 | 0.40341702245104232684… | yes (~16 sig figs; last-digit float noise) |

No engine disagreement. Outputs are fresh: sympy `.py` mtime 11:56:51 < `.txt` 12:47:15; mathematica `.wl` 11:56:53 < `.txt` 13:19:29. Both report exit code 0.

## Verdict justification

The math is correct end to end. I independently re-solved the 4×4 log-linear system and confirmed all of drho = ½dTheta, da = ½dKs − ¼dTheta, dcs = ½dKs − ¼dTheta + ⅕dP, dZ = dKq − ⅖dP, the §3 bundle rewrite, the §4 frozen-wall forms, and the numeric ρ_w^(χ) = √(4.06863235008162/25) ≈ 0.403417 — all match the paper's boxed forms exactly. Attacks that failed: sign and factor-of-2 checks on every coefficient (all consistent with the power-law origins in §1); the bundle substitution `dP→dN0−dD0` is applied correctly; the frozen limit is the genuine `dTheta=0` specialization; the engines agree to float precision. Two real but non-fatal gaps survive: (F1) the headline §2 slopes of `drho` and `da` on `dTheta` are printed but never asserted — the only checks touching them are the near-tautological re-substitution of solve's own output and the `dTheta=0` frozen limit, which deletes exactly those terms; and (F2) the Mathematica script is a line-by-line transliteration of the SymPy one, sharing even the wrong "STAGE 149" banner, so a common setup error would not be caught. Neither is a `paper_misalignment` (every printed value matches the paper) and neither propagates downstream, so verdict is `findings`, `stop_cold: null`, alignment `aligned`.

## Self-test notes

I checked: (1) Variable independence — no `diff` in either script; the inversion is purely algebraic via `solve`/`Solve`, so no identically-zero-derivative trap. (2) Symmetry/parity — no integrals; n/a. (3) Trivial-case pre-check on the proposed F1 assertions — `sol[drho] - dTheta/2`: sympy returns `dTheta/2 - dTheta/2 = 0` (assert_zero passes correctly); `sol[da] - (dKs/2 - dTheta/4)`: returns `(dKs/2 - dTheta/4) - (dKs/2 - dTheta/4) = 0` (passes correctly); both reduce to 0 only when the solver gives the paper form, so they genuinely fail on a wrong slope. (4) Paths — `.py` fix in `scripts/`, `.wl` fix in `mathematica/`, both pre-existing files. (5) Paper round-trip — the F1 targets `dTheta/2` and `dKs/2 - dTheta/4` are copied verbatim from notes §2 boxed equations; no new constant introduced; F2's matrix-inverse cross-check reproduces the same §2 coefficients, introducing no new value.
