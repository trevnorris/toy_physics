---
unit_id: 247
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 5
paper_alignment: misaligned
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 247 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_247.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows 40, 55, 92, 117, 138, 239-249 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.txt` (mtime 2026-05-11 12:52, newer than script 11:58 → fresh)
- mathematica output: (missing)

## What the paper claims

Stage 247 assembles the first stationary lowered-barrier law from a carried one-port short-range baseline plus three imported "relaxed" packets (Stage-244 leakage/work, Stage-245 U/V drain, Stage-246 compensated source). The paper's `Output` (notes Section 6, appendix theorem item 6) is the compiler `V_eff^(247)(r) = V_short^(1p)(r) - lambda_L S_leak - lambda_W W_w^sess - DeltaE_UV - M_sigma` together with the exact lowering identity `V_short^(1p) - V_eff^(247) = lambda_L S_leak + lambda_W W_w^sess + DeltaE_UV + M_sigma`, and the consequent inequality `V_eff^(247) <= V_short^(1p)` whenever the imported packet scalars and weights are nonnegative. Distinct deliverables: (1) the one-port baseline `V_short^(1p)` with `A_6 = 3 alpha_6 + C_6/2`, `A_4 = C_4`, `A_2 = alpha_2 + C_2/2` and the susceptibilities as entries of the inverse reduced stiffness matrix; (2) the Stage-244 leakage/work scalars; (3) the Stage-245 weighted drain `DeltaE_UV = eta_UV D_UV` with `D_UV >= 0`; (4) the Stage-246 source-response `M_sigma = xi_R[R_inf - R(r)]` with the supporting facts `g_inf = 2/pi`, `g(r) >= 2/pi`, and `M_sigma >= 0` on the Session-I branch; (5) the lowering identity AND the inequality; (6) the Session-I benchmark at `r_soft = 0.18` decomposing `3.74163698 - 0.26971918*0.31069599 - 1.51632107 - 0.21064278 - 0.18386120 = 1.74701126`. The notes/card additionally state the benchmark intermediate `Delta = 210.17750000` and `D0 = 3.76481862` (stage_247.tex:407-409).

## What the script claims to verify

The docstring enumerates six items mirroring the paper deliverables. Concretely the assertions test: (Sec 1) six susceptibility expressions equal the corresponding entries of `K_red.inv()`; (Sec 2) `S_leak` and `W_sess` match a factored form through `E0`, and `Lvar` is recoverable by inverting `W_sess`; (Sec 4) `g_inf == 2/pi`; (Sec 5) the lowering gap `expand(V_short - V_eff)` equals the symbolic packet sum; (Sec 6) the numeric `Delta = 142.1775`, `D0 ~ 3.76481862`, the rebuilt `V_eff` matches the recorded `1.74701126`, and `lambda_L > 0`. Sections 3 (`DeltaE_UV`) and most of Section 4 (`R`, `M_sigma`, `sigma_min`) are print-only with no assertion.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Susceptibilities = inverse reduced stiffness | asserts chi_** == K_red.inv()[i,j] (py:60-65) | match |
| (1) A_6/A_4/A_2 collection, V_short form | built but never asserted (py:71-80) | partial (constructed, not independently checked; flows into Sec 5/6) |
| (2) S_leak, W_sess scalars | asserts factoring through E0 (py:120-121) | partial (weak; see F2) |
| (2) Lvar recoverability | asserts Lvar_from_W == Lvar (py:122) | mismatch — tautological round-trip (F2) |
| (3) DeltaE_UV = eta_UV D_UV, D_UV >= 0 | print only, no assert (py:130-135) | missing |
| (4) g_inf = 2/pi | asserts g_inf == 2/pi (py:158) | match |
| (4) g(r) >= 2/pi, M_sigma >= 0 on branch | print only, no assert (py:147,155) | missing (F3) |
| (5) lowering identity (gap = packet sum) | asserts expand(gap - expected) == 0 (py:172) | match (but algebraically guaranteed; see F4 note) |
| (5) inequality V_eff <= V_short | not checked | missing (F4) |
| (6) Delta = 210.17750000 (paper) | asserts Delta_num == 142.1775 (py:235) | mismatch — paper value disagrees with script+formula (F1) |
| (6) D0 = 3.76481862 | asserts D0_num ~ 3.76481862 (py:236) | match |
| (6) benchmark decomposition = 1.74701126 | asserts Vrebuild == Veff_obs (py:237) | mismatch — tautological by construction (F5) |
| extra: sigma_min | print only (py:148) | extra (harmless scaffolding) |

`paper_alignment: misaligned` — driven by the F1 numeric Delta contradiction plus several uncovered deliverables.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1-A6 | sympy | 60-65 | `simplify(K_red_inv[i,j] - chi) == 0` | claim 1 (susceptibilities) | yes |
| A7 | sympy | 120 | `simplify(S_leak - S_expected) == 0` | claim 2 | partial (factor-through-E0, weak) |
| A8 | sympy | 121 | `simplify(W_sess - W_expected) == 0` | claim 2 | partial (factor-through-E0, weak) |
| A9 | sympy | 122 | `simplify(Lvar_from_W - Lvar) == 0` | claim 2 | no — tautological inverse round-trip |
| A10 | sympy | 158 | `simplify(g_inf - 2/pi) == 0` | claim 4 (g far-field) | yes |
| A11 | sympy | 172 | `expand(lowering_gap - lowering_expected) == 0` | claim 5 (identity) | partial — true by construction, see F4 |
| A12 | sympy | 235 | `abs(Delta_num - 142.1775) < 1e-8` | claim 6 (benchmark Delta) | yes, but contradicts paper value (F1) |
| A13 | sympy | 236 | `abs(D0_num - 3.76481862) < 1e-7` | claim 6 (benchmark D0) | yes |
| A14 | sympy | 237 | `abs(Vrebuild_soft - Veff_obs) < 1e-10` | claim 6 (decomposition) | no — tautological by construction (F5) |
| A15 | sympy | 238 | `lambda_L_soft > 0` | claim 6 (lambda_L sign) | yes |

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** value_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_247.tex:407-409`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_sympy_audit.md:407`
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage247...sympy_audit.py:235`

**What's wrong:**
The paper card states the benchmark intermediate `\Delta=210.17750000` (stage_247.tex:407; same value in the notes Section 4 block "Delta=210.17750000"). The script computes `Delta = OmU2*OmW2 - Rmix**2` (py:40) and asserts `abs(float(Delta_num) - 142.1775) < 1e-8` (py:235). With the recorded session parameters `OmU2 = 9.0`, `OmW2 = 16.0`, `Rmix = 1.35`, the formula gives `9*16 - 1.35^2 = 144 - 1.8225 = 142.1775` — the script's value. The paper's `210.17750000` is inconsistent with the very formula the paper itself states for `Delta` (`\Delta=\Omega_U^2\Omega_W^2-R_{\rm mix}^{2}`, stage_247.tex:11). Moreover the paper's own next line, `D_0=3.76481862` (stage_247.tex:408), is internally consistent ONLY with `Delta=142.1775`: `D0 = K_* - Q/Delta = 4 - 33.4375/142.1775 = 3.76483` (matches), whereas `Delta=210.1775` would give `D0 = 4 - 33.4375/210.1775 = 3.84091` (does NOT match the paper's stated D0). So the paper is internally self-contradicting and the `210.17750000` figure appears to be a typo; every downstream number in the paper benchmark (D0, V_short=3.74163698, the final 1.74701126) is consistent with the script's 142.1775.

**Why this matters:**
A published stage card carries a numeric value that contradicts both its own stated formula and its own subsequent numbers. This is a paper-side defect, not a script defect — but per the red-team contract a paper↔script numeric disagreement must be routed to the user, not auto-edited.

**Required change:**
None applied by Codex. See `## Resolve before fix_loop` in the directive. The likely resolution is a paper-card typo fix (`210.17750000` → `142.17750000`), but the direction is the user's call.

**Verification:**
After user resolution: if paper is corrected, stage_247.tex:407 reads `\Delta=142.17750000` and matches py:235. No script change.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage247...sympy_audit.py:113,122`

**What's wrong:**
`Lvar_from_W = sqrt(W_sess * pi**4 * lam**2 / (512 * eta_leak**2 * mu_w * q * rho0))` (py:113) is the algebraic inverse of the line that defined `W_sess = 512 * eta_leak**2 * mu_w * q * rho0 * Lvar**2 / (pi**4 * lam**2)` (py:110). Substituting gives `sqrt(Lvar**2) = Lvar` (Lvar is declared `positive`). The subsequent `assert sp.simplify(Lvar_from_W - Lvar) == 0` (py:122) is therefore a round-trip through an inverse: it can never fail and proves nothing about the physics — it only re-confirms that `sqrt` undoes squaring. (The A7/A8 factor-through-`E0` asserts at py:120-121 are weak but non-tautological, since `E0 = 16*eta_leak*Lvar/pi^2` is a distinct upstream intermediate; they are not flagged.)

**Why this matters:**
A tautological assertion masquerades as a verification of the Stage-244 inversion that recovers `Lvar` from the tabulated `W_sess`. The genuine claim that should be tested — that the *numeric* `Lvar(session) = 20.01677473` is the value that reproduces the recorded `W_sess(session) = 1.51632107` — is the one used in Section 6 (py:213-217) but is itself never asserted against a fixed expected number.

**Required change:**
Replace the tautological symbolic assert at py:122 with a substantive numeric check that the inverted `Lvar` reproduces the recorded benchmark `W_sess`: assert that `W_sess` evaluated at `Lvar = Lvar_soft` and the session params equals the recorded `Wsess_obs = 1.51632107` to ~1e-7, and that `Lvar_soft` matches the paper's stated `20.01677473` to ~1e-6. See directive for exact form.

**Verification:**
New numeric asserts appear in Section 6 tying `Lvar_soft` to `Wsess_obs` and to the paper figure 20.01677473; script still exits 0.

### F3 — script_missing_paper_claim

**Severity:** medium
**Subtype:** script_missing_paper_claim (a paper_misalignment subtype, but resolvable script-side with no paper edit)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage247...sympy_audit.py:130-135,147,155`
- paper: notes Section 2.2-2.3 (D_UV >= 0; g(r) >= 2/pi; M_sigma >= 0), appendix item 4 (`\mathcal D_{UV}\ge 0`) and item 5

**What's wrong:**
Three positivity facts that the paper makes load-bearing for the lowering theorem are present in the script only as printed expressions, with no assertion:
- `D_UV >= 0` (Stage-245 drain nonnegativity; appendix `eq:app-part08-main-uv-drain`). Script defines `D_UV` (py:130) and `E_UV` (py:131) but never checks the sign.
- `g(r) >= 2/pi` on the Session-I branch (notes 2.3). Script computes `g_r` (py:143) and `g_inf` (py:144) but only asserts the far-field limit, not the branch lower bound.
- `M_sigma >= 0` on the Session-I compensated branch (notes 2.3 boxed; the source-side half of the lowering inequality). Script computes `M_sigma` (py:147) but never asserts its sign on the branch.

This is not `paper_missing_script_claim` — the paper is correct and the script under-covers it. Resolvable purely script-side, no paper edit.

**Why this matters:**
The stage's headline result is the lowering inequality `V_eff <= V_short`, which holds *only because each imported packet scalar is nonnegative on the branch*. The script verifies the algebraic identity but not the nonnegativity that gives it physical content. `M_sigma >= 0` is the nontrivial one (it depends on `g(r) < rF1` and `g(r) >= 2/pi`); leaving it unchecked means the script could pass even if the source term had the wrong sign.

**Required change:**
Add asserts: (a) `D_UV >= 0` symbolically given the positive-symbol assumptions (it is `a_V*chi_UV**2*f_U**2/Delta_UV**2`, manifestly nonneg for real `Delta_UV != 0` and `a_V>0`); (b) a numeric check on the Session-I slice that `g(r_soft) >= 2/pi` and `g(r_soft) < rF1`; (c) a numeric check `M_sigma(r_soft) >= 0` (the printed value is 0.18386120 > 0). See directive for forms and the variable-independence self-test note.

**Verification:**
New asserts for `D_UV >= 0`, `g(r_soft)` bounds, and `M_sigma(r_soft) >= 0` appear; script exits 0.

### F4 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage247...sympy_audit.py:166-172`

**What's wrong:**
The lowering-identity assert `assert sp.expand(lowering_gap - lowering_expected) == 0` (py:172) is algebraically guaranteed: `V_eff` is *defined* as `V_short - lambda_L*S_leak - lambda_W*W_sess - E_UV - M_sigma` (py:166), so `lowering_gap = expand(V_short - V_eff)` is identically `lambda_L*S_leak + lambda_W*W_sess + E_UV + M_sigma`, which is exactly `lowering_expected` (py:168). The assert restates the definition of `V_eff`; it cannot fail. It documents the structure but does not exercise the *physical* deliverable of claim (5), which is the inequality `V_eff <= V_short` (appendix item 6, "Hence ... V_eff <= V_short pointwise"). That inequality is never checked.

**Why this matters:**
The script gives the appearance of verifying the lowering theorem while only re-expanding a definition. The substantive content (sign of the gap) is the part that is missing, and overlaps with F3 (the packet-scalar nonnegativities that make the gap nonnegative).

**Required change:**
Keep the identity assert (it is fine as a structural sanity check) but add the inequality content: assert that the symbolic `lowering_gap` is a sum of terms each known nonnegative under the stated assumptions (lambda_L, lambda_W nonnegative; S_leak, W_sess positive; E_UV nonneg via F3; M_sigma nonneg via F3 on the branch). Practically: after F3 lands, add a numeric Session-I assert that `V_eff(r_soft) <= V_short(r_soft)` and that the full numeric `lowering_gap(session) >= 0`. See directive.

**Verification:**
New numeric assert `Vshort_num - Veff_session >= 0` (equivalently `lowering_gap >= 0`) appears on the benchmark slice; script exits 0.

### F5 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/...stage247...sympy_audit.py:219-220,237`

**What's wrong:**
`lambda_L_soft = (Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs) / S_soft` (py:219) solves for `lambda_L` so that the rebuilt potential equals `Veff_obs`. Then `Vrebuild_soft = Vshort_num - lambda_L_soft*S_soft - Wsess_obs - UVdrop_obs - M_sigma_num` (py:220). Substituting the line-219 definition of `lambda_L_soft` makes the `lambda_L_soft*S_soft` term cancel to exactly `(Vshort_num - Wsess_obs - UVdrop_obs - M_sigma_num - Veff_obs)`, leaving `Vrebuild_soft = Veff_obs` identically. Therefore `assert abs(float(Vrebuild_soft - Veff_obs)) < 1e-10` (py:237) is true by construction for ANY input numbers — it is a self-solved-root resubstitution and cannot fail. The genuine benchmark claim (notes 4.2) is that the *independently derived* numbers `Vshort = 3.74163698`, `Wsess = 1.51632107`, `UVdrop = 0.21064278`, `M_sigma = 0.18386120`, `S_leak = 0.31069599`, and the *paper-stated* `lambda_L = 0.26971918` together reproduce `1.74701126` — i.e. `lambda_L` should be checked against its paper value, not solved away.

**Why this matters:**
The headline benchmark decomposition (eq:app-part08-stage247-benchmark-decomposition) is the stage's most concrete result, and the script's check of it cannot fail no matter what the inputs are. The check is vacuous: it verifies that solving a linear equation and back-substituting recovers the right-hand side. The real assertions that would catch an error — `Vshort_num ~ 3.74163698`, `M_sigma_num ~ 0.18386120`, `S_soft ~ 0.31069599`, and `lambda_L_soft ~ 0.26971918` (the paper value) — are computed but only `Delta`/`D0` (py:235-236) and `lambda_L > 0` (py:238) are asserted.

**Required change:**
Add substantive numeric asserts pinning each independently derived benchmark quantity to its paper-stated value: `Vshort_num ~ 3.74163698`, `M_sigma_num ~ 0.18386120`, `S_soft ~ 0.31069599`, `Lvar_soft ~ 20.01677473`, and `lambda_L_soft ~ 0.26971918`. Then keep the closure check by asserting the *forward* decomposition with the paper's `lambda_L` literal (not the solved one): `abs((Vshort_num - 0.26971918*S_soft - Wsess_obs - UVdrop_obs - M_sigma_num) - Veff_obs) < 1e-6`. That makes the closure falsifiable. See directive.

**Verification:**
New asserts pin V_short, M_sigma, S_leak, Lvar, and lambda_L to their paper figures; the closure assert uses the paper's literal lambda_L = 0.26971918 rather than the solved value; script exits 0.

### F6 — missing_verification_script

**Severity:** high
**Subtype:** missing_mathematica
**Files:**
- target: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage247_relaxed_stationary_barrier_compiler_from_one_port_short_range_leakage_uv_and_compensated_source_packets_mathematica_audit.wl` (does not exist)

**What's wrong:**
Stage 247 is checkpoint:`False`, `is_status_only_candidate:False` → both engines required. No `.wl` exists (`mathematica/moving_throat_pde_stage247_*.wl` not found). Every claim in this stage is independently verifiable in native Mathematica: the susceptibilities are entries of `Inverse[K_red]` (Mathematica `Inverse` + `Simplify`); the leakage/work factoring is `Simplify`; the source-response `g_inf` is `Limit[g[r], r -> Infinity]`; `M_sigma >= 0` / `g >= 2/pi` are `Reduce`/`Simplify[..., assumptions]` or numeric; the lowering identity is `Simplify[Vshort - Veff - packetSum]`; and the Session-I benchmark is straight `N[...]` arithmetic. There is no genuine single-engine obstruction.

**Why this matters:**
The dual-engine policy requires a second, independently decomposed engine wherever one is possible. Here it is clearly possible, so the absence is a finding.

**Required change:**
Create the `.wl` as a native, independently structured re-derivation (NOT a transliteration of the `.py`): build `K_red` and take `Inverse`, extract susceptibility entries, and confirm against the paper's closed forms via `Simplify`; derive the leakage/work scalars and confirm the factoring; take the `Limit` for `g_inf`; verify `D_UV >= 0`, `g(r_soft) >= 2/pi`, `M_sigma(r_soft) >= 0`; verify the lowering identity by `Simplify[lowering_gap - packetSum] == 0`; and reproduce the Session-I benchmark numerically with the paper's literal `lambda_L = 0.26971918` (forward direction, not solved). Use `expectZero`/`expectTrue`/`expectApprox` helpers. See directive claim manifest M1-M8.

**Verification:**
`redteam exec-mathematica 247` runs the new `.wl`, all `expect*` checks pass, exit 0; the susceptibility/identity routes use Mathematica primitives (`Inverse`, `Limit`, `Reduce`/`Simplify`) rather than echoing the SymPy choreography.

## Independent-derivation check (Mathematica)

No `.wl` exists, so no transliteration to assess. `mathematica_transliteration` is NOT flagged (nothing to compare). The directive prescribes an independent route to avoid transliteration when the file is created.

## Engine cross-check

Only SymPy is present; `engines_agree: n/a`. Cross-check deferred to after the `.wl` is created (F6).

## Verdict justification

`verdict: findings`, `stop_cold: null`. The susceptibility-vs-inverse checks (A1-A6) and the `g_inf = 2/pi` limit (A10) are genuine and hold up. But the stage has one paper-side numeric contradiction (F1: paper `Delta=210.17750000` vs the formula-correct `142.1775`, and the paper's own D0 confirms 142.1775 — a paper typo, routed to the user), two vacuous-by-construction assertions (F2 `Lvar` inverse round-trip; F5 the benchmark `lambda_L` self-solve/resubstitution — the headline decomposition check literally cannot fail), missing positivity checks for the packet scalars the lowering theorem rests on (F3), an identity assert that only re-expands a definition while the actual inequality goes unchecked (F4), and a missing required second engine (F6 missing_mathematica). Attacks that succeeded: the `Lvar_from_W` and `Vrebuild` round-trips collapse to `sqrt(x^2)=x` and `solve-then-resubstitute` respectively. Attacks that failed: the matrix-inverse susceptibility asserts are real and non-tautological; the `g_inf` limit is real; the symbolic lowering identity, while definitional, is at least correctly stated. I confirm I read the paper card, the notes, and the part-08 appendix rows before opening the script.

## Self-test notes

Variable-independence trap: the prescribed new asserts use no `sp.diff`/`D` against a variable the expression lacks; the only derivative-free positivity/limit checks are over variables the expressions genuinely contain (`g_r`, `M_sigma`, `D_UV` all depend on their checked variables), so no identically-zero vacuity is introduced. Trivial-case pre-check: `M_sigma(r_soft)` evaluates to the printed `0.18386120 > 0`, `D_UV` is a square-over-square times `a_V>0` (nonneg), and the forward benchmark with the paper's literal `lambda_L=0.26971918` reproduces `1.74701126` (hand-checked: `0.26971918*0.31069599 = 0.0838007`, `3.74163698 - 0.0838007 - 1.51632107 - 0.21064278 - 0.18386120 = 1.7470113`), so the proposed forward closure assert is satisfiable and non-vacuous. Paper round-trip: F2/F3/F4/F5 fixes pin to the paper-stated figures (20.01677473, 0.18386120, 0.31069599, 0.26971918, 3.74163698) and do not introduce any constant absent from the card/notes; F1 is left for the user and not auto-fixed.
