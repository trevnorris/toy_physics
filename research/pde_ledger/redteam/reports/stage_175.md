---
unit_id: 175
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-28T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage175_wall_normalized_load_shape.md]
  paper_appendix: present
---

# Audit unit 175 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_175.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage175_wall_normalized_load_shape.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 81, 611-616)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` states verbatim: "Factors the remaining defect into wall-normalized shape/load variables and proves conservative-shape preservation leaves only outgoing load slippage." The notes supply the detail: the three Stage-242 static bundles factor exactly by the wall baseline `K` as `B0 = K chi^2`, `Z0 = K Upsilon`, `N0 = Lambda^2` (homogeneity laws `Delta_r = K^2 hat_Delta_r`, `Q_r = K^3 hat_Q_r`, `P_r = K^2 hat_P_r`), so the three self-similarity defects become `Sigma_B = d ln chi^2`, `Sigma_Z = d ln Upsilon`, and `Sigma_N = 2 d ln Lambda - delta_K = d ln(Lambda^2/K)`. The stage's two headline theorems are: (a) on conservative-shape-preserving branches (`Sigma_B = Sigma_Z = 0`) the remaining defect is carried only by the outgoing load factor, `Xi_load = sum_r rho_r^(N)(2 d ln Lambda_r - delta_K)`; and (b) the no-go: if all wall-normalized shapes are frozen, `Xi_load = -delta_K` (using `sum_r rho_r^(N) = 1`), so naive common self-similarity does not kill the defect unless `delta_K = 0`. Distinct deliverables: D1 factorization (B0/Z0/N0 + the three K-homogeneity laws), D2 the three differential defect identities, D3 the conservative-shape reduction, D4 the no-go `Xi_load = -delta_K`.

## What the script claims to verify

Both engines verify (1) the exact homogeneity laws `Delta = K^2 Delta_hat`, `Q = K^3 Q_hat`, `P = K^2 P_hat`, and the bundle factorizations `B0 = K chi^2`, `Z0 = K Upsilon`, `N0 = Lambda^2`; (2) the differential identities `Sigma_B = d ln chi^2`, `Sigma_Z = d ln Upsilon`, `Sigma_N = d ln(Lambda^2/K)` and `Sigma_N = 2 d ln Lambda - delta_K`, computed through a one-parameter `eps`-flow with log-slopes `kappa` (=`delta_K`), `schi`, `su,sw,sr,sgu,sgw`; (3) the conservative-shape branch (`schi=0`, all port slopes `=0`) gives `Sigma_B = Sigma_Z = 0` and `Sigma_N = -kappa`. The docstring (lines 16-18) further claims items 3-4 as `Xi_load = <2 d ln Lambda - dK>_N` and `Xi_load = -dK`, i.e. the rho-weighted aggregate forms. The verdict applies to a single port `r` (no summation, no `rho` weights).

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 homogeneity `Delta=K^2 hat`, `Q=K^3 hat`, `P=K^2 hat` | sympy 67-78 / math 51-56 | match |
| D1 factorization `B0=K chi^2`, `Z0=K Upsilon`, `N0=Lambda^2` | sympy 50,82,83 / math 38,60,61 | match |
| D2 `Sigma_B = d ln chi^2` | sympy 118 / math 83 | match |
| D2 `Sigma_Z = d ln Upsilon` | sympy 125 / math 89 | match |
| D2 `Sigma_N = d ln(Lambda^2/K)` | sympy 132 / math 95 | match (weak/log-algebra; substantive K-homogeneity sits in D1 line 83) |
| D2 `Sigma_N = 2 d ln Lambda - delta_K` | sympy 133 / math 96 | mismatch (tautological — see F1) |
| D3 conservative branch `Sigma_B=Sigma_Z=0` | sympy 145,146 / math 102,103 | match |
| D4 no-go `Xi_load = -delta_K` (with `sum rho^(N)=1`) | sympy 147 / math 104 (per-port `Sigma_N=-kappa` only) | partial — no `rho`-weighting / `sum rho=1` step (see F2) |

The card's bottom-line `\stagefield{Output}` (factorization + "conservative-shape preservation leaves only outgoing load slippage") is faithfully exercised. `paper_alignment: aligned`. The two findings are script-internal quality issues (a tautological assertion and a missing aggregation step relative to the script's own docstring), not a misstatement of the paper's claim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `B0 - K*chi^2 == 0` | D1 | yes |
| A2 | sympy | 67-70 | `Delta - K^2*Delta_hat == 0` | D1 | yes |
| A3 | sympy | 71-74 | `Q - K^3*Q_hat == 0` | D1 | yes |
| A4 | sympy | 75-78 | `P - K^2*P_hat == 0` | D1 | yes |
| A5 | sympy | 82 | `Z0 - K*Upsilon == 0` | D1 | yes |
| A6 | sympy | 83 | `N0 - Lambda^2 == 0` | D1 | yes |
| A7 | sympy | 118 | `Sigma_B_direct - dln(chi^2) == 0` | D2 | yes |
| A8 | sympy | 125 | `Sigma_Z_direct - dln(Upsilon) == 0` | D2 | yes |
| A9 | sympy | 132 | `Sigma_N_direct - dln(Lambda^2/K) == 0` | D2 | partial (log-algebra; K-homogeneity already at A6) |
| A10 | sympy | 133 | `Sigma_N_direct - (2 dln Lambda - kappa) == 0` | D2 | **no (tautological, F1)** |
| A11 | sympy | 145 | `Sigma_B_cons == 0` | D3 | yes |
| A12 | sympy | 146 | `Sigma_Z_cons == 0` | D3 | yes |
| A13 | sympy | 147 | `Sigma_N_common + kappa == 0` | D4 (per-port) | partial (no rho-weighting, F2) |
| B1-B13 | mathematica | 38,51,56,60,61,83,89,95,96,102,103,104 | identical forms | same as A1-A13 | identical (transliteration, F3) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:127-133`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:91-96`

**What's wrong:**
The assertion `expect_zero("Sigma_N - (2 dln Lambda - dK)", Sigma_N_direct - (2 * dlog(expr_L) - kappa))` (sympy line 133) cannot fail by construction. `Sigma_N_direct` is defined at line 130 as `simplify(2*dlog(expr_ratio) - kappa)` with `expr_ratio = (P / Delta).subs(subs_hat).subs(subs_eps)` (line 128). The RHS uses `expr_L = Lambda.subs(subs_eps)` (line 129) with `Lambda = simplify((P / Delta).subs(subs_hat))` (line 81). Since `simplify` is value-preserving, `expr_L` and `expr_ratio` are the *same expression* — both are `(P/Delta)` after the same `subs_hat` then `subs_eps`. The assertion is therefore `[2 dlog(expr_ratio) - kappa] - [2 dlog(expr_ratio) - kappa] == 0`, i.e. `X - X == 0`. It holds independent of the actual forms of `P`, `Delta`, or any physics. The Mathematica line 96 (`expectZero["Sigma_N - (2 dln Lambda - dK)", sigmaNDirect - (2*dlog[exprL] - kappa)]`) has the identical defect (`exprRatio` line 91 vs `exprL` line 92 are the same object). Note also that line 132 (`Sigma_N - dln(Lambda^2/K)`) is only weakly non-tautological: because `(P/Delta).subs(subs_hat)` cancels its `K^2` and is degree-0 in `K`, that check reduces to the pure log-algebra identity `dlog(X^2/K) = 2 dlog(X) - kappa`; the substantive scale-invariance content of `Sigma_N` lives in the homogeneity check A6 (`N0 - Lambda^2 == 0`, line 83), which is genuine.

**Why this matters:**
Line 133/96 contributes a guaranteed PASS to the transcript, giving false confidence that the identity `Sigma_N = 2 d ln Lambda - delta_K` was independently exercised. The paper's `Sigma_N` claim (notes line 223) IS covered by the genuine homogeneity check A6 plus the weak A9, so this is not a correctness gap in the audit's conclusion — but the assertion as written is dead weight that cannot detect a regression.

**Required change:**
Make line 133/96 exercise the paper's `Sigma_N = 2 d ln Lambda_r - delta_K` against an *independent* computation of `2 d ln(P_r/Delta_r) - delta_K` that is not literally re-using `Lambda`. Compute the left ratio directly from the un-substituted physical primitives perturbed in `eps` (i.e. perturb `Ou, Ow, R, gu, gw` and `K` via the physical scaling `subs_hat ∘ subs_eps` so the `K` power is visibly present before cancellation), so that the `-delta_K` subtraction is actually load-bearing. Concretely: build `expr_PD_phys = (P/Delta)` then apply `subs_hat` then `subs_eps` WITHOUT a prior `simplify` to `Lambda`, and compare `2*dlog(expr_PD_phys) - kappa` against `2*dlog(Lambda.subs(subs_eps)) - kappa` only if a distinct route is used — otherwise replace line 133 with a check that the `kappa` term genuinely vanishes, e.g. `expect_zero("Sigma_N has no residual kappa beyond -dK", (Sigma_N_direct - (2*dlog(expr_L))) + kappa)` rephrased so it does not collapse to `X-X`. See directive for the exact safe replacement.

**Verification:**
After the fix, sympy line 133 (and math line 96) must no longer be of the form `Sigma_N_direct - (2*dlog(<same expr as in Sigma_N_direct>) - kappa)`. The verifier confirms the new assertion references a differently-constructed expression (the perturbed physical primitives rather than the cached `Lambda`/`expr_ratio`) and that both scripts still exit 0 with the line printing `= 0`.

### F2 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:16-18, 144-147`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:101-104`

**What's wrong:**
The script docstring (lines 16-18) states the deliverables as "3. Conservative-shape-preserving reduction `Xi_load = <2 d ln Lambda - dK>_N`. 4. Naive common-shape branch gives `Xi_load = -dK`," where `<...>_N` denotes the `rho_r^(N)`-weighted port average and the no-go `Xi_load = -dK` depends on `sum_r rho_r^(N) = 1` (notes line 330). The assertions only verify the *per-port* quantity `Sigma_N_common + kappa == 0` (line 147 / math 104), i.e. `Sigma_N = -kappa` for a single port. No `rho_r^(N)` weight, no summation, and no `sum_r rho_r^(N) = 1` step ever appears. The headline no-go result that the notes box (`Xi_load = -delta_K`, notes line 333) is therefore not assembled from the per-port building block; the closing aggregation `sum_r rho_r^(N)(-delta_K) = -delta_K` is asserted only in the printed "Conclusions" text, not checked.

**Why this matters:**
The script's own stated claim 4 (`Xi_load = -dK`) is stronger than what is verified (`Sigma_N = -dK` per port). The missing piece (`sum rho = 1` weighting) is mathematically trivial, but leaving it as printed prose rather than an assertion means a future change to the weighting convention would not be caught by this audit. This is the stage's sharpest theorem, so it deserves a real check.

**Required change:**
Add an assertion that exercises the weighted aggregation explicitly. Introduce a small fixed set of port weights `rho` (e.g. two symbolic non-negative weights `rho1, rho2` with the constraint `rho1 + rho2 = 1`, or a symbolic list summing to 1) and confirm `sum_r rho_r * Sigma_N_common` simplifies to `-kappa` under `sum rho = 1`. Equivalently, verify the identity `Xi_load_frozen = (rho1 + rho2)*(-kappa)` and substitute `rho1 + rho2 -> 1` to get `-kappa`, asserting `expect_zero("Xi_load (all shapes frozen) + dK", Xi_load_frozen.subs({rho1+rho2: 1}) + kappa)`. This anchors the `sum rho = 1` step that the no-go theorem relies on.

**Verification:**
After the fix, a new assertion referencing `rho` weights and `sum rho = 1` appears (sympy after line 147, math after line 104), prints `= 0`, and both scripts exit 0. The directive must keep the weight constraint as `sum rho = 1` to match notes line 330.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:36-104`
- (compare) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:48-147`

**What's wrong:**
The `.wl` script is a mechanical line-by-line port of the `.py` script, not an independent re-derivation. Corresponding excerpts:
- Bundle defs — sympy 55-57: `Delta = Ou * Ow - R**2` / `Q = gu**2 * Ow + 2 * gu * gw * R + gw**2 * Ou` / `P = Ou * gw + R * gu`; math 45-47: `delta = ou*ow - r^2` / `q = gu^2*ow + 2*gu*gw*r + gw^2*ou` / `p = ou*gw + r*gu`. Identical term order and structure.
- Substitutions — sympy `subs_hat` (59-65) and `subs_eps` (103-112) are reproduced verbatim as `subsHat` (49) and `subsEps` (69-78), same variable mapping, same `Exp[slope*eps]` choreography.
- Differential block — sympy 128-133 and math 91-96 build `expr_ratio`/`exprRatio`, `expr_L`/`exprL`, `Sigma_N_direct`/`sigmaNDirect` in the same order with the same `2*dlog(...) - kappa` forms (and therefore inherit the identical F1 tautology).
- Both even reproduce the same wrong banner string `"STAGE 158 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION"` (sympy 39, math 31), confirming copy-paste origin.

The second-engine policy requires the Mathematica audit to derive the result independently from the physical premises so that a transcription error in one engine cannot pass undetected in both. Here both engines share the same algebra and the same F1 defect, so they cannot cross-check each other.

**Why this matters:**
An identical bug (such as F1) is present in both engines and passes in both, defeating the purpose of a dual-engine cross-check. The cosmetic "STAGE 158" banner also misleads a reader of the transcript about which stage is being audited.

**Required change:**
Re-derive the Mathematica checks from the physical premises along a structurally different route rather than echoing the SymPy variable choreography. At minimum: (a) build the homogeneity laws by an explicit degree count (e.g. confirm `delta /. {ou->t^2*ou,...}` scales as `t^4` — using the physical interpretation that `ou=Omega_U^2` etc.) rather than the `subsHat` substitution copied from sympy; (b) compute the differential identities via `Series[..., {eps,0,1}]` coefficient extraction instead of the copied `dlog` helper; (c) fix the banner to `"STAGE 175 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION"`. If a full re-derivation is out of scope for the fix loop, at minimum apply (c) the banner correction and the F1 fix independently in each engine.

**Verification:**
The verifier confirms the `.wl` no longer mirrors the `.py` line-by-line for the differential block (different intermediate construction), the banner reads `STAGE 175`, and the script exits 0.

## Independent-derivation check (Mathematica)

Not independent. The `.wl` is a transliteration of the `.py`: same `delta = ou*ow - r^2`, `q`, `p` definitions in identical term order (math 45-47 ↔ sympy 55-57), identical `subsHat`/`subsEps` maps, identical `dlog` helper (math 26-29 ↔ sympy 100-101), identical assertion order and identical `Sigma_N` tautology (math 91-96 ↔ sympy 128-133), and even the identical copy-paste banner typo "STAGE 158". See F3.

## Engine cross-check

Both transcripts agree at the level claimed: every line prints `... = 0` / `PASS`, both exit 0. The sympy output (12:47) and mathematica output (13:21) are both newer than the script files (11:56), so outputs are fresh. Because the engines are transliterations (F3), this agreement is weak evidence — it confirms the shared algebra is internally consistent but cannot independently corroborate it. No `engine_disagreement`.

## Verdict justification

The factorization and homogeneity content (D1: `Delta=K^2`, `Q=K^3`, `P=K^2`, `B0=K chi^2`, `Z0=K Upsilon`, `N0=Lambda^2`) is genuine, non-tautological, and matches the notes exactly; I attacked the variable-renaming (script symbol `Ou` = paper `Omega_U^2`, scaled `Ou=K*Ou_hat` ⇔ paper `Omega_U^2 = K*hat_Omega_U^2`) and found the scaling self-consistent and the K-degrees (2,3,2 ⇒ Z0~K, N0~K^0) correct. The conservative-branch `Sigma_B=Sigma_Z=0` and the no-go `Sigma_N=-kappa` (D3, D4 per-port) are genuine and anchor the stage's headline theorem. The over-strong `positive=True` assumptions are harmless (identities are sign-independent; log-derivatives are branch-independent). The verdict is `findings` for three script-quality issues: F1 a clean tautology at sympy 133 / math 96 (`expr_L === expr_ratio`), F2 the docstring's `Xi_load = -dK` claim verified only per-port without the `sum rho=1` aggregation, and F3 the Mathematica script being a transliteration of the SymPy script (sharing the F1 defect and a wrong banner). None rise to stop-cold: the paper's bottom-line Output is faithfully covered (`paper_alignment: aligned`), and all three fixes are local and non-propagating.

## Self-test notes

Variable-independence: checked that `expr_ratio`/`exprL` and `expr_Z`/`exprU` are the load-bearing pairs — `expr_Z = K*expr_U` (differ by a K, so Sigma_Z is real) but `expr_ratio = expr_L` (identical, so Sigma_N line 133 is the tautology, F1). Trivial-case: for F2's proposed fix, with `rho1+rho2=1` and per-port `Sigma_N=-kappa`, the weighted sum is `(rho1+rho2)(-kappa) = -kappa`, residual `+kappa = 0` confirmed. Symbol domains: `positive=True` over-strong but harmless since every check is a polynomial/log-derivative identity independent of sign and branch. Paper round-trip: the F1/F2 fixes reuse only constants already in the notes (`delta_K=kappa`, `sum rho=1` from notes line 330) and introduce no constant the paper lacks.
