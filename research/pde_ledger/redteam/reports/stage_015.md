---
unit_id: 015
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 015 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_015.tex`
- notes: `(none)` (no `notes/stages/moving_throat_pde_stage015_*.md` files exist)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 52 + `\input{stages/stage_015}` at line 109)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.txt`

## What the paper claims

Stage 015 promotes the throat to an autonomous graph field `R(Omega, w, t)` and posits the parent action `S_Sigma[R] = int dt dw dOmega L_Sigma` with density `L_Sigma = (1/2) mu_Sigma R_t^2 - (1/2) T_{w,Sigma} R_w^2 - (1/2) T_{Omega,Sigma} |grad_Omega R|^2 - U_Sigma`. The total parent-level model is `S_total = S_psi[psi, A; Sigma] + S_EM[A] + S_Sigma[R]`. The stage's `\stagefield{Output}` is quoted verbatim: "Stage~015 exports the parent throat action promotion and the exact quadratic recovery formula \eqref{eq:stage015-keta}", i.e., `K_eta = U_{Sigma,RR}(R_0,w) - partial_w(T_{w,Sigma,R}(R_0,w) R_0') + (1/2) T_{w,Sigma,RR}(R_0,w) (R_0')^2`. The card also asserts that "The audit carries the boundary densities explicitly and tests both zero and nonzero boundary-discharge probes." The part-appendix row for stage 015 calls it the "parent throat-action packet and projection/reduction status boundary" with status `\StatusExact{}/\StatusReduced{}`. No notes file exists for the stage.

## What the script claims to verify

The sympy script tests: (a) a generic quadratic IBP product-rule identity and a concrete Gaussian instantiation, plus a `boundary_value`-operator sanity check on `atan(w)`; (b) the K_eta canonical-quadratic-form recovery by symbolically expanding `L` to O(eps^2), peeling the eta-eta_w cross term, and matching the IBPed result against the textbook K_eta plus a sign-mutation guard; (c) a wall-only specialization of full even-channel gates K1, H_even (with B/Z modes set to zero), Jacobian determinant `1/27`, perturbed-gate solves, a Gaussian-overlap closed-form check with a `5*delta_TO` vs `6*delta_TO` mutation guard; (d) three real-Y20 overlap ratios (`{1, 1/2, -1}`) plus same-sign-cross-term vanishing; (e) a grouped trace/anomaly identity `xbar = x0` and `bx = 3*ax`. The mathematica script tests the same nine M1-M9 blocks in the same order, with minor extras (closed-form values `dMoverlap = sqrt(pi/3)`, `dKoverlap = 23*sqrt(pi)/(3*sqrt(3))`, `wallDetShift = 2*eps/3`).

## Paper <-> script cross-check

| Paper-side deliverable | Script-side coverage | Status |
|---|---|---|
| Lagrangian density form `L_Sigma` (eq:stage015-parent-density) | Built into `lagrangian = mu0*Rt^2/2 - Tw*Rw^2/2 - TO0*eps^2*grad2/2 - U` (used as input to K_eta derivation; not separately asserted) | partial |
| Total action `S_total = S_psi + S_EM + S_Sigma` (eq:stage015-total-action) | No script-side check (additive declaration only - not algebraically testable in isolation) | missing |
| K_eta exact quadratic recovery formula (eq:stage015-keta) | sympy line 88-97 / mathematica M3 (`L2afterIBP` vs `canonicalL2`), plus dTwRR0p sign mutation | match |
| Boundary densities + zero-discharge probe | sympy line 59-72 / mathematica M2 (Gaussian profiles) | partial (degenerate - see F2) |
| Boundary densities + nonzero-discharge probe | sympy line 57-58 (atan probe); mathematica omits this entirely | partial |
| Wall-only K1/H_even gates, Jacobian det 1/27, Gaussian overlaps | sympy lines 123-193 / mathematica M4-M7 | extra (no paper-side counterpart - see F1) |
| Y20 overlap ratios `{1, 1/2, -1}` | sympy lines 195-200 / mathematica M8 | extra (see F1) |
| Grouped trace/line `xbar=x0, bx=3*ax` | sympy lines 202-208 / mathematica M9 | extra (see F1) |

`paper_alignment: partial`. The central K_eta deliverable is matched. The boundary-density story is partly matched but with a degenerate concrete profile that does not actually exercise the asymmetric-discharge case. The Lagrangian-form and total-action structural claims are declarative (carried by construction, not independently asserted). A large block of the script (M4-M9, six asserts in sympy plus their mathematica mirrors) verifies wall-only gates, Y20 ratios, and grouped-trace identities that the stage card and appendix row do not mention.

## Assertion inventory

| #  | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|----|--------|------|------|-------------------------------|---------------------|
| A1 | sympy | 47-50 | `assert_zero(generic quadratic IBP identity, ...)` | "boundary densities" (general IBP) | yes |
| A2 | sympy | 51-54 | `assert_nonzero(mutated quadratic IBP sign should fail)` | guard for A1 | yes |
| A3 | sympy | 58 | `assert_nonzero(boundary operator detects nonzero endpoint discharge)` (atan probe) | "nonzero boundary-discharge probe" (literal but decorative) | partial |
| A4 | sympy | 68 | `assert_zero(concrete quadratic IBP boundary discharge)` (Gaussian) | "zero boundary-discharge probe" | partial (degenerate, see F2) |
| A5 | sympy | 69-72 | `assert_zero(concrete quadratic IBP with boundary)` | IBP identity at Gaussian | partial (degenerate, see F2) |
| A6 | sympy | 93 | `assert_zero(raw eta eta_w cross term)` | K_eta derivation step | yes |
| A7 | sympy | 97 | `assert_zero(quadratic wall action after integration by parts)` | K_eta exact formula | yes |
| A8 | sympy | 101 | `assert_nonzero(mutated K_eta sign should fail)` | guard for A7 | yes |
| A9 | sympy | 126 | `assert_zero(wall-only K1 specialization)` | none (extra) | tautological - see F3 |
| A10 | sympy | 127 | `assert_zero(wall-only H_even specialization)` | none (extra) | tautological - see F3 |
| A11 | sympy | 146-150 | `assert_zero(wall-only K1 from Gaussian overlaps)` | none (extra) | yes (within extra scope) |
| A12 | sympy | 151-155 | `assert_zero(wall-only H_even from Gaussian overlaps)` | none (extra) | yes (within extra scope) |
| A13 | sympy | 163-167 | `assert_nonzero(wall-only K1 detects mutated 6*delta_TO)` | none (extra) | yes (within extra scope) |
| A14 | sympy | 174 | `assert_zero(wall-only even-gate determinant - 1/27)` | none (extra) | yes (within extra scope) |
| A15 | sympy | 176-177 | `assert_zero(wall-only dK closed form)`, `assert_zero(... dM ...)` | none (extra) | yes (within extra scope) |
| A16 | sympy | 179-180 | `assert_nonzero(wall-only solve detects perturbed K1 gate ...)` | none (extra) | yes (within extra scope) |
| A17 | sympy | 193 | `assert_nonzero(wall-only determinant detects perturbed K1 coefficient)` | none (extra) | yes (within extra scope) |
| A18 | sympy | 198-200 | `assert_zero(Y20 overlap lane 20/21/22)` | none (extra) | yes (within extra scope) |
| A19 | sympy | 207-208 | `assert_zero(grouped trace, grouped line b=3a)` | none (extra) | yes (within extra scope) |
| M1 | math | 60 | `expectZero[M1 generic IBP product-rule identity]` | mirrors A1 | yes |
| M2 | math | 65 | `expectNonzero[M1 mutated IBP boundary sign]` | mirrors A2 | yes |
| M3 | math | 75 | `expectZero[M2 Gaussian IBP boundary discharge]` | mirrors A4 | partial |
| M4 | math | 81 | `expectZero[M2 Gaussian IBP cross equals bulk]` | mirrors A5 | partial |
| M5 | math | 90 | `expectZero[M3 K_eta raw eta etaw cross coefficient]` | mirrors A6 | yes |
| M6 | math | 96 | `expectZero[M3 K_eta canonical quadratic form]` | mirrors A7 | yes |
| M7 | math | 102 | `expectNonzero[M3 K_eta dTwRR0p sign mutation]` | mirrors A8 | yes |
| M8-M9 | math | 112-113 | wall-only K1/H_even specializations | mirrors A9-A10 | tautological - see F3 |
| M10-M11 | math | 126-127 | `expectZero[M5 Gaussian dM/dK overlap closed form]` | none (extra, but useful checkpoint of integrals) | yes (within extra scope) |
| M12-M16 | math | 131-148 | K1wallNum / HevenwallNum / mutated guard | mirrors A11-A13 | yes (within extra scope) |
| M17 | math | 151 | `expectZero[M6 wall-only Jacobian determinant - 1/27]` | mirrors A14 | yes (within extra scope) |
| M18-M19 | math | 164-165 | wall det perturbation = `2*eps/3`, nonzero | extra to sympy A17 (closed-form check) | yes |
| M20-M22 | math | 168-171 | wall solve + perturbed solve | mirrors A15-A16 | yes (within extra scope) |
| M23-M25 | math | 175-177 | Y20 ratios | mirrors A18 | yes (within extra scope) |
| M26 | math | 178-184 | `Do[expectZero[same-sign cross term m=...]]` | mirrors `real_y20_square_ratio` interior `same_sign` guard | yes (within extra scope) |
| M27-M28 | math | 195-196 | grouped trace / line | mirrors A19 | yes (within extra scope) |

The mathematica script does not include an analogue of A3 (the `atan` nonzero-discharge sanity probe). Otherwise, the mathematica spine is a structurally identical port of the sympy spine, with two extras (closed-form `sqrt(pi/3)`, `23 sqrt(pi)/(3 sqrt(3))` for the overlap integrals, and `2*eps/3` for the determinant perturbation).

## Findings

### F1 - paper_misalignment

**Severity:** high
**Subtype:** paper_missing_script_claim
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_015.tex:44-46` (the entire `\paragraph{Output.}` block)
- `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex:52` (the appendix row for stage 015)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:103-208` (wall-only gates, Y20 ratios, grouped trace)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-196` (mirroring blocks M4-M9)

**What's wrong:**
The paper card's `\stagefield{Output}` says only: "Stage~015 exports the parent throat action promotion and the exact quadratic recovery formula \eqref{eq:stage015-keta}." The card body discusses promotion to `R(Omega,w,t)`, the Lagrangian density `L_Sigma`, the total action sum, and the K_eta formula - and nothing else. The appendix row reinforces this with "Parent throat-action packet and projection/reduction status boundary." Yet the scripts dedicate ~6 of 12 sympy asserts (A9-A19) and ~21 of 28 mathematica checks (M8-M28) to:
- wall-only specializations of full even-channel gates K1 = -dM + dK/9 and H_even = (2/3)dM - dK/27, the `wall_only_specialization` substitution `{B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}` (sympy:123 / math:109), the Jacobian determinant `1/27` (sympy:174 / math:151), perturbed-solve diagnostics (sympy:176-180 / math:168-171), Gaussian overlap integrals `dMoverlap = sqrt(pi/3)` / `dKoverlap = 23 sqrt(pi)/(3 sqrt(3))` (math:126-127);
- real-Y20 overlap ratios `{1, 1/2, -1}` for `m = 0, 1, 2` (sympy:195-200 / math:175-184);
- grouped trace/line identities `xbar = x0`, `bx = 3*ax` with `lam20=1, lam21=1/2, lam22=-1` (sympy:202-208 / math:186-196).

Nothing in the stage card, the eq:stage015-keta formula, the appendix row, or any (nonexistent) notes file references B/Z modes, K1/H_even gate algebra, Y20 STF coefficients, or a grouped trace-anomaly identity. The sympy docstring (line 2) even names the source as `step_13_parent_throat_action_master_notes.md` - a note that does not exist in `notes/stages/`. The script appears written against an unindexed master-note that may have been split or absorbed into other stages; without that source available, the wall-only / Y20 / grouped blocks are unmotivated relative to the published stage card.

**Why this matters:**
A reader of the paper card has no way to know what the wall-only K1/H_even, Y20 ratios, or grouped trace identities are doing in this audit, or what claim they substantiate. Either the audit over-promises (the script verifies more than the stage exports) or the stage card under-promises (the gate/Y20/grouped scaffolding is genuinely a deliverable that the card forgot to mention). Direction of resolution is the user's call - Codex must not silently trim either side.

**Required change:**
`## Resolve before fix_loop` (see directive). The directive asks the user whether (a) the paper card should be expanded to mention the wall-only gates, Y20 ratios, and grouped trace as additional deliverables (with reference to the missing `step_13_parent_throat_action_master_notes.md` source), or (b) the script blocks M4-M9 (sympy lines 103-208; mathematica lines 104-196) should be trimmed because they were imported from a stale master-note draft that no longer belongs in stage 015.

**Verification:**
After user picks a direction: if (a), `paper/stages/stage_015.tex` acquires new paragraphs and/or a notes file is added under `notes/stages/moving_throat_pde_stage015_*.md`; if (b), the script is trimmed to A1-A8 (sympy) and M1-M7 (mathematica), and the stale outputs are regenerated.

### F2 - insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:55-72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:67-81`

**What's wrong:**
The "concrete quadratic IBP with boundary" block uses `A_concrete = exp(-w^2)` and `eta_concrete = exp(-w^2/2)`. Both factors are even in `w`, so the cross integrand `-A * eta * eta_w = -exp(-w^2) * exp(-w^2/2) * (-w * exp(-w^2/2)) = w * exp(-2 w^2)` is odd and vanishes by parity over `(-oo, oo)`. The bulk integrand `(1/2) A_w * eta^2 = (1/2)(-2w) exp(-w^2) * exp(-w^2) = -w * exp(-2 w^2)` is also odd and vanishes. The boundary discharge `[-A * eta^2/2]` from `-oo` to `+oo` is 0 - 0 = 0. So the assertion `quad_cross_concrete - (quad_boundary_concrete + quad_bulk_concrete) = 0 - (0 + 0) = 0` is `0 = 0` - it does not exercise the IBP identity on a profile that produces nonzero pieces on either side. The paper card promises "both zero and nonzero boundary-discharge probes"; the nonzero probe (sympy line 57-58, `atan(w)`) is wired only to the `boundary_value` operator sanity check, never to the IBP identity itself. So the concrete IBP test never sees a case where boundary discharge is actually nonzero and IBP cancellation is doing real work.

**Why this matters:**
The IBP check at the concrete level is decorative; only the generic symbolic identity A1 carries the IBP claim, and a wrong sign in the concrete instantiation would still pass because all three pieces are zero independently. The paper-promised "nonzero discharge probe" is satisfied only by an unrelated `atan` test that has nothing to do with the action density.

**Required change:**
Add a second concrete pair whose profiles break the global parity so the integrals are nonzero and the IBP cancellation is meaningful. Concrete suggestion: keep the existing Gaussian pair as a degenerate-baseline check, and add new asserts using `A_extra = exp(-w^2)`, `eta_extra = w * exp(-w^2/2)`. Then the cross integrand is `-A_extra * eta_extra * D[eta_extra, w]`. With `eta_extra = w * exp(-w^2/2)`, `D[eta_extra, w] = exp(-w^2/2) - w^2 * exp(-w^2/2) = (1 - w^2) exp(-w^2/2)`. So `eta_extra * D[eta_extra, w] = w (1 - w^2) exp(-w^2)`; `A_extra * eta_extra * D[eta_extra, w] = w (1 - w^2) exp(-2 w^2)`, odd in `w`, still vanishes. Better profile: `A_extra = 1/(1 + w^2)` (even), `eta_extra = exp(-w^2/2)/(1+w^2)` - still even-times-even, same problem. The robust fix is to break the symmetry of `A`: use `A_extra(w) = w * exp(-w^2)` (odd) and keep `eta_extra = exp(-w^2/2)` (even). Then `A_extra * eta * eta_w = w * exp(-w^2) * exp(-w^2/2) * (-w * exp(-w^2/2)) = -w^2 * exp(-2 w^2)` is even, integral is nonzero (`-sqrt(pi/2)/4`). The bulk integrand `(1/2) A_extra' * eta^2` where `A_extra' = exp(-w^2) - 2 w^2 exp(-w^2) = (1 - 2 w^2) exp(-w^2)`, giving `(1/2)(1 - 2 w^2) exp(-2 w^2)`, also even and nonzero (`sqrt(pi/2)/4`). Boundary discharge `[-A_extra * eta_extra^2 / 2]_{-oo}^{+oo} = 0 - 0 = 0`. So the cross now equals bulk (both nonzero with opposite signs from the IBP sign), and the assertion `cross - (boundary + bulk) = 0` is non-trivial. Apply the analogous change in mathematica.

**Verification:**
After the fix, the residual print lines for the new concrete IBP block should display nonzero `quad_cross_concrete_extra` and `quad_bulk_concrete_extra` summands (with the residual still 0 after cancellation). Mathematica's M2 output likewise. The existing Gaussian-pair asserts may remain as a parity baseline but must be augmented, not replaced silently.

### F3 - tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:118-127`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl:104-113`

**What's wrong:**
The "wall-only K1/H_even specialization" asserts are tautological by construction. The script defines:
```
D01_full = dK - B01 - Z01
D21_full = -(dM + B21 + Z21)
D41_full = -(B41 + Z41)
K1_full  = D21_full + D01_full / 9
H_even_full = D41_full - (2/3)*D21_full - D01_full / 27
```
then substitutes `{B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}` to produce `K1_wall = -dM + dK/9` and `H_even_wall = (2/3) dM - dK/27`, and immediately asserts these match `-dM + dK/9` and `(2/3) dM - dK/27`. The right-hand sides are not derived from anything physical; they are simply the same substitution applied by inspection. The assert cannot fail under any choice of coefficient because the LHS *is* the substitution. The same pattern is mirrored in mathematica M4. (This contrasts with A11-A12 which then substitute Gaussian-overlap closed forms - those are non-tautological *within* the extra scope, but the substitution identities A9-A10 / M4 are not.)

**Why this matters:**
These two asserts spend script lines without testing anything. They look load-bearing in the verdict ("PASS") but cannot detect any error. They should either be removed, or replaced with a derivation of the K1/H_even reduction from the actual `\StatusExactClosure` reduction paper claim (which would lift the whole block out of `paper_missing_script_claim` territory if such a paper claim were available).

**Required change:**
If F1 resolves toward keeping the wall-only block (paper-side expansion), replace A9-A10 / M4 with a derivation that traces back to a non-trivial paper-side coefficient claim (the values `1/9`, `-1`, `2/3`, `-1/27` for the linear combinations defining K1 and H_even must come from somewhere paper-side - quote that source and derive the wall-only forms by independent substitution). If F1 resolves toward removing the wall-only block, delete A9-A10 (sympy 126-127) and M4 (mathematica 112-113) along with the rest of the wall-only scaffolding. Either way, do not leave a tautological substitution-against-itself check.

This finding is therefore **blocked on F1** - the right action depends on F1's user resolution.

**Verification:**
The assertion must, after fix, be capable of failing under at least one wrong choice of coefficient. If the derivation paths to K1/H_even are independent (one symbolic substitution, one derivation from the paper-stated invariant), the assert is non-tautological.

### F4 - mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl`

**What's wrong:**
The Mathematica script's algebraic spine is a structurally identical re-encoding of the SymPy script rather than an independent re-derivation. Paired examples:

- SymPy 81-86 vs Mathematica 83-87 - same `Tw = Tw0 + eps*TwR0*eta + eps^2*TwRR0*eta^2/2`, same `U`, same `Rt = eps*eta_t`, same `Rw = R0p + eps*eta_w`, same `L` decomposition.
- SymPy 87 `L2_raw = sp.diff(sp.expand(L), eps, 2).subs(eps, 0) / 2` vs Mathematica 88 `L2raw = Coefficient[Series[lagrangian, {eps, 0, 2}] // Normal, eps, 2]` - different SDK function calls, same operation (extract eps^2 coefficient).
- SymPy 96 `L2_after_ibp_derived = sp.expand(L2_raw - cross_term + cross_after_ibp)` vs Mathematica 95 `L2afterIBP = Expand[L2raw - (-TwR0*R0p*eta*etaw) + dTwRR0p*eta^2/2]` - character-for-character algebraic mirror.
- SymPy 118-125 vs Mathematica 104-111 - same `D01_full`, `D21_full`, `D41_full`, same `K1_full = D21_full + D01_full/9`, same `wall_only_specialization`/`wallSpec`. The variable choreography is line-for-line preserved.
- SymPy `real_y20_square_ratio` (25-31) vs Mathematica `realY20Ratio` (174) - identical normalization-by-`gauntBase`.
- SymPy 202-208 vs Mathematica 189-196 - identical `xbar = (x20 + 2 x21 + 2 x22)/5`, `ax = (2 x20 - x21 - x22)/10`, `bx = (x21 - x22)/2`.

The two scripts do not derive the K_eta formula or the wall-only K1/H_even combinations from independent physical premises; they share the same intermediate algebraic representation. The Mathematica does add `Sqrt[Pi/3]`, `23 Sqrt[Pi]/(3 Sqrt[3])`, and `2*eps/3` as extra closed-form checkpoints (M5, M6), which is genuine independent content - but the algebraic backbone is a port.

**Why this matters:**
A second engine is meant to catch silent algebra/typo errors in the first. When both engines walk the same intermediate steps with the same variable names and decompositions, an error in the derivation strategy propagates into both engines identically and is invisible to the cross-check. The closed-form checkpoints added in M5/M6 are good but do not redeem the K_eta or wall-only blocks.

**Required change:**
Restructure the Mathematica K_eta block (M3) to derive K_eta independently, *without* pre-expanding `lagrangian` to eps^2 with `Series[...] // Normal` and pattern-matching against a pre-written `canonicalL2`. Concrete prescription: compute the Euler-Lagrange operator on `lagrangian` for the field `R` directly, linearize around `R = R0(w)` by writing `R = R0[w] + eps*eta[t,w,Omega]` (treat eta as a separate field), expand the EL equation to O(eps), and read off the linearized PDE's mass coefficient. Specifically, define `RFull[t_, w_, om_] := R0[w] + eps*eta[t,w,om]` and `LSig[R_, Rt_, Rw_, gO_, w_] := mu0*Rt^2/2 - (Tw0 + (R - R0[w])*TwR0 + (R - R0[w])^2*TwRR0/2)*Rw^2/2 - TO0*gO/2 - (U0 + (R - R0[w])*UR0 + (R - R0[w])^2*URR0/2)`. Then the EL coefficient of `eta` linearized in eps gives the K_eta combination via `D[LSig, R] - D[D[LSig, Rt], t] - D[D[LSig, Rw], w]` evaluated at R = R0[w] with appropriate variable replacements. Compute that combination, extract the linear-in-eta coefficient (the mass term), and compare to `URR0 - dTwRR0p + TwRR0*R0p^2/2`. This derivation path does not go through `Series` and does not share variable names with the sympy script's intermediate `L2raw`/`canonicalL2` decomposition.

For the wall-only block (M4): the wall-only K1/H_even combinations are tautological (F3) - resolve F3 first; the transliteration call applies only conditionally if the wall-only block survives F1 resolution.

**Verification:**
After restructuring, the Mathematica script should reach the same K_eta formula but via Euler-Lagrange linearization, not via series-expansion-and-pattern-match. The diff between the SymPy and Mathematica scripts at the K_eta-block level should not be a literal `s/sp.diff/D/`, `s/sp.expand/Expand/` rewrite. Search-and-replace transliteration should be impossible to construct.

## Independent-derivation check (Mathematica)

See F4. The Mathematica script's spine is a transliteration of the SymPy script's algebra for the K_eta block and the wall-only block. Two-paired examples justifying the call:

SymPy 95-96:
```
cross_term = -TwR0 * R0p * eta * eta_w
cross_after_ibp = d_TwR_R0p * eta**2 / 2
L2_after_ibp_derived = sp.expand(L2_raw - cross_term + cross_after_ibp)
```
Mathematica 95:
```
L2afterIBP = Expand[L2raw - (-TwR0*R0p*eta*etaw) + dTwRR0p*eta^2/2];
```

SymPy 118-127:
```
D01_full = dK - B01 - Z01
D21_full = -(dM + B21 + Z21)
D41_full = -(B41 + Z41)
K1_full = sp.expand(D21_full + D01_full / 9)
H_even_full = sp.expand(D41_full - sp.Rational(2, 3) * D21_full - D01_full / 27)
wall_only_specialization = {B01: 0, B21: 0, B41: 0, Z01: 0, Z21: 0, Z41: 0}
K1_wall = sp.expand(K1_full.subs(wall_only_specialization))
H_even_wall = sp.expand(H_even_full.subs(wall_only_specialization))
```
Mathematica 104-111:
```
D01full = dK - b01 - z01;
D21full = -(dM + b21 + z21);
D41full = -(b41 + z41);
K1full = D21full + D01full/9;
Hevenfull = D41full - (2/3)*D21full - D01full/27;
wallSpec = {b01 -> 0, b21 -> 0, b41 -> 0, z01 -> 0, z21 -> 0, z41 -> 0};
K1wall = Expand[K1full /. wallSpec];
Hevenwall = Expand[Hevenfull /. wallSpec];
```
This is a port, not a re-derivation. The closed-form Gaussian checkpoints in M5 and the perturbation value `2*eps/3` in M6 are independent content; the K_eta and wall-only spines are not.

## Engine cross-check

Both engines produce identical PASS verdicts across all asserts. SymPy output (7 lines) is a summary; Mathematica output (56 lines) reports per-check residuals, all "= 0" or matching the expected nonzero residual. No engine disagreement. The Mathematica script also independently confirms two closed forms the SymPy script does not state explicitly: `dMoverlap = Sqrt[Pi/3]` and `dKoverlap = 23 Sqrt[Pi]/(3 Sqrt[3])`. Outputs are fresh (sympy txt mtime 1779397234 > script 1779390705; mathematica txt 1779397279 > script 1779390746).

## Verdict justification

`findings`. The K_eta exact-quadratic-recovery deliverable that the paper card calls out is genuinely verified (sympy A6-A8 / mathematica M5-M7) with an honest sign-mutation guard. However: (1) about half of each script is dedicated to wall-only gate algebra, Y20 overlap ratios, and grouped-trace identities that the paper card never mentions - this is paper_misalignment requiring user resolution, not a Codex fix. (2) The concrete Gaussian boundary-discharge probe is degenerate (cross, bulk, and boundary are all individually zero by parity), so it does not exercise the IBP cancellation. (3) The wall-only K1/H_even specialization asserts are tautological substitutions of the same coefficients against themselves (and depend on F1's resolution for the right action). (4) The Mathematica script's K_eta and wall-only spines are transliterations of the SymPy spine, weakening the second-engine cross-check (independent closed-form checkpoints in M5/M6 are a partial mitigation only). I am not flagging `CRITICAL_DOWNSTREAM` because the K_eta formula itself does match the paper, and the questionable extra blocks (wall-only, Y20, grouped) would only propagate downstream if downstream stages depend on them - which is exactly what the user-side resolution of F1 will determine.

## Self-test notes

For F2 I confirmed by symmetry argument that `A_concrete = exp(-w^2)` and `eta_concrete = exp(-w^2/2)` produce odd integrands on both the cross and bulk integrals over the symmetric domain `(-oo, oo)`, so both individually vanish - the IBP identity at the concrete level is `0 = 0 + 0`. I then mentally pre-checked the F2 fix profile `A_extra = w * exp(-w^2)`, `eta_extra = exp(-w^2/2)`: cross integrand `w * exp(-w^2) * exp(-w^2/2) * (-w) * exp(-w^2/2)` is even (w^2 factor), integral nonzero; A_extra' = (1 - 2w^2) exp(-w^2), bulk integrand `(1 - 2w^2) exp(-2w^2)/2` is even, integral nonzero; boundary `[-A_extra * eta_extra^2 / 2]` at +-oo is 0; so the IBP cancellation cross - bulk = 0 is a real (nonzero=nonzero) identity. For F3 I traced both `K1_wall` and `(-dM + dK/9)` back to the same substitution and confirmed the assert is provably 0 by construction. For F4 I paired four spine blocks line-by-line across the two scripts to confirm the transliteration call. For F1, I confirmed by file listing that no `notes/stages/moving_throat_pde_stage015_*.md` file exists, and that the paper card and appendix row contain no mention of wall-only gates, Y20 ratios, or grouped trace identities; therefore the routing of F1 to user resolution (not Codex auto-fix) is correct. I also verified the `assert_zero("wall-only dK closed form", wall_even_solve[dK])` substantively: solving `-dM + dK/9 = 0, (2/3) dM - dK/27 = 0` yields the unique trivial kernel `dK = 0, dM = 0`, so the assert is non-tautological even though the assertion's label "closed form" is misleading (there is no nontrivial closed form, only the trivial kernel).
