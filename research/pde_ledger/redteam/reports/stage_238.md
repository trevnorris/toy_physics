---
unit_id: 238
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 4
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 238 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_238.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 88, 1067-1099)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.txt`
- mathematica output: `(missing)`

## What the paper claims

Stage 238 rewrites the rigid-mouth same-charge packet in physical transfer-shape variables `(R_tr, T^2, eps_eta)`. The card `\stagefield{Output}`: "Physical transfer-shape compiler: the finite rigid-mouth packet is completely described by tracking, $\mathcal T^2$, and $\epsilon_\eta$." The card's derivation ledger and notes enumerate these distinct deliverables: (1) the exact coherent identity `R_target·T^2 = Lambda0(1-eps_eta)`; (2) the finite packet factorization `q_nt + (B*/C*)q_tr = ln(T^2/T_ref^2)`; (3) the rigid-mouth finite reduction `q_nt = ln(T^2/T_ref^2)`, `q_eta = ln(eps_eta/eps_eta_ref)`; (4) the exact first-order drift compiler for `δln R_tr`, `δln T^2`, `δln R_target` (three explicit closed forms in notes §3); (5) the first-order packet relation `q_nt + (B*/C*)q_tr = δln T^2`, `q_eta = d ln eps_eta`, and its rigid-mouth slice `q_nt = δln T^2`; (6) support-blindness of all three q's w.r.t. `zeta` and `M_mix`, with the physical content that the support sector enters only through `M_tr = M_mix[1 + zeta(1-eps)/(1-zeta·eps)]` (notes §4) yet the direct observables do not; (7) the three-gate theorem (tracking / transfer-shape / dressing).

## What the script claims to verify

The docstring (lines 4-14) claims the script verifies all six clusters above. It defines the coherent observables `Rtr`, `T2`, `Rtarget` explicitly (lines 57-59), then asserts the transfer-shape identity (62), the finite factorization via an exp-comparison (78-82), the rigid-mouth finite reduction (85-90), the three first-order compilers via `h`-perturbation and `diff(...,h)|_{h=0}` (96-137), the first-order packet relation (143), two rigid-mouth first-order lines (147, 151), twelve support-blindness `diff` checks (157-171), the tracking-gate compiler (178-184), and the two finite gates (187-188). All checks are exact symbolic (`assert_zero` via `simplify(factor(...))`).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `R_target·T^2 = Lambda0(1-eps_eta)` | line 62 | match (observables defined explicitly; non-tautological) |
| (2) finite factorization `q_nt+(B/C)q_tr = ln(T^2/T_ref^2)` | lines 78-82 | match (B-terms cancel correctly; exp-comparison reduces to identity (1)) |
| (3) rigid finite `q_nt = ln(T^2/T_ref^2)`, `q_eta = ln(...)` | lines 85-90 | match |
| (4) first-order `δln R_tr`, `δln T^2`, `δln R_target` | lines 102-137 | match (perturbation `diff` recovers the three hardcoded `_expected` forms; genuine) |
| (5a) first-order packet `q_nt+(B/C)q_tr = δln T^2` | line 143 | match (B-terms cancel; reduces to first-order of identity (1); genuine) |
| (5b) rigid first-order `q_nt = δln T^2` | lines 146,150-151 | mismatch/missing — line 147 is the trivial `q_eta` identity; line 151 duplicates line 118; the rigid `q_nt` reduction is never asserted (F2) |
| (6) support-blindness of q's w.r.t. zeta, M_mix | lines 157-171 | mismatch — all 12 are tautological: `zeta`/`M_mix` never appear in the expressions (F1) |
| (7) three-gate theorem (tracking) | lines 178-184 | mismatch — tautological re-arrangement of the hardcoded `dlnRtr_expected` (F3) |
| (7) three-gate (transfer-shape, dressing finite gates) | lines 187-188 | match (line 187 = identity (1); line 188 is the `eps_eta` ratio) |

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 62 | `simplify(Rtarget*T2 - Lambda0(1-eps_eta)) == 0` | claim 1 | yes |
| A2 | sympy | 78-82 | `exp(q_nt+(B/C)q_tr|sub) - T2/T2_ref == 0` | claim 2 | yes |
| A3 | sympy | 85-89 | `exp(q_nt_rigid) - T2/T2_ref == 0` | claim 3 | yes |
| A4 | sympy | 90 | `exp(q_eta) - eps_eta/eps_eta_ref == 0` | claim 3 (eta) | yes |
| A5 | sympy | 106 | `dlnRtr_calc - dlnRtr_expected == 0` | claim 4 (R_tr) | yes |
| A6 | sympy | 118 | `dlnT2_calc - dlnT2_expected == 0` | claim 4 (T^2) | yes |
| A7 | sympy | 137 | `dlnRtarget_calc - dlnRtarget_expected == 0` | claim 4 (R_target) | yes |
| A8 | sympy | 143 | `q_nt_first_direct + (B/C)q_tr_first - dlnT2_calc == 0` | claim 5a | yes (B cancels; reduces to 1st-order of id 1) |
| A9 | sympy | 147 | `q_eta_first - dlnepseta == 0` | claim 5 (eta) | no — `q_eta_first := dlnepseta` (tautology, F2) |
| A10 | sympy | 150-151 | `q_nt_factorized_rigid - dlnT2_expected == 0` | claim 5b (rigid q_nt) | no — alias of A6; does not touch `q_nt` (F2) |
| A11 | sympy | 157-162 | `diff(Rtr/T2/eps_eta, zeta or M_mix) == 0` | claim 6 | no — vars absent (tautology, F1) |
| A12 | sympy | 166-171 | `diff(q_tr/q_nt_factored/q_eta, zeta or M_mix) == 0` | claim 6 | no — vars absent (tautology, F1) |
| A13 | sympy | 178-184 | `dlnRtr_expected + prefactor*tracking_condition == 0` | claim 7 (tracking) | no — re-arranges hardcoded `dlnRtr_expected` (tautology, F3) |
| A14 | sympy | 187 | `exp(q_nt_rigid) - T2/T2_ref == 0` | claim 7 (transfer-shape gate) | yes (= A3) |
| A15 | sympy | 188 | `exp(q_eta) - eps_eta/eps_eta_ref == 0` | claim 7 (dressing gate) | yes (= A4) |

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py:157-171`

**What's wrong:**
All twelve support-blindness assertions differentiate expressions with respect to `zeta` and `M_mix`, but neither symbol ever appears in any of the differentiated expressions. `Rtr` (line 57) depends only on `chi0, deltaU`; `T2` (line 58) on `ZW, chi0, OmegaW2, eps`; `eps_eta` is a bare symbol; `q_tr`, `q_nt_factored`, `q_eta` are built only from those plus `B, C, *_ref`. `zeta` and `M_mix` are declared at line 46 but are dead symbols. Therefore every `sp.diff(EXPR, zeta)` and `sp.diff(EXPR, M_mix)` is identically zero by construction and the assertions cannot fail regardless of the physics.

The paper's actual support-blindness claim (notes §4, lines 239-267) has real content: the support sector enters through `M_tr = M_mix[1 + zeta(1-eps)/(1-zeta·eps)]`, and the nontrivial statement is that the direct observables, when expressed through the channel that *does* see `M_tr`, nonetheless have zero derivative. The script never introduces that `M_tr` dependence, so it proves nothing about whether support enhancement cancels.

**Why this matters:**
This is the central self-test trap-1 failure mode. Six paper deliverables (support-blindness of three observables and three q's) are reported "[ok]" in the output while being vacuously true. A reader trusts that support-blindness was verified when it was not exercised at all.

**Required change:**
Replace the dead-symbol checks with ones that route `zeta` and `M_mix` through the support channel the notes define, so the derivative is a genuine cancellation rather than a derivative-of-a-constant. Codex must (a) define `M_tr = M_mix*(1 + zeta*(1 - eps)/(1 - zeta*eps))` as the notes §4 prescribe, (b) construct the observables in the form where that factor would appear if the reduction had not already cancelled it (i.e., reintroduce the baseline `M_tr` factor that the reduced observables claim to have dropped), and (c) assert that `diff(observable, zeta)` and `diff(observable, M_mix)` vanish only after the cancellation, so a sign/factor error in the construction would make the assert fail. If the only honest statement is that the reduced observables are literally independent of the support symbols, the check must at least feed `zeta`/`M_mix` into the observables through a non-trivial substitution path so the vanishing is a derived fact, not an artifact of absent symbols. See claim manifest M1-M2 in the directive.

**Verification:**
After the fix, each of the (rewritten) lines 157-171 must contain at least one expression in which `zeta` (resp. `M_mix`) literally appears before differentiation; the verifier confirms the substituted expression is non-constant in the support symbol and that the residual still simplifies to zero.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py:142-151`

**What's wrong:**
The rigid-mouth first-order reduction (paper claim 5b / notes §3 boxed `q_nt = δln T^2`) is not actually asserted. Line 142 sets `q_eta_first = dlnepseta`; line 147 then asserts `q_eta_first - dlnepseta == 0`, which is `dlnepseta - dlnepseta == 0`, a pure tautology whose label ("rigid-mouth first-order dressing compiler") oversells it. Line 146 builds `q_nt_first_rigid = q_nt_first_direct.subs(dlnRtr_calc, 0)` but never uses it in any assertion. Line 150 sets `q_nt_factorized_rigid = dlnT2_calc` and line 151 asserts `dlnT2_calc - dlnT2_expected == 0` — identical to the already-made line 118 check, so it is a duplicate of A6 and has nothing structurally to do with `q_nt`. Net: the load-bearing rigid first-order statement `q_nt_first_rigid == dlnT2_calc` is never tested.

**Why this matters:**
A reader sees "[ok] rigid-mouth first-order transfer-shape compiler" and believes the rigid `q_nt = δln T^2` reduction is verified; it is not. The genuine non-trivial check (I confirmed by hand that `q_nt_first_direct|_{dlnRtr->0} = dlnT2_calc`, via `dlnRtarget_calc = -dlnT2_calc - eps_eta/(1-eps_eta)·dlnepseta`) is omitted.

**Required change:**
At line 147, replace the tautological `assert_zero(q_eta_first - dlnepseta, ...)` with the genuine dressing check tied to the perturbation result, e.g. assert `q_eta_first - dlnepseta` where `q_eta_first` is recomputed from a perturbed `eps_eta` (not defined as `dlnepseta`) — or, more simply, drop line 147's duplicate. At line 151, replace the duplicate-of-A6 assertion with the actual rigid-mouth claim: assert `q_nt_first_rigid - dlnT2_calc == 0` using the `q_nt_first_rigid` already built at line 146. See claim manifest M3.

**Verification:**
The new line must reference `q_nt_first_rigid` (line 146) and `dlnT2_calc` (line 116) and simplify to zero; the verifier confirms the assertion now exercises `q_nt`, not just `dlnT2`.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py:177-184`

**What's wrong:**
The "tracking gate compiler" assertion is `dlnRtr_expected + chi0·deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU))·tracking_condition == 0`, with `tracking_condition = (1+deltaU)dlnchi + (1+chi0)dlndelta` (line 177). But `dlnRtr_expected` is defined at lines 103-105 as exactly `-chi0·deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU))·[(1+deltaU)dlnchi + (1+chi0)dlndelta]`, i.e. `-prefactor·tracking_condition`. So the assertion is `(-prefactor·tracking_condition) + (prefactor·tracking_condition) == 0` — algebraically guaranteed by construction. It re-arranges the hardcoded `_expected` form against itself and never touches the independently `diff`-derived `dlnRtr_calc`. The tracking gate's physical content (`R_tr = R_tr_ref ⟺ tracking_condition = 0` at first order, notes §5.1) is not exercised.

**Why this matters:**
The tracking gate is one of the three pillars of the headline three-gate theorem (claim 7). The check that purports to verify it is vacuous.

**Required change:**
Tie the tracking gate to the independently derived `dlnRtr_calc` (line 102), not to `dlnRtr_expected`. For instance, assert that `dlnRtr_calc` factors as `-prefactor·tracking_condition` by checking `dlnRtr_calc + prefactor·tracking_condition == 0` (using `dlnRtr_calc`, the `diff`-derived quantity, so a perturbation error would surface), and additionally assert the gate statement itself: `dlnRtr_calc.subs(tracking_condition -> 0)`-style — i.e. substitute the gate condition into the perturbation symbols and confirm `dlnRtr_calc` vanishes. See claim manifest M4.

**Verification:**
The rewritten assertion must reference `dlnRtr_calc` (line 102), not `dlnRtr_expected`; the verifier confirms the check now depends on the perturbation derivation.

### F4 — missing_verification_script (missing_mathematica)

**Severity:** high
**Files:**
- `(missing)` — no `mathematica/moving_throat_pde_stage238_*.wl`

**What's wrong:**
This unit is `is_status_only_candidate: False` and `is_checkpoint: False` but has only a SymPy script. Every Stage 238 claim is exact symbolic algebra over a handful of rational/log expressions and first-order log-derivatives — Mathematica can verify all of it independently (native `D`/`Series` for the first-order compilers, `FullSimplify`/`Solve` for the transfer-shape identity and the exp-form factorizations, genuine support-channel substitution for support-blindness). The dual-engine rule's test is "is it possible," which is clearly yes. No `.wl` exists.

**Why this matters:**
Single-engine coverage on a non-checkpoint, non-status stage violates the dual-engine policy and leaves all claims unconfirmed by an independent route — especially important here because the SymPy script's support-blindness and tracking-gate checks (F1, F3) turned out tautological, so there is currently zero independent confirmation of those deliverables.

**Required change:**
Codex writes a NEW independent-route Mathematica script `mathematica/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.wl` that derives the claims via a different decomposition than the `.py` (NOT a transliteration). See the claim manifest M1-M7 in the directive.

**Verification:**
`redteam exec-mathematica 238` runs the new `.wl`, all in-file checks pass, exit 0; the support-blindness and tracking-gate checks must route the support symbols / gate condition through a non-trivial substitution (mirroring the F1/F3 fixes) so they are not vacuous.

## Independent-derivation check (Mathematica)

No `.wl` exists; `mathematica_transliteration` does not apply. The required new script must be an independent re-derivation (see F4 / directive manifest).

## Engine cross-check

Only one engine present; `engine_disagreement` does not apply. The SymPy output (exit 0, all 25 lines "[ok]") is internally consistent but, per F1/F2/F3, several "[ok]" lines are vacuous.

## Verdict justification

Verdict is `findings`. The substantive algebraic core holds up under attack: the transfer-shape identity (A1), the finite factorization (A2, with the B/C tracking feed-through correctly cancelling), the rigid finite reduction (A3-A4), the three first-order drift compilers (A5-A7, genuinely recovered by `h`-perturbation), and the first-order packet relation (A8, which I verified reduces to the first-order of the transfer-shape identity once the B-terms cancel) are all real and paper-aligned. But four deliverables are not honestly exercised: the twelve support-blindness checks are tautological because `zeta`/`M_mix` never enter the expressions (F1, high); the rigid-mouth first-order `q_nt` reduction is replaced by a trivial `q_eta` identity and a duplicate `dlnT2` check (F2, medium); the tracking-gate check re-arranges the hardcoded `dlnRtr_expected` against itself (F3, medium); and a required independent Mathematica engine is absent (F4, high). No `paper_misalignment`: the paper card, notes, and appendix are mutually consistent and agree with what the genuine checks verify; the issue is script-side rigor, not a paper↔script value disagreement. Not `stop_cold`: every finding is mechanically fixable within scope, and no downstream constant/sign would change (the math is already correct; the fixes make the checks actually test it).

## Self-test notes

Trap 1 (variable independence): this is exactly F1 — `sp.diff(Rtr/T2/eps_eta/q_*, zeta or M_mix)` is identically zero because those symbols are absent; flagged. Trap 3 (trivial-case pre-check): I hand-verified A8 reduces to `dlnRtarget_calc + dlnT2_calc + eps_eta/(1-eps_eta)dlnepseta == 0` (B-terms cancel) and that A2's exp-comparison reduces to the transfer-shape identity, both genuinely nonzero before cancellation — so A8/A2 are sound, not tautological. Trap 5 (paper round-trip): the prescribed F1/F2/F3 fixes reuse the exact constants and forms in notes §3-§5 (`M_tr = M_mix[1+zeta(1-eps)/(1-zeta·eps)]`, `q_nt_first_rigid = q_nt_first_direct|_{dlnRtr->0}`, the tracking prefactor), introducing no new paper_misalignment. Output freshness (mtime): script 11:58 < output 12:51, so no `stale_output`.
