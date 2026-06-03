---
unit_id: 238
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-03T08:25:16-06:00
verification_status: pending
needs_user_resolution: false
findings_applied: 4
findings_blocked: 0
---

# Codex directive — unit 238

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

SymPy file:
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`

---

## ITER 2 — orchestrator reframe of F1 + F4 (2026-06-03)

F2 and F3 are ALREADY APPLIED (their `## Applied` blocks below stand — do NOT redo them). F1 and F4 were correctly Blocked on iter 1: the notes define `M_tr` separately (§4 lines 241–242) but give NO pre-reduction observable form where `M_tr` cancels, and inventing one was forbidden. F1 and F4 are now REFRAMED to a faithful, non-vacuous form (negative control + leak detector + exclusion) that does NOT route `M_tr` into the observables. Apply the reframed F1 (SymPy) and write the reframed F4 `.wl`; the old `## Blocked` blocks are SUPERSEDED — replace each with an `## Applied: F<n>` block when done.

## F1 — tautological_check (support-blindness)

**Target:** `scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py:157-171`

**Issue:** All twelve support-blindness `sp.diff(..., zeta)` / `sp.diff(..., M_mix)` assertions are identically zero by construction: `zeta` and `M_mix` (declared line 46) never appear in `Rtr` (line 57), `T2` (line 58), the bare symbol `eps_eta`, or in `q_tr` / `q_nt_factored` / `q_eta`. The checks cannot fail and do not exercise the paper's support-blindness claim (notes §4), whose content is that the support sector enters through `M_tr = M_mix[1 + zeta(1-eps)/(1-zeta·eps)]` yet the direct observables drop out of it.

**Required change (orchestrator-reframed, iter 2):**
Do NOT route `M_tr` into the observables (the notes give only the reduced observables and define `M_tr` separately at notes §4 lines 241–242; there is no in-stage pre-reduction form where `M_tr` cancels, and inventing one is forbidden). Instead make the existing `∂_ζ(observable)=0` / `∂_{M_mix}(observable)=0` checks NON-VACUOUS by adding two live guards that prove the test discriminates:

1. **Support channel definition.** Define `M_tr = M_mix*(1 + zeta*(1 - eps)/(1 - zeta*eps))` exactly as notes §4.
2. **Negative control (support channel is live).** Assert `M_tr` genuinely depends on the support symbols: `sp.diff(M_tr, zeta)` is NOT identically zero (it equals `M_mix*(1 - eps)/(1 - zeta*eps)**2`) and `sp.diff(M_tr, M_mix)` is NOT zero. Use an `assert sp.simplify(...) != 0` form, or assert each derivative equals its known nonzero closed form. This proves `zeta`, `M_mix` are real, in-scope symbols — so an observable's `∂_ζ = 0` pass cannot be dismissed as "the symbol was never defined."
3. **Leak detector (the test catches contamination).** Construct a hypothetical support-contaminated observable, e.g. `Rtr_leak = Rtr * M_tr / M_mix`, and assert `sp.diff(Rtr_leak, zeta) != 0`. This demonstrates that IF the support sector had leaked into an observable, the `∂_ζ` check WOULD fail — i.e. the check is discriminating, not vacuous.
4. **The claim (kept, now load-bearing).** Retain the six observable/packet `∂_ζ = 0` and `∂_{M_mix} = 0` assertions (lines 156–170); they now carry weight because (2)–(3) prove the differentiation is live and the symbols are present. As belt-and-suspenders, also assert the structural exclusion (`zeta not in obs.free_symbols and M_mix not in obs.free_symbols`) for each of `R_tr, T^2, eps_eta, q_tr, q_nt_factored, q_eta`.

This is faithful to notes §4 (support enters ONLY through `M_tr`; the direct observables are `M_tr`-free) and invents no placement. A future edit that let support leak into an observable trips step 4; a typo that made `zeta` absent from `M_tr` trips step 2.

**Self-test reminder (must hold before you assert):** confirm `sp.simplify(sp.diff(M_tr, zeta) - M_mix*(1 - eps)/(1 - zeta*eps)**2) == 0` (and it is nonzero), `sp.diff(Rtr_leak, zeta) != 0`, and each true-observable `sp.diff(obs, zeta) == 0`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
- summary: Added the live `M_tr` support-channel negative controls, a support-leak detector, and structural support-exclusion guards before the reduced observable derivative checks.
- deviation: none

## F2 — insufficient_verification (rigid first-order q_nt)

**Target:** `scripts/...stage238..._sympy_audit.py:142-151`

**Issue:** The rigid-mouth first-order reduction `q_nt = δln T^2` (notes §3 boxed) is never asserted. Line 142 defines `q_eta_first = dlnepseta`; line 147 then asserts `q_eta_first - dlnepseta == 0` (a pure tautology). Line 146 builds `q_nt_first_rigid = q_nt_first_direct.subs(dlnRtr_calc, 0)` but never uses it. Line 150-151 asserts `dlnT2_calc - dlnT2_expected == 0`, a duplicate of the line-118 check, structurally unrelated to `q_nt`.

**Required change:**
1. At line 151, replace `assert_zero(q_nt_factorized_rigid - dlnT2_expected, "rigid-mouth first-order transfer-shape compiler")` with the genuine rigid `q_nt` claim using the already-built `q_nt_first_rigid` (line 146):
   assert that `q_nt_first_rigid - dlnT2_calc` simplifies to zero. (I confirmed by hand this holds: `q_nt_first_rigid = -dlnRtarget_calc - eps_eta/(1-eps_eta)·dlnepseta` and `dlnRtarget_calc = -dlnT2_calc - eps_eta/(1-eps_eta)·dlnepseta`, so the difference is 0.) Remove the now-unused `q_nt_factorized_rigid = dlnT2_calc` assignment at line 150.
2. At line 147, the `q_eta_first - dlnepseta` check is tautological. Either delete it (the finite dressing ratio is already covered by line 90/188), or replace it with a non-trivial dressing check: recompute the first-order dressing coordinate from a perturbed `eps_eta` (perturb `eps_eta -> eps_eta*exp(h*dlnepseta)`, take `diff(log(eps_eta_pert), h)|_{h=0}`) and assert that equals `dlnepseta`. Prefer the recompute (it exercises the perturbation machinery) over deletion.

**Self-test reminder:** confirm `simplify(q_nt_first_rigid - dlnT2_calc) == 0` with `q_nt_first_rigid` as defined at line 146 (subs `dlnRtr_calc -> 0`). Do not substitute into `dlnRtr_expected`; substitute into the `_calc` chain so the assert depends on the perturbation derivation.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
- summary: Recomputed `q_eta^(1)` from an exponential `eps_eta` perturbation and replaced the rigid first-order `q_nt` check with `q_nt_first_rigid - dlnT2_calc`.
- deviation: none

## F3 — tautological_check (tracking gate)

**Target:** `scripts/...stage238..._sympy_audit.py:177-184`

**Issue:** The "tracking gate compiler" assertion adds `dlnRtr_expected` (lines 103-105, `= -prefactor·tracking_condition`) to `+prefactor·tracking_condition`, so it is `(-X)+(X)==0` by construction and never references the independently `diff`-derived `dlnRtr_calc` (line 102). The gate's physical content (notes §5.1: `R_tr = R_tr_ref ⟺ (1+deltaU)dlnchi + (1+chi0)dlndelta = 0` at first order) is not exercised.

**Required change:**
Re-anchor the assertion to the `diff`-derived `dlnRtr_calc`, not `dlnRtr_expected`:
1. Assert the factorization on the derived quantity: `assert_zero(dlnRtr_calc + chi0*deltaU/((1+chi0)*(1+deltaU)*(1+chi0+deltaU))*tracking_condition, "tracking gate factorization")`. (Equivalent residual, but now a perturbation error in `dlnRtr_calc` would surface.)
2. Additionally assert the gate statement: substitute the gate condition into the *perturbation symbols* and confirm the drift vanishes. Concretely, solve `tracking_condition == 0` for one drift (e.g. `dlndelta = -(1+deltaU)/(1+chi0)*dlnchi`), substitute into `dlnRtr_calc`, and assert the result simplifies to zero — i.e. `q_tr_first` (= `-C*dlnRtr_calc`) vanishes exactly on the gate locus.

**Self-test reminder:** verify `dlnRtr_calc` (the `diff`-derived line-102 value) — not `dlnRtr_expected` — is what enters both new asserts. Substitute `dlndelta -> -(1+deltaU)/(1+chi0)*dlnchi` into `dlnRtr_calc` and confirm it reduces to 0 symbolically.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_sympy_audit.py`
- summary: Re-anchored the tracking gate factorization to `dlnRtr_calc` and added the gate-locus substitution check on the derived drift.
- deviation: none

## F4 — missing_verification_script (missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.wl`

(`.wl` lives in `mathematica/`. Filename suffix `_mathematica_audit.wl` is mandatory.)

**Issue:** Non-checkpoint, non-status-only stage with only a SymPy engine. Mathematica can independently verify every claim (the test is "is it possible," and it clearly is). No `.wl` exists.

**Required change:**
Write a NEW Mathematica script that derives the Stage 238 claims via a DIFFERENT decomposition than the `.py` — not a line-by-line transliteration. Independent-route guidance (Codex designs the actual code):
- Use native Mathematica primitives. For the first-order compilers, prefer `Series[Log[obs /. {chi0 -> chi0 Exp[h dlnchi], ...}], {h, 0, 1}]` and read the `h^1` coefficient, or `D[...,h] /. h->0` — but contrast against the `.py`'s `sp.diff(sp.log(simplify(...)), h)` path by, e.g., expanding the observables' logs additively first.
- For the transfer-shape identity and the exp-form factorizations, use `FullSimplify` / `Solve` directly on `Rtarget*T2 - Lambda0*(1-eps_eta)` and on the log-difference forms, rather than the `.py`'s `exp(...) - ratio` comparison.
- The support-blindness and tracking-gate checks MUST be written in the non-vacuous form prescribed in F1 and F3 (route `zeta`/`M_mix` through `M_tr = M_mix(1 + zeta(1-eps)/(1-zeta eps))`; anchor the tracking gate to the derived drift, not a hardcoded `_expected`). Do not reproduce the SymPy script's tautological versions.
- Provide an `expectZero[expr_, label_]` helper that strips `ConditionalExpression[0, ...]` and `Exit[1]`s on nonzero; print `[ok] label` on success.

**Claim manifest** (the new `.wl` must independently verify all of these; the SymPy fixes F1-F3 apply equally to the `.wl`):

- **M1 — Coherent transfer-shape identity:**
  `R_target · T^2 = Lambda0 (1 - eps_eta)`, with
  `R_tr = (1 + chi0/(1+deltaU))/(1+chi0)`,
  `T^2 = ZW (1+chi0)^2 / (OmegaW2 (1-eps)^2)`,
  `R_target = Lambda0 OmegaW2 (1-eps_eta)(1-eps)^2 / (ZW (1+chi0)^2)`.

- **M2 — Support-blindness (NON-VACUOUS, per reframed F1):** define `Mtr = Mmix(1 + zeta(1-eps)/(1-zeta eps))` (notes §4). (a) NEGATIVE CONTROL: show `D[Mtr, zeta]` and `D[Mtr, Mmix]` are NOT zero (support channel is live). (b) LEAK DETECTOR: form a contaminated observable `RtrLeak = Rtr Mtr/Mmix` and show `D[RtrLeak, zeta] =!= 0` (the check discriminates). (c) CLAIM: show `D[obs, zeta] = D[obs, Mmix] = 0` for `obs ∈ {R_tr, T^2, eps_eta, q_tr, q_nt, q_eta}`, and/or `FreeQ[obs, zeta] && FreeQ[obs, Mmix]`. Do NOT route `Mtr` into the observables (the notes define only the reduced observables; no in-stage cancellation form exists).

- **M3 — Finite packet factorization + rigid reduction:**
  `q_nt + (B*/C*) q_tr = ln(T^2 / T_ref^2)`, using `R_target_ref · T_ref^2 = Lambda0(1 - eps_eta_ref)`; and on `q_tr = 0`: `q_nt = ln(T^2/T_ref^2)`, `q_eta = ln(eps_eta/eps_eta_ref)`. Here
  `q_tr = -C* ln(R_tr/R_tr_ref)`,
  `q_nt = B* ln(R_tr/R_tr_ref) + ln((1-eps_eta)/(1-eps_eta_ref)) - ln(R_target/R_target_ref)`,
  `q_eta = ln(eps_eta/eps_eta_ref)`.

- **M4 — First-order drift compilers** (three closed forms, notes §3):
  `δln R_tr = -chi0 deltaU/((1+chi0)(1+deltaU)(1+chi0+deltaU)) [(1+deltaU) dlnchi + (1+chi0) dlndelta]`,
  `δln T^2 = dlnZW - dlnOm2 + 2chi0/(1+chi0) dlnchi + 2eps/(1-eps) dlneps`,
  `δln R_target = dlnOm2 - dlnZW - 2chi0/(1+chi0) dlnchi - 2eps/(1-eps) dlneps - eps_eta/(1-eps_eta) dlnepseta`.

- **M5 — First-order packet relation + rigid slice (per F2):**
  `q_nt^(1) + (B*/C*) q_tr^(1) = δln T^2` with `q_tr^(1) = -C* δln R_tr`, `q_nt^(1) = B* δln R_tr - δln R_target - eps_eta/(1-eps_eta) dlnepseta`; and on `δln R_tr = 0`: `q_nt^(1) = δln T^2`, `q_eta^(1) = dlnepseta`.

- **M6 — Tracking gate (NON-VACUOUS, per F3):** anchored to the derived `δln R_tr`, show `δln R_tr` factors as `-prefactor · tracking_condition` with `tracking_condition = (1+deltaU) dlnchi + (1+chi0) dlndelta`, and that substituting `tracking_condition = 0` into the derived `δln R_tr` gives 0.

- **M7 — Finite gates:** `q_nt = 0 ⟺ T^2 = T_ref^2` (via `exp(q_nt_rigid) = T^2/T_ref^2`) and `q_eta = 0 ⟺ eps_eta = eps_eta_ref` (via `exp(q_eta) = eps_eta/eps_eta_ref`).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 238` and confirms each manifest item M1-M7 has a corresponding `[ok]` check, that the support-blindness (M2) and tracking-gate (M6) checks route the support symbols / gate condition through a non-trivial substitution (not vacuous), and that the script exits 0. Also `redteam exec-sympy 238` confirms F1-F3 landed (rewritten lines reference live symbols / `_calc` quantities) and the SymPy script still exits 0.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage238_physical_branch_transfer_shape_compiler_packet_factorization_and_post_static_dressing_invariance_theorem_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M7 with reframed support-blindness negative controls, leak detection, and derived tracking-gate checks.
- deviation: none
