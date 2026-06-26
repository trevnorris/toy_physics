# Directive pathA_31 — Gate 2: Scalar Breathing (η₀₀) → legacy (a,L) collective closure

**STATUS:** ✅ DONE — executed 2026-06-25, verdict **`BREATHING_CALIBRATED`** (tri-review CLEAN after a
remediation: v1 was REJECTED for HF `x−x` + hardcoded counterfactual + a chosen-to-pass truncation threshold;
v2/v3 fixed those, the Codex review also caught a Parseval-muddled measure; a final surgical fix made the
STRUCTURE gate computed-not-hardcoded). The distributed lift's ℓ=0 reduction genuinely reproduces the legacy
`(a,L)` collective closure: profiles derived (`L₀`-harmonic lifting), `M_AB/K_AB` from operator projection
(not `∂²E_geom`), two genuinely-independent HF routes match, structure gate computed + able-to-fail probe
fires, dual-engine agreed (max Δ 5.7e-13). **CALIBRATED** because the wall coefficients `μ_η,T_w,K_η` are
calibration inputs. **Verified by:** the generalized-eigenproblem `V₂`-overlap truncation is clean at the
physical `β` (`β_L0=1.85 → o_1=0.993`) and across the whole order-unity stiffness window (`K_η/T_w≲2.6`),
failing only for sharp walls (the `β_L0=18` v1 silently used → `o_1=0.466`). **Caveat (adversarial):** the
overlap measure does not guard profile-correctness (it passes any mouth-BC shape) — the anti-wrong-profile
teeth come from the recomputed HF mismatch; the clean truncation is demonstrated *for an assumed order-unity
wall stiffness*. Review trail: Codex→GLM→Codex (v1) → reject → v2/v3 remediation (Codex) → execute →
tri-review (arbiter + fidelity + adversarial) → structure-gate fix → re-verify.
**Push:** the full nonlinear moving-throat PDE closure (memory `project-moving-throat-pde-push`;
`research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`). This is **Gate 2**, building on
Gate 1 (`pathA_30`, frozen-wall D/N unit test = `DN_UNITTEST_BC_DEPENDENT`).

**One-line goal.** Switch on **only the scalar breathing channel** `η₀₀(w,t)` (the ℓ=0 sector) of the
distributed wall action, and show that its **two-mode collective truncation** reproduces the legacy
finite-dimensional `(a,L)` closure — with the breathing kept **genuinely dynamical** — so the distributed
Hellmann–Feynman force reduces to the legacy `F_a, F_L`. (Scaffold §11.2; handoff §8.1–8.2, §5.3, §6.2.)

**VERDICT LADDER (the script must compute, not assert, which rung fires; `BREATHING_CALIBRATED` is the DEFAULT
— `PASS` must be *earned* per §2.5):**
- `BREATHING_PASS` — the ℓ=0 modal operator is derived and **dynamical** (inertia term kept); the two-mode
  truncation is a **consistent** reduction (modes separated from the continuum / projection closes to stated
  error); the reduced `M_AB, K_AB` (Q=(δa,δL)) come out with legacy **structure** (M symmetric pos-def, K
  symmetric, the `(a,L)` Hessian structure of `E_geom`); the distributed HF force reduces to `F_a, F_L`;
  **AND** the mode profiles `α_a(w), α_L(w)` are **DERIVED from the Gate-1 straight reference throat**
  `R₀(w)`, not free-calibrated; counterfactual rejects a wrong truncation.
- `BREATHING_CALIBRATED` *(default)* — everything above **except** the profile derivation: `α_a, α_L` (and/or
  the stiffnesses `μ_η, T_w, K_η`) are **calibration inputs**, not derived from `R₀(w)`. Honest toy-model
  outcome: the `(a,L)` recovery is structural + calibrated, not first-principles. Records which inputs are
  calibration. Does **not** block the push.
- `BREATHING_FAIL_<reason>` — `FAIL_NOT_DYNAMICAL` (M and K are not both produced as projections of the one
  ℓ=0 wall operator — i.e. K was built from the static `E_geom` Hessian and inertia merely appended, §2.1);
  `FAIL_TRUNCATION_INCONSISTENT` (the two lowest physical modes' `V₂`-overlap `o_k < 1−ε_trunc` for all
  physical `β`, or `min(ω₁²,ω₂²)≤0` (unstable), so the 2-mode `(a,L)` truncation is uncontrolled, §2.2); `FAIL_STRUCTURE` (M not
  pos-def, or K's Hessian *pattern* — off-diagonal sign / rank — incompatible with legacy `E_geom`, §2.3);
  `FAIL_HF_FORCE` (the *independently* computed distributed HF force ≠ legacy `F_a, F_L`, §2.4);
  `FAIL_COUNTERFACTUAL` (a wrong profile still "matches" with all calibrations frozen, §3);
  `FAIL_ENGINE_DISAGREE` (the two engines' `FullSimplify` disagree on a derived expression — distinct from
  ENV_BLOCKED).
- `BREATHING_ENV_BLOCKED` — missing Mathematica / license / timeout / infrastructure (NOT an algebra
  disagreement). **Not a physics verdict** — re-run.

---

## §0. Scope, framing, and what is able-to-fail (vs what is definitional)

**Toy-model contract (standing).** Calibrations allowed but every claim *earned*, not asserted; goal = a
legitimate brane/bulk/soliton structure supporting EM+GR. (memories `project-calibrated-pde-goal`,
`feedback-calibrate-predict-methodology`.)

**Nothing is static (keystone) — the sharpest trap for THIS gate.** The legacy `(a,L)` sector has a *static*
reading (`E_geom(a,L)` with equilibria). That reading is **forbidden** here. `(a,L)` must emerge as the **slow
/ lowest collective modes of the DYNAMICAL ℓ=0 wave equation on the throat**
`μ_η q_{0,tt} − ∂_w(T_w ∂_w q_0) + K_η q_0 = S_0^{(ψ)} + f_0^{ext}` (handoff §8.1) — i.e. inertia (`μ_η q_tt`)
is retained and the `(a,L)` closure is its lowest-mode truncation, NOT a frozen energy-minimization. Dropping
the time-derivative / posing `(a,L)` as static equilibria ⇒ `FAIL_NOT_DYNAMICAL`. (memory
`project-model-mechanics-corrections` §1.)

**The tautology to avoid.** `(a,L)` are *defined* as the lowest moments of `η` (handoff §5.2, §6.2), so
"recovering `(a,L)`" is **partly definitional** and is NOT by itself the test. The genuinely able-to-fail
content is everything that is NOT guaranteed by the definitions:
1. the ℓ=0 reduction keeps **inertia** (dynamical, not static);
2. the two-mode truncation is **consistent** (the two collective modes separate from the continuum / the
   projection closes — a wrong reference throat would give no clean separation);
3. the reduced `M_AB, K_AB` have legacy **structure** (M pos-def, K symmetric, matching the `E_geom` Hessian
   shape) — not just any 2×2;
4. the distributed Hellmann–Feynman force `f_0^{(ψ)} = −δH_ψ/δΣ|_0` genuinely **reduces** to the legacy
   `F_a=−∫ρ∂_aV_conf`, `F_L=−∫ρ∂_LV_conf` (handoff §2.3, scaffold §5.3);
5. (for `PASS`) the profiles `α_a, α_L` are **derived from `R₀(w)`**, not chosen to fit.

**Necessary-not-sufficient.** This gate validates the monopole/geometry-on coupling before any grouped-`P2`
(ℓ=2) work (Gate 3). It does not touch the quadrupole, Maxwell, or leak/return sectors (§7).

---

## §1. The object to build

Start from the **already-frozen** distributed wall action (handoff §8.1) — do not re-postulate it:
```
S_η^(2) = ½ ∫ dt dw dΩ √γ₀ [ μ_η(∂_t η)² − T_w(∂_w η)² − T_Ω η(−Δ_S²)η − K_η η² ]
```
Restrict to the **ℓ=0 (axisymmetric breathing) channel**: `η = η₀(w,t) Y₀₀(Ω)`, `Y₀₀=1/(2√π)`
(the `T_Ω` term drops since `l(l+1)=0`). This gives the ℓ=0 modal equation (handoff §8.1) — **a genuine wave
equation in `(w,t)`**, with the wall-drive source `S_0^{(ψ)}` from `δV_conf = −(V'_wall(Σ₀/ℓ_c)/ℓ_c) η`
(handoff §5.3).

**Two-mode collective truncation** (handoff §8.2, §6.2):
```
η₀(w,t) = 2√π [ α_a(w) δa(t) + α_L(w) δL(t) ],     Q^A = (δa, δL)
```
yielding the reduced Lagrangian `L_red^(0) = ½ M_AB Q̇^A Q̇^B − ½ K_AB Q^A Q^B` with
```
M_AB = 4π ∫ dw μ_η α_A α_B ,    K_AB = 4π ∫ dw [ T_w α_A' α_B' + K_η α_A α_B ].
```
(Notation: handoff §8.2 writes `K_0` for the ℓ=0 stiffness; `K_0 ≡ K_η` here, since the `T_Ω` angular term
drops at `l(l+1)=0`.)

**Legacy target** (handoff §2.3): `L_geom = ½ M_a ȧ² + ½ M_L L̇² − E_geom(a,L)` with force laws
`F_a=−∂H/∂a, F_L=−∂H/∂L`. The reduced `(M_AB, K_AB)` must reproduce this **structure**; the legacy magnitudes
`(M_a, M_L, ∂²E_geom)` are calibration unless `α_a, α_L` are derived from `R₀(w)`.

**Emit (don't assume):** the ℓ=0 modal operator, the truncation integrals `M_AB, K_AB`, the HF force
reduction, and the provenance of `α_a, α_L`.

---

## §2. Able-to-fail sub-checks (each independently computed)

1. **DYNAMICAL (solve-order enforced — the anti-static-relapse core).** Both `M_AB` **and** `K_AB` must be
   **produced as projections of the one ℓ=0 wall operator** onto the same profiles
   (`M_AB=4π∫μ_η α_Aα_B`, `K_AB=4π∫[T_w α_A'α_B' + K_η α_Aα_B]`), so inertia and stiffness come from the SAME
   operator and cannot be chosen independently. **Mandatory ORDER (emit each step, in this order):** (i) the
   ℓ=0 modal operator + BCs + the `μ_η`-weighted inner product → (ii) the profiles `α_A` + the projection
   integrals → (iii) `M_AB, K_AB` + the dynamical EOM `M_AB Q̈ + K_AB Q = (HF force)` → (iv) **only THEN**
   compare to legacy `E_geom`. Building `K_AB` from the static `E_geom` Hessian and merely *appending*
   `M_AB Q̈` is the static relapse ⇒ `FAIL_NOT_DYNAMICAL`. Emit the EOM with `Q̈` present **and** the
   provenance of `K_AB` (it must be the operator projection of step ii–iii, NOT `∂²E_geom`).
2. **TRUNCATION CONSISTENCY (decidable — separate the two BC problems).** Define the ℓ=0 operator
   `L₀ = μ_η⁻¹[−∂_w(T_w∂_w·) + K_η·]` on `w∈[0,L₀]`. **Two distinct BC problems (state both explicitly):**
   - the **`g`-lane eigenvalue problem** uses the *homogeneous* wall BCs (mouth `g(0)=0`, regular cap
     `T_w g'(L₀)=0`); `L₀` is self-adjoint w.r.t. `⟨f,g⟩ = 4π∫μ_η f g dw` on this homogeneous domain; its
     spectrum `{λ_n}` is **discrete**;
   - the **collective modes `α_a, α_L`** carry the *inhomogeneous* boundary data (`α_a(0)=1, α_L(0)=0`, regular
     cap) and are therefore **NOT** eigenmodes of the homogeneous problem — they span the 2D subspace `V₂` that
     the homogeneous eigenbasis cannot represent (this is the boundary-carrying / lifting-function part).
   The discarded `g(w,t)` lane is the component of `η₀` orthogonal to `V₂`.
   **Use the rigorous Galerkin truncation-validity measure (NOT a sum of raw basis overlaps — those are
   Parseval-complete representation coefficients, ≈1 for *any* function, and are NOT leakage; that conflation
   was the round-1-remediation error).** Build the **combined-basis** mass and stiffness matrices
   `𝕄_{IJ}=⟨b_I,b_J⟩_μ`, `𝕂_{IJ}=4π∫[T_w b_I'b_J' + K_η b_I b_J]dw` over `B={α_a, α_L, g_1, …, g_N}` (the
   2D `V₂` plus the lowest `N` homogeneous `g`-modes; emit `N` and show convergence in `N`). Solve the
   **generalized eigenproblem** `𝕂 v = ω² 𝕄 v`. Take the two lowest *true* eigenvectors `v₁, v₂` and compute
   their **`𝕄`-metric overlap with `V₂`**: `o_k = ‖P_{V₂} v_k‖_𝕄 / ‖v_k‖_𝕄` (`P_{V₂}` = the `𝕄`-orthogonal
   projector onto `V₂`). The 2-mode `(a,L)` truncation is consistent **iff the two lowest physical modes live
   in `V₂`**: `o₁, o₂ ≥ 1 − ε_trunc` with **FIXED `ε_trunc = 0.1`** (low modes ≥90% captured). This is
   Parseval-clean (it asks whether the low physical modes are in `V₂`, not whether `α_A` expands in a complete
   basis) and genuinely able-to-fail (if the low modes are mostly `g`-lane, `o_k < 0.9` ⇒
   `FAIL_TRUNCATION_INCONSISTENT`). Context: require **`min(ω₁², ω₂²) > 0`** — i.e. BOTH lowest physical
   eigenvalues positive (throat stability; ascending ordering makes a bare `λ₂>0` insufficient since `λ₁`
   could be ≤0) — and report the gap to the next mode `(ω₃²−ω₂²)/ω₂²`.
   - **`β` is NOT a free pass-knob (anti-gaming, mandatory).** The profile sharpness is set by
     `β = √(K_η/T_w)` — **a scalar only on the constant-coefficient straight reference** (state this
     assumption; for `w`-dependent `T_w,K_η` it is a parameter family). Large `β` sharply localizes `α_A`. The
     script must **SWEEP `β`** over a **predeclared** range — anchored to `beta_from_R0` (the `β` implied by
     the Gate-1 reference throat geometry) if derivable, else a broad fixed log sweep — and **emit the full
     `o_k(β)` table**, reporting whether a defensible `β`-window gives `o_k ≥ 1−ε_trunc` (disclose it) and the
     value at `beta_from_R0`. The verdict may **NOT** be obtained by silently picking the single `β` that
     passes. If `o_k < 1−ε_trunc` for all physical `β`, the honest verdict is `FAIL_TRUNCATION_INCONSISTENT` —
     a **legitimate** falsification-first result, not something to rescue.
   **Honest scope:** this is a quasi-static modal-overlap check — necessary, not a guarantee that the source
   `S₀^{(ψ)}` doesn't dynamically excite `g` (deferred, §7).
3. **STRUCTURE (specify the invariant pattern + allowed freedom).** Emit `M_AB, K_AB`; verify `M_AB` symmetric
   **positive-definite** and `K_AB` symmetric. **Emit the FULL legacy `E_geom` Hessian** (not only the
   `E_ratio` penalty part — `E_geom` may carry additional structure) and compare `K_AB` to it by its
   **invariant PATTERN**, not element values: the off-diagonal **sign** and the rank / zero-structure. Concrete
   reference for the penalty contribution: `E_ratio=½κ(L−αa)²` has Hessian `[[κα², −κα],[−κα, κ]]` — a
   rank-1-plus structure with a **negative** off-diagonal (for `α>0`); require the derived `K_AB` to match the
   full-Hessian sign/rank pattern. **Allowed calibration freedom = the overall scales** of `μ_η, T_w, K_η` (and the legacy
   `M_a, M_L, κ, α`). If matching requires fitting the **full** `K_AB` element-by-element (more freedom than
   those scales), that **downgrades the verdict to `CALIBRATED`**; a genuine sign/rank-pattern mismatch ⇒
   `FAIL_STRUCTURE`.
4. **HF FORCE REDUCTION (TWO genuinely independent computations — no `x−x`).** Compute the two sides by routes
   that share **NO intermediate expression**, then compare:
   - (a) **distributed:** `f_Σ = −δH_ψ/δΣ|_0` projected onto `(δa,δL)` via the **action / configuration-space
     measure** `∫ f_Σ α_A √γ₀ dw dΩ` — the force pairs with the displacement in the *action* measure, **NOT**
     the `μ_η`-weighted mass inner product (unless the projected source is explicitly `μ_η⁻¹ S`); state the
     measure;
   - (b) **legacy parametric:** `F_a = −∫ρ ∂_a V_conf(X;a,L)`, `F_L = −∫ρ ∂_L V_conf(X;a,L)` computed from the
     `V_conf(X;a,L)` parametrization **directly** — its own `∂_a, ∂_L` and its own integral — **NOT** by
     reusing the distributed integrand.
   **FORBIDDEN (the round-1 defect):** assigning the *same* expression to both sides (round-1 had `F_dist` and
   `F_legacy` byte-identical, so `differences=0` was a tautological `x−x`). The script must **emit BOTH
   unsimplified expressions** so the audit can confirm two different routes; the equality must be a genuine
   simplification of two distinct forms. (`α_A ≡ ∂η/∂Q` chain-rule identity also remains forbidden.) A genuine
   mismatch in sign, measure, or magnitude ⇒ `FAIL_HF_FORCE`.
5. **PROFILE PROVENANCE (selects `PASS` vs `CALIBRATED`; dynamical provenance, decidable from an artifact).**
   For `PASS`, the script must emit a **derivation formula** for `α_a(w), α_L(w)` as the **quasi-static
   response of the ℓ=0 operator `L₀`** to the inhomogeneous boundary drive — the harmonic lifting function
   solving `L₀ α_A = 0` in the interior with the inhomogeneous BCs (`α_a(0)=1, α_L(0)=0`, regular cap tied to
   `R₀(L₀)=0`, from the §6.2 moment definitions `q_00(0)=2√π δa`) — **NOT** an arbitrary static-deformation
   shape disconnected from `L₀`. Tying the profile to the *operator* is what makes §2.2's spectral check
   meaningful and is the honest sense in which `(a,L)` is the **slow limit of the dynamical system**: the
   profile is the operator's `ω→0` response, while the dynamical content lives in the genuine `μ_η` inertia
   `M_AB`. (`K_AB` equalling the static Hessian is the *desired* legacy match, **not** a relapse — the relapse
   would be losing the `μ_η` inertia, or using a profile unrelated to `L₀`.) A bare label "derived", or a
   static-deformation shape not derived from `L₀`, ⇒ default `CALIBRATED`. Likewise flag `μ_η, T_w, K_η` as
   `derived|calibration`.
6. **STATIC↔DYNAMIC consistency (procedural).** Confirm the static `(a,L)` equilibria (`∂E_geom/∂Q=0`) are the
   `Q̈→0` limit of the dynamical EOM of check 1 — the one-object-two-limits principle — emitted as a limit of
   the *same* reduced system, not a separate static solve. (Honest scope: like Gate 1's static control, this
   verifies limit-consistency, not an independent able-to-fail surface.)

---

## §3. Mandatory counterfactual self-test

The script **must** prove a wrong reduction **fails its own gate** (the pathA_29/pathA_30 lesson), with **all
flags COMPUTED by re-running the projection + the independent HF computation on each mutated profile — NOT
typed literals** (the round-1 defect: `wrong_M_det`, `hf_force_mismatch` were hardcoded `"0"`/`True`).
**Freeze ALL baseline calibrations** (`μ_η, T_w, K_η` and the legacy targets `M_a, M_L, κ, α`) and `α_L`, and
**forbid any re-fit**. Test **TWO** mutated `α_a`:
- (i) **degenerate** `α̃_a = 0` (trivially breaks `M_det→0`);
- (ii) **non-trivial** wrong shape that still satisfies the mouth BC `α_a(0)=1` but is wrong (e.g. `α̃_a = 1`
  constant, or a wrong-decay profile). This one will **NOT** lose `M_AB` positive-definiteness, so the gate
  must catch it via the **recomputed** HF-force mismatch and/or a drop in the truncation overlap `o_k` —
  proving the gate is not toothless against realistic wrong profiles.
For BOTH mutations emit the **recomputed** quantities (from the mutated profile, not typed): `wrong_M_det`,
`wrong_o_k`, and `wrong_F_a` (distributed, recomputed) vs the frozen legacy `F_a`. Emit
`counterfactual_guard: {degenerate: {M_det, fails: true}, nontrivial: {profile, M_posdef, o_k, F_a_dist,
F_a_legacy_frozen, hf_mismatch, fails: true}, calibrations_frozen: true, alpha_L_frozen: true,
refit_forbidden: true}`. If the **non-trivial** wrong profile still "matches" (recomputed HF agrees AND
`o_k ≥ 1−ε_trunc`) ⇒ `FAIL_COUNTERFACTUAL`.

**Tell to watch for:** the genuine method (the `∫dw` projection) must run on the headline path, not only in a
side check; if `M_AB, K_AB` are produced from a real projection somewhere but the *verdict* compares typed-in
legacy values, that is the hardcode one level deeper.

**Tell to watch for:** if `M_AB, K_AB` are typed-in to equal `(M_a, M_L, ∂²E_geom)` and read back, that is the
hardcode. They must be **produced** from the `∫dw` integrals over the actual profiles, with the integrals
emitted as artifacts.

---

## §4. Anti-tautology firewall + reduction certificate

- **FORBIDDEN:** writing the legacy `(M_a, M_L, E_geom)` values into `M_AB, K_AB` by hand and reading them
  back; asserting the `(a,L)` recovery without producing the projection integrals.
- **EMIT SOLVE ARTIFACTS:** the ℓ=0 modal operator, the eigenmodes/profiles `α_a, α_L` (with their
  provenance), the `∫dw` integrands and the resulting `M_AB, K_AB`, the HF-force reduction chain, the
  dynamical EOM with the `Q̈` term, and the counterfactual's broken match.
- **Reduction certificate (mandatory, emitted):** the ℓ=0 restriction (show the `T_Ω` term drops at `l=0`);
  the background (Gate-1 straight reference `R₀(w)`, `ρ₀=ρ_*`); the source `δV_conf` linearization; and the
  **provenance of every input** (`α_a, α_L, μ_η, T_w, K_η`) marked `derived|calibration`. Carry the **same
  phonon-limit caveat as Gate 1** if any matter-sector quantity (`c_s`, BdG `k⁴`) enters: deferred under
  `kξ≪1`, `ξ=ħ/(mc_s)`.
- **Controls must fire from computed inputs** — no literal `match=True`, no `x−x`.

---

## §5. Dimensional consistency

Units-restored homogeneity check (sympy.physics.units **or** explicit `{M,L,T}` dimension-tuple arithmetic;
must FAIL if a dimension is perturbed): confirm `M_AB` has dimension of (mass-like inertia) consistent with
`½M_AB Q̇Q̇` being an energy, `K_AB` a stiffness, and the HF force `F_a` a force. (memory
`feedback-dimensional-consistency-check`.)

## §6. Deliverables

- `software/stage1_solver/tools/pathA_31_scalar_breathing_sympy.py` + `..._scalar_breathing.wl` (dual-engine,
  genuinely independent routes; `timeout 600`; exit 0; `engine_agreement` via `FullSimplify`).
- `software/stage1_solver/reports/pathA_31_scalar_breathing.md` — LINE 1 = verdict rung; then the modal
  operator, the profiles + provenance, `M_AB/K_AB` with integrands, the HF-force reduction, the dynamical EOM,
  the static-limit consistency, the counterfactual, `engine_agreement`, `dim_check`.
- `software/stage1_solver/reports/pathA_31_results.yaml` — machine-readable: `verdict`;
  `reduction_certificate`; `dynamical: true`; `truncation`: `{o_1, o_2, N_modes, convergence_in_N,
  o_k_of_beta` (the full β-sweep table)`, beta_window, beta_from_R0, lambda2, gap, epsilon_trunc: 0.1}`;
  `M_AB`, `K_AB` (derived); `M_posdef`; `K_structure_ok`; `hf_force_reduces` (+ both
  unsimplified HF expressions `F_dist`, `F_legacy`); `profile_provenance: derived|calibration`;
  `stiffness_provenance`; `static_limit_consistent`; `counterfactual_guard` (degenerate + non-trivial, all
  recomputed); `engine_agreement`; `dim_check: pass`. **YAML/markdown only — no JSON.**
- A short feed note under `research/pde_ledger/notes/stages/` against scaffold §11.2 (do **not** edit audited
  ledger stages).

## §7. Deferred / out of scope (record, do not silently drop)

- the grouped real `P2` (ℓ=2) sector + isotropy gate (Gate 3); the quadrupole normalization (Gate 4);
- the localized Maxwell mixed channels `A_w/J^w/F_{μw}` (keep alive as ontology, don't excite here);
- the brane↔bulk leak/return + the radiative (imaginary) part of the response (later gates);
- a full first-principles derivation of `μ_η, T_w, K_η` from the matter sector (calibration inputs here unless
  the certificate derives them).

## §8. Process

- **Codex codes + runs both engines; Claude reviews only.** Directive states requirements + acceptance; does
  not pre-design the script route beyond "independent routes" + "emit artifacts".
- **Dual-engine** (Mathematica: Codex needs `--sandbox danger-full-access`, or orchestrator runs `math` as
  arbiter; ≤2 seats).
- **Review ordering:** iterate Codex (xhigh) to GREEN → one GLM pass → fold → Codex confirm → execute.
  **Run the review gauntlet via a gauntlet-runner agent** to protect orchestrator context
  (`feedback-offload-review-gauntlet`). Execution backgrounded, never `timeout`-wrapped.
- **Post-exec tri-review** (clean agents): orchestrator arbiter re-run + transliteration-fidelity +
  adversarial ("is the recovery EARNED, or definitional/hardcoded?"). Then user gate.
