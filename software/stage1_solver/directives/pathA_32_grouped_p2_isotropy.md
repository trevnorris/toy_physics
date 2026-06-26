# Directive pathA_32 — Gate 3: Grouped-P2 (ℓ=2) sector → isotropy gate a₂=b₂=a₄=b₄=0

**STATUS:** ✅ DONE — executed 2026-06-26, verdict **`ISOTROPY_CALIBRATED`** (full tri-review CLEAN after a
remediation: v1 was REJECTED for two pass-by-construction defects the adversarial-with-ablation leg caught and
fidelity missed — (a) dual-engine gaming (the Mathematica baseline lane-collapse was a vacuous `x−x` on three
byte-identical copies, honest only on the probes) and (b) 5 of 8 counterfactual probes typed their FAIL booleans
plus a tautological `able_to_fail` aggregate; v2 fixed both — honest per-harmonic `.wl` assembly + per-lane raw-`D`
cross-engine comparison, each probe computed-from-mutated-input with a self-ablation — and re-verified). The
distributed lift's ℓ=2 grouped-`P2` sector satisfies the isotropy / lane-degeneracy theorem at the linearized
isotropic reference: the three grouped lanes `{20,21,22}` collapse to common conservative coefficients (raw-`D`
lane defects = 0, the **PRIMARY** gate), the reduction is SO(3)-covariant (angular Gram = `I₅`; computed `−Δ_S²`
eigenvalue `λ_m=6` per harmonic; the angular stiffness `K₂=∫[T_w β₂'² + (K_η+6T_Ω)β₂²]` — the term that dropped at
ℓ=0 — is now alive), and the gate is genuinely able-to-fail (8 counterfactual probes, each flips under ablation).
**This is the gate where the quadrupole sector first appears.** **CALIBRATED** because the wall constants
(`μ_η, T_w, K_η, T_Ω, β₂`) and the symbolic radial scalars are calibration inputs (toy-model contract).
**Two-tier gate:** the raw-`D` lane collapse is the verdict-bearing PRIMARY test; the published
`a₂=b₂=a₄=b₄=0` (normalized `u`-defects) is a necessary-but-not-sufficient cross-check (a per-lane
order-independent prefactor cancels in the normalized response). **Verified by:** arbiter byte-for-byte re-run;
fidelity + adversarial re-ablation EARNED. Review trail: Codex→GLM→Codex (v1) → tri-review REJECT → v2 remediation
(Codex) → execute → tri-review (arbiter + fidelity + adversarial-with-ablation). Gate 3 of the moving-throat PDE
completion ladder (`research/pde_ledger/notes/stages/moving_throat_pde_completion_ladder.md`), building on Gate 1
(`pathA_30`, frozen-wall D/N unit test = `DN_UNITTEST_BC_DEPENDENT`) and Gate 2 (`pathA_31`, scalar breathing =
`BREATHING_CALIBRATED`). **NEXT = Gate 4** = the `54/5` quadrupole normalization `m̂₀²P₀ = 54Gc_s⁵/5a⁵c⁵` (handoff
§8.3/§12; `G` is CALIBRATED — the PDE delivers the FORM/branch + outgoing `P₀`, not `G`). Push memory:
`project-moving-throat-pde-push`.

**One-line goal.** Switch on the **grouped real-P2 channels** `{20, 21c, 21s, 22c, 22s}` (the ℓ=2 sector) of
the distributed wall action, reduce them to the **three grouped lanes** `{20, 21, 22}`, and show that on an
**isotropic reference throat** the lanes **collapse to common conservative coefficients**
(`D_{20,n}=D_{21,n}=D_{22,n}` — the **primary** raw-`D` gate; the normalized-defect form `a₂=b₂=a₄=b₄=0` is the
published cross-check), AND that this gate is **genuinely able-to-fail** (an injected anisotropy in the reference
**splits the raw-`D` lanes**; under an order-dependent anisotropy the normalized `a₂/b₂` also move).
(Scaffold §8.2; handoff §6.1, §8.3, §11.)

**What this gate is NOT.** It is *not* the `54/5` quadrupole normalization (`m̂₀²P₀ = 54Gc_s⁵/5a⁵c⁵`, handoff
§8.3/§12 = **Gate 4**), *not* the outgoing odd-coefficient `N_{A,n}` extraction (Gate 4/5), and *not* the
fully-solved nonlinear branch data `(K_A,M_A,B_{A,n},Z_{A,n},N_{A,n})` (handoff §13.2 = **Gate 6, the wall**).
Gate 3 is the **SO(3)-covariance / lane-degeneracy theorem at the linearized isotropic reference level** — the
exact analogue of Gate 2 reproducing the legacy `(a,L)` closure, one harmonic family up.

**VERDICT LADDER (the script must COMPUTE, not assert, which rung fires; `ISOTROPY_CALIBRATED` is the DEFAULT —
`PASS` must be *earned* per §2.5):**
- `ISOTROPY_PASS` — lane collapse holds (`D_{20,n}=D_{21,n}=D_{22,n}` at every order `n`, at the finest
  resolved granularity per §2.2) ⇒ `a₂=b₂=a₄=b₄=0` to tolerance; the reduction is **SO(3)-covariant** (§2.1)
  and **dynamical** (§2.3); the lanes are **stable/well-posed** (§2.4); **all** counterfactuals fire (§3);
  **AND** the radial profile `β₂(w)` is **DERIVED** from the Gate-1 straight reference throat `R₀(w)` (the ℓ=2
  support equation) and the angular stiffness `T_Ω` is derived — not free-calibrated.
- `ISOTROPY_CALIBRATED` *(default)* — everything above **except** the derivation: `β₂(w)`, `T_Ω`, and/or the
  wall stiffnesses `μ_η, T_w, K_η` are **calibration inputs**, not derived from `R₀(w)`. Honest toy-model
  outcome: the isotropy/lane-collapse is **structural + covariant + able-to-fail**, on a calibrated reference.
  Records which inputs are calibration. **Does not block the push.**
- `ISOTROPY_FAIL_<reason>`:
  - `FAIL_NOT_COVARIANT` — the reduction is not SO(3)-covariant on the isotropic reference: the angular
    Laplacian eigenvalue is **not** `l(l+1)=6` for all five real ℓ=2 harmonics in the stated normalization, or
    the angular self-overlaps `∫|Y^real_2m|²dΩ` differ across `m` on the isotropic measure (§2.1) — the lanes
    split even on the isotropic reference for a *machinery* reason.
  - `FAIL_ANISOTROPIC_BRANCH` — the lanes do **not** collapse on the intended isotropic reference despite a
    covariant reduction: any **raw-`D` defect** `a_{D,n}`/`b_{D,n}` ≠ 0 (primary), or an ungrouped c/s-pair split
    `D_{21c,n}≠D_{21s,n}` / `D_{22c,n}≠D_{22s,n}`, or a per-port split, beyond tolerance ⇒ the reference branch is
    genuinely anisotropic (a real physics result that blocks Gate 4). (Note: the normalized `u`-defects can stay
    zero under a pure per-lane prefactor split — §1c — so raw-`D` is the signal here.)
  - `FAIL_TAUTOLOGICAL` — the three lanes were **not** assembled from genuinely distinct angular projections
    (shared/copied lane datum, hardcoded equality, byte-identical lane assembly) ⇒ collapse is a construction
    artifact, not a result (§4 firewall fires).
  - `FAIL_NOT_ABLE_TO_FAIL` — an injected-anisotropy or m-dependent-profile counterfactual (§3) does **not**
    flip the gate to FAIL ⇒ the isotropy control is toothless.
  - `FAIL_STABILITY` — `M₂ ≤ 0` or `K₂ ≤ 0` (including the `6T_Ω` contribution) anywhere in the physical
    calibration window ⇒ unphysical lane (§2.4).
  - `FAIL_SINGULAR_RESPONSE` — `D_{A,0} = K_A − B_{A,0} − Z_{A,0} = 0` (or `u_2/u_4` unbounded) for any lane in
    the physical calibration window despite `M₂,K₂>0` ⇒ the normalized response is ill-defined (§2.4).
  - `FAIL_ENGINE_DISAGREE` — the two engines' `FullSimplify`/simplify disagree on a derived expression
    (distinct from `ENV_BLOCKED`).
- `ISOTROPY_ENV_BLOCKED` — missing Mathematica / license / timeout / infrastructure (NOT an algebra
  disagreement). **Not a physics verdict** — re-run.

---

## §0. Scope, framing, and what is able-to-fail (vs what is definitional)

**Toy-model contract (standing).** Calibrations allowed but every claim *earned*, not asserted; goal = a
legitimate brane/bulk/soliton structure supporting EM+GR (memories `project-calibrated-pde-goal`,
`feedback-calibrate-predict-methodology`). `β₂`/`T_Ω`/wall-stiffness being calibration inputs is fine and is
the *expected* `ISOTROPY_CALIBRATED` outcome — but the **covariance, lane-collapse, and able-to-fail content
must be genuinely computed**.

**Nothing is static (keystone) — the trap for THIS gate.** The conservative isotropy data `D_{A,n}` are the
**even low-frequency Taylor coefficients of the genuinely ω-dependent grouped response** `D_A^cons(ω) =
D_{A,0} + D_{A,2}ω² + D_{A,4}ω⁴ + O(ω⁶)`, with inertia retained in `D_{A,2} = −(M_A + B_{A,2} + Z_{A,2})`. The
isotropy result is the **ω→0 sector of that dynamic object**, NOT a frozen angular-energy Hessian. Posing the
ℓ=2 sector as a static energy minimization (dropping the `ω`-dependence / `M_A`) ⇒ `FAIL` (and is conceptually
the recurring relapse, memory `project-model-mechanics-corrections` §1; `feedback-dimensional-consistency-check`
companion). Add the static↔dynamic consistency reading: `D_{A,0}` (the ω→0 conservative datum) and the finite-ω
expansion must be the **two limits of one object**.

**The tautology to avoid (load-bearing for this gate).** "The ℓ=2 multiplet is degenerate under rotations" is a
**symmetry truth**, so naively it can be *built in* by feeding three identical computations and reading back
raw-`D` lane equality (and its `a₂=b₂=0` cross-check). That is theatre. The genuinely able-to-fail content is
that the **reduction machinery actually respects the symmetry through honest, independently-evaluated angular
integrals** of *different* functions (`Y^real_20 ∝ P₂(cosθ)`, `Y^real_22c ∝ sin²θ cos2φ`, …) against the
reference angular weight — and that the **same machinery breaks** (the **raw-`D` lanes split**; under an
order-dependent anisotropy `a₂` or `b₂` also moves) the moment the reference is made anisotropic. The collapse is
a *result of integrating distinct functions against an isotropic weight*, not three copies of one number (§4
enforces this).

**Definitional vs earned.** The grouping `{21}={21c,21s}`, `{22}={22c,22s}` and the weighted-trace/defect
algebra (`x̄=(x₂₀+2x₂₁+2x₂₂)/5`, `a_x=(2x₂₀−x₂₁−x₂₂)/10`, `b_x=(x₂₁−x₂₂)/2`) are **definitions** (handoff §11)
— not the test. The earned content is §2.1–§2.5 + §3.

---

## §1. The object to build

Work at the **linearized isotropic reference** level (the regime of Gates 1–2): a stationary straight finite
throat `R₀(w)`, wall motion linearized, the conservative low-frequency response assembled per harmonic channel.

**(1a) Grouped-P2 distributed wall reduction (handoff §6.1, §8.3).** Decompose the wall displacement
`η(Ω,w,t) = η₀(w,t)Y₀₀ + Σ_{m∈P2(real)} q_{2m}(w,t) Y^real_{2m}(Ω) + η_{≥3}`, with the grouped real-P2 set
`{20, 21c, 21s, 22c, 22s}`. For one grouped real-P2 component `η_{2m} = β₂(w) q_{2m}(t) Y^real_{2m}(Ω)`, the
reduced Lagrangian is `L_{2m} = ½ M₂ q̇_{2m}² − ½ K₂ q_{2m}²`, with
- `M₂ = ∫ dw μ_η β₂²`,
- `K₂ = ∫ dw [ T_w (β₂')² + (K_η + l(l+1) T_Ω) β₂² ]`, with **`l(l+1)=6` for ℓ=2** (the angular `−Δ_S²`
  term that **dropped at ℓ=0** in Gate 2; this is the new physics).

**(1b) Per-lane conservative response (handoff §11; scaffold §8.2).** For each lane `A ∈ {20, 21, 22}` build
`D_A^cons(ω) = D_{A,0} + D_{A,2}ω² + D_{A,4}ω⁴`, with
- `D_{A,0} = K_A − B_{A,0} − Z_{A,0}`,
- `D_{A,2} = −(M_A + B_{A,2} + Z_{A,2})`,
- `D_{A,4} = −(B_{A,4} + Z_{A,4})`,

where `B_{A,n}` (BdG support moments) and `Z_{A,n}` (Maxwell/mixed moments) are the linearized support data.
**On the isotropic linearized reference, every lane-dependent quantity factorizes** as `(angular self-overlap of
lane A) × (common radial/spectral scalar)`: the support moments `B_{A,n}=C_A·B̃_n`, `Z_{A,n}=C_A·Z̃_n`, **and the
wall sector** `M_A=C_A·M₂`, `K_A=C_A·K₂` — where `C_A=∫ w(Ω)|Y^real_{2m}|² dΩ` is the lane's angular self-overlap
(the `K_η` and `T_Ω`/`λ_m` parts of `K_A` share the **same** `C_A`, since `∫ w Y_{2m}(−Δ_{S²})Y_{2m}dΩ = λ_m C_A`).
**On the isotropic reference with orthonormal harmonics (Gram `=I₅`, §2.1), `C_A=1` for all lanes**, so
`M_A=M₂`, `K_A=K₂`, `B_{A,n}=B̃_n`, `Z_{A,n}=Z̃_n` — lane-independent ⇒ raw lane collapse. The **common radial
scalars `B̃_n, Z̃_n` (and `K̃, M̃`) are carried symbolically** (they are the still-unsolved nonlinear-branch
content, deferred to Gate 6, §7) — Gate 3 proves the **lane-equality of the prefactors**, not their numerical
value. The **angular projections MUST be computed as honest integrals** of the explicit real harmonics against
the reference angular weight `w(Ω)` (isotropic: `w≡1`; anisotropic probe: §3a — where a nonzero `ε` makes
`C_A(ε)` lane-dependent and the collapse breaks).

**(1c) Normalized response + isotropy variables (handoff §11).**
`Y_A(ω) = D_{A,0}/D_A^cons(ω) = 1 + u_2^{(A)}ω² + u_4^{(A)}ω⁴ + O(ω⁶)`, with
- `u_2^{(A)} = −D_{A,2}/D_{A,0}`,
- `u_4^{(A)} = (D_{A,2}² − D_{A,0}D_{A,4})/D_{A,0}²`.

Apply the defect algebra to the triples `{u_2^{(20)},u_2^{(21)},u_2^{(22)}}` and `{u_4^{(A)}}`:
- `a₂ = (2u_2^{(20)} − u_2^{(21)} − u_2^{(22)})/10`, `b₂ = (u_2^{(21)} − u_2^{(22)})/2`,
- `a₄, b₄` analogously from `u_4^{(A)}`.

**The isotropy gate (two tiers — NOT equivalent).** The **PRIMARY, strong** gate is **raw lane collapse**
`D_{20,n}=D_{21,n}=D_{22,n}=D_n` at every order `n∈{0,2,4}` — apply the defect algebra **directly to the raw
`D_{A,n}` triples**: `a_{D,n}=(2D_{20,n}−D_{21,n}−D_{22,n})/10`, `b_{D,n}=(D_{21,n}−D_{22,n})/2`, require all
`=0`. The **SECONDARY** gate is the handoff §11 published form on the normalized response,
`a₂=b₂=a₄=b₄=0` (from the `u`-defects). ⚠️ **These are NOT equivalent:** raw collapse ⟹ `u`-defects vanish, but
**not conversely** — a per-lane **order-independent** prefactor `g_A` (i.e. `D_{A,n}=g_A d_n` with the *same*
`g_A` for all `n`) **cancels** in `Y_A=D_{A,0}/D_A(ω)`, leaving `u₂,u₄` lane-independent (defects zero) while the
raw `D` lanes split by `g_A`. So the `u`-defects are **necessary but not sufficient**; the **raw-`D` collapse is
the verdict-bearing isotropy gate**, the `u`-defects are the published cross-check. Both are emitted; the §3a
probe must move the **raw-`D`** defects (and an order-dependent variant must additionally move the `u`-defects —
§3a).

---

## §2. Able-to-fail sub-checks (each independently computed)

**§2.1 COVARIANCE (the new ℓ=2 physics).** Compute, do not assert:
(i) the angular Laplacian eigenvalue **computed per harmonic** from `−Δ_{S²} Y^real_{2m} = λ_m Y^real_{2m}`
(symbolically, not asserted) equals `λ_m=6` for **all five** real ℓ=2 harmonics; any `λ_m≠6` ⇒
`FAIL_NOT_COVARIANT`. **The `K₂` angular coefficient MUST be this computed `λ_m`** — the `(K_η+λ_m T_Ω)` term
uses the per-harmonic computed eigenvalue, **NOT a typed `6`** — and the covariance check **asserts the
coefficient used in `K₂` equals the computed `−Δ_{S²}` eigenvalue per harmonic**; a mismatch ⇒
`FAIL_NOT_COVARIANT`. (This is what gives §3d teeth: a wrong typed value disagrees with the computed eigenvalue
and **flips the verdict**, not merely 'shifts diagnostics'.)
(ii) the **5×5 angular Gram matrix** `G_{mm'}=∫ Y^real_{2m} Y^real_{2m'} dΩ` on the isotropic measure equals the
**identity `I₅` exactly** — i.e. the real ℓ=2 basis is genuinely ortho**normal**, not merely orthogonal with a
common `C≠1`. (Equal self-overlaps alone are insufficient: the `M₂/K₂` formulas assume unit normalization, so a
common `C` would silently rescale `M_A,K_A,D_A`.) Emit the full `G` and assert `G=I₅`; **or**, if a non-unit
normalization is deliberately used, **carry the common Gram factor explicitly** through `M_A,K_A,D_A` and emit
it. These are the machinery-level statement that the reduction is SO(3)-covariant; failure here is distinct from
a genuinely anisotropic branch (§2.2).

**§2.2 LANE COLLAPSE (the isotropy gate).** Assemble each lane's `D_{A,n}` **independently** from its own
honest angular projection (§1b) against the isotropic reference weight, then verify `D_{20,n}=D_{21,n}=D_{22,n}`
at **every order `n∈{0,2,4}`**. **Ungrouped-first rule (anti-masking by premature grouping):** FIRST emit and
verdict-test the **five ungrouped** per-channel `D_{2m,n}` for `m∈{20,21c,21s,22c,22s}` — in particular the
intra-lane c/s pairs must agree (`D_{21c,n}=D_{21s,n}` and `D_{22c,n}=D_{22s,n}` at every order `n`) on the
isotropic reference; **only then** form the grouped lanes `{21}`,`{22}`. A c/s split that is averaged away by
grouping ⇒ `FAIL_ANISOTROPIC_BRANCH` (the grouped `D_{A,n}` must NOT be the first place collapse is tested).
**Granularity rule (anti-masking):** if `B_{A,n}/Z_{A,n}` are assembled as a
sum over ports/modes/sectors `r` (handoff §10), verify lane-equality **per-port (before summation)**, not only
on the summed `D_{A,n}` — summation can hide compensating port-level anisotropy (flag as `FAIL_ANISOTROPIC_BRANCH`
if per-port equality fails even when the sum is equal). **Primary verdict-bearing quantities = the RAW-`D` defects** `a_{D,n}=(2D_{20,n}−D_{21,n}−D_{22,n})/10`,
`b_{D,n}=(D_{21,n}−D_{22,n})/2` for each `n∈{0,2,4}` — these detect a pure per-lane prefactor split that the
normalized `u`-defects miss (§1c). Emit: the raw `D_{A,n}` triples + the six raw-`D` defects
`{a_{D,0},b_{D,0},a_{D,2},b_{D,2},a_{D,4},b_{D,4}}` (**PRIMARY**), and the normalized `{a₂,b₂,a₄,b₄}` + the
triples `{u_2^{(A)}}`, `{u_4^{(A)}}` (published §11 cross-check), all with magnitudes. Tolerance: prefer **exact
symbolic 0**; if numeric, the magnitudes must be `< 1e-10` AND `≥ 10` orders below the same quantities under the
§3a anisotropy probe.

**§2.3 NOT-STATIC / DYNAMICAL.** Confirm `D_A^cons(ω)` is the genuine ω-expansion with `M_A` retained in
`D_{A,2}` (not appended to a static Hessian); confirm `D_{A,0}` is the ω→0 limit of that one object (static↔
dynamic consistency). Building `K₂`/`D_{A,0}` from a frozen angular energy without the `M₂ q̇²` inertia ⇒ FAIL.

**§2.4 STABILITY / WELL-POSEDNESS.** `M₂ > 0` and `K₂ > 0` (with the `+6T_Ω` term) across the physical
calibration window; emit the lane frequency `Ω_{2m} = √(K₂/M₂)` and confirm it is **m-degenerate** on the
isotropic reference (a second, independent witness of lane collapse). `M₂≤0` or `K₂≤0` ⇒ `FAIL_STABILITY`.
**Denominator guard:** the normalized variables `u_2^{(A)}, u_4^{(A)}` require `D_{A,0} ≠ 0`; a non-degenerate
wall (`M₂,K₂>0`) can still have `D_{A,0} = K_A − B_{A,0} − Z_{A,0} = 0`. Emit `D_{A,0}` per lane and confirm it
is nonzero/bounded across the calibration window; `D_{A,0}=0` (or unbounded `u`) ⇒ `FAIL_SINGULAR_RESPONSE` (a
distinct well-posedness rung — the `0/0` in `u_2` must be caught, not silently swallowed).

**§2.5 PASS-vs-CALIBRATED gate.** Emit `ISOTROPY_PASS` **only if ALL non-deferred wall inputs are derived** —
i.e. `β₂(w)` is **derived** from the Gate-1 reference `R₀(w)` (the ℓ=2 radial support equation, the
`(K_η+6T_Ω)` analogue of the Gate-2 `β`-profile), `T_Ω` is derived, **and** the wall stiffnesses `μ_η, T_w, K_η`
are derived (not free-calibrated). If **any** of these (`β₂`, `T_Ω`, `μ_η`, `T_w`, `K_η`) remains a calibration
input, emit `ISOTROPY_CALIBRATED` and **list every calibration input** (including the symbolic radial scalars).
**Expectation for this gate:** `ISOTROPY_CALIBRATED` is the realistic Gate-3 outcome with `μ_η, T_w, K_η` (and
likely `β₂, T_Ω`) as calibration inputs — consistent with the verdict-ladder default. The which-rung decision
must be computed from a `derived_inputs` / `calibration_inputs` partition in the script, not typed.

---

## §3. Mandatory counterfactual self-test (able-to-fail probes baked into the script)

Each probe must show the gate's verdict **changes**; a probe that leaves the verdict unchanged ⇒
`FAIL_NOT_ABLE_TO_FAIL`. (Lesson: Gate-2 v1's counterfactual was a hardcoded `True` on a degenerate mutation —
memory `feedback-negative-verdict-short-circuit`. These must compute from real perturbed inputs.)

**§3a ANISOTROPY INJECTION (the decisive probe — TWO sub-probes that together demonstrate the two-tier gate).**

*§3a-i — pure-prefactor anisotropy (shows WHY raw-`D` must be primary).* Apply a common quadrupole angular weight
`w(Ω)=1 + ε·P₂(cosθ_n̂)` (axis n̂) **uniformly to every angular self-overlap** — the wall sector (`M_A,K_A`) AND
the support sectors (`B_{A,n},Z_{A,n}`) — so each lane scales by its **single** self-overlap
`C_A(ε)=∫ w(Ω)|Y^real_{2m}|² dΩ = 1 + ε c_A`. Then `D_{A,n}=C_A(ε)·d_n` (one prefactor, all orders). Assert:
(a) the **raw-`D` defects** `a_{D,n}`/`b_{D,n}` go nonzero and **scale linearly with `ε`** (the PRIMARY gate flips
to FAIL — emit ≥3 values of `ε`); **and** (b) the **normalized `u`-defects** `a₂`/`b₂` **stay zero** (the `C_A`
cancels in `Y_A`). This is the positive demonstration that the `u`-defects are necessary-but-not-sufficient and
the **raw-`D` defects are the load-bearing isotropy signal** (§1c).

*§3a-ii — order-dependent (sector-selective) anisotropy (shows the `u`-defects are not dead).* Apply the angular
weight `w(Ω)=1+ε·P₂(cosθ_n̂)` to the **support sectors only** (`B_{A,n},Z_{A,n}` → `C_A(ε)·B̃_n, C_A(ε)·Z̃_n`)
while leaving the **wall sector isotropic** (`M_A=M₂, K_A=K₂`, unscaled) — physically, an anisotropic
matter/EM background on an intrinsically-isotropic wall. Then `D_{A,0}=K₂−C_A(ε)(B̃_0+Z̃_0)` and
`D_{A,2}=−(M₂+C_A(ε)(B̃_2+Z̃_2))` acquire **different** lane dependence (the `C_A(ε)` no longer factors out of the
ratio), defeating the single-`g_A` cancellation **without perturbing the angular eigenvalue** (so this does NOT
trip the §2.1/§3d covariance check — it is a genuine anisotropic-branch probe, distinct from a covariance break).
Assert **both** the raw-`D` defects **and** the normalized `u`-defects `a₂`/`b₂` go nonzero and scale with `ε` —
proving the `u`-defect cross-check is genuinely live under non-prefactor anisotropy.

A probe that leaves its asserted defect family unchanged ⇒ `FAIL_NOT_ABLE_TO_FAIL`. (§3a-i is the one that
makes the raw-`D`/`u` distinction concrete; §3a-ii prevents the `u`-defect cross-check from being vacuous.)

**§3b M-DEPENDENT PROFILE.** Set `β₂ → β_{2,m}` differing across `m` (e.g. `β_{22} ≠ β_{20}`), recompute; assert
the **raw-`D` lanes split** (`D_{20,n}≠D_{21,n}` or `D_{21,n}≠D_{22,n}` at some order `n`) ⇒ FAIL. (A pure
per-lane profile *scale* can give `D_{A,n}=g_A d_n`, so the normalized `u`-defects `a₂/b₂` may stay zero even
though raw `D` splits — the raw-`D` defects are therefore the primary verdict-bearing signal here, consistent
with §1c/§2.2.) Only additionally require `a₂` or `b₂` ≠ 0 if the `β` perturbation is explicitly
order-dependent / shape-changing (not a single per-lane scale). Guards against a reduction accidentally blind to
the channel.

**§3c DEGENERATE.** `β₂ ≡ 0` ⇒ `M₂=K₂=0` ⇒ the lane is ill-posed; assert this is caught as
degenerate/`FAIL_STABILITY`, **not** a spurious pass (the `0/0` in `u_2` must be detected, not silently
swallowed).

**§3d WRONG-EIGENVALUE.** Force a mismatch between the `K₂` angular coefficient and the computed `−Δ_{S²}`
eigenvalue (e.g. type `2`, the ℓ=1 value, while §2.1 computes `λ_m=6`); assert the §2.1 covariance check **flips
the verdict to `FAIL_NOT_COVARIANT`** — NOT merely "diagnostics shift." (A *uniform* `6→2` preserves lane
degeneracy and may preserve stability, so the teeth come from the coefficient-vs-computed-eigenvalue **equality**
check in §2.1, not from any stability/degeneracy shift.)

---

## §4. Anti-tautology firewall + reduction certificate

**Firewall (the lanes must be genuinely distinct computations).** The script MUST:
1. **Emit the explicit real harmonics** `Y^real_{20}, Y^real_{21c}, Y^real_{21s}, Y^real_{22c}, Y^real_{22s}`
   (closed form in θ,φ) and the **full angular overlap matrix** (`5×5` Gram against the reference weight, and
   the `3×3` grouped-lane reduction) so the reviewer sees the collapse arises from integrating **different
   functions** against an isotropic weight.
2. Assemble each lane's `D_{A,n}` from its **own** harmonic projection — **forbidden:** hardcoding
   `D_{20,n}=D_{21,n}`; computing one lane and copying it to the others; feeding the lanes from identical input
   provenance (this is exactly the Gate-2 v1 `x−x` defect, one harmonic up). A **shared assembler routine called
   with genuinely different harmonic inputs is correct and desirable** — the forbidden thing is copied *data* or
   identical input provenance, NOT shared *code machinery*. To prove distinctness, emit a distinct
   harmonic/input hash per lane plus the per-lane integral traces; if the lane inputs are not provably distinct
   (matching hashes / copied data), emit `FAIL_TAUTOLOGICAL`.
3. The lane-collapse booleans and `{a₂,b₂,a₄,b₄}` must **compute from the assembled matrices**, carrying the
   §3 able-to-fail probes (not typed literals).

**Reduction certificate (state precisely).** FROZEN-INPUT: the linearized isotropic reference background
`R₀(w)`; `β₂(w)`; the wall stiffnesses `μ_η,T_w,K_η,T_Ω`; the **symbolic** common radial scalars
`B̃_n, Z̃_n, K̃, M̃`. COMPUTED: the explicit harmonics + the `5×5` angular Gram (asserted `=I₅`); the per-lane
`M₂,K₂,D_{A,0/2/4}`; the **raw-`D` defects `{a_{D,n},b_{D,n}}` (primary)** and the normalized `{u_2^{(A)}},
{u_4^{(A)}}, {a₂,b₂,a₄,b₄}` (cross-check); the §3 perturbed responses; the §2.4 frequencies.
DEFERRED (§7): numerical `B_{A,n}/Z_{A,n}` from the solved nonlinear branch (Gate 6); the `54/5` normalization
(Gate 4); the outgoing `N_{A,n}` (Gate 4/5). **The conservative claim:** *given* an isotropic linearized
reference, the grouped-P2 reduction is SO(3)-covariant ⇒ **raw lane collapse `D_{20,n}=D_{21,n}=D_{22,n}`** (⟹
`a₂=b₂=a₄=b₄=0`), **and** this is able-to-fail under anisotropy (§3a). No claim is made about the *value* of the
common `D_n` (that needs the nonlinear solve).

---

## §5. Dimensional consistency (standing pre-numbers step)

**Restore units** — do not trust `c_s=ħ=m=a=1` pins (natural units hide dropped Jacobians/`a⁵` factors —
memory `feedback-dimensional-consistency-check`). With SymPy dimensional homogeneity confirm: `[M₂]` from
`μ_η·length`, `[K₂]` consistent including the `6T_Ω` term (so `T_Ω` carries the right dimension to add to `K_η`),
`[D_{A,n}]` consistent across `n` with the `ω^{2}`/`ω^{4}` weighting, `[u_2]=ω^{-2}`, `[u_4]=ω^{-4}`, and
`a₂,b₂` sharing `[u_2]`, `a₄,b₄` sharing `[u_4]` (so the defect algebra is dimensionally legal). The angular
integrals are dimensionless (solid angle). Emit a homogeneity table.

---

## §6. Deliverables

**Dual-engine required** (memory `feedback-dual-engine-required`: a `.wl` wherever Mathematica *can*
independently verify — here it can: angular integrals, the eigenvalue check, the ω-expansion, the defect
algebra).
- `software/stage1_solver/tools/pathA_32_grouped_p2_isotropy_sympy.py`
- `software/stage1_solver/tools/pathA_32_grouped_p2_isotropy.wl`
- `software/stage1_solver/reports/pathA_32_grouped_p2_isotropy.md` (verdict + which-rung + narrative)
- `software/stage1_solver/reports/pathA_32_results.yaml` (machine-readable; **YAML not JSON** —
  `feedback-no-json-for-llm-io`) emitting, at minimum:
  - the five explicit `Y^real_{2m}` and the `5×5` angular Gram matrix with the `G=I₅` assertion boolean (§2.1),
    plus the `3×3` grouped reduction (isotropic + §3a perturbed);
  - the computed `−Δ_{S²}` eigenvalue `λ_m` per harmonic AND the boolean that the `K₂` angular coefficient equals
    the computed `λ_m` (§2.1, the §3d hook);
  - the **five ungrouped** per-channel `D_{2m,n}` for `m∈{20,21c,21s,22c,22s}` and the intra-lane c/s-pair
    equality booleans (`D_{21c,n}=D_{21s,n}`, `D_{22c,n}=D_{22s,n}`) — emitted BEFORE grouping (§2.2);
  - per-lane `M₂, K₂, D_{A,0}, D_{A,2}, D_{A,4}` (with the `D_{A,0}≠0` denominator-guard boolean, §2.4);
  - the **raw `D_{A,n}` triples + the six raw-`D` defects `{a_{D,0},b_{D,0},a_{D,2},b_{D,2},a_{D,4},b_{D,4}}`
    (PRIMARY verdict-bearing)** + magnitudes (§2.2);
  - the triples `{u_2^{(A)}}`, `{u_4^{(A)}}`;
  - `{a₂, b₂, a₄, b₄}` + magnitudes (published §11 cross-check); the lane-collapse booleans (computed,
    able-to-fail);
  - `Ω_{2m} = √(K₂/M₂)` + m-degeneracy flag (§2.4);
  - the §3 counterfactual results: §3a perturbed **raw-`D` defects** vs `ε` (≥3 points, scaling confirmed) AND
    the order-dependent variant's perturbed `{a₂,b₂}` vs `ε`; §3b m-dependent-profile split; §3c degenerate
    catch; §3d wrong-eigenvalue → `FAIL_NOT_COVARIANT` verdict flip — each with the verdict it produced;
  - `derived_inputs` / `calibration_inputs` partition (§2.5);
  - `verdict`, `which_rung`, `engine_agreement` (max symbolic/numeric delta; tolerance **`< 1e-10` for
    symbolic-exact checks, `< 1e-8` for numeric** — `FAIL_ENGINE_DISAGREE` if exceeded), `dim_homogeneity_table`.

A pde_ledger feed note `research/pde_ledger/notes/stages/moving_throat_pde_pathA_32_grouped_p2_isotropy_result.md`
(Codex-authored, Claude-verified) summarizing the verdict for the ledger.

---

## §7. Deferred / out of scope (record, do not silently drop)

- **The `54/5` quadrupole normalization** `m̂₀²P₀ = 54Gc_s⁵/5a⁵c⁵` and the outgoing prefactors `P₀,P₂,P₄`
  (handoff §8.3/§12) — **Gate 4**. (`G` is calibrated; the PDE delivers the FORM/branch — memory
  `project-pathA-build`.)
- **The outgoing odd-coefficient `N_{A,n}` extraction** and the `Ŷ₂^out(ω)` fingerprint (handoff §12.1) —
  **Gate 4/5**.
- **The solved nonlinear branch data** `(K_A,M_A,B_{A,n},Z_{A,n},N_{A,n})` on the true nonlinear solution
  (handoff §13.2) — **Gate 6 (the wall)**. Gate 3 carries these as symbolic common scalars.
- **Cross-ℓ consistency** (ℓ=0/1 return ↔ ℓ=2 quadrupole, the reconciliation's new gate) — **Gate 5**.
- **Higher harmonics** `η_{≥3}` (ℓ≥3) — out of scope.
- **BC provenance** — inherits Gate 1's `BC_DEPENDENT` (the D/N assignment is a labeled calibration input).
- **PN match-back** (`research/4d_*pn*`) — the decisive falsifier, downstream (Phase 4); do not re-derive the
  audited PN ladder (memory `project-pn-gravity-ladder`).

---

## §8. Process

1. **Codex design-review** of THIS directive (`codex exec -c model_reasoning_effort=xhigh`, backgrounded, never
   `timeout`-wrapped) — leading-question / can't-fail-gate / tautology trap hunt (memory
   `feedback-directive-design-review`). Fold fixes; Codex confirm-pass to GREEN.
2. **GLM** single fresh-perspective pass (user-run, out-of-band). Fold; Codex iterate to GREEN again
   (`feedback-review-ordering-codex-then-glm`).
3. **Codex executes** dual-engine (Claude reviews, Codex codes — `feedback-claude-reviews-codex-codes`); Codex
   iterates each script to exit 0 (`feedback-codex-iterates-until-clean`). Mathematica run needs
   `--sandbox danger-full-access` (`feedback-codex-mathematica-sandbox`); ≤2 concurrent `math -script`
   (`feedback-mathematica-single-seat`); scripts wrapped `timeout 600` (`feedback-script-timeout-policy`).
4. **Tri-review** (clean agents): (a) orchestrator arbiter re-run refreshing committed outputs; (b)
   transliteration-fidelity audit (code-vs-equations, term-by-term); (c) **adversarial-with-ablation** — force
   each verdict-bearing flag/threshold and confirm the verdict actually changes; verify the FIX too, not just
   the claim (both directions: false-PASS *and* false-FAIL — memory `feedback-negative-verdict-short-circuit`).
5. **Run the gauntlet via per-gate AGENTS**, not inline (`feedback-offload-review-gauntlet`); orchestrator keeps
   two carve-outs: verify flagged math folds, and one final directive read pre-exec.
6. User gate before commit (`feedback-sequential-audit-chunks`; commit only when the user asks).
