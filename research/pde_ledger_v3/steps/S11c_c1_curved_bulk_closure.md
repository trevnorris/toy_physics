# S11c-c1 — curved-interface bulk closure: DtN/impedance operator + permeable face response (step record)

Slug `S11c_c1_curved_bulk_closure`. Step record for **S11c-c1**, the first half of the S11c-c curved-interface
bulk closure (the S11c-c decision-list row was split c1/c2 by user choice 2026-09-03). c1 solves the perturbed
curved two-face outgoing bulk acoustic problem for the nonlocal **DtN/impedance** operator + the **permeable
face response** `(δp_s,J_s,t_s)(V_s,μ_θ)`, the Fredholm noninvertibility loci, and the three-object dissipation
audit; it **exports** the closed face response + the DtN operator/kernel for c2. Physics authority
`directives/S11c_c1_SHARED_PHYSICS.md` (spec committed `db5cbf88`). It generalizes S11b's flat-face B0b (impedance
`Z`, three regimes) and B0c (permeable response) to the tilted faces S11c-a shape-differentiated. The self-energy
fold back into the S11c-b slab operator is **c2** (not this record).

⭐ This record is the **interpretation** layer. The two engines PRINT objects and state no conclusions; every result
names the computed object and the commit/record behind it, and every quantitative claim about a run carries the
command in `directives/_measurements/` (rule 2).

> ⭐⭐ **STATUS OF THE CLOSE (2026-09-05; scoped after 2 step-record review legs — Grok + Codex sol — corrected an
> earlier "closure-wide AGREE" overreach).** Unlike S11c-b, S11c-c1's **cross-engine residual WAS RUN** (family-scoped,
> 46 of 50 families on 30 GB). **What is EARNED cross-engine:** the **core two-momentum DtN kernel AGREES** — all 4
> (anchoring×face) joined cases collapse to **EXACT ZERO** off-diagonal under physically-justified identifications,
> and the collapse is load-bearing (a wrong jet sign or freezing the two legs into one leaves it nonzero); the
> **face-response coefficients AGREE** for LAB_HELD (off-diagonal) and MATERIAL_ADVECTED on the physical diagonal
> `k=k′`. **No cross-engine DISAGREEMENT was found among anything measured.** ⛔ **What is NOT earned — a closure-wide
> AGREE:** seal 5 (background density) stays a **surfaced rule-17 freeze (UNDECIDED)**; the **raw DTN_OPERATOR** form,
> the **ENERGY** audit, and an off-diagonal **flat-resolvent leg-labeling** (MATERIAL) are **UNDECIDED**; the **4 giant
> families** (PERMEABLE_PORT_HERMITIAN, PERMEABLE_DISSIPATION_VS_OMEGA_TAU, UNIFORM_LIMIT_S11CC1_OPERAND,
> UNIFORM_LIMIT_RESIDUAL) + the **full per-family symbolic residual** are **UNMEASURED — DEFERRED** (≥64 GB,
> `DEFERRED_HEAVY_RUNS.md`), ⛔ not pre-adjudicated. Reconcile record + reproducible bridge:
> `_measurements/S11c_c1_comparator_reconcile.md`, `_measurements/S11c_c1_reconcile_reproduce.py`; leg reports
> `_measurements/S11c_c1_step_record_review_{grok,codex_sol}.txt`. ⇒ **the kernel is cross-engine closed; the rest is
> UNDECIDED/deferred.**

## What the step computes
On the S11c curved background (in-plane-varying thickness `W_bg(y)=W̄₀[1+η w₁(ξ)]`, two anchorings
LAB_HELD/MATERIAL_ADVECTED, two face-normal signs ±1, two density representatives ρ4/ρbr, three radiation regimes
PROPAGATING/EVANESCENT/GRAZING per leg): (1) the **nonlocal DtN/impedance operator** `Z_s(ω;k,k′)` — composition
`N₀∘M_{h_s}∘N₀ + Div(h∇) + κ²h` — AND its **two-momentum kernel** with BOTH branch legs `q_out(k), q_out(k′)`
explicit, `Z_1 ∝ Ŵ_bg(k−k′)·[q(k)q(k′)+k·k′−ω²/c_s0²]/(q(k)q(k′))` (⛔ no single-`k` multiplier, ⛔ no one-leg
left-quantization — both delete the mode mixing, the rule-17 freeze); (2) the **permeable curved face response**
`(δp_s,J_s,t_s)(V_s,μ_θ)` as an operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹` (not a scalar division), μ_θ a separate
opaque constitutive operand; (3) the **Fredholm noninvertibility** condition (formal) + the algebraic locus
protocol (flat/finite-dim); (4) the **three-object dissipation audit** — bulk-radiation Hermitian part `H_a[Z]`,
the two-port permeable-port Hermitian form, and the INDEPENDENT traction-vs-far-field-Poynting energy balance.

## The two blind engines + the comparator
- **SymPy engine** `scripts/S11c_c1_bulk_closure_sympy_audit.py` (`PY_S11CC1_*`). Imports the frozen base
  `scripts/S11c_b_exports.py` via `load_model` over the fold (44-root `IMPORT_KEYS`), binding only its declared
  keys. Migration `f90e7630`; build baseline `65afa1cd`; repair (5 emit-only controls now bite) CLEARED `d6e16471`.
- **Wolfram engine** `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` (`WL_S11CC1_*`), **blind** — imports
  nothing, re-derives the §§1–2 setup + S11c-a substrate + the S11c-b μ_θ operand from the sibling specs
  (`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control). Build `e139bc61`; repair R1–R4
  `13f0bd2c`; repair-2 (Fredholm re-freeze + dead parity axis) `dd34d564`. Blindness is the cross-engine control:
  an agreement is independent construction, not a copy.
- **Comparator** `scripts/S11c_c1_cross_engine_comparator.py` — the frozen-T7 join (contract N8): joins by object
  name on the S11c-a axis-typed keys, pairs residual operands, three-valued, rejects a native boolean, PRINTS and
  decides nothing (rule 2). Gated directive `84686a54` (2 decision legs, 12 folds); reviewed baseline `7141e6ad`
  (2 re-review legs CLEAN); scoped repair `704308af` (astra, +14 lines: R1 `cS0←c_s0` fold, R2 held-parse the WL
  4-arg triple-range energy integrals; 2 re-review legs CLEAN). Canonical `.out` committed `4a14100a` (GIN).

## The arc (each result: commit/record + how verified)
- **Both engines per-engine SOUND (2-leg each, form-ablated).** SymPy: core physics 2-leg-confirmed (two-momentum
  DtN kernel — Grok's "tangential freeze" adjudicated a FALSE POSITIVE; operator-inverse response; Λ_X-only
  traction; opaque μ_θ; delta topology) + 5 emit-only control defects repaired to bite. Wolfram: core 2-leg-SOUND
  (two-momentum kernel both legs, rigid-shift, operator inverse, T-a re-derivation, §5 controls bite, blindness,
  names, μ_θ) + 3 MUST-FIX repaired (composition momentum labels → both legs; real far-field Poynting; response
  `t_s` binding) + repair-2 (Fredholm two-leg re-freeze; dead parity axis) — all 2-leg CLEARED. ⚠ The WL repair
  directive INITIALLY skipped its rule-7 decision legs (a real gap, remediated: run retroactively it found the
  directive not sound; repair-2 got its legs BEFORE the build and caught a first-draft R2 that would have corrupted
  the correct PERMEABLE_PORT_HERMITIAN congruence). Records `_measurements/S11c_c1_{build_review,repair_directive,
  wl_build_review,wl_repair_directive,wl_remediation_plan}.md`.
- **Comparator BUILT + re-reviewed SOUND** (`7141e6ad` → repair `704308af`). 50 shared non-LOCAL families; every
  SEAL proven load-bearing by name-map ablation (two-momentum `qOut[...]` stays applied, μ_θ opaque, ω-assumption,
  bg-density field, regime/parity all UNRECONCILED — surfaced, not folded); raw-operand whitelist exactly
  `{DTN_OPERATOR.WHOLE_OBJECT, NONINVERTIBILITY_CONDITION.OPERATOR}`; three-valued residual preserved; no false
  agreement, no hidden disagreement. Records `_measurements/S11c_c1_comparator_{build,repair}_directive.md`.
- **CROSS-ENGINE RECONCILE — the two engines AGREE** (`_measurements/S11c_c1_comparator_reconcile.md`, this session).
  46 of 50 families run family-scoped (light 42 + mid 4), peak RSS ~317 MB, EXIT=0. The comparator prints the RAW
  residual (no identification applied, rule 2) so joined residuals are nonzero; OUR adjudication applies the known,
  independently-justified representational identifications and tests collapse (reproducible from the committed
  `.out`: `_measurements/S11c_c1_reconcile_reproduce.py`). Verdicts on the five pre-identified seals: **(1)
  two-momentum q_out — AGREE** (DtN kernel → EXACT ZERO, all 4 (face×anchoring) cases; the Stage-2 residual
  collapsed to the radiation-branch dispersion relation itself, →0 on-shell; adversarial corruptions bite); **(2) ω
  real-assumption — AGREE** (unified, contributes nothing); **(3) μ_θ — AGREE, scoped** (opaque/supplied both sides;
  response coefficients `ε·A≡B` to machine zero for LAB_HELD and for MATERIAL on the physical diagonal `k=k′`; ε
  sits in PY's μ_θ `.py:840` vs WL's coefficient, the coeff constructor `closed_coefficients` `.py:723-755` omits
  ε); **(4) regime/parity naming — AGREE on content** (`OUTPUT_X`/`INPUT_X` vs bare is a bridgeable keying
  convention; PARITY_MATRIX joins and IS the proven kernel); **(5) background density (rule 17) — UNDECIDED,
  SURFACED FREEZE** ⛔ NOT AGREE (the functional dependence agrees and neither engine differentiates ρ — 0
  derivatives — but a PY bare constant is not globally a WL live field `rhoBrBgRho4Constant(x)=(ρbr/W₀)WBg(x)`;
  folding one to the other is the rule-17 collapse this rebuild exists to catch — kept surfaced; c2 re-adjudication
  MANDATORY). **Also surfaced UNDECIDED:** the raw DTN_OPERATOR face-form; the off-diagonal flat-resolvent
  leg-labeling (PY output-leg `q_out` vs WL input-leg for the flat half-space symbol — diagonal-exact, differs
  off-diagonal in the MATERIAL response coefficients); the ENERGY audit (PY closed-form vs WL far-field integral).

## Established (cross-engine AGREE) vs owed (surfaced/deferred)
- **ESTABLISHED — per-engine SOUND (2-leg each) AND cross-engine AGREE:** the **two-momentum DtN KERNEL** (all 4
  cases, exact zero, load-bearing) and the **permeable face-response COEFFICIENTS** `(δp_s,J_s,t_s)` (LAB_HELD
  off-diagonal; MATERIAL diagonal), the flat symbol, the parity matrix (= the kernel), the degenerate loci, the
  dimensions. ⛔ **The raw DtN OPERATOR (whole-object noncommutative form) is NOT in this list — it is UNDECIDED**
  (§ below): kernel-AGREE does not extend to the operator. The Hermitian/reactive dissipation parts AGREE only at
  c1's first order (WL's true-area-weighted adjoint `BoundaryMeasure` → PY's plain conjugate-swap at first order;
  O(η²)→S11c-e; full symbolic DEFERRED). **For c2's consume-set** (`dtn_operator`, `dtn_flat_symbol`, `dtn_kernel`,
  face response — per-engine-verified 44-row delta over `S11c_b_exports.py`): the **kernel + flat symbol + response
  are cross-engine AGREE**; the **`dtn_operator` whole-form and seal-5 density are UNDECIDED — c2 must NOT treat them
  as cross-engine-closed** (re-adjudicate in c2, see What's next).
- **OWED / SURFACED (UNDECIDED, on this box) — ⛔ not AGREE, ⛔ no predicted residual:** seal 5 (background density
  field-vs-constant freeze); the raw **DTN_OPERATOR** face-form (PY per-face, WL per-branch); the **ENERGY** audit
  (PY closed-form propagating energy vs WL unevaluated far-field flux integral); the off-diagonal **flat-resolvent
  leg-labeling** (MATERIAL). `HOMOGENEITY_*` is a keying gap (j=0, PY `HEIGHT_NORMAL_NORMAL` vs WL `HOMOGENEITY`),
  ⛔ NOT AGREE-by-inheritance.
- **DEFERRED to a ≥64 GB box (`DEFERRED_HEAVY_RUNS.md`):** (a) the **4 giant families** by payload
  (PERMEABLE_PORT_HERMITIAN 62 MB, UNIFORM_LIMIT_RESIDUAL 37 MB, UNIFORM_LIMIT_S11CC1_OPERAND 33 MB,
  PERMEABLE_DISSIPATION_VS_OMEGA_TAU 12 MB) — UNMEASURED, the only place the R1 `cS0` fold activates a joined
  residual; (b) the **full per-family symbolic residual** as belt-and-suspenders for the inheritance-argument
  families (REP_INVARIANCE / CONTROL_* / DEGENERATE / DIMENSIONS — each plausibly an identical symbolic transform of
  the proven kernel, ⚠ argued not per-family-computed). The UNDECIDED items above (DTN_OPERATOR, ENERGY, flat-leg,
  seal 5) also close there.

## Method notes
- ⭐ **A nonzero cross-engine residual is not a disagreement** — the comparator prints raw (rule 2), so the reconcile
  is the staged representational bridge: apply the JUSTIFIED identities and watch what remains. A residual collapsing
  to the object's DEFINING relation (here the dispersion `q²=ω²/c_s0²−|k|²`) is the strong confirmation, not a
  coincidence. [[feedback_reconcile_representational_bridge]] [[feedback_matching_number_is_not_evidence]]
- ⛔ **Rule 17 (bg density) — kept a SURFACED freeze, ⛔ NOT waved through as AGREE.** WL keeps `ρ4,constant` as a
  LIVE applied field bound to `(ρbr/W_0)·W_bg(x)`; PY froze it to a bare constant. Neither engine differentiates ρ
  (0 derivatives) and both use it as a local `μ_s=μ_θ/ρ`, so the *functional dependence* agrees — but a bare scalar
  is not globally a live field, and folding one to the other to declare AGREE is precisely the rule-17 collapse this
  rebuild exists to catch (a step-record review leg caught the earlier over-claim). ⇒ seal 5 stays UNDECIDED; the
  c2 re-adjudication is **MANDATORY** (c2's fold sums over the face, where ρ(x)'s variation is load-bearing).
  [[feedback_never_freeze_a_varying_field]] [[feedback_basis_independence_must_not_freeze_spurion]]
  [[feedback_handcode_comparison_never_blanket_collapse]]
- ⛔ **Never blanket-collapse the cross-engine bridge.** Each identification was read from the two engine SOURCES and
  justified independently of the outcome (FT-of-derivative jet identity; σ_W=η_bg·W_0/L_W binding; on-shell
  dispersion; ε-placement) — ⛔ never tuned to force a zero. [[feedback_handcode_comparison_never_blanket_collapse]]
- ⛔ **Serialize CAS jobs; watch RSS.** c1 measured LIGHT (comparator peak ~317 MB on 30 GB); the 4 giants + full
  residual are the ≥64 GB work. Detached launch (harness reaps `run_in_background` ~87 s).

## What's next
**c2** (self-energy fold): fold the closed response into `S11CB_SLAB_OPERATOR`, re-extract the closed off-diagonal
kernel from the CLOSED full operator → the coupled nonlocal self-energy operator. Held folds: extract-then-close
ordering (close FIRST then re-extract — extract/eliminate don't commute), θ-row `Λ_X`(traction)/`J_s`(mass)
routing, substitution-increment emit (c2's self-energy = closed − open-symbolic, so S11c-b's deferred cross-engine
residual cancels out of c2's residual). c2 = a CODE build (gpt-6-astra) with its own gated directive (2 decision
legs BEFORE the build). ⚠ NO per-substep card — N1 specifies ONE S11c roll-up card, produced only after S11c-e.
The ≥64 GB c1 cross-engine re-run (above) is the outstanding c1 work; when it lands it confirms the
inheritance-argument families symbolically and closes the UNDECIDED items (DTN_OPERATOR, ENERGY, flat-leg, and —
with a valid field-vs-field comparison — seal 5). ⛔ Until then c2 must treat `dtn_operator` and the density as
NOT cross-engine-closed.
