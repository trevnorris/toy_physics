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
> **δp_s (pressure) and J_s (relative-flux) response coefficients AGREE** (all leaves collapse to exact zero at
> physical kinematics — sweep `AGREE=54`; `_measurements/S11c_c1_reconcile_coeff_sweep.py`). **No cross-engine
> DISAGREEMENT was found among anything measured.** ⛔ **What is NOT earned — a closure-wide AGREE** (two review legs
> narrowed my claims, rule 15): seal 5 (background density) stays a **surfaced rule-17 freeze (UNDECIDED)**; the **t_s
> (traction) response leaf** (WL zero-padded 4-vector vs PY scalar), the **raw DTN_OPERATOR** form, the **ENERGY**
> audit, and the off-diagonal **flat-resolvent leg-labeling** are **UNDECIDED**; the **4 giant
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
  an agreement is independent construction, not a copy (⚠ SCOPED — see *Method notes › Independence is SCOPED*: the
  spec supplied composition + expected values for some objects, so agreement there is partly fidelity-to-supplied
  structure; the two-momentum kernel is the independently-confirmed object).
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
  real-assumption — AGREE** (unified, contributes nothing); **(3) μ_θ — AGREE for δp_s/J_s** (opaque/supplied both
  sides; the pressure + relative-flux coefficients `ε·A≡B` to machine zero at physical kinematics — sweep `AGREE=54`,
  `_measurements/S11c_c1_reconcile_coeff_sweep.py`; ε sits in PY's μ_θ `.py:840` vs WL's coefficient, `closed_coefficients`
  `.py:723-755` omits ε; ⚠ the **t_s traction leaf is UNRESOLVED** — WL zero-padded 4-vector vs PY scalar, below);
  **(4) regime/parity naming — AGREE on content** (`OUTPUT_X`/`INPUT_X` vs bare is a bridgeable keying
  convention; PARITY_MATRIX joins and IS the proven kernel); **(5) background density (rule 17) — UNDECIDED,
  SURFACED FREEZE** ⛔ NOT AGREE (the functional dependence agrees and neither engine differentiates ρ — 0
  derivatives — but a PY bare constant is not globally a WL live field `rhoBrBgRho4Constant(x)=(ρbr/W₀)WBg(x)`;
  folding one to the other is the rule-17 collapse this rebuild exists to catch — kept surfaced; c2 re-adjudication
  MANDATORY). **Also surfaced UNDECIDED:** the raw DTN_OPERATOR face-form; the off-diagonal flat-resolvent
  leg-labeling (PY output-leg `q_out` vs WL input-leg for the flat half-space symbol — diagonal-exact, differs
  off-diagonal in the MATERIAL response coefficients); the ENERGY audit (PY closed-form vs WL far-field integral).

## Established (cross-engine AGREE) vs owed (surfaced/deferred)
- **ESTABLISHED — per-engine SOUND (2-leg each) AND cross-engine AGREE:** the **two-momentum DtN KERNEL** (all 4
  cases, exact zero, load-bearing), the **δp_s (pressure) + J_s (relative-flux) response COEFFICIENTS** (all leaves
  collapse at physical kinematics, sweep `AGREE=54`), the flat symbol, the parity matrix (= the kernel), the
  degenerate loci, the dimensions. ⛔ **NOT in this list — UNDECIDED (§ below):** the raw DtN OPERATOR (whole-object
  noncommutative form — kernel-AGREE does NOT extend to it), and the **t_s (traction) response leaf** (WL zero-padded
  4-vector `(0,0,0,scalar)` vs PY scalar; the flat-traction cases do not cleanly reduce). The Hermitian/reactive
  dissipation parts AGREE only at c1's first order (WL's true-area-weighted adjoint `BoundaryMeasure` → PY's plain
  conjugate-swap at first order; O(η²)→S11c-e; full symbolic DEFERRED). **For c2's consume-set** (`dtn_operator`,
  `dtn_flat_symbol`, `dtn_kernel`, face response — per-engine-verified 44-row delta over `S11c_b_exports.py`): the
  **kernel + flat symbol + δp_s/J_s response are cross-engine AGREE**; the **`dtn_operator` whole-form, the t_s
  traction, and the seal-5 density are UNDECIDED — c2 must NOT treat them as cross-engine-closed** (re-adjudicate in
  c2, see What's next).
- **OWED / SURFACED (UNDECIDED, on this box) — ⛔ not AGREE, ⛔ no predicted residual:** seal 5 (background density
  field-vs-constant freeze); the **t_s (traction) response leaf** (WL zero-padded 4-vector vs PY scalar — a
  scalar-vs-vector representation); the raw **DTN_OPERATOR** face-form (PY per-face, WL per-branch); the **ENERGY**
  audit (PY closed-form propagating energy vs WL unevaluated far-field flux integral); the off-diagonal **flat-resolvent
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
  rebuild exists to catch (a step-record review leg caught the earlier over-claim). ⚠ **The channel is O(εη), ⛔ NOT
  the "harmless O(η²)" an earlier reconcile draft claimed** (corrected by the retrospective spec review; both retro
  legs computed the derivative): since `μ_θ=O(ε)`, re-binding ρ to `background_density_map` gives
  `d(μ_s)/dη|₀ = −μ_θ·w₁/ρ_br` — **FIRST-order in the shape η**. ⛔ It is **not lost**: PY carries ρ **opaquely
  (0 derivatives, no `w1_profile` baked in — verified against the real export by both retro legs + a Codex-sol verify
  pass, residual 0)**, so c2 **recovers** the exact O(εη) channel by that re-bind. ⇒ seal 5 stays UNDECIDED (a
  field-vs-constant representational question, ⛔ NOT a dropped physics channel and ⛔ NOT harmless-because-higher-order —
  c1's engines/exports STAND, no reopen); the c2 re-adjudication is **MANDATORY** (c2's fold sums over the face, where
  ρ(x)'s variation is load-bearing; c2 v2 §3d.1 mandates the re-bind).
  [[feedback_never_freeze_a_varying_field]] [[feedback_basis_independence_must_not_freeze_spurion]]
  [[feedback_handcode_comparison_never_blanket_collapse]]
- ⚠ **Exact grazing — claim SCOPED: `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`** (retrospective spec review). The
  **DtN inverse `N⁻¹`** (hence `Z=iρ_mω·N⁻¹`) is **nonanalytic as `q_out→0`** (the grazing limit; at exact double
  grazing `N=η·B` so `N⁻¹` — and thus `Z` — carries a `~1/η` Laurent pole, ⛔ not a Taylor series in η), so the
  first-shape coefficient `Z₁` is a valid **non-grazing asymptotic** coefficient ONLY — `‖N₀⁻¹N₁‖≪1` imposed on both
  legs (both momentum legs bounded away from grazing). ⚠ It is `N⁻¹`/`Z` that is singular, ⛔ **NOT the permeable
  face-closure resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹`**, which for generic `Λ_A≠0` is **regular** at grazing (`O(η)`, and
  `[I+(Λ_A/ρ_m²)Z]⁻¹·Z → ρ_m²/Λ_A` finite — both retro legs printed "D⁻¹ 1/η-pole = False, Z=C·N⁻¹ 1/η-pole =
  True"). The GRAZING regime label is carried as a keying axis (reconcile §4 `DTN_GRAZING_BEHAVIOUR`) and the
  strict-`v_dr=0` rest-frame qualification is correct, but the exact-grazing *threshold response* is ⛔ **not
  claimed** at c1's first shape order — it needs a degenerate threshold solve (deferred).
  [[feedback_never_freeze_a_varying_field]]
- ⚠ **Independence is SCOPED, ⛔ not blanket** (retrospective spec review). The c1 spec supplied the composition
  recipe and some expected structural values — rigid-shift cancellation, the flat `Z₀=ρ_m ω/q_out`, the zero-jet
  outcome (a rule-5/rule-3 leak). So for THOSE objects part of the cross-engine "agreement" is **fidelity to the
  supplied structure**, ⛔ not independent discovery. What IS independently confirmed is the **two-momentum DtN
  kernel**: both blind engines AND both retrospective spec-review legs re-derived it from first principles (Grok
  Route A/B, Codex boundary solve — residual 0). Read the "agreement is independent construction" claim (in *The two
  blind engines* above) with that scope: the kernel earns it; the supplied objects earn it only as fidelity.
- ⛔ **Never blanket-collapse the cross-engine bridge.** Each identification was read from the two engine SOURCES and
  justified independently of the outcome (FT-of-derivative jet identity; σ_W=η_bg·W_0/L_W binding; on-shell
  dispersion; ε-placement) — ⛔ never tuned to force a zero. [[feedback_handcode_comparison_never_blanket_collapse]]
- ⛔ **Serialize CAS jobs; watch RSS.** c1 measured LIGHT (comparator peak ~317 MB on 30 GB); the 4 giants + full
  residual are the ≥64 GB work. Detached launch (harness reaps `run_in_background` ~87 s).

## Carry-forward corrections (lower-severity — from the retrospective spec review, 2026-09-05)
These change what may be CLAIMED, ⛔ not any computed object (c1's engines/exports STAND — no reopen). Full list +
both leg reports: `_measurements/S11c_c1_spec_retro_review_adjudication.md`,
`_measurements/S11c_c1_spec_retro_review_{grok,codex_sol}.txt`.
- **Energy-residual orientation.** The independent energy identity is `P_face + P_∞ = 0` (traction work + **positive
  outgoing** far-field Poynting = 0, with `P_∞` the outgoing flux), ⛔ NOT the literal `A−B` the spec wrote — a
  literal `A−B` on the written operands is `−2·δp·V*` for the supplied `t_s` and vanishes only for the WRONG `t`
  sign. PY correctly compared to `−outgoing_flux` (`.py:1407,1459`), so no recomputation is owed; the fix is the
  spec/tag SEMANTICS (to write it as a **vanishing `A−B` residual**, define the bulk subtraction operand `B` as
  **minus** outgoing Poynting — `A−B = P_face−(−P_∞) = 0` — or emit `A+B`). Carry to c2's energy control. (both legs)
- **`h_s` graph-vs-outward + `N`/`Z` terminology.** `h_s=s·W_bg/2+ζ_s` is the signed lab GRAPH height (face-odd); the
  shape displacement that enters the kernel is the face-EVEN outward displacement `a_s=(W_bg−W₀)/2`. Both engines used
  the correct outward height (PY `shape_source:592`, WL `:177`), so the kernel AGREES — notation only. Separately: `Z`
  is the impedance (Neumann-to-Dirichlet response), the mathematical DtN is `N`, related by `Z=iρ_mω·N⁻¹` (⛔ not
  synonyms); `N₁=−(N₀M_aN₀+Div(a∇)+κ²M_a)` vs `Z₁=iρ_mω·(−N₀⁻¹N₁N₀⁻¹)`. Terminology fix; no recomputation. (both legs)
- **Density as a multiplication operator.** Name the live `1/ρ_br,bg` as a multiplication operator `β=M_{1/ρ_br,bg}`
  so the O(εη) channel (Method notes › rule 17) cannot be emitted as a bare constant. ⇒ **the c2 build directive must
  require this naming** (it makes the freeze un-emittable). Carry to c2. (both legs)
- **`K_a` is Hermitian, ⛔ not anti-Hermitian.** The reactive object `K_a=(Z−Z†ₐ)/(2i)` is **Hermitian**
  (`K_a†−K_a=0`); the anti-Hermitian part is `(Z−Z†ₐ)/2`. The spec's "anti-Hermitian" label for `K_a` is wrong — fix
  the label wherever the dissipation objects are reused (c2 / S11c-e). (Codex NIT)
- **Evanescent caveat covers all second-shape grades.** The evanescent-nullspace completion is **second-shape order
  = η², η·σ_W, σ_W²** (η and σ_W are INDEPENDENT grades, spec §2c), ⛔ not only `O(η²)`. Physical conclusion (no
  first-order passivity violation) unchanged. (Codex NIT; Grok concurs the caveat is physically right)
- **Drain-projection `O(σ_W²)` wording.** The `O(σ_W²)` applies to the geometric **drain-tilt** projection
  (`n̂·ŵ−s=O(σ_W²)`), ⛔ NOT to convection — convection is dropped as an inheritance of N11a's standing rest-frame
  limit (`d(rel)/dη=0`, has no η), ⛔ not because "convection is O(σ_W²)". Wording fix. (Grok NIT)
- **Flat `Z₀` / rigid-shift expected-value leakage.** `§5d` typed `Z₀=ρ_m ω/q_out` and `§3a` typed the rigid-shift
  "must cancel" residual — rule-5 leaks in the (builder-facing) c1 spec. Physics correct; folds into the *Independence
  is SCOPED* note above (agreement on `Z₀`/rigid-shift/zero-jet is partly fidelity-to-supplied-values). Spec hygiene;
  the kernel remains independently confirmed. (both legs)

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
