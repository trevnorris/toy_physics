# S11c-c2 N6 RECONCILE — build directive (a diagnostic SCRIPT; astra-authored)

**Role.** Adjudicate the per-engine (SymPy) N6 representation-invariance finding: uncorrupted
`REP_INVARIANCE_RESIDUAL` (`R_N6 = I_E − I_{M→E}`) is certified-nonzero in ~18 forward-block columns
(`THETA`/`E_W`/`TRANSVERSE_TO_THICKNESS`), both densities, `LAB_HELD`. ⛔ This is **NOT** pre-adjudicated: a raw
nonzero residual is a computed measurement, not a verdict. This build produces the instrument that tests whether that
residual **vanishes modulo the fixed material↔Eulerian defining relations** (⇒ invariance holds) or leaves an
unexplained retained-order discrepancy (⇒ a real finding). The instrument PRINTS operands and residuals; ⛔ it decides
nothing and asserts nothing — the disposition is adjudicated on our side.

**Author = Codex `gpt-6-astra` (high).** ⛔ The orchestrator does not author this instrument (E1/CAS-authorship).
This directive had its E1 question-vet (`_legs/S11c_c2_N6_reconcile_question_vet.md`, Codex-sol) and gets its 2
decision legs before you build. Your job ends at **build → verify the deliverable → report**. ⛔ Do NOT run any review
legs, ⛔ do NOT read `.claude/skills/`, ⛔ do NOT create scratch/analysis files beyond the deliverable + its `.out`.

## The corrected reconcile question (this is what the instrument answers — do NOT restate it as R = J)
For each fixed `(α,ρ)`: does the independently-constructed **native-material route** (full prescribed material scalar
pullback for μ, material face velocity, and the face builder's internal covector map) **reduce, under ONLY the
predeclared material↔Eulerian field relations + their derivative prolongations + the allowed weak/on-shell
normalizations, to the SAME retained `(η^{≤1},σ_W^{≤1})` operator as the Eulerian route?** Because the cleared
construction already maps BOTH operands into the common Eulerian basis, there is **no nonzero "representation offset"
to subtract away**: N6 requires `[I_E − I_{M→E}]_retained = 0` in the quotient by those defining relations. ⛔ Do NOT
construct a "justified channel J" and call `R_N6 = J` invariance — that is author-freedom leakage (it lets an author
absorb whatever `R_N6` contains). The evidence is that the residual **factors through the defining relations**,
established by testing the two intermediate bridges directly, and cross-checked end-to-end.

## Supplied (verified; unfalsifiable within this build) vs withheld
- **Supplied:** the N6 companion `scripts/S11c_c2_N6_diagnostic_sympy.py` and its imports (`S11c_a_*`, `S11c_b_*`,
  the c2 `Compiler`/`build_increment`/`face_factory`/`constitutive`/`pit` machinery); the cleared route-2
  construction `_measurements/S11c_c2_N6_route2_spec_astra.md`; §5c of `directives/S11c_c2_SHARED_PHYSICS.md`; the
  defining relations and formulas below. Reuse the diagnostic's machinery by import — ⛔ do NOT re-derive Z, the
  resolvent, or the DtN; ⛔ do NOT introduce a new background-density law.
- **Withheld:** there is **no supplied expected value** for any bridge residual or for `R_N6`. ⛔ No acceptance
  criterion referencing an expected value; ⛔ residual-zero is NEVER an exit/assert/`while`-until predicate. A
  certified-nonzero modular numerator is a one-sided certificate of discrepancy; an all-zero sample table means only
  "no nonzero found" at the stated conditional false-negative bound — ⛔ never "certified zero."

## The frozen defining bridge relations (predeclare BEFORE inspecting any residual)
Fixed anchoring `α`, density `ρ`; both operands already carry Eulerian covectors/rows/fields/tests/pressure
identities. The fundamental relation is the **fixed-anchoring material pullback** (route-2 spec §2, :39-70):
- `g_i = D_iρ₄/ρ₄`, `a_ρ = u_i g_i` (θ-advection); `h_α = u_i D_iW_bg/W_bg` at `LAB_HELD`, `0` at
  `MATERIAL_ADVECTED` (the `e_W` thickness shift). Field maps `θ ↦ θ + a_ρ`, `e_W ↦ e_W + h_α`, with their spatial
  **derivative prolongations** `D_i(θ+a_ρ)`, `D_i(e_W+h_α)`, the density Jacobian (quadratic-projection only), and the
  face-covector map `material_inverse_transpose`. Equivalently `Δρ − δρ_E − u·∇ρ⁰ = 0` at fixed anchoring (⛔ NEVER
  used to bridge `LAB_HELD ↔ MATERIAL_ADVECTED` — that is the §5c category error).
- Representative laws (route-2 spec :72-83): `RHO4_CONSTANT` ⇒ `ρ₄=ρ_br/W₀`, `g_i=0` (advection structurally absent);
  `RHOBR_CONSTANT` ⇒ `ρ₄=ρ_br/W_bg`, `g_i=−D_iW_bg/W_bg` (advection live).
- **Jet-vocabulary bridge (sanctioned; route-2 spec :128; diagnostic `N6_JET_BRIDGE` :351-356):** `a.grad_theta[i]`
  (`theta_d1`…) and `b.grad_theta[i]` (`grad_theta_1`…) are derivatives of the **same field** at the jet/physical-field
  boundary. Add this exact identity to the predeclared **and emitted** frozen-relation set; ⛔ forbid any broader
  symbol renaming (⛔ NOT an indiscriminate rename — the explicit table only).

These induce the two **discriminating c2-interface bridges** (the strong evidence — test each DIRECTLY, not only
through the end-to-end increment):
1. **Carrier bridge** `C_E = C_{M→E}` (route-2 spec :171). `C_{r,p}` are the pressure-slot carrier coefficients
   (`S_{P,r} = Σ_{p∈P} C_{r,p} p`). Because `m_coeff` is produced AFTER the material builder's covector map and uses
   the SAME pressure-slot identities `P`, `C_E − C_M` must reduce to zero under the fixed face-map identities. ⚠ Test
   the **coefficient residual directly**, row/face/slot/grade-wise — `Increment(C_E−C_M, S_E) = 0` alone is weaker
   (could arise from a kernel/source nullspace). ⚠ The diagnostic's existing `CARRIER_RECONSTRUCTION_RESIDUAL`
   (diagnostic:832-836) is imported-Eulerian vs Eulerian-**factory** — NOT this imported-Eulerian-vs-material bridge,
   which is separately needed.
2. **Combined-source bridge** `b_{E,s} = b_{M,s}` (route-2 spec §5, :152-178). The c1 source amplitude is
   `b_{r,s} = (1 + Λ_V/ρ_m)·V̄_{r,s} + (Λ_A/(ρ_m·ρ_br^bg))·μ̄_r`, `V̄_{r,s}=V_{r,s}/ε`. Test `b_{E,s} − b_{M,s}` per
   face/grade directly. (The opaque c1 response `p_{r,s}=𝓡 Z[b_{r,s}]` is coordinate-map-opaque — ⛔ never substitute
   θ / multiply a Jacobian into `DELTA_P`/`Z`/resolvent; the bridge lives at the source amplitude `b`, not inside Z.)

## The exact THREE-WAY affine split (⛔ the two-way split hides the cross-term; ⛔ `build_increment(m_coeff, es−ms)` is WRONG)
`build_increment(C,S)` is **affine, not linear, in S**: signature 0 always contributes the source-independent bare
term `−C·p` (diagnostic:491-506), i.e. `I(C,S) = −C·p + B(C,S)` with `B` the closed-response contraction (signatures
6/9/12, bilinear in `C` and `S`). With `ΔC=C_E−C_M`, `ΔS=S_E−S_M`, the exact split of `R_N6 = I(C_E,S_E) − I(C_M,S_M)`
telescopes to **three** channels that separately localize the discrepancy:
```
R_N6 = I(ΔC, S_M)  +  B(C_M, ΔS)  +  B(ΔC, ΔS)
       └ CARRIER ┘    └ SOURCE ┘     └ CROSS ┘
```
- **CARRIER channel** `= I(ΔC, S_M)` — carrier difference at the **material** source `S_M`; carries the bare `−ΔC·p`.
  Build with the diagnostic's `build_increment` on `(ΔC, ms)`. ⚠ Use `S_M` (not `S_E`): `I(ΔC,S_E)=I(ΔC,S_M)+B(ΔC,ΔS)`
  would fold the cross-term into the carrier channel and destroy localization.
- **SOURCE channel** `= B(C_M, ΔS)` and **CROSS channel** `= B(ΔC, ΔS)` — **closed-response-only contractions over
  signatures 6/9/12 ONLY** (no bare `−C·p` term; the bare term is source-independent and cancels at fixed carrier).
  ⛔⛔ Do NOT call `build_increment(·, ΔS)` — it re-adds a spurious `−C·p` signature-0 term. Build a
  signature-{6,9,12}-only bilinear contraction (reuse the `template`/kernel path, omit signature 0).
- **Three-way localization:** carrier-only failure ⇒ CARRIER moves; source-only ⇒ SOURCE moves; only the CROSS
  channel carries the joint `B(ΔC,ΔS)` term (so a carrier or source discrepancy is not misattributed).
- **SPLIT_CHECK (a required emitted guard, ⛔ not an assert):** emit `CARRIER + SOURCE + CROSS − R_N6` where `R_N6` is
  recomputed in-process (below). ⚠ The guard is satisfied by **exact samplewise-zero numerators on the shared PIT
  samples** (every draw: numerator ≡ 0 mod each prime) — ⛔ NOT by a structural zero node (`plus`/`minus` at :137,:163
  drop only literal-zero operands; they do not reduce `x−x`, so a faithful independent check stays an `add` node), and
  ⛔ NOT by emitting `number(0)` directly. PRINT the residual node AND its per-sample numerators; ⛔ do not assert.

## What the instrument computes + emits (per fixed α,ρ; PRINT operands then residual)
Build every stage **independently of `R_N6`** (⛔ never by subtracting from `R_N6`). ⚠ **Recompute `E`, `M`, `R_N6`
IN-PROCESS** (reuse the diagnostic's `build_increment`/`residual` on the same imports) and run **ONE `pit()` over the
whole reconcile object dict** so every object — bridges, channels, `SPLIT_CHECK`, `R_N6` — shares identical sample
points. ⛔ Do NOT join stored diagnostic PIT tables from the N6 `.out` (a cross-run join makes `SPLIT_CHECK` a vacuous
cross-run comparison; the diagnostic re-samples per case). Emit (tag prefix `S11CC2_N6RC_`), each as columns +
`numerator_denominator` + `nonzero_modular_numerator`:
- `CARRIER_EULERIAN`, `CARRIER_MATERIAL`, `CARRIER_BRIDGE_RESIDUAL` = `C_E − C_M` — the **coefficient-level** carrier
  residual, keyed row/face/slot/grade (the strong direct test; ⛔ not `Increment(ΔC,S)`).
- `SOURCE_EULERIAN`, `SOURCE_MATERIAL`, `SOURCE_BRIDGE_RESIDUAL` = `es − ms`. ⚠ `SOURCE_{EULERIAN,MATERIAL}` **ARE the
  diagnostic's `es`/`ms` circuits** (`source_terms` :378-389, :841-843): `es` = imported μ_E (`inputs.mu/ε`) + imported
  `inputs.geometry['face_velocity']`; `ms` = material `mu_m`(t=1) + material `m_v`. The `b_{r,s}` formula
  `(1+Λ_V/ρ_m)V̄_{r,s}+(Λ_A/(ρ_m ρ_br^bg))μ̄_r` (route-2 :157-159, `V̄=V/ε` already inside `source_terms`) is the
  **slot IDENTIFICATION** of what `source_terms` computes — ⛔ NOT a re-coded second expression (that would double-count
  ε or miss imported slot factors, and violate corollary 1). `ΔS := es − ms`.
- `CARRIER_CHANNEL` = `I(ΔC, ms)`, `SOURCE_CHANNEL` = `B(C_M, ΔS)`, `CROSS_CHANNEL` = `B(ΔC, ΔS)`, and
  `SPLIT_CHECK` = `CARRIER_CHANNEL + SOURCE_CHANNEL + CROSS_CHANNEL − R_N6` (samplewise-zero numerators, per above).
- Re-emit `R_N6` (recomputed in-process) so the reconcile `.out` is self-contained.
- Provenance: fingerprints of `C_E`,`C_M`,`es`,`ms` sources; PIT primes/draws/δ (the diagnostic's honest
  `family·max(per_prime)` bound); a census that each stage's inputs are the material vs Eulerian builders (⛔ not one
  derived from the other); the frozen defining-relation list (incl. the jet-vocabulary bridge) as an emitted record.

## Controls — able-to-fail, one-sided, SEPARATED by which stage each corruption can move, ⛔ never A−A
⚠ `a_ρ` lives ONLY in the constitutive material μ (`constitutive` :328-342); the open pressure-slot coefficients carry
neither θ nor μ (route-2 §4 :130-150), so `a_ρ` moves the **SOURCE**, ⛔ NOT the carrier. Use **separate** one-sided
corruptions (in a `/tmp` COPY only — ⛔ never here; this is the build legs' ablation; the shipped instrument emits the
uncorrupted stages only, structured so each bites):
- **Carrier bite:** corrupt the **material** covector/normal-map that actually feeds `C_M` (`material_inverse_transpose`
  / the material normal, route-2 :113), holding the material source `ms` at its **uncorrupted** value (⛔ do NOT rebind
  `V`/`m_v` from the corrupted map) ⇒ must MOVE `CARRIER_BRIDGE_RESIDUAL` (`C_E − C_M`) / `CARRIER_CHANNEL`, leaving
  `C_E`, the Eulerian operand, and `SOURCE_BRIDGE_RESIDUAL` unchanged. ⛔ NOT the Eulerian-factory tilt (diagnostic
  :808, `tilt_coeff`): tilt touches neither `C_E` (imported `e_coeff` :797) nor `C_M` (material `m_coeff` :809), so it
  cannot move this bridge — the Eulerian-factory tilt is the SEPARATE reconstruction/independence probe
  (`CARRIER_RECONSTRUCTION_*` :832-836), which must itself be **gated on the unmodified factory matching imported
  `C_E`** (a baseline nonzero is reconstruction drift, ⛔ not tilt evidence; route-2 :195,:209).
- **Source bite (`RHOBR_CONSTANT` only):** `a_ρ → 0` on ONE material constitutive route ⇒ must MOVE
  `SOURCE_BRIDGE_RESIDUAL` / `SOURCE_CHANNEL`, and leave `CARRIER_BRIDGE_RESIDUAL` and the **Eulerian operand
  byte-identical** (computed carrier independence). ⛔ Do NOT require `a_ρ` to move the carrier. **OBSERVE** the
  end-to-end `R_N6` response (route-2 :209) — ⛔ do not require it to move.
- **RHO4 computed absence:** for `RHO4_CONSTANT` (`g_i=0 ⇒ a_ρ=0`) the advection corruption is structurally A−A —
  emit the computed absence, ⛔ never an A−A residual.
- **Split identity:** `SPLIT_CHECK` samplewise-zero (per the affine-split section) regardless of the physics
  disposition — a structural guard on the split, independent of whether `R_N6` vanishes.

## ⭐⭐⭐ THE THREE SCRIPT CLAUSES — non-negotiable (verbatim)
> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit` payload is a CAS object (an
> expression, a PIT numerator table, a symbolic-difference node) — ⛔ never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** Compute → emit → (only then, if at all) guard. ⛔ No `assert residual==0`,
> ⛔ no residual-zero exit/`while` predicate. A certified-nonzero numerator is the one-sided certificate.
> **3. Interpretation belongs to the reconcile record.** ⛔ The script does not editorialise or label a disposition.

## ⭐⭐ FOUR COROLLARIES (verbatim)
1. ⛔ A hand-typed CAS object is still hand-typed — the carrier/source/channel objects must be **reached by
   computation** from the material and Eulerian builders. The ONLY place physical symbols are combined by hand is in
   constructing the action/ansatz; every control re-enters the chain at the builder, ⛔ never at a result.
2. ⛔ The tag NAME names the object, ⛔ never its value/sign/shape (`*_BRIDGE_RESIDUAL` names the object, carries no
   claim about whether it is zero).
3. ⛔ No tautological residual — the two operands of every emitted residual come from **independent routes** (material
   builders vs imported Eulerian); verify independence by one-sided corruption (build-leg ablation). ⛔ Never difference
   an object against its own substitution.
4. ⛔ Emission is conditional only on which package/quantity a tag belongs to, ⛔ never on a payload's value; identical
   payloads across the two routes are the INVARIANCE finding and must both appear.

## Deliverable + verification
- Deliverable: `scripts/S11c_c2_N6_reconcile_sympy.py` (a companion; imports the N6 diagnostic's machinery). Run it,
  writing `.out` to an ABSOLUTE path OUTSIDE the repo, per case with `timeout` guards; the finite-field PIT is
  per-case compute-heavy — run ONE case at a time, ⛔ never a full-symbolic zero-test over all grades.
- Verify: the deliverable exists, is non-empty, imports cleanly, and its `.out` shows the emitted tags with PIT tables
  and `SPLIT_CHECK` present as an independent residual node whose **shared-PIT numerators are samplewise 0** (print the
  node + its per-sample numerators; ⛔ NOT `number(0)`, ⛔ NOT a structural node-identity). Report token usage. ⛔ Do NOT
  open/interpret the residual disposition — that is the orchestrator's, after the build legs.
