# S11c-c2 — SHARED PHYSICS (the self-energy fold: closing the slab operator with the curved-bulk response)

**S11c-c2** is the second half of the S11c-c curved-interface bulk closure (the S11c-c decision-list row was split
c1/c2 by user choice 2026-09-03; `directives/S11c_decisions.md` N2). c1 solved the curved two-face outgoing bulk
problem and **exported the closed permeable face response** `(δp_s, J_s, t_s)(V_s, μ_θ)`, the nonlocal DtN/impedance
operator, its flat symbol, and its two-momentum kernel (`steps/S11c_c1_curved_bulk_closure.md`; exports
`scripts/S11c_c1_exports.py`). **c2 folds that closed response into the S11c-b variable-coefficient slab operator
`S11CB_SLAB_OPERATOR`** — whose θ-row and mechanical rows still carry the face pressure `δp_s` **symbolically** — and
**re-extracts the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the CLOSED full operator**, yielding the
coupled **nonlocal self-energy operator**: the transverse sector now carries a self-energy threaded through the bulk
DtN. This document is the physics authority for the two blind c2 engines and their comparator. Tag prefix `S11CC2_`.

The SymPy engine reads the inherited model through `ledger_fold.load_model` over the atomic frozen base
`scripts/S11c_b_exports.py` **with the c1 delta `scripts/S11c_c1_exports.py`** folded on top (§7), binding only its
declared `IMPORT_KEYS`; the Wolfram engine imports nothing and re-derives every consumed object from the sibling
specs (`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control). Blindness is the control: an
agreement is independent construction, not a copy.

⭐ This is a **CODE build** authored against c1's and S11c-b's reviewed exports. Rule 7: it gets **two decision legs
(Codex `gpt-5.6-sol` xhigh + Grok) before any builder**; the build directive that follows gets its own two decision
legs before the build. ⚠ **This is spec v2** — folded once from the v1 two-leg decision gate
(`directives/_measurements/S11c_c2_spec_review.md`), which found the v1 isolation claim false, the fold operation
under-specified against the real θ-row, and the re-adjudication set incomplete.

---

## 0 · Scope

**In scope.** (1) **Close** the S11c-b slab operator by substituting the c1 closed face response into its symbolic
`δp_s`-slots, per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, summing the two faces `s∈{+,−}`
that the slab EOM already sums (§1c). (2) **Re-extract** the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the
**closed** full operator by the S11c-b §3c weak variational restriction — the **close-then-extract** ordering (§2).
(3) Emit the resulting **nonlocal self-energy** as the **substitution increment** (§3c), both operands re-extracted
from `SLAB_OPERATOR` with the same extract. (4) **Re-adjudicate the c1 items the fold makes load-bearing** (§3d):
the background-density field-vs-constant freeze (c1 seal 5, rule 17), the `t_s` traction representation, the DtN
whole-form-vs-kernel, the **traction-vs-slab mechanical-power pairing** c1 assigned to c2, the flat-resolvent
leg-labeling, and the `μ_R,bg` form control c1 reserved for c2.

**Out of scope (named, not solved).** The profile-conditioned spectrum/scattering and leakage rates (S11c-d, `N5`);
the leakage observable and falsification (S11c-e, `N7`); the nonlinear-light program (`N10`); the O(η²) evanescent-
leakage sign on the flat nullspace (S11c-e, c1 §3b caveat). c2 introduces **no** convective bulk operator (`N11a`).
The full per-family cross-engine residual and the four c1 giant families remain **deferred to a ≥64 GB box**
(`DEFERRED_HEAVY_RUNS.md`); c2 must be constructible and cross-engine-testable on this box for its own increment
(§7), and name — not silently absorb — anything it cannot close here.

---

## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build

Everything in §1 is an input. The substitutions, orderings, extractions, self-energy structure, and the six
re-adjudication verdicts of §§3–5 are **outputs**; ⛔ none is stated here.

### 1a · Inheritance and the two imported models

The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12`), and admissibility are exactly S11c-a
§§1–2 / S11c-b §§1–2, inherited by pointer. c2 consumes two already-built, per-engine-reviewed models:

- **S11c-b** (`scripts/S11c_b_exports.py`, committed `af560257`; step record `steps/S11c_b_variable_coefficient_operator.md`):
  the variable-coefficient slab operator `slab_operator` over `{u,θ,e_W}` (write-key `slab_operator`, tag
  `S11CB_SLAB_OPERATOR`), its term provenance `slab_operator_term_origins`, the off-diagonal coupling kernel
  `coupling_kernel`, the constitutive operand `mu_theta_operator`, and — inherited from S11c-a — the θ-row flux
  closure carrier `closure_shape_deriv`. ⚠ **S11c-b's full cross-engine residual is DEFERRED (≥64 GB); its coupling
  content is per-engine leg-verified only, and two whole-row sign conventions + #90's two flags are
  cross-engine-UNVALIDATED** (`steps/S11c_b_variable_coefficient_operator.md:112-115`): the **kinetic** convention
  `−K` PY vs `+K` WL (a bulk term, independent of the response slots), the **face generalized-force** convention PY
  `+diff` vs WL `−linearVirtualVariation`, and the **#90 closure-fold sign**. ⛔ The face-force and closure-fold
  conventions are the coefficients of the slots c2 substitutes into (§3c) — they do **not** cancel from c2's residual.
- **S11c-c1** (`scripts/S11c_c1_exports.py`, per-engine-reviewed 44-row delta; step record `steps/S11c_c1_curved_bulk_closure.md`):
  the closed face response, write-keys **`s11c_c1_face_response`** / **`s11c_c1_face_response_coeffs`** (⛔ NOT the
  bare S11b `face_response`/`face_response_coeffs`, which c1 imported as its uniform-limit regression operand);
  the DtN operator `dtn_operator` (per anchoring/face: `s11cc1_dtn_operator_{lab_held,material_advected}_{plus,minus}`),
  the flat symbol `dtn_flat_symbol`, the two-momentum kernel `dtn_kernel`; the per-case response resolvents
  `s11cc1_response_resolvent_{α}_{s}_{ρ}_constant`; the momentum symbols `s11cc1_{k,q}_{in,out}put*`; the FT-of-profile
  objects `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_{1,2,3}`.

### 1b · What is cross-engine AGREE vs UNDECIDED in the c1 import — SUPPLIED HONESTLY (rule 6/16)

⭐ **AGREE (cross-engine, load-bearing; c1 reconcile `directives/_measurements/S11c_c1_comparator_reconcile.md`):**
the **two-momentum DtN kernel** `dtn_kernel` (all four anchoring×face cases collapse to exact zero off-diagonal; a
wrong jet sign or a one-leg freeze leaves it nonzero) and the **`δp_s` (pressure) + `J_s` (relative-flux) response
coefficients** (all leaves collapse at physical kinematics; sweep `AGREE=54`). The flat symbol on the diagonal and
the parity matrix (= the kernel) AGREE.

⛔ **UNDECIDED — c2 must NOT treat these as cross-engine-closed (c1 step record §"Established vs owed"):** (i) the
**background density** `rho_br_bg_rho4_constant` — a **surfaced rule-17 freeze**: WL keeps it a live applied field,
PY froze it to a bare constant; (ii) the **`t_s` traction** — WL a zero-padded 4-vector `(0,0,0,scalar)`, PY a
scalar; (iii) the **raw `dtn_operator` whole-form** — kernel-AGREE does **not** extend to the whole noncommutative
object; (iv) the **off-diagonal flat-resolvent leg-labeling** — PY output-leg `q_out` vs WL input-leg, equal on the
`k=k′` diagonal, differing off-diagonal in the MATERIAL_ADVECTED response coefficients; (v) the **ENERGY** audit (PY
closed-form vs WL far-field integral). Each of (i)–(v) enters c2's fold load-bearingly (§3d); ⛔ folding any of them
to force AGREE is the exact defect this rebuild exists to catch.

### 1c · The S11c-b slab operator and the symbolic slots the fold closes — SUPPLIED (grounded in the real row)

`slab_operator` is the divergence-form variable-coefficient operator over `{u,θ,e_W}`, per `(anchoring, density)`,
with every background coefficient (`μ_R,bg`, `W_bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, the `Σ_E⁰` map) and its first jet
retained live. **It is a two-face object** — the mass row is the sourced evolution `∂_tΣ + ∇_x·(Σv) = −(J₊+J₋)`
(both faces summed) and the mechanical rows sum both per-face tractions. The rows are: `U_BODY_BALANCE`,
`E_W_BALANCE` (the constraint-reduced mechanical rows, pin B), and the θ/mass row `evolution_mass_balance − Σ_s
closure_shape_deriv_s` (the #90 fold). ⛔ **The explicit relative-flux `J_s` carrier has already been eliminated
from the θ-row by the #90 subtraction** — what remains is the flux **closure** written out as
`Λ_A𝒜_s + Λ_V V_s` with the face pressure `δp_s`/`d_w_delta_p_s` **symbolic**. In the real row (verified against
`closure_shape_deriv`), the θ-row carries `−4I·Λ_A(−δp_s/ρ_m + μ_θ/ρ_br)/(ωτ_A+I)`, `−2I·Λ_V W_0 e_{W,t}/(ωτ_V+I)`,
and the pressure-jet terms in `d_w_delta_p_s`. ⇒ **the fold operation is "substitute the closed `δp_s(V_s,μ_θ)` and
its w-jets `d_w_delta_p_s` into the symbolic `delta_p_±`/`d_w_delta_p_±` slots"**, ⛔ **NOT "substitute a closed
`J_s`" (there is no `J_s` slot; adding one would double-count the already-folded flux closure).**

`mu_theta_operator` (`S11CB_MU_THETA_OPERATOR`) is the separate held-fixed constitutive operand (neither row).
⛔ The face-force and #90 closure-fold sign conventions above are **not** normalized by c2; the §3c increment carries
each engine's own convention, and the §3d.4 mechanical-power pairing adjudicates the face-force sign (rule 1/6).

### 1d · The face closure laws and the Λ-channel routing — SUPPLIED (this routing is c2's task)

Per anchoring `α` and face `s`, on the curved-face objects of the S11c-a substrate (⛔ not flat Cartesian objects):

```text
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br,bg⁰ ,   μ_θ = S11CB_MU_THETA_OPERATOR (held-fixed operand) ,
J_s = Λ_A(ω) 𝒜_s + Λ_V(ω) V_s ,               ⇐ the FLUX closure carries only Λ_A, Λ_V ,
t_s = −(δp_s + Λ_X(ω) 𝒜_s) n̂_s ,              ⇐ the TRACTION carries Λ_X (a reciprocal-traction channel) ,
n̂_s·v_bulk,s = V_s + J_s/ρ_m ,                interfacial mass balance (kinematics, not a result) ,
Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,  the three τ_I independent .
```

⚠ **`Λ_X` appears ONLY in the traction `t_s`** — ⛔ **not** in the flux closure `J_s`, and ⛔ **not** in the S11c-a
T-i shape derivative `closure_shape_deriv` (which is the shape derivative of `J_s − Λ_A𝒜_s − Λ_V V_s = 0` only; c1
§1d:157-160). The mechanical routing of `t_s` (which carries `Λ_X` and `δp_s`) is c2's, via the S11c-a traction
`traction` (T-d) / virtual work. c1 supplied the closed `(δp_s, J_s, t_s)(V_s, μ_θ)` by solving these with the bulk
relation `δp_s = Z·v_bulk,s` (the DtN operator), keeping `Λ_I(ω)`/`μ_θ` symbolic. **The routing is c2's** (c1
§1d:160): the closed `δp_s` (+jets) closes the θ-row's already-folded flux terms and enters `t_s`; the closed `t_s`
enters the mechanical rows. The equations of motion are obtained by S11b's balance-law/virtual-work method, ⛔ not by
putting an irreversible response kernel in an ordinary action.

---

## 2 · The three fold objects and the non-commuting ordering — SUPPLIED framing

`coupling_kernel` was obtained by **weak variational restriction of the OPEN (δp_s-symbolic) §3b operator** (S11c-b
§3c). Closing the face response threads the nonlocal bulk DtN operator `Z` through the operator, including into the
transverse diagonal block, so:

```text
extract( close( SLAB ) )  ≠  close( extract( SLAB ) )      (extract and close do NOT commute) .
```

⭐ **The correct object is `extract(close(·))` — close the FULL operator first, then re-extract.** The counterexample
fixing the ordering: for `R_x = x + p` with the closure `p = α y`, extracting the `x`-block first (`∂R_x/∂x = 1`, `p`
dropped) then closing gives a zero `x→y` coupling, whereas closing first (`R_x = x + α y`) then extracting gives `α`;
the closure introduces coupling the extract-first route discards, and that induced coupling is the self-energy.

⚠ **Three DISTINCT objects — each gets its own name (⛔ do not call all three "self-energy"):**

```text
S11CC2_CLOSED_COUPLING_KERNEL  = extract(close(SLAB))                  (the re-extracted closed off-diagonal block) ,
S11CC2_ORDERING_COMMUTATOR     = extract(close(SLAB)) − close(extract(SLAB))   (the §5a form-ablation control) ,
S11CC2_SELF_ENERGY_INCREMENT   = extract(close(SLAB)) − extract(SLAB)  (c2's OWN-ROWS object, §3c) .
```

⛔ `close(extract(SLAB))` (close the already-extracted OPEN kernel) is the §5a ablation, ⛔ not a construction route.

---

## 3 · The self-energy construction (OUTPUTS)

Every object below is computed for both anchorings `α` and both density representatives `ρ`; it carries its computed
`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component value, sign, order, parity, or grade
in this document.

### 3a · Close the full operator

Substitute the c1 closed face pressure `δp_s(V_s,μ_θ)` and its w-jets `d_w_delta_p_s` (from `s11c_c1_face_response`)
into the symbolic `delta_p_±`/`d_w_delta_p_±` slots of `slab_operator` — the θ-row's already-folded flux terms
(`closure_shape_deriv`) and the mechanical-row traction (`traction`, carrying `Λ_X`). ⛔ Do **not** add a closed `J_s`
term (the θ-row's flux is already `Λ_A𝒜_s + Λ_V V_s`; the fold closes its symbolic `δp_s`, §1c). Because
`δp_s = Z·v_bulk,s` with `Z` a nonlocal two-momentum operator, the closed rows carry `Z` — the response elimination
is the c1 operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹`, ⛔ not a scalar division. Map c1's `V_s`, opaque `μ_θ`, and
per-case resolvent symbols onto the slab-row fields explicitly (the build directive freezes the exact symbol map).
Retain the curved-face measure/normal and the live background coefficients; ⛔ do not freeze any coefficient at its
constant binding before the fold (rule 17, §3d.1). **Emit the closed operator as the assembled two-face object per
`(α,ρ)`** (both faces summed, as the slab EOM is), the per-face contributions as provenance, and the face-parity
blocks.

```text
⇒ S11CC2_CLOSED_SLAB_OPERATOR (assembled, per (α,ρ)) , S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS (per face) ,
  S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS .
```

### 3b · Re-extract the off-diagonal coupling from the CLOSED operator

Apply the S11c-b §3c **weak variational restriction verbatim** (independent divergence-free transverse and
curl-free longitudinal trial/test displacements, independent `θ`/`e_W` trial/test fields, compact support in the
in-plane interior so the IBP boundary term is fixed to zero; the inherited face conditions still apply) to the
**closed** full operator of §3a, extracting both off-diagonal blocks → `S11CC2_CLOSED_COUPLING_KERNEL`. ⛔ Do not
implement the split by zeroing only undifferentiated field occurrences (S11c-b §3c). Emit both blocks and, **only if
the pairing-based adjointness residual is a genuine independent route**, that residual; ⛔ it is not the mixed second
derivative of a scalar energy (rule 2 corollary 3) — if the two blocks are adjoint by construction, emit them and say
there is no independent second route rather than dressing a structural zero as a check.

```text
⇒ S11CC2_CLOSED_COUPLING_KERNEL , S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS ,
  S11CC2_SELF_ENERGY_ADJOINTNESS_RESIDUAL  (emit only if an independent route exists) .
```

### 3c · The self-energy as the substitution increment — the export representation

Emit the self-energy as the **substitution increment**, both operands re-extracted from `slab_operator` **with the
same weak restriction** (⛔ NOT the imported open `coupling_kernel`, which the two-leg-gated
`directives/export_ledger_bind_closure_design.md:148-153` bars as a construction operand — c2 re-extracts):

```text
S11CC2_SELF_ENERGY_INCREMENT  =  extract( close(SLAB) )  −  extract( SLAB )    (per α, ρ) .
```

Because the weak restriction is linear, this equals `extract(close(SLAB) − SLAB)` — the restriction of the pure
substitution increment, supported only where the closed `δp_s` replaced the symbolic slot. ⚠ **This is an export
representation, ⛔ not a residual check** (rule 2 corollary 3). ⭐ **What it does and does NOT isolate (corrected from
v1):** the increment **drops the S11c-b bulk/kinetic base** (every term with no `δp_s` slot — including the deferred
≥64 GB slab-operator content and the kinetic `−K/+K` convention — appears identically in `close(SLAB)` and `SLAB`
within each engine and cancels). It does **NOT** drop the **face-generalized-force and #90 closure-fold sign
conventions**: those are the coefficients multiplying the `δp_s` slots, so they multiply the increment and — being
cross-engine-UNVALIDATED — can leak into c2's cross-engine residual `(≈ 2·carrier·increment)`. ⇒ the comparator
**SURFACES** them (rule 1/6, §7) and the §3d.4 mechanical-power pairing **adjudicates the face-force sign**. Emit both
operands and the increment; ⛔ do not claim the increment isolates c2 from those two conventions.

### 3d · The six re-adjudications (rule 17 / c1 UNDECIDED carry-ins the fold makes load-bearing)

Each is a **computation whose cross-engine disposition is a finding**, ⛔ not a value to match:

1. **Background density, field-vs-field (rule 17, c1 seal 5).** Before the fold, bind `rho_br_bg_rho4_constant` to
   its live relation from `background_density_map` (`rho_br_bg_rho4_constant = W_bg·ρ_br/W_0`, `W_bg=W_0(1+η w₁)`),
   ⛔ not to a bare constant — else the fold silently repeats the c1 PY freeze while passing the manifest guards. The
   **two density representatives are the field-vs-constant pair**: `RHO4_CONSTANT` carries the live `ρ_br,bg⁰(x)`,
   `RHOBR_CONSTANT` the frozen `ρ_br`. Emit the self-energy for both and the field-vs-field object
   `S11CC2_DENSITY_LIVE_MINUS_FROZEN = SelfEnergy[ρ_br,bg⁰(x)] − SelfEnergy[ρ_br]`. ⛔ This is a live-field difference,
   ⛔ NOT a `∇ρ` jet (neither c1 engine differentiates ρ; at first shape order it is an O(η) multiplicative
   difference) and ⛔ NOT `ρ_m` (the bulk fluid density, a different object). The comparator joins the two engines'
   live-field objects **field-against-field**; ⛔ never a PY constant folded onto a WL field. ⛔ `∇ρ→0` is not an
   accepted corruption (`N6`).
2. **The `t_s` traction representation.** Carry `t_s` in its **native covector form** `t_s = −(δp_s + Λ_X𝒜_s) n̂_s`
   (a covector along the curved outward normal), ⛔ not a pre-collapsed scalar, so the mechanical-row contribution is
   the genuine covector pairing; the scalar-vs-4-vector difference c1 left UNDECIDED is carried into the comparator as
   a representation to reconcile, ⛔ not silently collapsed on one side. ⇒ `S11CC2_TRACTION_MECHANICAL_CONTRIB`.
3. **The DtN kernel vs whole-form.** Structure the fold so the load-bearing bulk object is the AGREE'd two-momentum
   **kernel** `dtn_kernel` (its `q_out(k),q_out(k′)` legs + on-shell dispersion). Emit — separately, so a leg can
   ablate it — **which self-energy terms depend on the raw whole-form `dtn_operator`** beyond the kernel
   (`S11CC2_DTN_WHOLEFORM_DEPENDENCE`). ⛔ Do not assume whole-form AGREE; ⛔ do not fold the whole-form into the
   kernel.
4. **The traction-vs-slab mechanical-power pairing (c1-assigned to c2).** c1's energy audit "uses no slab EOM row and
   does not import `S11CB_SLAB_OPERATOR` (that pairing is S11c-c2's, after the fold)" (c1 §3b:328-330); parent N2
   gives c2 the folded operator and S11c-e only the conversion FORM. Emit the post-fold mechanical-power / traction
   pairing against the slab kinetic-and-stored variation: `S11CC2_TRACTION_SLAB_POWER_PAIRING` and its residual. ⭐
   This is the control that **settles the load-bearing face-generalized-force sign** the increment (§3c) does not
   cancel — a one-sided `t_s`-sign corruption must move the mechanical-power residual against the slab kinetic term.
5. **The flat-resolvent leg-labeling.** c1 left the PY-output-leg vs WL-input-leg convention UNDECIDED, equal on the
   `k=k′` diagonal, differing off-diagonal in the MATERIAL_ADVECTED response coefficients. If the fold uses
   `dtn_flat_symbol` **only** as the uniform-limit regression operand (diagonal `k=k′`), pin it so and say so; if any
   MATERIAL off-diagonal self-energy term consumes it, it is a **sixth UNDECIDED item to re-adjudicate**, ⛔ not
   inherited as AGREE. Emit `S11CC2_FLAT_SYMBOL_USAGE` naming where and how it enters.
6. **The `μ_R,bg` form control (c1-reserved for c2).** c1 reserved the `μ_R,bg` form ablation for c2, "where `μ_θ` is
   composed with the slab variables" (c1 §5a). Since the fold composes the opaque `μ_θ` operand into the slab rows,
   emit the `μ_R,bg`-form ablation of the self-energy (§5) — the c1 reservation is discharged here.

---

## 4 · Objects to compute and emit

Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, multigraded and dimensioned:

- The **closed full slab operator** (assembled two-face, per-face provenance, parity blocks) — §3a.
- The **re-extracted closed off-diagonal coupling kernel** (both blocks), adjointness residual only if independent — §3b.
- The **self-energy increment** and its two same-extract operands — §3c.
- The **six re-adjudication objects** — §3d: `S11CC2_DENSITY_LIVE_MINUS_FROZEN`, `S11CC2_TRACTION_MECHANICAL_CONTRIB`,
  `S11CC2_DTN_WHOLEFORM_DEPENDENCE`, `S11CC2_TRACTION_SLAB_POWER_PAIRING`, `S11CC2_FLAT_SYMBOL_USAGE`, and the
  `μ_R,bg`-form ablation output.
- The **control outputs** of §5, each emitted as the object and its literal residual.

Every result carries its `(ε,η,σ_W)` order (`N12`; the transverse↔thickness coupling is the inherited `O(εη)`, the
self-energy the increment threaded through the bulk `Z`) and its restored `[L,T,M]` dimension. ⛔ No result is
reported without both.

---

## 5 · Independent routes and controls

⭐ Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result. Each emits the
object and its literal residual; ⛔ none asserts a target value. A **coefficient** rescale tests arithmetic; only a
**form** change tests physics.

### 5a · The ordering ablation (the non-commutation, §2)

Emit `S11CC2_CLOSED_COUPLING_KERNEL = extract(close(·))` AND `close(extract(·))` (close the already-extracted OPEN
`coupling_kernel` — a legitimate **regression** use of the imported open kernel) and their difference
`S11CC2_ORDERING_COMMUTATOR`. The commutator is the closure-induced self-energy and must be **nonzero** in general; a
byte-identical result is the tell that the closure was not threaded through the full operator. ⛔ Report the literal
difference; ⛔ state no expected form.

### 5b · The routing ablation (§1d)

Re-run the fold with the routing corrupted one direction at a time: (i) route the closed `δp_s` only into the
mechanical rows, dropping it from the θ-row flux terms; (ii) drop `Λ_X` from `t_s`. Each corruption must move the
self-energy increment (nonzero residual vs §3c). ⛔ Tests that the §1d routing is load-bearing, ⛔ not that any
channel is nonzero.

### 5c · The N6 independent-route control — the two COORDINATE constructions at a FIXED anchoring

`N6`/`N4` require **two independently-constructed routes of the same object, compared after the exact
Eulerian↔material field redefinition, then a one-sided independence corruption** — ⛔ not a bare one-sided
corruption. Per the parent pattern (`S11c_a_SHARED_PHYSICS.md` §5a "Representation-invariance routes N4/N6" +
`S11c_b_SHARED_PHYSICS.md` §5a; `S11c_decisions.md` N4/N6), the two routes are the **Eulerian and
material-coordinate constructions of the SAME self-energy increment**, with the **anchoring `α` AND the density
representative `ρ` HELD FIXED across the two routes** (⛔ the routes are the *representation* axis, ⛔ never the
anchoring axis):

```text
route 1 (Eulerian):    I_E^{α,ρ}     = extract( close(SLAB) ) − extract(SLAB) ,
route 2 (material→E):   I_{M→E}^{α,ρ} = T_{M→E}[ extract( close(SLAB_M) ) − extract(SLAB_M) ] .
```

where, at this SAME fixed anchoring `α` (and density `ρ`):
- **route 1** is c2's §3c increment built in **Eulerian** variables — the §3a substitution of the imported same-`α`
  c1 closed face response into the imported Eulerian `slab_operator` (`SLAB`), then the §3c weak extraction. ⛔
  "Eulerian graph/level-set" is S11c-a **substrate** language, ⛔ not c2's construction: c2 **folds**, it does not
  re-differentiate the interface.
- **`SLAB_M`** is S11c-b's **material-coordinate** construction of `slab_operator` at this same `α` (`S11c-b §5a`:
  `x=x(X,t)` and the S11c-a §5a face-flattening coordinate that belongs to that `α`). It is **constructed in this
  control** (the direct sibling of route 1's `SLAB`), ⛔ not an already-imported operand. ⛔ `SLAB_M` is that
  material-coordinate operator **BEFORE** the `N4` map-back; `T_{M→E}` is applied **once, to the increment**, in the
  route-2 formula above — ⛔ it is NOT already folded into `SLAB_M` (else the increment would be double-mapped).
- **`close(SLAB_M)`** is the §3a substitution of the **imported same-`α`** c1 closed face response
  `s11c_c1_face_response` into `SLAB_M` — the close operand **IS that imported response**, ⛔ NOT a second DtN
  construction, ⛔ NOT c1's Hanzawa operand standing in for "material."
- **`T_{M→E}`** is the S11c-a/S11b map + Jacobian (`N4`: `Δρ = δρ_E + u·∇ρ⁰`).
- ⛔ Do NOT reconstruct the bulk DtN by S11c-a's global scaling `w′=[w−ζ_c]/[W_bg+δW]` (`c1 §5a`: secular at
  infinity). ⛔ `Δρ` still does not bridge `LAB_HELD ↔ MATERIAL_ADVECTED`.

Emit both operands and their residual, keyed by object and anchoring (and density):

```text
S11CC2_REP_INVARIANCE_EULERIAN_OPERAND[α,ρ] = I_E^{α,ρ} ,
S11CC2_REP_INVARIANCE_MATERIAL_OPERAND[α,ρ] = I_{M→E}^{α,ρ} ,
S11CC2_REP_INVARIANCE_RESIDUAL[α,ρ]         = I_E^{α,ρ} − I_{M→E}^{α,ρ} .
```

`N6` is the physics requirement that these two are the **same operator in two representations**; the uncorrupted
residual is that representation-invariance measurement — its **computed value is the finding**, ⛔ no target value
is supplied, and the diff is adjudicated on our side (⛔ never a builder exit condition). ⛔⛔ The field redefinition
`Δρ = δρ_E + u·∇ρ⁰` relates **two descriptions of one perturbation of a fixed background**, ⛔ NEVER the two
anchorings: using it to bridge `LAB_HELD ↔ MATERIAL_ADVECTED` is a **category error** — a same-state rewrite applied
to two distinct physical setups, the exact freeze/identification `N4` exists to prevent. The two anchorings are
**distinct physical anchorings, not alternate representations of one branch** (`S11c_a_SHARED_PHYSICS.md` §2c); `c1`
§5a makes this same cut and ⛔ forbids treating `MATERIAL_ADVECTED` as the second N6 route.

**Then** the `N3`/`N4` independence probe, still at fixed `α`, `ρ`: mutate **one representation route only at its
source** and recompute it, leaving the other route unchanged. Two probes: (i) a **tilt** probe — reverse one face's
first-jet slope term in `n̂_s` on the **Eulerian** route; (ii) an **N4 advection** probe — omit the advective
`u·∇Σ_E⁰` term only from the **material→Eulerian map-back** (⛔ NOT from the common `δp_s`-independent slab base,
which cancels in the increment per §3c and so would not move it). Emit the package
`S11CC2_CONTROL_INDEPENDENCE_{BASE,CORRUPTED,RESIDUAL}[α,ρ,probe]`: **PRINT** the operands and residual, ⛔ do not
assert; the disposition — the corrupted route moved while the **uncorrupted** route did not — is adjudicated on our
side. ⛔ Do not emit an `A−A` control where emitted provenance shows the mutated source is structurally absent.
⛔ Corrupting one *anchoring* is NOT this test (it only shows two distinct setups differ, which was never in doubt).

**The cross-anchoring difference is a SEPARATE, MANDATORY contract, ⛔ not N6.** Emit it — this is the reclassified
former object, ⛔ not optional:

```text
S11CC2_ANCHORING_L_MINUS_M[ρ] = S11CC2_SELF_ENERGY_INCREMENT[LAB_HELD, ρ] − S11CC2_SELF_ENERGY_INCREMENT[MATERIAL_ADVECTED, ρ]
```

on the same footing as `S11CC2_DENSITY_LIVE_MINUS_FROZEN` (§5d): it measures whether two distinct physical setups
differ, its residual is a **computed outcome with ⛔ no prescribed zero target**, and it must ⛔ never be labeled a
representation-invariance residual.

### 5d · The background-density field-vs-field re-adjudication (rule 17, §3d.1)

Emit the self-energy increment for both density representatives (`RHO4_CONSTANT` live, `RHOBR_CONSTANT` frozen) and
their difference `S11CC2_DENSITY_LIVE_MINUS_FROZEN`. ⛔ `∇ρ→0` is not an accepted corruption (`N6`); the genuine test
is that the live-field representative differs from the frozen one wherever `ρ_br,bg⁰(x)` enters, and that the two
engines join field-against-field. ⛔ Do not manufacture a `∇ρ` jet; ⛔ do not substitute `ρ_m`.

### 5e · Three DISTINCT reduction limits (⛔ not conflated) + the μ_R,bg form ablation

⛔ `Z→0`, `Λ_A⁰=0`, and impermeability are **three different limits** (from the c1 face solve
`δp_s=Z·v_bulk,s`, `J_s=Λ_A𝒜_s+Λ_V V_s`): `Z→0` ⇒ `δp_s=0`; `Λ_A⁰=0` ⇒ `δp_s=Z V(Λ_V+ρ)/ρ ≠ 0`; **impermeable**
requires `Λ_A⁰=Λ_V⁰=0` (S11b step `:89-90`). Emit each as its own regression with its literal residual:
- **Uniform limit** (`W_bg→W̄₀`, `η→0`): the off-diagonal self-energy increment must vanish (S11b decoupling, `N6`);
  a **secondary** smoke test only (cannot see coefficient/sign/parity).
- **Zero-DtN** (`Z→0`): the increment's **bulk-nonlocal** (Z-dependent) part vanishes; ⛔ do not equate this with
  `Λ_A⁰=0` or impermeability.
- **μ_R,bg form ablation** (§3d.6): perturb the FORM of `μ_R,bg(x)` in the composed `μ_θ` operand and require the
  self-energy to move (⛔ a coefficient rescale is insufficient — only a form change tests the coupling).

---

## 6 · Method, dimensions, and script obligations

- **Method.** Balance laws + the binding material virtual-displacement rule + variational derivatives with held-fixed
  fields named + prescribed external virtual work (S11b), ⛔ never an irreversible kernel in an ordinary action. The
  weak restriction is S11c-b §3c verbatim; the operator inverse is c1 §3a/§3b verbatim.
- **Dimensions.** Restore `[L,T,M]` on every emitted object, dimensional consistency able-to-fail
  ([[feedback_dimensional_consistency_check]]); `(ε,η,σ_W)` multigrade on every object (`N12`).
- **Rest-frame limit.** Inherit `N11a` inert; c2 constructs **no** convective operator. Every c2 result inherits the
  c1/S11b smallness domain (`|q_out·v_bulk_normal_0/ω|≪1` + boundary-layer/subsonic), ⛔ never aliasing
  `v_bulk_normal_0` to `v_0`.
- **Script obligations.** The three build-skill clauses bind the build directive (`.claude/skills/build/SKILL.md`):
  a script PRINTS computed objects and never states conclusions; PRINT the residual, do not assert it; interpretation
  is the step record. ⛔ No hand-typed CAS object standing in for a computed one; every control re-enters at the
  ACTION/imported operands. ⛔ No tautological residual (rule 2 corollary 3): the §3c increment is an export
  representation, ⛔ not a check; the §3b adjointness and §3d.4 pairing residuals are emitted only when a genuine
  independent route exists.
- **Serialize CAS jobs; watch RSS.** c1 measured LIGHT (comparator peak ~317 MB on 30 GB); c2's fold threads the
  nonlocal `Z`, so the full increment may be heavy — measure the process that runs, defer heavy controls
  in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`), ⛔ never two memory-heavy CAS jobs concurrently. Detached launch
  (harness reaps `run_in_background` ~87 s). Mathematica: 2-seat licence, `--sandbox danger-full-access`.

---

## 7 · Names, F9 reservations, chain output, and export schema

**F9 / `N14` reservations.** Every new object gets a **fresh** injective `mechanical_lower_camel` name; ⛔ never reuse
an imported S11c-b/c1/S11b key for a new object. ⛔ Never reuse `slab_operator`, `coupling_kernel`,
`closure_shape_deriv`, `dtn_operator`, `dtn_kernel`, `face_response`, or any imported constant key for a c2 object.

**Chain output (`N1`/`N8`; topology = the two-leg-gated `directives/export_ledger_bind_closure_design.md` §D1–§D3).**
The SymPy engine reads the inherited model via the **positional** call
`load_model("scripts/S11c_b_exports.py", "scripts/S11c_c1_exports.py")` (signature `load_model(base_path,
*delta_paths)`, `scripts/ledger_fold.py:102`; ⛔ NOT `load_model(base=…, deltas=[…])` — that TypeErrors), binding only
its declared `IMPORT_KEYS`, and writes `scripts/S11c_c2_exports.py` as its **own-rows delta** (§D2, ⛔ not the
accumulated whole-model file). Verified on the real files: base 2441 rows + c1 delta 44 rows → fold 2485, exact-key
intersection empty, no overwrite, `check_consumer` closure resolves. What the step **binds** (§D1) is the §3
consume-set from both parents, named by their **F9 write-keys**:

- from **S11c-b**: `slab_operator`, `slab_operator_term_origins`, `mu_theta_operator`, `closure_shape_deriv`
  (S11c-a-authored, in the fold), the S11c-a T-a..T-i face substrate (`face_normal`, `conormal_deriv`,
  `face_measure_shape_deriv`, `face_velocity`, `relative_flux`, `kinematic_balance`, `traction`, `face_shift`),
  `background_density_map`, and the constants/kernels those rows reference (`Lambda_A_0`, `Lambda_V_0`, `Lambda_X_0`,
  `tau_A`, `tau_V`, `tau_X`, `rho_m`, `rho_br`, `W_0`, `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`, `mu_R`,
  `epsilon_shape`, …). ⛔ The open `coupling_kernel` is **emit-only unless c2 declares it a REGRESSION operand** for
  §5a (`export_ledger_bind_closure_design.md:148-153`) — it is ⛔ **not** a construction operand for §3c.
- from **S11c-c1**: `s11c_c1_face_response`, `s11c_c1_face_response_coeffs`, `dtn_operator`, `dtn_flat_symbol`,
  `dtn_kernel`, the per-case `s11cc1_response_resolvent_*`, the momentum symbols `s11cc1_{k,q}_{in,out}put*`,
  `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_*`. ⛔ **NOT the bare `face_response`/`face_response_coeffs`**
  (those are the S11b open/flat rows c1 imported as its regression operand; binding them folds the wrong physics —
  the N14/F9 false-equal).

⛔ The exact `IMPORT_KEYS` **root set** (minimal roots whose recursive closure covers the consume-set) is fixed at the
c2 build directive against the two real export files, ⛔ not enumerated-then-frozen here; its two decision legs verify
it, and that the guard (`check_consumer`/`assert_lookups_equal_manifest`/`assert_delta_is_minimal`) passes on the
two-parent fold — ⚠ noting the guard passes on **key existence**, so it will **not** catch a wrong-provenance
`face_response` binding; that check is the directive's + legs' responsibility. `BUILD_INPUT_DIGESTS` pins, per §D3,
`{this sub-step's SymPy audit, scripts/S11c_b_exports.py, scripts/S11c_c1_exports.py, this spec, scripts/ledger_fold.py}`.
⛔ Never `git add -f` a big `.out`; ⛔ never annex an `*_exports.py`.

**The comparator (`N8`, frozen `T7` contract).** The c2 comparator joins the two blind engines' emitted objects by
name, pairs residual operands, is three-valued, rejects a native boolean, and PRINTS/decides nothing (rule 2). Its
load-bearing residual is on the **self-energy increment** (§3c) — which drops the deferred ≥64 GB S11c-b bulk residual
but **carries the face-force/#90 closure-fold sign conventions** (§3c); the comparator **SURFACES** those (rule 1/6),
⛔ does not normalize them, and the §3d.4 mechanical-power pairing adjudicates the face-force sign. It also surfaces
the §3d representation questions (density field-vs-field, `t_s` scalar-vs-4-vector, DtN whole-form-vs-kernel,
flat-symbol usage); the reconcile is the staged representational bridge
([[feedback_reconcile_representational_bridge]]), ⛔ never a blanket collapse. The four c1 giant families + the full
per-family symbolic residual remain deferred (`DEFERRED_HEAVY_RUNS.md`); c2 names, ⛔ does not pre-adjudicate,
whatever it cannot close on this box.

**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
rows it folds into, and the c1 closed response — importing nothing (the only cross-engine control). ⛔ The denylist
stays cut (`N9`/rule 12); blindness is enforced by absence.

---

## 8 · Supplied versus computed; builder report

**SUPPLIED (unfalsifiable in this build):** all of §1 (the two imported models and their AGREE/UNDECIDED disposition,
the slab-operator rows and the real θ-row structure, the face closure laws and Λ-channel routing), the §2
non-commuting-ordering framing and its counterexample and the three named objects, the S11c-b §3c weak-restriction
convention, the c1 §3a/§3b operator inverse, `background_density_map`, `N11a`, `N12`.

**COMPUTED (outputs, ⛔ none stated here):** the closed two-face operator (§3a); the re-extracted closed coupling
kernel and its adjointness disposition (§3b); the self-energy increment and its operands (§3c); the six
re-adjudication objects (§3d); every control residual (§5); the `(ε,η,σ_W)` orders and `[L,T,M]` dimensions.

**Builder report.** The build directive states, per emitted object, which line computed it (`.claude/skills/build`);
declares the fold's symbol map (§3a), routing (§1d), ordering (§2), and increment representation (§3c) it implemented;
and reports the literal residuals of §5 — ⛔ never a prose conclusion. The disposition of the six §3d re-adjudications
and every §5 residual is read on **our** side, in the step record, ⛔ not asserted by the script (rule 5).
