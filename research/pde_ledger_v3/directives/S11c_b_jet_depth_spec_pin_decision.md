# Spec-pin decision list — what IS the S11c-b slab U-momentum row (θ-independent EL, or θ-eliminated)?

## Why this exists

The cross-engine `scripts/S11c_b_row_residual.py` will show a disagreement on the bulk U-momentum row: WL carries
order-3 background jets there, PY order-2. A three-way reconciliation (orchestrator + Codex + Grok; record
`_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`, advisor scripts+stdout under
`_measurements/s11c_b_jet_depth_consult_{codex,grok}/`) established the mechanism with three consistent
computations. It is NOT a jet-depth freeze and NOT a physics error in either engine. It is a REPRESENTATION
mismatch, and the ledger must PIN which representation the slab U-row is before either engine is changed and before
`row_residual` is read. This list states the decision, the settled facts, the two objects, a proposed resolution,
and the consequence of each branch. Per rule 7 it gets two decision legs (Codex + Grok) before any engine change.

## Settled facts (measured on BOTH engines — GIVEN; do not re-derive, but you may spot-check)

1. The §3a energy basis carries only FIRST background jets (`∂W_bg`, `∂μ_R,bg` spurions); no Hessian is an energy
   invariant. ⇒ the `∂u`-flux coefficient `∂L/∂(∂_i u)` is background-jet order 1.
2. The RAW held-fixed variational derivative `δU/δu` (θ, e_W held fixed) — one EL divergence on an order-1 flux —
   is background-jet order **2** on BOTH engines. Retained-grade depth 3/4 generate nothing new (PY measured
   identical at depths 2,3,4; WL's raw energy EL measured order 2).
3. The order-3 in WL's emitted `U_INTERNAL` is ENTIRELY the θ-constraint reaction: WL's
   `constrainedRowsWithLiveEnergyEL` eliminates the virtual θ through the material virtual constraint `δ_vΣ_mat=0` (linearized)
   (`θ + (W0/W) e_W + u·∇W/W + Div u`), leaving `Inactive[Div]` of a flux that IS `μ_θ`; activating that outer Div
   with `W_bg` live adds one derivative → `∇μ_θ`, background-jet order **3** (10 atoms). PY reproduces the SAME 10
   order-3 atoms when the same constraint is applied — so the engines agree on BOTH objects (⚠ measured as matching
   max jet order + the same 10 order-3 atom NAMES on the probed case/term — one PY term, WL representative 16; NOT
   verified as full coefficient equality across all branches/densities); they merely emit DIFFERENT ones as the U-row.

## The two candidate objects (they differ by exactly the order-3 `∇μ_θ` constraint reaction)

- **(A) θ-INDEPENDENT.** The slab U-momentum row is the held-fixed variational derivative `δU/δu|_{θ,e_W fixed}`
  (PY `operator_from_density`). θ is an independent DOF carrying its OWN row `μ_θ`; the constraint is NOT folded
  into the U-row. Bulk U-row max background-jet order 2.
- **(B) θ-ELIMINATED.** The slab U-momentum row is the constraint-reduced momentum equation: S11b's binding
  virtual-displacement rule eliminates virtual θ via the material virtual constraint `δ_vΣ_mat=0` (linearized), folding `∇μ_θ` into the U-row
  (WL `constrainedRowsWithLiveEnergyEL`). Bulk U-row max background-jet order 3.

## Proposed resolution — (A), θ-INDEPENDENT — TO BE ADVERSARIALLY VERIFIED, not assumed

The proposal is (A): the slab U-momentum row is `δU/δu` with θ held fixed, and WL should stop folding `∇μ_θ` into
`BULK_ENERGY`/`U_INTERNAL`. Reasoning to CHECK against the sources (⛔ do not take these as established — read the
governing text and confirm or OVERTURN):
- §3b (`S11c_b_SHARED_PHYSICS.md:272-288`) defines the slab operator as "the first-order equations of motion for
  the slab DOFs `{u,θ,e_W}`" and emits `S11CB_MU_THETA_OPERATOR (μ_θ = δU/δθ … kept as a named operand)`
  SEPARATELY — i.e. θ appears to have its own row, and folding `∇μ_θ` into the U-row would double-represent it.
- §1c (`…:126-131`): `μ_θ ≡ (δU/δθ)|_{u, e_W and all other fields fixed}` and "θ may not be eliminated through a
  constraint before this derivative."
- BUT §1c (`…:134,146-148`) also says the equations of motion use "S11b's method — balance laws, the binding
  virtual-displacement rule, … and prescribed external virtual work," and in S11b `μ_θ` drives the FACE mass-flux
  affinity `𝒜_± = μ_s − δp_±/ρ_m`, `μ_s = μ_θ/ρ_br⁰` (`S11b_SHARED_PHYSICS.md:~228`). If that binding
  virtual-displacement rule constrains δθ (the mass balance is a constraint on the slab-DOF virtual
  displacements), the reduced U-equation legitimately carries `∇μ_θ` — that is object (B).

The crux the legs must settle FROM THE SOURCES: **does S11b's binding virtual-displacement rule leave the slab-DOF
virtual displacements `{δu, δθ, δe_W}` INDEPENDENT (→ one held-fixed EL row each, θ's row = μ_θ; object A), or does
it CONSTRAIN δθ via the mass balance so the momentum equation absorbs `∇μ_θ` (→ object B)?** Also: is the
divergence-form requirement of §3b ("variable coefficients sit inside the in-plane divergences") satisfied by (A),
or does it REQUIRE the `∇·(flux)` structure WL produces?

## Consequences (spell out so the follow-on builder has an unambiguous target)

- **If (A) is pinned:** WL is the engine to change — its `U_INTERNAL`/`BULK_ENERGY` slot must NOT activate the
  constraint `Div` that produces `∇μ_θ`; that content belongs to the separately-emitted `μ_θ` row. PY unchanged.
  PY's `STRONG_ROW_JET_DEPTH = 2` STAYS (raising it is a measured no-op). After the fix, both engines emit an
  order-2 bulk U-row and `row_residual` should agree there.
- **If (B) is pinned:** PY is the engine to change — it must fold the same mass-balance constraint into its U-row
  and raise the strong-row background-jet depth to 3 (then, and only then, PY depth-2 is a genuine freeze). §3b
  must be amended to state the U-row is the constraint-reduced equation and to reconcile that with the separate
  `μ_θ` operand (avoid double-representation).
- **If the sources are genuinely AMBIGUOUS:** the pin is itself a SPEC finding — §3b/§1c must be amended to state
  explicitly whether the slab momentum row is θ-independent or θ-eliminated. Name the exact sentence to add. (A
  closed spec is not a correct spec.)

## What each leg must do (document/spec adjudication — read the sources FIRST, quote both sides)

Read and quote the governing text: `S11c_b_SHARED_PHYSICS.md` §0/§1a/§1c/§3a/§3b (esp. `:52-64,126-148,272-288`);
`S11b_SHARED_PHYSICS.md` (the balance-law method, the virtual-displacement/virtual-constraint rule, the face
closure/affinity, ~`:195-232` and the virtual-work sections); `steps/S11b_interface_coupling_law.md` (the coupling
law RESULT — is θ an independent DOF or slaved?); the #88 KINETIC+θ treatment
(`directives/S11c_b_88_blast_radius_build_directive.md`). Then:
1. State whether the settled facts above are accurate (spot-check permitted; the advisor scripts are under
   `_measurements/s11c_b_jet_depth_consult_{codex,grok}/`).
2. Determine, WITH CITATIONS, whether §3b + S11b's method specify object (A) or (B) for the slab U-momentum row —
   or whether the spec is genuinely ambiguous. ⛔ Do not defer to the proposed resolution; overturn it if the text
   says (B). ⭐ The decisive question is the independence/constraint status of the slab-DOF virtual displacements
   in S11b's binding rule — ground your answer in that rule's text, not in the jet-order measurements.
3. State the consequence for the engines under your verdict, and flag any conflict with #88 or the S11b coupling
   law that the pin must stay consistent with.

Output: your verdict (A / B / AMBIGUOUS-needs-spec-addition) with the governing citations quoted, whether the
settled facts hold, the engine consequence, and any consistency conflict. A prose claim without a quoted spec
citation (or, for a spot-check, a script + literal stdout) is discarded.

---

## FOLDED VERDICT — PINNED (B). The proposed (A) is OVERTURNED. (both legs + orchestrator verification)

Both decision legs (Codex + Grok, raw verdict transcripts `~/.s11_build/S11c_b_jet_depth/{codex,grok}_leg_verdict.log`;
`.log` is gitignored, so the distilled verdict + citations live in this section) independently
returned **B** and overturned the proposed (A), citing the SAME governing S11b text. The orchestrator verified the
decisive citations verbatim against the source (rule 13). **PIN: (B) — the slab U-momentum row is the
constraint-reduced in-plane equation carrying `−∇μ_θ`; θ's virtual displacement is NOT independent.**

**Decisive governing text (orchestrator-verified verbatim):**
- `S11b_SHARED_PHYSICS.md:341-342` (constraint eq at `:337`): "Thus `δ_vθ, δ_v(δW), δ_vu` are not independent. ⛔ **Do NOT vary `U` with `θ`
  held fixed.** Impose the VIRTUAL CONSTRAINT either by eliminating one virtual variation or by a Lagrange
  multiplier… The **same** multiplier supplies the in-plane restoring force and the thickness term."
- `S11b_SHARED_PHYSICS.md:426` (convention check a): "**The in-plane equation your variation produces must carry
  the restoring force `−∇(δU/δθ)`. This single check selects the convention uniquely.** … **Wrong derivation
  caught:** … or **varying at fixed `θ`, removes this contribution.**"
- S11b's OWN engine encodes it: `S11b_interface_coupling_law_mathematica_audit.wl:280`
  `constrainedUL = Expand[explicitUL + I k muTheta]` — the in-plane EOM carries the `μ_θ` gradient (verified).

**Why (A)'s citations do not select (A) (both legs agree):** §1c's "θ may not be eliminated through a constraint
before this derivative" scopes the CONSTITUTIVE `μ_θ` (compute it held-fixed first), not the U-row EOM. §3b's
separate `S11CB_MU_THETA_OPERATOR` is a scalar constitutive/face-affinity operand (`𝒜_± = μ_s − δp/ρ_m`), not an
independent `θ = 0` equation — so folding the VECTOR `∇μ_θ` into U is not double-representation. θ's actual EOM is
the sourced MASS-EVOLUTION row (restored after the variation), a different object. §3c's independent θ trials probe
the already-built §3b operator; they do not unbind the forming variations. Divergence form fits either object.

**Additional consistency (Grok):** pinning (A) would make the S11c-b U-row a different object from S11b's
`INPLANE_EOM`, so the §5c uniform-limit residual would carry the constraint reaction even at `(η,σ_W)=0`, breaking
the inherited method and the §5c regression. ⇒ (A) is inconsistent with §5c; (B) is consistent.

**ENGINE CONSEQUENCE (the follow-on builder's target — needs its OWN build legs, rule 7/9):**
- **PY changes.** `operator_from_density` (`sympy_audit.py:1968-2062`) currently emits the held-fixed EL (object A)
  — correct for the CONSTITUTIVE `μ_θ` operand, WRONG for the slab U-row. PY must (i) fold S11b's virtual
  constraint `δ_vθ + δ_ve_W + ∇_x·δ_vu = 0` into `U_MOMENTUM_ROWS` and the thickness row (same multiplier),
  (ii) raise `STRONG_ROW_JET_DEPTH` 2→3 (only WITH the constraint reaction present is depth-2 a genuine freeze;
  PY already reproduces the 10 order-3 atoms once the constraint is applied), (iii) keep `μ_θ` as the separate
  constitutive operand, (iv) treat the θ-row as MASS-EVOLUTION, not `μ_θ = 0`.
- **WL is already correct** (`constrainedRowsWithLiveEnergyEL`); no WL change. #89b's un-freeze was RIGHT.
- **§3b amendment** (add after `:276`, both legs proposed near-identical text): the U-momentum row is S11b's
  in-plane momentum balance under the binding virtual-displacement rule (virtual θ eliminated via `δ_vΣ_mat=0`
  linearized), carrying `∇μ_θ`; it is NOT `δU/δu|_{θ,e_W fixed}`; `S11CB_MU_THETA_OPERATOR` remains the
  constitutive operand, not that reaction and not the mass-evolution equation.
- **#88 re-adjudication:** #88 identified the strong U-rows with the raw held-fixed EL (`operator_from_density`);
  that identification conflicts with this pin. #88's energy-basis-completion disturbance measurement STANDS (it
  witnesses the unconstrained-EL ingredients), but the KINETIC-family/full-row adjudication must be REDONE once PY
  carries the constraint reaction, and #88's θ result is a disturbance of `MU_THETA_OPERATOR`, not of an
  independent θ equation.
- **row_residual** then compares the two constraint-reduced (order-3) U-rows; emit the raw EL under a distinct
  name (`RAW_BULK_U_EL`) only if still useful.

⛔ No engine change yet — the PY constraint-fold is the follow-on BUILDER and gets its own decision list + build legs.
