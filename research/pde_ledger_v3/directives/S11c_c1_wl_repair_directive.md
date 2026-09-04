# S11c-c1 Wolfram engine — repair directive

## Role and single deliverable

Modify `research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` in place (baseline
`e139bc61`). Its only product is the flushed stdout tag stream (regenerated). ⛔ Import nothing (the blindness
control stands — no `Get`/`Import`/`<<`/`ReadString`/`OpenRead`, no absolute repo path, no SymPy engine/export/
`.out`); the copy-to-empty-dir + structural-scan acceptance from the build directive still holds. This repair
fixes three defects (R1–R3, code-verified) and folds four biting NITs (R4). ⛔ Do **not** touch the sound core.

`CLAUDE.md` binds. The sole physics authority is `directives/S11c_c1_SHARED_PHYSICS.md` (§§0–8; HEAD `f90e7630`);
the build directive `directives/S11c_c1_wl_build_directive.md` governs naming/blindness/idioms unchanged. Add no
expected value or acceptance criterion (rule 5): each fix names the **object** the emit must be and the
**construction invariant** its control must satisfy, never a computed value, sign, or residual.

## The three defects — SUPPLIED, code-verified against baseline `e139bc61`; each fix re-enters at the CONSTRUCTION

### R1 — `DTN_OPERATOR` composition must carry BOTH momentum legs (spec §3a `SHARED_PHYSICS.md:247-261`; rule 17)

**Defect** (`operatorCompositionFromDerivation`, `.wl:542-560`): `gZero = FourierMultiplier[I qOut[omega,
momentumOutput]]` is used in BOTH positions of `OperatorComposition[gZero, multiplication, gZero]`, and `nZero`
likewise in `OperatorComposition[nZero, dtnVariation, nZero]` — every `N_0`/`N_0^{-1}` factor carries
`momentumOutput`; `momentumInput` appears nowhere in the composition. In the two-sided composition `N_0 ∘ M_h ∘
N_0`, the **rightmost** factor acts on the INPUT field first, so it must carry the **input** momentum; the
leftmost acts on the OUTPUT. As emitted, the input-leg `N_0` is frozen to the output momentum — the WKB/
left-quantized freeze spec §3a forbids, relocated into the operator emit, and inconsistent with the verified
two-leg `DTN_KERNEL` (`.wl:478-479` `qOutputLive`/`qInputLive`).

**The object the fix must produce:** the composition form must **manifestly carry both legs** — the momentum a
`FourierMultiplier` applies is the momentum of the field it acts on, so the rightmost `N_0`/`N_0^{-1}` factor of
each `OperatorComposition` carries `momentumInput` and the leftmost carries `momentumOutput` (and any
intermediate `M_h`/`Div(h∇)`/`κ²h` carries the transfer `Ŵ_bg(momentumOutput−momentumInput)`), so that a
consumer expanding the composition to a kernel recovers `DTN_KERNEL`. **Construction invariant (emit a probe as a
computed object, ⛔ not a bare check):** a probe over the emitted `DTN_OPERATOR` returns
`COMPOSITION_HAS_Q_INPUT = True` **and** `COMPOSITION_HAS_Q_OUTPUT = True`; and the composition, kernel-expanded
on-shell, minus `DTN_KERNEL` is a computed residual that is zero for the correct labelling and **moves** if the
rightmost factor is (as a control) re-frozen to the output leg. ⛔ Do not alter `DTN_KERNEL` (verified two-leg).

### R2 — energy BULK operand must be the outgoing far-field acoustic Poynting flux of `φ` at `|w|→∞` (spec §3b `:321-330`)

**Defect** (`deriveEnergyOperands`, `.wl:994-1028`): `bulkAmplitude` is solved from the FACE Neumann law
(`I qOut·A == energyVelocity`, `:1008-1010`), and `farFieldPressure = pressureAtFace·farFieldPhase` (`:1015`),
`farFieldVelocity = energyVelocity·farFieldPhase` (`:1017`) — so `OUTGOING_FARFIELD_POYNTING_FLUX =
measure·Re[pressureAtFace·Conj[energyVelocity]·|farFieldPhase|²]/2`, with `farFieldPhase` a free undefined
placeholder. At `farFieldPhase→1` the bulk operand equals the face pairing byte-identical (residual 0): it is the
`½Re(δp·V*)` face quantity spec §3b **forbids**, ⛔ not a Poynting flux at `|w|→∞`. The two "independent routes"
are one.

**The object the fix must produce:** the bulk operand is the **outgoing far-field acoustic Poynting (energy)
flux of the bulk potential `φ`**, computed from the **half-space outgoing solution** `φ` (the §1b
radiation-selected outgoing wave, carrying `q_out` and the first-shape curvature correction) evaluated on a
control surface **in the half-space / at `|w|→∞`** — `S = ½ Re[δp · v_bulk,n*]` with `δp = −ρ_m ∂_tφ`,
`v_bulk,n = ∂_nφ`, evaluated on that far surface, ⛔ never `δp`/`v` at the face and ⛔ with no free `farFieldPhase`
factor. Remove `farFieldPhase` entirely. **Construction invariant:** the bulk operand and the face traction
pairing are two constructions from independent data, so **a one-sided corruption of ONLY the bulk route** (e.g. a
`q_out` branch flip in the far-field outgoing solution) moves ONLY the bulk operand and drives `ENERGY_RESIDUAL`
nonzero, while the face operand is unchanged; and the residual is **not** a structural zero killed by any
`farFieldPhase→1` substitution. (Same fix the SymPy engine's R1 applied.) Emit operands + residual; ⛔ assert
nothing.

### R3 — energy FACE operand must be built from the response's actual `t_s` (spec §3b `:322-330`)

**Defect** (`.wl:1011-1012`): `tractionVector = tractionOrientation·pressureAtFace·normal` uses a fresh
impermeable flat solve and carries **no** `Λ_X`; it never binds `FACE_RESPONSE`'s `t_s`
(`FACE_HAS_TS_FROM_RESPONSE=False`, `RESPONSE_TRACTION_HAS_LAMBDAX=True`). So flipping the reconstruction's own
orientation moves the pairing (a self-test), but flipping the **emitted** response `t_s` cannot — a `t_s`-sign or
`Λ_X`-placement error in the closed response is invisible to the audit, which spec §3b says must catch it.

**The object the fix must produce:** the face operand is the true-area traction pairing built from the **§3b
`t_s` object of the closed permeable response** (`t_s = −(δp_s + Λ_X 𝒜_s) n̂_s`), so it **binds** the same `t_s`
the response emits (in the real-`ω`/propagating/impermeable/`Λ_X⁰=0` slice it reduces to `−½Re Σ_s∫a_s^0 δp_s
V_s*`, but the operand must be that response `t_s`, ⛔ not a fresh solve). **Construction invariant:** a one-sided
sign flip of the **response** `t_s` (at its source in the response construction) moves ONLY the face operand and
drives `ENERGY_RESIDUAL` nonzero (`FACE_HAS_TS_FROM_RESPONSE=True`); the bulk operand (R2) is unchanged. Emit a
computed probe object showing the face operand references the response traction.

## R4 — fold the four biting NITs (both legs)

- **`SOURCE_EQUATIONS`** must emit the **real** boundary/closure equations as re-parseable CAS relations, ⛔ not
  `HoldForm[{flatBulk$NNNNN,…}]` Module-local gensyms (the `Solve` already uses the real equations — keep them in
  the payload).
- **`PERMEABLE_PORT_HERMITIAN`** parity keys `DELTA_W`/`ZETA_C` must carry the **actually computed** `δW`- and
  `ζ_c`-combination blocks (spec §3a: "whether it couples the `δW` and `ζ_c` combinations — a computed block"),
  ⛔ not the same per-face matrix duplicated under both keys. If the two combinations genuinely coincide at this
  order, emit that as a computed equality (operands + residual), ⛔ not as an unused loop.
- **`§5a` layer-potential**: remove the decorative typed head `RadiationPreservingLayerPotential[…]` that carries
  no data dependence; the emitted Hanzawa/layer-potential operand must be the **computed** second-route kernel
  (the one whose one-sided corruption already moves `REP_INVARIANCE_RESIDUAL`), so the operand and the residual
  reference the same computed object.
- **Locus `REAL_ADMISSIBLE`**: the `STATUS_TOKEN∈{ADMISSIBLE·EXCLUDED·UNDECIDED}` must come from a **genuine**
  real-admissibility test of each solution branch (per spec §6), ⛔ not from `Head[reduced]=!=Reduce` labelling
  any non-`True`/`False` result `ADMISSIBLE`. An undecidable branch is `UNDECIDED` (a coverage finding), ⛔ never
  silently `ADMISSIBLE`. Keep `TEST_OBJECT` the live re-parseable CAS boolean.

## What must stay byte-identical (⛔ do not touch — the 2-leg-sound core)
`DTN_KERNEL`, `DTN_FLAT_SYMBOL`, `DTN_RIGID_SHIFT_*`, `DTN_BY_REGIME_PAIR`, `DTN_BY_PARITY`,
`DTN_HERMITIAN_PART`/`DTN_REACTIVE_PART` (their true-area adjoint is correct), `DTN_GRAZING_BEHAVIOUR`, the
permeable `FACE_RESPONSE`/`FACE_RESPONSE_COEFFS` operator inverse, the `NONINVERTIBILITY_CONDITION`, the §5b/§5d/
§5e controls and their liveness, the T-a..T-i re-derivation, the reserved-name spellings, and `μ_θ` opacity.
R3 binds the existing `FACE_RESPONSE` `t_s` — it does not re-derive it.

## Method and script clauses (spec §6 — unchanged, binding)
Every fix re-enters at the construction (the composition operator, the bulk `φ` solution, the response `t_s`),
⛔ never at a result. Each control emits operand A, operand B, and `A−B` before any guard; a physics disagreement
emits and continues (exit 0); nonzero exit is operational only. ⛔ No `VERDICT`/`PASS`/`FAIL`; a boolean is a
typed CAS object (unevaluated relational / `STATUS_TOKEN`), ⛔ never a native `True`/`False` residual operand. ⛔
No hand-typed CAS object with no data dependence (a probe like `COMPOSITION_HAS_Q_INPUT` is computed from the
emitted operator, ⛔ not asserted). Run discipline unchanged (detached; one kernel; `danger-full-access`;
600s-no-output/RSS kill recorded-not-narrowed; per-case streaming).

## Builder report (under 25 lines)
The diff vs `e139bc61` (functions changed); confirmation `DTN_KERNEL`/flat symbol/core rows are byte-identical;
the three repaired controls' probe objects (`COMPOSITION_HAS_Q_INPUT`, the bulk-only and face-only one-sided
corruptions moving only their operand, `FACE_HAS_TS_FROM_RESPONSE`); the four NIT fixes; tag count unchanged (51);
tasks run; runtime; any ambiguity. ⛔ No computed physics value or conclusion.
