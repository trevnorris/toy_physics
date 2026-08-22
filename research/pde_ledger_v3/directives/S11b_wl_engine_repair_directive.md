# S11b WL engine — repair directive (folded once after two directive legs)

## Authority and boundary
Repair `research/pde_ledger_v3/mathematica/S11b_interface_coupling_law_mathematica_audit.wl` IN PLACE
(committed baseline `ec89f9df`). `CLAUDE.md` binds; `directives/S11b_SHARED_PHYSICS.md` is the sole physics
authority (⛔ do not re-open). The engine is BLIND (imports nothing, re-derives). ⛔ Add no expected value or
acceptance criterion (rule 5): the orchestrator holds the withheld references and diffs on our side. ⛔ Do
NOT state or bless the current non-target values as correct (that both leaks them as references and can
protect an undetected defect); confine "leave it unchanged" to the NAMED non-target objects listed under
Scope. The derived physics at the non-target sites must stay behaviourally unchanged.

## F-WL-1 — `ZPERM_SLICE_MAP` (B2c): wrong extraction (sign) and wrong level (dynamic, not static)
Site: `emitShared["ZPERM_SLICE_MAP", ...]` (~L865); `equationZeroForm` L367; the `μ_s=0` slice of
`generalFace["EQUATIONS"]`.
1. **Extraction (sign).** The current `rawPressureCoefficient` is `Coefficient[equationZeroForm[flux_eq],
   facePressure]` — the coefficient of `facePressure` in `faceFlux − RHS`. Fix: extract the coefficient of
   `facePressure` from the **flux expression itself (the RHS of the constitutive closure)** — the
   coefficient the supplied §4 affinity forces on the raw-pressure channel — ⛔ not from the residual
   `faceFlux − RHS`. ⛔ Do not state the resulting sign or value.
2. **Level (static).** That RHS coefficient is still DYNAMIC while the RHS carries `Λ_A(ω)`. B2c asks for the
   STATIC `Λ_p⁰ ↔ Λ_A⁰` map — the coefficient on the ω-independent `Λ_I⁰` law (equivalently the `ω→0` / the
   `Λ_I(ω)→Λ_I⁰` reduction of the RHS coefficient). Emit that ω-independent object as `S11B_ZPERM_SLICE_MAP`.
   Any dynamic relation, if emitted, goes under a DISTINCT tag. ⛔ Do not conflate the static map with a
   frequency-dependent coefficient.
⚠ `S11B_ZPERM_SLICE` (the specialized face response) is already correct — leave it; only the MAP is wrong.

## F-WL-2 — the §6 energy checks are tautological / detached; make them genuine with NAMED independent routes
Sites: `PRESSURE_WORK_SIGN_CHECK` (~L1176; `slabPressurePower = −½Re[X]`, `outgoingBulkPower = +½Re[X]` from
the SAME `X` ⇒ residual ≡ 0), and the two-port / energy-sink/source / unattributed checks consuming the
hand-written `energyEquationRules` (~L1118), not the assembled system.
Fix: build the power balance from **two INDEPENDENT routes named by §6**, ⛔ not both from `fullSystem`
(which is the SOLVED system — pairing two of its objects re-creates a tautology). Route A (slab side): the
§6 discriminator-2 virtual-displacement pairing on the **PRE-elimination / off-shell** EOM and face laws
(in-plane eq with `(∂_t u)*`, thickness eq with `(∂_tδW)*`, mass balance with `μ_s`). Route B (bulk side):
the outgoing bulk acoustic power from the supplied §2 acoustic ansatz / B0b impedance `Z`, derived
independently. The residual of A + B must be able to FAIL under a sign/form change to the EOM or the closure.
⛔ `NOT_ESTABLISHED` is FORBIDDEN for §6 discriminators 2–3 (§6 supplies both routes); it is permitted ONLY
for a NAMED, demonstrably-absent operand, emitted as available-operands + `MISSING_OPERAND`, with ⛔ no
residual and ⛔ no `TEST_OBJECT`.

## F-WL-3 — three distinct checks (split)
- **3a `KERNEL_ORIENTATION_IDENTITIES`** (~L1102) extracts kernels from `closureLawForExtraction` (a shadow);
  and `solveFace` (L298–299) has already EXPANDED `𝒜`, so a naive `Coefficient[assembled RHS, 𝒜]` is 0 (false
  pass) and a naive `kernelReplacement` moves non-target `KERNEL_PROPAGATION_RESIDUALS`. Fix: introduce ONE
  raw face-law association (with `𝒜` restored as an un-expanded symbol on the PRE-elimination closure)
  consumed by the solve, the traction assembly, AND all three kernel extractions; extract the kernel as the
  §6 `coeff_𝒜(J)|_{V=0}` from that pre-elimination closure so the identity bites when a kernel is
  mis-oriented. State explicitly whether `KERNEL_PROPAGATION_RESIDUALS` and any artifact tag legitimately
  change as a result (that is expected, not a regression).
- **3b `CAUSALITY_CHECK`** (~L1105) is unconditionally true: its `Cases` returns `{}` and `And[] = True`. Fix:
  use explicit association lookups; require the expected NUMBER of nonempty records to be present before
  aggregating (a positive presence check), so a missing/renamed record FAILS rather than silently passing.
- **3c `GRAZING_MODE_CLASSIFICATION`** (~L1368): `Limit[q·fullSystem["MATRIX"], q→0]` zeros a finite matrix
  (rank 0). Fix: classify from the **leading `q⁰` block** of the matrix pencil as `q→0` (which is
  non-degenerate (⚠ SUPERSEDED by fix-2 a6b0b684: the actual sound-cone block is rank 3 = non-grazing, NOT rank 2)), via a Laurent/nullspace analysis of that leading-order block, and map the
  leading-order stratum to the mode categories. Require a non-degenerate leading-order payload. ⛔
  `NOT_ESTABLISHED` here is bounded to "no leading-order object exists," ⛔ never "cannot classify."

## The three script clauses + corollaries (verbatim, non-negotiable)
1. PRINT computed objects; ⛔ never a prose conclusion. 2. PRINT the residual; ⛔ do not `assert` it zero.
3. Interpretation → the step record. ⛔ A hand-typed CAS object is still hand-typed: every emitted object is
REACHED BY COMPUTATION from the assembled system; a control re-enters at the ACTION/closure, ⛔ never at a
result. ⛔ No tautological residual — the two operands come from INDEPENDENT routes, verified by one-sided
corruption (break route A only; route B must not move).

## Run discipline (Mathematica)
One kernel at a time (2-seat licence); orchestrator writes the committed `.out` after review; demonstration
runs to scratch. Kill a demo kernel at 600 s no-new-output or RSS > 6 GB (⚠ this engine runs ~0.5 GB; 6 GB is
a runaway tripwire, not a budget).

## Acceptance — executable, DECISIVE, no expected values (rule 5)
- **F-WL-1 (representation-invariance, not "moves under the affinity form").** ⚠ An affinity-form change moves
  BOTH the wrong residual extraction and the correct flux extraction, so that test is NOT decisive.
  Demonstrate instead: (i) corrupt `equationZeroForm` / swap the closure's LHS↔RHS — the OLD residual
  extractor's map MOVES, the NEW flux-RHS extractor's map does NOT (representation invariance); and (ii) the
  emitted `S11B_ZPERM_SLICE_MAP` is independent of `ω` and `τ` (static). ⛔ Do not state the withheld relation.
- **F-WL-2/3:** each repaired check must move under a one-sided corruption of exactly its own load-bearing
  object (route A only for F-WL-2; the assembled kernel for 3a; a removed record for 3b; the leading-order
  block for 3c), and must NOT move under a corruption of an unrelated object.
- The named non-target objects (Scope) stay behaviourally unchanged; the orchestrator diffs the withheld
  criteria on our side.

## Scope — objects that must NOT change behaviourally (name them; ⛔ do not bless their VALUES as correct)
Impedance / regimes, added mass, `S11B_ZPERM_SLICE`, the transverse mode, the breathing slice, the
longitudinal dispersion, and the mode roots. (F-WL-3a may legitimately move `KERNEL_PROPAGATION_RESIDUALS`.)

## Report (§13) — under 20 lines
The edits (tag:line), which checks are now genuine vs honestly `MISSING_OPERAND`, the one-sided-corruption /
representation-invariance demonstration for each, and confirmation no Scope object changed and no expected
value was stated.
