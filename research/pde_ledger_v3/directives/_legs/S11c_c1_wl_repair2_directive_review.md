# Directive design review — S11c-c1 Wolfram engine REPAIR-2 directive (rule-7 gate, BEFORE the build)

## What you are reviewing

`research/pde_ledger_v3/directives/S11c_c1_wl_repair2_directive.md` — an **orchestrator-written** repair directive
for the blind Wolfram engine of S11c-c1. It fixes **two emit-only defects** that a retroactive gate found had
propagated into the built engine `research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl`
(baseline `13f0bd2c`). **This review runs BEFORE the repair-2 build** — nothing has been built yet against this
directive. Your job: judge whether the **directive itself** is sound — the two defects it names are real in the
code, the fix it prescribes is the correct physics at the correct construction site, it leaks no expected value,
it mis-specifies no control, and it correctly bounds what stays byte-identical. A defect you catch here is caught
before it can propagate; that is the entire point of this gate.

## Sources of truth — form your own view FIRST, then read the directive

- The current engine code at the two defect sites (read these ranges and judge for yourself whether the defect is
  real): `mathematica/S11c_c1_bulk_closure_mathematica_audit.wl`
  - **R1 site:** `fredholmFunctionSpaceOperator` at `.wl:580-597`, and its consumer `NONINVERTIBILITY_CONDITION`
    at `.wl:1537-1549`. Compare against the already-repaired two-leg operator
    `operatorCompositionRecordFromDerivation` at `.wl:542-572` (note `gZeroOutput`/`gZeroInput`,
    `nZeroOutput`/`nZeroInput`). Is every `N_0` factor in `fredholmFunctionSpaceOperator` really frozen to
    `momentumOutput`, with `momentumInput` (and any genuine dummy) absent? Does the `.wl:577-579` comment claim a
    "dummy coordinate … rather than either external kernel leg" while the code uses `momentumOutput`?
  - **R2 site:** `portParityCombination` at `.wl:1608-1650` and the emission at `.wl:1651-1706`. Is
    `DELTA_W_PORT_MATRIX` really `(plus+minus)/4`, `ZETA_C_PORT_MATRIX` really `(plus+minus)`, and the difference
    `(plus−minus)/2` really placed in the off-diagonal coupling block? `portCaseForAxes` (`.wl:1588-1593`) gives
    the per-face objects for face ∈ {+1,−1}.
- The sole physics authority: `directives/S11c_c1_SHARED_PHYSICS.md`.
  - §1 `:76-83` — the face DOF definitions: `δW ≡ ζ₊ − ζ₋`, `ζ_c ≡ (ζ₊ + ζ₋)/2`, `ζ_s = ζ_c + s·δW/2`, `s∈{+1,−1}`.
  - §3a `:274-277` — emit the face-parity structure, whether the first-shape correction couples the `δW` and `ζ_c`
    combinations, as a **computed** block.
  - §3b `:299-306` — the Fredholm/noninvertibility condition is failure of `[I + (Λ_A/ρ_m²) Z]` to be invertible
    (Z the two-momentum operator, §3a); emit "the operator whose invertibility is in question, **and** its symbol
    where it is diagonal." §3a `:247-278` — the two-momentum DtN operator/kernel carries **both** legs
    `q_out(k), q_out(k′)`; the single-`k`/left-quantized single-leg form is the rule-17 freeze the step forbids.
  - Siblings `S11c_a_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`, `S9_export_chain_rebuild_directive.md:16-18`
    (the blindness control) as needed.
- `CLAUDE.md` rules 2, 3 (name the object not the recipe), 5 (no expected value), 6, 11, 16, 17.

⛔ Do **not** read the SymPy engine (`scripts/S11c_c1_bulk_closure_sympy_audit.py`), the step records, or any
prior-art script computing the same object — check the directive on its merits and against the code + spec above,
never against a remembered answer.

## The questions this review must answer (report a finding for any "no")

1. **The two defects are real.** From the code ranges above, independently confirm (or refute) that (R1) the
   Fredholm operator is frozen to a single leg while the operator whose invertibility is in question is the
   two-leg DtN operator, and (R2) the `DELTA_W`/`ZETA_C` labels are attached to the wrong parity combination given
   the spec `:81-82` definitions. If either defect is NOT real as described, that is a MUST-FIX finding against the
   directive.
2. **Fidelity + fix correctness (rule 3, at the construction).** Does R1 correctly demand the emitted
   `NONINVERTIBILITY_CONDITION` operator carry BOTH legs (matching `DTN_OPERATOR`), with the diagonal-symbol
   alternative permitted only via a GENUINE fresh dummy momentum (never `momentumOutput`)? Does R2 correctly bind
   `DELTA_W`↔(ζ₊−ζ₋) and `ZETA_C`↔(ζ₊+ζ₋)/2 via a computed forward/inverse binding residual, WITHOUT pre-registering
   which combination is even/odd or whether they couple/coincide? Does each fix re-enter at the construction (the
   operator's leg labelling; the per-face→parity map), never at a result?
3. **Leak / rule 5.** Does the directive state, imply, or pre-register ANY computed value, sign, parity, residual
   magnitude, or expected result? Scrutinize the construction-invariant probes (`FREDHOLM_Z_HAS_Q_INPUT/_OUTPUT`,
   the parity binding residual) and the re-freeze / label-swap controls: are these structural invariants described
   only by **what they MOVE** (legitimate), or does any of them encode an expected physics value or a "residual is
   zero for the correct X" (the specific rule-5 leak this remediation exists to remove)?
4. **Controls able to fail.** For each control, is it able to MOVE under the named corruption (re-freeze the input
   leg; swap the two parity labels), and is the binding residual a genuine two-route check rather than a tautology
   that is zero for any input? If a control cannot discriminate the fixed object from the defective one, that is a
   MUST-FIX.
5. **Scope / byte-identical bound.** Does the directive correctly EXCLUDE `NONINVERTIBILITY_CONDITION` and
   `PERMEABLE_PORT_HERMITIAN` from the byte-identical-protected core (they are what is being fixed), while
   correctly protecting the genuinely 2-leg-sound core (kernel, flat symbol, the two-leg `DTN_OPERATOR`, the
   repaired energy audit, `FACE_RESPONSE`, T-a..T-i, μ_θ)? Does it smuggle any new physics or exceed fixing these
   two defects?
6. **Completeness.** Are these two the only emit-only objects that inherited the single-leg freeze or a parity
   mislabel? In particular, does any OTHER emitted object build its own operator with both `N_0` factors on
   `momentumOutput` (like `fredholmFunctionSpaceOperator`), or attach a `δW`/`ζ_c` label to the wrong combination?
   Name any you find.
7. **The laziest passing rewrite.** Imagine the laziest edit that satisfies each fix as worded and say whether it
   leaves the emitted object correct — e.g. could the builder satisfy R1 by renaming `momentumOutput` to a symbol
   that is still not a genuine dummy, or satisfy R2 by adding a residual that is structurally zero regardless of
   the labelling?

## Method and output

This is a DOCUMENT review of the directive, cross-checked against the concrete code + spec cited above. Quote the
directive text and the code/spec text it must honor (file:line) for every claim (rule 2). State severity:
MUST-FIX (the directive would produce a wrong/weak fix, leak an answer, or misbound the byte-identical core) vs
NIT. If the directive is sound, say so and name the two or three things you checked most closely (with the
code/spec lines). Report only findings that catch a way the physics or a control could be wrong. ⛔ A prose
re-derivation is worth nothing on its own — ground every claim in a file:line you actually read.
