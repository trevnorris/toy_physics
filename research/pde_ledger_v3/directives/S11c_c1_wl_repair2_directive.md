# S11c-c1 Wolfram engine — repair-2 directive

## Role and single deliverable

Modify `research/pde_ledger_v3/mathematica/S11c_c1_bulk_closure_mathematica_audit.wl` in place (baseline
`13f0bd2c`). Its only product is the flushed stdout tag stream (regenerated). ⛔ Import nothing (the blindness
control stands — no `Get`/`Import`/`<<`/`ReadString`/`OpenRead`, no absolute repo path, no SymPy engine/export/
`.out`); the copy-to-empty-dir + structural-scan acceptance from the build directive still holds. This repair
fixes **two propagated emit-only defects** (R1, R2 below), both confirmed against baseline `13f0bd2c` by the two
decision legs on this directive. ⛔ Do **not** touch the sound core (§ "byte-identical").

`CLAUDE.md` binds. The sole physics authority is `directives/S11c_c1_SHARED_PHYSICS.md` (§§0–8); the build
directive `directives/S11c_c1_wl_build_directive.md` governs naming/blindness/idioms unchanged. Add no expected
value or acceptance criterion (rule 5): each fix names the **object** the emit must be and the **construction
invariant** its control must satisfy — ⛔ never a computed value, sign, residual magnitude, block-arithmetic
formula, or the phrase "is zero for the correct …". A control is described by **what it MOVES**, never by what a
residual equals.

## Why this directive exists (context, not a task)

The prior repair directive `S11c_c1_wl_repair_directive.md` was built against without its rule-7 decision legs.
Run retroactively, the two legs flagged two propagated emit-only defects. This directive's own two decision legs
(Codex + Grok, both in) then **corrected the second finding's location**: the retroactive legs had reported a
"parity swap" in `PERMEABLE_PORT_HERMITIAN`, but that object is **correct** (see the note under R2); the real
propagated parity defect is a **dead parity axis** elsewhere. This directive fixes exactly R1 and the corrected
R2, both emit-only (neither is exported; c2 consumes only the DtN operator / two-momentum kernel / face response,
all separately 2-leg-verified). The full record is `_measurements/S11c_c1_wl_remediation_plan.md`.

## The two defects — SUPPLIED, code-verified against baseline `13f0bd2c`; each fix re-enters at the CONSTRUCTION

### R1 — the Fredholm / `NONINVERTIBILITY_CONDITION` operator must be the two-leg DtN operator `Z` (spec §3b `SHARED_PHYSICS.md:299-306`, §3a `:247-258`; rule 17)

**Defect** (`fredholmFunctionSpaceOperator`, `.wl:580-597`; consumed by `NONINVERTIBILITY_CONDITION` at
`.wl:1537-1549`): this function builds its own operator with `nZero = FourierMultiplier[1/(I qOut[omega,
momentumOutput])]` and `gZero = FourierMultiplier[I qOut[omega, momentumOutput]]`, then composes
`OperatorComposition[gZero, multiplication, gZero]` and `OperatorComposition[nZero, dtnVariation, nZero]` — every
`N_0`/`N_0^{-1}` factor carries `momentumOutput`; `momentumInput` (and any genuine dummy) appears nowhere, and the
`multiplication` even drops the `FourierKernelData` profile source the real operator carries (`.wl:585-586` vs
`.wl:556-558`). This is the identical single-leg WKB/left-quantized freeze spec §3a forbids and that R1 of the
prior repair removed from `operatorCompositionRecordFromDerivation` (`.wl:542-572`, now correctly two-leg via
`gZeroOutput`/`gZeroInput`, `nZeroOutput`/`nZeroInput`, `rightLegMomentum` defaulting to `momentumInput`) — but
`fredholmFunctionSpaceOperator` is a **separate** function that was never updated and was explicitly listed
"byte-identical, do not touch," so the freeze persists in the operator whose invertibility `NONINVERTIBILITY_
CONDITION` reports. The `.wl:577-579` comment claims the Fourier momentum is "a dummy coordinate, rather than
either external kernel leg," but the code uses `momentumOutput`, which **is** the external output leg (used as
such everywhere else in the file) — so the comment asserts a property the code does not have.

**The object the fix must produce:** the operator inside `NONINVERTIBILITY_CONDITION` whose invertibility is in
question is `[I + (Λ_A/ρ_m²) Z]` with `Z` the **two-momentum DtN operator** — the SAME operator, carrying **both**
legs `q_out(momentumOutput)` and `q_out(momentumInput)` and the profile source `Ŵ_bg(momentumOutput−
momentumInput)`, as `DTN_OPERATOR`/`operatorCompositionRecordFromDerivation` (the rightmost `N_0`/`N_0^{-1}` of
each `OperatorComposition` carrying the input leg, the leftmost the output leg). The cleanest realization is for
`fredholmFunctionSpaceOperator` to **be** that two-leg `Z` (e.g. reuse `operatorCompositionFromDerivation
[anchoring]`), so `NONINVERTIBILITY_CONDITION` wraps the genuine `[I + (Λ_A/ρ_m²) Z]`. Spec §3b asks for "the
operator whose invertibility is in question, **and** its symbol where it is diagonal": the diagonal symbol is the
**already-emitted flat Fourier-diagonal symbol** `FLAT_DIAGONAL_SYMBOL_RELATION`/`DEGENERATE_LOCI_*` — keep those
**unchanged**. ⛔ Do **NOT** add a "diagonal-symbol reduction of the full operator at `k=k′`": the nonuniform
operator is not diagonal, `k=k′` merely selects a diagonal kernel slice, and a fresh dummy there only **renames
the freeze** — there is no legitimate second single-momentum object to emit, so tag count stays **51**. Make the
`.wl:577-579` comment true (the operator carries both external legs) or remove it. **Name the object; ⛔ do not
re-protect it byte-identical.**

**Construction invariant** (emit computed probe objects, ⛔ not a bare check, ⛔ not a leaked-zero residual):
over the emitted `NONINVERTIBILITY_CONDITION` operator (its `Z`-part), a probe returns `FREDHOLM_Z_HAS_Q_INPUT`
**and** `FREDHOLM_Z_HAS_Q_OUTPUT` — both computed from the emitted operator by the same containment instrument the
file already uses (`.wl:619-626`), ⛔ not asserted. The **re-freeze control must be constructor-level**: re-invoke
the operator constructor with the input-leg momentum set to the output leg (`rightLegMomentum → momentumOutput`,
exactly as `DTN_OPERATOR` re-freezes at `.wl:1329-1332`), computed on the **same** `Z` the operator field emits —
⛔ NOT a post-hoc `momentumInput → momentumOutput` substitution over the whole expression (a syntactic `Count`
probe can be gamed by sprinkling `qOut[omega, momentumInput]` without a genuine two-leg `N_0`; the control must
change the construction, not the surface text). That constructor-level re-freeze **MOVES** `FREDHOLM_Z_HAS_Q_INPUT`
(both-legs-present → input-absent). ⛔ Do **not** emit a "`operator − DTN_OPERATOR` is zero" residual and ⛔ do not
state any residual value — that is the rule-5 leak the prior R1 carried. Describe the control only by the probe it
moves.

### R2 — remove the DEAD parity axis in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`: emit the parity-COMBINATION form (spec §3b `:308-320`)

**Defect** (`PERMEABLE_DISSIPATION_VS_OMEGA_TAU`, emission `.wl:1713-1736`): the emit iterates over both `face ∈
{+1,−1}` **and** `parity ∈ {DELTA_W, ZETA_C}`, but the payload is `"HERMITIAN_FORM" -> portCase["HERMITIAN_FORM"]`
with `portCase = portCaseForAxes[anchoring, face, density, outputRegime, inputRegime]` (`.wl:1717`) — a
**per-face** object that does **not** depend on `parity`. The `parity` iterator affects only the key
(`portCaseKey[…, parity]`); `INDEPENDENT_MEMORY_LIMITS` likewise comes from `memoryLimitPackage[anchoring, face,
…]` (`.wl:1081-1090`), per single face. So `PARITY_DELTA_W` and `PARITY_ZETA_C` receive the **identical
untransformed per-face** memory Hermitian form — a fake parity axis. Spec §3b `:308-320` requires the dissipation
objects "per regime pair and **parity combination**," alongside `PERMEABLE_PORT_HERMITIAN`.

**The object the fix must produce:** for each parity combination (`δW`, `ζ_c`), the memory/dissipation Hermitian
form is the **parity-combination of the two per-face** memory Hermitian forms, formed by the **same change of
basis** `PERMEABLE_PORT_HERMITIAN` already uses (`portParityCombination` `.wl:1608-1650`, `faceToParityMatrix`
`.wl:696-698`) — ⛔ **not** the per-face form duplicated under both keys — and its ωτ_I `ZERO_LIMIT`/`INFINITE_LIMIT`
memory limits computed on **that combined form**. Emit under `PARITY_DELTA_W` the `δW`-combination and under
`PARITY_ZETA_C` the `ζ_c`-combination. ⛔ Do **not** pre-register which combination is even/odd, whether they
couple, or that they differ — those are outputs; ⛔ do not hand-write any block-arithmetic combination of the
matrices (`(plus−minus)` etc.). Reuse the existing, already-protected change-of-basis construction.

**Construction invariant:** the `PARITY_DELTA_W`-key and `PARITY_ZETA_C`-key payloads are **distinct computed
objects** wherever the two per-face memory forms differ (the curved case) — so a **one-sided corruption of only
the `+`-face** memory Hermitian form (at its construction) **MOVES** the two parity-combination payloads
**differently**, and the two keys' payloads are ⛔ not byte-identical in the curved case. ⛔ Do not state any
residual value or say a combination "is zero"; describe the control only by the corruption it moves.

> **NIT (same class, lighter — a redundant axis, not wrong values):** `UNIFORM_LIMIT_S11CC1_OPERAND` /
> `UNIFORM_LIMIT_S11B_OPERAND` (`.wl:2017-2039`) carry the same `PARITY_DELTA_W`/`PARITY_ZETA_C` keys on a **flat**
> (parity-independent, equal-face) object, so both keys duplicate one value. Emit the flat object as
> parity-independent, or emit the parity-coincidence as a **computed** statement — ⛔ not a silently duplicated
> dead axis — **while keeping `UNIFORM_LIMIT_RESIDUAL` a valid C1-vs-S11b comparison** (`.wl:2043-2044`). If the
> cleanest fix would alter that residual's operands, leave the object and flag it for the step record instead.

**⭐ Note — why `PERMEABLE_PORT_HERMITIAN` is NOT touched (it is correct; the retroactive finding mislocated the
defect).** `portParityCombination` (`.wl:1623-1636`) sets `deltaWHermitian = (plusHermitian+minusHermitian)/4`,
`zetaCHermitian = plusHermitian+minusHermitian`, coupling `(plus−minus)/2`. These are the **congruence** blocks
`Aᵀ·diag(P₊,P₋)·A` of the power-conjugate form under `A = faceToParityMatrix = [[1/2,1],[1/2,-1]]` — correct given
the outward orientation `V_s = s∂_tζ_s` (`S11c_a_SHARED_PHYSICS.md:58`): both diagonal blocks are even
combinations, the odd combination is the coupling. The naive reading (that `DELTA_W` should be `(ζ₊−ζ₋)`) is a
category error (port **matrices** are not face **values**) and would make the thickness port vanish when the faces
share a map. `PERMEABLE_PORT_HERMITIAN` therefore stays byte-identical. (A separate control-strength observation —
that these blocks are hand-written to match the congruence rather than computed by applying the stored
`portParityTransformation` — is deferred to the step record / T7 cross-check, ⛔ not folded here.)

## What must stay byte-identical (⛔ do not touch — the 2-leg-sound core)
`DTN_KERNEL`, `DTN_FLAT_SYMBOL`, `DTN_OPERATOR` (the two-leg `operatorCompositionRecordFromDerivation`, already
repaired), `DTN_RIGID_SHIFT_*`, `DTN_BY_REGIME_PAIR`, `DTN_BY_PARITY`, `DTN_HERMITIAN_PART`/`DTN_REACTIVE_PART`,
`DTN_INERTIAL_LOADING`, `DTN_GRAZING_BEHAVIOUR`, the permeable `FACE_RESPONSE`/`FACE_RESPONSE_COEFFS` operator
inverse and its `t_s` (and its similarity-transform `PARITY|…` view, `.wl:1507-1528`), **`PERMEABLE_PORT_HERMITIAN`**
(correct — see the note under R2), the repaired energy audit (`ENERGY_*`: far-field Poynting bulk operand +
response-`t_s` face operand), the flat diagonal-symbol loci (`FLAT_DIAGONAL_SYMBOL_RELATION`, `DEGENERATE_LOCI_*`),
the §5b/§5d/§5e controls, the T-a..T-i re-derivation, the reserved-name spellings, and `μ_θ` opacity. ⛔
**`NONINVERTIBILITY_CONDITION`/`fredholmFunctionSpaceOperator` and `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` are NOT in
this list** — they are what R1 and R2 fix.

## Method and script clauses (spec §6 — unchanged, binding)
Every fix re-enters at the construction (the Fredholm operator's leg labelling; the per-face→parity combination of
the memory form), ⛔ never at a result. Each control emits its operands and the computed probe/residual **before**
any guard; a physics disagreement emits and continues (exit 0); nonzero exit is operational only. ⛔ No `VERDICT`/
`PASS`/`FAIL`; a boolean is a typed CAS object (unevaluated relational / `STATUS_TOKEN`), ⛔ never a native `True`/
`False` residual operand. ⛔ No hand-typed CAS object with no data dependence — every probe/combination is
**computed** from the emitted operator/per-face forms via the existing constructions, ⛔ not asserted and ⛔ not
hand-written block arithmetic. ⛔ No leaked residual value and no "is zero for the correct labelling" (rule 5; the
specific defect being remediated). Run discipline unchanged (detached; one kernel; `danger-full-access`;
600s-no-output/RSS kill recorded-not-narrowed; per-case streaming). Tag count unchanged (**51**).

## Builder report (under 25 lines)
The diff vs `13f0bd2c` (functions changed: expect `fredholmFunctionSpaceOperator` and the
`PERMEABLE_DISSIPATION_VS_OMEGA_TAU` emission plus their probe emits; and the `UNIFORM_LIMIT` NIT if taken);
confirmation that the byte-identical core rows above — including `PERMEABLE_PORT_HERMITIAN`, `DTN_KERNEL`, the flat
symbol/loci, `FACE_RESPONSE` — are unchanged; the two repaired controls' computed probe objects
(`FREDHOLM_Z_HAS_Q_INPUT`/`_OUTPUT` and their constructor-level re-freeze move; the two parity-combination payloads
being distinct and their one-sided per-face-corruption move); tag count (51, or state what changed and why); tasks
run; runtime; any ambiguity. ⛔ No computed physics value or conclusion.
