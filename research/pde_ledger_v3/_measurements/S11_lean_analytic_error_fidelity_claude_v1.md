# Independent fidelity review: S11 tail/Abel and inverse-stability contract (T1–T4)

## Verdict: **CLEAR**

This covers only the bounded conditional contract (T1–T4). It is not a full S11c error certificate. I found no blocking problems. The four recommendations below are not blocking; the first two are wording or count fixes, the last two are optional stronger results.

**Packet identifier (as supplied, not recomputed):** revision `S11 tail/Abel and inverse stability T1–T4 fidelity contract v1`, `aggregate_sha256 79a8697aebca1ba1c1ac29e694083daed6c1df512259d07fc59444b27e2b8a48`, 23 files.

## Limits of my own checks
- **Nothing was run.** I had only read and search tools, so I ran no Lean, Python or SymPy, and computed no hashes.
- **Author evidence only:** all build results, axiom output, the fact that each tactic succeeds (e.g. `change` in `perturbEquiv_eq`, `nlinarith` in `abel_truncation_bound`), mutant diagnostics and hash matches.
- **What I did:** read every T1–T4 source, instrument and report in the packet, and checked the mathematics by hand.
- **Mathlib from memory:** the lemma meanings I relied on (`Units.oneSub`, `setIntegral_compl`, `Integrable.bdd_mul`, `setIntegral_dirac`, `Integral.transform`) come from my knowledge of Mathlib and SymPy. The pinned Mathlib source is not in the packet.

## T1: actual integrals (`Tail.lean`, `Abel.lean`)
**Definitions**
- `tailMass` (Tail.lean:11) is the real Bochner integral of ‖g‖ over sᶜ. `firstMoment` (Tail.lean:14) is ∫|x|·‖g‖.
- Both use the complex norm; there is no real-part projection.
- `abelFactor a x = exp(-a*|x|)` (Abel.lean:8): fixed origin, absolute position.

**Integrability.**
- `truncation_error_bound` (Tail.lean:37) requires integrable g, an a.e.-strongly-measurable b with ‖b x‖ ≤ B for every x, and a measurable s. From these it derives that b·g is integrable before applying `setIntegral_compl`.
- `abel_truncation_bound` (Abel.lean:62) also requires that |x|·‖g‖ is integrable. So neither `tailMass` nor `firstMoment` can silently be Bochner's zero value for a non-integrable function.
- `omitted_integral_bound` is true but on its own certifies nothing, as the docs say.

**Combining the two errors.**
- The regulated truncation step uses c = abelFactor·b. It is measurable (continuous × a.e.-measurable), and ‖c‖ ≤ B needs a ≥ 0 through `abelFactor_bounds`.
- The pairing step uses |1 − e^{−a|x|}| ≤ a|x|, which follows from `add_one_le_exp`.
- The triangle inequality then gives B·(tailMass + a·firstMoment), the intended bound.
- The integral statements hold for any measure on ℝ. `measurable_cutoff` holds for every real R, so a negative R simply gives an empty kept set; no hidden premise is needed.
- With a = 0 the result reduces to plain truncation (`abel_zero`).
- `hB0` is redundant (it already follows from `hB`) and is unused in `abel_pairing_bound` (`_hB0`). It is consistent and harmless.

## T2: actual inverse (`Stability.lean`)
**Construction.**
- `perturbEquiv` forms `Units.oneSub(−A⁻¹E)`, i.e. the unit 1 + A⁻¹E, in the complete ring X →L X, then composes with A. The forward map is x ↦ A x + E x, which is what `perturbEquiv_eq` states.
- `CompleteSpace X` appears only where the construction needs it: `perturbEquiv`, `perturbEquiv_eq` and `exists_controlled_inverse`.
- The other estimates hold on any normed spaces over a general `NontriviallyNormedField`. No finite-dimensional, self-adjoint, real-frequency or L² premise is imported.

**Estimates.** I checked each by hand; all are correct.
- `relative_error_small`: ‖A⁻¹E‖ ≤ κε < 1.
- `coercive_bound`: (1 − κε)‖x‖ ≤ κ‖(A+E)x‖.
- `inverse_norm_bound`: ≤ κ/(1 − κε). The denominator is positive because κε < 1, and κ ≥ 0 because ‖A⁻¹‖ ≤ κ.
- `resolvent_identity`: B⁻¹ − A⁻¹ = −B⁻¹EA⁻¹.
- `inverse_difference_bound`: ≤ κ²ε/(1 − κε).
- `solution_error_identity` / `solution_error_bound`: ≤ κ/(1 − κε)·(‖δf‖ + ε‖A⁻¹f‖).
- `observation_error_bound`: a fixed map C multiplies the bound by ‖C‖.

**Existence.** `exists_controlled_inverse` combines existence with both operator bounds.
- The supplied-B lemmas (solution and observation error) apply to that B, because an equivalence's inverse is fixed by its forward map.
- Optionally, the theorem could also bundle the solution and observation bounds, or a comment could note this uniqueness.

**Scope.** The docs correctly leave the graph/outgoing realization, the common operator norm and approximate channel extraction as application obligations.

## T3: native correspondence (`S11_lean_analytic_error_source_check.py`)
**Selection.** The AST selection runs only the `xi`/`profiles` declarations, the Gaussian normalization (audit.py:1618–1623) and the half-line block (1624–1647). It neither imports the production module nor builds the closed operator.

**Hand check.**
- The Gaussian Fourier mass is √(π/a)·√(4πa) = 2π.
- With phase e^{−isξ}: the positive half-line gives 1/(a+is) and the negative half-line gives 1/(a−is).
- Their sum is 2a/(a²+s²). The even and odd parts are a/(a²+s²) and −is/(a²+s²).
- Setting s = L_W·q, multiplying by L_W and dividing by 2π gives:
  - the constant kernel P_h, with h = a/L_W;
  - the step kernel (P_h − iQ_h)/2.
- As h → 0 these tend to δ(q) and δ/2 − (i/2π)·PV(1/q), matching the spec's f̂_red convention (SHARED_PHYSICS:233–235).
- `L_W = 10` is asserted from the input file. `delta_mass` is recorded but not accepted, which is correct since its value is unclear.

**Recommendation 1 (wording, not blocking).**
- The L_W Jacobian, the substitution s → L_W·q and the division by the Fourier mass are written by the checker itself (source_check.py:60–61). They are not executed from native code.
- By reading `hat` (audit.py:1739–1751: `transform(y₃, (ℓξ, ξ))` includes the Jacobian ℓ, then divides by `fourier_mass`) and `prescribe` (1676, s = ℓ·q₃), I confirmed that the native code performs these same operations.
- However, `wrong_width_a_times_ell` checks the checker's own map, not a native operand.
- ANALYTIC_ERROR_FIDELITY.md §"Compact source connection" should say that this part of the scale/width link is checked by inspection.
- I did not trace `node.args[2]` to k_out − k_in beyond the spec's definition of Q.

**Partition.** The docs split contributions into local differential, ordinary integrable and distributional constant/step terms, all marked open or conditional (FIDELITY.md table). They claim no operator-wide estimate and no rate for the whole profile class. The assessment's Fourier/Sobolev bounds are correctly labelled deferred.

## T4: controls (contract check script and its JSON report)
**Rejections.**
- All 12 mutants have exit status 1 and exactly one diagnostic, `unsolved goals ⊢ False` in `contract_control`, with nothing else in the output. There were no timeouts or import errors.
- Note that the instrument itself only requires at least one matching diagnostic. The "exactly one" claim rests on the recorded outputs, which I inspected.

**Witnesses.**
- The tail and Abel witnesses are real Dirac integrals: the mass at 2 lies outside `Icc (-1) 1`, and 1 − e^{−1} at 1.
- The inverse witnesses are tight exact scalars:
  - inverse norm: 2 = κ/(1 − κε);
  - inverse difference: 1 = κ²ε/(1 − κε);
  - source error: 3/2;
  - observation error: 9/2.
- They are refuted as the docs describe. The `absolute_position` and `physical_width` pairs would catch a regulator without |x| or a width of a·ℓ.

**Recommendation 2 (count accuracy, not blocking).** `inverse_difference_positive` and `operator_error_omission_positive` are the same statement. The "sixteen positives" are therefore 15 distinct statements.

**Recommendation 3 (optional, stronger).**
- The scalar controls do not instantiate `exists_controlled_inverse` or `abel_truncation_bound`, and the positives never satisfy all of either theorem's premises at once.
- The premises are plainly satisfiable (e.g. g = 0 or A = id), and the limitation is accurately stated in FIDELITY.md:97–103, VERIFICATION.txt:64–66 and the validation report's `limitations`.
- An ℝ instance with A = id and E = −½·id would tie the controls directly to these theorems.

**Recommendation 4 (optional).** `real_phase_norm` is unused. The docs could say explicitly that the Fourier phase is absorbed into b, so B is uniform in real transfer only; complex transfer is not covered, as the assessment already notes.

**Audit root.** It contains exactly 41 `#print axioms` lines, and the recorded output lists only `propext`, `Classical.choice` and `Quot.sound`. I found no `sorry`, `admit` or `axiom` in the four modules.

## Not reviewed (outside T1–T4)
The nonlinear-pole work, variable coefficients, limiting absorption, and uniform physical moments or operator estimates are open application obligations. None of them blocks this conditional contract.
