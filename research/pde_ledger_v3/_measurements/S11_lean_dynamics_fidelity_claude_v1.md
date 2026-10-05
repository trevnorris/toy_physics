# Independent fidelity review: S11 D2 odd-invariant dynamics, E1–E4

**Packet:** `MANIFEST.json`, revision "S11 D2 odd-invariant dynamics E1–E4 fidelity contract v1"
**Aggregate SHA256 (as stated in the manifest):** `16a1dd5784c344f032716d86974b52a58c79110464e2cb1353b915c17d2581c0`

**Verdict: CLEAR.** No blocker remains for this bounded contract.

**Limits of this review:**
- **Nothing recomputed or rebuilt.** This session had no shell, so I did not recompute any file hash or the aggregate. I did not rebuild Lean or rerun either instrument.
- **Hash evidence is string comparison only.** Every source, instrument and native hash in the three records matches the manifest text. The four rebuilt S10 `.olean` hashes and the three new `S11OddDynamics` `.olean` hashes match between `DYNAMICS_VERIFICATION.txt`, `_validation.json` and `_contract_checks.json`.
- **Scope.** The review is by reading and by hand calculation.

## E1 — Density and object identity
- **Gradient indexing.** `spatialGradient J j i = J j.succ i` gives `G_ij = ∂_i u_j`, with Lean coordinates (t, x1, x2).
- **P matches `oddPairing`.** `S11Invariants.oddPairing = invariantForm ![0,0,0,1]` evaluates to `coords0·coords1 = (G00+G11)(G01−G10)`. That is `divergence J · skew J` with `d = ∂1u1+∂2u2` and `c = ∂1u2−∂2u1`.
- **`density_identity`.** It states `lagrangian β J = −β/2 · oddPairing`. `lagrangian` is defined directly as `−β/2·d·c`, and β is an arbitrary real constant with no premises.
- **Native P_D.** It comes from the V6 basis in `compute_q9(2)`. The sum of the RREF rows has leading monomial `g1g2` with coefficient +1. The recorded value `G11G12 − G11G21 + G12G22 − G21G22` equals `(G11+G22)(G12−G21)`, and I checked it by hand.
- **Package difference.** The spec (`S11_SHARED_PHYSICS.md` §7: `L = T − W`, `W_XFORM_EXTRA = … + (β/2)P_D`) and `package_build` (lines 449–474) both give XFORM_EXTRA − MAIN = `−β/2·P_D`. The kinetic, curl and div terms cancel exactly, so the instrument's replacement symbols for ρ, μ_R and B_comp do not matter.
- **Result.** The recorded `actual_delta` has the right sign and the half factor.

## E2 — First variation
- **Momenta are real derivatives.** `momentum` is `deriv` of the density along `basisJet`. `lagrangian_variation` gives `HasDerivAt`, so `deriv` is not a default value.
- **`momentum_eq` checked by hand.** Row t is (0, 0). Row x1 is (−βc/2, −βd/2). Row x2 is (βd/2, −βc/2).
- **`eulerLagrange_eq` checked by hand.** `E = −Σ_j ∂_j p_{j·}` gives `β/2·(∂1c − ∂2d, ∂1d + ∂2c)`. No mixed derivatives are commuted.
- **Relative action.**
  - `relativeAction` integrates the pointwise density change, and `u` is only assumed `SmoothField`.
  - `relative_density_integrable` works through the exact quadratic expansion `lagrangian_change`. The linear term is integrable because each term is a continuous momentum times a compactly supported `∂h`. The s² term is integrable because it is half of `linearDensity h h`.
  - The derivative at 0 comes from the proved polynomial expansion in s; no dominated convergence is needed.
  - Integration by parts uses S10's `coord_integration_by_parts` with a smooth momentum and a compact test field.
- **All compact tests ⇔ pointwise EL.** `actionStationary_iff_eulerLagrange` uses Mathlib's `ae_eq_zero_of_integral_contDiff_smul_eq_zero` with `single_testField`. `Measure.eq_of_ae_eq` then uses continuity and the open-positive Lebesgue measure. The converse direction is immediate.
- **Witness.** `witnessField = planeWave 0 ![1,0] ![1,0] = (cos x1, 0)`. By hand, d = −sin x1 and c = 0, so E(0) = (0, −β/2), which matches `witness_eulerLagrange`.
- **Nonzero first variation and nullness.**
  - `exists_nonzero_firstVariation` gives, for every β ≠ 0, an admissible compact test with nonzero first variation. The proof is classical and non-constructive, which is acceptable for an existence claim.
  - `variationallyNull_iff` quantifies over all smooth fields and gives nullness exactly when β = 0.
  - Nothing claims that every background is non-stationary, or classifies null Lagrangians in general.

## E3 — Modal operator and mixing
- **Amplitude derivative.** `modalAction_eq` gives `−β/2 (k·a)(r·a)` with r = (−k2, k1), and `modal_variation` shows that `M a` is its gradient.
- **Independent PDE route.** `eulerLagrange_planeWave` holds for every real ω, k, a and point x. It goes through `fieldJet_planeWave` (the jet is −sin φ · modeJet(−ω)), `momentum_smul`, `partial_const_sin` and `momentum_contraction`. I checked the contraction by hand for both components. The phase is `waveCovector = (−ω, k)`, i.e. k·x − ωt.
- **Polarization.**
  - `polarization_decomposition` covers every amplitude when k ≠ 0; `dot_turn` and `normSq_turn` supply orthogonality and equal norms.
  - The directional actions are `M k = −β|k|²/2·r` and `M r = −β|k|²/2·k`, and both cross elements equal `−β|k|⁴/2`.
- **Mixing.**
  - `mixing_iff` shows mixing exactly when β ≠ 0 and k ≠ 0.
  - `mixing_cases` covers β = 0, then (β ≠ 0, k = 0), then (β ≠ 0, k ≠ 0). The three cases exclude each other by their guards.
  - Nothing is claimed about the full XFORM_EXTRA spectrum.

## E4 — Native link and controls
**Compact source instrument.**
- **What it runs.** It executes only the selected native helpers, taken from the source by AST.
- **What it compares.**
  - P_D is taken from `compute_q9(2)`.
  - The action difference is the real `package_build` difference.
  - The odd V5 combination is selected by its RREF pivot coordinates [0, 1, 0, 0], with the reconstruction asserted.
- **Hard-coded formula.** The only hard-coded formula is the Lean-side comparison target; nothing is inserted into the engine.
- **Recorded objects, checked by hand:**
  - native EL = −(Lean EL), because the native code uses the positive-divergence convention (as in the spec's Q1 and in Wolfram line 2183);
  - native EL = −β/2 · V5_odd;
  - route A = −M, with the cos(phase) factor stripped;
  - route B = M/2;
  - the period average is modalAction/2;
  - at β = −2, k = (3, 4), the cross value is 625.
- **Production wiring.** The production driver (lines 2463–2477, 2807) routes `compute_q9` output into `package_build`, the EL and both route helpers. Those helpers are linear in the action, so checking the difference is enough.
- **Limit.** This is a tested SymPy translation, not a kernel-certified link. The Wolfram evidence is limited to the four source anchors; I confirmed that each exists and fits this reading.

**Rejected controls (10).**
- **Source mutations.** I inspected the primary failure of each; all three are genuine coefficient mismatches in the named declaration:
  - `density_sign`: `density_identity` is left with +β/2 against −β/2 on every monomial.
  - `density_factor`: the coefficients are ±β against ±β/2.
  - `local_EL_sign`: both component goals of `eulerLagrange_eq` have opposite coefficients.
  - The secondary failures (`change`, `rewrite`) are not counted.
- **Statement mutants.** All seven end in `⊢ False`.
- **Environment.** There are no import, syntax or timeout failures (the `invalid` regex guards against them), and exit status is 1 in each case.

**Passing controls (8).** All eight are admissible:
- (1, ![1,0]) mixes;
- (−2, ![3,4]) mixes;
- zero β and zero k do not mix;
- the null / non-null loci;
- the witness EL value is −1/2;
- a compact test with nonzero first variation exists at β = 1.

**Audit.** The axiom audit covers all 47 theorems in the three modules, and each uses only `propext`, `Classical.choice` and `Quot.sound`. There is no `sorry`, `admit` or custom axiom.

## Blockers
None.

## Optional suggestions (not conditions of closure)
1. **Lean-side modal control.** No Lean-side mutation targets `modalOperator` itself. The modal sign and factor are tested by the source instrument (`missing_action_half` fails the route-A identity) and indirectly by `witness_EL_sign`. A one-line sign mutation of `modalOperator` would make this control explicit.
2. **Adjudication wording.** For `density_sign` and `density_factor`, `DYNAMICS_VERIFICATION.txt` calls the secondary failures "change-tactic" failures. There is also a mathematical unsolved goal in `lagrangian_increment`, so the wording could be more precise.
3. **"Instantiates" overstates the reuse.** `DYNAMICS_FIDELITY.md` says the S10 lemmas are "instantiated". The Analytic lemmas are reused, but the Variation-level chain (about 13 lemmas) is copied nearly verbatim for the odd density. A density-generic S10 variation theorem would fit L1 better; this does not affect correctness.
4. **Instrument locus checks.** The five locus checks evaluate the Python copy of the Lean M, not the native routes. They are covered by the symbolic native identities anyway, and it would help to label them as such.
5. **Object hashes are tied to this run.** Dependency objects were written with `lake env lean -o`, not `lake build`. The four changed S10 `.olean` hashes are therefore bound to this run's invocation, as the records already state, and a later `lake build` may produce different objects again.

**Final verdict: CLEAR** for the bounded D2 odd-invariant dynamics contract E1–E4. This covers one of the two independent review legs, with the limits stated at the top: no hash recomputation and no independent rebuild.
