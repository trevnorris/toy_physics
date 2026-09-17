# D4 odd density: supplied object and evidence boundary

Fidelity record for [D4_ODD_COVERAGE.md](D4_ODD_COVERAGE.md).
**D4B.1–D4B.4 is complete.** Local verification and current correspondence
pass, and both independent fidelity reviews are CLEAR. Review provenance,
optional findings and their dispositions are recorded in
[D4_ODD_FIDELITY_REVIEW.md](D4_ODD_FIDELITY_REVIEW.md).

The existing D4 classification at `d6119f75` fixes the density
`P=(G01-G10)(G23-G32)-(G02-G20)(G13-G31)+(G03-G30)(G12-G21)`.
Native Q9's actual emitted `P_D` equals P; the fully summed epsilon expression
is `2P`. The current increment uses `L_beta=-(beta/2)P`, matching the odd term
in the supplied stiffness action. Beta is any constant real coefficient.

`S11D4Odd.Action` reuses `odd_classification`, links the density to
`invariantForm ![0,0,0,beta]`, and defines momenta by actual differentiation of
L. With `G_ij=partial_i u_j`, those momenta are proved to be `-(beta/2)M_ij` for the
explicit complementary matrix in the coverage contract. The time momentum
is zero. No field equation, wave ansatz, positivity or nonzero-beta premise
enters this definition.

`S11D4Odd.Calculus` supplies the finite-sum/product and commuting-partial rules
on S10's existing real spacetime. `Boundary` uses smoothness to establish
`sum_i partial_i M_ij=0`. The current `K_i=(1/2)sum_j u_j M_ij` then satisfies
`div K=P`. The factor follows from `sum_ij G_ij M_ij=2P` and is not arbitrary.

`Variation` differentiates the actual finite relative-action integral with
smooth compact variations, reusing S10's integration-by-parts machinery. A
smooth background need not have integrable total density. The proved result
is zero local Euler–Lagrange expression and zero first variation for every
such background and variation. This is a bulk statement: a divergence may
contribute when a physical boundary is present. The nonzero density and
momentum witnesses distinguish that statement from pointwise vanishing.

The compact native instrument extracts only the selected original Q9 and
coordinate/EL helpers. It checks the actual computed P, all momenta, the
explicit current and the odd combination of native V5. Native EL has the
opposite overall sign to Lean's variational expression; both vanish here.
Consequently sign/factor controls use nonzero momenta. An added `G00²` term
is a deliberate non-null control, also used to catch a native helper mutated
to return zero everywhere. These are exact Python/SymPy translation checks,
not kernel-certified CAS code. The native V5 comparison concerns the actual
odd combination in its emitted V1 basis, before multiplying by `-beta/2`;
linearity preserves its zero result. Wolfram evidence is source inspection only.

The supplied specification is `directives/S11_SHARED_PHYSICS.md`: §Q1 fixes
the native EL sign, §7 gives `W_XFORM_EXTRA` the additive term `(beta/2)P_D`,
and the §7 `P_D` rule uses the actual computed V6 basis without rescaling.
With `L=T-W`, the isolated odd term is exactly `-(beta/2)P_D`. The current
native check recomputes this actual P through selected original helpers;
it does not substitute a hand-chosen P into the production script.
Lean numbers time as coordinate 0, followed by the four spatial coordinates.
The native functions use argument order `(x1,x2,x3,x4,t)`; the correspondence
uses these named spatial derivatives and t, not the position of t in that list.

The recorded suite built fourteen unchanged local imports in dependency
order in run1, then reused those guarded objects in run2. Run2 freshly built
the four new modules and their 41-declaration audit root. Twelve mathematical
mutations were rejected and twelve positives passed, covering density/momentum/current
normalization, the derivative-index cancellation, the bulk conclusion and
admissible zero/negative coefficients. A timeout or instrument failure is not
accepted as mathematical rejection. Old proof and evidence files remain
historical if rebuilding unchanged imports changes object hashes.

The first recorded run passed the compact native check and all fourteen
unchanged imports. It stopped in `Action` on opaque matrix-entry reduction
and an unreachable trailing tactic; no mutation ran. The repairs explicitly
reduce the finite entries before using smoothness or polynomial tactics and
use the current finite-sum calculus lemma name. They do not change the
definitions or theorem statements. The failed run and focused diagnostics
are retained under `_scratch/S11_lean_d4_odd/`. All fourteen dependency
objects stayed byte-identical on this rebuild. Run2 rebuilt the new modules
and ran every control; focused checks were not accepted as recorded reuse.

Before run2, all four new modules passed focused builds with strict warnings.
The explicit divergence and first-variation statements compiled. Repairs also
preserve the finite sum used for the EL cancellation and explicitly reduce
the concrete density, momentum and current witnesses. The admissible affine
positive control uses the same explicit finite-entry reduction. The completed
run2 supplies the recorded build, axiom audit and all twelve mathematical
rejections and twelve positives. Current source, dependency, instrument,
generator, input/output object and log hashes have been validated, as have
the fifteen clean tracked package checkouts at their pinned commits.

`D4_ODD_VERIFICATION.txt` records each control group and distinguishes the
intended false mathematical identities from secondary diagnostics. The
mutation replacing every divergence derivative by `partial_1` is false for
`u=(0,0,0,x1*x2)`; doubling the current is false for `u=(0,x1,0,x3)`.
The mathematical failures occur in the required named declarations. Later
rewrite, simp or concrete-conversion errors are not accepted as evidence.

This verifies the existing constant-coefficient D4 odd cancellation and its
normalization exhaustively over smooth fields. It is not a new physical
discovery, a proof of the complete `XFORM_EXTRA` spectrum, or a claim that the
odd density vanishes pointwise. In particular, `nonzero_density` gives P=1,
and `nonzero_momentum` gives a momentum of -1 at beta=2.

The domain and exclusions in the coverage contract are binding. No full D4
even-family bulk classification, variable beta, D5, spectra, interfaces,
S11c work, production rerun or systematic CAS bridge is included. The user
approved this fixed 38-file packet; both independent reviewers cleared it
without required corrections. Closure changes only documentation and review
state: all proof, instrument, native-source and object bytes retain their
reviewed correspondence. The approved snapshot and historical evidence remain
unchanged. See `S11_lean_d4_odd_closure_validation.json` for the explicit
post-review documentation changes. No substantive re-review or rebuild is
needed, and no further proof expansion belongs to this completed contract.
