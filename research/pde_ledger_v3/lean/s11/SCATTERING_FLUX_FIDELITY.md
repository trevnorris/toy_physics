# F1–F4 fidelity record — bounded contract complete

The object is a supplied finite complex matrix J acting on full amplitude
vectors, with `pair J x y = star x dot (J mulVec y)`. There is no diagonal,
positive, transverse-only or nondegenerate restriction. The physical application
must supply all eligible channel directions and any needed closed-mode matching
data. A matrix type cannot certify that a numerical basis spans the physical
mode space. Empty finite types and indefinite forms remain in the theorem domain.

`flux J x` takes the real part of this pairing. `hermitian_flux_real` proves
the imaginary part vanishes when `J.conjTranspose = J`. Arbitrary non-Hermitian
physical-current residuals are not thereby declared acceptable: the nonzero
imaginary witness deliberately distinguishes the raw pairing from its real part.
The specification permits a possibly non-Hermitian mode pencil; that is not a
proof that any supplied energy-current matrix is Hermitian or correctly real.

`flux_add` retains both cross contractions for every matrix. Their real sum must
vanish for flux additivity (`flux_add_iff_cross_zero`). This applies to a proposed
sector split or incoming/outgoing boundary split only after the actual amplitude
decomposition is supplied. There is no automatic conservation or current
orthogonality assertion. In particular, closed matching modes cannot simply be
deleted from a finite-boundary interference identity.

`pullback J C = C^H J C` and `pair_pullback` identify the whole form under a
coordinate map. Rectangular C is allowed, but selecting fewer coordinates is
not an invertible change of the complete mode basis. The scattering-coordinate
identity supplies `Cout * Dout = 1` and transforms S to `Dout * S * Cin`; it
preserves the flux for the represented input. It does not assert that an
arbitrary S solves the PDE. The current changes together with coordinates;
amplitudes alone are not flux-normalized observables.

The normal coordinate convention follows shared physics section 3a: outward
left/right signs are -1/+1, while the incoming current has the opposite sign.
The two-end outgoing value is `-left + right`, not a sum of absolute values.
`fraction` returns `none` at exactly zero denominator, including zero flux at a
nonzero vector of an indefinite form. Negative denominators are included
algebraically but not relabelled as positive incident physical channels.
Nonnegative fractions require nonnegative numerator and positive denominator;
the upper bound one requires numerator no larger than denominator.

`conditional_balance` starts with the explicit physical application premise
`converted + survived + defect = incoming`, with nonzero incident flux. It gives
the normalized sum with the defect retained. The converse identifies zero
defect as necessary for that sum to be one. This supplies no conservation law,
no assumption S^H S=I, and no missing bound-capture observable.

## Compact native identification

The instrument extracts only the original `multiply`/`quadratic` functions from
`S11c_d_continuum_currents.py`, `adjoint`/`gram` from
`S11c_d_continuum_response.py`, and `RectangularModeJets.multiply` from the
existing engine. They execute on small integer/Gaussian-integer fixtures at
grade (0,0). Their sources are not imported as modules; no saved scientific
packet, physical constructor or production driver is executed. The tested
operations are the actual complete conjugate contraction and congruence. Source
inspection separately anchors the outgoing s and incoming -s factors in
`open_metrics`; that physical map is not executed by this check.

The fixtures test complex off-diagonal interference, basis covariance, zero
amplitudes, indefinite null amplitudes and an empty amplitude space. Wrong
transpose, discarded off-diagonals, transposed current and stale basis metric
formulas are compared. A separate isolated mutation removes `.conj()` from the
original native adjoint function. These are tested translations outside the
Lean kernel, not certified numerical scattering results. The selected NumPy
fixture values and operations are exact in binary arithmetic; no approximate
equality threshold is used here. The report pins full source and selected AST
hashes. Python/NumPy runtime versions were not captured by this native instrument;
this record does not claim a separately reproduced environment or a certified
numerical scattering computation.

The Lean `fraction` has no executed native counterpart in this increment.
The native `quotient` uses a complex diagonal denominator without a zero guard;
its physical reality and denominator obligations belong to the subsequent
bookkeeping assessment. No numerical defect is inferred here. The selected
engine `RectangularModeJets.multiply` is tested only at grade (0,0); its
componentwise grade cutoff is not the unrestricted current-series product.

## Verification and review boundary

Recorded run8 freshly built Current, Balance, Controls and the audit root.
All 34 selected axiom lists contain only `propext`, `Classical.choice` and/or
`Quot.sound`. Eleven paired false statements each reached exactly one unsolved
`False` in `contract_control`, without warnings or other errors. Their eleven
true counterparts and four further positives passed. These are fresh statements
decided with canonical witnesses (or exact arithmetic), not independent
rederivations. The conservation mutant is arithmetic; separate admissible
balance positives exercise the theorem and its nonvacuous premises. The compact native check
passed twelve checks (seven identities/examples and five wrong-formula controls).
There are 31 records: one native check, four builds and 26 formal controls.

Development repairs preserved the intended definitions and claims. Initially
external compiled caches were absent; the exact pinned cache was restored from
2989 local archives with no download or historical ledger rebuild. A native
stale-metric fixture was refused because its basis map fixed the chosen vector;
the corrected map changes that vector. Broad tactic imports hit the internal
allocator budget; specific imports plus a 2048 MiB Lean budget passed inside
the unchanged mandatory 2 GiB whole-job/no-swap/one-CPU guard. Run8's resource
receipt records the observed cgroup peak and zero max/OOM/swap events.
Explicit pointwise star-add rewriting,
unused binder removal and the concrete matrix/finite-sum imports completed the
canonical builds. Run1's supervisor correctly recorded ERROR but did not return
a nonzero process status; subsequent supervisors propagate child failures.

Runs6–7 refused the zero-denominator mutant because its goal remained
`none = some 0`, not the required literal False. The pinned NormNum implementation
omits constructor simprocs, so run8 supplies the proved Option constructor
rewrite facts explicitly; all controls then passed their intended checks. No
canonical proof or control statement changed in that repair. No
syntax/import/resource failure or unreduced control diagnostic
is accepted as a rejected mathematical statement. Inspect both the scientific
report and the separate containment receipt; neither a process exit alone nor
an empty log is a proof of success.

All new proof objects are isolated. No previous proof object or historical
report is replaced. `S11_lean_flux_preserved_inputs.json` records 528 historical
files and 232 shared objects. The new record was taken after the first failed
native attempt; it is a preservation baseline, not evidence of that failed job.

`S11_lean_flux_validation.json` rechecks all recorded logs and controls, live and
snapshot source bytes, four output objects, fifteen clean package pins and seven
direct Mathlib source/object pairs. The external transitive cache remains the
pinned baseline; it is not a fresh rebuild of Mathlib. All 528 historical files
and 232 shared local objects match the preservation manifest. The separate D5B
portable registration passed ten tooling regressions; it is not a new execution
of the D5B proof. Original installation validation records remain unchanged.

Claude and Grok independently returned substantive CLEAR verdicts with no
blocking findings on fixed packet
`5945c5b4b9d980613e7ff1c31261547a50f559799c14250173e325d0e02584ae`.
Both read sources and recorded evidence; neither reran Lean. Optional findings
and documentation clarifications are disposed in SCATTERING_FLUX_FIDELITY_REVIEW.md.
The approved snapshot and historical author validation remain unchanged; the
closure record explicitly records these documentation-only deltas. F1–F4 is
complete at its stated algebraic scope. No physical scattering solve, numerical
conversion value or new physical discovery is inferred.
