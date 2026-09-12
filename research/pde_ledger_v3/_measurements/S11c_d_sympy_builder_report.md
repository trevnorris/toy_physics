# S11c-d SymPy builder checkpoint

The original mechanical-load mismatch is resolved in the fresh reduced
LAB_HELD/RHO4_CONSTANT reference. All five mechanical coefficients now agree,
the five mass coefficients still agree, and all 30 recorded current/energy
residual scalars are zero. The independent stiffness anchor also agrees.
The existing d formulas pass with the repaired exports. The
[current report](S11c_d_nonlocal_current_report.md) records the scope and
source-pinned evidence; the full repair history is in the
[repair record](S11c_mechanical_repair_report.md).

S11c-b now obtains the face-work row multiplier from the action and stored
stiffness coefficient, preserving the physical generalized force separately.
The native four-case rebuild changes exactly four exported roots. All 188
physical comparison scalars are zero, including kinetic orientation, mass,
chemical, non-face preservation, and separate uniform S11b face-load checks.
The c1 refresh leaves all 44 exported values unchanged.

The rebuilt c2 export changes its closed slab and coupling roots; its
self-energy increment is unchanged in the value census. All 431 serialization
checks completed. The strengthened power reference uses separately computed
constrained energy-term variations and a kinetic action variation, retaining
its dependency on b's energy construction. Four canonical power residuals and
16 kinetic entries are zero; both source controls are detected at all three
PIT samples in each case. This is sampled control detection, not global
nonvanishing. The scoped legacy dependency triage retains its documented
boundary: four trials ran, two encountered the existing inequality error, and
the invariant dimension-binding tag leaves the literal verdict `FAIL`.

The new reference reduction checks 39 scalar reconstruction digests against
literal zero, with resolved dimensions. Its transcript and the repaired
current transcript are published under `scripts/out/`, together with the
upstream and repair-check outputs. The final native four-case d rebuild and
all nine inventory stages are complete. The [full inventory](S11c_mechanical_repair_d_full_checks.json)
records 432 isolated root/lift candidates across 24 packets, complete sampled
basis/pairing ranks, 480 transported paths and 24 unresolved branch-locus paths,
432 constant-end Laurent residues, and 24 right/24 left threshold chain spaces.
Point, slice, bank, contour and exceptional-domain limitations remain explicit;
these momentum poles do not supply the profile-frequency poles of §3b.

The native run took 9,989.25 seconds with 1,734,840 KiB peak RSS and empty
stderr. Source pins are stable. Fourier reconstruction residuals/projections
are zero, and metadata/coverage inventories have no gaps. Lossless round trips
preserve all 229,636 tags and 229,634 indexed source-line assignments. The
regenerated 83,848,737-byte production transcript is published under
`scripts/out/`, with its previous annex payload preserved. The repair is
complete; historical counts were not imposed as expected answers. Next is the
one-case nonlocal current and mode/flux normalization work.

All ten broad engine TODOs remain, including variable-profile nonlocal current,
mode/flux normalization, complete two-ended scattering, profile-frequency
bound poles, survival, bookkeeping, weak coefficients, controls, and own-row
export. Existing exceptional-domain and upstream-debt boundaries remain open.
The supplied physics and retained contract below are unchanged. No authority
change, S10/Lean edit, review leg, comparator, Wolfram run, incomplete export,
commit, or push occurred in this repair turn.

## Retained user-approved solver/export contract


1. Preserve `EdgeReduction`, the positional three-parent fold, and the exact
   direct-lookup manifest. All numerical assembly must consume the computed
   reduced rows, including the full nonlocal terms and full coupling vertex.
2. Accept explicit, independent dimensionless profile functions w(xi), m(xi),
   their derivatives, asymptotic limits and tail information. A selected smooth
   step with an independently adjustable localized modulus bump is a numerical
   instance; it does not replace the interface class. Store the profile formula
   and digest in every case record. The current preflight input is recorded in `S11c_d_channel_preflight_input.json`.
3. Require a complete parameter map in a declared L/T/M unit frame, real
   continuum frequency and tangential momentum, small contrast, and
   sigma_W = eta_bg W_0/L_W for evaluations on the physical homotopy. Retain
   independent eta/sigma grades in the symbolic calculation. Test actual
   reference/end channel availability before attempting flux normalization.
   Do not manufacture an incident channel by assigning a sector label.
4. Compute both-end modes, left/right normalization, the S11b-derived current,
   and the variable-profile matching problem. Re-expand the continuum response
   to the retained rectangle; do not present a finite-contrast numerical
   solution as a higher-order continuum prediction. Retain evanescent matching
   modes, channel degeneracies and domain failures explicitly.
5. Numerical pole searches have an explicit profile, parameter map, sheet,
   bounded search region and isolating contours. Evaluate the retained operator
   without the continuum re-expansion. Record boundary/quadrature resolution,
   domain size, precision, root residuals, contour-count evidence, and changes
   under refinement. A bounded search does not establish a global pole set.
   An unsuccessful or inconclusive search is unresolved, not an empty pole set.
   Compute residues/projectors and sheet/decay/width/closure tests only for
   actually resolved candidates; emit spectral overlap, not capture probability.
6. Separate transparent symbolic expressions from evaluated numerical records.
   Symbolic operator/continuum/weak-coefficient exports remain differentiable
   SymPy expressions, compacted with algebraic equivalence checks. Numerical
   mode and pole datasets retain their input bindings, domain, convergence
   evidence, dimensions and truncated-model status. They are not stand-ins for
   a generic symbolic profile-dependent root function. This is the export
   distinction motivating the user-approved contract; the downstream consumer
   will need to bind the appropriate representation explicitly.
7. Fingerprints summarize already constructed/evaluated objects. They do not
   evaluate nonlocal integrals or replace a spectral solve. Algebraic PIT and
   physical numerical evaluation are separate records. All completed roots use
   fresh lowerCamel write-keys and the existing bind-closure/minimal-delta guards.
8. Finish the one-case path, then implement controls and bookkeeping, then run
   all four cases once and write the complete export. No review legs,
   comparator, Wolfram engine, downstream stage, or commit belongs to this lane.
