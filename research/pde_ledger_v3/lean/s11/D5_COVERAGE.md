# D5 quadratic invariant completeness contract

Authorized after D4C checkpoint `47d4c249` by "great. let's get to work on D5".
Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Status:
local verification PASS; independent Claude and Grok fidelity reviews CLEAR.
The bounded D5.1–D5.4 contract is complete. See
[D5_FIDELITY_REVIEW.md](D5_FIDELITY_REVIEW.md) for dispositions and limits.

## Claim, domain and bounded obligations

The object is the space of real quadratic forms on **all** real 5×5 matrices,
with conjugation `G -> R G Rᵀ`. Native coordinates are row-major
`G_ij = partial_i u_j`, with no symmetry, transverse, positivity or field-equation
restriction on G. This is a pointwise density classification, before any
Euler–Lagrange or total-divergence quotient. Dimension five is supplied.

| Item | Required conclusion |
|---|---|
| D5.1 — complete density space | Every SO(5)-invariant quadratic form is uniquely `a (tr G)^2 + b tr(G^2) + c tr(G Gᵀ)`, and every such form is invariant under the full O(5) group. The theorem covers all real coefficients, including zero and negative values. |
| D5.2 — parity and census | The actual SO and O subspaces coincide and each has dimension three. The reflection-odd subspace inside them is zero, with dimension zero. Connect the empty odd basis to native `P_D=0`; this is not a claim about bulk equivalence of the three even densities. |
| D5.3 — compact native identification | Identify the complete computed Q9 V1/V2/V6 spans and reflection action against the proved forms, using the actual row-major monomial ordering. Dimensions alone are insufficient. Record the generator orientation and the normalization of the three forms. |
| D5.4 — controls and review | Audit standard axioms and reject meaningful wrong counts, incomplete spans, nonzero odd forms and normalization claims, with admissible positives. Validate source/dependency/instrument/object correspondence; obtain two independent non-author fidelity reviews and resolve findings, then stop. |

## Proof method and evidence boundary

Reuse the finite-coordinate completeness method and trace identities of the
completed D3/D4 contracts without changing their reviewed sources. A complete
quadratic representation has 325 coefficients. Explicit rotations impose
necessary equations; Lean checks the coefficient representation, every retained
equation and reconstruction. Trace identities establish sufficiency for every
orthogonal matrix. The generator selects rational certificates independently
of the native Q9 computation and is not a trusted theorem oracle.

In odd dimension, `R` and `-R` have the same conjugation action and opposite
determinants. Use this compact identity to connect SO(5) and O(5); no new
five-dimensional determinant expansion is needed. It does not alone establish
the completeness of the three-dimensional invariant space.

The native check may execute selected original Q9 helpers for D5, including a
wrong-orientation control in memory. It must not import the production module,
run its driver or regenerate exports. Full-span comparison and non-invariant
same-count controls identify the space. Wolfram anchors are source inspection
unless an actual engine run is explicitly recorded. No kernel-certified CAS
translation is claimed.

Keep one Lean worker, `-j1 -M4096`, strict warnings and a 600-second process limit
with whole-process-group cleanup. Split certificates into bounded compilation
units as needed. A timeout or instrument error is never a rejected mathematical
mutation. Recorded reuse requires transitive sources, pins, generator, commands
and input/output object hashes; focused builds are diagnostic only.

## Review, preservation and completion

Prior completed proof sources/evidence, current D4C objects and other-session
files remain unchanged. Record a preservation manifest before starting. Shared
lakefile/README additions belong to D5 and do not silently revise old packets.
No S11c source, calculation, pinned export or production job is modified.

Local verification passed in recorded run11: 45 guarded canonical objects,
51 standard-axiom audits, thirteen fresh mathematical rejections and seventeen
fresh positives. Compact native full-span/reflection checks and live-source,
dependency, object and preservation validation passed. See
[D5_VERIFICATION.txt](D5_VERIFICATION.txt) and
[D5_FIDELITY.md](D5_FIDELITY.md).

The user approved the fixed 64-file packet. Both independent reviews returned
substantive CLEAR reports with no required fix. Author validation confirmed
archive/snapshot/transport correspondence and unchanged proofs, instruments,
objects and protected inputs. Optional suggestions are dispositioned in
[D5_FIDELITY_REVIEW.md](D5_FIDELITY_REVIEW.md); only closure documentation changed.

Stop after D5.1–D5.4. No D5 bulk variation, general-dimensional classification,
general null-Lagrangian theorem, new spectrum/stability/interface/scattering/pole
work, production reruns or systematic CAS bridge. No commit is requested.
