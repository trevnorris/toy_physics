# D3 bulk variation: object and fidelity map

Scope: [K1–K4](D3_BULK_COVERAGE.md), authorized after `0ecf30a8`.
Status: **bounded K1–K4 complete**. Local verification and current hashes pass;
Claude and Grok independently returned CLEAR with no blockers. See
[D3_BULK_FIDELITY_REVIEW.md](D3_BULK_FIDELITY_REVIEW.md) for reports and dispositions.

## Supplied object and conventions

The completed J1–J4 theorem classifies every real quadratic form on all real
3×3 gradients under SO(3), which here equals O(3) invariance. This increment
uses that same `S11D3Invariants.invariantForm v`, with `v=![a,b,c]`:

`Q(G)=a(tr G)^2+b tr(G^2)+c tr(G G^T)` and `L=-Q/2`.

`density_identity` identifies the new density with this classified object;
`all_invariant_densities` reuses its unique-coefficient theorem. Coefficients
are arbitrary constant real stiffness coefficients. No positivity, nonzero
coefficient, symmetric-gradient or transverse assumption is introduced. The
coefficients share the supplied stiffness units of `mu` and `B`; Lean models
them as real scalars and does not prove dimensional consistency.

The row convention remains `G_ij=partial_i u_j`. Lean spacetime is
`(t,x1,x2,x3)`, with `fieldJet u x j i=partial_j u_i`. The native functions use
arguments `(x1,x2,x3,t)`. The temporal jet does not enter this density: its
momentum is zero. Kinetic energy is not introduced by this increment.

## Actual variation and finite contract

`momentum` is defined by differentiating L with respect to each jet entry.
The proved spatial momentum is
`P_ri=-(a delta_ri div u+b partial_i u_r+c partial_r u_i)`.
The local variational derivative is `-div P`, hence

`EL=(a+b) grad(div u)+c Delta u`.

Smoothness supplies commuting mixed derivatives through mathlib's symmetric
second Fréchet derivative theorem; it is not a separate field hypothesis.
The existing S10 integration-by-parts and compact-test machinery supports the
finite relative-action expansion and its first derivative. Background fields
need not have integrable total action.

Within this full invariant family, the proved exhaustive null locus is
`c=0, a+b=0`. The corresponding current is
`F_i=sum_j(u_i partial_j u_j-u_j partial_j u_i)`:
`div F=(div u)^2-tr(G^2)`. This identifies a divergence, not a pointwise zero
density. Physical boundary conditions and variations that reach a boundary
remain outside the bulk statement.

`responseMap v=![c,a+b]` is a linear map with a one-dimensional kernel and
two-dimensional range. `bulkEquivalent_iff` proves equality on every smooth
field and point, not a test at generic wavevectors. Smooth plane waves serve
as admissible separating examples for necessity, while the full local formula
supplies sufficiency. The action/EL equivalence supplies nonzero compact first
variations outside the null locus.

For phase `k·x-omega*t`, the density-derived modal operator is
`M a_amp=-c|k|^2 a_amp-(a+b)(k·a_amp)k`. It agrees with
`S11Homogeneous.modalOperator 0 c (a+b+c) omega k a_amp` for all k, including
zero. The longitudinal stiffness is `a+b+c`, the transverse stiffness is c.
No positivity or new frequency/root conclusion follows from this identity.

## Native boundary and controls

`_measurements/S11_lean_d3_bulk_source_check.py` executes only ten selected
original native helper definitions for D3. It reads the actual Q9 V1 basis,
solves its invertible basis change to the three trace forms and compares every
native V5 vector, including arbitrary real linear combinations. A basis index
is not assumed to identify a particular form. The native EL helper uses
positive divergence of momentum, the negative of the Lean variational sign.
V5 for Q therefore equals twice the displayed Lean EL for L=-Q/2. The
actual native basis rows in trace-form coordinates are
`[[0,1,0],[1/2,-1/2,0],[0,-1,1]]`, so general native weights are
`(a+b+c,2a,c)`. The symbolic check calls `.doit()` on mixed derivatives; its
smooth-field interpretation is matched by the explicit Lean commutation proof.

The check also compares the plane-wave operator with the zero-inertia
homogeneous stiffness map. Native sign/factor mutations and target-side false
null/wrong-map controls test this compact correspondence. Wolfram evidence is
source inspection only. No production audit, package sweep or export is run.
This is tested symbolic translation; no kernel-certified CAS parser is claimed.

The formal suite covers density sign and half-factor, divergence-current sign,
the coefficient map, false null/equivalence claims and the exact dimensions.
Passing controls include a nonzero null density, negative coefficients, a
nonzero compact first variation and both polarization directions. Instrument
errors and timeouts cannot count as mathematical rejection.

## Provenance and completed review

Five new modules and their axiom root passed sequential verification. Fifteen
unchanged local imports were built in dependency order in run1 and reused with
hash guards in run2; earlier H/I/E/J records stay historical. Reuse requires the local
transitive source closure, pins, generator, command, input objects and output
hash to agree. Keep one Lean worker, with allocator limit 4096 MiB and a
600-second per-process timeout; the allocator cap is not an OS RSS cap.

The local evidence is recorded in `D3_BULK_VERIFICATION.txt` and
`../../_measurements/S11_lean_d3_bulk_validation.json`: 49 standard-axiom
audits, twelve mathematical rejections and eleven positives, with live source,
object, instrument, log and package-pin correspondence validated.

The user explicitly authorized the fixed 41-file packet. Both independent
non-author reviews returned CLEAR. They inspected statement fidelity and
mathematical identities; neither recompiled Lean, reran native instruments or
recomputed hashes. The author revalidated snapshot/archive/transport and live
evidence correspondence. Optional additions were assessed against the contract
and are not closure conditions. No proof or instrument changed after review.
The frozen packet is retained; current closure evidence is in
`../../_measurements/S11_lean_d3_bulk_closure_validation.json`. Completed
H/I/E/J evidence and all S11c files and pinned exports are preserved.
