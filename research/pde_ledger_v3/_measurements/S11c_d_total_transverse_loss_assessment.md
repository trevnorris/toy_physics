# Total leading transverse loss: conditional feasibility assessment

2026-09-28. **Preparation only; no scientific launch authority.** This assessment
used no new science, payload restoration, worker repair or root continuation.
The subsequent user instruction authorizes amendment reviews, then only the
transverse/face repair and its reviews, followed by a costed outgoing-field/power
plan and a go/no-go stop. No scientific runs are authorized in that sequence.
The accompanying [short draft](../directives/S11c_d_TOTAL_TRANSVERSE_LOSS_AMENDMENT_DRAFT.md)
is a proposed scope change, not a cleared replacement for existing directives.

**Assessment: viable in principle, conditional in this model.** It can remove
the need to label every damped thickness root before stating a total loss.
It does not remove the need for an actual first-order outgoing response and
its source-derived power balance. Those are not currently supplied by the
saved separated-point Green candidate or the formal forcing inventory.

## What the headline means

Use loss **out of the transverse sector**, normalized by physical incident
current. Count all reflected and transmitted transverse polarizations as
surviving light. Thus this differs from forward-beam extinction, which includes
elastic reflection. The proposed total includes non-transverse escape, physical
interface absorption and acoustic bulk radiation, counted once on one control
volume. It is a stationary flux fraction, not a temporal lifetime or photon
capture probability.

Thickness conversion followed by absorption must not be counted twice. On a
finite volume, count non-transverse power that crosses its boundary and
dissipation inside it; extending the volume can move power between those two
entries. A supported total may be invariant even when that split is not.
Actual lateral terms, bulk tails, external work and limit conventions matter.

## When first-order fields suffice

The small expansion parameter is the existing contrast path, not the physical
permeability. Keep the latter fixed. Let the complete state and each relevant
loss-producing trace/amplitude have expansions

```
Psi(lambda) = Psi0 + lambda Psi1 + ...
a(lambda) = C(lambda) Psi(lambda)
          = a0 + lambda a1 + ...
a0 = C0 Psi0,     a1 = C0 Psi1 + C1 Psi0.
```

If the actual baseline satisfies `a0=0` for every contributing loss channel,
and its supported power form is regular,

```
B(lambda)[a(lambda),a(lambda)]
    = lambda^2 B0[a1,a1] + higher orders.
```

No second-order field occurs in that leading term. This is a conditional
order-counting argument, not a computed slab result. It includes first-order
face shape, normal, closure and end-projection changes through `C1 Psi0`;
using only a fixed trace map on `Psi1` is insufficient. Measures and incident
normalization must be expanded consistently. Keep eta and sigma_W independent
until the declared physical path is applied.

If only a *net* baseline power cancels while individual drives remain nonzero,
or a required form/limit is singular, the conclusion does not follow. In the
general expansion, `B0[a0,a2]`, its conjugate, and map/measure variations enter
the quadratic coefficient. Our governing
[retained-order rule](../directives/S11c_d_SCATTERING_FORM_AMENDMENT.md:269)
already makes this distinction.

The same obstruction occurs if one computes `1 - outgoing transverse power`
from only a first-order transmitted field: its nonzero baseline interferes
with the omitted second-order transverse correction. An optical-theorem
identity can determine that loss-side coefficient without constructing the
second-order field, but cannot justify dropping its interference. This order
distinction is illustrated by the first-Born/forward-amplitude discussion in
[UT Austin's optical-theorem notes](https://web2.ph.utexas.edu/~vadim/Classes/2024f-qft/Optical.pdf),
p. 5. The familiar absorption-plus-scattering/source-work relation is discussed
in [Miller et al., Optics Express 24 (2016), §2](https://millergroup.yale.edu/sites/default/files/files/MillerPol16.pdf).
Neither external electromagnetic/quantum result is a proof for this slab.

## Premises still needed, and a possible shorter calculation

- A real, regular incoming transverse channel with nonzero physical current,
  and lossless transverse transport at both ends. Test actual face/affinity
  drives, signed power and physical current; zero TH/HT blocks alone do not
  establish losslessness. This preserves the transverse-first/face repair.
- Baseline-free loss amplitudes/port drives and the relevant current
  orthogonality, not merely a zero net loss. Include the actual first-order
  geometry and closure maps and their source joins.
- A passive, stationary balance at the claimed settings: no unaccounted power
  from the held background or other external source, no secular accumulation,
  and disjoint boundary/dissipation accounting. Otherwise retain signed
  transfer and the unresolved work; do not force a positive loss.
- A supported outgoing first-order field, including its bulk/face state, and
  a regular weak-contrast domain away from unresolved thresholds, resonances
  and divergent tails. No classification label can replace this response.

There is a concrete source-work route worth reviewing: evaluate the power
supplied by the first-order forcing to the non-transverse receiving field,
using the actual row-power map, and establish equality to its boundary flux
plus interface dissipation. This need not individually normalize every
damped eigenmode. If the forcing also scatters transverse light, that elastic
part must be separated; total source work is not automatically transverse loss.
Do not insert a generic Euclidean `Im(f†Gf)` with guessed sign or normalization.

The native source already constructs `PLUS_ROW_POWER_MAP`, separate interface
and bulk power, and a finite-volume balance
(`S11c_d_mixing_scattering_sympy_audit.py:3307-3420`). The inherited
[two-port identity](../directives/S11b_SHARED_PHYSICS.md:704) includes transferred-
mass pressure work. These are reusable ingredients, not an established identity
for the new spatially forced field. A finite-point or source-specific check is
the intended scope, not a global operator-certification campaign.

## Current availability and the secondary continuation

Saved symbolic baseline-zero evidence supports trying the route. The
[outgoing candidate](S11c_d_outgoing_prescription_candidate_report.md) still has
unevaluated integrals, separated-point and evanescent-slice restrictions. The
[forcing inspection](S11c_d_two_asymptote_saved_inspection_report.md) did not
compute its outgoing action. The localized-thickness attempt failed before an
action candidate; the old finite omega=1 solve is not already a physical total-
loss result. A general coincident Green formula cannot be substituted for the
missing source pairing. No numerical saving or days-to-result is established;
the benefit is removing a root-label prerequisite, not eliminating the solve.

For optional branch interpretation retain the user's diagnostic
`Lambda_A(omega;s)=s Lambda_A(omega)`, `0<=s<=1`, with the memory time fixed and
physical results at `s=1`. Below the *actual* bulk-radiation boundary, use an
established lossless `s=0` problem; collisions/singularities stay unresolved.
Above bulk availability, **complex-at-s=0 is not alone a leaky label**: it can
also be evanescent or on a wrong sheet. Require outgoing-sheet, nonzero face-
drive and bulk-radiation evidence, or keep the label unresolved. Diagnostic
lineage does not by itself assign a flux or percentage of total loss.

**Recommendation:** review the proposed primary observable and its finite
premise list before implementation. Make the thickness/root continuation
secondary so it cannot hold the headline hostage. If the required forced-field
power pairing needs a new unsupported method, return that specific blocker;
do not expand into a full radiating Green construction automatically. This
assessment supplies no new CLEAR verdict, A11/A12 clearance or launch authority.
