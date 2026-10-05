# Variable-coefficient and flat-interface fidelity boundary

Fidelity record for [VC1–VC4](VARIABLE_COEFFICIENT_COVERAGE.md), governed by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). T1–T4 is complete at
`ad365b5f`. The canonical VC proofs passed run1 and the complete control rerun
passed run2. Both independent non-author fidelity reviews are CLEAR; bounded
VC1–VC4 is complete. See [review and closure](VARIABLE_COEFFICIENT_FIDELITY_REVIEW.md). This
increment concerns the reviewed local D3 family and D4 odd density with
prescribed coefficient profiles, plus an actual normal-slice integral identity. It does not identify
the complete S11c reduced operator or prove a multidimensional transmission
problem.

## Object, derivatives and domain

The coordinate convention is G_ij=partial_i u_j (spatial derivative rows).
`Point D` contains time followed by D spatial coordinates; Lean spatial indices
use `i.succ`. `Vec D` contains real components. The prescribed profiles and
background fields in the local identities are smooth; no sign, nonzero,
transverse, field-equation or stationarity premise is imposed. Spacetime
profiles are permitted, including spatial profiles. Profiles are held fixed
under variation. Since these densities contain no time derivative, their time
momentum is zero even for time-dependent fields and prescribed coefficients.

`D3.lagrangian` is the reviewed action

```
L(x,G) = -[a(x)(tr G)^2 + b(x)tr(G^2) + c(x)sum_ij G_ij^2]/2.
```

`D4.lagrangian` is the reviewed odd action `-beta(x)P(G)/2`, with
`P=F01 F23-F02 F13+F03 F12`, `F=G-G^T`, and `M_ij=partial P/partial G_ij`.
P uses the native computed odd-basis normalization, not an assumed epsilon
prefactor. The preceding D4 classification proves this is the one-dimensional
odd density family. The fully summed Levi-Civita contraction on G is 2P; this
is prose/hand algebra here and is recorded in the earlier native D4 check.
Lean proves the derivative-defined dual matrix and its contraction G:M=2P,
not a new theorem about an explicit epsilon tensor in this increment.

In both modules, momenta are defined as actual derivatives of the jet action.
`momentum_identity` connects them to the old jet derivative at the local
coefficient value, and `pointwise_variation` reuses the actual directional
variation theorem. The local EL is defined as minus the coordinate divergence
of these momenta, consistently with the negative sign of L. These are local
EL/pointwise-variation statements; a new integrated first-variation theorem
with variable coefficients is not claimed.

## Product-rule corrections and bulk equivalence

The D3 momentum is
`p_ij=-(a delta_ij div u+b partial_j u_i+c partial_i u_j)`.
`eulerLagrange_eq` proves, for each component j,

```
EL_j = (a+b)partial_j(div u) + c Delta u_j
       + (partial_j a)div u
       + sum_i (partial_i b)partial_j u_i
       + sum_i (partial_i c)partial_i u_j.
```

The old bulk-null direction a=-b,c=0 retains the residual
`(partial_j a)div u-sum_i(partial_i a)partial_j u_i`; it need not vanish.
Constants recover the old local operator. This proves the stated extension
for every smooth profile in that family; it is not a classification of every
profile for which a residual happens to vanish.

For D4, `eulerLagrange_eq` proves
`EL_j=(1/2)sum_i(partial_i beta)M_ij`.
It uses the reviewed divergence identity `sum_i partial_i M_ij=0`.
A constant beta recovers zero EL; variable beta generally does not.

The dimension-independent `weighted_divergence` theorem proves
`div(aJ)=grad a dot J+a div J`. The reviewed D3 current
`J_i=sum_j[u_i partial_j u_j-u_j partial_j u_i]` and D4 current
`K_i=(1/2)sum_j u_j M_ij` then give

```
L_D3(null profile) = -(1/2)div(aJ) + (1/2)grad a dot J
L_D4(odd profile)  = -(1/2)div(beta K) + (1/2)grad beta dot K.
```

Thus replacing a density by a divergence after inserting a profile requires
retaining the displayed correction. Bulk equivalence at constant coefficients
does not remove boundary effects or identify tractions.

## Interface term and precise application boundary

`split_integration_by_parts` is an actual interval-integral theorem. It assumes
separate real functions p_minus and p_plus with specified derivatives on the
closed intervals, a common test function h with specified derivative, and
interval integrability of all three derivative data on the appropriate sides.
It proves

```
integral_a^c p_minus h' + integral_c^b p_plus h'
 = p_plus(b)h(b)-p_minus(a)h(a)
   +[p_minus(c)-p_plus(c)]h(c)
   -integral_a^c p_minus' h-integral_c^b p_plus' h.
```

The interval integrals are oriented, so the theorem requires no order on a,c,b;
the adjacent physical-interval interpretation uses a<=c<=b. The two flux
functions are separately defined representatives on the real line with actual
derivatives at the closed-interval endpoints. This is a sufficient extension
hypothesis, stronger than a theorem using only one-sided weak traces. Their
values at c need not agree. The same h supplies a common trace.
`compact_endpoint_split` assumes only h(a)=h(b)=0; its name does not assert a
compact-support theorem. An explicit polynomial h supplies a nonzero example.

`normalFlux` contracts a real normal vector with momentum, and D3/D4
`traction_eq` identifies it with the same derivative-defined action. Unit-normal
normalization is an application convention, not a premise of the algebraic
identity. With the normal oriented from minus to plus, the interface pairing
is `(p_minus-p_plus) dot h`. `jumpPair_zero_iff` proves that this pairing vanishes
for every finite-dimensional real test value iff the two flux vectors agree.
Physical stiffness traction is the negative of the momentum flux used here.

The integral theorem is scalar and applies componentwise. The finite-dimensional
trace criterion is separate. This increment does not prove a tangential Fubini
lifting, a general multidimensional variational principle, existence of weak
traces, continuity of the field, or that arbitrary traces arise from solutions.
No independent surface action/source or product of discontinuous distributions
is assumed. Such application hypotheses must be supplied before interpreting
this as a full transmission law.

## Compact native connection and controls

`S11_lean_variable_source_check.py` uses the existing selected-source extractor
from the D4B instrument, the original Q9 D3/D4 basis computations and original
coordinate helpers. These native checks use spatial profiles; Lean's local
identities also allow time dependence, with zero time momentum. The instrument
identifies the same D3 trace basis and D4 P_D=P exactly,
then computes the prescribed-profile momentum divergence itself. The new
variable-profile calculation is not attributed to the old constant-coefficient
native EL helper. The native engine, historical reports and extraction helper
are hashed read-only inputs. No production driver, S11c calculation/export or
Wolfram execution is performed. This is exact symbolic translation checking,
not kernel certification of the instrument.

The seven native identities concern the D3 density, EL and weighted current;
the D4 density, EL and weighted current; and a nonzero normal-slice integral.
Ten native controls cover omitted gradient/current corrections, D3 gradient
sign and derivative-index transpose, D4 factor/sign/false zero response, and
omitted/reversed interface terms. The wrong-index witness uses b=x1, a=c=0,
u=(x2,0,0): actual EL=(0,1,0), while the transposed b-gradient formula gives zero.

The formal suite has twelve paired false statements and sixteen positive
executions. They are concrete statement mutations, not automatic replacements
inside canonical general proofs. Every accepted false statement must reduce to
exactly one `False` goal in `contract_control` for this accepted evidence.
The runner alone accepts at least one intended diagnostic; author validation
separately checks that every recorded rejection has exactly one. They use actual smooth-field
responses, nonzero momentum tractions, an actual interval integral and finite
trace pairings. No syntax, import, timeout or resource failure counts as a
mathematical rejection. Constant profiles and equal traces are explicit passing
cases. Native derivative-index sensitivity is recorded separately from these
twelve Lean pairs. Lean's D3/D4 response witnesses prove the specified nonzero
component; the native instrument additionally checks their full vector values.

| Obligation | Principal declarations | Sensitivity and admissible examples |
|---|---|---|
| VC1 | `D3.momentum_identity`, `pointwise_variation`, `eulerLagrange_eq`, `null_profile_residual`, `constant_profile` | Actual u=(0,x2,0),a=x1,b=-x1,c=0 gives EL component 1; omitted/reversed correction rejected; native transpose control; constant null profile passes |
| VC2 | `weighted_divergence`, `weighted_density`, `D3.weighted_null_density`, `D4.eulerLagrange_eq`, `weighted_odd_density`, `constant_profile` | u=(0,0,0,x3),beta=x1 gives EL component 1/2; missing/sign/factor controls; nonzero weighted correction; constant beta passes |
| VC3 | `split_integration_by_parts`, `compact_endpoint_split`, `jumpPair_zero_iff`, D3/D4 `traction_eq` | p_minus=2,p_plus=5,h=1-x^2 on [-1,1] gives integral -3; momentum signs/factors, omission and reversed jump rejected; equal traces pass |
| VC4 | `Controls.lean`, audit root and compact native check | Full recorded verification passes; both independent fidelity reviews CLEAR |

## Verification history and remaining work

Focused builds are debugging evidence only. The general D3 product expansion,
its full local EL and weighted current passed after routine concrete-vector and
sum-normalization repairs. The D4 gradient response and weighted current passed.
The normal-slice interface module also passed. Concrete witness proofs needed
explicit finite-coordinate reduction, typed derivative function normalization
and explicit interval-integrability proofs. The last focused Controls build
reported only an unnecessary tactic-sequencing linter; that sequencing is now
removed. Run1 subsequently passed all 34 unchanged imports, the five new modules
and audit root, all 31 standard-axiom audits, all twelve paired mathematical
rejections and their twelve passing partners. It stopped at the extra D3
constant-profile positive: the final entry of the literal vector [1,-1,0] was
not reduced. This is a control-proof failure, not a rejected mathematical
mutation. The repair adds only `Matrix.cons_val_two` to that control's
normalization. Canonical proof sources, theorem statements and all forty
compiled objects are unchanged and have been hash-validated, along with the
historical evidence and five preserved T1 objects. No resource limit was raised.

The initial native check passed seven identities and nine controls. To discharge
VC4's explicitly promised derivative-index sensitivity, a tenth native control
compares the correct b-gradient contraction with its transpose. Its smooth
witness is b=x1, a=c=0 and u=(x2,0,0): the correct EL is (0,1,0), while the
transposed formula gives zero. This is a targeted check of an existing coverage
obligation, with no new mathematical scope or native source change. It passed
in run2 together with all seven native identities and the other nine controls.

Run2 passed with all forty canonical builds reused from recorded run1 under
the full transitive source, generator, dependency-pin, command and input/output
object guards. All twelve paired mathematical rejections and sixteen positives
ran freshly. Each rejection has exactly one diagnostic: unsolved `False` in
`contract_control`. The 31 selected axiom declarations match by name and use
only standard axioms or none. There are 68 check records and forty output
objects. No focused object is accepted as recorded build evidence.

Author validation checks current source/instrument/report/log hashes, every
input/output object, all fifteen clean pinned package revisions and all thirteen
direct Mathlib source/object pairs. The transitive Mathlib compiled cache
remains the existing pinned dependency trust baseline; the full library was
not rebuilt. Historical proof bytes/manifests and the five T1 objects match.
Earlier D3/D4 object records remain historical after unchanged import rebuilds;
their reviewed source/evidence was not rewritten. Run2 has no timeout or
instrument failure. Both independent fidelity reviews are CLEAR with no required
mathematical corrections. Non-blocking observations are resolved in the review
record. Supervisor state/logs are in `_scratch/S11_lean_variable/verification_run2`;
per-check sources/logs are in `lean/s11/_scratch/variable_verification`.
The user authorized the completed VC checkpoint and then question 3.
