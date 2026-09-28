# S11c-d two-asymptote response: concrete method plan

Status: **source-grounded design; implementation and independent method review
remain required**. The user said **Continue** after the completed outgoing
prescription and the stated next step of two-asymptote planning. This document
continues that design. It does not authorize an external submission or a new
scientific launch. No completed calculation is a reconstruction queue.

## Proposed route and first deliverable

Use the validated full reference inverse to construct exact asymptotic field
corrections at the two inherited real pole blocks, then lift those corrections
smoothly across the interface. Apply the original unbounded reduced operator
to that lift, retaining its local derivatives and all ordered nonlocal terms.
The resulting defect is the forcing for the next response construction.

The first bounded deliverable is **the two-end first-grade lift and complete
forcing defect**, with exact end cancellation and a source-specific disposition
of its admissible convolution domain. It covers the saved LAB_HELD /
RHO4_CONSTANT input and independent grades `(1,0)` and `(0,1)`, at both ends
and both inherited real blocks. It does not yet publish a response, scattering
matrix, mixed response coefficient or FORM. This is one construction task,
not another input-reader campaign or separate jobs for each caveat.

The [first-stage contract](S11c_d_two_asymptote_first_stage_contract.md) fixes
its operands, output/check requirements and stopping rules. The
[source receipt](../_measurements/S11c_d_two_asymptote_plan_sources.json)
pins the definitions and selected saved routes without importing scientific
modules or restoring a pickle.

## What source inspection established

| Existing object | Reuse and concrete limitation |
|---|---|
| Validated reference prescription | Full 5x5 inverse, exact real-pole residues, symmetric PV plus signed-current delta assembly, units and branch evidence are saved. Acceptance is for the fixed-input separated-point candidate; it does not yet define its action on every nondecaying forcing. |
| Native reduced action and assembly | `actions.pickle` retains five-field probe actions; `assembly.pickle` retains actual local coefficient matrices and original ordered nonlocal integrals. These are the whole-line source operands. Do not call `EdgeReduction` or reconstruct the accepted assembly. |
| `EdgeReduction.subtraction/prescribe` | The source already carries the canonical constant/Heaviside/localized profile split and its regulated half-line transforms. The Abel weak limit occurs after the complete convolution. Reuse this convention, including its tail premise and fixed origin. |
| `ConstantEndPencil.background` | Ends are simultaneous translations of the full kernel coordinates, including the profile subtraction. A one-coordinate coefficient limit is not the nonlocal end operator. All three full strong symbols are already saved. |
| `continuum-grades.pickle` | Actual symbolic coefficient grades are reusable. Its `termJoins` inherit a finite-domain source factorization, so their integral bounds do not define a whole-line action. |
| `BoundedSourceFourierAssembly` | `bounded` replaces each infinite limit by a cutoff; `construct` stores both `ORIGINAL` and `BOUNDED` and reorders the bounded source integral. Its accepted certificate explicitly excludes an unbounded interchange. Rejoin original integrals and preserve their order; do not send saved cutoffs to infinity and call that proved. |
| Saved end pencils and invariant pairs | `left/right-pencil.pickle` contain symbolic pencil/radical coefficient tables. Saved `R` and `K` mode variations are numerical arrays with checked numerical residuals. They support comparison and gauge bookkeeping, but a small nonzero tail residual cannot be declared exact cancellation. |
| Saved finite response/current | Supplies the correct dependency pattern, phase order, current forms and incident normalization. Its fixed matrices and endpoint at ±64 do not define the whole-line inverse or its asymptotic extraction. |

The input remains `w(xi)=(1+tanh(xi))/2`,
`m(xi)=(1-tanh(xi)^2)/3`, with its saved parameters and units. The actual
reference grades remain zero; retain eta and sigma independently until the
specified physical homotopy. The accepted rest-acoustic real axis is still
evanescent. No radiating witness is obtained by this construction.

## Exact end corrections without a new root search

The following is a proposed algebraic route to review and implement, not a
newly evaluated scientific identity for the saved objects.

Let `P0(k)` and `R0(k)=P0(k)^(-1)` be the saved full strong reference symbol and
inverse. Write the saved full end symbols as `P_e(k;eta,sigma)` for
`e=LEFT,RIGHT`. Their coefficient at zero grade must join `P0` using the
physical branch, units and normal-momentum coordinate. Use the existing
symbolic coefficient tables where they match these operands; do not call the
old end-mode constructor. Let `P_e,h` denote a first-grade coefficient.

For each already accepted real pole `k_j`, use the local Laurent series of

```text
H_e,00(k) = R0(k)
H_e,h(k)  = -R0(k) P_e,h(k) R0(k),     h = (1,0),(0,1).
```

This is an ordinary local meromorphic matrix product on the inherited branch,
not multiplication of the PV/delta distributions. It introduces no new poles
or root search. The local square-root expansion must join the saved physical
branch value and relation. Do not select a different sign at a nearby point.

If `A_j` is the saved residue and the actually computed first-grade principal
part is `C_e,h,-2/(k-k_j)^2 + C_e,h,-1/(k-k_j)`, define

```text
E_j,00(z)  = exp(i k_j z) A_j
E_e,j,h(z) = exp(i k_j z) [C_e,j,h,-1 + i z C_e,j,h,-2].
```

Equivalently these are coefficients of `(k-k_j)^(-1)` in
`exp(i k z) H_e,g(k)`. This notation does not require a contour integration
or a complex-frequency prescription. At zero grade the residue columns form
an overcomplete five-column frame of the saved nullity-two block. They are
not five independent physical channels or flux-normalized amplitudes. Retain
that redundancy and the source-derived column units explicitly. A physical
channel-coordinate/current-normalization map is a later construction duty.

Verify the new coefficients against the actual end pencil, not merely the
inverse-variation formula used to build them. For a polynomial vector `p(z)`,
the convention `exp(+i k z)` gives the finite polynomial action

```text
P_e,g(-i d/dz)[exp(i k_j z) p(z)]
 = exp(i k_j z) sum_a [(-i)^a/a!] (d_k^a P_e,g)(k_j) d_z^a p(z).
```

Using that separate action, save and check `P0 E_j,00=0` and
`P0 E_e,j,h + P_e,h E_j,00=0`, entry by entry, with complete operands.
Branch transport derivatives and off-diagonal blocks must remain. A
numerical residual cannot certify an exactly vanishing nondecaying tail.
If exact cancellation cannot be established, stop at the saved coefficients
and name that failure; do not substitute the existing approximate mode jets.

This construction reuses known poles/residues but computes genuinely new
exact field corrections. It is not a replay of accepted numerical mode
variations, a new modal census, or clearance of physical end channels.

## Lift and retain the complete forcing

Choose a dimensionless auxiliary smooth partition at the existing origin,
`chi_+(z)=(1+tanh(z/L_W))/2`, `chi_-=1-chi_+`. It is held fixed under eta/sigma
differentiation and changes no physical profile or input. It is separate from
the source's canonical Heaviside subtraction and its profile Abel regulator.

For each reference residue frame and first grade define

```text
T_j,h = chi_- E_LEFT,j,h + chi_+ E_RIGHT,j,h
F_j,h = -L_h E_j,00 - L0 T_j,h
U_j,h = T_j,h + V_j,h,        L0 V_j,h = F_j,h.
```

Only the first two lines are constructed in the first bounded stage. Every
`L` here is the full native reduced action at that grade. In particular,

```text
[L0,chi] E = L0(chi E) - chi L0(E)
```

contains both local product derivatives and the nonlocal difference. If a
source term has a plain kernel representation, its latter part is
`Integral K0(z,z') [chi(z')-chi(z)] E(z') dz'`. The actual implementation must
substitute at the field slot inside each saved integral and retain derivative
terms there; this shorthand must not replace nested or differentiated source
operands by an invented simple kernel.

Save both the direct defect and its commutator decomposition,

```text
F_j,h = -[L_h - sum_e chi_e L_e,h] E_j,00
        - sum_e [L0,chi_e] E_e,j,h.
```

Their equality uses the checked exact end equations. The equality of these
two full expressions is a source reconstruction check. Neither it nor a
coefficient fingerprint alone proves decay/integrability of the remaining
nonlocal action. Preserve both cross-interface directions and the original
nested integration order before any further interchange is justified.

## Domain needed before applying the prescription

The saved 50 tail summaries include powers 1, 2, 3 and 4 of reciprocal
momentum. Thus decay alone cannot be used to claim absolute Fourier
integrability of every entry at coincident positions. The next action should
be defined by the saved momentum PV/delta distribution paired with the
actual forcing transform, with its measure and sign inherited unchanged.
There is no need to invent a pointwise value at `z=z'` to do this, but there
is a need to establish that this particular pairing is defined.

First specify the action on smooth rapidly decreasing forcing. Then determine
from the constructed source expressions whether each lifted defect belongs
to that class, or to a precisely stated weaker class sufficient at the real
poles and at infinity. Polynomial factors from phase derivatives require the
corresponding weighted tails. The general supplied L1 profile-tail premise
must not silently become a stronger weighted-tail premise for every profile.
For this pilot, analyze only the actual saved tanh/sech profile and full
nonlocal kernels; no global profile theorem is a prerequisite.

Keep the profile Abel regulator, real-pole symmetric exclusion and any
auxiliary test regularization distinct. Combine the complete source action
before its prescribed weak Abel limit. Specify any subsequent order of
pairing/limits and show the required existence for these operands; do not
multiply singular distributions formally or exchange limits just because
finite-cutoff checks passed. If the needed bound or pairing is unavailable,
persist the exact offending term and its limits as the next dependency.

Standard test-function/distribution definitions and Fourier transforms are
summarized in [NIST DLMF §1.16](https://dlmf.nist.gov/1.16). They support the
language of this domain requirement, not the physical outgoing choice or
its validity for this coupled source. Native phase signs and Fourier mass,
not the reference's normalization, govern the implementation.

## Following construction, after this ingredient passes

Construct the action `V_h=G_out F_h` only on the established source class,
prove its equation and outgoing asymptotic extraction on that class, and
join the residue frames to physical end coordinates, phases and incoming
normalization. Lift-gauge changes can move homogeneous components between
`T` and `V`; compare the reconstructed field and its fixed incoming state,
not the lift alone. Retain the full continuous Fourier contribution and
closed-channel response, even though only nondecaying real-pole tails are
lifted explicitly.

For the mixed grade, the source remains
`-L10 U01 - L01 U10 - L11 U00`, with both ordered cross terms. Its tail lift
must include end-field variations contracted with lower-grade scattering
amplitudes and the mixed end variation. A formal extension of the local
inverse coefficients is
`R0 P10 R0 P01 R0 + R0 P01 R0 P10 R0 - R0 P11 R0`; this is a dependency
formula, not a calculated mixed response or permission to drop either order.

Only after those constructions should the other three cases, general
parameter/profile bindings, current/incident-denominator coefficients and
own-row FORM exports be completed under their applicable gates. Existing
fixed-input numerical end variations do not supply general differentiable
parameter dependence. A11/A12 and supported bulk-depth flux retain their
separate domain/face-map obligations; the closed development slice cannot
clear them. No new centre mechanics, physical input, pole campaign or
whole-amendment review is proposed.

## Cost and execution boundary

The completed reference constructor cost 19.314 seconds; the last outgoing
continuation cost 3394.435 seconds, dominated by one exact derivative. These
measurements warn against predicting this new algebra's cost from packet size.
Use source-bound local coefficient arithmetic and persistent per-entry
returns; do not replay determinant/cofactor construction or the earlier
expensive derivatives to obtain new local coefficients.

The proposed first job retains the normal 900-second whole-job / 840-second
native cap, 2 GiB, zero swap, one CPU, nice 15, 32 tasks and one native thread,
under the shared guard and existing supervisor. These are execution limits,
not a runtime estimate. The prior indefinite-duration approval applied only
to the completed outgoing continuation. It is not inherited by this stage.

Before any launch, author the actual worker and fresh input manifest, complete
the independent implementation/method review on a fixed explicitly authorized
packet, resolve substantive findings, and pin a fresh gate to the exact final
bytes and applicable user scope. No reviewer submission or worker launch is
part of this planning turn. Preserve original failures, accepted artifacts,
incident history, Lean/shared guard and the protected builder suffix.
