# D4C.1–D4C.4 fidelity record

Recorded local verification passes; both independent fidelity reviews are CLEAR
and bounded D4C.1–D4C.4 is complete. Governed by `../FORMALIZATION_POLICY.md`
and `D4_BULK_COVERAGE.md`. See
[D4_BULK_FIDELITY_REVIEW.md](D4_BULK_FIDELITY_REVIEW.md) for the fixed packet,
review identities, limits and dispositions.

The supplied full D4 quadratic density is the already classified SO(4) family
`a (tr G)^2+b tr(G^2)+c tr(G G^T)+beta P`, with `G_ij=partial_i u_j` and
`P=F01 F23-F02 F13+F03 F12`, `F=G-G^T`. The action is `L=-Q/2`.
The coefficients are arbitrary constant real parameters. This is a spatial
action on the existing smooth spacetime fields, with zero time momentum and
no inertia. Coordinates and coefficients carry the existing dimensionless
convention; no extra dimensional rescaling is introduced. Lean spacetime index
zero is time and spatial derivative row `i.succ` corresponds to native `x_(i+1)`.
Native coordinate functions list spatial arguments before time; this explicit
index identification is retained, and the native time momentum is zero.

`Even.lean` instantiates the previously reviewed D3 even-family argument in
four spatial dimensions. `Action.lean` adds the unchanged reviewed D4 odd
Lagrangian and proves equality to the full classified density. Its momentum
is defined by an actual jet derivative. The first variation in `Variation`
is the derivative of the actual integral of the density change for a smooth
background and smooth compact variation. No convergence of a separate absolute
background action is assumed. Existing compact-support integration by parts
and the fundamental lemma justify stationarity iff the pointwise local EL
vanishes.

`Bulk` uses the reviewed odd EL cancellation, and proves the even result:
`EL=(a+b)grad(div u)+c Delta u`. No coefficient sign is excluded.
`Census` uses this identity for sufficiency of bulk equivalence on all smooth
fields. Transverse and longitudinal smooth waves establish necessity; they
are witnesses for a universal operator theorem, not a restriction of its domain.
The linear response map is `(a,b,c,beta)->(c,a+b)`. Its image and kernel both
have dimension two; the null family is exactly `(t,-t,0,beta)`.

The explicit even current has divergence `(div u)^2-tr(G^2)`; the unchanged
odd current has divergence `P`. Their weighted combination represents the
whole null family. Density, momentum and current witnesses remain nonzero.
The new even current witness is an actual smooth affine field with divergence
2. The reused odd current normalization is a value/jet witness with current
component 1; no separate odd affine-field witness is claimed in the new module.
The odd divergence theorem itself holds for every smooth field. This clarifies
the frozen review prompt's overly broad plural wording about affine witnesses.
Thus the result concerns bulk variation with compact tests, not arbitrary
boundary-value problems. No general null-Lagrangian classification or variable
coefficient result is asserted here. The homogeneous parameter identity
`rho=0,mu=c,B=a+b+c` connects the modal operator without a new spectrum analysis.

The native instrument selects existing SymPy Q9 and coordinate/EL helper
functions using AST extraction. It does not import or run the production
module. An exact invertible basis change identifies the full computed D4
basis with the four trace/orientation forms. It compares every basis vector's
V5 response and the combined actual `L=-Q/2` operator. Native helper EL uses
`+div(momentum)`; Lean's first-variation EL uses `-div(momentum)`, so the native
expression equals negative Lean EL. The instrument also checks the full
momentum, zero time row, response rank/kernel and even current identity.
Wolfram anchors are source inspections only, not a Wolfram execution.
This compact translation is outside the Lean kernel.

Canonical source mutations are not claimed. The paired mathematical
controls change concrete true values/propositions for coefficients, nullness,
operator equivalence, dimensions, action/momentum/current signs and factors,
and longitudinal/transverse responses. Each accepted false control has
exactly one diagnostic in `contract_control`, reducing to `False`, with no
warning, import, syntax, resource or timeout failure. General null-family,
negative-coefficient, zero-wavevector and nonzero-first-variation positives
supplement the paired positives.

Historical proof/evidence is preserved. The preservation manifest records
read-only native/generator inputs, prior Lean sources/evidence and completed
NP/T1/VC objects that need no rebuild. Imported older objects may change during
fresh recorded compilation; historical reports retain their original hashes.
The shared lakefile and this README's new entry belong to this increment.

The first focused builds required explicit reduction of the equal gradient
definitions and pointwise function addition when splitting the full EL into
its even and odd parts. These are proof-code repairs, with definitions and
mathematical statements unchanged. Focused builds are diagnostic only and
were not accepted as recorded reuse evidence.

Run1's supervisor stopped before the recorded suite on the concrete census
example: the third entry of `[1,-1,0,0]` had not reduced. Explicit
`Matrix.cons_val_two` closes that proof. The Controls witness proofs also
explicitly reduce `Fin.cases` and finite-vector entries. Definitions and
theorem statements are unchanged. No mutation evidence was produced by this
failed preliminary build, and run2 subsequently compiled every module freshly.

## Recorded verification and compact native map

Run2 freshly built all 33 objects (26 unchanged imports, six new modules and
audit root), passed 67 standard-axiom lists and all 36 control executions.
The sixteen false statements each have one diagnostic, `unsolved goals` with
`False` in `contract_control`, and no warnings or instrument errors. The twenty
positive executions contain eighteen distinct sources because the even and odd
momentum normalization positives are each paired with both sign and factor
mutations. No canonical source-replacement controls are claimed. The tested
process-group timeout cleanup is unchanged; no timeout occurred in run2.

In the native computed V1 order, the exact basis-to-target rows are
`(0,1,0,0)`, `(1/2,-1/2,0,0)`, `(0,-1,1,0)`, `(0,0,0,1)`.
Their determinant is `-1/2`. Thus target coefficients `(a,b,c,beta)` correspond
to native basis weights `(a+b+c,2a,c,beta)`. The second native basis vector is
half the even null density; the fourth is the odd density. Native `P_D=P`
and the full action momenta match, so a vanishing odd EL expression does not
hide the odd normalization. Native V5 is the helper EL of `Q`, hence twice the
Lean EL of `L=-Q/2`. Mixed derivatives are simplified under the same smooth
field interpretation. All four actual basis vectors are compared.

The native instrument passes twelve identities, catches independent sign and
factor mutations of the original native EL helper, and passes seven target
checks: a wrong homogeneous coefficient map, two wrong odd momentum formulas,
and four admissible null/non-null/nonzero witnesses. Those seven entries are
not seven negative controls. No Wolfram execution or production run is claimed.

`_measurements/S11_lean_d4_bulk_validation.json` validates live sources, pins,
all input/output objects, saved logs and instruments. Historical proof/evidence
and the protected completed NP/T1/VC objects match their preservation manifest.
The external baseline is fifteen clean pinned package repositories and twelve
direct Mathlib source/object pairs. The complete transitive Mathlib cache is
not rebuilt by this increment. The existing D4 generator passed `--check`
without changing its certificates. These validations are author evidence;
reviewers must state whether they independently checked them.

Both completed reviews checked statement fidelity and supplied algebra by
source reading, without independently running Lean, native instruments or hash
checks. Neither identified a required mathematical correction. The optional
stronger controls/theorems were dispositioned within the existing contract;
they do not enlarge it. The author revalidated all reviewed sources, objects,
logs, package pins and protected historical inputs at closure. Only status and
evidence-clarification documents changed after review, so no proof build or
mutation rerun was needed. The approved snapshot/archive/transport and all
canonical sources and instruments are unchanged.
