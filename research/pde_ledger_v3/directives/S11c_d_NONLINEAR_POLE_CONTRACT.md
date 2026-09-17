# S11c-d nonlinear-pole contract correction — nonlinearPoleV2

**Authority and scope.** This is the user-approved correction to the nonlinear
pole contract in `S11c_d_SHARED_PHYSICS.md` v10, approved after the exact issue
checkpoint `d9d4787e`. It governs new SymPy and blind Wolfram work. It supersedes
the unrestricted residue/projector and pole-promotion prescriptions in §3b and
their repeated requirements in §§4, 7, 8 and the build instructions. All other
physics, supplied premises, reduced-representation requirements and output
obligations remain in force. This approval is not a claim of independent review
or of a completed Lean proof.

The active finite-action production pins the original v10 file, SHA256
`fd76447db9c6f2b31d7a72219f4705052023c45788392ec70d761f03aea626c5`.
That baseline stays byte-identical during the run. This addendum is the effective
correction for new work, not a silent relabeling of old artifacts. Future pole
construction and exports must pin **both** the baseline and this addendum until
an explicitly joined consolidation replaces them.

## 1. Analytic domain and computed existence

For the retained §1c-reduced pencil L(omega): X -> Y, fix the profile, parameter
map, tangential momentum, analytic sheet and frequency domain. Name fixed complex
domain/codomain spaces and norms. For an unbounded realization use a justified
common domain and topology; frequency-dependent boundary conditions require an
explicit analytic identification with fixed spaces. An outgoing continuum
prescription alone does not establish a bounded inverse on an unspecified L2
space.

For finite matrices require an analytic regular pencil (determinant not
identically zero). For the nonlocal operator state and justify an analytic
Fredholm domain of index zero, an invertible point in its connected analytic
domain, and the resulting finite-type meromorphic singular structure. Declare
thresholds, branch points/cuts, denominator and domain failures explicitly.
An unresolved analytic premise or singularity classification remains unresolved.
Generic witnesses do not establish these premises on exceptional loci.

Compute the bounded search region, its boundary and candidate isolating contours.
Retain the original pole-search precision, quadrature/domain refinement and
sheet/normalizability/width/all-channel-closure requirements. A resolved empty
set is limited to its certified search domain; failed or unexecuted searches
do not produce empty sets. Constant-end normal-momentum poles remain distinct
from profile-dependent frequency poles.

## 2. Canonical singular data for every resolved pole

At an isolated resolved frequency omega_*, compute the **entire principal part**:

```text
delta = omega - omega_* ,
L(omega)^-1 = sum_{j=1}^p C_{-j} delta^-j + H(omega) ,
H holomorphic near omega_* ,   C_{-p} nonzero ,
R_* = C_{-1}: Y -> X .
```

Emit p, every coefficient, geometric multiplicity, algebraic multiplicity and
partial multiplicities when established, plus the full computed nullspaces and
their completeness evidence. Keep inverse pole order, determinant/characteristic
multiplicity and nullity as distinct quantities. A zero residue C_{-1} does not
imply absence of a pole or of singular coupling. Verify Laurent reconstruction,
left/right operator identities and contour moments with actual operands.

The logarithmic-derivative integral is a separately typed object:

```text
J_Gamma = (1/(2 pi i)) integral_Gamma L(omega)^-1 L'(omega) d omega : X -> X .
```

It is **not a projector in general**. For one isolated pole its residue is the
ordered sum `sum_{j=1}^p C_{-j} L^(j)(omega_*)/(j-1)!`; derive and check this from
the computed Laurent/Taylor operands. In finite regular analytic matrices its
trace counts determinant zeros with algebraic multiplicity inside an admissible
contour. In operator problems use a justified finite-dimensional reduction or
finite-type characteristic-multiplicity theorem. Do not assume a full operator
determinant, a trace-class integrand or a pointwise operator trace exists.
Count, rank and idempotency are distinct checks.

## 3. Conditional semisimple modal projections

For **one isolated semisimple zero**, let V: C^m -> X span the entire right
nullspace and W: Y -> C^m span the left annihilator of the range of L(omega_*).
Under §1's hypotheses compute

```text
D = W L'(omega_*) V ,                 D invertible ,
R_* = V D^-1 W : Y -> X ,
P_X = R_* L'(omega_*) : X -> X ,
P_Y = L'(omega_*) R_* : Y -> Y .
```

Here invertibility of the full pairing supplies the semisimplicity condition
in the stated finite-type setting. Emit D, its rank/inverse evidence, and full
nullspace residuals before using this construction. A numerical rank decision
retains its precision/tolerance and separation domain. A singular or unresolved
D cannot be made invertible by dropping basis directions or assigning a
normalization.

Both P_X and P_Y are conditional local modal projections; compute their
idempotency, rank, range and basis/coordinate covariance residuals. For a contour
enclosing only this semisimple zero, J_Gamma = P_X. This does not identify sums
over arbitrary nonlinear zeros with a Riesz projector on X. Simple modes are
the m=1 case with the derivative normalization retained from v10. Modal
normalization is separate from the independently derived energy-current and
flux normalization.

The semisimple residue formula and full-pairing criterion are stated in
[Schumacher, Corollary 7.5 and Proposition 7.6](https://arxiv.org/html/2412.15985v1#S7).
Their use here requires the stated analytic hypotheses, not just a fitted matrix.

## 4. Higher-order poles and genuine spectral projections

If the semisimple condition fails, compute the necessary left/right root chains
and partial multiplicities or an equivalent finite-type singular construction,
including completeness and obstructions. The right-chain equations are the
coefficients of L(omega) v(omega); extract them from the actual Taylor operands,
and do likewise on the left. They must reconstruct the full principal part
in §2. A single kernel vector is insufficient for a chain or higher-dimensional
space. Keep unresolved chain/completeness questions explicit; do not exclude
defective poles from the required construction.

Use the term **Riesz projector** only with a justified spectral realization.
For example, for a frequency-independent closed operator A on a named state
space Z and a contour in its resolvent set, the Riesz projector is

```text
P_Z = (1/(2 pi i)) integral_Gamma (omega I - A)^-1 d omega : Z -> Z .
```

Provide the actual linearization/realization, domain, spectral correspondence
and analytic physical reconstruction maps, such as a verified local identity
`L^-1 = E(omega) (omega I - A)^-1 F(omega) + H_0(omega)` with H_0 holomorphic,
E: Z -> X and F: Y -> Z. Show that the realization introduces no unexplained
spurious or canceled poles in the claimed correspondence. A descriptor pencil
requires its own justified projection formula; do not silently omit its mass
matrix. Retain the nilpotent/root-chain data needed for higher-order response.
The compressed map E P_Z F is not automatically idempotent or even an X -> X
map. Different realizations must be compared through their verified physical
maps; raw state-space matrices are not canonical cross-engine equality keys.

The canonical common data are the physical inverse principal part and its
typed response, with conditional modal or realized Riesz data clearly labeled.
See [Beyn's nonlinear contour formulation](https://arxiv.org/abs/1003.1580) for
the distinction between nonlinear inverse principal parts and linear spectral
data. An unavailable realization does not authorize inventing P_Z, but it also
does not erase already computed physical principal-part data.

## 5. Spectral overlap has typed forcing and observation maps

Name the transverse incident amplitude space U, the forcing injection
B_T(omega): U -> Y and the requested field/channel observation O(omega): X -> V_o.
Compute spectral overlap from the singular part of

```text
O(omega) L(omega)^-1 B_T(omega) : U -> V_o .
```

For analytic B_T and O at a semisimple pole, its residue is
`O(omega_*) R_* B_T(omega_*)`. At a higher-order pole their Taylor derivatives
also enter the Laurent coefficients and must be retained. State any different
regularity/domain if these maps are not analytic. Applying an X -> X projection
directly to a Y-valued forcing is not licensed by square matrix sizes; a source
lift must be given and justified. Emit basis changes, dimensions, the operands
and reconstruction residuals. These remain **spectral overlaps**, never capture
probabilities/rates without a specified physical protocol.

The bound/physical-sheet, normalizability, zero-width and all-channel-closure
tests remain separate computed objects. Principal-part rank or algebraic count
alone establishes none of them.

## 6. Truncation and perturbative promotion

Poles and singular data are initially results of the retained first-shape-order
operator. For an analytic omitted remainder DeltaL: X -> Y on and inside an
admissible isolating contour, use compatible fixed spaces and norms and compute
a bound

```text
q = sup_{omega in Gamma} || L_ret(omega)^-1 DeltaL(omega) ||_{X -> X} < 1 .
```

Together with analytic regularity/Fredholm premises throughout the homotopy
L_ret + t DeltaL, 0 <= t <= 1, this gives contour invertibility and preservation
of the enclosed **total algebraic characteristic multiplicity**. It does not
preserve the number of distinct zeros, individual nullities, semisimplicity or
chain lengths. Rank preservation of a genuine Riesz projector requires a
justified continuous spectral realization and an admissible contour there.

Quantitative displacement or residue/subspace error bounds require additional
separation, conditioning and derivative estimates. The inequality alone is not
a universal linear-in-error pole-motion bound. Refinement at a few parameters
does not establish uniform contour/operator bounds. Without these premises and
actual remainder bounds, keep the results explicitly at truncated-model status.

## 7. Emission, export, controls and provenance

Retain the tag family `S11CD_BOUND_POLE_SET_AND_RIESZ_DATA` as a compatibility
container, with explicit `nonlinearPoleV2` schema and typed records. Its name
does not assert every candidate has a Riesz projector. Its required computed
content is: search/premise domains; pole locations and multiplicities; all Laurent
coefficients; conditional full-pairing modal projections; root-chain/completeness
data; a Riesz realization and its maps when claimed; logarithmic-derivative/count
diagnostics; and the separate physical-bound tests. Unresolved constructions have
explicit status and evidence, not zero or empty replacement payloads. The
spectral overlap family binds the physical response of §5.

The same typed physical data replace every v10 shorthand for "normalized Riesz
residues/projectors" in the emit/export lists, comparator description and build
brief. Export only the required downstream roots and their recursive bind
closure, using fresh injective write-keys. Retain transparent symbolic operators
and scoped numerical records as required by the retained solver/export contract.
New pole builds and exports include this addendum in their input digests.

Join pole records by profile, parameters, sheet, contour and multiplicities,
with compatible source/field frames and full principal-part/response operands.
Join realized projections only after an explicit state-space correspondence.
Do not compare a multiplicity count with a projector rank or a determinant value
with a canonical response. No comparator or blind engine is launched by this
repair.

Every physical emitted component retains its restored [L,T,M] unit and computed
(eps,eta,sigma_W)/lambda order or explicit nonanalytic/retained-operator grade
status. In particular `[C_-j] = [omega]^j [L^-1]` as a map Y -> X, while modal
projection entries have the relevant field or source unit ratios. Do not label
every matrix entry dimensionless merely because the whole map is an endomorphism.

Controls re-enter actual pencil/forcing operands. Include simple, full semisimple,
higher-order and multi-root contours; basis/coordinate changes; singular-pairing
rejection; source/field type controls; full-principal-part and overlap checks;
and splitting under an admissible small perturbation. Exact synthetic controls
validate mathematical distinctions, not the physical S11c-d pole set. Preserve
exceptional strata, numerical conditioning and all original supplied/debt limits.
