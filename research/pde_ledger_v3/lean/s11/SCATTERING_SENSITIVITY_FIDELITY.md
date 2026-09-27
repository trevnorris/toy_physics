# S1–S4 sensitivity: fidelity record

Status: bounded S1–S4 COMPLETE; local PASS and independent Claude/Grok CLEAR.
Recorded run6 freshly built eight isolated objects and passed 51 selected standard-axiom audits,
eighteen paired mathematical rejections and twenty-four positive executions
(twenty-three distinct). Author validation re-adjudicated the recorded logs,
source/snapshot/object bindings, native lineage, thirteen direct Mathlib import
source/object pairs, fifteen clean package pins and historical preservation.
The existing transitive external cache remains a pinned baseline; it was not
wholly rebuilt. Both reviewers found no blocking defect. Optional dispositions
and exact documentation-only deltas are recorded separately; proofs, instruments,
approved snapshot/archive/transport and historical evidence remain unchanged.

Earlier attempts are retained as failures. Offline restoration validated 2999
pinned local archives and unpacked nine without downloading or rebuilding old
ledger objects. Run1 passed the compact native checks, then the unchanged
Stability import hit Lean's internal 2048 MiB accounting limit. Actual cgroup
peak was 1,146,425,344 bytes, with no max/OOM/swap events. The internal setting
became 4096 MiB, matching Stability's historical verification; the hard whole-job
cap remained 2 GiB. Bounded proof repairs supplied the complex normed-field
import, made two fixed-term inequality additions explicit, omitted an unused
finite-index instance from a pointwise theorem, and removed unnecessary tactic
sequencing. Formulas and finite-space applications were preserved. Every repair
followed complete source/snapshot/log/object/native/dependency/preservation
validation. No compiler, resource, linter or syntax failure counts as a rejected
mathematical control. No focused or historical object substitutes for a fresh
recorded build. Run5's final extra-positive lexer failure was fixed by whitespace
only in the instrument, without changing canonical sources. Run6's actual guard
peak was 418,717,696 bytes, with no max/OOM/swap events.

## Finite object and normalization

`S11c_d_finite_scattering.py:construct` assembles M with its modal boundary rows,
forms row scales and column scales, solves the balanced system with least
squares, and converts back to coefficients. If R is reciprocal row scaling and
D is column scaling, K=R M D⁻¹, g=Rb, z=Dc. The native scaled equation residual
is Kz-g. A bound for K⁻¹ controls errors in z; conversion to c requires D⁻¹.
The separately saved `unreplacedOperator` residual belongs to a different
matrix and cannot replace the residual of the solved boundary system.

The Lean residual theorem takes an actual continuous linear equivalence
and a bound on its inverse. It does not turn `lstsq` rank, singular values,
condition estimates or an empirical convergence study into that premise.
`Residual` uses the already reviewed T2 inverse and solution perturbation bounds
for a supplied operator change E, rather than proving a second inverse theorem.
The new perturbed-residual estimate supplies both continuous linear equivalences
A and B and assumes B=A+E. It does not construct B from numerical data. Reviewed
T2 inverse existence has its own completeness premise; no such premise is
silently discharged by the new finite-solve instrument.

For fixed incident data, trace evaluation, modal inversion and open selection
give a=Cz+d. The offset d includes the subtracted incoming trace. Norm estimates
apply to the difference C(z_hat-z), while the flux estimate is based on the
complete amplitude a, including d. Origin phase conventions require the
simultaneously transformed current matrix as in completed F1–F4.

`Current` reuses `S11ScatteringFlux.pair`, `flux`, their full coordinate sum and
interference identity. Amplitude mass is explicitly ∑|a_i|. The observation
bound is the sum of the component-map operator norms. With |J_ij|≤β, the current
bound is β mass(x) mass(y), not a spectral-matrix-norm assertion. This choice
keeps all off-diagonal entries and permits non-Hermitian/indefinite forms and
empty finite spaces. A supplied change H in J receives its own error term.

`Fraction` uses positive absolute-denominator margins and therefore permits
signed denominators and numerators. The physically oriented native incoming
denominator is separately positive. An error smaller than the reference
denominator magnitude gives a nonzero perturbed denominator. A zero margin
does not. Lean's totalized division is not a physical zero-denominator policy.

## Compact translation boundary

The native instrument parses only the balancing/solve/unscaling, residual and
observation/current statements from the existing `construct` AST. It executes
them on a synthetic ten-unknown/two-incident system, with nontrivial scaling,
complex modal mixing, incoming subtraction and off-diagonal current entries.
The complete constructor, boundary classifier, symbolic imports, source operands
and saved physical system are never executed or loaded. Reference inspector
and phase-map files are read-only provenance anchors.


The native fixture checks scaling with K and D, but its numerical inverse-bound
examples use the unscaled M and coefficient sup norm, with the explicitly
known synthetic inverse bounded by 0.625. Its observed-amplitude comparison
uses C_c=P V⁻¹T and c; in balanced variables C_z=C_c D⁻¹. The Lean scaling and
unscaling identities describe the conversion. No numerical bound on the saved
physical K inverse, no norm equivalence constant and no complete physical
residual-to-flux certificate is silently supplied by this split test.

The selected observation AST ends at construction of incoming_flux. The native
positive-incident-current guard and saved outgoingFluxRatio expression immediately
following it are source-inspected only. Synthetic fraction examples execute
separate arithmetic on the tested currents and explicit signed scalars. Thus
no executed native zero-denominator-policy or physical ratio-validation claim
follows from the generic Lean fraction theorem.

The native J selects open outgoing channels and fills only same-end entries;
cross-end blocks are zero. The theorem accepts arbitrary finite J without
certifying that physical channel selection. Existing complex fixtures are
Hermitian; no non-Hermitian or indefinite fixture execution is claimed. The
Lean fullJ witness is singular, so it must not be called positive definite.

For fixed incident modes/current, incoming_flux is independent of the numerical
solution coefficients: residual-only perturbations use denominator error η=0.
Nonzero denominator error is separately supplied uncertainty; the signed scalar
examples exercise that algebra. Pipeline composes the fixed-operator/fixed-J
residual estimate into a fraction. Its perturbed-operator composition stops at
amplitude, and H/unscaling bounds remain separate application steps. No combined
perturbed/H fraction theorem is claimed. Existing denominator_margin can supply
a reference-based margin, but the pipeline still explicitly asks for both bounds.

The independent arithmetic checks compare the residual identity, scaling maps,
affine observation and full quadratic current. Floating-point synthetic solves
use a declared 1e-12 scaled absolute tolerance; they are not exact arithmetic
proofs or an interval-certified inverse. Omission/sign controls must separate
by more than 1e-8. Native run1 passed thirteen identities/examples, fourteen
wrong-formula controls and nine conditional bound examples, recording Python
3.10.12, NumPy 1.26.4 and source/AST hashes. These synthetic checks do not
certify the supplied bounds for any saved physical operator.

Eighteen paired statement rejections and twenty-four positive executions passed,
along with the selected 51-declaration axiom audit. The cross-term and
quadratic-error omissions have distinct scalar values: with a=2,e=1, the
actual change is 5, omission of interference gives 1, and omission of e² gives 4.
Controls use canonical witnesses or arithmetic, not independent derivations or
canonical source replacement. Compiler/resource/import/linter failures do not
count as mathematical rejections. Every accepted false statement reached one
contract_control diagnostic and exactly one unsolved False. No warning, import
error or resource failure was accepted.

All 605 earlier source/evidence files and 243 objects are pinned by
`_measurements/S11_lean_sensitivity_preserved_inputs.json`. F1–F4/P1–P4 closures
are historical and remain unchanged after the user's new resume instruction.
Scientific work stays under the default 2 GiB/no-swap/one-CPU/32-task guard;
the previous Grok resource exception does not apply to proof jobs.

Canonical verification, controls, author source/object/pin validation and two
independent source-fidelity reviews passed. Neither reviewer reran Lean/NumPy or
rehashed the external filesystem; those bindings are author-validated. The user
authorized this fixed-packet transfer. No commit or new increment is authorized.
This conditional finite-model result does not establish physical inverse bounds, full-operator approximation error,
continuum convergence, channel completeness or conservation.
