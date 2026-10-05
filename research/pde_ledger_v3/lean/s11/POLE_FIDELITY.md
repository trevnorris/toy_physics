# NP1–NP4 source and statement fidelity

Governed by [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md) and
[POLE_COVERAGE.md](POLE_COVERAGE.md). Local verification passes and Claude and
Grok independently returned CLEAR. The bounded contract is complete; review
provenance, limits and optional findings are recorded in
[POLE_FIDELITY_REVIEW.md](POLE_FIDELITY_REVIEW.md).
The mathematical/application assessment is [POLE_ASSESSMENT.md](POLE_ASSESSMENT.md).

## Authority and compact object identification

The calculation session already adopted `nonlinearPoleV2` in
`directives/S11c_d_NONLINEAR_POLE_CONTRACT.md`. Its original shared v10 source
and all production inputs remain unchanged. This increment does not discover
new physical poles or rerun the synthetic repair or production scripts.

`_measurements/S11_lean_pole_source_check.py` parses selected assignments in the
original diagnostic sources without importing or executing them. It identifies
L=z², L=z²−1, L=zI−N with N=[[0,1],[0,0]], the second-coordinate forcing,
first-coordinate observation, and affine maps B=1+3z, O=2+5z in their existing
records. Four controls reject a same-count wrong Jordan sign, a same-kernel
rescaling, a different forcing coordinate, and a frozen observation. Selected
existing native identity/rejection evidence is read and marked historical.
This compact translation is tested and independently reviewed. It is not
kernel-certified; the reviewers read sources and evidence without executing
the instrument or recomputing hashes themselves.
No claim is made that the old diagnostic's complete 89-check computation was
rerun by this instrument.

In the basic modal algebra, A and W are supplied linear maps; the name
`derivative` for A does not prove that it is the derivative of a pencil. That
identification belongs to the application. W is required to annihilate L0 only
in the separate inverse-coefficient identification theorem. Its factors retain
their recorded order; no frequency-dependent conjugation is introduced.
Coordinates and supplied operators are dimensionless. The complex coordinate
is z and every circle is centered at zero, including the two-root example. The circle
integral has positive orientation and normalization 1/(2πi). Positive radii
exclude the zero-centered poles from the path. The z²−1 example uses radius 2,
so both roots are enclosed and neither is on the contour. Field inversion is
totalized at zero in Lean; inverse identities require z≠0 and the contour never
uses the singular center as an inverse value.

## What Lean proves

| Module | Mathematical statement and boundary |
|---|---|
| `Modal` | Over a field, independent vector-space types X,Y,K and maps A:X→Y,V:K→X,W:Y→K with a supplied linear equivalence D=WAV. R=VD⁻¹W has RAR=R; RA and AR are idempotent with their full ranges and modal ranks. A separate theorem identifies a supplied inverse coefficient using range V=ker L0, W L0=0, L0 R=0 and L0 H+A R=I. These are coefficient identities from a supplied simple inverse expansion; existence of such an expansion is not proved. |
| `Moments` | On a complete complex normed vector space, actual normalized circle integrals of every finite sum of integer-power monomials extract exactly exponent −1. Polynomial terms have zero moment by integration, not by definition. |
| `Laurent` | Actual integrals of finite double-principal-part operands. The ordered product with the supplied affine operand A+zB retains C1 A+C2 B. Interpreting that operand as an actual pencil derivative is an application premise. The response with rectangular observation and forcing retains O0 C1 B0+O0 C2 B1+O1 C2 B0. No unspecified holomorphic remainder is discarded. |
| `Scalar` | Actual scalar derivatives and inverses, circle integrals for z² and z²−1, nonidempotency, the nonzero higher coefficient of z⁻², and the actual response residue 11. |
| `Jordan` | Two-sided inverse of zI−N for z≠0, determinant z², the one-axis kernel at zero and explicit Jordan chain, actual inverse/logarithmic contour I and nonzero higher coefficient N, idempotency and trace 2, and the faithful physical transfer z⁻² with zero residue. |
| `Controls` | Admissible scaled and full two-dimensional modal pairings, a simple zero-map noninjectivity witness, noncommuting product-order witness and frozen-map contrast. The exclusion of singular pairings from `PairingData` is structural and is not itself tested by the zero-map control. |

Matrices use Mathlib's finite elementwise norm for their vector-valued circle
integrals. This is an explicit choice in the finite examples, not a physical
outgoing-operator norm or domain assertion. The Jordan state space is C²; its
physical input and output are scalar coordinates selected from the actual
proved inverse. The modal algebra keeps X and Y separate even when a numerical
example has equal dimensions.

The abstract modal theorem permits a zero-dimensional modal space. It does not
assert that a pole exists; a nonzero mode/residue or an actual singular inverse
is an additional existence obligation.

The coverage is universal for the modal algebra under its hypotheses and for
the finite Laurent formulas. The examples discriminate the named cases; they
are not an exhaustive classification of arbitrary nonlinear pencils. In
particular, the finite core does not prove the analytic Fredholm/Keldysh
existence theorem, a general argument principle, a general Riesz calculus,
all root chains/partial multiplicities, or a physical S11c realization.
No formal scattering, bound-pole existence, normalizability or channel claim
follows from these synthetic examples.

## Verification and controls

The recorded suite `_measurements/S11_lean_pole_contract_check.py` passes in
`_scratch/S11_lean_poles/verification_run3`: six modules plus an audit root,
55 standard-axiom audits, 17 paired mathematical rejections and 21 positives,
45 check records and seven output objects. The seven canonical builds were
fresh in recorded run1 and reused in run3 under the full guards; all 38 controls
ran fresh. The compact native check passes 19 checks and four translation
controls. These are verification statistics, not a percentage of general pole
theory. [POLE_VERIFICATION.txt](POLE_VERIFICATION.txt) records the final hashes.

Each paired false statement yielded exactly one diagnostic in
`contract_control`, containing exactly one unsolved `False`. Syntax, imports,
warnings, timeouts, memory limits and other instrument failures are not
mathematical rejection. The scaled pairing, full rank, product order, contour
sign, higher coefficient, nonlinear cluster, Jordan residue and response-map
omissions have explicit controls. There are no canonical source-replacement
mutations in this increment.

The singular control proves noninjectivity of multiplication by zero on C;
it does not mutate `PairingData`. A supplied `LinearEquiv` and `pairing_eq`
exclude singular pairings in the general theorem. The zero-residue/response
control checks the nonzero value at z=1; the retained singular coefficient is
proved by `square_higher_coefficient` together with `jordan_transfer_exact`,
and is not inferred from that point value alone. The scalar response-residue
proof uses its own exact expansion, not an instantiation of the generic matrix
theorem. The rejected targets 6 and 5 match the omitted terms by hand algebra;
there are no separate Lean theorems constructing those truncated responses.

Focused builds are development diagnostics, not recorded reuse evidence.
The full suite uses one worker, -j1 -M4096 and 600 seconds per process with the
previously tested whole-process-group timeout cleanup. Reuse requires identical
transitive local sources, package pins, commands and input/output objects.
Historical reports remain historical; all forty VC and five T1 objects are
protected and are not rebuilt by this increment.

The initial read-only preflight is
`_scratch/S11_lean_poles/initial_preflight.json`, copied unchanged to durable
`_measurements/S11_lean_pole_preserved_inputs.json`. The object-preservation
guard is a workspace safeguard, not an extra mathematical premise.
`_measurements/S11_lean_pole_validation.json` records author validation of
all live local sources, formal/native instruments and reports, dependency pins,
direct Mathlib source/object pairs, input/output objects and preserved inputs.
It also checks the literal saved controls and their diagnostics against the
current instrument, and guarded build records against both archived prior
runs. The existing transitive Mathlib cache is the pinned trust baseline;
its entire dependency graph was not rebuilt by this increment. The user
explicitly approved the fixed packet's transfer, and two independent non-author
reviews now clear it. Author validation separately confirms archive, snapshot,
transport and live source correspondence before closure documentation edits.

Run1 compiled all seven objects and printed 55 standard-axiom lists, then stopped
in the audit parser: it removed spaces but not embedded newlines in a wrapped
list. No controls ran. The bounded instrument repair strips all token whitespace;
a regression check rejects nonstandard axioms and accepts wrapped/empty lists.
Paired proof code now uses one normalization call, avoiding an unreachable
follow-up when a true rewrite already closes its goal. Statements and values
are unchanged. Canonical proofs, pins, native evidence and all seven NP objects
matched run1, and all 45 historical objects remained unchanged. Run2 used guarded
canonical reuse and restarted the control suite.

Run2 reused all seven guarded canonical objects, passed all 55 audits and four
paired positives/mathematical rejections, then stopped at the noncommuting
product-order positive control. Broad simplification erased the zero summand
before the named identity could match. Only that pair now uses `norm_num only`
with the named identity; statements and values are unchanged. Run3 passed all
controls after the strict paired probe, with the same guarded canonical reuse.
No failed instrument or positive control is counted as a rejected mutation.

The axiom-parser regression, ordered-control probes and repair preflights remain
in `_scratch/S11_lean_poles/`. The final author validation confirms every
canonical proof, pin, native input and object retained its recorded bytes.
The timeout cleanup is unchanged from the previously tested process-group
implementation. No further proof build is needed for these documentation-only
status and scope clarifications. Both reviewers assessed statement fidelity
and disclosed their verification limits; their source reviews are independent
review legs, while the build/hash validation remains author evidence.

The unchanged `Laurent.lean` docstring about B as a second derivative is read
conditionally: if the supplied affine operand is the actual pencil derivative,
B has that interpretation. The formal theorem only quantifies over four
matrices and an explicit finite operand. No derivative identification was added
by the comment. This interpretation is explicit in the reviewed coverage,
fidelity document and prompt; the comment was left unchanged to preserve the
reviewed proof bytes. General analytic existence/remainder theorems and extra
control families suggested by the reviewers are optional stronger results,
not outstanding obligations of NP1–NP4.
