# Lean scope, fidelity and completion policy

Effective 2026-09-11, at the user's direction. This is the operating policy for
the v3 PDE ledger's Lean work, including associated generators, verification
instruments and reports. It codifies the user's refocus and the subsequent
discussion of fidelity, coverage and mutation controls. It supersedes older
plans to expand the S10 CAS bridge; those plans remain historical evidence.

Read this policy before starting or resuming Lean work and before expanding its
scope. State the particular obligation being addressed in the working update.
The user can revise the scope; do not infer a revision from an instruction such
as "keep going" or from the existence of more unverified transcript fields.

**Lean's deliverable is a mathematical proof, a compact identification of the
intended object, an exhaustive coverage contract, and mutation controls.**
Statement fidelity must be independently reviewed. Exhaustive transcription of
CAS output is not a completion requirement.

## Division of responsibility

| Layer | Responsibility | Limit of the evidence |
|---|---|---|
| Lean | Deduction from explicit definitions and assumptions; complete mathematical classification; compact fidelity link and coverage contract | A valid proof does not establish that its supplied physics or its translation is the intended one. |
| SymPy and Wolfram | Derive the objects required by the CAS specification, including exceptional cases and controls | Agreement can retain a common modeling or convention error. |
| Comparator and numeric/PIT checks | Check correspondence of the actual outputs, including parameter maps, normalization, domains and deliberately selected strata | Report the actual symbolic or sampled guarantee. Sampling does not prove a global classification or discover every exceptional locus. |
| Independent review | Check statement fidelity, assumptions, retained order, conventions and interpretation | Compilation and declaration counts do not substitute for this review. |

Do not claim that the current comparator or PIT already covers a surface without
checking its evidence. Exact symbolic identities can give stronger assurance
than sampling. Their value does not require checking every emitted expression
in Lean. A bridge creates translation and maintenance obligations; it does not
enlarge Lean's kernel. No layer's unfinished obligation is waived by this policy.

## L1 — Bound the formal proof

Before adding proofs, identify the mathematical conclusion or fidelity gap they
will close. Reuse a general theorem across dimensions, parameter values and
engines wherever it applies. A generated artifact is justified by its contribution
to the declared contract, not by the number of transcript tags it covers.

Do not build or expand a systematic per-output CAS bridge: no default program
to translate every expression, repeated field, solver status, filter trace,
basis presentation, minor or dimensional metadata record into a Lean theorem.
Do not replicate that program across every engine, action package or dimension.

A small exact action/operator identity, a mathematical calculation needed by
the proof, or a targeted certificate needed to close a named fidelity gap is
allowed. Keep the reason and boundary explicit. "Another output is unchecked"
is not such a gap. If an allegedly compact connection starts expanding into a
transcript-wide project, stop the expansion and redesign the connection.

## L2 — Identify the object and prove the coverage contract

The contract must identify:

1. The supplied action or operator, parameter mapping, coordinate conventions,
   units, normalization and mathematical notion of equivalence being used.
2. The domain: dimensions, coefficient signs, nonzero denominators, parameter
   exclusions, regularity and any retained-order or ansatz restrictions.
3. The conclusion and its quantifiers, including the distinction between a
   candidate root and a root with a nonzero mode, algebraic multiplicity and
   kernel dimension, and full versus transverse amplitude spaces where relevant.
4. An exhaustive set of precisely defined cases, with a proof of coverage and
   the required disjointness or explicit treatment of overlaps. State which
   cases can be empty under the chosen dimension or parameter assumptions.
5. The roots, loci, dimensions/counts and other necessary invariants on each
   case, and the corresponding CAS audit obligations.

Matching a few counts or eigenvalues does not identify an operator: different
operators can have those invariants while acting on different polarization
spaces. Retain a compact identification of the actual action/operator, or a
suitable equivalence for the stated claim, in addition to invariant agreement.
For example, multiplication by an everywhere-nonzero scalar preserves kernels;
that fact alone does not preserve a resolvent's normalization or its residues.

The fidelity record must distinguish the identity proved inside Lean from the
tested translation or comparator connection to the actual CAS artifacts. Keep
source/provenance references at this boundary. Do not describe unformalized
translation software as kernel-certified.

"Exactly N strata" means an exhaustive specified classification on a stated
domain, not a claim that there is a unique possible partition. Generic samples
cannot discharge exceptional-case coverage. The CAS audit must deliberately
cover each required case/component and check full subspaces when the claim
concerns a kernel; one residual-zero vector does not establish a complete basis.

## L3 — Keep mutation controls meaningful

Every load-bearing contract claim needs an identified mutation control. Organize
controls around mathematical claims and assumptions; do not require a separate
copy for every helper lemma or repeated presentation of the same claim.

- Use deliberately wrong coefficients, signs, counts, branches, multiplicities,
  denominator assumptions or subspace-completeness claims as appropriate.
- Demonstrate rejection for the intended mathematical reason. A syntax error,
  missing import, timeout or broken test environment does not count.
- Include passing controls. Inspect assumptions for accidental inconsistency
  and include admissible examples or nonemptiness evidence where needed.
- Record the mutation, expected failure, observed result and source provenance.
  Keep canonical proofs unchanged by isolated mutation runs.

Mutation controls test sensitivity and fidelity; they do not replace proof or
guarantee detection of every possible mistranslation. Run the relevant controls
after substantive changes. Do not repeat expensive builds or mutation suites
for documentation-only edits when the verified sources still match.

## L4 — Review fidelity and match effort to the step

For new or materially changed physics-bearing formal claims, obtain two
independent non-author reviews and resolve their findings before declaring
fidelity review complete. Review definitions, hypotheses, parameter maps,
normalizations, excluded cases, quantifiers and interpretation. Record reviewer
identity, reviewed revision and the disposition of findings. A build or mutation
run is not a review leg. Missing review is a named outstanding review obligation;
it is not a reason to keep adding unrelated formalization.

For a closed or calibration result such as S10, keep the deliverable to the
proof, compact fidelity link, coverage contract and mutation controls. For a
novel result, strengthen the invariants and object identification where a
specific claim needs it. Novelty does not license exhaustive transcript bridging.
Cost is not a reason to omit fidelity, exceptional cases or mutation controls.

## Required work contract and stopping rule

Before substantive implementation, write or identify a short per-step contract
in the step's coverage document. Reuse an existing contract rather than creating
a new planning document on every turn. It must record:

| Field | Required content |
|---|---|
| Claim | The named mathematical result and its domain. |
| Fidelity link | The action/operator identity or equivalence, parameter map, and connection to actual CAS artifacts. |
| Coverage | The exhaustive cases and invariant/count obligations. |
| Evidence | Existing theorem names/files, essential remaining proofs, mutation controls and provenance. |
| Review | Independent reviewers, reviewed revision and unresolved findings. |
| Exclusions | What the theorem does not claim and which obligations belong to CAS, comparator/export or physics. |
| Completion | A finite list of deliverables whose completion ends the Lean task. |

Before each proposed increment, name the open contract item it closes. If it
does not close one, do not add it merely because it is possible to formalize.
An expansion must identify a concrete mathematical or fidelity gap and be
recorded as a scope change. Work beyond the user's agreed deliverable requires
agreement on that changed deliverable; work within it proceeds autonomously.

The Lean task is complete when its stated theorem and exhaustive contract are
proved, its compact artifact connection is documented and checked at the stated
level, relevant mutation controls pass, the build/axiom audit passes without
admissions or custom physics axioms, and fidelity reviews are resolved. At that
point, stop adding Lean coverage. Report CAS production or paper work separately.
If a CAS discrepancy invalidates the fidelity link, resolve or explicitly retain
that blocker; do not declare the link established.

Report progress by obligations closed and gaps remaining. Numbers of theorems,
generated lines, transcript fields and axiom audits are verification statistics,
not goals or percentages of mathematical completion.

## Binding S10 disposition

The existing bridge in commit `f21e6459` is retained as evidence. Its full build,
provenance and mutation results remain useful. Preserve it without expanding it
or deleting it as part of this refocus. Reuse the mathematical and fidelity
theorems it already supplies.

The next Lean work is to consolidate the existing six-family mathematical
results into the compact contract above, identify any essential gaps in that
contract, map the relevant mutation controls and complete fidelity review.
In particular, preserve the exhaustive anisotropic directional classification,
root coincidence and full/transverse kernel counts. Reuse the general proofs
for D4 and other covered dimensions.

The former plan to translate the remaining D3 reality-filter traces,
Q5/Q6/period-average metadata, and then all D4/other-package outputs into Lean
is superseded. These are not automatic Lean completion requirements. Production
Q7 construction, CAS stratum handling, the broad comparator/export pipeline and
ledger/paper reconciliation remain separate S10 work with their own acceptance
criteria. This policy does not declare those tasks completed.

An explicit exponentially growing spacetime solution is only a new Lean
obligation if the agreed contract claims that construction. The existing negative
squared-frequency result must retain its narrower stated meaning. Deriving the
supplied action or physical dimension is outside S10's conditional theorem.

Use [s10/COVERAGE.md](s10/COVERAGE.md) for current work. Earlier checkpoint and
bridge reports document the evidence at their recorded revisions; their former
expansion suggestions do not override this policy.
