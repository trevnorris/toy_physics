# Next S11 contract: invariant-space fidelity

Authorized after checkpoint `d80d02d2`, 2026-09-16 UTC: the user said
“Go ahead and start” to this D2 I1–I4 contract. **I1–I4 complete; both independent fidelity reviews CLEAR.** The completed H1–H4 contract
remains in [COVERAGE.md](COVERAGE.md). Its exclusions explicitly reserved the
invariant census for a separately agreed contract.

## Why this is the next useful gap

Q9 of [the shared specification](../../directives/S11_SHARED_PHYSICS.md)
defines quadratic forms in the entries of a real matrix G that are invariant
under every conjugation `G -> R G Rᵀ`, for either SO(D) or O(D). Counts alone
cannot identify this space. This is also where the historical claim of three
SO(D) invariants in every dimension failed.

Read-only reconnaissance found a more concrete fidelity discrepancy in the
checkpoint SymPy source. `compute_q9` (lines 497–515) places the images of basis
monomials in **rows**, then takes the right nullspace of those rows. A column of
polynomial coefficients instead needs the transpose of each generator block.
The Wolfram constructor explicitly uses `Transpose[actionRows]` (line 910).
This is a source-level distinction, not a basis-presentation difference.

The [small D2 probe](../../_measurements/S11_lean_q9_orientation_probe.py) and
[its evidence](../../_measurements/S11_lean_q9_orientation_probe.json) execute
five original helpers with the original real coordinate names. They do not run
the production audit or initialize the registry. One emitted basis polynomial is

`p(G) = G11² + G12 G21 + G22²`.

For `G=diag(1,0)` and the proper rotation
`R=[[3/5,-4/5],[4/5,3/5]]`, it gives `p(G)=1` but
`p(R G Rᵀ)=481/625`. Thus the emitted space contains a non-invariant polynomial.
It also omits `(tr G)²`: adding that invariant raises its rank from 4 to 5.
Nevertheless, the emitted SO/O dimensions are **4/3**, the expected counts.
The transposed D2 control has dimension 4 and passes this rotation test. That
positive test is not a proof of full group invariance.

This was not found in the searched defect-register prose. It is not yet a
survey of D3–D5 or downstream outputs. H1–H4 explicitly excluded Q9 and checked
that the selected MAIN action is independent of its input; that closure stands.

## Authorized bounded increment: D2

| Item | Deliverable |
|---|---|
| I1 — Object and coefficient action | Define actual real quadratic forms on 2×2 matrices and invariance under all SO(2)/O(2) conjugations. Prove the coefficient-action orientation; do not substitute a few sampled rotations for the quantified definition. |
| I2 — Complete spaces | Prove spanning and independence: SO(2) dimension 4, O(2) dimension 3, and a one-dimensional reflection-odd complement. One possible basis uses `(tr G)²`, `(G12-G21)²`, `(G11-G22)²+(G12+G21)²`, and `(tr G)(G12-G21)`; the last is reflection-odd. The basis must be proved complete, not supplied as an assumption. |
| I3 — Compact native identification | Identify the actual polynomial convention and compare the full invariant spans, including the explicit counterexample above. Any minimal native Q9 correction must be named and reviewed; an engine disagreement cannot be hidden behind matching counts. Broader production reruns and downstream export reconciliation remain separate work. |
| I4 — Controls and closure | Reject the wrong coefficient-action orientation, omission of the odd pairing and the false SO(2) count 3; include admissible positive controls. Sequential builds and axiom audits, compact source evidence and two independent non-author fidelity reviews, then stop. |

This increment deliberately starts with the smallest exceptional dimension,
where the concrete source discrepancy can be proved without a large group-theory
or transcript-translation project. D3–D5 completeness, Euler–Lagrange/total-
divergence classification, S11b and S11c are not part of this increment.
Use one Lean worker with the existing memory settings. The fixed 28-file review transfer was explicitly authorized; both independent
reviews returned CLEAR. New scope or a new transfer packet retains its own
authorization requirements.

I1–I4 are now agreed scope. Stop once the stated proofs, compact identification,
controls and two independent fidelity reviews are complete. H1–H4 remains a
separate completed contract.


## Implementation and verification handoff

The bounded implementation is in
[Quadratic.lean](S11Invariants/Quadratic.lean),
[Rotation.lean](S11Invariants/Rotation.lean),
[Classification.lean](S11Invariants/Classification.lean),
[Census.lean](S11Invariants/Census.lean), and
[Controls.lean](S11Invariants/Controls.lean), with
[S11Invariants.lean](S11Invariants.lean) as the axiom-audit root.
The object conventions, exact span connection and production boundary are in
[INVARIANT_FIDELITY.md](INVARIANT_FIDELITY.md).

The native correction is limited to transposing each generator's monomial-image
block before stacking in `compute_q9`. The compact source check passes for the
D2 SO/O/odd spans and rejects the wrong-orientation control despite matching
counts. Wolfram remains unchanged. These are local source checks; production
results have not been regenerated or cleared.

The recorded one-worker suite is
[`S11_lean_invariant_contract_check.py`](../../_measurements/S11_lean_invariant_contract_check.py).
It builds the five modules and audit root, checks standard axioms, tests nine
mathematical mutations and eight positive controls, and validates source and
dependency hashes. Reuse requires matching commands, source/dependency hashes
and compiled-module hashes; mutations always run afresh. Tool, syntax,
resource and import failures do not count as mutation rejection.

The recorded suite passed: five modules and audit root, 46 standard-axiom
audits, nine mathematical rejections and eight positive controls. All live
source, output and dependency checks pass; see
[INVARIANT_VERIFICATION.txt](INVARIANT_VERIFICATION.txt).

Completion: both independent non-author reviews returned CLEAR with no
blocking findings. All local verification and live correspondence checks pass.
See [INVARIANT_FIDELITY_REVIEW.md](INVARIANT_FIDELITY_REVIEW.md) for the fixed
revision, reviewer identities and optional-finding dispositions. I1–I4 are
complete; the stopping rule applies. Downstream production work remains
separate, and no commit is made by this closure.
