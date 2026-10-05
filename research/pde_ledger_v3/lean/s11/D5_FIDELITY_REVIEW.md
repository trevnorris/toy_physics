# D5.1–D5.4 independent fidelity review

Both independent non-author reviews are **CLEAR**, with no required finding
remaining. The bounded D5 quadratic density classification is complete under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). This closes D5.1–D5.4,
not D5 bulk variation, S11 as a whole or the separate CAS production work.

## Reviewed revision and reviewers

The user authorized transfer of the fixed 64-file packet by replying
"I approve" to the exact-packet request. Authorization, immutable packet and
run identities are recorded in
`S11_lean_d5_review_state.json` (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_d5_review_state.json`).

- Packet: `S11 D5 quadratic invariant D5.1–D5.4 fidelity contract v1`.
- Aggregate SHA256: `ff8aa6e60d1403128615212ccbf601fc0c038ef82d0f64690144d9a9457c081c`.
- Archive SHA256: `e716fd68b6f6f0875925efdb5816e6e23e25dfac7aa31a10b34c56d341622cb5`.
- Claude session: `0cdefacf-bdcc-4371-8840-862f65657d85`, completed
  2026-09-18 01:33:57 UTC. Reported model usage includes
  `claude-opus-5[1m]` and `claude-haiku-4-5-20251001`.
- Grok session: `e31a560b-073b-4a76-bd6f-bcd9e2f5f425`, completed
  2026-09-18 01:43:33 UTC. Reported model: `grok-4.6-build`.

These were separate, sequential source-only sessions given the same isolated
packet. Neither session received the other review. Durable terminal reports:

- [Claude report](../../_measurements/S11_lean_d5_fidelity_claude_v1.md)
- [Grok report](../../_measurements/S11_lean_d5_fidelity_grok_v1.md)

Both inspected the object, domain, group action, coefficient normalization,
classification/census, compact native identification, controls and execution
records. They read the core proofs and sampled the generated constraint and
reconstruction blocks. Neither independently compiled Lean, executed the
native computation, recomputed hashes or exhaustively checked every generated
linear combination by hand. Compilation and preservation remain author-validated
execution evidence. The full certificate is checked by Lean, not by treating
reviewers' samples as a completeness proof.

Raw JSON, run records and stderr are retained alongside the extracted reports.
Claude ended with `success/end_turn/completed`, exit 0 and empty stderr.
Grok ended with `end_turn` and exit 0. Its stderr includes startup warnings and
one refused read request whose line range exceeded the tool's token limit.
The local session trace confirms subsequent smaller reads of the same
`BilinearExpansion.lean` file completed, including its public final tail lemma.
The error was recovered before its substantive terminal report; it is not
accepted as review evidence. Only terminal text after the CLEAR heading is
extracted, excluding earlier commentary and the separate thought field.
Extraction offsets, raw hashes and the error disposition are recorded in the
closure validation.

## Findings and dispositions

Both reviewers confirm that the supplied object is an arbitrary real quadratic
form on all real 5×5 gradients, under full conjugation. The 325-coefficient
representation is exhaustive; necessary equations and rational reconstruction
are kernel checked; trace identities supply full-group sufficiency. They
confirm SO/O equality in odd dimension, the unique three-coefficient
representation, actual subspace dimensions 3/3/0 and the zero odd subspace.
Neither requested a mathematical repair.

| Review note | Disposition |
|---|---|
| Claude: definition-level Lean mutations could strengthen sensitivity. | Retain the agreed thirteen paired statement controls and seventeen positives. Canonical source replacement is not claimed. Definitions were included in independent fidelity review; no broader mutation campaign is needed to close D5.4. |
| Claude: the same-count wrong-span control could exercise the exact RREF comparison path. | Its rank-three wrong space has stacked rank greater than three against the target, which directly demonstrates unequal spans. The real comparison checks both stacked rank and RREF. Current evidence meets D5.3; no claim of exhaustive instrument mutation coverage is made. |
| Claude: transpose the native reflection block to remove reliance on a diagonal reflection. | The actual reflection is diagonal, so its monomial action is diagonal and equals its transpose. The current source matches the proved object. Generalizing native reflection handling is optional CAS maintenance, not a D5 fidelity fix; engine sources remain unchanged. |
| Claude: add an explicit `O_unique` corollary. | O uniqueness already follows from `SO_iff_O` and `SO_unique`. The stated mathematical conclusion is covered without adding a redundant declaration. |
| Grok: general dimension, EL/divergence quotient, D5 bulk, kernel-certified CAS, source mutations or further production/interface/spectral work would be stronger. | These remain exclusions, not unfinished D5.1–D5.4 obligations. |

Reports are retained verbatim. Claude's sentence counting 320 quarter-turn
instances is a prose slip: the actual source contains **321** quarter-turn
instances plus **one** rational-rotation instance, totaling 322 necessary
equations. `ConstraintBlock16` contains the final quarter-turn and the rational
instance. This correction changes no certificate, theorem, test or native
computation. The claim that every retained equation is used is generator/source
coverage evidence; the theorem's completeness rests on checked reconstruction
and full-group sufficiency, not that bookkeeping claim alone.

The dispositions retain the agreed scope. Only closure/status documentation
changed; no substantive re-review, build or mutation rerun is required under
the policy when all reviewed proof and instrument hashes still match.

## Closure evidence and meaning

Before closure edits, author validation matched the archive, snapshot,
transport and all 64 live packet files; all 45 canonical objects; fifteen
clean package pins and five direct Mathlib source/object pairs; 383 historical
proof/evidence files; two read-only native inputs; and 51 protected D4C/NP/T1/VC
objects. No shared-source drift was found. The transitive Mathlib cache remains
the existing pinned baseline, not a fresh rebuild.

Accepted run11 retains 51 standard-axiom audits, thirteen fresh mathematical
rejections and seventeen fresh positives (75 records). Canonical build lineage
is 43 modules built in recorded run9, then Controls/root in run10; all 45 were
reused in run11 under the complete guards. No focused build, timeout, warning
or compiler failure counts as recorded success or mathematical rejection.

Native evidence still identifies the full actual V1/V2/V6 spans and reflection
with 3/3/0 and `P_D=0`. Both wrong-span controls were detected; historical wrong
orientation retains the same counts while failing span agreement. Only selected
native helpers executed; Wolfram evidence is source inspection. This tested
translation remains outside Lean's kernel.

For dimension five, every invariant quadratic density has a unique coefficient
representation in `(tr G)^2`, `tr(G^2)`, `tr(G Gᵀ)`; no nonzero reflection-odd
invariant remains. This is exhaustive certification of the supplied density
family, not a new physical discovery or a statement that three distinct bulk
responses or particular boundary conditions follow.

`S11_lean_d5_closure_validation.json` (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_d5_closure_validation.json`)
records packet/report correspondence, preservation and documentation-only
differences. The approved snapshot, archive, transport and mathematical sources
remain unchanged. Installation/reproducibility fixes retain their separate
completed validation. No commit was requested or performed for this closure.
