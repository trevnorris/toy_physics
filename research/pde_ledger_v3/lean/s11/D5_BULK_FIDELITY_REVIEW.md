# D5B.1–D5B.4 independent fidelity review

Both independent non-author reviews are **CLEAR**, with no required correction
remaining. The bounded full D5 constant-coefficient bulk contract is complete
under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). This does not close
S11 as a whole or the separate CAS production and inhomogeneous-operator work.

## Fixed revision, authorization and reviewers

The user explicitly authorized the current Claude/Grok review and a commit
after both sign off: “If you're going to ask for approval for grok and Claude
to review, you have it. Wrap this up and if they sign off then commit.”
That advance authorization was bound to the frozen verified 80-file packet.

- Revision: S11 D5 constant-coefficient bulk D5B.1–D5B.4 fidelity contract v1.
- Packet SHA256: `d1c84f4d2ebb5d5bfc2933f41ec2c70e17e6059b0cfcd74dfe580beb92759b5e`.
- Archive SHA256: `f4cef6ae169e0c34c07bdacdf8b496ee40180139c2285fc74189818e5f055dde`.
- Claude session: `8bd7f971-28dc-410a-9dd6-8b9c2853da0b`, completed
  2026-09-18 04:21:58 UTC. Reported model usage includes
  `claude-opus-5[1m]` and `claude-haiku-4-5-20251001`; no subagents spawned.
- Grok session: `9cfd3052-ca2c-4b13-9a51-151f14e9804a`, completed
  2026-09-18 04:30:27 UTC. Reported model `grok-4.6-build`.

The sessions were separate and sequential, with read-only tools and the same
isolated copy of the packet outside the repository. Neither received the other
review. Durable substantive terminal reports:

- [Claude report](../../_measurements/S11_lean_d5_bulk_fidelity_claude_v1.md)
- [Grok report](../../_measurements/S11_lean_d5_bulk_fidelity_grok_v1.md)

Both read source and supplied execution records. Neither rebuilt Lean, reran
native computations or controls, or independently computed file hashes. Their
clearance is statement-fidelity review; execution and preservation are separately
author-validated. Raw JSON, launch/run records and stderr are retained. Original
run records preserve the pre-adjudication `NOT_ADJUDICATED` label; the review
validation, closure record and review state supply the subsequent CLEAR verdicts.

Claude ended `success/end_turn/completed`, exit zero, with empty stderr. Grok
ended `end_turn`, exit zero, with a substantive completed report. Its stderr
contains startup configuration/repository warnings and one read error for the
optional duplicate manifest path. The session trace confirms a successful read
of the authoritative `MANIFEST.json`; the prompt and final report explicitly
use it. That failed lookup is not review evidence and did not leave a required
proof file unread. Grok's durable report starts at its final review heading;
earlier progress text and its separate thought field are not extracted.
Raw hashes, extraction offsets and error disposition are recorded in
`_measurements/S11_lean_d5_bulk_review_validation.json`.

## Findings and dispositions

Both reviewers checked the three-coefficient/five-component distinction,
full gradient and derivative-row convention, actual `L=-Q/2` and jet momentum,
smooth-background/compact-variation hypotheses, actual relative-action derivative,
universal EL/equivalence/nullness statements, response/kernel dimensions 2/1,
explicit current and normalization, native basis/sign/factor map, controls and
scope limits. Neither required a mathematical or implementation repair.

| Optional note | Disposition |
|---|---|
| Claude, echoed by Grok: native `current_factor_rejected` compares `2X` with `X` directly. | Accepted as a limitation. This is a generic polynomial sanity check, not an independent mutation of the computed current. The nine target checks comprise two wrong computed-formula rejections, this one sanity check and six positive witnesses. Actual native `current_divergence`, the Lean boundary identity, nonzero current/density witnesses and the Lean current-factor pair supply the meaningful current-normalization evidence. Reports retain their original bytes; closure prose now states the distinction. |
| Claude: the homogeneous parameter map can be misread as density equality. | Explicitly identify `rho=0,mu=c,B=a+b+c` as a bulk/modal operator map. No pointwise action-density equality is claimed. Divergence terms and their possible boundary effects remain. |
| Claude and Grok: the nonzero-first-variation witness is nonconstructive. | `exists_nonzero_firstVariation` is a classical existence theorem obtained from non-nullness, not a formula for an explicit test field. This satisfies its stated claim; a constructive witness is outside the bounded completion requirement. |

The optional notes need only these fidelity clarifications. No reviewed proof,
definition, theorem statement, generator, instrument, native source or test
changed; no substantive re-review or verification rerun is required by policy.

## Closure evidence

Before closure edits, author checks matched all eighty live packet files,
snapshot, archive and transport; sixty isolated current objects; fifteen clean
package pins and thirteen direct Mathlib source/object pairs; all source/log/
build-lineage records; 499 historical source/evidence files, 232 shared local
objects, two native inputs and the installation records. No shared-source drift
was found. The transitive Mathlib cache remains the existing pinned baseline.

Recorded run4 has 92 records: sixty canonical object records, 54 standard-axiom
audits, fourteen paired mathematical rejections and eighteen fresh positives
(seventeen distinct positive statements). Forty-two classification objects
come from the fully validated fresh portable D5 replay. Twelve other imports
and Action/Variation were built in run1; Bulk in run2; Census/Controls/root in
run3. Run4 reuses them under full guards. The 54 root audit lists are revalidated
run3 output; only the 32 controls execute freshly in run4. Focused builds,
resource failures, warnings and failed positives never count as mathematical
rejections or recorded reuse evidence.

Native evidence retains twelve exact identities and two mutations of original
EL helper implementations, plus the nine target checks distinguished above.
The full three-dimensional basis and actual momenta match; native EL uses the
opposite sign and V5(Q) equals twice Lean EL for L=-Q/2. Only selected helpers
executed; Wolfram evidence is source inspection. Translation stays outside the
Lean kernel.

The separately authorized portable D5 density integration is also complete:
ten tooling regressions, 45 fresh objects, 51 axiom audits, thirteen rejections,
seventeen positives, generator and native checks. `INSTALL_D5_VALIDATION.json`
remains the immutable record of that execution. After review, only installation
status prose changes to say D5B is now complete but has not been added to the
portable catalog. Its earlier document hashes remain historical; the closure
record lists these precise documentation changes. Runner/setup/test sources,
pins and all mathematical evidence are unchanged. No full `all` replay or
fresh OS bootstrap is claimed.

## Meaning and stopping rule

Three classified D5 densities have two independent constant-coefficient bulk
responses and one null direction `(t,-t,0)`. The null density is a divergence
and can retain nonzero momentum and boundary effects. This is exhaustive
verification of the specified family, not a new physical prediction or the
absence of interface effects. No variable-coefficient/general-dimension/null-
Lagrangian expansion, spectrum/scattering/pole work, production rerun or
systematic CAS bridge was added. The stopping rule for D5B.1–D5B.4 now applies.

Closure: `_measurements/S11_lean_d5_bulk_closure_validation.json`.
Review state and authorization: `_measurements/S11_lean_d5_bulk_review_state.json`
and `_measurements/S11_lean_d5_bulk_transfer_authorization.json`.
