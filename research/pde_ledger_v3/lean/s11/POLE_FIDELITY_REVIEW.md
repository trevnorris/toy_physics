# NP1–NP4 independent fidelity review and closure

The bounded nonlinear-pencil finite-core contract is complete under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Claude and Grok returned
substantive terminal CLEAR reports with no required mathematical corrections.
Local verification and current source/dependency/object correspondence pass.
This closes the conditional algebra, finite contour identities and explicit
examples in [POLE_COVERAGE.md](POLE_COVERAGE.md), not the analytic existence
theory or a physical S11c pole/scattering calculation.

The user explicitly authorized sending the fixed 31-file packet with
“I authorize sending those to grok and Claude”. Two earlier transfer requests
were rejected before execution; the subsequent explicit authorization resolved
the approval block. Both reviewers received the same isolated read-only packet,
in separate sequential sessions, without the other review's report.

- Revision: `S11 nonlinear-pencil finite-core NP1–NP4 fidelity contract v1`.
- Aggregate SHA256: `9c7e8e41dd1820647c3dc787ad171043944038c70d37b873699d414f49396e57`.
- Archive SHA256: `2c354358720bb6fd03afd8298a237131cf326de785222598ee15bef84b9bb9d7`.
- Manifest SHA256: `a85ac31521dcf002a4baaa4d4a5244c6f98f00da917d470287003cd9f2ec0e6e`.
- Manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_pole_review_packet.json`),
  archive (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_pole_review_packet_v1.tar.gz`),
  authorization/state (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_pole_review_state.json`).

| Reviewer | Session | Result |
|---|---|---|
| Claude `claude-opus-5[1m]` (run also records Haiku usage) | `13c6554e-2dc7-41d7-aa84-735448bc3fc3` | CLEAR; terminal success/completed, end_turn, exit 0 |
| Grok `grok-4.6-build` | `e3941bb7-1d01-49bb-82aa-c0d2016f7547` | CLEAR; terminal end_turn, exit 0 |

Durable reports: [Claude](../../_measurements/S11_lean_pole_fidelity_claude_v1.md)
and [Grok](../../_measurements/S11_lean_pole_fidelity_grok_v1.md).
The raw JSON, stderr and run records are preserved. Grok's substantive report
begins at character offset 848 in `text`; preceding progress text and provider
metadata are not clearance. Claude's report is the full `result` field.
The raw run records retain their original NOT_ADJUDICATED launch disposition;
the final adjudication is recorded separately in the review state and closure
validation rather than rewriting the execution evidence.

Both reviewers performed source reading and hand algebra. Neither executed
Lean, native Python or hash audits; the pinned Mathlib sources and compiled
objects were not supplied in the packet. Their reports are independent fidelity
reviews, not independent builds. Claude stderr is empty. Grok's stderr contains
startup plugin/hook/repository warnings, but it completed with a substantive
report, end_turn and exit 0. No partial, cancelled, errored or plan-only response
was accepted.

Author validation separately checked archive/snapshot/transport and all 31 live
files before closure edits, seven current objects, all 45 source/log records,
input objects and transitive source hashes, fifteen clean package pins and four
direct Mathlib source/object pairs. All seven read-only native/addendum inputs
and forty VC plus five T1 objects match their preservation record. No shared
source drift occurred. The other session's application-mapping commit
`8e8bd430` is preserved and is not treated as a change to the approved packet.
See closure validation (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_pole_closure_validation.json`).

| Finding | Disposition |
|---|---|
| Claude: the singular-pairing control only proves noninjectivity of the zero map; it does not mention `PairingData`. | Clarified in the live fidelity record. General singular-pairing exclusion follows structurally from the supplied `LinearEquiv` and `pairing_eq`; the mutation tests only the concrete zero-map fact. No claim of a structural mutation is retained. |
| Claude: no explicit positive/control instantiates all hypotheses of residue uniqueness/inverse-coefficient identification; basis covariance is not separately stated. Grok also suggests optional uniqueness/identification mutants. | Optional strengthening, deferred under the stopping rule. The general theorems have the explicit reviewed hypotheses and no analytic-existence conclusion. A hand consistency check is `scaledPair`, L0=0, R=(1/2)id, H=0: range V=ker L0 is the whole scalar space and AR=id. This is hand algebra, not a new Lean positive control. |
| Claude: the scalar residue-11 example is not an application of the generic matrix response theorem. Grok suggests separate theorems for truncated values 6 and 5. | Documented the independent scalar expansion and hand identification of the omitted terms. Both generic and scalar identities are proved; no formal instantiation or separately constructed truncated-response theorem is claimed. Additional proofs are optional. |
| Claude: the nonzero-response control tests only the value at z=1. | Clarified that nonzero singular behavior is carried by `square_higher_coefficient` and `jordan_transfer_exact`, together with the retained Jordan coefficient; it is not inferred from the point-value control alone. A direct transfer-weighted-moment mutation would be optional strengthening. |
| Claude: the `Laurent.lean` comment saying B is the second derivative should read conditionally. | Its application-premise meaning is explicit in the already reviewed coverage, fidelity record and prompt, and is reinforced in the live fidelity record. The theorem quantifies over a supplied affine operand only. The canonical comment and proof bytes remain unchanged; no new derivative-identification claim or rebuild is introduced. |
| Grok: general holomorphic-remainder vanishing and analytic-existence theorems would be stronger results. | Explicitly excluded from NP1–NP4. No unspecified remainder is discarded by the finite formulas. The application must supply any additional analytic theorem it needs. |

Two minor phrases in the verbatim reviews do not describe the source literally:
Claude's “Real calculus throughout” means actual differentiation here; the
`Scalar` and `Jordan` derivatives are over C. Grok calls the logarithmic
expansion “eight monomials” while correctly listing four: the logarithmic
product has four terms, and the observation/forcing response has eight.
Both sources and reports' displayed formulas agree on the mathematics.

Only closure and clarification documents changed. Definitions, hypotheses,
proofs, instruments, native inputs and recorded local verification evidence
remain unchanged. No rebuild, control rerun or substantive re-review is needed.
The approved snapshot, archive and transport remain frozen with their original
pending-review wording; the closure validation records the exact live document
differences. The requested `POLE_HANDOFF.md` remains outside the fixed packet.

What this adds to the ledger is a reviewed machine-checked separation of the
objects: R need not be a projection; the full-pairing products RA and AR are
conditional projections; nonlinear logarithmic contour values can fail
idempotency; an affine Jordan contour can still give I; and a physical transfer
can have zero residue with a retained higher pole. Frequency-dependent forcing
and observation can change the actual residue from the frozen value 0 to 11.
These are established mathematical distinctions, not new physical discoveries.

The user authorized the checkpoint commit with “After they approve go ahead
and commit”; both required clearances and finding dispositions now satisfy
that condition. See commit authorization (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_pole_commit_authorization.json`).
This supersedes the older cached watcher text saying no NP commit was requested.
The commit uses a `lean:` prefix and preserves unrelated/other-session work.
No new mathematical increment is authorized by this completion.
