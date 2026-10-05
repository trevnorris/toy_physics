# VC1–VC4 independent fidelity review and closure

The bounded variable-coefficient and flat-interface normal-slice/trace contract
is complete under [FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).
Claude and Grok returned substantive terminal CLEAR reports. Local verification
and current source/dependency/object correspondence pass. This closes the
stated local D3/D4 contract, not the full S11c operator or a multidimensional
transmission theorem.

The user approved the fixed 61-file packet with “Approved.” Both reviewers
received the same isolated read-only snapshot, in separate sequential sessions,
without the other's report.

- Revision: `S11 variable-coefficient and flat-interface VC1–VC4 fidelity contract v1`.
- Aggregate SHA256: `8152df86a31ecaabe89713a3fb178c7c469885f905b9e486207bbc695c02cf45`.
- Archive SHA256: `32e7b1acf2caaf3c329d3f8a31c53b5255d95e4c356a6d72ccefcd66c77ca8f6`.
- Manifest SHA256: `1abf833f0806dcc518836ded7775a7808d73846f84d12b0dbc4462fab8b359bf`.
- Manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_variable_review_packet.json`),
  archive (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_variable_review_packet_v1.tar.gz`),
  authorization/state (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_variable_review_state.json`).

| Reviewer | Session | Result |
|---|---|---|
| Claude `claude-opus-5[1m]` (run also records Haiku usage) | `0cb30c1f-e1ca-4234-9877-f1a0d194e63b` | CLEAR; terminal success, exit 0 |
| Grok `grok-4.6-build` | `d60445f6-38e7-466e-adc3-04cc6a82ad45` | CLEAR; terminal end_turn, exit 0 |

Durable reports: [Claude](../../_measurements/S11_lean_variable_fidelity_claude_v1.md)
and [Grok](../../_measurements/S11_lean_variable_fidelity_grok_v1.md).
Raw JSON, stderr and run records are preserved. Grok's terminal report begins
at character offset 655 of `text`; preceding plans and other provider metadata
are not the substantive terminal report. Both reviewers explicitly performed
source reading and hand algebra, not Lean/native execution or hash audits.
Their clearances are independent fidelity reviews, not independent builds.

Claude stderr is empty. Grok recorded startup plugin/hook/repository warnings;
it nevertheless completed with a substantive report, terminal end_turn and
exit 0. No partial, cancelled, errored or plan-only response was accepted.

Author validation separately checked archive/snapshot/transport and all 61 live
files before closure edits, all forty output objects, all 68 source/log records,
input object and transitive source hashes, generators, fifteen clean package
pins and thirteen direct Mathlib source/object pairs. Historical proof bytes
and manifests and all five T1 objects match. There was no shared-source drift.
See closure validation (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_variable_closure_validation.json`).

| Finding | Disposition |
|---|---|
| Claude: the runner accepts at least one intended False diagnostic, although the report says exactly one. | Clarified the distinction: exactly one was observed in every recorded mutant and independently enforced by author validation; the unchanged runner's acceptance predicate alone is weaker. No additional errors are present in the accepted evidence. |
| Claude: the epsilon expression is not itself a new Lean theorem. | Clarified that the Levi-Civita contraction on G is prose/hand algebra at this layer, with the earlier native D4 check separately recorded. Lean proves the derivative-defined dual matrix and its contraction G:M=2P. No kernel-certified epsilon translation is claimed. |
| Claude: outer run directory differs from the instrument's log directory. | Both are intentional: `_scratch/S11_lean_variable/verification_run2` holds supervisor state/logs; `lean/s11/_scratch/variable_verification` holds per-check source/log files. Documented both. |
| Claude: the generic weighted-current witness, weighted ±1/2 factors and D3 index could have additional Lean mutants. | Optional strengthening. General identities are proved; the selected native checks cover D3/D4 omission and derivative-index sensitivity, and the formal controls retain their disclosed scope. No additional theorem/control added under the stopping rule. |
| Grok: full-vector witnesses, integrated variable-profile variation, weak traces, profile traction or a classification of all vanishing profiles would be stronger results. | All optional and beyond the agreed contract. Lean proves the specified nonzero components, native checks the full vectors, and the local/slice application limits remain explicit. |
| Grok: the name `compact_endpoint_split` is broader than its endpoint-zero hypotheses. | Already documented: the theorem assumes h(a)=h(b)=0, not compact support. No rename or proof change needed. |

Only closure/clarification documentation changed. Definitions, hypotheses,
proofs, instruments, native sources and recorded verification evidence remain
unchanged; no rebuild, control rerun or substantive re-review is needed.
The approved snapshot, archive and transport remain frozen with their original
pending-review wording.

The completed result retains D3 coefficient-gradient terms, D4 beta-gradient
response, weighted-current corrections and the signed interface momentum-flux
jump under explicit regularity assumptions. It makes no claim about every
S11c sector, general weak traces, multidimensional stationarity, surface sources
or distribution products. Constant-coefficient bulk-null densities can retain
nonzero profile responses and boundary momentum.

The user subsequently authorized committing this completed work with a `lean:`
message and then proceeding to question 3. That instruction is recorded in
continuation authorization (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_variable_continuation_authorization.json`)
and supersedes the older cached watcher's no-commit/no-next-scope wording after
these clearances. Question 3 requires its own bounded contract; no future
review-packet transfer is authorized by this approval. S11c calculation files
and pinned exports remain untouched.
