# E1–E4 independent fidelity review and closure

**The bounded D2 odd-invariant dynamics contract E1–E4 is complete.** Both
independent non-author reviews returned CLEAR, with no blocking findings.
Local builds, axioms, mathematical controls, compact native correspondence and
current source/object hashes passed. This closes only
[DYNAMICS_COVERAGE.md](DYNAMICS_COVERAGE.md), under
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md).

## Fixed reviewed revision

The user explicitly approved transfer of the fixed 34-file packet by answering
“Yes” to the packet-specific request. Both reviewers read the same snapshot in
an isolated transport directory, sequentially, with read-only tools. Neither
reviewer received the other's response or acted as an author of the proofs.

- Packet: `S11 D2 odd-invariant dynamics E1–E4 fidelity contract v1`.
- Aggregate SHA256:
  `16a1dd5784c344f032716d86974b52a58c79110464e2cb1353b915c17d2581c0`.
- Aggregate method: SHA256 of UTF-8 `json.dumps(files, sort_keys=True)`.
- Archive SHA256:
  `d91f74171af7d5fcfbb3bc623d58bc07ccc5de30a5f96a2fd53138789dfacd58`.
- Manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_dynamics_review_packet.json`),
  archive (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_dynamics_review_packet_v1.tar.gz`),
  authorization and review state (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_dynamics_review_state.json`).

| Independent reviewer | Session | Result |
|---|---|---|
| Claude, Opus 5 (run also records Haiku usage) | `88a2fb0e-0d69-4ed3-8882-c738bb29e852` | CLEAR; successful terminal response, exit 0 |
| Grok, `grok-4.6-build` | `7b511736-39b0-4aef-a03b-3a6cf27d4582` | CLEAR; terminal `end_turn`, exit 0 |

Durable reports: [Claude](../../_measurements/S11_lean_dynamics_fidelity_claude_v1.md)
and [Grok](../../_measurements/S11_lean_dynamics_fidelity_grok_v1.md). Each has its
unchanged raw JSON, stderr and run metadata alongside it. Grok's raw text
includes preliminary progress messages; the extracted report starts at its
unique final review heading (character offset 1026). Only the substantive final
report is the verdict, and the entire raw response is retained.

Both are source-fidelity reviews by reading and mathematical analysis. Neither
reviewer independently rebuilt Lean, reran the instruments or recomputed hashes.
The author separately validated the actual snapshot/archive/transport hashes,
live sources, instrument and terminal-report hashes, 31 check logs, 13 compiled
objects, and all 15 clean dependency checkouts against their pins. Reviewer
statements about hash correspondence mean comparison of the supplied records,
not independent cryptographic validation.

Claude stderr was empty. Grok emitted startup warnings about plugin precedence,
hook configuration and repository discovery under `/tmp`. Its source review
completed with the correct packet identifier and substantive final verdict;
these warnings did not truncate or replace the review. No plan-only, cancelled,
errored or partial response was accepted.

## Findings and dispositions

Neither reviewer identified a blocking defect. The following optional points
were considered without expanding the completed contract.

| Point | Disposition |
|---|---|
| Direct Lean mutation of `modalOperator` would make its sign control more explicit. | No extra mutation added. Existing exact action-gradient/PDE identities, native sign/factor mutations, false-mixing controls and witness-sign control cover the load-bearing claims. Both reviewers accepted this coverage. |
| Density mutations also fail `lagrangian_increment` algebraically, beyond the required `density_identity` failure. | Clarified the verification wording. The later change-tactic failures remain secondary; no instrument failure is counted. |
| S10 analytic lemmas are reused, but the action-specific variation chain is repeated for the odd density. | Clarified `DYNAMICS_FIDELITY.md`. No density-generic calculus refactor is required for correctness or this bounded closure. |
| The five Python locus controls evaluate the Lean-side target matrix, not the native routes directly. | Clarified their role. Exact symbolic route identities establish native correspondence; the finite evaluations are target-side instrument controls, not exhaustive coverage evidence. |
| Object hashes depend on the recorded `lake env lean -o` build. | Retained that explicit boundary. Current objects match this run; historical reports were not rewritten, and a later build must validate its own objects. |
| An explicit compact bump would be more constructive than the first-variation existence proof. | Not added. The agreed conclusion is existence; the proved fundamental-lemma argument supplies an admissible compact test. |
| Native period averaging is a `sin²→1/2` rewrite. | Made this tested-translation limit explicit. No new formalized CAS integration claim is made. |
| The native V6 RREF normalization selects coefficient +1 for the odd pairing. | Already checked against actual P_D and the package difference. No additional formalization is needed. |

## Closure evidence and limits

[DYNAMICS_VERIFICATION.txt](DYNAMICS_VERIFICATION.txt) records three new module
builds and the audit root, nine unchanged local dependency builds, all 47
standard-axiom audits, ten mathematical rejections, eight passing controls and
the compact native checks. The proofs establish:

- The exact supplied density `-beta P/2`, its actual first variation and local
  variational derivative on smooth fields with compact test variations.
- All-fields variational nullness exactly at beta=0, and an admissible nonzero
  first-variation witness for every beta≠0.
- The odd contribution's modal action on both longitudinal and transverse
  directions, with nonzero cross-sector action exactly when beta≠0 and k≠0.

The full live packet matched the reviewed snapshot before closure documentation.
Only the four current documentation files (coverage, fidelity, verification and
README) were then updated, and this closure record was added. Proofs, native
sources, dependencies, instruments and terminal check records remain the reviewed
bytes. The approved archive and snapshot remain unchanged. No substantive
re-review, proof rebuild or mutation rerun is needed for these documentation
clarifications. See the
closure validation (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S11_lean_dynamics_closure_validation.json`).

This does not complete all of S11, the full XFORM_EXTRA spectrum, D3–D5 invariant
or divergence classes, interface/radiation physics, production exports or
comparator checks. No S11c calculations or pinned S11c-b/c1/c2 exports were
modified. No commit was requested or made for this closure. The E1–E4 stopping
rule now applies.
