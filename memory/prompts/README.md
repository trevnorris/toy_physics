# Semantic extraction prompt suite

These prompts turn a prepared memory transaction into staged Markdown. They
are instructions for a semantic writer; they do not discover sources, read the
working tree, publish pages, or update synchronization state.

Use `00-snapshot-contract.md` with every task, then select one task prompt:

| Task | Prompt |
|---|---|
| Monolithic or componentized paper capsule | `source-capsule-paper.md` |
| Research charter, plan, defect register, or other governance capsule | `source-capsule-governance.md` |
| Step record with engines, comparators, and outputs | `source-capsule-step.md` |
| Software project, verification project, or result-family capsule | `source-capsule-software.md` |
| Deleted/retired source capsule | `source-capsule-lifecycle.md` |
| Cross-source topic page | `topic-synthesis.md` |
| Domain script catalog | `script-catalog.md` |
| Candidate-conflict refresh | `conflict-register.md` |
| Navigation/index refresh | `index-refresh.md` |

The caller must supply the transaction manifest path, task record or unit ID,
staged output path, target commit, refresh date, and any staged base page whose
human-owned text must be preserved. A prompt must stop with a structured
missing-input report when one of those inputs is absent; it must not inspect
live repository sources to repair the packet.

The normal order is source capsules, topics and script catalogs, conflicts,
then the index. All outputs remain `memory_review: ai_draft` until a separate
recorded review checks the exact wording and citations.
