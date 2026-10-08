# O2 record directive: review dispositions (orchestrator)

**Artifact:** `directives/O2_record_directive.md`, the pre-builder decision list for O2 sub-step 7 (scoping §7 row
7). The reviewed version is sha `ddcb5ac0…`, frozen at `_scratch/s9b_build/O2_record_directive_reviewed_v0.md`. It is
orchestrator-written, so it gets Codex + Grok (G2): one pass, findings verified, folded once.

**Legs.** Both used the identical prompt `_scratch/s9b_build/o2_record_directive_review_prompt.md` and reported
before adjudication.
- **Codex (gpt-6.1-sol, xhigh):** NEEDS REVISION, 1 finding (`_scratch/s9b_build/o2_record_directive_review_codex_final.txt`).
- **Grok:** NEEDS REVISION, 1 finding (`_scratch/s9b_build/o2_record_directive_review_grok.txt`). Grok also checked
  the following and found them sound:
  - the cited commits resolve to the named sources;
  - item 2 leaves a printed difference open unless a cited source settles it;
  - the knife-triple sentence matches the comparator's scope.

Each verification is a mechanical lookup. Commands and literal output are in
`O2_record_directive_review_disposition_lookups.md`.

| # | Finding (leg) | Disposition | What must be true after the fold |
|---|---|---|---|
| D1 | The record's source packet includes the comparator build directive. `CLAUDE.md`'s step-record row forbids a build directive in the packet. (Codex 1) | **ACCEPT.** Lookups D1: v0 line 28; `CLAUDE.md:58` ("⛔ no build directive in the packet"). | No build directive is in the record author's packet. The comparator's scope reaches the record through its printed output, its source and its acceptance disposition, which item 3 already requires. The record legs' packet also carries no build directive. |
| D2 | Item 6 has the record state the user's 2026-10-07 premise that the flowing brane's momentum density is `ρ_br V`. The O2 sources leave that map open, so the identification would enter the O2 handoff without support from its sources. It belongs in an S9b premise decision. (Grok 1) | **ACCEPT.** Lookups D2. Spec `4680e251` L133 keeps `ℐ_br^live` OPEN with no `ρ_br V` identification. The engines' acceptance keeps the momentum density OPEN (`O2_build_r3_review_disposition.md:23`). The user's premise is real, but it governs the S9b repair. It is folded into the S9b repair decision list, which gets its own two legs. | Item 6 asks only which O2 objects and conditional inputs Part D may use, and under which conditions. Item 1 keeps every input as its sources state it, so the momentum map stays OPEN in the record. |

**Fold.** D1 and D2 are folded once. The folded version is sha `17fa6706…`; the lookups' "Fold applied" section
shows 0 build-directive mentions and 0 `ρ_br V` mentions. Under G2 there is no second pass. The record author may
start.
