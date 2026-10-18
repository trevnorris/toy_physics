# The 64 GB run plan in DEFERRED_HEAVY_RUNS.md: review dispositions (orchestrator)

The mechanical lookups behind each disposition are in `DEFERRED_HEAVY_RUNS_plan_lookups.md` beside this file
(generator `_scratch/backlog_plan/gen_lookups.sh`).

## The round
- **Artifact.** The section "Run plan for the 64 GB box (user-approved 2026-10-18)" in
  `research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md`. The orchestrator wrote it, and the user approved the order on
  2026-10-18 and asked for it to be committed before the RAM upgrade.
- **Role.** The plan orders runs that the file's entries already define. It does not change what any run computes,
  any premise, any check, or what may be claimed, so it is not physics-bearing. The two-report commit gate still
  applies, so it got one pass of two legs, with the findings folded once.
- **Prompt and legs.** The prompt is `_scratch/backlog_plan/review_prompt.md`. The legs were Codex (`gpt-6.1-sol`,
  xhigh, read-only sandbox) and Grok (`grok-4.7`, effort high), as the authorship table gives for an
  orchestrator-written artifact. Their reports are `_scratch/backlog_plan/review_codex.txt` and `…/review_grok.txt`.
  Both returned **NEEDS REVISION**.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Codex 1, Grok 2 | **Wave 2 can exceed the guard's ceiling.** If the Wave 1 SymPy job is still running, its 30 GiB plus Wave 2's 24 and 16 GiB is 70 GiB. The guard refuses a launch when the reserved caps plus the new cap pass its ceiling (guard lines 251–253). | **Confirmed.** Folded: a new pre-launch line says the guard reserves each cap until the job exits and that a refused launch waits. Wave 2 now says only one of its two jobs fits while the SymPy job holds 30 GiB. |
| 2 | Codex 2 | **The comparator run does not settle the c1 entry's open reconciliations.** The comparator compares `DTN_OPERATOR` by its outer operator signature (line 607) and defers exact canonicalization of oversized leaves (line 1948). | **Confirmed.** Folded: the c1 item now calls the run one measurement and lists the entry's open reconciliations as further work, which may need a comparator change that Codex builds and two legs review. |
| 3 | Codex 3 | **"After each run" put the comparison before the review legs.** G3 says to launch legs before inspecting the result (`CLAUDE.md:161`). | **Confirmed.** Folded: check only that the run finished and its output exists, launch the legs, record the comparison after both report, then commit. |
| 4 | Codex 4 | **"A job killed at its cap gives its real peak" is false.** The cap is both the cgroup memory limit and the address-space limit (guard lines 369–370), so a killed job shows only that it needed more than the cap. | **Confirmed.** Folded: a job killed at its cap needed more than the cap, its peak is still unknown, and it reruns with a larger cap once admitted. |
| 5 | Grok 1 | **The S9b item's 13 GiB caps are below SymPy's recorded round-9 peak, and Wave 1 should name the builders' reruns at 18 and 8 GiB.** | **Partly confirmed.** The 13 GiB caps are for the two review legs, which are told not to rerun the full SymPy harness and to retry a killed SymPy job at 18 GiB (review prompt lines 122 and 132). The builders' runs come before the review and are not Wave 1 jobs. What is right: the plan did not say the review waits for the SymPy builder, or that a leg may retry at 18 GiB. Folded both, with the resume brief for the case where the shutdown interrupts the builder. |

## Outcome
All five findings are folded once; nothing is outstanding that changes what a run computes or what may be claimed.
The plan is committed as folded.
