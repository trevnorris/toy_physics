# Guarded-run launcher, shared machine-wide pool: review dispositions (orchestrator)

The change is orchestrator-written, so the legs are Codex and Grok (`CLAUDE.md` G1). The launcher constrains every
CAS run, but it does not change what any run computes or what may be claimed from it. It is reviewed until clear
anyway, because a defect here can let jobs overcommit the machine.

**The change, round 0.** The user (2026-10-15) said the machine is dedicated to this work: use every resource, with
no host-memory padding. In `scripts/s11c_guarded_run.py`, round 0 made four changes:
- one machine-wide reservation registry (`~/.local/state/s11c-guard/resource-pool`) in place of one per checkout;
- a pooled budget equal to `MemTotal`, rounded down to whole GiB, in place of 16 GiB;
- no host-memory reserve (minimum 0 in place of 4 GiB);
- all logical CPUs, physical cores first.

`AGENTS.md` was changed to match. Packet: `_scratch/s9b_build/guard_share_review/`:
- before: `s11c_guarded_run_before.py`, sha256 `6a636942…`, the committed guard;
- after: `r0_s11c_guarded_run_after.py`, sha256 `984c4256…`;
- `r0_guard.diff` (`17f7b844…`) and `r0_AGENTS.diff` (`a909f669…`);
- `r0_baseline.sha256`, written before the files were renamed with the `r0_` prefix;
- prompt `r0_review_prompt.md`, sha256 `98973aab…`.

**Legs.**
- Codex (gpt-6.1-sol, xhigh), NEEDS REVISION, seven findings. Report: `r0_review_codex_final.txt`; evidence:
  `r0_codex_evidence/` (scripts 01–04 and their stdout).
- Grok, NEEDS REVISION, three findings. Report: `r0_review_grok.txt`; evidence: `r0_grok_evidence/exercise.{py,stdout}`.

Both found the four intended changes implemented. Both found these refusals intact: a live unaccounted unit, an
inactive reservation whose owner is dead, a malformed or over-budget record, a duplicate CPU, a live limit mismatch,
and a release while the unit is active. Both also found the startup limit verification intact. No leg started a unit
or touched the live registry.

| # | Finding (leg) | Disposition |
|---|---|---|
| X1 | With the reserve at 0, admission credits a live unit's reclaimable file cache twice: once as `memory.current` subtracted from its cap, and again inside `MemAvailable`. Codex's model: 16 GiB machine; a live 8 GiB cap holding 2 GiB anonymous and 4 GiB cache; 6 GiB outside the pool; 4 GiB free. A new 6 GiB cap is admitted, though filling both caps needs 20 GiB. The committed guard refused this only because of its 4 GiB reserve. (Codex 1; Grok saw the same arithmetic and judged it within `MemTotal`) | **ACCEPT.** Without a reserve, admission has to count each byte once. Otherwise jobs plus processes outside the pool can exceed RAM, and the kernel's OOM killer may then pick an agent session. Required: a live unit's credited use is only memory the kernel cannot reclaim from it. Its file-LRU pages are already counted in `MemAvailable`. |
| X2 | Switch-over: a reservation written by an old-guard copy is invisible to the new registry until its unit is live, and the reverse holds too. Two jobs can then share a CPU and exceed the budget. (Codex 2, Grok 2) | **ACCEPT, procedure.** No code can see every checkout's old registry. Required: install the new guard in every checkout at a moment when no guarded unit is live and no launch is in flight. Then confirm that every old per-checkout registry is empty. Install record: below. |
| X3 | Legacy and pooled launches from different checkouts can both pass. The legacy overlap check runs under the per-checkout execution lock, not the shared admission lock. A scheduled interleaving admits both on CPU 7, with caps summing past `MemTotal`. (Codex 3, Grok 1; present in the committed guard too) | **ACCEPT.** Required: a legacy launch cannot interleave with or overlap a pooled admission from any checkout. |
| X4 | A record whose CPU is the string `"7"` validates, because limits are compared through `str(cpu)`. It does not occupy integer CPU 7, so a second job is admitted on CPU 7. (Codex 4; present before) | **ACCEPT.** Required: a reservation's CPU is an integer, or the record is refused as malformed. |
| X5 | A reservation whose unit is live but whose launcher has died is admitted beside. (Codex 5; present before) | **NO CHANGE.** The live unit's cap is still counted, and its actual limits are verified against its record at every admission, so admission beside it cannot overcommit. Once it ends, the dead owner makes the record stale, and admission refuses until it is reconciled. My prompt overstated the property as "a dead owner refuses admission". |
| X6 | `release_pool` compares the stored owner with the owner field of the spec it is given, never with the calling process. A second process carrying the spec released a pending reservation, and a replacement was admitted on the same CPU. (Codex 6; Grok noted the same comparison; present before) | **ACCEPT.** Cheap to close. Required: only the process recorded as owner can release. |
| X7 | The docstring says "Legacy exclusive launches remain unchanged", but the legacy reserve is now 0. The docstring and `AGENTS.md` also say admission needs only that available memory covers the caps. It also needs a free CPU and slot, a clean registry, the budget, and the unused caps of running jobs. The budget is `MemTotal` rounded down to whole GiB. (Codex 7, Grok 3) | **ACCEPT.** Required: both texts state what the code does. |

**Route.** The orchestrator, as the author, repairs X1, X3, X4, X6 and X7. Round 1 follows, with fresh Codex and
Grok legs. X2 is carried out at install time and recorded here.

## Round 1

**The change, round 1.** The round-0 change was repaired by its author, for X1, X3, X4, X6 and X7. Packet: the same
directory.
- after: `r1_s11c_guarded_run_after.py`, sha256 `c2c9ef01…`;
- diffs: `r1_guard.diff` (`8179eea4…`), `r1_AGENTS.diff` (`541fd8d3…`), and `r1_repair_only.diff` (`66b8607a…`, the
  diff against round 0, committed at `b9876bc3`);
- baseline `r1_baseline.sha256`;
- prompt `r1_review_prompt.md`, sha256 `b106429d…`: round 0's prompt plus these dispositions and checks on X1–X7.

**Legs (fresh).**
- Codex (gpt-6.1-sol, xhigh), NEEDS REVISION, four findings. Report: `r1_review_codex_final.txt`; evidence:
  `r1_codex_evidence/`.
- Grok, NEEDS REVISION, one finding. Report: `r1_review_grok.txt`; evidence: `r1_grok_evidence/` (scripts 01–04 and
  their stdout).

Both found X3, X4 and X6 now hold, the file-cache part of X1 holds, and the X5 reasoning is right. On X2, both found
the switch-over procedure enough, provided launches stay stopped through the whole install and the empty-registry
check. Both found these refusals intact: malformed, over-budget, duplicate-CPU, mismatched-limit and out-of-range
file-cache records.

| # | Finding (leg) | Disposition |
|---|---|---|
| Y1 | X1 is incomplete. Reclaimable slab is inside a unit's `memory.current` and also inside `MemAvailable`, but only the file LRU is subtracted. Grok's injected case admits caps 1 GiB past RAM; Codex's admits 3 GiB past. The committed guard refused both, because of its reserve. (Codex 1, Grok 1) | **ACCEPT.** This is the X1 property, not yet met. Required: the credited use excludes every part of a unit's memory that `MemAvailable` already counts. |
| Y2 | Admission reads every live unit's memory, then `MemAvailable`. A worker that frees memory between the two reads has its earlier use credited while the freed memory is counted as available. Codex's trace admits 7 GiB past RAM. Reordering the reads alone does not close it. (Codex 2; present before, hidden by the reserve) | **ACCEPT.** Required: the credited use is a bound that cannot credit memory freed, or count memory allocated, between the samples. |
| Y3 | "Distinct physical cores first" fails after releases. The CPU order is static and admission takes the first free entry, so once CPUs 7 and 6 are released CPU 7 is chosen beside its running sibling 3 while another core is idle. (Codex 3; introduced in round 0, when sibling threads were added) | **ACCEPT.** Required: a job is placed on an idle physical core whenever one exists. |
| Y4 | The docstring and `AGENTS.md` say pooled launches wait while a legacy run holds the admission lock. A pooled launch from the legacy run's own checkout refuses at once, on that checkout's execution lock. (Codex 4) | **ACCEPT.** Required: the texts distinguish refusal in the same checkout from waiting in another. |

**Route.** The author repairs Y1–Y4, and round 2 follows with fresh Codex and Grok legs. Y3 and Y4 are defects in this
change's own repairs. If round 2 finds a defect in what round 2 changes, the author changes (G4).

## Round 2

**The change, round 2.** The author repaired Y1–Y4. Packet: the same directory.
- after: `r2_s11c_guarded_run_after.py`, sha256 `93e61d5f…`;
- diffs against `40dfecc9`: `r2_guard.diff` (`84a49377…`) and `r2_AGENTS.diff` (`7010516f…`);
- `r2_repair_only.diff` (`cc35e036…`), the diff against round 1, which is committed at `0dfecf6e`;
- baseline `r2_baseline.sha256`;
- prompt `r2_review_prompt.md`, sha256 `18cdfdd6…`.

**Legs (fresh).**
- Codex (gpt-6.1-sol, xhigh), NEEDS REVISION, two findings. Report: `r2_review_codex_final.txt`; evidence:
  `r2_codex_evidence/`.
- Grok, NEEDS REVISION, one finding. Report: `r2_review_grok.txt`; evidence: `r2_grok_evidence/`.

Both found Y1, Y3 and Y4 now hold, and X3, X4 and X6 still hold. Both found X2's procedure enough, provided no launch
runs from the idle check through the empty-registry confirmation. Both found X5's reasoning right.

| # | Finding (leg) | Disposition |
|---|---|---|
| Z1 | Y2 is not met. Round 2 reads every unit, then `MemAvailable`, then every unit again, and keeps the smaller of the two samples. A unit that frees memory and reallocates it between the two walks shows high use in both samples, while `MemAvailable` was read during the free. Codex's injected case admits caps 6 GiB past RAM; Grok's, 4 GiB past. A docstring caveat does not meet Y2, and `AGENTS.md` states no caveat. (Codex 1, Grok 1) | **ACCEPT.** No atomic snapshot of every unit and `MemAvailable` exists. Required: each `MemAvailable` reading is compared only with live usage read next to it, and admission refuses unless every such comparison passes. A free and reallocation between non-adjacent reads then cannot raise admission. Whatever window remains is stated exactly in the docstring and in `AGENTS.md`. |
| Z2 | Round 2 keeps one sample per unit and drops the other before the limit and accounting checks. A second sample that shows an actual cap other than the reserved one, or inconsistent reclaimable accounting, is ignored, and an extra job is admitted. The one-sided branch overwrites the reclaimable figure, which also hides a malformed value. (Codex 2) | **ACCEPT.** Required: every sample of every live unit is validated against its reservation before any usage figure from it is used. The credit is carried separately from the raw accounting, which is never overwritten. |

**Route.** Both findings are in material round 2 changed: the sampling repair for Y2. Under G4 the author changes:
Codex writes the round-3 repair, and fresh Codex and Grok legs review it. Until it is accepted and installed, every
checkout keeps the committed guard. S9b work continues under it, one guarded job at a time.

## Round 3

**The change, round 3.** Under G4 the author changed: Codex (gpt-6.1-sol, xhigh, session `01a14268…`) repaired Z1
and Z2. Its brief is `r3_author_brief.md`, and its report, verification scripts and stdout are in
`r3_author_workspace/`. Packet: the same directory.
- after: `r3_s11c_guarded_run_after.py`, sha256 `1db855a3…`;
- diffs against `40dfecc9`: `r3_guard.diff` (`9ab64725…`) and `r3_AGENTS.diff` (`ab42037e…`);
- `r3_repair_only.diff` (`705120b4…`), the diff against round 2, which is committed at `f487979c`;
- baseline `r3_baseline.sha256`;
- prompt `r3_review_prompt.md`, sha256 `7a512cf7…`.

**Legs.** The round-3 repair is Codex-written, so the legs are a fresh Claude agent and Grok (G1). The orchestrator's
parts were reviewed by Codex legs in rounds 0–2.
- Fresh Claude (opus), NEEDS REVISION, one finding. Report: `r3_review_claude.md`; evidence: `r3_claude_evidence/`.
- Grok, CLEAR. Report: `r3_review_grok.txt`; evidence: `r3_grok_evidence/`.

Both found that the four intended changes, and X1, X3, X4, X6, X7, Y1–Y4, Z1 and Z2, now hold. Both found the X2
procedure enough and X5's reasoning right. Both found the remaining window stated in the docstring and in `AGENTS.md`:
memory that changes inside both read pairs. Claude also checked the author's report against the code.

| # | Finding (leg) | Disposition |
|---|---|---|
| W1 | The registry path is `Path.home()`, which follows `$HOME`, while `systemctl --user` reaches the same per-user manager whatever `HOME` is. Two launchers with different `HOME` therefore keep separate registries: in a real-process test, two pending 20 GiB jobs both got CPU 7 against a 31 GiB budget. The docstring and `AGENTS.md` say "one machine-wide registry". (Claude 1; the line dates from round 0) | **ACCEPT.** Lookup: every agent process running now has `HOME=/home/trevnorris` (`/proc/<pid>/environ`, below), so nothing has hit this. But the stated property is one registry per guarded-unit namespace, and `HOME` does not identify that namespace. Required: the registry path depends only on the user, matching `systemctl --user`. |
| — | Claude, outside the filter: the `AGENTS.md` diff deletes "Memory budgets and parallel pooling are the user's call … Use it only when the user authorizes parallel jobs". | **NO CHANGE.** The user (2026-10-15) said that line was set by Codex, not by them, and that the machine is dedicated: "You can use every last resource available". |
| — | Claude: a later checkout of a commit from before this change brings back the old guard with its own registry. | **PROCEDURE.** X2 covers it. Installation updates every checkout on this machine (11 copies of the launcher, below), and new builder exports come from commits after this change. |

Lookups (mechanical):
```
$ for p in <codex wrapper/session pids> <this shell>; do grep HOME /proc/$p/environ; done
433891 bash: HOME=/home/trevnorris
433893 codex: HOME=/home/trevnorris
467282 bash: HOME=/home/trevnorris
$ ls -d /var/projects/toy_physics*/scripts/s11c_guarded_run.py | wc -l
11
```

**Route.** W1 is in round-0 material, not in what round 3 changed, so the round-3 author repairs it. Round 4 is
scoped to that repair, with a fresh Claude agent and Grok.
