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
