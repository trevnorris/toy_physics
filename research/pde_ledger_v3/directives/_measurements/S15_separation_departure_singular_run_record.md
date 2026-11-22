# S15 Singular build — unfinished, November 12, 2026

The user stopped both object-125 workers to discuss a simulation in which the
drains determine the orbits. At 23:14:55 UTC both guards had recorded the stop
and released their 24 GiB reservations. The four continuation hooks exited.
The 144-object survey never started. Final object dispositions remain open.
Governing commits are `d4cd6d09` and `d544ef73`. No new commit was made.

**Trial results.** The compact construction gave Singular's derivative rows
as equations on a finite grid of reference separations. Thirteen of the
fourteen representative objects returned trial conditions. Object 125, the
both-drain contribution to the regional x-momentum rate, produced 41 of
115 rows. Direction 1 contains orders 0 through 20. Direction 5 contains
orders 0 through 19. Both workers stopped during the next order. Neither
worker wrote an engine completion receipt. The two guard exit codes were 1.

Both stopped transcripts were archived with their physical definitions.
Restoration reproduced the full captured transcripts and matched the raw
Singular bytes. Each completed equation has an individual polynomial hash.
All 41 rows matched an independently differentiable shell profile for their
direction at one reference separation. The first four rows in each direction
also matched direct scalar differentiation as polynomial identities. On the
separate region-momentum object 124, all 16 equations saved from the expanded
baseline had matched the compact representation.

**Open physics.** The trial measures the fluid momentum in a shell. Determining
the acceleration of a throat still requires its inertia and the net momentum
delivered to it, including exchange through its boundaries. See the
[S15 open questions](../S15_decision_list.md#retained-open).
The [unfinished run manifest](S15_separation_departure_singular_run_manifest.json)
contains the attempt receipts, benchmark evidence and archive hashes.

<details>
<summary>Finite-grid trial attempts, recorded UTC times, and resource use</summary>

| Trial attempt | Start, UTC | Supervisor end, UTC | Cap (GiB) | Outcome | Guard wall seconds | Recorded peak bytes |
|---|---|---|---:|---|---:|---:|
| finite_grid01_object002 | 2026-11-10T17:02:09.025233+00:00 | 2026-11-10T17:02:21.749132+00:00 | 54 | trial conditions returned | 12.700002054916695 | 176066560 |
| finite_grid01_object124 | 2026-11-10T17:04:29.648678+00:00 | 2026-11-11T17:25:07.533912+00:00 | 54 | explicit user stop | 87637.85839541396 | 3941429248 |
| compact_grid01_first_001 | 2026-11-11T17:26:37.897959+00:00 | 2026-11-11T17:26:39.304761+00:00 | 53 | trial conditions returned | 1.3854872800875455 | 59658240 |
| compact_grid01_object002 | 2026-11-11T17:27:25.572888+00:00 | 2026-11-11T17:27:31.472559+00:00 | 1 | trial conditions returned | 5.8656953000463545 | 113577984 |
| compact_grid01_object124 | 2026-11-11T17:27:25.574249+00:00 | 2026-11-11T23:08:49.762096+00:00 | 16 | trial conditions returned | 20484.105882315896 | 9598545920 |
| compact_grid01_object120 | 2026-11-11T23:44:03.804639+00:00 | 2026-11-11T23:44:15.444193+00:00 | 8 | trial conditions returned | 11.61175618506968 | 207888384 |
| compact_grid01_object121 | 2026-11-11T23:44:03.805917+00:00 | 2026-11-11T23:48:29.242654+00:00 | 8 | trial conditions returned | 265.380098619964 | 848453632 |
| compact_grid01_object012 | 2026-11-11T23:44:03.808088+00:00 | 2026-11-11T23:44:04.851656+00:00 | 2 | trial conditions returned | 0.8850700319744647 | 56836096 |
| compact_grid01_object013 | 2026-11-11T23:44:03.808093+00:00 | 2026-11-11T23:44:04.726837+00:00 | 2 | trial conditions returned | 0.8147568251006305 | 56844288 |
| compact_grid01_object125 | 2026-11-11T23:44:03.809469+00:00 | 2026-11-12T23:14:55.596869+00:00 | 24 | explicit user stop | 84651.51047008205 | 13794443264 |
| compact_grid01_object003 | 2026-11-11T23:44:03.817322+00:00 | 2026-11-11T23:44:05.991987+00:00 | 2 | trial conditions returned | 1.8326413279864937 | 58748928 |
| compact_grid01_object014 | 2026-11-11T23:44:03.820088+00:00 | 2026-11-11T23:44:05.159066+00:00 | 2 | trial conditions returned | 0.8820763928815722 | 56840192 |
| compact_grid01_object007 | 2026-11-11T23:44:03.824035+00:00 | 2026-11-11T23:44:05.528030+00:00 | 2 | trial conditions returned | 1.0878107189200819 | 57303040 |
| compact_grid01_object001 | 2026-11-11T23:44:04.850510+00:00 | 2026-11-11T23:44:05.971219+00:00 | 2 | trial conditions returned | 0.9833936791401356 | 56840192 |
| compact_grid01_object018 | 2026-11-11T23:44:04.963741+00:00 | 2026-11-11T23:44:05.981155+00:00 | 2 | trial conditions returned | 0.8348734930623323 | 56881152 |
| compact_grid01_object060 | 2026-11-11T23:44:05.275382+00:00 | 2026-11-11T23:44:07.099407+00:00 | 2 | trial conditions returned | 1.668251679977402 | 59670528 |
| compact125_dirs54 | 2026-11-12T09:09:05.995161+00:00 | 2026-11-12T23:14:55.612824+00:00 | 24 | explicit user stop | 50749.558699744055 | 11292332032 |

</details>

<details>
<summary>Historical unfinished report of November 10</summary>

### S15 Singular build — UNFINISHED (2026-11-10)

**This is an unfinished-build report, not the directive's final report.** No full
run took place; no object is closed or assigned a final label, and no trial output
is presented as the final transcript. The two completed trial outputs are
unreviewed. Governing commits: directive `d4cd6d09`, resume note `d544ef73`.
No new commits, agents or review rounds were made.

**Construction tried.** Python formed the exact surface/region objects and their
linear analytic state systems. Singular extracted direction coefficients and
closed a polynomial row module under the state-system derivative, dividing only
whole rows by known nonzero geometric factors. This replaced attempt01's
polynomial-ideal closure. Exponents and other parameters remained polynomial
variables. The construction and its stated premises are in
[`CONSTRUCTION.md`](../../singular/CONSTRUCTION.md).

**Attempts.** All times below are UTC. Start times are the guard invocation's;
new end times are the supervisor's observation of guard return. Wall seconds are
the guard's recorded value, not subtraction of these two timestamps. Peaks are
the maxima recorded by the guard's resource sampler, in bytes.

| Attempt | Start → end | Cap (GiB) | Outcome | Wall seconds | Peak bytes |
|---|---|---:|---|---:|---:|
| `attempt01` | Nov 9 03:42:21.409292 → stop reported 18:32:56 | 56 | Historical user stop | Not recorded | 12057948160 |
| `module_trial01_object000` | Nov 9 20:18:28.590900 → 20:18:30.690566 | 54 | Code error | 2.038331691874191 | 105451520 |
| `module_trial02_object000` | Nov 9 20:20:20.698930 → 20:20:22.512545 | 54 | Trial conditions returned | 1.7882936710957438 | 77996032 |
| `module_trial02_object001` | Nov 9 20:21:33.865702 → 20:21:34.483440 | 54 | Trial conditions returned | 0.5955363628454506 | 63672320 |
| `module_trial02_object002` | Nov 9 20:21:34.785268 → Nov 10 16:44:18.291612 | 54 | Explicit user stop | 73363.48406163091 | 10287329280 |

**What returned and what did not.** Indices 0, 1 and 2 are the inner-surface,
component-1 carried-current groups labelled `{drain1}`, `{drain2}` and
`{drain1,drain2}` respectively. The first two returned trial conditions (8 and
0 printed equations). The third printed its item-3 object and reached
`WORKING_ON: differential module round 11`, `BASIS_SIZE: 149`; its stdout is
1,705,568,406 bytes. It returned no completed item-4 condition set. Of 14 selected
trial objects, two completed, one was stopped, and eleven were not attempted.
No region trial completed. The full 144-object run did not start.

**Stop and missing evidence.** The user said, “Stop and write unfinished-build
report.” The unit was stopped with its guard launcher alive. The launcher
recorded exit 1 and released its reservation; the batch ended without starting
another object. No guarded unit or owned computation remained active at the
post-stop check. No out-of-memory event was recorded. The stop is a user decision,
not a measured proof that the method cannot finish. The stopped child's
`child-outcome.json` and engine `run_data.json` are absent because unit termination
ended those writers. Attempt01 also lacks its launcher outcome because its
launcher had already ended; its accepted note reports about 14 h 50 min, but no
completed guard wall-time measurement exists. Missing values are not inferred.

**Evidence and open question.** The adjacent
[`run manifest`](S15_separation_departure_singular_run_manifest.json) embeds the
quoted output, guard facts, attempt metadata and source hashes. Its reproduction
commands replace temporary path arguments with archived source/work roles;
original command hashes are retained. It is an unfinished-trial evidence
manifest, not a complete final-run manifest. Sources for the historical and
syntax-error attempts are preserved under `singular/preserved/`. The current
engine and evidence collector remain uncommitted. The open question is whether
a smaller exact state representation and explicit parameter-case treatment can
complete the both-drain and region objects. No further method or run was started.

</details>
