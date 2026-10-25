# S15 build run record — H1 repair

Date: 2026-10-24. Decision-list contract: `91248b67`. Reviewed baseline: `efe5ce73`, as supplied by the orchestrator. Repair scope: H1 in `S15_build_review_r0_disposition.md`. No git, git-annex or datalad command was run; no new commit was made.

Each engine replaces the unsupported `balance_construction` claim with a list of printed objects and adds a `region/component_i/balance/reduction` label for each component. Those labels explicitly state that no route in that engine reduces the residual with the supplied local balance. This uses H1's permitted disclosure alternative; no computation of the asserted implication or closure is claimed. The independent constructions, expressions, live fields, domains and groups are unchanged. Wolfram imports nothing from SymPy. No cross-engine output comparison was performed.

## SymPy

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/sympy_h1_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_two_drain_region.py
```

- Source: `research/pde_ledger_v3/scripts/S15_two_drain_region.py`; SHA256 `55f1716333c04b64ce0ea2745463cf180ceef3d692f2b01ddf769570611986dc`.
- Transcript: `research/pde_ledger_v3/scripts/out/S15_two_drain_region.out` (24641539 bytes).
- Transcript SHA256: `4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c`.
- Guard and child exit codes: `0`; `limitsVerified: true`; `guardReason: null`; stderr bytes: `0`.
- Child runtime: `387.7250224899908` seconds; guard runtime: `387.798950522003` seconds; peak cgroup memory: `145035264` bytes.
- Requested printed objects not produced: none. The balance residual is printed without a reduction route under the supplied local balance, as explicitly labelled.

## Wolfram

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/wolfram_h1_01 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_two_drain_region.wl
```

- Source: `research/pde_ledger_v3/mathematica/S15_two_drain_region.wl`; SHA256 `37fa7ec54b9c91c4a59b364f60f629c36180ac34f5836e8a3b7e0b9996d030ec`.
- Transcript: `research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out` (9797390 bytes).
- Transcript SHA256: `811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6`.
- Guard and child exit codes: `0`; `limitsVerified: true`; `guardReason: null`; stderr bytes: `0`.
- Child runtime: `10.22031184399384` seconds; guard runtime: `10.292321781016653` seconds; peak cgroup memory: `224870400` bytes.
- Requested printed objects not produced: none. The balance residual is printed without a reduction route under the supplied local balance, as explicitly labelled.

## Verification and measurement scope

The following non-CAS checks compare each transcript only to its own preserved baseline, excluding exactly `balance_construction` and the three new `balance/reduction` lines:

```sh
python3 _scratch/s15/build_work/h1/verify_transcript.py sympy
python3 _scratch/s15/build_work/h1/verify_transcript.py wolfram
```

For each command, the literal output includes `"other_lines_byte_identical": true`, `"completion_marker_present": true`, and `"transcript_is_symlink": false`. Thus the mathematical objects and their printed generality are unchanged. Source diffs contain only the two reporting statements per engine. All required object labels are present. Wolfram guard stdout contains progress lines and no diagnostic messages.

Each annex link was removed before its run; the new transcripts are plain files. Nothing under `.git/` was written, removed or chmodded. Runs were sequential, with 50 GiB, zero swap, no wall/CPU/inactivity deadline and one Wolfram kernel. Runtime is the guard's measured child wall time; peak memory is the maximum recorded kernel `memory.peak`, including the guard supervisor, not process RSS. Each command's guard directory contains `outcome.json`, limit/duration validation, resource samples, stdout and stderr. Working copies and mechanical check output remain under `_scratch/s15/build_work/h1/`.

Stopped after this repair and record. No computational closure claim is made; repair review remains with the orchestrator.
