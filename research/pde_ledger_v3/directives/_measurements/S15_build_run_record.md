# S15 build run record

Contract: user-supplied sections of the S15 decision list at `91248b67`. Date: 2026-10-24. Checkout HEAD: `a9acfb1d8ff85a8b33893be0c2af4624efb9abe8`. No new commits; review remains with the orchestrator.

Both independent final sources produced their transcripts in sequential guarded runs. Wolfram imports no SymPy source or output. No cross-engine comparison was performed. The scripts retain live fields and exact symbolic surface/volume integrals; pressure expansion conditions and order are printed with the objects. Transcript label inventories include the required surface pieces, all seven label groups per grouped piece, pressure objects, region momentum and balance, and dimensional objects. Neither engine omitted a requested object. This records execution and object coverage, not mathematical review clearance.

## SymPy

Command:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/sympy_20261024_03 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_two_drain_region.py
```

- Source: `research/pde_ledger_v3/scripts/S15_two_drain_region.py`.
- Source SHA256: `faf6ad225bc9ac5afd03953ecd0388de4dc6e88333888faecd19f6e7bfbff152`.
- Transcript: `research/pde_ledger_v3/scripts/out/S15_two_drain_region.out` (24641288 bytes).
- Transcript SHA256: `85d99aed185464f9a9e68b2d0cc7b64193c073fdc5c321c0c60825069b2e8259`.
- Guard: exit code `0`; child exit code `0`; `limitsVerified: true`; `guardReason: null`; stderr bytes `0`.
- Child wall runtime: `387.2189354400034` seconds. Guard wall runtime: `387.28407184098614` seconds.
- Peak cgroup memory: `144961536` bytes.

## Wolfram

Command:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/wolfram_20261024_02 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_two_drain_region.wl
```

- Source: `research/pde_ledger_v3/mathematica/S15_two_drain_region.wl`.
- Source SHA256: `6ea5638ab1184e3f5bd77fb86ea315052e7e1b6fe78c738135f3b9099aa2d0b7`.
- Transcript: `research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out` (9797152 bytes).
- Transcript SHA256: `8774462d38a926329ce15f524001972b9af2a5b0225ec19a2f48b7c3b2982fd1`.
- Guard: exit code `0`; child exit code `0`; `limitsVerified: true`; `guardReason: null`; stderr bytes `0`; no Wolfram diagnostic messages in guard stdout.
- Child wall runtime: `10.320520654000575` seconds. Guard wall runtime: `10.391379517997848` seconds.
- Peak cgroup memory: `225038336` bytes.

## Measurement scope and stop

Runtime is the guard's `childOutcome.wallSeconds`; peak memory is the maximum recorded kernel `memory.peak`, including the guard supervisor and CAS in their cgroup, not process RSS. Both runs used 50 GiB, zero swap, one job at a time, and no wall/CPU/inactivity deadline. Guard evidence is in each command's log directory (`outcome.json`, `child-outcome.json`, `limit-validation.json`, `duration-validation.json`, `resource-samples.jsonl`, `stdout`, `stderr`); the `_scratch/s15/guard` link resolves under `_scratch/s11c/s15`.

Earlier launch-path rejection and engine code-error attempts are not the reported final runs; their working copies and evidence remain under `_scratch/s15/`. Completed final jobs were not replayed. Missing requested objects: none. Open work: the later, separate comparison and orchestrator review. No force law, coefficient, normalization, mass, throat-force identification, or connection to J, Q, Theta_Q or N_infinity_3 is claimed. Stopped after producing this run record.
