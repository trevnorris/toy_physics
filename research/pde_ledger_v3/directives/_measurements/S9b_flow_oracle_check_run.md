# S9b flow-only oracle check: run record

Produced the transcript-sourced comparison instrument, its final guarded-run output, and the identification listing. No physical conclusion or review clearance is claimed. No commits were made.

Files to read (relative to `research/pde_ledger_v3/`):
- `scripts/S9b_flow_oracle_check.py`: mechanical reader and single-kernel launcher.
- `scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl`: computation from both engines' operands.
- `scripts/S9b_flow_oracle_check_data/oracle_quotes.json` and `identification_map.json`: cited oracle transcriptions and map metadata.
- `scripts/out/S9b_flow_oracle_check.out`: complete final stdout.
- `directives/_measurements/S9b_flow_oracle_check_identifications.md`: every printed `UNMATCHED` map entry and its literal reason/source.

Command:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 12 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_005 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py
```

The driver launches `math -script /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl` inside the same guard, one kernel. Guard outcome: **admitted, exit 0, limitsVerified=true, guardReason=null**. Verified memory cap **12884901888 bytes**, swap cap **0**, with no time limit. Workload wall time **15.142501781003375 seconds**; launcher wall time **15.217535862000659 seconds**. Final cgroup peak memory **493305856 bytes**; swap use and OOM-kill counters remained zero. Evidence: `_scratch/s11c/s9b-oracle-check/run_005/` at repository root.

Transcript SHA256: `ba22d74705e4aa2704f750b8dc8bdd84b28ba009d5734117dac4ca60a3c3428b`.

Validation: final files and source hashes match the emitted provenance; annex targets match accepted revision `503736ff`. The transcript contains four general/retained profile packages, all 18 source exchange cases, both supplied-dispersion comparisons, and 14 line-element `UNMATCHED` records. It reaches its computed structural guards and completion count, with no Wolfram messages. This is execution validation, not an independent review.

Oracle copies, extracted text, page images and working inputs remain under repository-root `_scratch/s9b_oracle_check/`. Local PDF versions: gr-qc/0411060v2 (31 Aug 2006), gr-qc/9712010v2 (15 Dec 1997); PDF hashes and per-equation citations appear in the transcript.

Open items: the **18** literal `UNMATCHED` identification-map entries in the accompanying listing; the absence of an HL continuity operand and an S9b line element. Unruh full-text verification, Cahill/Sbitnev comparisons, and the register's historical second-route question remain outside this check. Build review remains with the orchestrator; none was launched here. The S9b record, register, unrelated working files and other jobs were untouched.
