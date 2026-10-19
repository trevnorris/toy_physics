# S9b flow-only oracle check: final exact-zero-test run record

Deliverable: restore exact zero tests, remove sampling and the ineffective per-test bound, drop the Hamilton–Lisle FORM control, and regenerate the baseline, Visser FORM control and Visser-only value comparison.

**All three final runs completed.** Exact tests use the same `FullSimplify[residual == 0, Assumptions -> printedDomain, TimeConstraint -> Infinity]` method as `run_019`, with internal cache clearing and memoization of identical residual/domain pairs in the main instrument. There is no per-test memory bound or containment claim. The sampler and sampling metadata were removed. Source hashes remained fixed across all three runs. Nothing was committed; comparison base `2caec139`.

The **2806** original operand/domain/residual triples in each main transcript match its corresponding `2caec139` transcript exactly. All **2806** baseline zero-test records also match `run_019` exactly. Each main run completes all **242** structural checks, with no Wolfram messages. The value checker completes all **18** Visser balance comparisons, with equal printed domains and no messages. Exact zero tests may return a conditional zero locus; zero `NOT_ESTABLISHED` entries is not a claim that every residual vanishes identically.

Every final run used a **53687091200-byte** cgroup and native address-space cap (50 GiB), zero swap, one kernel, and no time limit. Admission succeeded and limits were verified throughout; no run was refused, killed or exhausted memory.

**baseline (run_032)**

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_032 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py
```

Admitted; exit **0**, `limitsVerified=true`, `guardReason=null`. Workload **291.42767466499936 seconds**; launcher **291.4934450340079 seconds**; recorded cgroup peak **3617349632 bytes**. Transcript `scripts/out/S9b_flow_oracle_check.out`; SHA256 `e66ff63512e4f044bab2a60eb327665db865179531fe6077622de8fab553a25e`. Exact zero tests: **2806**. Printed `NOT_ESTABLISHED`: **0**; reason counts: **none**.

**visser_form (run_033)**

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_033 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py --control visser_form
```

Admitted; exit **0**, `limitsVerified=true`, `guardReason=null`. Workload **291.220714487994 seconds**; launcher **291.2882461309928 seconds**; recorded cgroup peak **3616137216 bytes**. Transcript `scripts/out/S9b_flow_oracle_check_visser_form.out`; SHA256 `e194229c3ea66c4329f6f2cc02fd34307d9aaf29889c39bae30ae33f53352009`. Exact zero tests: **2806**. Printed `NOT_ESTABLISHED`: **0**; reason counts: **none**.

**values (run_034)**

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_034 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check_values.wl
```

Admitted; exit **0**, `limitsVerified=true`, `guardReason=null`. Workload **5.871138028000132 seconds**; launcher **5.937171835001209 seconds**; recorded cgroup peak **276946944 bytes**. Transcript `scripts/out/S9b_flow_oracle_check_values.out`; SHA256 `bedcbfe9e24e6cd0a090549a1e0d5ff41ccff86da938e7749561d74d797c8a4a`. Exact zero tests: **18**. Printed `NOT_ESTABLISHED`: **0**; reason counts: **none**.

The **Hamilton–Lisle FORM control is dropped by user instruction**. Its exact-test attempts did not complete at 50 GiB (`run_023` and `run_025`); the latter's per-test `MemoryConstrained` wrapper did not contain the failure. The commands, partial stdout hashes and measurements are cited in `S9b_flow_oracle_check_controls.md`. No alternate-method transcript is presented as that control.

Files changed: `scripts/S9b_flow_oracle_check.py` removes the HL control option and sampler provenance; `scripts/S9b_flow_oracle_check_data/S9b_flow_oracle_check.wl` restores exact tests and removes the HL FORM mutation; `S9b_flow_oracle_check_values.wl` restores exact tests and retains only Visser comparisons. Removed `S9b_flow_oracle_check_numeric.wl` and `scripts/out/S9b_flow_oracle_check_hl_spatial_form.out`. Regenerated the three transcripts named above, this record and the controls record; updated only the baseline hash in the identification listing. Identification statuses and the quoted oracle equations are unchanged. No new physical equivalence or review clearance is claimed.
