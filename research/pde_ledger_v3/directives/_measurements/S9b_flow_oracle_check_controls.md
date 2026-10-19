# S9b flow-only oracle check: final FORM controls record

**Visser FORM sensitivity is reported as a change in computed values.** The value checker reads the baseline and control residuals, prints both operands and their domains, computes their difference, then prints its exact zero test under that same printed domain. It does not count textual inequality.

Of **18** balance comparisons, **10** computed deltas print exactly `0` and their zero tests print `True`. In the other **8** (cases 1 and 2 of each engine's deflection and radar families), the computed delta vanishes exactly when the live density's radial derivative at the comparison radius is zero. Thus the FORM change moves those eight values away from that derivative-zero locus. This is a conditional value-sensitivity result; no constant-density premise or oracle density law was introduced. All **18** printed-domain identity guards pass. Literal output examples from the value-checker transcript:

```text
FORM_VALUE_DELTA: Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FORM_VALUE_DELTA", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> Derivative[1][py`rho$br$profile][py`s9b$r] == 0|>
FORM_VALUE_DELTA: Sqrt[2]*Sqrt[GM/r]*Derivative[1][rhoBr][r]
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FORM_VALUE_DELTA", "domain" -> "WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/1", "value" -> Derivative[1][rhoBr][r] == 0|>
FORM_VALUE_DELTA: 0
RESIDUAL_ZERO_TEST: <|"residual_tag" -> "FORM_VALUE_DELTA", "domain" -> "PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_DEFLECTION/case/3", "value" -> True|>
```

**The Hamilton–Lisle FORM control is dropped.** With exact zero tests it did not complete at 50 GiB in `run_023` and `run_025`. These are historical commands for the pre-removal code, not commands supported by the final driver:

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_023 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py --control hl_spatial_form
```

`run_023`: admitted; exit **1**, `limitsVerified=true`, `guardReason=null`. Workload **1011.3028902869992 seconds**; launcher **1011.3708832169941 seconds**; recorded peak **48076103680 bytes**. Partial stdout `_scratch/s11c/s9b-oracle-check/run_023/stdout`; SHA256 `b7ffb349d587891780970f605c08248803cb182b231aede507b89d4aa49b3b80`. The kernel prints `No more memory available.` and shuts down without instrument completion. The driver reports kernel exit 1 and `completion_record=false`.

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_025 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py --control hl_spatial_form
```

`run_025`: admitted; exit **1**, `limitsVerified=true`, `guardReason=null`. Workload **1015.0558314240043 seconds**; launcher **1015.1223335620016 seconds**; recorded peak **48079503360 bytes**. Partial stdout `_scratch/s11c/s9b-oracle-check/run_025/stdout`; SHA256 `b584901e8485f8383d732a13dfd335e539d8d66a86c88077daaf044357f1c9a7`. The kernel prints `No more memory available.` and shuts down without instrument completion. The driver reports kernel exit 1 and `completion_record=false`.

`run_023` used exact tests without a per-test bound. `run_025` used the attempted 1 GiB `MemoryConstrained` wrapper around exact tests; it printed no memory-bound fallback before kernel exhaustion. That wrapper did not contain the tests, and it has been removed. Both runs used a 50 GiB whole-run cap, zero swap and no time limit; sampled cgroup max/oom/oom_kill counters remained zero and the guard reported no kill. These partial logs document the resource stops, not a completed HL control. The sampled HL output artifact has been removed, and no transcript produced by another method substitutes for this dropped control.

The final baseline, Visser control and Visser-only value checker were then run one at a time from unchanged final sources. Their commands, transcripts and hashes follow. All were admitted with verified 50 GiB limits, zero swap, one kernel and no time limit. Nothing was committed.

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
