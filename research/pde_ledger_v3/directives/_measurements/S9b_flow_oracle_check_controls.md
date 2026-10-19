# S9b flow-only oracle check: FORM controls

Execution evidence only; no physical conclusion or review clearance. Controls alter computed copies after printing the unchanged quoted oracle equations. Both retain source fields and domain handling. All commands launch a single `math -script` child inside the pooled guard, with a 12 GiB cap, zero swap and no time limit.

- `hl_spatial_form`: admitted, exit 0, limits verified, no guard reason; workload 126.56015280199063 seconds, launcher 126.62687389799976 seconds, cgroup peak 918228992 bytes. Transcript `scripts/out/S9b_flow_oracle_check_hl_spatial_form.out`; SHA256 `ea08717e0e3087f7038273119e69fc0c7b1bc7f1cf151fe98d9ea7d0a55fab46`.

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 12 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_011 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py --control hl_spatial_form
```

- `visser_form`: admitted, exit 0, limits verified, no guard reason; workload 113.61128813499818 seconds, launcher 113.67690052300168 seconds, cgroup peak 374734848 bytes. Transcript `scripts/out/S9b_flow_oracle_check_visser_form.out`; SHA256 `3c8a1dffcb95314b05cd571cae0dbebc411e0512c525fd0faebc1ccadc611217`.

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s9b-oracle --memory-gib 12 --log-directory /var/projects/toy_physics/_scratch/s11c/s9b-oracle-check/run_012 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check.py --control visser_form
```

Visser operator FORM control moves the density factor outside the divergence. The full operator is evaluated with live S9b density in each engine. The first differing native balance-residual pair is printed below (baseline, then control); these are literal computed outputs, not hand-derived expressions.

```text
CONTINUITY_ORACLE_EVALUATED_EQUATION: (2*py`rho$br$profile[py`s9b$r]*py`s9b$V[py`s9b$r])/py`s9b$r + py`s9b$V[py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r] + py`rho$br$profile[py`s9b$r]*Derivative[1][py`s9b$V][py`s9b$r] == 0
CONTINUITY_ORACLE_EVALUATED_EQUATION: py`rho$br$profile[py`s9b$r]*((2*py`s9b$V[py`s9b$r])/py`s9b$r + Derivative[1][py`s9b$V][py`s9b$r]) == 0
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: (Sqrt[2]*Sqrt[py`GM]*py`s9b$r*Derivative[1][py`rho$br$profile][py`s9b$r] - Sqrt[2]*Sqrt[py`GM/py`s9b$r]*py`s9b$r^(3/2)*Derivative[1][py`rho$br$profile][py`s9b$r])/py`s9b$r^(3/2)
EXCHANGE_BALANCE_OPERATOR_RESIDUAL: Sqrt[2]*Sqrt[py`GM/py`s9b$r]*Derivative[1][py`rho$br$profile][py`s9b$r]
```

Textually differing balance residuals: 8 of 18. Text inequality is the FORM sensitivity observation, not a claim that all unevaluated radical forms are canonical. Original source-exchange residuals and density zero sets remain printed.

Hamilton–Lisle spatial FORM control adds the squared quoted (2) profile times the radial differential squared to (1), by computation from the transcription. It changes the extracted rr coefficient and optical residuals. The coordinate coefficient and the construction-normalized unit coefficient are separately labelled; the latter remains fixed by construction. The final control reaches all 242 structural checks and completion without Wolfram messages.

The earlier local control's cgroup recorded `pids.events: max 3`, `pids.current: 31`, `pids.max: 32`; its kernel threads were all waiting. That blocked run was manually stopped, not killed by the guard. The check now sets `ParallelThreadNumber` and `MKLThreadNumber` to 1 and disables threaded inequality disjunctions, preserving the cap. Wolfram documents the internal thread control in [Parallel Computation](https://reference.wolfram.com/language/Compile/tutorial/Parallel.html). The reproduced wait and successful repaired control are evidence for this instrument; the historical reviewer run was not diagnosed directly.

Debug-run guard records remain in scratch. No guard refusal or guard-triggered kill occurred; incomplete runs 006 (boundary-evaluation bug) and 009 (blocked internal threading) were manually stopped. Their transcripts are not result evidence. No independent review was launched.
