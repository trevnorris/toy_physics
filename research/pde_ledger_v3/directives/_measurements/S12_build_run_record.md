# S12 build run record — flux face repair

2026-10-20. Contract: `../S12_decision_list.md`; repair authorized against the reviewed baseline at `63d99406`. No new review, engine-output comparison, or commit was performed.

Both engines now represent J independently in all three open regions, with free one-sided traces, and retain arbitrary distributions supported at either face. The supplied balance is differentiated to obtain the jump contributions. The additional supported-flux derivatives remain formal because their forms are OPEN. No continuity, exterior vanishing, boundary condition, or delta profile is imposed. The previous n, u, faces, general interior order, conditional two-offset G specialization, and dimensional computations retain their generality.

Wolfram isolates supported-flux contributions with algebraic markers after substitution into the balance. SymPy extracts the flux differential operator from the supplied balance and applies it to the supported distributions. The implementations import neither each other nor each other's outputs.

Both final runs used the pooled guard serially, Wolfram first: 50 GiB memory cap, zero swap, 32 tasks, verified limits, no execution deadlines. Each final transcript equals its own guard stdout. Output-tag and field-identifier coverage checks found no omissions; they do not adjudicate the physics.

## Wolfram

Script: `research/pde_ledger_v3/mathematica/S12_order_source.wl` (repository relative). Source SHA256: `57aa7e0f35b62fbdb1cb5660393c8fb060e333b2accca0fa1408f0e071786033`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/wolfram_flux_faces_02 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S12_order_source.wl
```

Guard and child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Workload wall time: **4.317915863997769 s**; guard wall time: **4.38882244497654 s**. Peak guarded-cgroup memory: **222773248 bytes** (maximum recorded `memory.peak`).

Transcript: `research/pde_ledger_v3/mathematica/out/S12_order_source.out` (repository relative), **199694 bytes**. SHA256: `eab1149dd2b5cd0a05a982ab2c382710cdd23254032012b180627e9524935d84`.

## SymPy

Script: `research/pde_ledger_v3/scripts/S12_order_source.py` (repository relative). Source SHA256: `b052834f2d78909daa584de36761d3609adc73153c1ff017f7c178fbc308a0ee`.

```bash
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s12 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s12/guard/sympy_flux_faces_01 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S12_order_source.py
```

Guard and child exit **0**, `limitsVerified: true`, `guardReason: null`, stderr **0 bytes**. Workload wall time: **1.167489800020121 s**; guard wall time: **1.2381242029950954 s**. Peak guarded-cgroup memory: **55160832 bytes** (maximum recorded `memory.peak`).

Transcript: `research/pde_ledger_v3/scripts/out/S12_order_source.out` (repository relative), **204491 bytes**. SHA256: `440a207292b29c7ceb52180a01897bc42bf5a6258cf16148c55ce72261d7113d`.

## Coverage and run evidence

Both final transcripts contain the general source, both exterior regions, the interior, both face contributions, the labelled conditional G specialization, and dimensional operands/residuals. Face objects include unrestricted regional flux traces and formal derivatives of arbitrary face-supported flux distributions. Requested objects not produced: **none**. Contract OPEN items remain untouched; inventories were not undertaken.

An earlier admitted Wolfram attempt used the same Wolfram command above with log directory `/var/projects/toy_physics/_scratch/s12/guard/wolfram_flux_faces_01`: child/guard exit **1**, `limitsVerified: true`, `guardReason: null`, workload **4.36955198799842 s**, guard **4.440165937994607 s**, peak **223039488 bytes**. This was a pure-function binding bug in the repair, fixed before the final run; no guard refusal or kill occurred. Its partial transcript is `_scratch/s12/build_work/repair_flux_faces/wolfram_flux_faces_01.partial.out`.

Guard receipts, limit/duration validations, and resource samples are in the command log directories. `_scratch/s12/guard` remains the existing alias to `_scratch/s11c/s12_guard`; the guard was not changed. Working scripts and mechanical coverage/hash metadata are in `_scratch/s12/build_work/repair_flux_faces/`. The committed annex transcript objects were preserved; new transcripts replace their working-tree links. Commits: **none**.
