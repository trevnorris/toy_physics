# S15 separation departure — run record

Deliverable: two independent reconstruction scripts and final guarded transcripts, with the exact-object gaps below. No cross-engine object comparison or review was performed. Contract inputs are the task’s Part A (`91248b67`) and Part B (`7f52b4d2`).

## Commands and guard outcomes

**SymPy**

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_sympy_02 -- python3 /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
```

Guard exit `0`; child exit `0`; `limitsVerified: true`; `guardReason: null`; stderr `0` bytes. Guard runtime **1986.4803606170462 s**; cgroup peak **326676480 bytes**.

- [Script](../../scripts/S15_separation_departure.py); SHA-256 `72d694f5cf9de6480127e6bd50a87c95d086052df7c8e3cdfb919b563caa23b7`.
- [Transcript](../../scripts/out/S15_separation_departure.out); SHA-256 `3a0375e8ec2c8bbd4811cf32a26ecae5a090acb9796384a1edae36ba8e5d836b`.

The SymPy artifact inventory contains 336 Q sections, with no missing contract names or required section slots and with its end marker present.

**Wolfram**

```sh
python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/reading_wolfram_02 -- math -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
```

Guard exit `0`; child exit `0`; `limitsVerified: true`; `guardReason: null`; stderr `0` bytes. Guard runtime **133.07696074899286 s**; cgroup peak **3511361536 bytes**.

- [Script](../../mathematica/S15_separation_departure.wl); SHA-256 `5c094e8e8d2e459c3cf5c1de996d46a75ce5c6dcce3d590d3fde720205e9a980`.
- [Transcript](../../mathematica/out/S15_separation_departure.out); SHA-256 `651c200ffb36196fb3241e56dc9e9f85432d6f224e06c7664699577fbfcc56b7`.

The Wolfram artifact inventory contains 336 Q sections, with no missing contract names or required section slots and with its end marker present.

Runs were sequential, each with a 50 GiB cap, zero swap and no time limit. Wolfram used the two requested thread settings. Runtime is the guard’s `wallSeconds`; peak memory is the maximum recorded cgroup `memory.peak`, including the guard and descendants.

## Exact objects not obtained

In the following complete index, `s` ranges over `inner`, `outer`; `i` over `1,2,3`. Each pattern denotes every combination. The missing objects are item-3 closed forms of `E_geo`. Exact family integral expressions remain printed; the downstream item-4 conditions and Q evaluations under those conditions were not obtained for these entries.

**SymPy — 51 entries**, with `g` in `drain1`, `drain1+drain2`: `s/carried/component_i/labels_g`; `s/transport_centre/component_i/labels_g`; `s/transport_radius/component_i/labels_g`; `region/component_i/momentum_labels_g`; `region/component_i/rate_labels_g`; `region/component_i/rate_labels_drain2`.

Reason printed by this engine: NOT OBTAINED IN CLOSED FORM: angular integrand is outside this engine's exact polynomial integration route. The exact per-object names also appear under `EXACT_OBJECT_GAP` at the end of its transcript.

**Wolfram — 48 entries**, with `g` in `{drain1}`, `{drain1, drain2}`: `s/carried/component_i/labels_g`; `s/transport_centre/component_i/labels_g`; `s/transport_radius/component_i/labels_g`; `region/component_i/momentum_labels_g`; `region/component_i/rate_labels_g`.

Reason printed by this engine: NOT OBTAINED IN CLOSED FORM: non-polynomial angular factors remain outside this engine's exact integration route. The exact per-object names also appear under `EXACT_OBJECT_GAP` at the end of its transcript.

General-field integrals retained under item 3’s explicit exception are not counted as closed-form gaps. For the closed forms reached, item-4 conditions are printed as implicit equations in the named quantities with the engine’s exponent branches, alongside Q restricted to the condition.

Open work: the closed forms and dependent conditions indexed above. Commit IDs: **none created**; build review and the comparator are outside this task.
