# Runtime Falsification README

This directory now has a working chain for turning the current reduced
moving-throat / Branch-B work into concrete falsification checks.

The goal here is not to prove the whole program. The goal is to kill bad
physical identifications and bad simulation behavior as quickly as possible.

## Scope

The files added in the latest pass are centered on:

- the dimensional port map,
- the runtime monitor identities,
- snapshot adapters,
- a fail-fast classifier,
- and a partial front-end screen for the existing monopole solver output.

The relevant audited steps are:

- `step_32_parent_throat_action_dimensional_port_map_sympy.py`
- `step_33_parent_throat_action_runtime_monitor_falsifier_sympy.py`
- `step_34_parent_throat_action_cfd_runtime_postprocessor_sympy.py`
- `step_35_parent_throat_action_failfast_classifier_sympy.py`
- `step_36_parent_throat_action_snapshot_adapter_sympy.py`
- `step_37_parent_throat_action_monopole_jsonl_fastscreen_sympy.py`

The main runtime tools are:

- `cfd_runtime_monitor_postprocess.py`
- `cfd_runtime_failfast.py`
- `cfd_snapshot_adapters.py`
- `single_throat_monopole_jsonl_fastscreen.py`

## What Is Fixed

These constants are treated as frozen carry-forward values from the 4D summary
and 1PN bridge summary:

- `n = 5`
- `alpha_opt = (n-1)/2 = 2`
- `kappa_add = 1/2`
- `alpha^2 = 3/4`
- `K_vec = 2/pi^2`

Important interpretation points:

1. `c_s` and `c` are decoupled in this model.
2. `step_31` only falsified the naive direct-SI choice `S_port = 1`.
3. `step_32` makes the exact source-map bridge explicit:

   `mhat0^2 * S_port * P0 = 54 G c_s^5 / (5 a^5 c^5)`

   On the hardcoded Branch-B target sheet where `mhat0^2 P0 = 54/5`, this
   collapses to:

   `S_port = G c_s^5 / (a^5 c^5)`

So the reduced patch does not by itself determine a unique physical `c_s`.

## Files And Roles

### Symbolic / audited steps

- `step_32_*`: dimensional/source-map bridge
- `step_33_*`: exact runtime monitor identities and hard falsifiers
- `step_34_*`: CFD postprocessor self-test
- `step_35_*`: fail-fast classifier self-test
- `step_36_*`: snapshot adapter self-test
- `step_37_*`: existing monopole JSONL fast-screen self-test

Each has a matching `*_notes.md` file.

### Runtime tools

#### `cfd_runtime_monitor_postprocess.py`

Consumes a simulation snapshot and computes:

- `S_rho`
- `R_cont`
- `R_Pois_exact`
- `R_Pois_lin`
- shell-averaged `Q_r`
- shell-averaged / tail `mu_eff^2`
- shell-averaged / tail `alpha_fit` when optics fields are present

Supported snapshot schemas:

1. `full_4d`

   Required keys:

   - `rho`
   - `jx`
   - `jy`
   - `jz`
   - `jw`
   - `W`
   - `x`
   - `y`
   - `z`
   - `w`

2. `projected_3d`

   Required keys:

   - `rho_brane`
   - `Jx_brane`
   - `Jy_brane`
   - `Jz_brane`
   - `S_rho`
   - `x`
   - `y`
   - `z`

Optional keys for either path:

- `dWdw`
- `dt_rho`
- `phi3`
- `rho0`
- `N_probe`
- `Phi_eff`

#### `cfd_runtime_failfast.py`

Consumes the JSON summary produced by the postprocessor and emits:

- `PASS`
- `FAIL`
- `INCOMPLETE`

with explicit reasons.

Default thresholds:

- continuity relative RMS max: `0.05`
- exact Poisson relative RMS max: `0.05`
- `Q_r` tail coefficient of variation max: `0.05`
- `|mu_eff^2|` tail median max: `0.25`
- `|alpha_fit - 2|` tail mean max: `0.1`
- `std(alpha_fit)` tail max: `0.1`

Important verdict behavior:

- missing optics makes the verdict `INCOMPLETE` unless another channel already
  fails;
- near-zero `S_rho` makes projection residuals non-load-bearing and defaults to
  `INCOMPLETE`;
- source-free exterior synthetic tests may opt out with
  `require_projection_source = False`.

#### `cfd_snapshot_adapters.py`

Converts solver-side state dumps into the runtime-monitor snapshot schema.

Supported upstreams:

1. `wavefunction-4d`

   Input is a 4D complex wavefunction dump with:

   - `psi` or `psi_real`/`psi_imag`
   - `x,y,z,w`
   - `W` unless `--allow-uniform-W` is intentionally used

   Optional:

   - `dWdw`
   - `dt_rho`
   - `phi3`
   - `rho0`
   - `N_probe`
   - `Phi_eff`
   - `mass` or `m`
   - `hbar`

   It computes currents via:

   `j = (hbar/m) Im(conj(psi) grad psi)`

   `W` is required by default. Use `--allow-uniform-W` only when a flat
   projection weight is an intentional modeling choice.

2. `monopole-3d`

   Input is a projected 3D solver state with:

   - `rho`
   - `mx`
   - `my`
   - `mz`
   - `x`
   - `y`
   - `z`

   and either:

   - direct `S_rho`, or
   - `W`, `Mdot`, and `lambda` / `lambda_bulk`

   In the second case it reconstructs the monopole source law:

   `S_rho = -Mdot * W + lambda * Mdot / V_domain`

   If `V_domain` is present it is used directly; otherwise the adapter computes
   it from the uniform grid axes.

#### `single_throat_monopole_jsonl_fastscreen.py`

This is a partial early-kill screen for the current
`single_throat_monopole.py` solver logs.

It reads JSONL stdout and checks the last diagnostic event for:

- `dP_slope` near `-1`
- `geff_slope` near `-2`
- `mach_max <= 0.6`
- enough fit points

This is weaker than the full snapshot-based runtime monitor, but it works with
the solver as currently written. Malformed JSONL lines are reported as
`INCOMPLETE` rather than crashing the screen.

## Exact Runtime Monitors

The core projected source is:

`S_rho = -[W J_w] + ∫ W' J_w dw`

The primary exact residuals are:

- `R_cont = dt_rho + div(J_brane) - S_rho`
- `R_Pois_exact = rho_brane * lap(phi3) - S_rho + dt_rho + (grad rho_brane) · v_brane`
- `R_Pois_lin = rho0 * lap(phi3) - S_rho`

The main hard falsifiers are:

1. Exterior Yukawa-like screening:
   - `Q_r(r)` does not plateau
   - or `mu_eff^2(r)` tends to a nonzero exterior constant

2. Wrong weak-field optical coefficient:
   - `alpha_fit` does not sit near `2`

## How To Run

### 1. Re-run the audited symbolic steps

```bash
python verify_step_outputs.py \
  step_32_parent_throat_action_dimensional_port_map_sympy.py \
  step_33_parent_throat_action_runtime_monitor_falsifier_sympy.py \
  step_34_parent_throat_action_cfd_runtime_postprocessor_sympy.py \
  step_35_parent_throat_action_failfast_classifier_sympy.py \
  step_36_parent_throat_action_snapshot_adapter_sympy.py \
  step_37_parent_throat_action_monopole_jsonl_fastscreen_sympy.py
```

Logs go to:

```text
/tmp/em_projected_verify
```

### 2. Full snapshot-based falsification path

#### From a 4D wavefunction dump

```bash
python cfd_snapshot_adapters.py \
  wavefunction-4d raw_wave.npz runtime_snapshot.npz

python cfd_runtime_monitor_postprocess.py \
  runtime_snapshot.npz --output-json summary.json

python cfd_runtime_failfast.py \
  summary.json --output-json verdict.json
```

If the dump intentionally has no `W`, add `--allow-uniform-W`.

#### From a 3D monopole state dump

The current 4D research monopole solver can emit this dump directly:

```bash
python /var/projects/toy_physics/research/4d/scripts/single_throat_monopole.py \
  --N 256 \
  --steps 2000 \
  --diag_every 200 \
  --snapshot_out raw_state.npz \
  > monopole_N256_snapshot.log
```

Use the same physical/numerical flags as the run you want to reproduce; the
important addition is `--snapshot_out`.

```bash
python cfd_snapshot_adapters.py \
  monopole-3d raw_state.npz runtime_snapshot.npz

python cfd_runtime_monitor_postprocess.py \
  runtime_snapshot.npz --output-json summary.json

python cfd_runtime_failfast.py \
  summary.json --output-json verdict.json
```

### 3. Partial fast-screen for the current monopole solver

If you only have the stdout JSON log from `single_throat_monopole.py`:

```bash
python single_throat_monopole_jsonl_fastscreen.py \
  monopole.log --output-json monopole_verdict.json
```

This is the shortest current path to rejecting obviously bad `1/r` and `1/r^2`
behavior.

## Current Self-Test Results

### Step 34

- periodic exact-consistency snapshot:
  - `max |R_cont| = 1.4210854715202004e-14`
  - `max |R_Pois_exact| = 0.0`
  - `max |R_Pois_lin| = 0.0`
- Newton tail:
  - `Q_r` tail cv `= 0.00048311949183278685`
  - `mu_eff^2` tail median `= -0.003597967054167384`
- Yukawa tail:
  - `Q_r` tail cv `= 0.25427545643672567`
  - `mu_eff^2` tail median `= 1.954292666547206`

### Step 35

Classifier verdicts:

- Newton-like exterior: `PASS`
- Yukawa exterior: `FAIL`
- bad optics: `FAIL`
- projection-broken snapshot: `FAIL`
- missing optics snapshot: `INCOMPLETE`
- near-zero source snapshot: `INCOMPLETE`

### Step 36

Wavefunction adapter path:

- `rel jx error = 5.470539562816678e-05`
- `rel jy error = 5.470539562816674e-05`
- `rel jz error = 0.00020415416182304036`
- `max |R_cont| = 1.8245167713208398e-16`

Important caveat:

- the wavefunction adapter self-test validates current extraction and continuity
  closure
- it does **not** claim machine-zero exact Poisson closure
- that part depends on the chosen longitudinal reconstruction convention after
  current extraction

Projected 3D monopole adapter path:

- `max |R_cont| = 1.4210854715202004e-14`
- `max |R_Pois_exact| = 0.0`
- reconstructed source relative error `= 0.0`
- missing `W` without `--allow-uniform-W`: rejected
- reconstructed source without `lambda`: rejected

### Step 37

Monopole JSONL fast-screen:

- good log: `PASS`
- bad log: `FAIL`
- weak log: `INCOMPLETE`
- threshold-boundary log: `PASS`
- just-outside-threshold log: `FAIL`
- malformed JSONL log: `INCOMPLETE`

## What This Does Not Yet Do

1. It does not patch every live solver code path to emit snapshots. The
   `/var/projects/toy_physics/research/4d/scripts/single_throat_monopole.py`
   path now supports `--snapshot_out raw_state.npz`; the 4D wave path still
   needs its own export hook.
2. It does not prove the analog works.
3. It does not turn the partial monopole JSONL screen into an exact field-level
   falsifier.

## Next Practical Move

Add one export hook to whichever solver you want to test first:

- for the 4D GNLS path: emit `psi`, axes, `W`, and optionally `phi3`,
  `N_probe`, `Phi_eff`
- for the monopole path: emit `rho`, `mx`, `my`, `mz`, axes, and either
  direct `S_rho` or `W`, `Mdot`, `lambda`

Once that exists, the full runtime path is:

```bash
adapter -> postprocess -> failfast
```

and the simulation can be rejected or provisionally passed without manual note
reading.
