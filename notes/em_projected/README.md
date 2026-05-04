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
- `cfd_source_convention_sweep.py`

## Current Project State

We are pausing the computational-physics work here and carrying the derivation
side into the paper.

The derivation/audit chain is in good enough shape to use as the source for a
complete paper derivation, subject to normal paper-level cleanup and notation
unification. The computational tools below should be treated as a falsification
scaffold, not as established simulation evidence for the paper.

Current interpretation:

1. The reduced Branch-B derivation fixes the dimensionless target
   `mhat0^2 * P0 = 54/5`.
2. The physical dimensional bridge is still carried by `S_port`.
3. The direct-SI electron map with `S_port = 1` is conditionally falsified by
   `step_31`.
4. `c_s` and `c` are decoupled model speeds; there is no derivation-side rule
   forcing `c_s <= c`.
5. The first monopole runtime data are not admissible falsification data. The
   current simulator does not yet provide a stationary, non-reflecting,
   mass-balanced throat configuration from which clean exterior field readings
   can be taken.
6. Prior quiet-simulation work in `/var/projects/4d_1pn_sim/` also did not
   produce a usable clean-readout simulator, despite substantial effort.
   Simulation scripts under research folders should therefore be treated as
   legacy or diagnostic experiments, not as evidence-producing infrastructure.

So the paper should use the derivation results and the conditional electron-map
falsification. It should not claim that the CFD realization has passed or
failed.

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

Optional resolution-derived `Q_r` threshold:

- `--q-tail-cv-max-from-resolution` sets
  `q_tail_cv_max = sqrt(2 / max(tail_n_points, 4))`
- `tail_n_points` is emitted by the runtime monitor from the number of exterior
  shell samples used in the `Q_r` tail
- heuristic derivation: near-stationary tail shells are treated as independent
  samples, so the sampling noise of a mean/CV-type tail statistic scales like
  `1/sqrt(n)`; applying a factor-of-2 safety margin gives `sqrt(2/n)`, with
  `max(n, 4)` preventing overconfident thresholds for very short tails
- `--q-tail-cv-max` and `--q-tail-cv-max-from-resolution` are mutually
  exclusive; the CLI raises instead of silently choosing one

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

#### `cfd_source_convention_sweep.py`

This is a diagnostic companion to the runtime monitor. It does not change a
verdict. It asks whether a failing projected snapshot can be rescued by a simple
schema convention change:

- stored `Jx/Jy/Jz` as currents vs. stored velocities;
- source scale;
- source scale plus a uniform offset;
- `phi3` sign flip.

Run it on the adapted runtime snapshot:

```bash
python cfd_source_convention_sweep.py runtime_snapshot.npz
```

Reading guide:

- if the velocity interpretation sharply improves the residuals, the adapter
  probably mapped velocities as currents;
- if the best source scale is near `-1`, suspect a source sign convention;
- if the best source scale is a constant near but not equal to `1`, suspect a
  source-normalization mismatch;
- if no convention drops the residuals near the fail-fast threshold, treat the
  failure as a real solver/projection issue.

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

python cfd_source_convention_sweep.py \
  runtime_snapshot.npz
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

python cfd_source_convention_sweep.py \
  runtime_snapshot.npz
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
- default source-free exterior: `INCOMPLETE`
- fixed `Q_r` threshold boundary: `PASS`
- fixed `Q_r` just outside threshold: `FAIL`
- `mu_eff2` threshold boundary: `PASS`
- `mu_eff2` just outside threshold: `FAIL`
- resolution-derived `Q_r` threshold boundary at `tail_n_points=200`: `PASS`
- resolution-derived `Q_r` just outside threshold: `FAIL`

### Step 36

Wavefunction adapter path:

- `rel jx error = 5.470539562816678e-05`
- `rel jy error = 5.470539562816674e-05`
- `rel jz error = 0.00020415416182304036`
- `max |R_cont| = 1.8245167713208398e-16`
- uniform `W` fallback is checked on a non-unit `w` span:
  - integral error `= 0.0`
  - value error `= 0.0`

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

## Monopole Smoke-Test Output

These are contaminated simulator diagnostics only. They are included here so the
next simulation pass starts from the actual observed state, not because they are
valid model evidence.

The existing monopole setup is not a clean readout simulator:

1. inserting a throat drives large transient sloshing;
2. the throat drains the superfluid, so the background density drifts over time;
3. waves reflect off the domain edges and return into the particle region;
4. there is no settled measurement window where exterior field readings can be
   trusted.

Therefore the outputs below should be read as implementation smoke tests only.
They do not pass or falsify the analog model.

This is stronger than a single bad run. Prior attempts to build even a quiet
simulation in `/var/projects/4d_1pn_sim/` did not yield a reliable readout
setup. Any simulation scripts under research directories are therefore lower
confidence still: useful for inspecting implementation ideas, but not useful
for model evidence until replaced by a purpose-built clean-readout simulator.

### N256 / 20k JSONL fast-screen

For the `N=256`, `20000`-step monopole run, the stdout fast-screen returned:

- status: `PASS`
- `mach_max = 0.008957277052104473`
- `dP_slope = -0.9590923887813703`
- `geff_slope = -2.1394193895285403`
- `dP_npts = 93`
- `geff_npts = 76`

This only says the late diagnostic stream can sometimes fit approximate `1/r`
pressure and `1/r^2` acceleration slopes in a contaminated run. It is not
evidence that a clean stationary exterior field has formed.

### N256 snapshot monitor

On the adapted snapshot from the same family of runs, the stricter runtime
monitor returned:

- `max |R_cont| = 0.05680320603942618`
- `max |R_Pois_exact| = 0.09229334860573492`
- `max |R_Pois_lin| = 0.09053007745399924`
- `continuity_rel = 0.1303949052528477`
- `poisson_rel = 0.21620849389895325`
- `Q_r_tail_cv = 0.16723952034586026`
- `mu_eff2_tail_median = -0.0028201693218378662`
- status: `FAIL`

A second tail-window probe gave `Q_r_tail_cv = 0.2644081266725705`. In a clean
solver this would be a serious flux-plateau failure; in this run it is better
treated as evidence that the current simulator is not yet suitable for
readout-level falsification.

### Source convention sweep

The source-convention sweep on that adapted snapshot reported:

- schema: `projected_3d`
- source RMS: `0.0005039559150496455`
- mean density: `0.9999975516504094`
- `phi3` present: `False`
- periodic `x/y/z`: `True`

Treating stored `Jx/Jy/Jz` as brane currents gave:

- continuity relative residual: `0.1303949052528477`
- exact Poisson relative residual: `0.21620849389895322`
- linear Poisson relative residual: `0.21333244148302308`
- best continuity source scale: `0.8954031495849052`
- best continuity residual after that scale: `0.07786096710895239`
- best exact Poisson source scale: `0.8240759659735598`
- best exact Poisson residual after that scale: `0.12568550865520497`

Treating stored `Jx/Jy/Jz` as velocities made the residuals slightly worse.
Uniform offsets and `phi3` sign flips did not rescue the run.

Interpretation:

1. the adapter's current-vs-velocity convention is probably correct;
2. the source normalization appears roughly 10-20% mismatched, depending on the
   residual being fit;
3. that mismatch is not enough to explain the failure;
4. because `phi3` is absent, the monitor reconstructs a velocity potential, not
   a solver-emitted pressure/enthalpy potential.

The current computational state is therefore: **the available monopole run is a
contaminated smoke test, not a valid pass/fail falsification run**.

## What This Does Not Yet Do

1. It does not patch every live solver code path to emit snapshots. The
   `/var/projects/toy_physics/research/4d/scripts/single_throat_monopole.py`
   path now supports `--snapshot_out raw_state.npz`; the 4D wave path still
   needs its own export hook.
2. It does not prove the analog works.
3. It does not turn the partial monopole JSONL screen into an exact field-level
   falsifier.
4. It does not support arbitrary meshes. The runtime monitor assumes uniform
   Cartesian `x/y/z` axes.
5. It does not provide a nonperiodic Helmholtz reconstruction. If `--nonperiodic`
   is used, the snapshot must already include `phi3`.
6. It does not validate physical units. The tools assume the solver snapshot has
   already been nondimensionalized consistently with the chosen port map.
7. It does not infer a trustworthy background density in every case. If `rho0`
   is missing, the monitor uses the mean projected density as a fallback.
8. It does not make `R_Pois_lin` a decisive verdict channel. The fail-fast
   classifier uses the exact Poisson residual; the linearized residual is a
   diagnostic readout.
9. It does not convert per-cell optical checks into an independent ray-traced
   bending or Shapiro test. `alpha_fit` is a local coefficient diagnostic.
10. It does not guarantee a minimum grid resolution or tail-window quality. The
    default thresholds are fail-fast rules of thumb, not convergence criteria;
    the optional resolution-derived `Q_r` threshold only ties one tolerance to
    the number of tail-shell samples.
11. It does not use identical data-quality fields in every tool. The fail-fast
    classifier and JSONL screen both report incomplete reasons, but their input
    schemas differ.
12. It does not make a missing or flat projection weight physical. `--allow-uniform-W`
    is an explicit fallback for testing, not a derived brane projection.
13. It does not treat `1/L_w` or the flat `W` fallback as a physical projection
    law. Those are normalization placeholders until the solver emits a derived
    projection weight.
14. It does not claim exact Poisson closure for the wavefunction adapter path;
    that depends on a separate longitudinal reconstruction convention.
15. It does not yet test the monopole solver against a solver-emitted
    pressure/enthalpy potential. The current adapted monopole snapshot lacks
    `phi3`, so the monitor reconstructs a velocity potential.
16. It does not make a JSONL fast-screen `PASS` decisive. The N256 run currently
    passes the slope screen but is contaminated by sloshing, density drift, and
    boundary reflections.
17. It does not resolve whether the remaining continuity/source mismatch is a
    numerical-resolution issue, a boundary/projection issue, or a missing solver
    channel. The convention sweep only rules out the easiest schema mistakes.
18. It does not yet include the simulator features required for clean readout:
    mass balance, adiabatic throat ramp-in, absorbing boundaries, and a
    separated settling/readout window.
19. It does not inherit a working quiet-simulation baseline from
    `/var/projects/4d_1pn_sim/`; that effort did not reach usable readout
    quality.

## Next Practical Move

Computational work is paused. When it resumes, do not jump straight to a larger
GPU run or reuse research-folder simulation scripts as if they were validated.
First build a purpose-designed simulator that can produce clean readings:

1. add a mass-balanced throat/source configuration so the mean density does not
   drain away;
2. ramp the throat in adiabatically instead of inserting it impulsively;
3. add absorbing or sponge boundaries so reflected waves do not re-enter the
   readout region;
4. separate the transient settling phase from the measurement window;
5. emit the full fields needed for readout:
   - for the 4D GNLS path: `psi`, axes, `W`, and optionally `phi3`, `N_probe`,
     `Phi_eff`;
   - for the monopole path: `rho`, `mx`, `my`, `mz`, axes, direct `S_rho` or
     `W`, `Mdot`, `lambda`, plus the pressure or enthalpy channel needed to
     test the intended longitudinal field.

Once that exists, the full runtime path is:

```bash
adapter -> postprocess -> failfast
```

The derivation-side next move is separate: move the audited step chain into the
paper and keep the computational results clearly marked as preliminary
falsification scaffolding.
