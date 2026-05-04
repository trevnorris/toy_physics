# Parent Throat Action — Snapshot Adapters

## Purpose

The runtime monitor and fail-fast classifier are only useful if solver output
can be translated into their snapshot schema with minimal friction.

This step adds exactly that translation layer in
[cfd_snapshot_adapters.py](/home/trevnorris/Downloads/em_projected/cfd_snapshot_adapters.py).

---

## 1. Supported upstream formats

The adapter currently supports two realistic solver-side families:

1. **4D wavefunction snapshots**
   - input: `psi` or `psi_real`/`psi_imag`, plus `x,y,z,w`
   - required for production: `W`
   - optional: `dWdw`, `dt_rho`, `phi3`, `rho0`, `N_probe`, `Phi_eff`
   - output schema: full 4D
     ```text
     rho, jx, jy, jz, jw, W, x, y, z, w
     ```

2. **3D monopole state dumps**
   - input: `rho, mx, my, mz, x, y, z`
   - plus either direct `S_rho` or enough information to reconstruct it
     from the monopole source law (`W`, `Mdot`, `lambda`)
   - output schema: projected 3D
     ```text
     rho_brane, Jx_brane, Jy_brane, Jz_brane, S_rho, x, y, z
     ```

The monitor now accepts either schema directly.

---

## 2. 4D wavefunction conversion

For a complex wavefunction \(\psi\), the adapter computes

\[
\rho=|\psi|^2,
\qquad
j_i=\frac{\hbar}{m}\Im(\psi^*\partial_i\psi),
\qquad
j_w=\frac{\hbar}{m}\Im(\psi^*\partial_w\psi).
\]

So a 4D GNLS snapshot can be converted into the exact fields needed by the
projected leakage monitor without inventing a separate current format.

If no projection weight `W` is provided, the adapter now rejects the snapshot by
default. A uniform \(w\)-weight can still be requested explicitly with
`--allow-uniform-W`, but that fallback is treated as a deliberate modeling
choice rather than a silent default.

---

## 3. 3D monopole conversion

The single-throat monopole solver already lives on the brane, so the adapter
does not fabricate a fake \(w\)-dimension.

Instead it maps

\[
\rho_{\rm brane}\leftarrow \rho,
\qquad
\mathbf J_{\rm brane}\leftarrow (m_x,m_y,m_z),
\]

and uses either:

- direct `S_rho`, or
- the built-in monopole source law
  \[
  S_\rho=-\dot M\,W+\lambda\dot M/V_{\rm domain}.
\]

That keeps the monitor tied to the solver’s actual source term.

For reconstructed sources, `lambda` / `lambda_bulk` is required. The adapter can
also consume an explicit `V_domain`; otherwise it computes the uniform cell
volume from the axes.

---

## 4. Self-test

[step_36_parent_throat_action_snapshot_adapter_sympy.py](/home/trevnorris/Downloads/em_projected/step_36_parent_throat_action_snapshot_adapter_sympy.py)
checks both paths:

1. **4D wavefunction path**
   - adapts a synthetic complex wave snapshot,
   - verifies that the projected extracted current matches the known phase
     gradient to small relative error,
   - passes the adapted snapshot through the runtime monitor and verifies
     machine-zero continuity closure.

2. **3D monopole path**
   - adapts a projected state dump,
   - passes it through the same runtime monitor,
   - verifies the same residual closure,
   - separately exercises the reconstructed monopole source law.

The self-test also checks these guardrails:

- a wavefunction snapshot without `W` is rejected unless
  `--allow-uniform-W` is used;
- using `--allow-uniform-W` produces a projection weight whose trapezoidal
  integral over the \(w\) axis is `1.0`;
- the uniform-`W` fallback is tested on a non-unit \(w\)-span, so dropping the
  `1/span` normalization factor would no longer be invisible;
- a reconstructed monopole source without `lambda` is rejected.

The wavefunction path intentionally does **not** claim machine-zero exact
Poisson closure in this self-test, because that depends on the separate
longitudinal reconstruction convention used after current extraction.

So the adapter layer is now tested, not just sketched.

---

## 5. CLI

```bash
python cfd_snapshot_adapters.py wavefunction-4d raw_wave.npz runtime_snapshot.npz
python cfd_snapshot_adapters.py monopole-3d raw_state.npz runtime_snapshot.npz
python cfd_runtime_monitor_postprocess.py runtime_snapshot.npz --output-json summary.json
python cfd_runtime_failfast.py summary.json --output-json verdict.json
```

This is the first end-to-end path from realistic solver-side dumps to a direct
falsification verdict.
