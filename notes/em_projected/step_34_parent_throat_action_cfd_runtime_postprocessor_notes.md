# Parent Throat Action — CFD Runtime Postprocessor

## Purpose

The symbolic monitor identities are no longer enough.

This step turns them into a concrete CFD-side postprocessor that can read a
snapshot, project the \(w\)-direction, reconstruct the longitudinal brane
potential, and compute the hard falsifier diagnostics directly.

The implementation lives in
[cfd_runtime_monitor_postprocess.py](/home/trevnorris/Downloads/em_projected/cfd_runtime_monitor_postprocess.py),
and this step self-tests it on synthetic exact fields.

---

## 1. Snapshot schema

The postprocessor consumes a `.npz` snapshot with required keys

```text
rho, jx, jy, jz, jw, W, x, y, z, w
```

where:

- `rho`, `jx`, `jy`, `jz`, `jw` have shape `(nx, ny, nz, nw)`,
- `W` and `w` have length `nw`,
- `x`, `y`, `z` are the three brane axes.

Optional keys:

```text
dWdw, dt_rho, phi3, rho0, N_probe, Phi_eff
```

- `dWdw` overrides the numerical derivative of `W`.
- `dt_rho` lets the continuity residual use an explicit time derivative.
- `phi3` bypasses the FFT Helmholtz reconstruction if the simulation already
  provides the longitudinal potential.
- `rho0`, `N_probe`, `Phi_eff` enable the linearized monitor and optics fit.

---

## 2. Computed diagnostics

The tool projects the snapshot to the brane:

\[
\rho_{\rm brane}=\int W\rho\,dw,\qquad
\mathbf J_{\rm brane}=\int W\mathbf J_{xyz}\,dw,
\]

and computes the exact leakage source

\[
S_\rho
=
-[WJ_w]_{-\infty}^{+\infty}
+\int W'J_w\,dw.
\]

It then evaluates:

1. continuity residual
   \[
   R_{\rm cont}=\partial_t\rho_{\rm brane}+\nabla\!\cdot\mathbf J_{\rm brane}-S_\rho,
   \]
2. exact longitudinal / Poisson residual
   \[
   R_{\rm Pois}^{\rm exact}
   =
   \rho_{\rm brane}\nabla^2\phi_3-S_\rho+\partial_t\rho_{\rm brane}
   +(\nabla\rho_{\rm brane})\cdot\mathbf v_{\rm brane},
   \]
3. linearized quasi-static residual
   \[
   R_{\rm Pois}^{\rm lin}=\rho_0\nabla^2\phi_3-S_\rho,
   \]
4. exterior flux plateau
   \[
   Q_r(r)=4\pi r^2\partial_r\phi_3,
   \]
5. effective exterior mass scale
   \[
   \mu_{\rm eff}^2(r)=\frac{\nabla^2\phi_3-S_\rho/\rho_0}{\phi_3},
   \]
6. optics coefficient fit
   \[
   \alpha_{\rm fit}(r)=-\frac{c_{\rm probe}^2\,[N_{\rm probe}(r)-1]}{\Phi_{\rm eff}(r)}.
   \]

So the postprocessor computes exactly the runtime quantities needed for the fast
falsification suite.

---

## 3. Synthetic self-test

The self-test script
[step_34_parent_throat_action_cfd_runtime_postprocessor_sympy.py](/home/trevnorris/Downloads/em_projected/step_34_parent_throat_action_cfd_runtime_postprocessor_sympy.py)
runs three checks:

1. **Periodic exact-consistency snapshot**
   - constructs a periodic harmonic \(\phi_3\),
   - defines \(S_\rho=\rho_0\nabla^2\phi_3\),
   - encodes that source exactly through the boundary leakage term,
   - verifies the postprocessor returns
     `max |R_cont|`, `max |R_Pois_exact|`, and `max |R_Pois_lin|`
     at machine-zero scale.

2. **Newton exterior calibration**
   - feeds a softened \(1/r\) field,
   - checks that the exterior `Q_r` shell average is flatter than the Yukawa
     case,
   - checks that the exterior `mu_eff^2` tail stays near zero.

3. **Yukawa exterior calibration**
   - feeds a softened \(e^{-\mu r}/r\) field,
   - checks that `mu_eff^2` rises to a positive exterior value,
   - checks that the `Q_r` shell average is less plateau-like.

The optics side is tested with the exact \(n=5\) target profile
\[
N_{\rm probe}=1-\frac{2\Phi_{\rm eff}}{c_{\rm probe}^2},
\]
so the recovered \(\alpha_{\rm fit}\) sits near \(2\).

---

## 4. CLI

The real-simulation entry point is:

```bash
python cfd_runtime_monitor_postprocess.py snapshot.npz --output-json summary.json
```

That summary now gives the exact quantities needed to falsify the analog on
runtime data:

- continuity/projection failure,
- exact or linearized Poisson failure,
- Yukawa-like screening,
- wrong weak-field optics coefficient.

This is the first point where the Branch-B patch can be tested directly against
simulation output instead of against reduced symbolic surrogates.
