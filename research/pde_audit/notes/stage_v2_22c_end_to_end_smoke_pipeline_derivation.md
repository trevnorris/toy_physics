# Stage V2-22C — End-to-end branch-realization smoke pipeline and tolerance budget

## Purpose

V2-22C closes the executable handoff chain:

\[
\text{solver export}
\rightarrow
\text{V2-22B validation}
\rightarrow
\text{V2-22A profile adapter}
\rightarrow
\text{V2-21 grouped observable extraction}
\rightarrow
\text{residual/tolerance packet}.
\]

This stage does **not** claim a solved moving-throat branch. Its purpose is to prove that a future PDE branch can be tested without changing the definitions after seeing the residuals.

The resulting executable stage distinguishes three classes of failure:

1. **handoff/validation failure**, such as a hard-cap geometry or nonpositive mixed-block determinant;
2. **open-throat/stability failure**, such as `D0 <= 0` or `Delta <= 0`;
3. **target-realization failure**, such as nonzero `R_norm`, `R_pole`, `R_P2`, or `R_P4`.

A stable branch is allowed to fail the target packet. That is a useful falsification, not a pipeline error.

---

## Inputs

The stage reuses the existing executable tools:

- `stage_v2_22b_solver_handoff_validator.py`
- `stage_v2_22a_profile_to_coefficient_adapter.py`
- `stage_v2_21_branch_extraction_fixture.py`

The built-in solver export is the open finite-throat D/N prototype:

\[
R_{\rm exit}>0,
\qquad
\texttt{boundary\_class}=\texttt{open\_impedance}.
\]

The invalid control is the forbidden hard-cap packet:

\[
R_{\rm exit}=0,
\qquad
\texttt{boundary\_class}=\texttt{hard\_cap}.
\]

---

## Formula checks performed in V2-22C

V2-22C runs a lightweight orchestration audit. It intentionally does not rerun the heavier component SymPy audits, because this stage is a handoff pipeline rather than another symbolic compiler.

The checked identities are:

### 1. Grouped inverse map

For

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}2,
\]

the inverse is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

### 2. Weak-axisymmetric grouped signature

For

\[
\lambda=(1,\tfrac12,-1),
\]

one has

\[
\bar x=x_0,
\qquad
b_x=3a_x.
\]

### 3. Constant-prefactor branch

The constant-prefactor conditions are

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
\]

Substitution gives

\[
P_2=0,
\qquad
P_4=0.
\]

### 4. Quadrupole normalization equivalence

The target

\[
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}
\]

when

\[
\gamma_{\rm quad}^{\rm eff}
=\widehat m_0^{\,2}\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}.
\]

The orchestration audit passed all `8/8` formula checks.

---

## Pipeline gates

### V2-22B validation gates

The solver export must satisfy:

\[
R_{\rm exit}>0,
\]

\[
\texttt{boundary\_class}=\texttt{open\_impedance},
\]

\[
\texttt{pre\_target\_freeze}=\texttt{true},
\qquad
\texttt{target\_blind}=\texttt{true}.
\]

The mixed block must satisfy

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0.
\]

The valid open-throat sample passed. The invalid hard-cap sample was rejected.

### V2-21 extraction gates

For each grouped lane, the extracted moments are

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

The normalized response and outgoing prefactor are

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

The target residual packet is

\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}=\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,
\qquad
R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

---

## Smoke-branch output

The open D/N smoke branch passed validation and stability, but it did not hit the target surface:

```text
open_gate_pass: True
stability_gate_pass: True
target_packet_pass: False
D0_bar: 1.560684600911219
N0_bar: 0.0008248594058984478
P0_bar: 0.0005285240883499759
P0_target: 10.8
R_pole: -2.452467306741504
R_norm: -10.79947147591165
R_P2: 0.0006833291043178326
R_P4: 0.0006221198975181069
R_tail: 0
```

This is the intended behavior. The sample branch is a stable open-throat profile fixture, not a calibrated realization of the universal quadrupole target.

---

## Calibrated direct-coefficient control

The calibrated direct-coefficient control passed the full target packet:

```text
target_packet_pass: True
R_pole: 0
R_norm: 0
R_P2: 0
R_P4: 1.4210854715202e-14
```

The tiny nonzero `R_P4` is floating-point roundoff.

This control proves that the gate can pass when the target surface is actually satisfied. Therefore the smoke branch failure is not caused by a broken extraction pipeline.

---

## Tolerance budget

The default extraction tolerance is

\[
\varepsilon_{\rm extract}=10^{-9}.
\]

The validator profile normalization tolerance is

\[
\varepsilon_{\rm norm}=5\times10^{-3}.
\]

The mixed-block determinant tolerance is

\[
\varepsilon_\Delta=10^{-12}.
\]

The smoke branch’s large normalized residual is dominated by `R_norm`, not by numerical roundoff:

\[
\left|\frac{R_{\rm norm}}{P_0^{\rm target}}\right|
\approx
0.999951.
\]

So the branch is very far below the required outgoing normalization.

---

## Final V2-22C conclusion

The executable branch-realization chain is now closed:

\[
\boxed{
\text{frozen solver packet}
\to
\text{validated profile manifest}
\to
\text{V2-21 observable packet}
\to
\text{residual/tolerance report}.
}
\]

The mechanical pipeline passes:

```text
mechanical_pipeline_pass: True
```

But no PDE branch realization is claimed:

```text
branch_target_realization_claimed: False
```

The next actual science step is to feed this pipeline a real moving-throat PDE branch export and see whether the resulting branch, without post-target refitting, lands on the isotropic full-bundle target surface.
