# Step 1 — Gauge-invariant circulation and bookkeeping

## Purpose

This step fixes the loop bookkeeping before any force law is attempted. The throat-mouth circulation label is kept separate from the electric charge branch, but this step does not claim a dynamical force result.

The 3D brane coordinates are
\[
\mathbf x=(x,y,z),
\]
and the hidden coordinate is \(w\). A finite mouth \(A\) has center \(\mathbf X_A\), radius \(R_A\), and a local orientation axis \(\hat s_A\). A small loop \(C_A\) links the mouth/rim circulation.

## Electric branch

The electric branch is
\[
\eta_Q=\pm1,\qquad q_*=\eta_Q e_*,\qquad e_*>0,
\]
with brane-normalized coupling
\[
q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}}.
\]
This branch label is not the circulation label.

## Gauge-invariant circulation constraint

With
\[
D_i=\partial_i-\frac{i q_*}{\hbar}A_i,
\]
the gauge-invariant phase-gradient is
\[
P_i\equiv \partial_i\theta-\frac{q_*}{\hbar}A_i.
\]
Under
\[
A_i\mapsto A_i+\partial_i\chi,
\qquad
\theta\mapsto \theta+\frac{q_*}{\hbar}\chi,
\]
we get
\[
P_i\mapsto P_i.
\]

The integer winding comes from single-valuedness of
\[
\psi=\sqrt{\rho}\,e^{i\theta}.
\]
Around a closed loop,
\[
\psi(2\pi)=\psi(0)
\quad\Longrightarrow\quad
e^{i\Delta\theta}=1.
\]
Equivalently,
\[
\cos\Delta\theta=1,
\qquad
\sin\Delta\theta=0,
\]
so
\[
\Delta\theta=2\pi n_A,\qquad n_A\in\mathbb Z.
\]

The fluxoid/winding integral is therefore
\[
\oint_{C_A}\left(\partial_i\theta-\frac{q_*}{\hbar}A_i\right)dl^i=2\pi n_A.
\]
The corresponding gauge-invariant hydrodynamic circulation is
\[
\Gamma_A
=\oint_{C_A}v_i dl^i
=\frac{\hbar}{m}\oint_{C_A}P_i dl^i
=\frac{2\pi\hbar}{m}n_A
\equiv \kappa n_A.
\]

For a single-valued gauge function \(\chi\),
\[
\oint_{C_A}\partial_i\chi\,dl^i=0.
\]
The SymPy audit checks this explicitly using
\[
\chi(t)=\chi_0+\chi_c\cos t+\chi_s\sin t.
\]
It also shows that a non-single-valued \(\chi_w t\) would have winding \(2\pi\chi_w\), so that case is not silently included as an ordinary small gauge transformation.

## Fixed data choice for the first pass

For the first derivation pass I hold \(n_A\) fixed as a topological/fluxoid label. Any physical current, throat throughput, or mixed-sector transport coefficient is introduced later as an additional closure.

That choice is conservative: the fluxoid constraint alone tells us the loop holonomy, not the radial force or mouth-throughput law.

## SymPy audit

The script verifies:

1. the gauge-invariant phase-gradient is invariant;
2. a single-valued \(\chi(t)\) has zero loop winding;
3. \(\psi(2\pi)=\psi(0)\) implies \(\Delta\theta\in2\pi\mathbb Z\);
4. in a pure-gauge loop representative, the \(\eta_Q\)-dependent pieces cancel before integration;
5. \(\Gamma_n=(2\pi\hbar/m)n\).

Run:

```bash
python step_01_fluxoid_firewall_sympy.py
```
