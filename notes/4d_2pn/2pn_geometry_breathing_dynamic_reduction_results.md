# 2PN geometry-breathing dynamic reduction: current result

## What this step adds

The previous step showed that the missing raw monopole wall closure
\[
\Delta K_{00}=\frac{109}{280}
\]
can be generated from the reduced static geometry Hessian in \((a,L)\), but it still left the **pole**
\[
Y_{\rm mono}(\omega)=\frac{109/280}{1-\omega^2/\Omega_{\rm mono}^2}
\]
underdetermined.

This step supplies the missing **dynamic reduction**.

Using a minimal affine breathing ansatz for the throat geometry,
we derive an exact reduced inertia matrix, obtain the full **two-pole** monopole response, and then show that the familiar single-pole auxiliary is the controlled low-frequency Padé reduction of that exact 2DOF dynamics.

This step does **not** depend on the earlier corrected \(z^8\) isotropic 4D-ball quotient, because it uses only the already-passed static geometry data.

---

## 1. Minimal affine inertia model

Use dimensionless geometry coordinates
\[
q=\left(\frac{\delta a}{a_0},\frac{\delta L}{a_0}\right),
\]
and the affine displacement field inside the 4D cylinder-like throat:
\[
\boldsymbol\xi_\perp=\frac{\delta a}{a_0}\,\mathbf r, 
\qquad
\xi_w=\frac{\delta L}{L_0}\,w.
\]

For the 3-ball cross section,
\[
\int_{B^3_{a_0}} r^2\,d^3x = \frac{4\pi a_0^5}{5},
\]
and along the axial interval,
\[
\int_{-L_0/2}^{L_0/2} w^2\,dw = \frac{L_0^3}{12}.
\]

So the reduced kinetic energy is
\[
T^{(2)} = \frac12\rho_{\rm eff}V_0 a_0^2\,\dot q^T \widehat M\,\dot q,
\qquad
\widehat M=
\begin{pmatrix}
\frac35 & 0 \\
0 & \frac1{12}
\end{pmatrix},
\]
with
\[
V_0 = \frac{4\pi}{3}a_0^3L_0 = \frac{4\pi}{3}\Lambda a_0^4.
\]

This is the minimal entrained-fluid inertia metric for the \((a,L)\) geometry sector.

---

## 2. Exact dynamic monopole susceptibility

From the previous static geometry step,
\[
E^{(2)} = \frac12\Sigma\,q^T \widehat H\,q,
\]
with
\[
\widehat H = \frac{a_0^2}{\Sigma}H_0,
\qquad
\bar g = \frac{a_0\nabla V_0}{V_0} = \left(3,\frac1\Lambda\right).
\]

Define the scaled frequency variable
\[
s = \omega^2\,\frac{\rho_{\rm eff}V_0 a_0^2}{\Sigma}.
\]
Then the exact reduced dynamic monopole response is
\[
\boxed{
Y_{\rm geom}(s)=
\frac1\Sigma\,
\bar g^T(\widehat H-s\widehat M)^{-1}\bar g.
}
\]

The raw wall monopole completion is therefore
\[
\boxed{
K_{00}^{\rm raw}(s)
=
-\frac{757}{2520}+Y_{\rm geom}(s).
}
\]
At \(s=0\), this reproduces the static completion from the earlier step.

---

## 3. Exact low-frequency reduction

Expanding at small \(s\),
\[
Y_{\rm geom}(s)=\Delta_0 + \Delta_2 s + O(s^2),
\]
with
\[
\Delta_0
=
\frac1\Sigma\bar g^T\widehat H^{-1}\bar g,
\qquad
\Delta_2
=
\frac1\Sigma\bar g^T\widehat H^{-1}\widehat M\widehat H^{-1}\bar g.
\]

The associated \([1/1]\) Padé single-pole reduction is
\[
\boxed{
Y_{\rm geom}(s)\approx
\frac{\Delta_0}{1-s/\lambda_{\rm eff}},
\qquad
\lambda_{\rm eff}=\frac{\Delta_0}{\Delta_2}.
}
\]

So the old monopole auxiliary really is the low-frequency reduction of a genuine reduced geometry dynamics.

---

## 4. EM-branch worked point

Using the same worked point as the previous geometry-Hessian closure,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac1{10},
\qquad
\beta=12,
\qquad
\Sigma_*=0.2076143291835488854\ldots,
\]
we get
\[
\widehat H \approx
\begin{pmatrix}
114.33174685 & 19.35786733 \\
19.35786733 & 3.80598758
\end{pmatrix},
\qquad
\widehat M=
\begin{pmatrix}
0.6 & 0 \\
0 & 1/12
\end{pmatrix},
\qquad
\bar g=\left(3,\frac1\Lambda\right).
\]

The generalized eigenproblem
\[
\widehat H v_i = \lambda_i \widehat M v_i
\]
gives the exact dimensionless pole parameters
\[
\boxed{
\lambda_- \approx 5.23115613,
\qquad
\lambda_+ \approx 230.99360624.
}
\]
Hence the physical pole frequencies are
\[
\boxed{
\Omega_{\pm}^2
=
\frac{\Sigma_*}{\rho_{\rm eff}V_0 a_0^2}\,\lambda_{\pm}.
}
\]
For \(a_0=1\) and \(\rho_{\rm eff}=1\),
\[
\Omega_-^2\approx 0.14034117,
\qquad
\Omega_+^2\approx 6.19708399.
\]

---

## 5. Exact two-pole decomposition

At the worked point, the exact dynamic response can be written as
\[
\boxed{
Y_{\rm geom}(s)=
\frac{R_-}{1-s/\lambda_-}
+
\frac{R_+}{1-s/\lambda_+},
}
\]
with positive static residues
\[
\boxed{
R_- \approx 0.00327376153,
\qquad
R_+ \approx 0.38601195275.
}
\]
They sum exactly to
\[
R_-+R_+ = \frac{109}{280}.
\]

So the exact geometry breathing response is a **two-pole Stieltjes function with positive residues**.

The dominant residue fraction is
\[
\frac{R_+}{R_-+R_+} \approx 99.1590\%.
\]
So the old single-pole auxiliary is not just qualitatively justified; it is quantitatively a very good reduction of the exact 2DOF geometry dynamics.

---

## 6. Single-pole Padé reduction

At the same worked point,
\[
\boxed{
\lambda_{\rm eff}\approx 169.48205088,
\qquad
\Omega_{\rm eff}^2
=
\frac{\Sigma_*}{\rho_{\rm eff}V_0 a_0^2}\,\lambda_{\rm eff}.
}
\]
For \(a_0=1\) and \(\rho_{\rm eff}=1\),
\[
\Omega_{\rm eff}^2 \approx 4.54685531.
\]

On the low-frequency band
\[
0\le s\le 0.1\lambda_-,
\]
the maximum relative error of the one-pole Padé form against the exact two-pole response is
\[
\boxed{\max |\delta Y/Y| \approx 8.87\times 10^{-5}.}
\]

So for PN/local response purposes, the one-pole monopole auxiliary is extremely accurate on the natural low-frequency band.

---

## 7. Interpretation

This is the cleanest dynamic monopole picture so far.

The raw monopole wall sector is now
\[
\boxed{
K_{00}^{\rm raw}(s)
=
-\frac{757}{2520}
+
\frac{R_-}{1-s/\lambda_-}
+
\frac{R_+}{1-s/\lambda_+},
}
\]
with the controlled low-frequency reduction
\[
\boxed{
K_{00}^{\rm raw}(s)
\approx
-\frac{757}{2520}+
\frac{109/280}{1-s/\lambda_{\rm eff}}.
}
\]

So the earlier “global breathing auxiliary” is no longer just a bookkeeping closure. It is the low-frequency reduction of the same \((a,L)\) geometry sector that already produced the exact static \(109/280\) coefficient.

---

## 8. What remains

The remaining task is now narrow and concrete:

1. derive the overall inertial scale \(\rho_{\rm eff}\) (or its more accurate soft-wall analog) from the actual Family-1 confinement / traction PDE,
2. and, if desired, improve the affine inertia ansatz to a soft-wall boundary-layer inertia.

Once that is done, the monopole pole is fully fixed from throat-side physics rather than inferred from operator matching.
