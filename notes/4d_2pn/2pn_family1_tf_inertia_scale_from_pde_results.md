# 2PN Family-1 TF inertia scale from the parent PDE

## What this step adds

The previous geometry-breathing reduction showed that the missing monopole wall closure
\[
\Delta K_{00}=\frac{109}{280}
\]
is supplied by the reduced geometry Hessian in \((a,L)\), and that its dynamics is an exact two-pole response. The remaining open item was the **bulk inertia scale**
\[
\rho_{\rm eff}
\]
appearing in the dimensionless frequency variable
\[
s=\omega^2\frac{\rho_{\rm eff}V_0 a_0^2}{\Sigma}.
\]

This step derives that scale from the **Family-1 parent PDE** itself, using the already-frozen \(n=5\) barotrope in the harmonic interior trap.

---

## 1. Parent-PDE reduction

Use the Family-1 interior throat branch with:

- steep radial wall and endcaps,
- harmonic interior brane trap
  \[
  V_{\rm in}(w)=\frac12 m_\psi \Omega_{\rm in}^2 w^2,
  \]
- and the frozen barotrope
  \[
  P(\rho)=K_{\rm EOS}\rho^n,\qquad
  h(\rho)=\frac{nK_{\rm EOS}}{n-1}\rho^{n-1}.
  \]

In the Thomas-Fermi / hydrostatic limit,
\[
h(\rho)+\frac12 m_\psi \Omega_{\rm in}^2 w^2=\mu_{\rm TF}.
\]

On the **filled-to-endcap** branch, the support reaches \(|w|=L_0/2\), so
\[
\mu_{\rm TF}=\frac12 m_\psi\Omega_{\rm in}^2\left(\frac{L_0}{2}\right)^2.
\]

That gives the exact axial profile
\[
\rho_0(w)=\rho_c\left(1-\frac{4w^2}{L_0^2}\right)^{\!1/(n-1)},
\qquad |w|\le \frac{L_0}{2},
\]
with central density
\[
\boxed{
\rho_c=
\left[
\frac{(n-1)m_\psi\Omega_{\rm in}^2L_0^2}{8nK_{\rm EOS}}
\right]^{\!1/(n-1)}.
}
\]

---

## 2. Exact bulk inertia scale

Write \(u=2w/L_0\) and define
\[
c_0(n)=\frac12\int_{-1}^{1}(1-u^2)^{1/(n-1)}\,du.
\]

Then the effective bulk inertia density entering the monopole breathing channel is
\[
\boxed{
\rho_{\rm eff}^{\rm TF}(n)=c_0(n)\,\rho_c.
}
\]

The exact closed form is
\[
\boxed{
c_0(n)=
\frac{\sqrt{\pi}\,\Gamma\!\left(1+\frac{1}{n-1}\right)}
     {2\,\Gamma\!\left(\frac32+\frac{1}{n-1}\right)}.
}
\]

For the already-frozen \(n=5\) branch,
\[
\boxed{
\rho_{\rm eff}^{\rm TF}(5)=
\frac{\sqrt{\pi}\,\Gamma(1/4)}{6\,\Gamma(3/4)}
\left(\frac{m_\psi\Omega_{\rm in}^2L_0^2}{10K_{\rm EOS}}\right)^{1/4}.
}
\]

Numerically,
\[
\frac{\sqrt{\pi}\Gamma(1/4)}{6\Gamma(3/4)}
\approx 0.87401918476404.
\]

So the Family-1 parent PDE fixes the previously free bulk inertia scale up to the already-physical throat parameters \((m_\psi,K_{\rm EOS},\Omega_{\rm in},L_0)\).

---

## 3. Exact inertia metric on the \((a,L)\) geometry sector

Using the same affine breathing ansatz as before,
\[
\boldsymbol\xi_\perp=\frac{\delta a}{a_0}\mathbf r,
\qquad
\xi_w=\frac{\delta L}{L_0}w,
\]
the reduced kinetic energy is
\[
T^{(2)}=
\frac12\,\rho_{\rm eff}^{\rm TF}V_0a_0^2\,\dot q^T\widehat M_{\rm TF}\dot q,
\qquad
q=\left(\frac{\delta a}{a_0},\frac{\delta L}{a_0}\right).
\]

The radial piece stays unchanged:
\[
\widehat M_{aa}=\frac35.
\]

The axial piece is fixed by
\[
c_2(n)=\frac12\int_{-1}^{1}u^2(1-u^2)^{1/(n-1)}\,du,
\]
with exact ratio
\[
\frac{c_2(n)}{c_0(n)}=\frac{n-1}{3n-1}.
\]

Therefore
\[
\boxed{
\widehat M_{LL}^{\rm TF}(n)=\frac{c_2}{4c_0}
=\frac{n-1}{4(3n-1)}.
}
\]

At \(n=5\),
\[
\boxed{
\widehat M_{\rm TF}(5)=
\begin{pmatrix}
\frac35 & 0\\[4pt]
0 & \frac1{14}
\end{pmatrix}.
}
\]

So the Family-1 \(n=5\) parent PDE renormalizes the axial breathing moment from the uniform-slice value
\[
\frac1{12}\longrightarrow \frac1{14}.
\]

This is a clean EOS-sensitive result.

---

## 4. Dynamic monopole response with the TF inertia branch

Carry forward the same reduced geometry Hessian \(\widehat H\) and volume-coupling vector \(\bar g\) from the previous geometry-Hessian step. The exact reduced response is still
\[
\boxed{
Y_{\rm geom}^{\rm TF}(s)=
\frac1\Sigma\,\bar g^T\left(\widehat H-s\widehat M_{\rm TF}(5)\right)^{-1}\bar g,
\qquad
s=\omega^2\frac{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}{\Sigma}.
}
\]

So the only change from the earlier affine-uniform model is:

- the free scale \(\rho_{\rm eff}\) is replaced by the explicit TF value,
- and the axial inertia coefficient becomes \(1/14\).

---

## 5. EM-worked point

At the same worked point used in the passed geometry-Hessian closure,
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac{1}{10},
\qquad
\beta=12,
\qquad
\Sigma_* \approx 0.2076143291835488854,
\]
the TF inertia branch gives dimensionless poles
\[
\boxed{
\lambda_- \approx 5.92556258,
\qquad
\lambda_+ \approx 237.91117494.
}
\]

The exact positive residues are
\[
\boxed{
R_- \approx 0.00262800,
\qquad
R_+ \approx 0.38665771,
}
\qquad
R_-+R_+=\frac{109}{280}.
\]

So the static closure is unchanged, as it must be.

The dominant residue fraction is
\[
\boxed{
\frac{R_+}{R_-+R_+}\approx 0.99324917.
}
\]

Thus the monopole wall completion remains an exact two-pole Stieltjes response with a highly accurate one-pole low-frequency reduction.

The Padé pole is
\[
\boxed{
\lambda_{\rm eff}^{\rm TF}\approx 188.17695898.
}
\]

On the low-frequency band \(0\le s\le 0.1\lambda_-\),
the max relative error of the one-pole reduction is only
\[
\boxed{
\max {\rm rel.\ err.}\approx 7.10\times 10^{-5}.
}
\]

---

## 6. Physical pole scales

Once the TF inertia scale is inserted, the physical poles are no longer free:
\[
\boxed{
\Omega_{\pm}^2=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}\,\lambda_{\pm},
\qquad
\Omega_{\rm eff}^2=
\frac{\Sigma_*}{\rho_{\rm eff}^{\rm TF}(5)V_0a_0^2}\,\lambda_{\rm eff}^{\rm TF}.
}
\]

Since \(\rho_{\rm eff}^{\rm TF}(5)\) is now explicit, the monopole breathing scale is fixed by the same Family-1 throat microphysics already used to define the parent PDE.

---

## Bottom line

This step closes an important remaining gap:

1. the Family-1 parent PDE fixes the bulk inertia scale,
2. the frozen \(n=5\) EOS fixes the axial breathing moment to \(1/14\),
3. the dynamic monopole wall channel remains a positive two-pole response,
4. and the one-pole breathing auxiliary stays an excellent low-frequency reduction.

So the bulk part of the monopole dynamics is no longer phenomenological.

The main open dynamic ingredient now is **not** the bulk inertia scale. It is the genuine **soft-wall / surface inertia completion**, if you want to go beyond the bulk TF branch.
