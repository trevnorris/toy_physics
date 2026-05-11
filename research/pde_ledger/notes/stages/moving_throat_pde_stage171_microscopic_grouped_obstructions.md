# Moving-Throat PDE — Stage 171: Exact Microscopic Decomposition of the Linear Grouped Outlet Obstructions

## Purpose

Stage 238 proved that the entire **linear** grouped real `P2` outlet problem collapses to the pair
\[
\mathcal K_A:=\delta D_{A,2}+\frac19\,\delta D_{A,0},
\qquad
\mathcal G_A:=\delta N_{A,0}-P_0\,\delta D_{A,0},
\qquad
A\in\{20,21,22\}.
\]
Those two combinations feed the hidden outlet deformations
\[
\delta\kappa_W^{(A)}=
\frac{3(1-\sigma_*)}{\sigma_* D_0}\,\mathcal K_A,
\qquad
\delta\gamma_W^{(A)}=
-\frac{1-\sigma_*}{9\sigma_* N_0}\,\mathcal G_A.
\]

So the next honest theorem gate is no longer to derive another abstract grouped map.
It is to compute the actual microscopic grouped-lane anisotropies that feed
\(\mathcal K_A\) and \(\mathcal G_A\).

This stage performs that decomposition.

The main result is that the full linear grouped-anisotropy problem collapses to **four microscopic defect bundles**:
\[
\boxed{
\mathcal K_A = \mathcal W_A-\mathcal B_A-\mathcal Z_A,
}
\]
with
\[
\mathcal W_A:=\frac19\,\delta K_A-\delta M_A,
\qquad
\mathcal B_A:=\delta B_{A,2}+\frac19\,\delta B_{A,0},
\qquad
\mathcal Z_A:=\delta Z_{A,2}+\frac19\,\delta Z_{A,0},
\]
and
\[
\boxed{
\mathcal G_A=
-P_0\,\delta K_A
+P_0\,\delta B_{A,0}
+\mathcal N_A,
}
\]
with
\[
\mathcal N_A:=\delta N_{A,0}+P_0\,\delta Z_{A,0}.
\]

So the linear grouped-even and grouped-odd outlet defects are not sourced by the full microscopic bundle independently. They are driven only by:

1. a wall baseline anisotropy \(\mathcal W_A\),
2. a BdG support anisotropy \(\mathcal B_A\),
3. a conservative Maxwell/mixed anisotropy \(\mathcal Z_A\),
4. and an outgoing-transfer bundle \(\mathcal N_A\).

That is a much sharper microscopic theorem gate than the raw grouped variables of Stage 238.

---

## 1. Carry-forward full-bundle coefficients

From the full grouped real `P2` bundle of Stage 6, each lane has low-frequency data
\[
D_A^{\rm(cons)}(\omega)=D_{A,0}+D_{A,2}\,\omega^2+D_{A,4}\,\omega^4+O(\omega^6),
\]
with
\[
D_{A,0}=K_A-B_{A,0}-Z_{A,0},
\]
\[
D_{A,2}=-\bigl(M_A+B_{A,2}+Z_{A,2}\bigr),
\]
\[
D_{A,4}=-\bigl(B_{A,4}+Z_{A,4}\bigr),
\]
and outgoing-transfer coefficient
\[
N_{A,0}=N_0+\delta N_{A,0}.
\]
Here:

- \(K_A,M_A\) are the grouped wall baseline coefficients,
- \(B_{A,n}\) are the grouped BdG support moments,
- \(Z_{A,n}\) are the grouped conservative Maxwell/mixed moments,
- \(N_{A,0}\) is the grouped static outgoing-transfer moment.

Linearizing around an isotropic compensated branch,
\[
D_{A,n}=D_n+\delta D_{A,n},
\qquad
N_{A,0}=N_0+\delta N_{A,0},
\]
Stage 238 already gave
\[
\mathcal K_A=\delta D_{A,2}+\frac19\,\delta D_{A,0},
\qquad
\mathcal G_A=\delta N_{A,0}-P_0\,\delta D_{A,0}.
\]

Substituting the full-bundle decomposition immediately gives
\[
\delta D_{A,0}=\delta K_A-\delta B_{A,0}-\delta Z_{A,0},
\]
\[
\delta D_{A,2}=-\delta M_A-\delta B_{A,2}-\delta Z_{A,2}.
\]
So
\[
\boxed{
\mathcal K_A
=
\left(\frac19\,\delta K_A-\delta M_A\right)
-
\left(\delta B_{A,2}+\frac19\,\delta B_{A,0}\right)
-
\left(\delta Z_{A,2}+\frac19\,\delta Z_{A,0}\right).
}
\]
This is the exact origin of the three pieces \(\mathcal W_A,\mathcal B_A,\mathcal Z_A\) quoted above.

Likewise,
\[
\boxed{
\mathcal G_A
=
\delta N_{A,0}
-P_0\bigl(\delta K_A-\delta B_{A,0}-\delta Z_{A,0}\bigr)
=
-P_0\,\delta K_A+P_0\,\delta B_{A,0}+\mathcal N_A.
}
\]
So the odd grouped-normalization defect has two distinct microscopic sources:

1. **induced odd anisotropy** coming from static conservative lane splitting, and
2. **direct transfer anisotropy** coming from the outgoing Maxwell/mixed bundle.

That split is one of the main conceptual results of the stage.

---

## 2. Exact BdG contribution to the outlet obstructions

For each grouped lane \(A\), the BdG moments are
\[
B_{A,0}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\qquad
B_{A,2}=\sum_\alpha \frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4}.
\]

Linearize around an isotropic branch
\[
c_{A\alpha}=c_\alpha+\delta c_{A\alpha},
\qquad
\varpi_{A\alpha}=\varpi_\alpha+\delta\varpi_{A\alpha}.
\]
Then
\[
\boxed{
\delta B_{A,0}
=
\sum_\alpha
\left(
\frac{2 c_\alpha}{\varpi_\alpha^2}\,\delta c_{A\alpha}
-
\frac{2 c_\alpha^2}{\varpi_\alpha^3}\,\delta\varpi_{A\alpha}
\right),
}
\]
\[
\boxed{
\delta B_{A,2}
=
\sum_\alpha
\left(
\frac{2 c_\alpha}{\varpi_\alpha^4}\,\delta c_{A\alpha}
-
\frac{4 c_\alpha^2}{\varpi_\alpha^5}\,\delta\varpi_{A\alpha}
\right).
}
\]
Therefore the exact BdG combination entering the even grouped outlet defect is
\[
\boxed{
\mathcal B_A
=
\sum_\alpha
\left[
2 c_\alpha\left(\frac{1}{\varpi_\alpha^4}+\frac{1}{9\varpi_\alpha^2}\right)\delta c_{A\alpha}
-
2 c_\alpha^2\left(\frac{2}{\varpi_\alpha^5}+\frac{1}{9\varpi_\alpha^3}\right)\delta\varpi_{A\alpha}
\right].
}
\]

So the BdG support sector feeds \(\mathcal K_A\) through one very specific weighted combination of lane-dependent coupling anisotropy and BdG-frequency anisotropy.

At the same time, the same static BdG anisotropy enters the odd grouped-normalization defect only through
\[
P_0\,\delta B_{A,0}.
\]
So the support sector can source \(\gamma_W\) even when the outgoing transfer is perfectly isotropic, simply because the normalization quotient reweights any static conservative lane splitting.

---

## 3. Exact Maxwell/mixed conservative bundle contribution

For each port-active internal pair \((U_{A,r},W_{A,r})\), define the isotropic baseline invariants
\[
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
\qquad
S_r=\Omega_{U,r}^2+\Omega_{W,r}^2,
\]
\[
Q_r=g_{U,r}^2\Omega_{W,r}^2+2 g_{U,r}g_{W,r}R_r+g_{W,r}^2\Omega_{U,r}^2,
\qquad
G_r=g_{U,r}^2+g_{W,r}^2,
\]
\[
P_r=\Omega_{U,r}^2 g_{W,r}+R_r g_{U,r}.
\]
The portwise conservative moments are
\[
Z_{A,0}^{(r)}=\frac{Q_{A,r}}{\Delta_{A,r}},
\qquad
Z_{A,2}^{(r)}=\frac{Q_{A,r}S_{A,r}-G_{A,r}\Delta_{A,r}}{\Delta_{A,r}^2}.
\]

Linearizing around the isotropic branch gives the exact variations
\[
\boxed{
\delta Z_{A,0}^{(r)}
=
\frac{\Delta_r\,\delta Q_{A,r}-Q_r\,\delta\Delta_{A,r}}{\Delta_r^2},
}
\]
\[
\boxed{
\delta Z_{A,2}^{(r)}
=
\frac{S_r}{\Delta_r^2}\,\delta Q_{A,r}
+
\frac{Q_r}{\Delta_r^2}\,\delta S_{A,r}
-
\frac{1}{\Delta_r}\,\delta G_{A,r}
+
\left(\frac{G_r}{\Delta_r^2}-\frac{2Q_rS_r}{\Delta_r^3}\right)\delta\Delta_{A,r}.
}
\]

The exact port contribution to the even outlet obstruction is therefore
\[
\mathcal Z_A^{(r)}:=\delta Z_{A,2}^{(r)}+\frac19\,\delta Z_{A,0}^{(r)},
\]
so
\[
\boxed{
\mathcal Z_A^{(r)}
=
\left(\frac{S_r}{\Delta_r^2}+\frac{1}{9\Delta_r}\right)\delta Q_{A,r}
+
\frac{Q_r}{\Delta_r^2}\,\delta S_{A,r}
-
\frac{1}{\Delta_r}\,\delta G_{A,r}
+
\left(\frac{G_r}{\Delta_r^2}-\frac{Q_r}{9\Delta_r^2}-\frac{2Q_rS_r}{\Delta_r^3}\right)\delta\Delta_{A,r}.
}
\]
Summing over ports gives
\[
\boxed{\mathcal Z_A=\sum_r \mathcal Z_A^{(r)}.}
\]

So the conservative Maxwell/mixed sector feeds \(\mathcal K_A\) through only four microscopic port variations:
\[
\delta Q_{A,r},\qquad \delta S_{A,r},\qquad \delta G_{A,r},\qquad \delta\Delta_{A,r}.
\]

---

## 4. Exact outgoing-transfer bundle contribution

The static outgoing-transfer moment is
\[
N_{A,0}^{(r)}=\frac{P_{A,r}^2}{\Delta_{A,r}^2}.
\]
So its first variation is
\[
\boxed{
\delta N_{A,0}^{(r)}
=
\frac{2P_r}{\Delta_r^2}\,\delta P_{A,r}
-
\frac{2P_r^2}{\Delta_r^3}\,\delta\Delta_{A,r}.
}
\]

The Maxwell/mixed contribution to the odd grouped-normalization obstruction is not \(\delta N_{A,0}^{(r)}\) alone, but
\[
\mathcal N_A^{(r)}:=\delta N_{A,0}^{(r)}+P_0\,\delta Z_{A,0}^{(r)}.
\]
Using the formula above for \(\delta Z_{A,0}^{(r)}\), one gets
\[
\boxed{
\mathcal N_A^{(r)}
=
\frac{P_0}{\Delta_r}\,\delta Q_{A,r}
+
\frac{2P_r}{\Delta_r^2}\,\delta P_{A,r}
-
\left(
\frac{P_0Q_r}{\Delta_r^2}
+
\frac{2P_r^2}{\Delta_r^3}
\right)\delta\Delta_{A,r}.
}
\]
Summing over ports gives
\[
\boxed{\mathcal N_A=\sum_r \mathcal N_A^{(r)}.}
\]

This formula is one of the key outputs of the stage.
It shows that the direct transfer sector feeds \(\mathcal G_A\) through only **three** microscopic port variations:
\[
\delta Q_{A,r},\qquad \delta P_{A,r},\qquad \delta\Delta_{A,r}.
\]
In particular, \(\delta S_{A,r}\) and \(\delta G_{A,r}\) drop out of the direct odd transfer bundle at linear order.

---

## 5. Primary microscopic variations for each Maxwell/mixed port

It is useful to unpack the port variations into the primitive grouped-lane data
\[
\delta\Omega_{U,A,r}^2,
\qquad
\delta\Omega_{W,A,r}^2,
\qquad
\delta R_{A,r},
\qquad
\delta g_{U,A,r},
\qquad
\delta g_{W,A,r}.
\]
Then the exact first variations are
\[
\boxed{
\delta\Delta_{A,r}
=
\Omega_{W,r}^2\,\delta\Omega_{U,A,r}^2
+
\Omega_{U,r}^2\,\delta\Omega_{W,A,r}^2
-
2R_r\,\delta R_{A,r},
}
\]
\[
\boxed{
\delta S_{A,r}
=
\delta\Omega_{U,A,r}^2+\delta\Omega_{W,A,r}^2,
}
\]
\[
\boxed{
\delta G_{A,r}
=
2 g_{U,r}\,\delta g_{U,A,r}+2 g_{W,r}\,\delta g_{W,A,r},
}
\]
\[
\boxed{
\delta P_{A,r}
=
 g_{W,r}\,\delta\Omega_{U,A,r}^2
+
 g_{U,r}\,\delta R_{A,r}
+
 R_r\,\delta g_{U,A,r}
+
 \Omega_{U,r}^2\,\delta g_{W,A,r},
}
\]
\[
\boxed{
\delta Q_{A,r}
=
 g_{W,r}^2\,\delta\Omega_{U,A,r}^2
+
 g_{U,r}^2\,\delta\Omega_{W,A,r}^2
+
 2 g_{U,r}g_{W,r}\,\delta R_{A,r}
}
\]
\[
\boxed{
\qquad
+
2\bigl(g_{U,r}\Omega_{W,r}^2+g_{W,r}R_r\bigr)\delta g_{U,A,r}
+
2\bigl(g_{W,r}\Omega_{U,r}^2+g_{U,r}R_r\bigr)\delta g_{W,A,r}.
}
\]

So once the physical moving-throat branch tells us how the grouped lanes perturb the microscopic port frequencies, port mixing, and wall–gauge couplings, the outlet obstructions are completely explicit.

---

## 6. Weak axisymmetric `Y_{20}` branch

On the weak axisymmetric branch write, for every microscopic grouped quantity,
\[
\delta X_A=\epsilon\,\lambda_A X^{(1)},
\qquad
(\lambda_{20},\lambda_{21},\lambda_{22})=\left(1,\frac12,-1\right).
\]
Then the two Stage 238 outlet obstructions collapse to two scalar microscopic amplitudes:
\[
\boxed{
\mathfrak K_1
:=
\frac19 K^{(1)}-M^{(1)}
-
\left(B_2^{(1)}+\frac19 B_0^{(1)}\right)
-
\left(Z_2^{(1)}+\frac19 Z_0^{(1)}\right),
}
\]
\[
\boxed{
\mathfrak G_1
:=
N_0^{(1)}-P_0 K^{(1)}+P_0 B_0^{(1)}+P_0 Z_0^{(1)}.
}
\]
So
\[
\mathcal K_A=\epsilon\,\lambda_A\,\mathfrak K_1,
\qquad
\mathcal G_A=\epsilon\,\lambda_A\,\mathfrak G_1.
\]
Therefore the direct outlet amplitudes are
\[
\boxed{
\kappa_1
=
\frac{3(1-\sigma_*)}{\sigma_* D_0}\,\mathfrak K_1,
}
\]
\[
\boxed{
\gamma_1
=
-\frac{1-\sigma_*}{9\sigma_* N_0}\,\mathfrak G_1.
}
\]

So the weak grouped-lane outlet problem has now collapsed one stage further.
It is no longer the whole bundle of grouped microscopic anisotropies.
It is just the pair
\[
(\mathfrak K_1,\mathfrak G_1).
\]

---

## 7. What Stage 239 changes

Stage 238 already reduced the linear grouped-anisotropy outlet problem to
\[
(\mathcal K_A,\mathcal G_A).
\]
But at that point those quantities were still phrased in terms of the grouped operator data
\[
(\delta D_{A,0},\delta D_{A,2},\delta D_{A,4},\delta N_{A,0}).
\]

Stage 239 now makes the microscopic meaning of those obstructions completely explicit.

The new theorem status is:

1. the hidden even outlet deformation is sourced only by
   \[
   \mathcal W_A,\quad \mathcal B_A,\quad \mathcal Z_A;
   \]
2. the hidden odd normalization deformation is sourced only by
   \[
   -P_0\delta K_A,
   \quad
   P_0\delta B_{A,0},
   \quad
   \mathcal N_A;
   \]
3. the direct transfer bundle \(\mathcal N_A\) depends only on
   \[
   \delta Q_{A,r},\quad \delta P_{A,r},\quad \delta\Delta_{A,r};
   \]
4. so the full weak grouped-anisotropy outlet bottleneck has collapsed to the scalar pair
   \[
   (\mathfrak K_1,\mathfrak G_1).
   \]

This is a substantial sharpening.
The next honest theorem gate is no longer

> “derive the grouped outlet anisotropies somehow.”

That part is now finished.

The next gate is:

> compute the actual moving-throat grouped-lane perturbations of
> \(
> K_A,\ M_A,\ c_{A\alpha},\ \varpi_{A\alpha},\ \Omega_{U/W,A,r},\ R_{A,r},\ g_{U/W,A,r}
> \)
> on the physical branch, and evaluate the two scalar amplitudes
> \(
> \mathfrak K_1,\ \mathfrak G_1
> \).

That is the direct continuation point.
