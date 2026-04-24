# Stage V2-12 — STF Angular Source-Map Theorem

## 0. Purpose

This stage audits whether the surviving orbital/worldtube STF quadrupole branch has an **angular normalization defect** when mapped into the moving-throat grouped real \(P_2\) ports.

The result is:

\[
\boxed{\widehat m_{\rm ang}=1}
\]

in the canonical normalized real-STF basis. Therefore the remaining quadrupole-normalization problem is not angular. It is radial/axial, source-amplitude, and outgoing-port normalization data.

This stage is deliberately narrow. It does not compute the radial overlap integrals, the Maxwell/mixed transfer prefactor, or the final product \( \widehat m_0^2P_0 \). It proves that the angular map itself is an identity.

---

## 1. Canonical real STF basis

Let \(n^i\) be the unit direction on the throat/mouth sphere. Use five real symmetric trace-free tensors \(E_A^{ij}\), with

\[
A\in\{20,21c,21s,22c,22s\}.
\]

The script uses

\[
E_{20}
=
\frac{1}{\sqrt6}
\begin{pmatrix}
-1&0&0\\
0&-1&0\\
0&0&2
\end{pmatrix},
\]

\[
E_{21c}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&0&1\\
0&0&0\\
1&0&0
\end{pmatrix},
\qquad
E_{21s}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&0&0\\
0&0&1\\
0&1&0
\end{pmatrix},
\]

\[
E_{22c}
=
\frac1{\sqrt2}
\begin{pmatrix}
1&0&0\\
0&-1&0\\
0&0&0
\end{pmatrix},
\qquad
E_{22s}
=
\frac1{\sqrt2}
\begin{pmatrix}
0&1&0\\
1&0&0\\
0&0&0
\end{pmatrix}.
\]

The audit verifies

\[
{\rm Tr}(E_A)=0,
\]

and

\[
{\rm Tr}(E_AE_B)=\delta_{AB}.
\]

So these five tensors are an orthonormal real STF basis.

---

## 2. Normalized real \(l=2\) angular functions

Define

\[
Y_A(n)
=
\sqrt{\frac{15}{8\pi}}\,
E_A^{ij}n_in_j.
\]

The fourth-moment identity on the unit sphere is

\[
\int_{S^2}n_in_jn_kn_l\,d\Omega
=
\frac{4\pi}{15}
\left(
\delta_{ij}\delta_{kl}
+\delta_{ik}\delta_{jl}
+\delta_{il}\delta_{jk}
\right).
\]

Because \(E_A\) is trace-free,

\[
\int_{S^2}
(E_A^{ij}n_in_j)(E_B^{kl}n_kn_l)\,d\Omega
=
\frac{8\pi}{15}
{\rm Tr}(E_AE_B).
\]

Therefore

\[
\int_{S^2}Y_A(n)Y_B(n)\,d\Omega
=
{\rm Tr}(E_AE_B)
=
\delta_{AB}.
\]

The script verifies this twice:

1. by the fourth-moment identity;
2. by direct \((\theta,\phi)\) integration of the explicit real harmonics.

Both checks return the \(5\times5\) identity matrix.

---

## 3. Source-map theorem

Let the angular part of the orbital/worldtube STF source be

\[
S(n)=\sum_B S_B Y_B(n).
\]

The port projection onto the same real harmonic basis is

\[
S_A^{\rm port}
=
\int_{S^2}Y_A(n)S(n)\,d\Omega.
\]

Using orthonormality,

\[
S_A^{\rm port}
=
\sum_B S_B
\int_{S^2}Y_A(n)Y_B(n)\,d\Omega
=
S_A.
\]

Thus the angular source-map matrix is exactly

\[
\boxed{
M_{AB}^{\rm ang}=\delta_{AB}.
}
\]

There is no angular rescaling, no angular kernel, and no missing angular component.

---

## 4. Orbital/worldtube STF quadrupole reconstruction

For a point/worldtube quadrupole

\[
Q_{ij}
=
\mu
\left(
x_ix_j-\frac{r^2}{3}\delta_{ij}
\right),
\qquad
r^2=x^2+y^2+z^2,
\]

the real STF coefficients are

\[
Q_A={\rm Tr}(E_AQ).
\]

The script derives

\[
Q_{20}
=
-\frac{\sqrt6\mu}{6}(x^2+y^2-2z^2),
\]

\[
Q_{21c}=\sqrt2\,\mu xz,
\qquad
Q_{21s}=\sqrt2\,\mu yz,
\]

\[
Q_{22c}=\frac{\sqrt2\mu}{2}(x^2-y^2),
\qquad
Q_{22s}=\sqrt2\,\mu xy.
\]

It then checks the reconstruction

\[
Q_{ij}=\sum_A Q_AE_{Aij}.
\]

The residual is exactly zero.

It also checks the angular identity

\[
\sqrt{\frac{15}{8\pi}}Q_{ij}n^in^j
=
\sum_A Q_A Y_A(n),
\]

with exact zero residual.

---

## 5. Norm and grouped metric

For the full five-mode real STF source,

\[
\int_{S^2}S(n)^2\,d\Omega
=
S_{20}^2+S_{21c}^2+S_{21s}^2+S_{22c}^2+S_{22s}^2.
\]

When the two \(21\) real components are grouped together and the two \(22\) real components are grouped together, a representative grouped vector

\[
x_{\rm grp}=(x_{20},x_{21},x_{22})
\]

corresponds to the full vector

\[
(x_{20},x_{21},x_{21},x_{22},x_{22}).
\]

Therefore the full five-mode norm becomes

\[
x_{20}^2+2x_{21}^2+2x_{22}^2
=
x_{\rm grp}^T
\operatorname{diag}(1,2,2)
x_{\rm grp}.
\]

This is the grouped metric

\[
\boxed{
G_{\rm grp}=\operatorname{diag}(1,2,2).
}
\]

So the grouped metric used in the \(20/21/22\) projector calculus is exactly the norm inherited from the five real STF modes.

---

## 6. Convention bridge to the older \(\Pi\) notation

The 2.5PN package uses an older real-STF component convention in which the normalized port variables are related by

\[
q_{20}=\sqrt{\frac23}\Pi_{20},
\]

\[
q_{21c}=\Pi_{21c},
\qquad
q_{21s}=\Pi_{21s},
\]

\[
q_{22c}=2\Pi_{22c},
\qquad
q_{22s}=2\Pi_{22s}.
\]

The conversion matrix is

\[
B_{\Pi\to q}
=
\operatorname{diag}
\left(
\sqrt{\frac23},1,1,2,2
\right).
\]

The script verifies

\[
\det B_{\Pi\to q}
=
\frac{4\sqrt6}{3}
\neq0,
\]

and

\[
{\rm rank}(B_{\Pi\to q})=5.
\]

So even in the older \(\Pi\) convention, the map is invertible. It introduces only a convention conversion, not an angular obstruction.

---

## 7. Interpretation

The audit proves the angular part of the quadrupole map is complete and exactly normalized in the canonical real-STF basis:

\[
\boxed{
\widehat m_{\rm ang}=1.
}
\]

Therefore the remaining source normalization can be factored as

\[
\widehat m_0
=
\widehat m_{\rm ang}\widehat m_{\rm rad}
=
\widehat m_{\rm rad}.
\]

The unresolved product is still

\[
\widehat m_0^2 P_0
=
\frac{54Gc_s^5}{5a^5c^5},
\]

but V2-12 shows that any miss in this equation must come from

1. radial/axial overlap amplitudes,
2. mouth/worldtube source amplitude normalization,
3. the Maxwell/mixed transfer factor \(N_0\),
4. the conservative denominator \(D_0\),
5. or the final moving-throat branch selection.

It cannot be blamed on a missing angular STF component or an angular normalization mismatch.

---

## 8. SymPy result

The script reports:

```text
all_STF_traces_zero: PASS
STF_basis_orthonormal: PASS
angular_gram_from_moments_identity: PASS
direct_angular_gram_identity: PASS
source_projection_identity: PASS
orbital_STF_reconstructs_exactly: PASS
angular_function_reconstructs_exactly: PASS
source_norm_identity: PASS
grouped_metric_matches_full_pair_norm: PASS
Pi_to_q_map_invertible: PASS
```

Final verdict:

\[
\boxed{
\text{The STF angular source-map has no angular normalization defect.}
}
\]

\[
\boxed{
\text{The remaining quadrupole normalization gap is radial/axial and dynamical.}
}
\]
