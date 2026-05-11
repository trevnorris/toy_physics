# Moving-Throat PDE — Stage 174: Static Self-Similarity and Exact Collapse of the Remaining Linear Grouped Loading Defect

## Purpose

Stage 241 reduced the remaining linear grouped `2.5`PN defect to one scalar,
\[
\Xi_{\rm load}
:=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
\]
That was already a strong collapse, but it still left the two static slopes
\(D_{01}/D_0\) and \(N_{01}/N_0\) looking like unrelated quantities.

The next honest step is therefore:

> rewrite the remaining defect directly as a **weighted failure of static self-similarity** between the wall baseline, the BdG support bundle, the conservative Maxwell/mixed bundle, and the outgoing-transfer bundle.

This stage does exactly that.

The main result is that the entire remaining linear grouped defect can be written as
\[
\boxed{
\Xi_{\rm load}
=
(\delta_N-\delta_K)
+
\omega_B(\delta_B-\delta_K)
+
\omega_Z(\delta_Z-\delta_K),
}
\]
where
\[
\delta_K:=\frac{K_1}{K},
\qquad
\delta_B:=\frac{B_0^{(1)}}{B_0},
\qquad
\delta_Z:=\frac{Z_0^{(1)}}{Z_0},
\qquad
\delta_N:=\frac{N_1}{N_0},
\]
and
\[
\omega_B:=\frac{B_0}{D_0},
\qquad
\omega_Z:=\frac{Z_0}{D_0},
\qquad
D_0=K-B_0-Z_0.
\]

So the last linear grouped normalization defect is not a generic mismatch of many coefficients.
It is exactly the weighted failure of the three support/transfer sectors to co-load with the static wall baseline.

That makes the next theorem gate much sharper than Stage 241:

> compute the **self-similarity defect fields** of the actual moving-throat branch.

If those vanish, then the linear grouped `2.5`PN defect vanishes automatically.

---

## 1. Exact static-slope decomposition of the conservative operator

From the full grouped bundle,
\[
D_0=K-B_0-Z_0,
\qquad
D_{01}=K_1-B_0^{(1)}-Z_0^{(1)}.
\]
Define the logarithmic static slopes
\[
\delta_D:=\frac{D_{01}}{D_0},
\qquad
\delta_K:=\frac{K_1}{K},
\qquad
\delta_B:=\frac{B_0^{(1)}}{B_0},
\qquad
\delta_Z:=\frac{Z_0^{(1)}}{Z_0},
\qquad
\delta_N:=\frac{N_1}{N_0},
\]
and the conservative static weights
\[
\omega_K:=\frac{K}{D_0},
\qquad
\omega_B:=\frac{B_0}{D_0},
\qquad
\omega_Z:=\frac{Z_0}{D_0}.
\]
Then
\[
\boxed{
\delta_D
=
\omega_K\,\delta_K
-
\omega_B\,\delta_B
-
\omega_Z\,\delta_Z.
}
\]
Because
\[
D_0=K-B_0-Z_0,
\]
one also has the exact weight identity
\[
\boxed{
\omega_K-\omega_B-\omega_Z=1.
}
\]
Therefore the Stage 241 defect becomes
\[
\Xi_{\rm load}
=
\delta_N-
\omega_K\delta_K
+
\omega_B\delta_B
+
\omega_Z\delta_Z,
\]
which can be rewritten, using the weight identity above, as
\[
\boxed{
\Xi_{\rm load}
=
(\delta_N-\delta_K)
+
\omega_B(\delta_B-\delta_K)
+
\omega_Z(\delta_Z-\delta_K).
}
\]

This is the first exact statement that the wall baseline slope \(\delta_K\) is the natural reference for the remaining grouped loading defect.

---

## 2. Microscopic logarithmic forms of the three support/transfer slopes

The point of the previous formula is that the three sector slopes \(\delta_B,\delta_Z,\delta_N\) are themselves exact weighted averages of microscopic logarithmic deformations.

### 2.1 BdG support bundle

For each BdG mode \(\alpha\),
\[
B_{0,\alpha}=\frac{c_\alpha^2}{\varpi_\alpha^2},
\qquad
B_0=\sum_\alpha B_{0,\alpha}.
\]
Under a weak-axisymmetric grouped perturbation,
\[
c_{A\alpha}=c_\alpha+\epsilon\lambda_A c_\alpha^{(1)},
\qquad
\varpi_{A\alpha}=\varpi_\alpha+\epsilon\lambda_A\varpi_\alpha^{(1)},
\]
one gets
\[
\boxed{
\delta_B
=
\sum_\alpha \rho_\alpha^{(B)}
\left(
2\frac{c_\alpha^{(1)}}{c_\alpha}
-
2\frac{\varpi_\alpha^{(1)}}{\varpi_\alpha}
\right),
}
\]
with normalized mode weights
\[
\rho_\alpha^{(B)}:=\frac{B_{0,\alpha}}{B_0},
\qquad
\sum_\alpha \rho_\alpha^{(B)}=1.
\]
So equivalently,
\[
\boxed{
\delta_B
=
\sum_\alpha \rho_\alpha^{(B)}
\,2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right).
}
\]

### 2.2 Conservative Maxwell/mixed static bundle

For each port-active pair \(r\),
\[
Z_0^{(r)}=\frac{Q_r}{\Delta_r},
\qquad
Z_0=\sum_r Z_0^{(r)}.
\]
So the weak-axisymmetric static slope is
\[
\boxed{
\delta_Z
=
\sum_r \rho_r^{(Z)}
\left(
\frac{Q_{1r}}{Q_r}
-
\frac{\Delta_{1r}}{\Delta_r}
\right),
}
\]
with normalized port weights
\[
\rho_r^{(Z)}:=\frac{Z_0^{(r)}}{Z_0},
\qquad
\sum_r \rho_r^{(Z)}=1.
\]
Equivalently,
\[
\boxed{
\delta_Z
=
\sum_r \rho_r^{(Z)}\,\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right).
}
\]

### 2.3 Outgoing-transfer bundle

For each outgoing port,
\[
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2},
\qquad
N_0=\sum_r N_0^{(r)}.
\]
So the weak-axisymmetric outgoing-transfer slope is
\[
\boxed{
\delta_N
=
\sum_r \rho_r^{(N)}
\left(
2\frac{P_{1r}}{P_r}
-
2\frac{\Delta_{1r}}{\Delta_r}
\right),
}
\]
with normalized outgoing weights
\[
\rho_r^{(N)}:=\frac{N_0^{(r)}}{N_0},
\qquad
\sum_r \rho_r^{(N)}=1.
\]
Equivalently,
\[
\boxed{
\delta_N
=
\sum_r \rho_r^{(N)}\,2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right).
}
\]

So the last grouped loading defect is already controlled by three explicit weighted logarithmic deformations of the microscopic moving-throat bundle.

---

## 3. Exact wall-referenced self-similarity defect fields

Because \(\delta_K\) is the natural wall-baseline slope, define the wall-referenced microscopic self-similarity defect fields
\[
\Sigma_\alpha^{(B)}
:=
2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right)-\delta_K,
\]
\[
\Sigma_r^{(Z)}
:=
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)-\delta_K,
\]
\[
\Sigma_r^{(N)}
:=
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)-\delta_K.
\]
Then
\[
\delta_B-\delta_K
=
\sum_\alpha \rho_\alpha^{(B)}\Sigma_\alpha^{(B)},
\]
\[
\delta_Z-\delta_K
=
\sum_r \rho_r^{(Z)}\Sigma_r^{(Z)},
\]
\[
\delta_N-\delta_K
=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}.
\]
So the remaining grouped loading defect becomes
\[
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}
+
\omega_B \sum_\alpha \rho_\alpha^{(B)}\Sigma_\alpha^{(B)}
+
\omega_Z \sum_r \rho_r^{(Z)}\Sigma_r^{(Z)}.
}
\]

This is the sharpest exact formula reached so far.

It says that the remaining linear grouped `2.5`PN defect is not sourced by the full microscopic bundle independently.
It is sourced only by the weighted failure of:

1. the BdG support bundle,
2. the conservative Maxwell/mixed bundle,
3. and the outgoing-transfer bundle,

to co-load with the static wall baseline.

---

## 4. Static self-similarity theorem

A particularly useful sufficient condition is now immediate.

Assume the weak-axisymmetric moving-throat branch is **statically self-similar** in the sense that, relative to the wall baseline slope \(\delta_K\), all three defect families vanish:
\[
\Sigma_\alpha^{(B)}=0
\quad \text{for every BdG mode } \alpha,
\]
\[
\Sigma_r^{(Z)}=0
\quad \text{for every conservative Maxwell/mixed port } r,
\]
\[
\Sigma_r^{(N)}=0
\quad \text{for every outgoing-transfer port } r.
\]
Equivalently,
\[
2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right)=\delta_K,
\]
\[
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)=\delta_K,
\]
\[
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)=\delta_K.
\]
Then the exact formula above gives
\[
\boxed{
\Xi_{\rm load}=0.
}
\]

So the remaining linear grouped normalization defect vanishes automatically under static self-similarity.

This is the grouped weak-axisymmetric analog of the Stage 252 tangent-compensation theorem: isotropic bundle drift did not push the system off the compensated Family-1 branch, and now weak grouped loading does not generate a linear grouped `2.5`PN defect if the full static support/transfer bundle co-loads self-similarly with the wall baseline.

---

## 5. Combined consequence on the even-preserving branch

Stage 241 already proved that on the even-preserving branch
\[
u_2^{(1)}=0
\quad\Longrightarrow\quad
D_{21}=-\frac{D_{01}}{9},
\qquad
D_{41}=-\frac{D_{01}}{27},
\]
and that the grouped defect pattern is
\[
\Delta_Q^{(20)}=\epsilon\,\Xi_{\rm load},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_{\rm load},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_{\rm load}.
\]
Therefore, on the even-preserving branch,
\[
\boxed{
\text{static self-similarity}
\quad\Longrightarrow\quad
\Delta_Q^{(20)}=\Delta_Q^{(21)}=\Delta_Q^{(22)}=0
\quad \text{at linear order.}
}
\]

So the remaining linear grouped `2.5`PN problem has narrowed again.
It is no longer to compute \(D_{01}/D_0\) and \(N_{01}/N_0\) separately.
It is to determine whether the actual moving-throat branch breaks static self-similarity, and if so through which defect family:

- BdG support,
- conservative Maxwell/mixed loading,
- or outgoing transfer.

---

## 6. What Stage 242 changes

Before this stage, the live theorem gate was:

> compute the static conservative slope \(D_{01}/D_0\) and the static outgoing-transfer slope \(N_{01}/N_0\).

After this stage, that is no longer the sharpest formulation.

The new theorem status is:

1. the remaining linear grouped defect is exactly
   \[
   \Xi_{\rm load}
   =
   (\delta_N-\delta_K)
   +
   \omega_B(\delta_B-\delta_K)
   +
   \omega_Z(\delta_Z-\delta_K);
   \]
2. each of the three support/transfer slopes is an exact weighted average of microscopic logarithmic deformations;
3. so the last linear grouped `2.5`PN defect is the weighted failure of static self-similarity relative to the wall baseline;
4. and if the moving-throat branch is statically self-similar on the even-preserving branch, the full linear grouped `2.5`PN defect vanishes.

That is the clean continuation point.

The next honest theorem gate is now:

> compute the microscopic wall-referenced defect fields
> \(\Sigma_\alpha^{(B)},\Sigma_r^{(Z)},\Sigma_r^{(N)}\)
> on the actual weak-axisymmetric moving-throat branch.

That will tell us whether the first linear grouped `2.5`PN defect comes from support non-self-similarity, conservative Maxwell/mixed non-self-similarity, or outgoing-transfer non-self-similarity.
