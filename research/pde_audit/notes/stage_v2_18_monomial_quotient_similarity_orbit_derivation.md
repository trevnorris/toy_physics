# Stage V2-18 — Monomial Quotient and Similarity-Orbit Audit

## 0. Purpose

This stage turns the weak-axisymmetric grouped-`P2` defect into a quotient problem.

Earlier stages showed that the grouped weak-axisymmetric branch has a single leading
outgoing-prefactor slope

\[
\Xi_1=\frac{P_1}{P_0}
=\frac{N_{01}}{N_0}-\frac{D_{01}}{D_0}.
\]

The purpose of this stage is to prove that this scalar is not an arbitrary
microscopic drift. It is one component of an exact quotient-residual packet

\[
q=(q_{\rm tr},q_{\rm nt},q_\eta)
=
\left(
\delta\ln\mathfrak C_{{\rm tr},*},
\delta\ln\mathfrak C_{{\rm nt},*},
\delta\ln\epsilon_\eta
\right),
\]

where the remaining microscopic drift directions are pure similarity-orbit
motion.

---

## 1. Original microscopic monomial map

Use the microscopic grouped drift vector

\[
\delta\mathbf x
=
(
\delta\ln\lambda_W,\,
\delta\ln c_{\eta U},\,
\delta\ln\gamma,\,
\delta\ln K_U,\,
\delta\ln K_\eta^{\rm eff},\,
\delta\ln K_W^{\rm eff},\,
\delta\ln\mu_W,\,
\delta\ln T_U
)^T.
\]

The exact monomials are

\[
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\left(
\frac{\lambda_W^2\mu_W}
{K_\eta^{\rm eff}(K_W^{\rm eff})^2}
\right)
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_UK_W^{\rm eff}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
\]

\[
\epsilon_\eta
=
\frac{c_{\eta U}^2}{K_UK_\eta^{\rm eff}}.
\]

Taking logarithmic differentials gives

\[
M_*\delta\mathbf x
=
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}.
\]

In the above variable order,

\[
M_*=
\begin{pmatrix}
0&1+\delta&1+\delta&-(2+\delta+\chi)&0&0&0&1+\chi\\
2+2E&0&2E&F-E&-1&-(2+E)&1&-F\\
0&2&0&-1&-1&0&0&0
\end{pmatrix},
\]

where \(\delta=\delta_{U,*}\), \(\chi=\chi_{0,*}\), \(E=E_*\), and \(F=F_*\).

The script verifies the quoted rank witness:

\[
\det M_*^{(T_U,K_\eta^{\rm eff},\mu_W)}
=
1+\chi.
\]

Therefore, on the physical branch \(1+\chi\ne0\),

\[
{\rm rank}(M_*)=3,
\qquad
\dim\ker M_*=5.
\]

Thus the microscopic drift space splits into a five-dimensional similarity
orbit plus three quotient directions.

---

## 2. Exact normal basis and similarity-orbit split

A convenient exact right inverse is

\[
n_{\rm tr}
=
\left(
0,0,0,0,0,0,
\frac{F_*}{1+\chi_{0,*}},
\frac{1}{1+\chi_{0,*}}
\right)^T,
\]

\[
n_{\rm nt}
=
(0,0,0,0,0,0,1,0)^T,
\]

\[
n_\eta
=
(0,0,0,0,-1,0,-1,0)^T.
\]

Let

\[
N=(n_{\rm tr},n_{\rm nt},n_\eta).
\]

The script verifies

\[
M_*N=I_3.
\]

Therefore any microscopic tangent has the exact split

\[
\delta\mathbf x
=
\delta\mathbf x_{\rm orbit}
+
Nq,
\]

with

\[
M_*\delta\mathbf x_{\rm orbit}=0.
\]

Equivalently,

\[
P_{\rm normal}=NM_*,
\qquad
P_{\rm orbit}=I-NM_*,
\]

are idempotent projectors, and

\[
M_*P_{\rm orbit}=0.
\]

This is the precise algebraic meaning of the similarity orbit: the orbit
directions preserve the three monomials, and only \(q\) moves the physical
quotient.

---

## 3. Normalized Stage-12/13 version

Using normalized variables

\[
(
d\ln G_W,\,
d\ln G_U,\,
d\ln R,\,
d\ln K,\,
d\ln M,\,
d\ln\Omega_U,\,
d\ln\Omega_W,\,
d\ln\delta_U
),
\]

the monomials become

\[
\mathfrak C_{{\rm tr},*}
=
\left(
\frac{RG_U}{\Omega_U^2G_W}
\right)^{1+\delta_U}
\delta_U^{1+\chi_0},
\]

\[
\mathfrak C_{{\rm nt},*}
=
\frac{MG_W^2}{K\Omega_W^4}
\left(
\frac{R^2\sigma}{\Omega_U^2\Omega_W^2}
\right)^{E_*}
\delta_U^{-F_*},
\]

\[
\epsilon_\eta
=
\frac{MG_U^2}{K\Omega_U^2}.
\]

The corresponding drift matrix is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta)&1+\delta&1+\delta&0&0&-2(1+\delta)&0&1+\chi\\
2&0&2E&-1&1&-2E&-(4+2E)&-F\\
0&2&0&-1&1&-2&0&0
\end{pmatrix}.
\]

The script verifies a nonzero rank witness

\[
\det M_{\rm norm}^{(\delta_U,M,\Omega_W)}
=
2(E+2)(1+\chi).
\]

So the normalized quotient map also has rank three on the physical branch.

The zero-defect equations solve triangularly as

\[
d\ln\delta_U
=
-\frac{1+\delta_U}{1+\chi_0}
\left(
d\ln R+d\ln G_U-d\ln G_W-2d\ln\Omega_U
\right),
\]

\[
d\ln M
=
d\ln K-2d\ln G_U+2d\ln\Omega_U,
\]

\[
d\ln\Omega_W
=
\frac{
d\ln G_W-d\ln G_U+(1-E_*)d\ln\Omega_U
+E_*d\ln R-\frac{F_*}{2}d\ln\delta_U
}
{E_*+2}.
\]

The script substitutes these formulas back into \(M_{\rm norm}v\) and verifies
that all three monomial drifts vanish.

---

## 4. Support-blind extension

If explicit BdG-support primitive drift directions

\[
d\ln\lambda_B,\qquad d\ln\varpi
\]

are appended, the monomial matrix gains two exact zero columns.

The script verifies that these two columns are zero and that the extended
monomial kernel dimension becomes

\[
2+5=7.
\]

This is an important boundary on what the quotient theorem proves. It governs
the weak-axisymmetric normalization defect \(\Xi_1\), but it does not by itself
close every conservative even gate involving BdG moments.

---

## 5. Exact physical defect compiler

The quotient residuals map into the physical first-order defect triplet by

\[
\Theta_1
=
-C_{\rm tr}q_{\rm tr},
\]

\[
\Xi_1
=
A_{\rm tr}q_{\rm tr}+q_{\rm nt},
\]

\[
\mathcal R_1
=
-\Xi_1
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
\]

In matrix form,

\[
\begin{pmatrix}
\Theta_1\\
\Xi_1\\
\mathcal R_1
\end{pmatrix}
=
\begin{pmatrix}
-C_{\rm tr}&0&0\\
A_{\rm tr}&1&0\\
-A_{\rm tr}&-1&-\dfrac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}
\end{pmatrix}
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}.
\]

The script verifies

\[
\det D_{\rm defect}
=
\frac{C_{\rm tr}\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}},
\]

so the compiler is invertible when

\[
C_{\rm tr}\ne0,\qquad
\epsilon_{\eta,*}\ne0,\qquad
\epsilon_{\eta,*}\ne1.
\]

The inverse is

\[
q_{\rm tr}
=
-\frac{\Theta_1}{C_{\rm tr}},
\]

\[
q_{\rm nt}
=
\Xi_1+\frac{A_{\rm tr}}{C_{\rm tr}}\Theta_1,
\]

\[
q_\eta
=
-\frac{1-\epsilon_{\eta,*}}{\epsilon_{\eta,*}}
(\mathcal R_1+\Xi_1).
\]

Thus

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\quad\Longleftrightarrow\quad
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

And since \(q=M_*\delta\mathbf x\),

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\quad\Longleftrightarrow\quad
\delta\mathbf x\in\ker M_*.
\]

That is the exact similarity-orbit zero-defect theorem.

---

## 6. Bridge back to the prefactor slope

The outgoing-prefactor slope is

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_1
=
\frac{N_{01}D_0-N_0D_{01}}{D_0^2}.
\]

Therefore

\[
\frac{P_1}{P_0}
=
\frac{N_{01}}{N_0}
-
\frac{D_{01}}{D_0}.
\]

The script verifies this identity exactly.

Combining this with the defect compiler gives

\[
\boxed{
\Xi_1
=
\frac{P_1}{P_0}
=
A_{\rm tr}q_{\rm tr}+q_{\rm nt}.
}
\]

So the weak-axisymmetric prefactor obstruction is not arbitrary. It is the
projection of the branch tangent onto the two quotient coordinates
\((q_{\rm tr},q_{\rm nt})\), while \(q_\eta\) enters the selected-branch residual
\(\mathcal R_1\).

---

## 7. Stage verdict

The audit passes.

The quotient variables

\[
(q_{\rm tr},q_{\rm nt},q_\eta)
\]

are exact physical coordinates on the defect quotient. The five remaining
microscopic directions are similarity-orbit directions. The physical
weak-axisymmetric first-order defects are exactly invertible functions of the
quotient residuals.

The immediate carry-forward theorem gate is therefore:

\[
\boxed{
\text{Extract the actual moving-throat branch tangent and test whether }
M_*\delta\mathbf x=0.
}
\]

Equivalently:

\[
\boxed{
\text{Compute }
\Xi_1=\frac{P_1}{P_0}
\text{ and the companion residuals from the quotient packet.}
}
\]

If \(M_*\delta\mathbf x=0\), the weak-axisymmetric prefactor problem is pure
similarity-gauge at first order. If not, the branch leaves the monomial-preserving
orbit, and the exact residual coordinates identify which physical defect
\((\Theta_1,\Xi_1,\mathcal R_1)\) failed.
