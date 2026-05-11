# Moving-Throat PDE — Stage 175: Wall-Normalized Port Shape Variables and the Outgoing-Load Theorem

## Purpose

Stage 242 reduced the remaining linear grouped `2.5`PN defect to the wall-referenced microscopic self-similarity fields
\[
\Sigma_\alpha^{(B)}
=
2\,\delta\ln\!\left(\frac{c_\alpha}{\varpi_\alpha}\right)-\delta_K,
\qquad
\Sigma_r^{(Z)}
=
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)-\delta_K,
\qquad
\Sigma_r^{(N)}
=
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)-\delta_K.
\]
That already showed that the last linear grouped loading defect is a weighted failure of static self-similarity relative to the wall baseline slope \(\delta_K\).

The next honest step is now smaller than “solve the whole branch.” It is:

> factor the microscopic support/port data by the wall baseline and determine whether the remaining grouped defect is a generic static-shape problem, or a much narrower **outgoing-load mismatch**.

This stage performs that factorization exactly.

The main result is that the three Stage 242 defect families collapse to **wall-normalized shape/load variables**:
\[
B_{0,\alpha}=K\,\chi_\alpha^2,
\qquad
Z_0^{(r)}=K\,\Upsilon_r,
\qquad
N_0^{(r)}=\Lambda_r^2,
\]
with
\[
\Sigma_\alpha^{(B)}=\delta\ln\chi_\alpha^2,
\qquad
\Sigma_r^{(Z)}=\delta\ln\Upsilon_r,
\qquad
\Sigma_r^{(N)}=\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right).
\]

So the remaining linear grouped defect is no longer a generic static mismatch of the whole microscopic bundle.
It is exactly:

1. wall-normalized BdG support-shape drift,
2. wall-normalized conservative port-shape drift,
3. and a **wall-loading law for the outgoing transfer factor**.

The sharpest new theorem is:

> on branches that preserve the wall-normalized BdG and conservative port shapes, the full remaining linear grouped `2.5`PN defect is carried only by the outgoing load factor \(\Lambda_r=P_r/\Delta_r\).

Even more sharply, if all wall-normalized shapes are frozen, then the defect does **not** vanish:
\[
\Xi_{\rm load}=-\delta_K.
\]
So a naive common self-similarity branch is not enough.
The outgoing transfer sector must actively load with the wall baseline.

---

## 1. Carry-forward microscopic static bundles

From Stage 242, the static grouped bundles were
\[
B_{0,\alpha}=\frac{c_\alpha^2}{\varpi_\alpha^2},
\qquad
Z_0^{(r)}=\frac{Q_r}{\Delta_r},
\qquad
N_0^{(r)}=\frac{P_r^2}{\Delta_r^2},
\]
with
\[
\Delta_r=\Omega_{U,r}^2\Omega_{W,r}^2-R_r^2,
\]
\[
Q_r=g_{U,r}^2\Omega_{W,r}^2+2g_{U,r}g_{W,r}R_r+g_{W,r}^2\Omega_{U,r}^2,
\]
\[
P_r=\Omega_{U,r}^2 g_{W,r}+R_r g_{U,r}.
\]

The wall baseline slope is
\[
\delta_K:=\frac{K_1}{K},
\]
and the Stage 242 grouped defect was
\[
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\Sigma_r^{(N)}
+
\omega_B \sum_\alpha \rho_\alpha^{(B)}\Sigma_\alpha^{(B)}
+
\omega_Z \sum_r \rho_r^{(Z)}\Sigma_r^{(Z)}.
\]

So everything now depends on how the three self-similarity fields are most naturally factorized by the wall baseline.

---

## 2. Exact wall-normalized factorization of the microscopic bundles

### 2.1 BdG support shape variable

Define the wall-normalized BdG support amplitude
\[
\boxed{
\chi_\alpha:=\frac{c_\alpha}{\sqrt K\,\varpi_\alpha}.
}
\]
Then
\[
\boxed{
B_{0,\alpha}=K\,\chi_\alpha^2.
}
\]

So the Stage 242 BdG self-similarity defect is simply
\[
\boxed{
\Sigma_\alpha^{(B)}
=
\delta\ln\chi_\alpha^2.
}
\]

Thus the BdG support sector contributes to the remaining grouped loading defect only through the drift of the wall-normalized support shape \(\chi_\alpha\).

### 2.2 Wall-normalized port variables

For each Maxwell/mixed port \(r\), define wall-normalized primitive variables
\[
\hat\Omega_{U,r}^2:=\frac{\Omega_{U,r}^2}{K},
\qquad
\hat\Omega_{W,r}^2:=\frac{\Omega_{W,r}^2}{K},
\qquad
\hat R_r:=\frac{R_r}{K},
\]
\[
\hat g_{U,r}:=\frac{g_{U,r}}{K},
\qquad
\hat g_{W,r}:=\frac{g_{W,r}}{K}.
\]
Then define the normalized port invariants
\[
\hat\Delta_r
=
\hat\Omega_{U,r}^2\hat\Omega_{W,r}^2-\hat R_r^2,
\]
\[
\hat Q_r
=
\hat g_{U,r}^2\hat\Omega_{W,r}^2
+
2\hat g_{U,r}\hat g_{W,r}\hat R_r
+
\hat g_{W,r}^2\hat\Omega_{U,r}^2,
\]
\[
\hat P_r
=
\hat\Omega_{U,r}^2 \hat g_{W,r}
+
\hat R_r \hat g_{U,r}.
\]

A direct substitution gives the exact homogeneity laws
\[
\boxed{
\Delta_r=K^2\hat\Delta_r,
\qquad
Q_r=K^3\hat Q_r,
\qquad
P_r=K^2\hat P_r.
}
\]

So the static conservative and outgoing-transfer bundles become
\[
\boxed{
Z_0^{(r)}=K\,\Upsilon_r,
\qquad
\Upsilon_r:=\frac{\hat Q_r}{\hat\Delta_r}=\frac{Q_r}{K\Delta_r},
}
\]
\[
\boxed{
N_0^{(r)}=\Lambda_r^2,
\qquad
\Lambda_r:=\frac{\hat P_r}{\hat\Delta_r}=\frac{P_r}{\Delta_r}.
}
\]

This is the central factorization of the stage.

- \(\Upsilon_r\) is the **wall-normalized conservative port shape**.
- \(\Lambda_r\) is the **dimensionless outgoing load factor**.

### 2.3 Exact rewrite of the Stage 242 port defects

Using the factorization above,
\[
\delta\ln\!\left(\frac{Q_r}{\Delta_r}\right)
=
\delta_K+\delta\ln\Upsilon_r,
\]
\[
2\,\delta\ln\!\left(\frac{P_r}{\Delta_r}\right)
=
2\,\delta\ln\Lambda_r.
\]
Therefore
\[
\boxed{
\Sigma_r^{(Z)}=\delta\ln\Upsilon_r,
}
\]
\[
\boxed{
\Sigma_r^{(N)}=2\,\delta\ln\Lambda_r-\delta_K
=\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right).
}
\]

So the conservative port self-similarity defect is pure **shape drift**, while the outgoing port defect is a **load mismatch** relative to the wall baseline.

That is the most important structural sharpening of the stage.

---

## 3. Exact wall-normalized load-shape formula for the remaining grouped defect

Define the weighted aggregate defect measures
\[
\Theta_B
:=
\sum_\alpha \rho_\alpha^{(B)}\,\delta\ln\chi_\alpha^2,
\]
\[
\Theta_Z
:=
\sum_r \rho_r^{(Z)}\,\delta\ln\Upsilon_r,
\]
\[
\Theta_N
:=
\sum_r \rho_r^{(N)}\,\delta\ln\!\left(\frac{\Lambda_r^2}{K}\right).
\]
Then the exact Stage 242 formula becomes
\[
\boxed{
\Xi_{\rm load}
=
\Theta_N+\omega_B\Theta_B+\omega_Z\Theta_Z.
}
\]

This is the cleanest exact form of the remaining linear grouped `2.5`PN defect reached so far.

It shows that the last grouped defect is not a generic mismatch of all microscopic variables.
It is exactly the weighted sum of:

1. wall-normalized support-shape drift,
2. wall-normalized conservative port-shape drift,
3. outgoing load mismatch.

So the remaining theorem gate has become smaller still.

---

## 4. Conservative-shape theorem

A particularly important special case is the branch on which both conservative shape families are preserved:
\[
\delta\ln\chi_\alpha^2=0
\qquad\text{for all BdG modes }\alpha,
\]
\[
\delta\ln\Upsilon_r=0
\qquad\text{for all conservative Maxwell/mixed ports }r.
\]
Then
\[
\Theta_B=0,
\qquad
\Theta_Z=0,
\]
and the remaining grouped defect collapses to
\[
\boxed{
\Xi_{\rm load}
=
\sum_r \rho_r^{(N)}\left(2\,\delta\ln\Lambda_r-\delta_K\right).
}
\]

So on conservative-shape-preserving branches, the full remaining linear grouped `2.5`PN defect is carried only by the outgoing load factor \(\Lambda_r\).

This is the stage’s main theorem.

---

## 5. Naive common-self-similarity no-go

The previous formula has an immediate and surprisingly sharp consequence.

Suppose the weak-axisymmetric moving-throat branch preserves **all** wall-normalized shapes:
\[
\delta\ln\chi_\alpha^2=0,
\qquad
\delta\ln\Upsilon_r=0,
\qquad
\delta\ln\Lambda_r=0.
\]
Then
\[
\Theta_B=0,
\qquad
\Theta_Z=0,
\qquad
\Theta_N
=
\sum_r \rho_r^{(N)}(-\delta_K)
=
-\delta_K,
\]
because \(\sum_r \rho_r^{(N)}=1\). Therefore
\[
\boxed{
\Xi_{\rm load}=-\delta_K.
}
\]

So a naive common self-similarity branch does **not** kill the remaining linear grouped defect unless the wall baseline itself does not load:
\[
\delta_K=0.
\]

This is a genuine no-go theorem.

It says that on any nontrivial weak-axisymmetric wall-loading branch, simply preserving all wall-normalized microscopic shapes is not enough.
The outgoing transfer sector must carry an **additional wall-loading law**.

That is a sharper statement than Stage 242.

---

## 6. Exact outgoing-load theorem

On conservative-shape-preserving branches, the remaining grouped defect vanishes if and only if
\[
\boxed{
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r
=
\delta_K.
}
\]

So the remaining linear grouped `2.5`PN bottleneck is now completely explicit:

> the weighted outgoing load factor must carry **half** of the wall-baseline logarithmic slope.

A stronger sufficient condition is the lane-by-lane port identity
\[
\boxed{
2\,\delta\ln\Lambda_r=\delta_K
\qquad\text{for every outgoing port }r,
}
\]
which implies \(\Sigma_r^{(N)}=0\) individually.

This is the exact replacement for the old vague statement “the outgoing bundle must be self-similar.”
The true requirement is not mere shape preservation of \(\Lambda_r\), but a very specific wall-loading law.

---

## 7. Direct consequence for the grouped weak-axisymmetric defect pattern

Stage 241 already proved that, on the even-preserving branch,
\[
\Delta_Q^{(20)}=\epsilon\,\Xi_{\rm load},
\qquad
\Delta_Q^{(21)}=\frac{\epsilon}{2}\,\Xi_{\rm load},
\qquad
\Delta_Q^{(22)}=-\epsilon\,\Xi_{\rm load}.
\]
So on conservative-shape-preserving branches this becomes
\[
\boxed{
\Delta_Q^{(20)}
=
\epsilon
\left[
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\delta_K
\right],
}
\]
\[
\boxed{
\Delta_Q^{(21)}
=
\frac{\epsilon}{2}
\left[
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\delta_K
\right],
}
\]
\[
\boxed{
\Delta_Q^{(22)}
=
-\epsilon
\left[
2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r-\delta_K
\right].
}
\]

So on the even-preserving branch, the whole remaining linear grouped `2.5`PN defect is a direct measurement of how the outgoing load factor fails to track the wall baseline.

That is the cleanest live theorem gate so far.

---

## 8. What Stage 243 changes

Before this step, Stage 242 said that the remaining linear grouped defect was the weighted failure of static self-similarity relative to the wall baseline.

After this step, the theorem status is much sharper:

1. the BdG support defect is exactly the drift of the wall-normalized support shape \(\chi_\alpha\);
2. the conservative Maxwell/mixed defect is exactly the drift of the wall-normalized port shape \(\Upsilon_r\);
3. the outgoing-transfer defect is exactly the wall-loading mismatch of the dimensionless load factor \(\Lambda_r=P_r/\Delta_r\);
4. therefore the final linear grouped bottleneck is **not** a generic static-shape problem;
5. on conservative-shape-preserving branches it is purely an **outgoing-load theorem**.

So the next honest derivation step is no longer to recompute \(D_{01}/D_0\) and \(N_{01}/N_0\), or even to re-derive the full static self-similarity story.

The next theorem gate is now:

> compute the grouped weak-axisymmetric drift of the outgoing load factor
> \[
> \Lambda_r=\frac{P_r}{\Delta_r}
> \]
> on the actual moving-throat branch, and test whether its weighted logarithmic slope satisfies
> \[
> 2\sum_r \rho_r^{(N)}\,\delta\ln\Lambda_r=\delta_K.
> \]

That is a real narrowing of the moving-throat PDE bottleneck.
