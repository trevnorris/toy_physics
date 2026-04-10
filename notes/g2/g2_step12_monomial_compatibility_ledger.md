
# Step 12 — Direct microscopic monomials and the exact compatibility ledger

## Goal

Step 11 compressed the coherent weak-axisymmetric problem to three exact
branch-adapted slippages,
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
```
The next clean move is to remove even that intermediate layer and identify the
**direct microscopic quantities** whose logarithmic drifts are those three scalars.

That is what this step does.

The main result is that the coherent branch is now controlled by three direct
microscopic monomials:
```math
\mathfrak C_{{\rm tr},*},\qquad \mathfrak C_{{\rm nt},*},\qquad \epsilon_\eta,
```
and the full zero-defect condition becomes an explicit three-equation
compatibility ledger for the microscopic grouped drifts
```math
(\lambda_1,c_1,\gamma_1,\kappa_U,\kappa_\eta,\kappa_W,\mu_1,\tau_1).
```

So the continuation point is now even smaller than Step 11 suggested:
it is no longer “track three abstract slippages,” but rather

> determine whether the actual moving-throat branch preserves three direct microscopic kernel monomials.

---

## Inputs carried forward

From Step 11, the coherent-kernel ratios are
```math
\chi_0=\frac{\gamma c_{\eta U}}{K_U},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U},
\qquad
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
```
```math
\epsilon_W=\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}},
\qquad
\frac{Z_W}{\Omega_W^2}
=
\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}.
```
The direct microscopic logarithmic drifts are
```math
\Sigma_\chi=\gamma_1+c_1-\kappa_U,
\qquad
\Sigma_\delta=\tau_1-\kappa_U,
```
```math
\Sigma_Z=2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W,
\qquad
\Sigma_\epsilon=2\gamma_1+2\lambda_1-\kappa_U-\kappa_W,
```
```math
\Sigma_\eta=2c_1-\kappa_U-\kappa_\eta.
```
And Step 11 already defined
```math
\Sigma_{\rm tr}
=
(1+\chi_{0,*})\Sigma_\delta+(1+\delta_{U,*})\Sigma_\chi,
```
```math
\Sigma_{\rm nt}
=
\Sigma_Z
+
E_*\Sigma_\epsilon
-
F_*\Sigma_\delta,
```
with
```math
\epsilon_* = \epsilon_{W,*}\!\left(1-\frac{2}{11}\frac{\delta_{U,*}}{1+\delta_{U,*}}\right).
```
And
```math
E_*=
\frac{2\epsilon_{W,*}}{1-\epsilon_*}\,
\frac{11+9\delta_{U,*}}{11(1+\delta_{U,*})},
```
```math
F_*=
\frac{2\chi_{0,*}}{1+\delta_{U,*}}
+
\frac{4\epsilon_{W,*}\delta_{U,*}}{11(1-\epsilon_*)(1+\delta_{U,*})^2}.
```

---

## Step 12A — Direct microscopic tracking monomial

The tracking scalar is already linear in the two logarithmic drifts
```math
\Sigma_\chi,\Sigma_\delta.
```
So freeze the reference-branch coefficients and define
```math
\boxed{
\mathfrak C_{{\rm tr},*}
:=
\chi_0^{\,1+\delta_{U,*}}
\delta_U^{\,1+\chi_{0,*}}.
}
```
Then
```math
\delta\ln\mathfrak C_{{\rm tr},*}
=
(1+\delta_{U,*})\,\delta\ln\chi_0
+
(1+\chi_{0,*})\,\delta\ln\delta_U
=
\Sigma_{\rm tr}.
```
So the tracking coordinate is not an abstract reduced object anymore:
```math
\boxed{
\delta\ln\mathfrak C_{{\rm tr},*}=\Sigma_{\rm tr}.
}
```

### Explicit microscopic form
```math
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}}.
}
```

---

## Step 12B — Direct microscopic nontracking monomial

The genuine nontracking transfer-shape scalar is also logarithmic in direct
microscopic ratios. Define
```math
\boxed{
\mathfrak C_{{\rm nt},*}
:=
\frac{Z_W}{\Omega_W^2}\,
\epsilon_W^{E_*}\,
\delta_U^{-F_*}.
}
```
Then
```math
\delta\ln\mathfrak C_{{\rm nt},*}
=
\delta\ln\!\left(\frac{Z_W}{\Omega_W^2}\right)
+
E_*\,\delta\ln\epsilon_W
-
F_*\,\delta\ln\delta_U
=
\Sigma_{\rm nt}.
```
So
```math
\boxed{
\delta\ln\mathfrak C_{{\rm nt},*}=\Sigma_{\rm nt}.
}
```

### Explicit microscopic form
```math
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}
{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}
{K_UK_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*}.
}
```

---

## Step 12C — Dressing invariant

The third coordinate is already direct:
```math
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}},
\qquad
\delta\ln\epsilon_\eta=\Sigma_\eta.
}
```

So the full coherent weak-axisymmetric zero-defect theorem now reads
```math
\boxed{
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln\mathfrak C_{{\rm tr},*}=0,
\quad
\delta\ln\mathfrak C_{{\rm nt},*}=0,
\quad
\delta\ln\epsilon_\eta=0.
}
```

That is a real sharpening. The direct weak-axisymmetric branch equations are now
microscopic monomial equations.

---

## Step 12D — Exact monomial-drift matrix

Now collect the microscopic grouped drift vector as
```math
\delta\mathbf x
=
(\lambda_1,\ c_1,\ \gamma_1,\ \kappa_U,\ \kappa_\eta,\ \kappa_W,\ \mu_1,\ \tau_1)^T.
```
Then the direct monomial drifts satisfy
```math
\begin{pmatrix}
\delta\ln\mathfrak C_{{\rm tr},*}\\
\delta\ln\mathfrak C_{{\rm nt},*}\\
\delta\ln\epsilon_\eta
\end{pmatrix}
=
M_*\,
\delta\mathbf x,
```
with
```math
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
```

The useful minor built from the columns
```math
(\tau_1,\kappa_\eta,\mu_1)
```
has determinant
```math
\boxed{
\det M_*^{(\tau_1,\kappa_\eta,\mu_1)}=1+\chi_{0,*}>0.
}
```
So on the physical branch
```math
\operatorname{rank}(M_*)=3,
\qquad
\dim\ker M_*=5.
```

That is the first exact sign that the monomial-rigid branch should be a
five-parameter similarity family rather than a fine-tuned isolated locus.

---

## Step 12E — Exact microscopic compatibility ledger

Setting the three monomial drifts to zero gives an explicit three-equation
compatibility system.

### Tracking compatibility
```math
\boxed{
(1+\chi_{0,*})(\tau_1-\kappa_U)
+
(1+\delta_{U,*})(\gamma_1+c_1-\kappa_U)
=0.
}
```

### Dressing compatibility
```math
\boxed{
2c_1-\kappa_U-\kappa_\eta=0.
}
```

### Nontracking compatibility
```math
\boxed{
2\lambda_1+\mu_1-\kappa_\eta-2\kappa_W
+
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*(\tau_1-\kappa_U)
=0.
}
```

The script solves this system exactly for the three dependent drifts
```math
(\tau_1,\kappa_\eta,\mu_1).
```

### Solved form
```math
\boxed{
\tau_1
=
\kappa_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U),
}
```
```math
\boxed{
\kappa_\eta=2c_1-\kappa_U,
}
```
```math
\boxed{
\mu_1
=
2c_1-\kappa_U+2\kappa_W-2\lambda_1
-
E_*(2\gamma_1+2\lambda_1-\kappa_U-\kappa_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
(\gamma_1+c_1-\kappa_U).
}
```

So the zero-defect branch is now an explicit microscopic rigidity ledger rather
than a broad slippage statement.

---

## Step 12F — Quartic anomaly gate in the monomial coordinates

For the anomaly problem, Step 11 gave
```math
\Xi_1=A_{\rm tr}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
A_{\rm tr}=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})}.
```
In the new direct coordinates that becomes
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,\delta\ln\mathfrak C_{{\rm tr},*}
+
\delta\ln\mathfrak C_{{\rm nt},*}.
}
```
So the carried quartic anomaly target
```math
\Xi_1=\Lambda_1
```
is now
```math
\boxed{
A_{\rm tr}\,\delta\ln\mathfrak C_{{\rm tr},*}
+
\delta\ln\mathfrak C_{{\rm nt},*}
=
\Lambda_1.
}
```
with the carried numerical target
```math
\Lambda_1\approx 0.279605891931464.
```

### Tracking-rigid specialization
If the branch preserves the tracking monomial exactly,
```math
\delta\ln\mathfrak C_{{\rm tr},*}=0,
```
then the entire quartic anomaly target collapses to
```math
\boxed{
\delta\ln\mathfrak C_{{\rm nt},*}=\Lambda_1.
}
```

That is the cleanest direct microscopic statement of the missing anomaly layer so far.

---

## What this changes

This step is important because it removes the last layer of abstraction before
the similarity-orbit geometry.

Before Step 12, the continuation point was still phrased in terms of three exact
branch-adapted slippages.

After Step 12, the continuation point is:

1. one tracking monomial,
2. one nontracking transfer-shape monomial,
3. one dressing ratio,
4. and one explicit `3 x 8` microscopic compatibility matrix.

That is the right base for the exact similarity-orbit closure.

---

## Continuation point

The next clean move is now immediate:

> identify the exact five-parameter monomial-preserving similarity orbit whose tangent space is cut out by the Step-12 compatibility matrix \(M_*\).

That will turn the present linear rigidity ledger into a full geometric closure statement.
