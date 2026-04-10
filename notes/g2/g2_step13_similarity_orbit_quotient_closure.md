# Step 13 — Exact similarity orbit and quotient closure

## Goal

Step 12 identified the three **direct microscopic monomials**
```math
\mathfrak C_{{\rm tr},*},\qquad \mathfrak C_{{\rm nt},*},\qquad \epsilon_\eta,
```
and proved that their logarithmic drifts are exactly
```math
\Sigma_{\rm tr},\qquad \Sigma_{\rm nt},\qquad \Sigma_\eta.
```

The remaining question was then unavoidable:

> what is the exact microscopic redundancy that leaves those three monomials unchanged?

This step answers that question completely.

It proves that the coherent weak-axisymmetric zero-defect branch is the tangent
space of an exact **five-parameter multiplicative similarity orbit**, and that the
three monomials above are the exact **quotient coordinates**.

For the g-2 derivation, that means something very concrete:

- the anomaly problem is no longer an 8-variable microscopic drift problem,
- it is exactly a 3-coordinate quotient problem,
- and the direct quartic anomaly gate itself lives only in a 2-dimensional
  quotient plane.

---

## Inputs carried forward

From Step 12, define the direct microscopic state vector
```math
x=(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U),
```
with all entries positive on the coherent branch.

The three exact direct monomials are
```math
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
```
```math
\boxed{
\mathfrak C_{{\rm nt},*}
=
\left(\frac{\lambda_W^2\mu_W}{K_\eta^{(\mathrm{eff})}(K_W^{(\mathrm{eff})})^2}\right)
\left(\frac{\gamma^2\lambda_W^2\sigma}{K_UK_W^{(\mathrm{eff})}}\right)^{E_*}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{-F_*},
}
```
```math
\boxed{
\epsilon_\eta=\frac{c_{\eta U}^2}{K_UK_\eta^{(\mathrm{eff})}}.
}
```
And the exact monomial-drift matrix from Step 12 is
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

---

## Step 13A — The finite fibre equations are controlled by the same matrix

Take two positive microscopic states `x` and `\widetilde x`, and define their
finite logarithmic ratio vector
```math
\Delta\mathbf x=
\begin{pmatrix}
\Delta_\lambda,
\Delta_c,
\Delta_\gamma,
\Delta_U,
\Delta_\eta,
\Delta_W,
\Delta_\mu,
\Delta_T
\end{pmatrix}^T
=
\begin{pmatrix}
\ln(\widetilde\lambda_W/\lambda_W),
\ln(\widetilde c_{\eta U}/c_{\eta U}),
\ln(\widetilde\gamma/\gamma),
\ln(\widetilde K_U/K_U),
\ln(\widetilde K_\eta/K_\eta),
\ln(\widetilde K_W/K_W),
\ln(\widetilde\mu_W/\mu_W),
\ln(\widetilde T_U/T_U)
\end{pmatrix}^T.
```
Then the exact invariant log-ratios are
```math
q_{\rm tr}:=\ln\frac{\widetilde{\mathfrak C}_{{\rm tr},*}}{\mathfrak C_{{\rm tr},*}},
\qquad
q_{\rm nt}:=\ln\frac{\widetilde{\mathfrak C}_{{\rm nt},*}}{\mathfrak C_{{\rm nt},*}},
\qquad
q_\eta:=\ln\frac{\widetilde\epsilon_\eta}{\epsilon_\eta}.
```
Because the three invariants are monomials, their exact finite log-ratios are
still linear in `\Delta\mathbf x`:
```math
\boxed{
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_*\,\Delta\mathbf x.
}
```

This is the first key result of the step.
The matrix that governed the infinitesimal compatibility ledger in Step 12 also
controls the **exact finite invariant fibres**.

---

## Step 13B — Exact five-parameter similarity orbit

Choose the five free finite logarithmic co-scalings
```math
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W)
```
for
```math
(\lambda_W,c_{\eta U},\gamma,K_U,K_W^{(\mathrm{eff})}).
```
Setting the invariant log-ratios to zero,
```math
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
```
gives the exact solved dependent shifts
```math
\boxed{
\Delta_\eta=2\Delta_c-\Delta_U,
}
```
```math
\boxed{
\Delta_T
=
\Delta_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U),
}
```
```math
\boxed{
\Delta_\mu
=
2\Delta_c-\Delta_U+2\Delta_W-2\Delta_\lambda
-
E_*(2\Delta_\gamma+2\Delta_\lambda-\Delta_U-\Delta_W)
-
F_*\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U).
}
```
Exponentiating gives the exact five-parameter multiplicative similarity orbit
`\mathcal G_*`.

So the zero-defect microscopic branch is not an isolated compatibility locus.
It is an exact codimension-3 similarity family.

---

## Step 13C — Tangent generator matrix

Writing the free orbit coordinates as
```math
g=(g_1,g_2,g_3,g_4,g_5)^T
```
for
```math
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W),
```
the exact tangent generator matrix is
```math
\boxed{
G=
\begin{pmatrix}
1 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 1 & 0\\
0 & 2 & 0 & -1 & 0\\
0 & 0 & 0 & 0 & 1\\
-2(1+E_*) & 2-F_*\alpha & -2E_*-F_*\alpha & -1+E_*+F_*\alpha & 2+E_*\\
0 & -\alpha & -\alpha & 1+\alpha & 0
\end{pmatrix},
\qquad
\alpha=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.
}
```
The script verifies
```math
\boxed{M_*G=0.}
```
So the tangent space of the similarity orbit is exactly the kernel of the
monomial-drift map.

That is the infinitesimal version of the orbit theorem.

---

## Step 13D — Canonical exact quotient section

Now solve the general finite fibre equation
```math
M_*\Delta\mathbf x=q,
\qquad
q=(q_{\rm tr},q_{\rm nt},q_\eta)^T,
```
using the same five orbit coordinates as free variables.
The dependent coordinates become
```math
\boxed{
\Delta_\eta=2\Delta_c-\Delta_U-q_\eta,
}
```
```math
\boxed{
\Delta_T
=
\Delta_U
-
\frac{1+\delta_{U,*}}{1+\chi_{0,*}}(\Delta_\gamma+\Delta_c-\Delta_U)
+
\frac{q_{\rm tr}}{1+\chi_{0,*}},
}
```
```math
\boxed{
\Delta_\mu
=
\Delta_\mu^{\rm orbit}
+
q_{\rm nt}-q_\eta+
\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}.
}
```
So the exact finite microscopic state splits into
```math
\boxed{
\Delta\mathbf x = Gg + Sq,
}
```
where the chosen right-inverse section is
```math
\boxed{
S=
\begin{pmatrix}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & -1\\
0 & 0 & 0\\
\dfrac{F_*}{1+\chi_{0,*}} & 1 & -1\\
\dfrac{1}{1+\chi_{0,*}} & 0 & 0
\end{pmatrix},
\qquad
M_*S=I_3.
}
```

This is the exact quotient closure in the form most useful for g-2.
It says that the eight microscopic log-ratios split into

- **five similarity-orbit directions** `g`, and
- **three exact quotient coordinates** `q`.

---

## Step 13E — Canonical finite quotient representative

The chosen gauge slice keeps the five similarity coordinates fixed and moves only
along the three quotient directions. In finite multiplicative form this means
```math
\boxed{
K_\eta^{(\mathrm{eff})}\mapsto e^{-q_\eta}K_\eta^{(\mathrm{eff})},
}
```
```math
\boxed{
T_U\mapsto e^{q_{\rm tr}/(1+\chi_{0,*})}T_U,
}
```
```math
\boxed{
\mu_W
\mapsto
\exp\!\left(q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}\right)\mu_W,
}
```
with the other five microscopic variables held fixed.

The script verifies directly that this representative satisfies
```math
\Delta\ln\mathfrak C_{{\rm tr},*}=q_{\rm tr},
\qquad
\Delta\ln\mathfrak C_{{\rm nt},*}=q_{\rm nt},
\qquad
\Delta\ln\epsilon_\eta=q_\eta.
```

So the quotient coordinates are not abstract. They have an explicit finite
microscopic representative.

---

## Step 13F — Observables in exact quotient coordinates

Now the coherent weak-axisymmetric observables become completely transparent.
Using the carried triangular normal form,
```math
\boxed{
\Theta_1
=
-
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})}
q_{\rm tr},
}
```
```math
\boxed{
\Xi_1
=
A_{\rm tr}\,q_{\rm tr}+q_{\rm nt},
\qquad
A_{\rm tr}=\frac{2\chi_{0,*}}{(1+\chi_{0,*})(1+\delta_{U,*})},
}
```
```math
\boxed{
\mathcal R_1+\Xi_1
=
-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
}
```

This is one of the most useful simplifications so far.
The direct quartic anomaly gate depends only on
```math
(q_{\rm tr},q_{\rm nt}),
```
while the third exact quotient coordinate `q_\eta` is orthogonal to the direct
`\Xi_1` match at this order.

So the g-2 problem has collapsed from

- 8 microscopic grouped drifts,

to

- a 3-dimensional exact quotient,

and then further, for the direct quartic match itself, to

- a 2-dimensional quotient plane.

---

## Step 13G — Canonical quartic anomaly representatives

The carried quartic anomaly target is
```math
\Xi_1=\Lambda_1.
```
In exact quotient coordinates this is simply
```math
\boxed{
A_{\rm tr}q_{\rm tr}+q_{\rm nt}=\Lambda_1.
}
```
So on the chosen quotient section the general canonical representative is
```math
\Delta\mathbf x
=
S\begin{pmatrix}
q_{\rm tr}\\
\Lambda_1-A_{\rm tr}q_{\rm tr}\\
q_\eta
\end{pmatrix}.
```
That is the exact gauge-fixed form of the whole quartic matching family.

### Simplest direct slice
Take the tracking-rigid and dressing-rigid slice
```math
q_{\rm tr}=0,
\qquad
q_\eta=0.
```
Then
```math
q_{\rm nt}=\Lambda_1,
```
and the canonical representative reduces to
```math
\boxed{
\Delta\mathbf x
=
(0,0,0,0,0,0,\Lambda_1,0)^T.
}
```
So in the chosen quotient gauge, the carried quartic anomaly target is
represented entirely by a **pure `\mu_W` drift**.

That does **not** mean the real moving-throat branch must physically vary only
`\mu_W`. It means that after modding out the exact five-parameter similarity
redundancy, every tracking-rigid, dressing-rigid quartic match can be represented
in that minimal way.

This is the cleanest microscopic simplification of the anomaly problem so far.

---

## What this changes

This step is the point where the moving-throat quotient geometry becomes directly
useful for g-2.

Before Step 13, the direct anomaly statement still lived inside an 8-variable
microscopic drift ledger.

After Step 13:

1. the exact microscopic redundancy is known,
2. the exact quotient coordinates are explicit,
3. the observables are linear in those quotient coordinates,
4. the direct quartic anomaly gate uses only two of the three,
5. and the chosen canonical slice gives an explicit minimal microscopic
   representative of the whole matching family.

So the next step no longer needs to “simplify the similarity-orbit algebra.”
That is done.
The next step should instead use the quotient closure to decide which exact
quotient trajectory the actual moving-throat branch is most naturally taking.

---

## Continuation point

The clean next move is now:

> use the exact quotient section to choose a physically meaningful branch condition — for example tracking-rigid, dressing-rigid, or minimum-norm in the quotient plane — and translate that into the simplest microscopic law for the missing quartic anomaly layer.

At this point the direct g-2 problem is genuinely small enough to test specific closures rather than reorganizing the algebra again.
