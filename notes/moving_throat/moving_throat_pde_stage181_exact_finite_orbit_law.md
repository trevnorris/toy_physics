# Moving-Throat PDE — Stage 181: Exact Finite Orbit Law for the Dependent Triple, Exact Mismatch Coordinates, and the Orbit-Restoration Theorem

## Status

**Exact within the carried Stage-174 Packet-B hierarchy, the Stage-175 orbit/quotient projector calculus, and the exact coherent-branch monomial definitions already frozen before Stage 180.**

This stage does **not** introduce a new constitutive law.
It is the finite orbit-side complement to the Packet-A retarded finish-line theorem of Stage 180.

---

## Purpose

Stage 180 completed the exact Packet-A retarded statement:
\[
\Delta_{\rm branch}=0
\iff
\chi_Q=1.
\]
But Stage 180 also left the second home-stretch packet untouched:
\[
\Delta_{\rm orbit}=(q_{\rm tr},q_{\rm nt},q_\eta),
\]
which must still vanish for full reduced closure.

Stage 174 had already shown that Packet B can be represented equivalently by
\[
(m_T,m_K,m_\mu),
\qquad
(\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta),
\]
and Stage 175 then upgraded the infinitesimal quotient side into an exact projector calculus on the full eight-dimensional microscopic drift space. What was still missing was the **finite microscopic orbit law itself**:

> given the five free microscopic coordinates and the invariant triple, what are the exact dependent microscopic coordinates on the same similarity orbit, and how does a general candidate branch fail to follow them?

This stage answers that completely.

The main outputs are:

1. the exact finite **single-orbit law** for the dependent triple
   \[
   (T_U, K_\eta^{(\mathrm{eff})}, \mu_W),
   \]
2. the exact three-dimensional **dependent residual mismatch triple**
   \[
   (m_T,m_K,m_\mu),
   \]
3. the exact logarithmic chart
   \[
   q_{\rm tr}=(1+\chi_{0,*})\ln m_T,
   \qquad
   q_\eta=-\ln m_K,
   \qquad
   q_{\rm nt}=\ln m_\mu-\ln m_K-F_*\ln m_T,
   \]
4. the exact restoration map that returns a candidate branch to the same similarity orbit by changing only the dependent triple,
5. and the sharp orbit-lock theorem
   \[
   \Delta_{\rm orbit}=0
   \iff
   m_T=m_K=m_\mu=1.
   \]

So Stage 181 turns the orbit side of the home stretch into a direct finite comparison problem.

---

## 1. Carry-forward split of microscopic coordinates and exact coherent monomials

Work on the positive microscopic state
\[
\boxed{
\mathbf x=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
}
\]
Keep the same five free microscopic coordinates used throughout the Stage-174/175 orbit package:
\[
\boxed{(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})})}
\]
and the same dependent microscopic triple:
\[
\boxed{(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W).}
\]

Let the carried coherent-branch constants be
\[
\chi_{0,*},\qquad \delta_{U,*},\qquad E_*,\qquad F_*,
\]
and let the invariant triple be fixed at
\[
\boxed{\bigl(\mathfrak C_{{\rm tr},*},\ \mathfrak C_{{\rm nt},*},\ \epsilon_{\eta,*}\bigr).}
\]

The exact coherent monomials are
\[
\boxed{
\mathfrak C_{{\rm tr},*}
=
\left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
\left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\mathfrak C_{{\rm nt},*}
=
\frac{\lambda_W^2\mu_W}
{K_\eta^{(\mathrm{eff})}\bigl(K_W^{(\mathrm{eff})}\bigr)^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{E_*}
\left(
\frac{\pi^2T_U}{L^2K_U}
\right)^{-F_*},
}
\]
\[
\boxed{
\epsilon_\eta=
\frac{c_{\eta U}^2}{K_U K_\eta^{(\mathrm{eff})}}.
}
\]

These are the same exact monomials whose infinitesimal drift compiler was encoded by `\(M_*\)` in Stage 175.

---

## 2. Exact finite orbit law for the dependent triple

Because the three monomials are triangular in the dependent triple, they can be solved in closed form.

### 2.1 Exact orbit value of `\(K_\eta^{(\mathrm{eff})}\)`

The third monomial gives immediately
\[
\boxed{
K_\eta^{(\mathrm{orbit})}
=
\frac{c_{\eta U}^2}{K_U\,\epsilon_{\eta,*}}.
}
\]

### 2.2 Exact orbit value of `\(T_U\)`

Substituting the free microscopic point into the tracking monomial gives
\[
\boxed{
T_U^{(\mathrm{orbit})}
=
\frac{L^2K_U}{\pi^2}
\left[
\frac{\mathfrak C_{{\rm tr},*}}
{(\gamma c_{\eta U}/K_U)^{1+\delta_{U,*}}}
\right]^{\!1/(1+\chi_{0,*})}.
}
\]

### 2.3 Exact orbit value of `\(\mu_W\)`

Finally the nontracking monomial fixes
\[
\boxed{
\mu_W^{(\mathrm{orbit})}
=
\frac{\mathfrak C_{{\rm nt},*}\,K_\eta^{(\mathrm{orbit})}\bigl(K_W^{(\mathrm{eff})}\bigr)^2}{\lambda_W^2}
\left(
\frac{\gamma^2\lambda_W^2\sigma}{K_U K_W^{(\mathrm{eff})}}
\right)^{-E_*}
\left(
\frac{\pi^2T_U^{(\mathrm{orbit})}}{L^2K_U}
\right)^{F_*}.
}
\]

So the dependent triple is not abstract any more. For fixed free microscopic point and fixed invariant triple,
\[
\boxed{
(T_U, K_\eta^{(\mathrm{eff})}, \mu_W)
=
\bigl(T_U^{(\mathrm{orbit})},K_\eta^{(\mathrm{orbit})},\mu_W^{(\mathrm{orbit})}\bigr)
}
\]
is the unique microscopic point on the same exact similarity orbit.

### 2.4 Exact finite orbit theorem

\[
\boxed{\textbf{Theorem (Stage 181 finite orbit law).}}
\]

**Given the five free microscopic coordinates and the invariant triple**
\[
\bigl(\mathfrak C_{{\rm tr},*},\mathfrak C_{{\rm nt},*},\epsilon_{\eta,*}\bigr),
\]
**the dependent microscopic triple is determined uniquely and exactly by the formulas above.**

This is the finite microscopic orbit law promised by Stage 175.

---

## 3. Exact dependent residual mismatch triple

Now let a general candidate branch share the same five free microscopic coordinates as the orbit point, but not necessarily the same dependent triple. Define the exact residual mismatch ratios by
\[
\boxed{
T_U=m_T\,T_U^{(\mathrm{orbit})},
\qquad
K_\eta^{(\mathrm{eff})}=m_K\,K_\eta^{(\mathrm{orbit})},
\qquad
\mu_W=m_\mu\,\mu_W^{(\mathrm{orbit})}.
}
\]
These are the Stage-174 residual mismatch ratios, now written as direct ratios to the explicit orbit point.

Then the invariant ratios collapse exactly to
\[
\boxed{
\frac{\mathfrak C_{{\rm tr},*}^{\rm(actual)}}{\mathfrak C_{{\rm tr},*}^{\rm(orbit)}}
= m_T^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\frac{\epsilon_\eta^{\rm(actual)}}{\epsilon_{\eta,*}}
=\frac{1}{m_K},
}
\]
\[
\boxed{
\frac{\mathfrak C_{{\rm nt},*}^{\rm(actual)}}{\mathfrak C_{{\rm nt},*}^{\rm(orbit)}}
=\frac{m_\mu}{m_K m_T^{F_*}}.
}
\]

So the entire finite orbit-side branch-selection problem is exactly three-dimensional:
\[
\boxed{(m_T,m_K,m_\mu).}
\]

If desired, the invariant-ratio packet of Stage 174 is therefore
\[
\boxed{
\mathfrak R_{\rm tr}=m_T^{1+\chi_{0,*}},
\qquad
\mathfrak R_{\rm nt}=\frac{m_\mu}{m_K m_T^{F_*}},
\qquad
\mathfrak R_\eta=\frac{1}{m_K}.
}
\]

---

## 4. Exact logarithmic chart and agreement with Stage 175

Define the logarithmic mismatch coordinates
\[
\tau:=\ln m_T,
\qquad
\kappa:=\ln m_K,
\qquad
\mu:=\ln m_\mu.
\]
Then the quotient coordinates are exactly
\[
\boxed{q_{\rm tr}=(1+\chi_{0,*})\tau,}
\]
\[
\boxed{q_\eta=-\kappa,}
\]
\[
\boxed{q_{\rm nt}=\mu-\kappa-F_*\tau.}
\]

So the Stage-174 Packet-B logarithmic chart is not merely a finite reparameterization by analogy. It is the exact logarithmic chart of the finite mismatch triple.

### 4.1 Direct connection to the Stage-175 drift compiler

Take the pure dependent mismatch drift vector in the ordered Stage-175 microscopic basis
\[
\Delta\mathbf x_{\rm mis}=
\begin{pmatrix}
0\\0\\0\\0\\ \kappa \\0\\ \mu \\ \tau
\end{pmatrix}.
\]
Then the Stage-175 quotient map gives
\[
M_*\Delta\mathbf x_{\rm mis}
=
\begin{pmatrix}
(1+\chi_{0,*})\tau\\[4pt]
\mu-\kappa-F_*\tau\\[4pt]
-\kappa
\end{pmatrix}
=
\begin{pmatrix}
q_{\rm tr}\\ q_{\rm nt}\\ q_\eta
\end{pmatrix}.
\]
So the Stage-175 first-order formulas are not just infinitesimal approximations of the finite mismatch language. They are the exact logarithmic chart of the finite mismatch ratios.

---

## 5. Exact restoration map

Given the quotient packet
\[
\Delta_{\rm orbit}=(q_{\rm tr},q_{\rm nt},q_\eta),
\]
restoration to the same exact similarity orbit is achieved by changing only the dependent triple:
\[
\boxed{
T_U^{(\mathrm{restore})}
=
T_U\,e^{-q_{\rm tr}/(1+\chi_{0,*})},
}
\]
\[
\boxed{
K_\eta^{(\mathrm{eff}),\,\mathrm{restore}}
=
K_\eta^{(\mathrm{eff})}\,e^{q_\eta},
}
\]
\[
\boxed{
\mu_W^{(\mathrm{restore})}
=
\mu_W\,e^{-q_{\rm nt}+q_\eta-F_*q_{\rm tr}/(1+\chi_{0,*})}.
}
\]
By construction,
\[
\boxed{
T_U^{(\mathrm{restore})}=T_U^{(\mathrm{orbit})},
\qquad
K_\eta^{(\mathrm{eff}),\,\mathrm{restore}}=K_\eta^{(\mathrm{orbit})},
\qquad
\mu_W^{(\mathrm{restore})}=\mu_W^{(\mathrm{orbit})}.
}
\]
So orbit restoration is exact and algebraic.

---

## 6. Exact finite orbit-lock theorem

The orbit-side criterion is now completely explicit.

\[
\boxed{\textbf{Theorem (Stage 181 finite orbit-lock theorem).}}
\]

**Within the carried coherent-branch hierarchy, a candidate branch with fixed free microscopic coordinates lies on the exact similarity orbit determined by the invariant triple if and only if any one of the following equivalent conditions holds:**
\[
\boxed{m_T=m_K=m_\mu=1,}
\]
\[
\boxed{\mathfrak R_{\rm tr}=\mathfrak R_{\rm nt}=\mathfrak R_\eta=1,}
\]
\[
\boxed{q_{\rm tr}=q_{\rm nt}=q_\eta=0.}
\]
Equivalently,
\[
\boxed{\Delta_{\rm orbit}=0.}
\]

So Stage 181 is the exact Packet-B complement of Stage 180:

- Stage 180 says the Packet-A retarded finish line is
  \[
  \chi_Q=1.
  \]
- Stage 181 says the Packet-B orbit-side finish line is
  \[
  \Delta_{\rm orbit}=0.
  \]

Together they sharpen the full reduced home stretch to
\[
\boxed{
\text{full reduced closure}
\iff
\chi_Q=1
\quad\text{and}\quad
\Delta_{\rm orbit}=0.
}
\]

---

## 7. What Stage 181 changes in the theorem problem

Stage 175 gave the exact infinitesimal orbit/failure split, but the actual finite orbit point was still implicit. Stage 181 removes that last ambiguity.

### 7.1 The finite similarity orbit is now explicit

The completed moving-throat PDE no longer needs an abstract “orbit-lock surface” to be compared by hand. The exact orbit point in the dependent triple is already known algebraically once the free coordinates and invariant triple are supplied.

### 7.2 The quotient coordinates are exact finite mismatch coordinates

The finite quotient packet
\[
(q_{\rm tr},q_{\rm nt},q_\eta)
\]
is now recognized as the exact logarithmic chart of the dependent mismatch triple
\[
(m_T,m_K,m_\mu).
\]

### 7.3 The remaining microscopic realization problem is localized completely

Once the PDE returns the actual microscopic state,
all residual orbit failure is localized exactly into the dependent triple:
\[
T_U,
\qquad
K_\eta^{(\mathrm{eff})},
\qquad
\mu_W.
\]
The free microscopic coordinates carry the similarity-orbit transport only; they do not themselves signal orbit-lock failure.

---

## 8. Immediate next derivation step

The next clean continuation is now fully finite and microscopic:

1. choose the actual free microscopic point returned by the moving-throat PDE,
2. compute the exact orbit-predicted dependent triple from the Stage-181 formulas,
3. compare the actual dependent triple to that orbit prediction,
4. and read off
   \[
   (m_T,m_K,m_\mu)
   \quad\text{or equivalently}\quad
   (q_{\rm tr},q_{\rm nt},q_\eta).
   \]

That is the sharpest direct orbit-lock test available after Stage 181.
