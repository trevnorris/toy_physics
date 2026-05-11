# Moving-Throat PDE — Stage 199: Exact Pairwise Orbit-Transport Law, Reference-Independent Mismatch Packet, and the Two-Point Orbit-Lock Theorem

## Status

**Exact within the carried Stage 243 orbit/quotient projector calculus, the Stage 249 finite orbit law for the dependent triple, and the exact coherent-branch monomial definitions frozen before Stage 248.**

This stage does **not** introduce a new constitutive law.
It removes the last orbit-side reference-point privilege from the reduced home stretch.

---

## Purpose

Stage 249 solved the dependent microscopic triple
\[
(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W)
\]
exactly once a free microscopic point and an invariant triple were fixed. That was the right finite complement to the Packet-A finish line of Stage 248, but it still privileged one orbit base point.

So the natural next question is already forced:

> can two arbitrary positive microscopic states be compared **directly**, without first choosing a preferred orbit representative?

This stage answers that completely.

The main outputs are:

1. the exact finite **two-point quotient packet**
   \[
   \Delta_{\rm orbit}^{(2\leftarrow1)}
   =
   \bigl(q_{\rm tr}^{(2\leftarrow1)},q_{\rm nt}^{(2\leftarrow1)},q_\eta^{(2\leftarrow1)}\bigr),
   \]
   obtained by applying the Stage 243 drift compiler to the logarithmic pairwise ratio vector,
2. the exact **pairwise orbit-transport factors**
   \[
   \Phi_T^{(2\leftarrow1)},\qquad
   \Phi_K^{(2\leftarrow1)},\qquad
   \Phi_\mu^{(2\leftarrow1)},
   \]
   which determine how the dependent triple must transform between any two points on the same similarity orbit once the five free coordinate ratios are fixed,
3. the exact **reference-independent mismatch triple**
   \[
   (m_T^{(2\leftarrow1)},m_K^{(2\leftarrow1)},m_\mu^{(2\leftarrow1)}),
   \]
4. the exact finite pairwise decomposition
   \[
   \Delta\mathbf x^{(2\leftarrow1)}
   =
   O_{\rm orb}\,\Delta\mathbf x^{(2\leftarrow1)}
   +
   Q_{\rm quot}\,\Delta\mathbf x^{(2\leftarrow1)},
   \]
5. the exact cocycle laws for transport, mismatch, and quotient packet,
6. and the sharp two-point theorem
   \[
   \mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)}
   \iff
   m_T^{(2\leftarrow1)}=m_K^{(2\leftarrow1)}=m_\mu^{(2\leftarrow1)}=1
   \iff
   q_{\rm tr}^{(2\leftarrow1)}=q_{\rm nt}^{(2\leftarrow1)}=q_\eta^{(2\leftarrow1)}=0.
   \]

So Stage 250 is the exact pairwise/reference-independent completion of the Packet-B side of the home stretch.

---

## 1. Two-point logarithmic ratio packet and exact pairwise monomial compiler

Work with two positive microscopic states
\[
\mathbf x^{(1)},\qquad \mathbf x^{(2)},
\]
in the same ordered microscopic coordinates as Stages 243 and 249:
\[
\boxed{
\mathbf x
=
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
}
\]

Define the exact finite pairwise logarithmic ratio vector
\[
\boxed{
\Delta\mathbf x^{(2\leftarrow1)}
:=
\ln\!\left(\frac{\mathbf x^{(2)}}{\mathbf x^{(1)}}\right)
=
\begin{pmatrix}
\ln \mathfrak r_\lambda \\
\ln \mathfrak r_c \\
\ln \mathfrak r_\gamma \\
\ln \mathfrak r_U \\
\ln \mathfrak r_K \\
\ln \mathfrak r_W \\
\ln \mathfrak r_\mu \\
\ln \mathfrak r_T
\end{pmatrix},
}
\]
where the raw pairwise ratios are
\[
\mathfrak r_\lambda:=\frac{\lambda_{W,2}}{\lambda_{W,1}},
\quad
\mathfrak r_c:=\frac{c_{\eta U,2}}{c_{\eta U,1}},
\quad
\mathfrak r_\gamma:=\frac{\gamma_2}{\gamma_1},
\quad
\mathfrak r_U:=\frac{K_{U,2}}{K_{U,1}},
\quad
\mathfrak r_W:=\frac{K_{W,2}}{K_{W,1}},
\]
\[
\mathfrak r_T:=\frac{T_{U,2}}{T_{U,1}},
\qquad
\mathfrak r_K:=\frac{K_{\eta,2}^{(\mathrm{eff})}}{K_{\eta,1}^{(\mathrm{eff})}},
\qquad
\mathfrak r_\mu:=\frac{\mu_{W,2}}{\mu_{W,1}}.
\]

Because the coherent monomials are multiplicative, the Stage 243 linear compiler becomes an **exact finite two-point compiler**:
\[
\boxed{
\mathbf q^{(2\leftarrow1)}
:=
\begin{pmatrix}
q_{\rm tr}^{(2\leftarrow1)}\\
q_{\rm nt}^{(2\leftarrow1)}\\
q_\eta^{(2\leftarrow1)}
\end{pmatrix}
=
M_*\,\Delta\mathbf x^{(2\leftarrow1)}.
}
\]
Equivalently,
\[
\boxed{
q_{\rm tr}^{(2\leftarrow1)}
=
\ln\!\frac{\mathfrak C_{{\rm tr},*}^{(2)}}{\mathfrak C_{{\rm tr},*}^{(1)}},
\qquad
q_{\rm nt}^{(2\leftarrow1)}
=
\ln\!\frac{\mathfrak C_{{\rm nt},*}^{(2)}}{\mathfrak C_{{\rm nt},*}^{(1)}},
\qquad
q_\eta^{(2\leftarrow1)}
=
\ln\!\frac{\epsilon_{\eta,2}}{\epsilon_{\eta,1}}.
}
\]

Written directly in the raw pairwise ratios,
\[
\boxed{
\frac{\mathfrak C_{{\rm tr},*}^{(2)}}{\mathfrak C_{{\rm tr},*}^{(1)}}
=
\left(\frac{\mathfrak r_\gamma\mathfrak r_c}{\mathfrak r_U}\right)^{1+\delta_{U,*}}
\left(\frac{\mathfrak r_T}{\mathfrak r_U}\right)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\frac{\epsilon_{\eta,2}}{\epsilon_{\eta,1}}
=
\frac{\mathfrak r_c^2}{\mathfrak r_U\mathfrak r_K},
}
\]
\[
\boxed{
\frac{\mathfrak C_{{\rm nt},*}^{(2)}}{\mathfrak C_{{\rm nt},*}^{(1)}}
=
\frac{\mathfrak r_\lambda^2\mathfrak r_\mu}{\mathfrak r_K\mathfrak r_W^2}
\left(\frac{\mathfrak r_\gamma^2\mathfrak r_\lambda^2}{\mathfrak r_U\mathfrak r_W}\right)^{E_*}
\left(\frac{\mathfrak r_T}{\mathfrak r_U}\right)^{-F_*}.
}
\]

So the two-point quotient packet is exact and requires no distinguished orbit representative.

---

## 2. Exact pairwise orbit-transport law

Suppose now that the two microscopic states lie on the **same** exact similarity orbit, so that their coherent monomials agree:
\[
\mathfrak C_{{\rm tr},*}^{(2)}=\mathfrak C_{{\rm tr},*}^{(1)},
\qquad
\mathfrak C_{{\rm nt},*}^{(2)}=\mathfrak C_{{\rm nt},*}^{(1)},
\qquad
\epsilon_{\eta,2}=\epsilon_{\eta,1}.
\]

Define the tracking exponent
\[
\boxed{
\alpha_*:=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.
}
\]
Then the same-orbit condition solves exactly for the dependent pairwise ratios.

### 2.1 Tracking transport factor

\[
\boxed{
\Phi_T^{(2\leftarrow1)}
:=
\mathfrak r_U
\left(\frac{\mathfrak r_U}{\mathfrak r_\gamma\mathfrak r_c}\right)^{\alpha_*}.
}
\]
So on the same exact orbit,
\[
\boxed{T_{U,2}=\Phi_T^{(2\leftarrow1)}\,T_{U,1}.}
\]

### 2.2 Dressing transport factor

\[
\boxed{
\Phi_K^{(2\leftarrow1)}
:=
\frac{\mathfrak r_c^2}{\mathfrak r_U}.
}
\]
So on the same exact orbit,
\[
\boxed{K_{\eta,2}^{(\mathrm{eff})}=\Phi_K^{(2\leftarrow1)}\,K_{\eta,1}^{(\mathrm{eff})}.}
\]

### 2.3 Nontracking transport factor

A factorized exact form aligned with the coherent monomials is
\[
\boxed{
\Phi_\mu^{(2\leftarrow1)}
:=
\frac{\Phi_K^{(2\leftarrow1)}(\mathfrak r_W)^2}{(\mathfrak r_\lambda)^2}
\left(\frac{\mathfrak r_\gamma^2\mathfrak r_\lambda^2}{\mathfrak r_U\mathfrak r_W}\right)^{-E_*}
\left(\frac{\Phi_T^{(2\leftarrow1)}}{\mathfrak r_U}\right)^{F_*}.
}
\]
Equivalently, after eliminating `\(\Phi_T\)` and `\(\Phi_K\)`,
\[
\boxed{
\Phi_\mu^{(2\leftarrow1)}
=
(\mathfrak r_c)^{2-F_*\alpha_*}
(\mathfrak r_\gamma)^{-2E_*-F_*\alpha_*}
(\mathfrak r_\lambda)^{-2(1+E_*)}
(\mathfrak r_U)^{-1+E_*+F_*\alpha_*}
(\mathfrak r_W)^{2+E_*}.
}
\]
So on the same exact orbit,
\[
\boxed{\mu_{W,2}=\Phi_\mu^{(2\leftarrow1)}\,\mu_{W,1}.}
\]

### 2.4 Exact pairwise transport theorem

\[
\boxed{\textbf{Theorem (Stage 250 exact pairwise orbit-transport law).}}
\]

**For any two positive microscopic states, once the five free pairwise ratios**
\[
(\mathfrak r_\lambda,\mathfrak r_c,\mathfrak r_\gamma,\mathfrak r_U,\mathfrak r_W)
\]
**are fixed, the exact same-orbit transport of the dependent triple is uniquely determined by**
\[
\Phi_T^{(2\leftarrow1)},\qquad
\Phi_K^{(2\leftarrow1)},\qquad
\Phi_\mu^{(2\leftarrow1)}.
\]

So Stage 249 no longer requires a privileged orbit point. Any microscopic state can serve as the transport source.

---

## 3. Exact reference-independent mismatch triple

For an arbitrary pair of positive microscopic states define the exact reference-independent mismatch ratios by dividing the raw dependent pairwise ratios by the exact same-orbit transport factors:
\[
\boxed{
m_T^{(2\leftarrow1)}:=\frac{\mathfrak r_T}{\Phi_T^{(2\leftarrow1)}},
\qquad
m_K^{(2\leftarrow1)}:=\frac{\mathfrak r_K}{\Phi_K^{(2\leftarrow1)}},
\qquad
m_\mu^{(2\leftarrow1)}:=\frac{\mathfrak r_\mu}{\Phi_\mu^{(2\leftarrow1)}}.
}
\]

Then the pairwise invariant-ratio packet collapses **exactly** to the same functional form found in Stage 249:
\[
\boxed{
\frac{\mathfrak C_{{\rm tr},*}^{(2)}}{\mathfrak C_{{\rm tr},*}^{(1)}}
=
\bigl(m_T^{(2\leftarrow1)}\bigr)^{1+\chi_{0,*}},
}
\]
\[
\boxed{
\frac{\epsilon_{\eta,2}}{\epsilon_{\eta,1}}
=
\frac{1}{m_K^{(2\leftarrow1)}},
}
\]
\[
\boxed{
\frac{\mathfrak C_{{\rm nt},*}^{(2)}}{\mathfrak C_{{\rm nt},*}^{(1)}}
=
\frac{m_\mu^{(2\leftarrow1)}}{m_K^{(2\leftarrow1)}\bigl(m_T^{(2\leftarrow1)}\bigr)^{F_*}}.
}
\]

So the same three packets remain exactly equivalent in the pairwise setting:
\[
\boxed{
(m_T,m_K,m_\mu)
\longleftrightarrow
\left(\frac{\mathfrak C_{{\rm tr},*}^{(2)}}{\mathfrak C_{{\rm tr},*}^{(1)}},
      \frac{\mathfrak C_{{\rm nt},*}^{(2)}}{\mathfrak C_{{\rm nt},*}^{(1)}},
      \frac{\epsilon_{\eta,2}}{\epsilon_{\eta,1}}
\right)
\longleftrightarrow
(q_{\rm tr},q_{\rm nt},q_\eta).
}
\]
But now the packet is completely base-point free.

---

## 4. Exact logarithmic chart and the finite two-point projector split

Write the finite pairwise logarithmic mismatch coordinates as
\[
\tau^{(2\leftarrow1)}:=\ln m_T^{(2\leftarrow1)},
\qquad
\kappa^{(2\leftarrow1)}:=\ln m_K^{(2\leftarrow1)},
\qquad
\mu^{(2\leftarrow1)}:=\ln m_\mu^{(2\leftarrow1)}.
\]
Then exactly as in Stage 249,
\[
\boxed{q_{\rm tr}^{(2\leftarrow1)}=(1+\chi_{0,*})\tau^{(2\leftarrow1)},}
\]
\[
\boxed{q_\eta^{(2\leftarrow1)}=-\kappa^{(2\leftarrow1)},}
\]
\[
\boxed{q_{\rm nt}^{(2\leftarrow1)}=\mu^{(2\leftarrow1)}-\kappa^{(2\leftarrow1)}-F_*\tau^{(2\leftarrow1)}.}
\]

So the Packet-B chart of Stage 249 is already the exact pairwise chart as well.

### 4.1 Exact finite use of the Stage 243 projectors

Because `\(\Delta\mathbf x^{(2\leftarrow1)}\)` is a logarithmic pairwise ratio vector, the Stage 243 projectors apply **exactly** to it:
\[
\boxed{
\Delta\mathbf x^{(2\leftarrow1)}
=
O_{\rm orb}\,\Delta\mathbf x^{(2\leftarrow1)}
+
Q_{\rm quot}\,\Delta\mathbf x^{(2\leftarrow1)}.
}
\]
The exact orbit part is
\[
\boxed{
O_{\rm orb}\,\Delta\mathbf x^{(2\leftarrow1)}
=
\begin{pmatrix}
\ln \mathfrak r_\lambda \\
\ln \mathfrak r_c \\
\ln \mathfrak r_\gamma \\
\ln \mathfrak r_U \\
\ln \Phi_K^{(2\leftarrow1)} \\
\ln \mathfrak r_W \\
\ln \Phi_\mu^{(2\leftarrow1)} \\
\ln \Phi_T^{(2\leftarrow1)}
\end{pmatrix},
}
\]
while the exact quotient-failure part has support only on the dependent triple:
\[
\boxed{
Q_{\rm quot}\,\Delta\mathbf x^{(2\leftarrow1)}
=
\begin{pmatrix}
0\\0\\0\\0\\
\ln m_K^{(2\leftarrow1)}\\
0\\
\ln m_\mu^{(2\leftarrow1)}\\
\ln m_T^{(2\leftarrow1)}
\end{pmatrix}.
}
\]

So Stage 243 is now recognized as an exact finite two-point orbit/failure decomposition.

### 4.2 Exact pairwise restoration map

Restoration of state `2` to the exact similarity orbit through state `1`, while holding the free pairwise ratios fixed, changes only the dependent triple:
\[
\boxed{
T_{U,2}^{\rm restore}
=
T_{U,2}\,e^{-q_{\rm tr}^{(2\leftarrow1)}/(1+\chi_{0,*})}
=
\Phi_T^{(2\leftarrow1)}T_{U,1},
}
\]
\[
\boxed{
K_{\eta,2}^{(\mathrm{eff}),\,\rm restore}
=
K_{\eta,2}^{(\mathrm{eff})}\,e^{q_\eta^{(2\leftarrow1)}}
=
\Phi_K^{(2\leftarrow1)}K_{\eta,1}^{(\mathrm{eff})},
}
\]
\[
\boxed{
\mu_{W,2}^{\rm restore}
=
\mu_{W,2}\,e^{-q_{\rm nt}^{(2\leftarrow1)}+q_\eta^{(2\leftarrow1)}-F_*q_{\rm tr}^{(2\leftarrow1)}/(1+\chi_{0,*})}
=
\Phi_\mu^{(2\leftarrow1)}\mu_{W,1}.
}
\]

So the pairwise orbit-restoration map is exact and algebraic.

---

## 5. Exact composition laws and reference independence

Take three positive microscopic states
\[
\mathbf x^{(1)},\qquad \mathbf x^{(2)},\qquad \mathbf x^{(3)}.
\]
Then the raw pairwise ratios compose multiplicatively, hence the logarithmic pairwise ratio vectors compose additively:
\[
\boxed{
\Delta\mathbf x^{(3\leftarrow1)}
=
\Delta\mathbf x^{(3\leftarrow2)}+\Delta\mathbf x^{(2\leftarrow1)}.
}
\]
Applying the exact pairwise transport factors gives
\[
\boxed{
\Phi_T^{(3\leftarrow1)}=
\Phi_T^{(3\leftarrow2)}\Phi_T^{(2\leftarrow1)},
\qquad
\Phi_K^{(3\leftarrow1)}=
\Phi_K^{(3\leftarrow2)}\Phi_K^{(2\leftarrow1)},
\qquad
\Phi_\mu^{(3\leftarrow1)}=
\Phi_\mu^{(3\leftarrow2)}\Phi_\mu^{(2\leftarrow1)}.
}
\]
Therefore the mismatch ratios compose multiplicatively as well:
\[
\boxed{
m_T^{(3\leftarrow1)}=m_T^{(3\leftarrow2)}m_T^{(2\leftarrow1)},
\qquad
m_K^{(3\leftarrow1)}=m_K^{(3\leftarrow2)}m_K^{(2\leftarrow1)},
\qquad
m_\mu^{(3\leftarrow1)}=m_\mu^{(3\leftarrow2)}m_\mu^{(2\leftarrow1)}.
}
\]
Since `\(M_*\)` is linear,
\[
\boxed{
\mathbf q^{(3\leftarrow1)}
=
\mathbf q^{(3\leftarrow2)}+\mathbf q^{(2\leftarrow1)}.
}
\]
So the Packet-B quotient packet is an exact additive two-point cocycle.

This is the precise sense in which the orbit verdict is now **reference independent**: any intermediate orbit representative gives the same final two-point mismatch packet after composition.

---

## 6. Exact two-point orbit-lock theorem

We can now remove the last orbit-side reference-point privilege completely.

\[
\boxed{\textbf{Theorem (Stage 250 two-point orbit-lock theorem).}}
\]

For any two positive microscopic states `\(\mathbf x^{(1)},\mathbf x^{(2)}\)`, the following are equivalent:

1. they lie on the same exact similarity orbit,
   \[
   \boxed{\mathbf x^{(2)}\in\mathcal G_*\!\cdot\mathbf x^{(1)};}
   \]
2. their dependent raw pairwise ratios equal the exact pairwise orbit-transport factors,
   \[
   \boxed{
   \mathfrak r_T=\Phi_T^{(2\leftarrow1)},
   \qquad
   \mathfrak r_K=\Phi_K^{(2\leftarrow1)},
   \qquad
   \mathfrak r_\mu=\Phi_\mu^{(2\leftarrow1)};
   }
   \]
3. the reference-independent mismatch triple is trivial,
   \[
   \boxed{
   m_T^{(2\leftarrow1)}=m_K^{(2\leftarrow1)}=m_\mu^{(2\leftarrow1)}=1;
   }
   \]
4. the pairwise invariant-ratio packet is trivial,
   \[
   \boxed{
   \frac{\mathfrak C_{{\rm tr},*}^{(2)}}{\mathfrak C_{{\rm tr},*}^{(1)}}
   =
   \frac{\mathfrak C_{{\rm nt},*}^{(2)}}{\mathfrak C_{{\rm nt},*}^{(1)}}
   =
   \frac{\epsilon_{\eta,2}}{\epsilon_{\eta,1}}
   =1;
   }
   \]
5. the pairwise quotient packet vanishes,
   \[
   \boxed{
   q_{\rm tr}^{(2\leftarrow1)}=q_{\rm nt}^{(2\leftarrow1)}=q_\eta^{(2\leftarrow1)}=0.
   }
   \]

Equivalently,
\[
\boxed{
Q_{\rm quot}\,\Delta\mathbf x^{(2\leftarrow1)}=0
\iff
M_*\,\Delta\mathbf x^{(2\leftarrow1)}=0.
}
\]

So the orbit-lock test is now an exact **two-point** statement.

---

## 7. Reduction to Stage 249

If the two states share the same free microscopic coordinates, then
\[
\mathfrak r_\lambda=\mathfrak r_c=\mathfrak r_\gamma=\mathfrak r_U=\mathfrak r_W=1.
\]
Hence
\[
\boxed{
\Phi_T^{(2\leftarrow1)}=\Phi_K^{(2\leftarrow1)}=\Phi_\mu^{(2\leftarrow1)}=1.
}
\]
Therefore
\[
\boxed{
m_T^{(2\leftarrow1)}=\mathfrak r_T,
\qquad
m_K^{(2\leftarrow1)}=\mathfrak r_K,
\qquad
m_\mu^{(2\leftarrow1)}=\mathfrak r_\mu.
}
\]
So if state `1` is chosen to be the explicit orbit point of Stage 249, Stage 250 reduces **exactly** to Stage 249.

This shows that Stage 249 is the fixed-base specialization of the present pairwise law.

---

## 8. What Stage 250 changes in the theorem problem

Stage 249 solved the dependent triple exactly but still privileged one orbit base point.
Stage 250 removes that remaining asymmetry.

### 8.1 The orbit side is now a direct two-point test

The completed moving-throat PDE no longer needs a distinguished orbit representative before the microscopic orbit-lock verdict can be read off. Any two positive microscopic states can be compared directly.

### 8.2 The Stage 243 projectors are now exact finite two-point objects

The orbit/quotient split is no longer merely infinitesimal language. It is the exact decomposition of the finite logarithmic pairwise ratio vector.

### 8.3 The full reduced home stretch is now cleanly split into one Packet-A scalar and one two-point Packet-B test

- Packet A (Stage 248):
  \[
  \chi_Q=1.
  \]
- Packet B (Stage 250):
  \[
  \Delta_{\rm orbit}^{(2\leftarrow1)}=0.
  \]

So the actual moving-throat realization problem is now:

1. realize the Packet-A outgoing finish line,
2. and realize the exact two-point Packet-B orbit-lock condition.

Nothing else survives in the reduced endgame algebra.

---

## 9. Immediate next derivation step

The natural continuation is now the fully reference-free home-stretch theorem:

1. combine the Stage 248 Packet-A scalar finish line with the Stage 250 two-point Packet-B orbit-lock theorem,
2. state the exact reduced closure criterion without privileging any orbit base point,
3. and then feed the actual PDE-selected branch directly into that two-packet compiler.

That is the sharpest next theorem gate after Stage 250.
