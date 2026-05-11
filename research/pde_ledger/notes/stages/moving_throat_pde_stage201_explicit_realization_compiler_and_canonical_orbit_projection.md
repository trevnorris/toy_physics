# Moving-Throat PDE — Stage 201: Exact Realization Compiler, Canonical Orbit Projection, and the First Explicit Four-Scalar Audit

## Status

**Exact within the carried Stage 243 orbit/quotient projector calculus, the Stage 249 finite restoration law, and the Stage 251 exact reference-free full home-stretch theorem.**

This stage does **not** introduce a new constitutive law.
It turns the already-carried home-stretch theorem into an explicit realizability compiler that can be applied directly to any candidate moving-throat branch.

---

## Purpose

Stage 251 completed the reduced home stretch in theorem form:
\[
\Delta_{\rm full}=0
\iff
\chi_Q=1
\quad\text{and}\quad
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
\]

That already answers the logical question.
But it still leaves the practical audit question open:

> given one candidate microscopic state from the completed moving-throat PDE, how do we test it against the target branch **without** choosing an arbitrary orbit representative, and how do we identify the unique microscopic correction that would return it to the target orbit if only the quotient packet failed?

This stage answers that exactly.

The main outputs are:

1. the **intrinsic** four-scalar realization packet written directly against the target monomial values,
2. the exact **canonical dependent-triple repair vector**
   \[
   \Delta\mathbf x_{\rm rep},
   \]
3. the exact **canonical orbit projection**
   \[
   \Pi^{\rm can}_{\mathcal O_*}(\mathbf x),
   \]
   which lands on the target similarity orbit while preserving the five free similarity coordinates,
4. the exact uniqueness theorem that this projected state is the **only** point on the target orbit with the same free quintuple,
5. and the explicit audit criterion
   \[
   \boxed{
   \mathbf x\in\mathcal Z_*
   \iff
   \chi_Q(\mathbf x)=1
   \ \text{and}\
   \Pi^{\rm can}_{\mathcal O_*}(\mathbf x)=\mathbf x,
   }
   \]
   where `\(\mathcal Z_*:=\{\mathbf x\in\mathcal O_*:\chi_Q(\mathbf x)=1\}\)` is the fully closed target set.

So Stage 252 is the first exact **operational** realization test of the completed branch.

---

## 1. Carry-forward target data and intrinsic branch ratios

Let the positive microscopic state be
\[
\boxed{
\mathbf x:=
\bigl(
\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U
\bigr).
}
\]

Carry forward the three exact coherent monomials from the Stage 238 / Stage 243 quotient hierarchy:
\[
\mathfrak C_{{\rm tr},*}(\mathbf x),
\qquad
\mathfrak C_{{\rm nt},*}(\mathbf x),
\qquad
\epsilon_\eta(\mathbf x).
\]

Fix the target orbit by the exact target values
\[
\boxed{
\mathfrak C_{{\rm tr},*}^{\rm target},
\qquad
\mathfrak C_{{\rm nt},*}^{\rm target},
\qquad
\epsilon_{\eta,*}^{\rm target}.
}
\]
Then
\[
\boxed{
\mathcal O_*=
\left\{
\mathbf x>0:
\mathfrak C_{{\rm tr},*}(\mathbf x)=\mathfrak C_{{\rm tr},*}^{\rm target},\ 
\mathfrak C_{{\rm nt},*}(\mathbf x)=\mathfrak C_{{\rm nt},*}^{\rm target},\ 
\epsilon_\eta(\mathbf x)=\epsilon_{\eta,*}^{\rm target}
\right\}.
}
\]

Define the exact intrinsic multiplicative ratios
\[
\boxed{
\mathfrak R_{\rm tr}(\mathbf x)
:=
\frac{\mathfrak C_{{\rm tr},*}(\mathbf x)}{\mathfrak C_{{\rm tr},*}^{\rm target}},
\qquad
\mathfrak R_{\rm nt}(\mathbf x)
:=
\frac{\mathfrak C_{{\rm nt},*}(\mathbf x)}{\mathfrak C_{{\rm nt},*}^{\rm target}},
\qquad
\mathfrak R_\eta(\mathbf x)
:=
\frac{\epsilon_\eta(\mathbf x)}{\epsilon_{\eta,*}^{\rm target}}.
}
\]

Their exact logarithmic versions are
\[
\boxed{
q_{\rm tr}(\mathbf x)=\ln \mathfrak R_{\rm tr}(\mathbf x),
\qquad
q_{\rm nt}(\mathbf x)=\ln \mathfrak R_{\rm nt}(\mathbf x),
\qquad
q_\eta(\mathbf x)=\ln \mathfrak R_\eta(\mathbf x).
}
\]

So the target orbit condition is already intrinsic:
\[
\boxed{
\mathbf x\in\mathcal O_*
\iff
\mathfrak R_{\rm tr}=\mathfrak R_{\rm nt}=\mathfrak R_\eta=1
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0.
}
\]

No orbit witness is needed.

---

## 2. Exact intrinsic four-scalar realization packet

The Packet-A scalar remains
\[
\boxed{
\Delta_Q(\mathbf x):=\chi_Q(\mathbf x)-1.
}
\]

Define the exact intrinsic full realization packet by
\[
\boxed{
\Delta_{\rm real}^{\rm int}(\mathbf x\mid \mathcal Z_*)
:=
\Bigl(
\Delta_Q(\mathbf x),\
q_{\rm tr}(\mathbf x),\
q_{\rm nt}(\mathbf x),\
q_\eta(\mathbf x)
\Bigr),
}
\]
where
\[
\boxed{
\mathcal Z_*:=\{\mathbf x\in\mathcal O_*:\chi_Q(\mathbf x)=1\}.
}
\]

Equivalently, the exact multiplicative chart is
\[
\boxed{
\mathcal V_{\rm real}^{\rm int}(\mathbf x\mid \mathcal Z_*)
:=
\Bigl(
\chi_Q(\mathbf x),\
\mathfrak R_{\rm tr}(\mathbf x),\
\mathfrak R_{\rm nt}(\mathbf x),\
\mathfrak R_\eta(\mathbf x)
\Bigr).
}
\]

By Stage 249/250, the exact mismatch chart is
\[
\boxed{
\mathcal M_{\rm real}^{\rm int}(\mathbf x\mid \mathcal Z_*)
:=
\Bigl(
\chi_Q(\mathbf x),\
m_T(\mathbf x),\
m_K(\mathbf x),\
m_\mu(\mathbf x)
\Bigr),
}
\]
with
\[
\boxed{
m_T=\exp\!\left(\frac{q_{\rm tr}}{1+\chi_{0,*}}\right)
=
\mathfrak R_{\rm tr}^{1/(1+\chi_{0,*})},
}
\]
\[
\boxed{
m_K=e^{-q_\eta}=\mathfrak R_\eta^{-1},
}
\]
\[
\boxed{
m_\mu=\exp\!\left(q_{\rm nt}-q_\eta+\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}\right)
=
\mathfrak R_{\rm nt}\,\mathfrak R_\eta^{-1}\,\mathfrak R_{\rm tr}^{F_* /(1+\chi_{0,*})}.
}
\]

So the exact closure criterion is already
\[
\boxed{
\Delta_{\rm real}^{\rm int}=0
\iff
\mathcal V_{\rm real}^{\rm int}=(1,1,1,1)
\iff
\mathcal M_{\rm real}^{\rm int}=(1,1,1,1)
\iff
\mathbf x\in\mathcal Z_*.
}
\]

This is the intrinsic version of Stage 251.

---

## 3. Exact canonical dependent-triple repair vector

Carry forward the Stage 243 exact quotient map
\[
\boxed{
\mathbf q=
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_*\,\Delta\mathbf x,
}
\]
on the ordered microscopic drift basis
\[
\Delta\mathbf x=
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_{K_\eta},\Delta_W,\Delta_\mu,\Delta_T)^T.
\]

Carry forward the exact canonical quotient section
\[
\boxed{
S_{(T,K_\eta,\mu)}.
}
\]
Then the exact canonical quotient representative is
\[
\Delta\mathbf x_{\rm quot}=S_{(T,K_\eta,\mu)}\,\mathbf q,
\]
and the exact canonical **repair** vector is
\[
\boxed{
\Delta\mathbf x_{\rm rep}:=-S_{(T,K_\eta,\mu)}\,\mathbf q.
}
\]

Using the explicit Stage 243 section, this becomes
\[
\boxed{
\Delta\mathbf x_{\rm rep}
=
\begin{pmatrix}
0\\
0\\
0\\
0\\
q_\eta\\
0\\
-\;q_{\rm nt}+q_\eta-\dfrac{F_*}{1+\chi_{0,*}}\,q_{\rm tr}\\[8pt]
-\dfrac{q_{\rm tr}}{1+\chi_{0,*}}
\end{pmatrix}.
}
\]

So the repair demand has support only on the dependent triple
\[
(T,\ K_\eta^{(\mathrm{eff})},\ \mu_W).
\]

### 3.1 Exact cancellation property

Because
\[
M_*S_{(T,K_\eta,\mu)}=I_3,
\]
the repair vector kills the quotient packet exactly:
\[
\boxed{
M_*\,\Delta\mathbf x_{\rm rep}
=
-\mathbf q.
}
\]
So for any candidate microscopic drift,
\[
\boxed{
M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm rep})=0.
}
\]

This is the exact microscopic meaning of “repair to the target orbit while keeping the free similarity coordinates fixed.”

---

## 4. Exact canonical orbit projection in multiplicative variables

Write the intrinsic Packet-B ratios
\[
\mathfrak R_{\rm tr}=e^{q_{\rm tr}},
\qquad
\mathfrak R_{\rm nt}=e^{q_{\rm nt}},
\qquad
\mathfrak R_\eta=e^{q_\eta}.
\]

Exponentiating the repair vector gives the exact canonical orbit projection
\[
\boxed{
\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)
=
\bigl(
\lambda_W,\ 
c_{\eta U},\
\gamma,\
K_U,\
K_\eta^{(\mathrm{eff})}\,\mathfrak R_\eta,\
K_W^{(\mathrm{eff})},\
\mu_W\,\mathfrak R_{\rm nt}^{-1}\mathfrak R_\eta\,\mathfrak R_{\rm tr}^{-F_* /(1+\chi_{0,*})},\
T_U\,\mathfrak R_{\rm tr}^{-1/(1+\chi_{0,*})}
\bigr).
}
\]

Equivalently,
\[
\boxed{
T_U^{(\mathrm{proj})}
=
T_U\,e^{-q_{\rm tr}/(1+\chi_{0,*})},
}
\]
\[
\boxed{
K_\eta^{(\mathrm{eff}),\,\mathrm{proj}}
=
K_\eta^{(\mathrm{eff})}\,e^{q_\eta},
}
\]
\[
\boxed{
\mu_W^{(\mathrm{proj})}
=
\mu_W\,e^{-q_{\rm nt}+q_\eta-F_*q_{\rm tr}/(1+\chi_{0,*})}.
}
\]

So Stage 249’s restoration law is exactly the canonical target-orbit projection of the Stage 251 intrinsic packet.

---

## 5. Exact target-orbit projection theorem

\[
\boxed{\textbf{Theorem (Stage 252 canonical orbit projection).}}
\]

For any positive microscopic state `\(\mathbf x\)` and any target orbit `\(\mathcal O_*\)` in the carried coherent quotient hierarchy:

1. the projected state lies on the target orbit,
   \[
   \boxed{
   \Pi^{\rm can}_{\mathcal O_*}(\mathbf x)\in\mathcal O_*,
   }
   \]
2. it preserves the five free similarity coordinates
   \[
   (\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}),
   \]
3. and it is the **unique** state on `\(\mathcal O_*\)` with those five free coordinates fixed.

### 5.1 Proof

Because
\[
M_*\,\Delta\mathbf x_{\rm rep}=-\mathbf q,
\]
the corrected drift satisfies
\[
M_*(\Delta\mathbf x+\Delta\mathbf x_{\rm rep})=0.
\]
So the corrected state has
\[
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]
which is equivalent to lying on `\(\mathcal O_*\)`.

The support of `\(\Delta\mathbf x_{\rm rep}\)` lies only on the dependent triple, so the free quintuple is unchanged.

Finally, uniqueness follows from the exact pivot block on `(T,K_\eta,\mu)`: once the five free coordinates are fixed, the condition `\(M_*\Delta\mathbf x=0\)` determines the dependent triple uniquely.

---

## 6. Fixed-point criterion and explicit realization test

The projection is a fixed point exactly on the target orbit:
\[
\boxed{
\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)=\mathbf x
\iff
q_{\rm tr}=q_{\rm nt}=q_\eta=0
\iff
\mathbf x\in\mathcal O_*.
}
\]

Therefore the Stage 251 full closure theorem can be rewritten in the most explicit audit form so far:
\[
\boxed{
\mathbf x\in\mathcal Z_*
\iff
\chi_Q(\mathbf x)=1
\ \text{and}\
\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)=\mathbf x.
}
\]

Equivalently,
\[
\boxed{
\mathbf x\in\mathcal Z_*
\iff
\Delta_{\rm real}^{\rm int}(\mathbf x\mid\mathcal Z_*)=0.
}
\]

So the first exact realizability audit of a completed PDE branch is now:

1. compute the Packet-A scalar `\(\chi_Q(\mathbf x)\)`,
2. compute the three intrinsic monomial ratios `\((\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta)\)`,
3. form the four-scalar packet
   \[
   (\chi_Q-1,\ln\mathfrak R_{\rm tr},\ln\mathfrak R_{\rm nt},\ln\mathfrak R_\eta),
   \]
4. or equivalently compute the canonical orbit projection and check whether it is a fixed point.

No orbit witness is required anywhere.

---

## 7. Exact same-free-quintuple interpretation

The projection theorem gives a very concrete microscopic reading of the Packet-B failure.

Take a candidate state `\(\mathbf x\)`. Freeze the five free coordinates
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_W^{(\mathrm{eff})}).
\]
Then:

- there exists exactly one point on the target orbit with those frozen free coordinates,
- that point is `\(\Pi^{\rm can}_{\mathcal O_*}(\mathbf x)\)`,
- and the three quotient coordinates
  \[
  (q_{\rm tr},q_{\rm nt},q_\eta)
  \]
  are exactly the logarithmic mismatch that tells us how far the dependent triple
  \[
  (T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W)
  \]
  must move to land on that orbit.

So Packet B is now not only a quotient theorem.
It is an explicit same-free-quintuple orbit locator.

---

## 8. First-order linearized version

At first order around a target-orbit witness,
\[
\Delta\mathbf x_{\rm rep}^{\rm lin}
=
-S_{(T,K_\eta,\mu)}\,M_*\,\Delta\mathbf x,
\]
and the full linearized realization packet is simply
\[
\boxed{
\Delta_{\rm real}^{\rm lin}
=
\Bigl(
\delta\chi_Q,\
M_*\,\Delta\mathbf x
\Bigr).
}
\]

So the fully nonlinear intrinsic packet and the exact projector calculus already collapse to the expected linear compiler without further assumptions.

---

## 9. What Stage 252 changes in the theorem problem

Stage 251 completed the reduced theorem logically.
Stage 252 makes it executable.

### 9.1 Packet B is now directly computable from target monomial values

No orbit representative and no pairwise witness are needed. The target orbit is specified by the three target monomial values, and the quotient packet is just their intrinsic log-ratio packet.

### 9.2 The unique microscopic repair is now explicit

The exact repair vector is carried only by
\[
(T_U,\ K_\eta^{(\mathrm{eff})},\ \mu_W),
\]
and it is given algebraically by the Stage 249 restoration law.

### 9.3 The realization audit has become fully procedural

The completed PDE now has to provide only:

- the Packet-A scalar `\(\chi_Q\)`,
- and the three target monomial ratios.

Everything else in the reduced realization audit is exact downstream algebra.

---

## 10. Immediate next derivation step

The next stage should use this compiler on the first actual reduced branch family returned by the moving-throat construction.

That means:

1. compute `\(\chi_Q\)` from the Packet-A grouped bundle,
2. compute the intrinsic Packet-B ratios
   \[
   \mathfrak R_{\rm tr},\ \mathfrak R_{\rm nt},\ \mathfrak R_\eta,
   \]
3. evaluate the exact four-scalar realization packet,
4. and, if the branch misses the target orbit, compute the exact canonical dependent-triple repair
   \[
   \Pi^{\rm can}_{\mathcal O_*}(\mathbf x)
   \]
   to locate the unique same-free-quintuple point on the target orbit.

That is the first truly explicit realization audit after Stage 252.
