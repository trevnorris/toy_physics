# Moving-Throat PDE — Stage 175: Exact Orbit/Quotient Projectors and the Microscopic Orbit-Lock Split

## Status

**Exact within the carried Stage-170 coherent quotient closure and the finite-packet hierarchy isolated in Stage 174.**

This stage does **not** introduce any new constitutive law or new physical closure.
It upgrades the finite quotient packet of Stage 174 into an exact **microscopic projector calculus** on the full eight-dimensional drift space.

---

## Purpose

Stage 174 showed that the reduced home-stretch verdict depends only on two exact finite packets:

- the grouped bundle packet `\(\mathcal P_A\)`, and
- the quotient / orbit-lock packet
  \[
  \Delta_{\rm orbit}=(q_{\rm tr},q_{\rm nt},q_\eta).
  \]

That result is already strong enough for the endgame compiler. But it still leaves one important microscopic question open:

> given an arbitrary candidate microscopic drift of the actual moving-throat branch, how do we split it **exactly** into
> 
> 1. the pure similarity-orbit motion that preserves the branch invariants, and
> 2. the true quotient-failure piece measured by
>    \((q_{\rm tr},q_{\rm nt},q_\eta)\)?

That split is the missing bridge between the finite packet theorem and the actual microscopic branch-selection problem.

The main outputs of this stage are:

1. the exact pivot block in the dependent coordinates
   \[
   (\Delta_T,\Delta_{K_\eta},\Delta_\mu),
   \]
2. the exact canonical quotient section
   \[
   S_{(T,K_\eta,\mu)},
   \]
3. the exact complementary projectors
   \[
   Q_{\rm quot},\qquad O_{\rm orb}=I-Q_{\rm quot},
   \]
4. the sharp result that the entire quotient-failure piece has support only on the dependent triple
   \[
   (\Delta_T,\Delta_{K_\eta},\Delta_\mu),
   \]
5. and the exact theorem that the orbit-lock condition is simply
   \[
   Q_{\rm quot}\Delta\mathbf x=0
   \iff
   M_*\Delta\mathbf x=0.
   \]

So Stage 175 turns the abstract finite quotient packet of Stage 174 into a direct microscopic orbit/failure decomposition.

---

## 1. Carry-forward microscopic drift space and finite quotient packet

Work in the ordered finite log-ratio drift vector
\[
\boxed{
\Delta\mathbf x:=
\begin{pmatrix}
\Delta_\lambda\\
\Delta_c\\
\Delta_\gamma\\
\Delta_U\\
\Delta_{K_\eta}\\
\Delta_W\\
\Delta_\mu\\
\Delta_T
\end{pmatrix},
}
\]
where the entries are the logarithmic ratios of
\[
(\lambda_W,\ c_{\eta U},\ \gamma,\ K_U,\ K_\eta^{(\mathrm{eff})},\ K_W^{(\mathrm{eff})},\ \mu_W,\ T_U).
\]

Carry forward the exact monomial-drift map
\[
\boxed{
\mathbf q:=
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=M_*\,\Delta\mathbf x,
}
\]
with
\[
\boxed{
M_*=
\begin{pmatrix}
0 & 1+\delta_{U,*} & 1+\delta_{U,*} & -(2+\chi_{0,*}+\delta_{U,*}) & 0 & 0 & 0 & 1+\chi_{0,*}\\[4pt]
2(1+E_*) & 0 & 2E_* & F_*-E_* & -1 & -(2+E_*) & 1 & -F_*\\[4pt]
0 & 2 & 0 & -1 & -1 & 0 & 0 & 0
\end{pmatrix}.
}
\]

So the finite Packet B of Stage 174 is not abstract. It is literally the quotient projection of the microscopic drift under `\(M_*\)`.

The five free similarity coordinates are
\[
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W),
\]
while the three dependent microscopic coordinates are
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu).
\]

---

## 2. Exact pivot block on the dependent triple

Select the columns of `\(M_*\)` corresponding to the dependent triple in the order
\[
(T,K_\eta,\mu).
\]
Then
\[
\boxed{
P_{(T,K_\eta,\mu)}:=M_*[:,(T,K_\eta,\mu)]
=
\begin{pmatrix}
1+\chi_{0,*} & 0 & 0\\[4pt]
-F_* & -1 & 1\\[4pt]
0 & -1 & 0
\end{pmatrix}.
}
\]
Its determinant is
\[
\boxed{
\det P_{(T,K_\eta,\mu)}=1+\chi_{0,*}>0.
}
\]
So the dependent triple is an exact pivot block on the constructive coherent branch.

The inverse is therefore exact:
\[
\boxed{
P_{(T,K_\eta,\mu)}^{-1}
=
\begin{pmatrix}
\dfrac{1}{1+\chi_{0,*}} & 0 & 0\\[8pt]
0 & 0 & -1\\[4pt]
\dfrac{F_*}{1+\chi_{0,*}} & 1 & -1
\end{pmatrix}.
}
\]

This already shows that the quotient coordinates can be reinserted into microscopic drift space using only the dependent triple.

---

## 3. Exact canonical quotient section

Let
\[
E_{(T,K_\eta,\mu)}
\]
be the ambient `\(8\times 3\)` embedding that inserts a three-vector into the dependent positions
\[
(T,K_\eta,\mu)
\]
in the ordered drift basis `\((\lambda,c,\gamma,U,K_\eta,W,\mu,T)\)`.

Define the canonical quotient section
\[
\boxed{
S_{(T,K_\eta,\mu)}:=E_{(T,K_\eta,\mu)}\,P_{(T,K_\eta,\mu)}^{-1}.
}
\]
Explicitly,
\[
\boxed{
S_{(T,K_\eta,\mu)}=
\begin{pmatrix}
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & 0\\
0 & 0 & -1\\
0 & 0 & 0\\
\dfrac{F_*}{1+\chi_{0,*}} & 1 & -1\\[8pt]
\dfrac{1}{1+\chi_{0,*}} & 0 & 0
\end{pmatrix}.
}
\]
By construction,
\[
\boxed{M_*S_{(T,K_\eta,\mu)}=I_3.}
\]
So `\(S_{(T,K_\eta,\mu)}\)` is an exact right inverse of the quotient map.

This means every quotient packet
\[
\mathbf q=(q_{\rm tr},q_{\rm nt},q_\eta)^T
\]
already has a canonical microscopic representative supported only on the dependent triple:
\[
\boxed{\Delta\mathbf x_{\rm quot}:=S_{(T,K_\eta,\mu)}\,\mathbf q.}
\]

---

## 4. Exact complementary orbit/quotient projectors

Define
\[
\boxed{Q_{\rm quot}:=S_{(T,K_\eta,\mu)}M_*,}
qquad
\boxed{O_{\rm orb}:=I_8-Q_{\rm quot}.}
\]
Then the projector identities are exact:
\[
\boxed{Q_{\rm quot}^2=Q_{\rm quot},}
qquad
\boxed{O_{\rm orb}^2=O_{\rm orb},}
qquad
\boxed{Q_{\rm quot}O_{\rm orb}=O_{\rm orb}Q_{\rm quot}=0.}
\]
Moreover,
\[
\boxed{M_*O_{\rm orb}=0,}
qquad
\boxed{M_*Q_{\rm quot}=M_*.}
\]
So for every microscopic drift,
\[
\boxed{
\Delta\mathbf x
=
\underbrace{O_{\rm orb}\Delta\mathbf x}_{\Delta\mathbf x_{\rm orbit}}
+
\underbrace{Q_{\rm quot}\Delta\mathbf x}_{\Delta\mathbf x_{\rm fail}}.
}
\]
This decomposition is exact and unique, with
\[
\boxed{\Delta\mathbf x_{\rm orbit}\in\ker M_*,}
qquad
\boxed{M_*\Delta\mathbf x_{\rm fail}=M_*\Delta\mathbf x=\mathbf q.}
\]

So `\(O_{\rm orb}\)` is the exact microscopic orbit projector, while `\(Q_{\rm quot}\)` is the exact quotient-failure projector.

---

## 5. Exact support of the quotient-failure piece

Let
\[
\mathbf q=M_*\Delta\mathbf x=
\begin{pmatrix}
q_{\rm tr}\\
q_{\rm nt}\\
q_\eta
\end{pmatrix}.
\]
Then the quotient-failure piece is
\[
\boxed{\Delta\mathbf x_{\rm fail}=Q_{\rm quot}\Delta\mathbf x=S_{(T,K_\eta,\mu)}\mathbf q.}
\]
So its only nonzero components are
\[
\boxed{(\Delta_T)_{\rm fail}=\frac{q_{\rm tr}}{1+\chi_{0,*}},}
\]
\[
\boxed{(\Delta_{K_\eta})_{\rm fail}=-q_\eta,}
\]
\[
\boxed{(\Delta_\mu)_{\rm fail}=\frac{F_*}{1+\chi_{0,*}}q_{\rm tr}+q_{\rm nt}-q_\eta.}
\]
All five free similarity directions vanish in `\(\Delta\mathbf x_{\rm fail}\)`.

This is the sharpest microscopic result of the stage:

> **the entire failure of single-orbit lock is carried only by the dependent triple**
> \((\Delta_T,\Delta_{K_\eta},\Delta_\mu)\).

So the five free similarity coordinates never directly carry quotient failure.

---

## 6. Exact orbit piece and the single-orbit law

Because
\[
\Delta\mathbf x_{\rm orbit}=O_{\rm orb}\Delta\mathbf x,
\]
the orbit projector leaves the five free coordinates unchanged:
\[
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W)_{\rm orbit}
=
(\Delta_\lambda,\Delta_c,\Delta_\gamma,\Delta_U,\Delta_W).
\]
It rewrites only the dependent triple, forcing it onto the exact Stage-170 single-orbit law.

Define
\[
\boxed{\alpha_*:=\frac{1+\delta_{U,*}}{1+\chi_{0,*}}.}
\]
Then the orbit piece satisfies
\[
\boxed{(\Delta_{K_\eta})_{\rm orbit}=2\Delta_c-\Delta_U,}
\]
\[
\boxed{(\Delta_T)_{\rm orbit}=\Delta_U-\alpha_*(\Delta_\gamma+\Delta_c-\Delta_U),}
\]
\[
\boxed{
(\Delta_\mu)_{\rm orbit}
=
2\Delta_c-\Delta_U+2\Delta_W-2\Delta_\lambda
-
E_*(2\Delta_\gamma+2\Delta_\lambda-\Delta_U-\Delta_W)
-
F_*\alpha_*(\Delta_\gamma+\Delta_c-\Delta_U).
}
\]
These are exactly the Stage-170 finite orbit formulas recovered by projection.

So the projector calculus does not invent a new orbit law.
It isolates the exact same law as the unique orbit representative with the same five free coordinates as the original microscopic drift.

---

## 7. Exact orbit-lock theorem in projector form

The exact orbit-lock condition now has several equivalent forms:
\[
\boxed{\mathbf q=0,}
\qquad
\boxed{M_*\Delta\mathbf x=0,}
\qquad
\boxed{Q_{\rm quot}\Delta\mathbf x=0,}
\qquad
\boxed{\Delta\mathbf x=O_{\rm orb}\Delta\mathbf x,}
\qquad
\boxed{\Delta\mathbf x\in\ker M_*.}
\]
So the Packet-B zero-set of Stage 174 is now a literal microscopic projector statement.

In words:

- Stage 174 said the reduced orbit packet must vanish.
- Stage 175 shows that this is exactly the statement that the candidate microscopic drift has no quotient-supported dependent-triple correction.

That is the clean bridge from the finite-packet theorem back to the actual moving-throat PDE branch.

---

## 8. What this stage fixes, and what remains for Stage 176

### Fixed here

This stage fixes the exact microscopic geometry of the orbit packet.
In particular:

1. Packet B is exactly the quotient projection `\(M_*\Delta\mathbf x\)`.
2. The orbit/failure split is exact and complementary.
3. The quotient failure lives only on the dependent triple.
4. The orbit projector leaves the five similarity coordinates untouched and restores the exact Stage-170 single-orbit dependent law.

### Still missing

What this stage does **not** yet do is compose the microscopic projector with the observable defect packet
\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]
That is the next step.

Stage 176 will take the exact observable inversion from the earlier compiler stages and feed it through the canonical quotient section
\[
S_{(T,K_\eta,\mu)}
\]
to obtain the exact observable-to-microscopic correction compiler.

---

## 9. Script-backed status

The accompanying SymPy audit verifies:

- the exact Stage-170 monomial-drift matrix `\(M_*\)`,
- the dependent pivot block `\(P_{(T,K_\eta,\mu)}\)` and its inverse,
- the exact right-inverse section `\(S_{(T,K_\eta,\mu)}\)`,
- the complementary projector identities for `\(Q_{\rm quot}\)` and `\(O_{\rm orb}\)`,
- the support of the quotient-failure projector only on the dependent triple,
- the exact formulas for the failure components in terms of `\((q_{\rm tr},q_{\rm nt},q_\eta)\)`,
- the exact orbit-law formulas for the projected dependent coordinates,
- and the equivalence
  \[
  Q_{\rm quot}\Delta\mathbf x=0
  \iff
  M_*\Delta\mathbf x=0.
  \]

Supporting files:

- `moving_throat_pde_stage175_orbit_quotient_projectors_sympy_audit.py`
- `moving_throat_pde_stage175_orbit_quotient_projectors_sympy_audit_output.txt`
