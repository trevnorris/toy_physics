# Stage V2-09 — Maxwell/mixed-sector kernel audit

## Status

**Verdict:** conditional pass inside the declared one-lane reduced Maxwell/mixed closure.

The audit checks the smallest reduced block that can carry an honest passive/outgoing channel:

\[
(Q,\ U,\ W),
\]

where \(Q\) is a wall/worldtube amplitude, \(U\) is a brane-like localized Maxwell coordinate, and \(W\) is the mixed-sector coordinate representing the \(A_w/F_{\mu w}/J^w\)-active block.

This stage does **not** prove that the full nonlinear moving-throat PDE realizes the needed branch. It proves that the reduced Maxwell/mixed mechanism has the right algebraic structure, stability gates, outgoing-transfer sign, and scalar-rescue condition.

---

## 1. Mixed-sector gauge invariants

Use the gauge convention

\[
A_0\mapsto A_0-\partial_t\chi,\qquad
A_a\mapsto A_a+\partial_a\chi,\qquad
A_w\mapsto A_w+\partial_w\chi.
\]

The mixed fields are

\[
E_w=-\partial_t A_w-\partial_w A_0,
\]

\[
C_a=\partial_a A_w-\partial_w A_a.
\]

The script verifies exactly that

\[
\delta_\chi E_w=0,
\qquad
\delta_\chi C_a=0.
\]

So the mixed variables used in this stage are gauge-invariant reduced observables, not gauge artifacts.

---

## 2. Reduced one-lane Maxwell/mixed model

The reduced quadratic Lagrangian is

\[
L
=
\frac12M\dot Q^2-\frac12KQ^2
+
\frac12\dot U^2-\frac12\Omega_U^2U^2
+
\frac12\dot W^2-\frac12\Omega_W^2W^2
+
R\,UW
+
g_UQU
+
g_WQW.
\]

Here:

- \(Q\) is the wall/worldtube amplitude,
- \(U\) is the localized brane-like Maxwell coordinate,
- \(W\) is the mixed \(A_w/F_{\mu w}/J^w\)-active coordinate,
- \(R\) mixes the two internal gauge-sector coordinates,
- \(g_U,g_W\) couple the wall to the internal gauge block.

With the convention \(e^{-i\omega t}\), define

\[
A(\omega)=\Omega_U^2-\omega^2,
\qquad
B(\omega)=\Omega_W^2-\omega^2,
\]

\[
\Delta(\omega)=A(\omega)B(\omega)-R^2.
\]

The internal equations are

\[
\begin{pmatrix}
A(\omega)&-R\\
-R&B(\omega)
\end{pmatrix}
\begin{pmatrix}
U\\W
\end{pmatrix}
=
\begin{pmatrix}
g_UQ\\g_WQ
\end{pmatrix}.
\]

The inverse internal block is

\[
\frac{1}{\Delta(\omega)}
\begin{pmatrix}
B(\omega)&R\\
R&A(\omega)
\end{pmatrix}.
\]

Therefore the conservative self-energy inherited by the wall is

\[
\boxed{
\Sigma_{\rm cons}(\omega)
=
\frac{
g_U^2B(\omega)+2g_Ug_WR+g_W^2A(\omega)
}{
\Delta(\omega)
}.
}
\]

The wall operator is

\[
\boxed{
D_{\rm cons}(\omega)=K-M\omega^2-\Sigma_{\rm cons}(\omega).
}
\]

The SymPy audit verifies the Schur-complement numerator and denominator exactly.

---

## 3. Conservative low-frequency expansion

Let

\[
\Delta_0=\Omega_U^2\Omega_W^2-R^2,
\]

\[
S_2=\Omega_U^2+\Omega_W^2,
\]

\[
Q_0=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
G_2=g_U^2+g_W^2.
\]

Then

\[
\Sigma_{\rm cons}(\omega)
=
z_0+z_2\omega^2+z_4\omega^4+O(\omega^6),
\]

with

\[
\boxed{
z_0=\frac{Q_0}{\Delta_0},
}
\]

\[
\boxed{
z_2=\frac{Q_0S_2-G_2\Delta_0}{\Delta_0^2},
}
\]

\[
\boxed{
z_4=
\frac{
Q_0(S_2^2-\Delta_0)-S_2G_2\Delta_0
}{
\Delta_0^3
}.
}
\]

The script verifies this by multiplying the truncated series by the denominator and checking equality through \(O(\omega^4)\).

---

## 4. Stability and positivity gates

The static internal Hessian is

\[
H_{\rm int}
=
\begin{pmatrix}
\Omega_U^2&-R\\
-R&\Omega_W^2
\end{pmatrix}.
\]

The internal block is positive if

\[
\boxed{
\Omega_U^2>0,\qquad
\Omega_W^2>0,\qquad
\Delta_0>0.
}
\]

The full static Hessian for \((Q,U,W)\) is

\[
H_{\rm full}
=
\begin{pmatrix}
K&-g_U&-g_W\\
-g_U&\Omega_U^2&-R\\
-g_W&-R&\Omega_W^2
\end{pmatrix}.
\]

The script verifies

\[
\det H_{\rm full}
=
K\Delta_0-Q_0
=
\Delta_0\left(K-\frac{Q_0}{\Delta_0}\right).
\]

Thus the full reduced one-lane stability gate is

\[
\boxed{
\Delta_0>0,
\qquad
K-\Sigma_{\rm cons}(0)>0.
}
\]

Equivalently,

\[
\boxed{
K\Delta_0-Q_0>0.
}
\]

This is the Maxwell/mixed analogue of the BdG softening gate: conservative support can soften the wall, but the branch must stay on the positive side of the Schur boundary.

---

## 5. Outgoing-port dressing

Attach a passive outgoing channel to the mixed coordinate by replacing

\[
B(\omega)\mapsto B(\omega)-\Pi_{\rm out}(\omega).
\]

Then

\[
\Sigma_{\rm full}(\omega)
=
\frac{
g_U^2[B(\omega)-\Pi_{\rm out}]+2g_Ug_WR+g_W^2A(\omega)
}{
A(\omega)[B(\omega)-\Pi_{\rm out}]-R^2
}.
\]

Expanding to first order in \(\Pi_{\rm out}\),

\[
\Sigma_{\rm full}(\omega)
=
\Sigma_{\rm cons}(\omega)
+
\Pi_{\rm out}(\omega)\,N(\omega)
+
O(\Pi_{\rm out}^2),
\]

where the script verifies the exact transfer factor

\[
\boxed{
N(\omega)
=
\frac{[A(\omega)g_W+Rg_U]^2}{\Delta(\omega)^2}.
}
\]

At zero frequency,

\[
\boxed{
N_0
=
\frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta_0^2}\ge0.
}
\]

So the reduced Maxwell/mixed block transfers the outgoing branch with a nonnegative static weight.

There is one important caveat:

\[
\Omega_U^2g_W+Rg_U=0
\]

is a **dark-port condition**. If this happens, then

\[
N_0=0,
\]

and the leading quadrupole outgoing channel is not transferred to the wall.

---

## 6. Transfer-factor low-frequency expansion

Let

\[
P=\Omega_U^2g_W+Rg_U.
\]

Then

\[
N(\omega)=N_0+N_2\omega^2+N_4\omega^4+O(\omega^6),
\]

with

\[
\boxed{
N_0=\frac{P^2}{\Delta_0^2},
}
\]

\[
\boxed{
N_2=
\frac{2P(PS_2-\Delta_0g_W)}{\Delta_0^3},
}
\]

\[
\boxed{
N_4=
\frac{
\Delta_0^2g_W^2
-2\Delta_0P^2
-4\Delta_0PS_2g_W
+3P^2S_2^2
}{
\Delta_0^4
}.
}
\]

These are the one-lane transfer moments that later grouped-\(P_2\) stages must generalize.

---

## 7. Outgoing \(l=2\) wall coefficient

The compact outgoing \(l=2\) branch has the normalized odd coefficient

\[
\Gamma_2^{\rm port}
=
\frac{a^5}{27c_s^5}.
\]

So

\[
\Pi_2^{\rm out}(\omega)
=
+i\frac{a^5}{27c_s^5}\omega^5+\cdots.
\]

Because

\[
D(\omega)=K-M\omega^2-\Sigma(\omega),
\]

the wall operator inherits

\[
\boxed{
\delta D_2^{\rm odd}(\omega)
=
-i\,N_0\,\frac{a^5}{27c_s^5}\omega^5
+
O(\omega^7).
}
\]

The negative imaginary sign is the passive damping sign in the \(e^{-i\omega t}\) wall-operator convention. Written in normalized response/admittance language, this corresponds to the positive outgoing sign used in the 2.5PN package.

---

## 8. Scalar derivative-coupling compatibility

A naive scalar outgoing port can generate an \(i\omega\) term. The reduced Maxwell/mixed block avoids that only if the scalar outlet is derivative-coupled.

Take

\[
g_U=0,\qquad
g_W(\omega)=\eta\omega.
\]

Then

\[
N_0^{\rm scalar}(\omega)
=
\frac{A(\omega)^2\eta^2\omega^2}{\Delta(\omega)^2}
=
\frac{\Omega_U^4\eta^2}{\Delta_0^2}\omega^2+O(\omega^4).
\]

If

\[
\Pi_0^{\rm out}(\omega)=i\gamma_1\omega+\cdots,
\]

then the wall-level correction is

\[
\delta D_0^{\rm odd}
=
-i\gamma_1
\frac{\Omega_U^4\eta^2}{\Delta_0^2}
\omega^3
+
O(\omega^5).
\]

So the scalar outlet is demoted from \(i\omega\) to \(i\omega^3\).

This is not a full scalar no-go theorem. It is a precise reduced compatibility condition:

\[
\boxed{
\text{direct non-derivative scalar port coupling is dangerous;}
}
\]

\[
\boxed{
\text{derivative scalar port coupling delays the leading scalar odd term.}
}
\]

---

## 9. Pass/fail ledger

| Gate | Result |
|---|---:|
| Mixed \(E_w,C_a\) gauge invariance | Pass |
| Conservative Schur self-energy | Pass |
| \(z_0,z_2,z_4\) low-frequency expansion | Pass |
| Outgoing transfer factor | Pass |
| \(N_0,N_2,N_4\) transfer expansion | Pass |
| Static stability determinant | Pass |
| \(l=2\) odd wall coefficient | Pass |
| Scalar derivative-coupling demotion | Pass |
| Full PDE branch realization | Open |

---

## 10. Carry-forward theorem gates

This stage leaves three concrete conditions for the full moving-throat branch:

1. **Stability gate**

   \[
   \Delta_0>0,
   \qquad
   K-\frac{Q_0}{\Delta_0}>0.
   \]

2. **Quadrupole transfer gate**

   \[
   N_0=
   \frac{(\Omega_U^2g_W+Rg_U)^2}{\Delta_0^2}
   \neq 0.
   \]

3. **Universal normalization gate**

   After grouped-\(P_2\) lifting,

   \[
   \widehat m_0^2\frac{N_0}{D_0}
   =
   \frac{54Gc_s^5}{5a^5c^5}.
   \]

This stage proves that the reduced Maxwell/mixed mechanism is algebraically capable of carrying the passive/outgoing \(l=2\) branch. It does not prove that the actual moving-throat PDE lands on the target value.
