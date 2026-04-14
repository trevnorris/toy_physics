# Same-Charge Barrier Audit — Stage 018: Rigid-Mouth Packet Projectors, the Static-Blind Dressing Line, and the Codimension-Two Orbit-Lock Point

## Status

**Exact within the carried Stage-017 direct branch-observable compiler and the later `5PN` orbit/quotient projector calculus.**

This stage does not add a new physical mechanism.
It upgrades the Stage-017 rigid-mouth strip picture into an exact **two-coordinate projector calculus** on the direct observable plane
\[
(R_1,E_1)
:=
(\delta\ln R_{\rm target},\,\delta\ln\epsilon_\eta).
\]

The main new result is sharper than the Stage-017 wording:

> on the rigid-mouth branch, the first static same-charge scalar sees only one quotient coordinate, `q_nt = Xi_1`, while a second exact dressing coordinate `q_eta` survives completely invisible to that static gate. So the static strip is a codimension-one test inside a codimension-two orbit-lock problem.

In other words, clearing the first static ceiling is necessary but not sufficient even after the mouth side is frozen.

---

## 0. Why this stage is needed

Stage 017 already proved that on the direct coherent branch
\[
\Theta_1=\delta\ln R_{\rm tr},
\qquad
\Xi_1=-\delta\ln R_{\rm target}-c_\eta\,\delta\ln\epsilon_\eta,
\qquad
\mathcal R_1=\delta\ln R_{\rm target},
\]
with
\[
c_\eta:=\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}.
\]
So under rigid-mouth lock
\[
\delta\ln R_{\rm tr}=0,
\]
the surviving static scalar is already
\[
\Xi_1=-R_1-c_\eta E_1.
\]

But the later orbit/quotient program also shows that full orbit lock is an exact quotient-packet statement: the quotient packet must vanish, not just one scalar projection of it. In the full microscopic language, orbit lock is exactly the zero-set condition
\[
Q_{\rm quot}\,\Delta x=0
\iff
M_*\Delta x=0,
\]
where the finite Packet-B coordinates are `q = (q_tr,q_nt,q_eta)`. So the right next step is to reduce that projector picture onto the rigid-mouth slice and see what the static strip is actually testing. This is exactly what Stage 175 isolates at the microscopic level. 

---

## 1. Exact rigid-mouth direct packet

On the rigid-mouth branch define the direct observable vector
\[
\boxed{
\mathbf x_{\rm rm}:=
\begin{pmatrix}
R_1\\
E_1
\end{pmatrix}
=
\begin{pmatrix}
\delta\ln R_{\rm target}\\
\delta\ln\epsilon_\eta
\end{pmatrix}.
}
\]
The surviving quotient packet is
\[
\boxed{
\mathbf q_{\rm rm}:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
\begin{pmatrix}
\Xi_1\\
E_1
\end{pmatrix}.
}
\]
Using the Stage-017 direct compiler,
\[
q_{\rm nt}=-R_1-c_\eta E_1,
\qquad
q_\eta=E_1,
\]
so
\[
\boxed{
\mathbf q_{\rm rm}=M_{\rm rm}\,\mathbf x_{\rm rm},
\qquad
M_{\rm rm}:=
\begin{pmatrix}
-1 & -c_\eta\\
0 & 1
\end{pmatrix}.
}
\]
This matrix is an involution:
\[
\boxed{M_{\rm rm}^2=I_2.}
\]
So the direct observable plane and the rigid-mouth quotient packet are related by an exact one-step involutive compiler.

The inverse map is therefore identical:
\[
\boxed{
\mathbf x_{\rm rm}=M_{\rm rm}\,\mathbf q_{\rm rm},
}
\]
that is,
\[
\boxed{R_1=-q_{\rm nt}-c_\eta q_\eta,}
qquad
\boxed{E_1=q_\eta.}
\]

---

## 2. Exact canonical packet projectors on the rigid-mouth plane

In packet space, the two obvious complementary projectors are
\[
Q_{\rm nt}:=
\begin{pmatrix}
1&0\\
0&0
\end{pmatrix},
\qquad
Q_\eta:=
\begin{pmatrix}
0&0\\
0&1
\end{pmatrix},
\qquad
Q_{\rm nt}+Q_\eta=I_2.
\]
Push them back to direct-observable space with the exact section `S_rm = M_rm`:
\[
\boxed{P_{\rm nt}:=S_{\rm rm}Q_{\rm nt}M_{\rm rm}},
\qquad
\boxed{P_\eta:=S_{\rm rm}Q_\eta M_{\rm rm}}.
\]
Because `M_rm` is its own inverse, these are explicit:
\[
\boxed{
P_{\rm nt}=
\begin{pmatrix}
1 & c_\eta\\
0 & 0
\end{pmatrix},
\qquad
P_\eta=
\begin{pmatrix}
0 & -c_\eta\\
0 & 1
\end{pmatrix}.
}
\]
They are exact complementary projectors:
\[
P_{\rm nt}^2=P_{\rm nt},
\qquad
P_\eta^2=P_\eta,
\qquad
P_{\rm nt}P_\eta=P_\eta P_{\rm nt}=0,
\qquad
P_{\rm nt}+P_\eta=I_2.
\]

So every rigid-mouth direct branch point splits uniquely as
\[
\boxed{
\mathbf x_{\rm rm}=\mathbf x_{\rm nt}+\mathbf x_\eta,
\qquad
\mathbf x_{\rm nt}:=P_{\rm nt}\mathbf x_{\rm rm},
\qquad
\mathbf x_\eta:=P_\eta\mathbf x_{\rm rm}.
}
\]
Explicitly,
\[
\boxed{
\mathbf x_{\rm nt}=
\begin{pmatrix}
R_1+c_\eta E_1\\
0
\end{pmatrix}
=
\begin{pmatrix}
-\Xi_1\\
0
\end{pmatrix},
}
\]
\[
\boxed{
\mathbf x_\eta=
\begin{pmatrix}
-c_\eta E_1\\
E_1
\end{pmatrix}
=
\begin{pmatrix}
-c_\eta q_\eta\\
q_\eta
\end{pmatrix}.
}
\]
So the rigid-mouth direct plane already contains two exact, complementary, physically meaningful pieces:

- `x_nt`: the piece seen by the first static scalar,
- `x_eta`: the pure dressing/selected-target drift that the static scalar cannot see.

---

## 3. Exact codimension-two orbit-lock theorem on the rigid-mouth branch

Because the rigid-mouth quotient packet is
\[
\mathbf q_{\rm rm}=(q_{\rm nt},q_\eta)^T=(\Xi_1,E_1)^T,
\]
full rigid-mouth orbit lock is exactly
\[
\boxed{
\mathbf q_{\rm rm}=0
\iff
q_{\rm nt}=0\ \text{and}\ q_\eta=0.
}
\]
Equivalently,
\[
\boxed{
R_1=0,
\qquad
E_1=0.
}
\]
Using the direct compiler again,
\[
\boxed{
\mathbf q_{\rm rm}=0
\iff
\Xi_1=0\ \text{and}\ R_1=0.
}
\]
Since `q_eta = E_1`, the condition `Xi_1 = 0` alone is not enough.

So on the rigid-mouth branch:

- the static strip `Xi_1 = 0` is codimension one,
- the true orbit-lock point is codimension two.

This is the sharpest version yet of “the static gate is not the whole orbit-lock problem.”

---

## 4. The static-blind dressing line

The compensated strip from Stage 017 is now simply the direct-space image of the packet projector `Q_eta`:
\[
\Xi_1=0
\iff
q_{\rm nt}=0
\iff
\mathbf x_{\rm rm}=\mathbf x_\eta=
\begin{pmatrix}
-c_\eta q_\eta\\
q_\eta
\end{pmatrix}.
\]
So the entire static-blind line is parameterized by the single dressing coordinate `q_eta`.

Its direct-space norm is exact:
\[
\boxed{
\|\mathbf x_\eta\|^2=(1+c_\eta^2)\,q_\eta^2.
}
\]
Therefore the static strip contains points arbitrarily far from the orbit-lock point.
For any prescribed size `L > 0`, choose
\[
q_\eta = \frac{L}{\sqrt{1+c_\eta^2}},
\]
then
\[
\Xi_1=0,
\qquad
\|\mathbf x_{\rm rm}\|=L.
\]

So the first static same-charge gate does **not** bound the true rigid-mouth orbit-failure amplitude.
It only bounds the `q_nt` component.

This is the exact static-blindness theorem on the rigid-mouth slice.

---

## 5. Canonical correction compilers

The projector formulas immediately give the two natural rigid-mouth corrections.

### 5.1 Static-only correction

Project to the static strip by removing only the `q_nt` component:
\[
\boxed{
\Delta\mathbf x_{\rm static}:=-\mathbf x_{\rm nt}
=
\begin{pmatrix}
\Xi_1\\
0
\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf x_{\rm rm}+\Delta\mathbf x_{\rm static}=\mathbf x_\eta,
\qquad
\Xi_1\to 0,
\qquad
q_\eta\to q_\eta.
\]
So the branch clears the first static ceiling but generally still fails orbit lock.

### 5.2 Full orbit-lock correction

Project all the way to the orbit-lock point by removing both packet components:
\[
\boxed{
\Delta\mathbf x_{\rm orbit}:=-\mathbf x_{\rm rm}
=
\begin{pmatrix}
q_{\rm nt}+c_\eta q_\eta\\
-q_\eta
\end{pmatrix}.
}
\]
This is exactly the sum of the static correction and the static-blind dressing correction:
\[
\boxed{
\Delta\mathbf x_{\rm orbit}
=
\Delta\mathbf x_{\rm static}+
\begin{pmatrix}
 c_\eta q_\eta\\
-q_\eta
\end{pmatrix}.
}
\]
So the extra step beyond the static gate is not mysterious. It is the exact removal of the packet component `q_eta`.

---

## 6. What changes physically after Stage 018

Stage 017 already showed that first-order same-charge survival on the rigid-mouth branch is governed by the strip
\[
|\epsilon\Xi_1|\le B,
\qquad
\Xi_1=-R_1-c_\eta E_1.
\]
Stage 018 now sharpens that into an exact orbit-packet statement:

1. the static gate constrains only `q_nt = Xi_1`,
2. the dressing coordinate `q_eta = E_1` survives completely outside that gate,
3. and the full rigid-mouth orbit-lock problem is exactly the vanishing of both packet coordinates.

So the next honest theorem gate is no longer “does the branch clear the static strip?”
It is:

> compute the actual rigid-mouth dressing coordinate `q_eta = \delta\ln\epsilon_\eta` (equivalently `R_1` once `Xi_1` is known), because that is the exact static-blind residue that still blocks orbit lock after the first static same-charge ceiling is cleared.
