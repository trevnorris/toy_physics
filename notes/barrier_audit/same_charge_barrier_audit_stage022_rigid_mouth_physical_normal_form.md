# Same-Charge Barrier Audit — Stage 022: Rigid-Mouth Physical Normal Form, Exact Physical-to-Microscopic Correction Compiler, and the Cartesian Orbit-Lock Theorem

## Status

**Exact within the carried Stage-019 rigid-mouth dependent-plane projector calculus and the Stage-021 physical-branch transfer-shape compiler.**

This stage does not solve the full moving-throat PDE branch.
It does the next clean reduction after Stage 021:

> it diagonalizes the rigid-mouth same-charge packet in the actual physical logarithmic variables 
> \((\mathcal T^2,\epsilon_\eta)\),
> converts that diagonal packet into an exact microscopic dependent-plane compiler,
> and shows that rigid-mouth orbit lock is already a Cartesian product of transfer-shape rigidity and dressing rigidity.

So after this stage the surviving rigid-mouth same-charge geometry is no longer triangular.
It is exactly a two-axis normal form:

1. a **pure transfer-shape axis**,
2. a **pure dressing axis**.

And the corresponding microscopic correction splits just as sharply:

- clearing the static transfer-shape defect changes only `\mu_W`,
- clearing the post-static dressing defect adds the equal `K_\eta^{(\mathrm{eff})}`–`\mu_W` shift.

---

## 0. Why this stage is needed

Stage 019 already proved that on the rigid-mouth slice the quotient-failure image in the dependent microscopic plane is
\[
\mathbf y_{\rm rm}
=
\mathbf y_{\rm nt}+\mathbf y_\eta,
\qquad
\mathbf y_{\rm nt}=\begin{pmatrix}0\\0\\q_{\rm nt}\end{pmatrix},
\qquad
\mathbf y_\eta=-q_\eta\begin{pmatrix}0\\1\\1\end{pmatrix}.
\]
So after the first static gate is cleared, the remaining obstruction is already the equal-drift `K_\eta`–`\mu` ray.

Stage 021 then rewrote the rigid-mouth finite packet directly in the physical branch observables and proved
\[
q_{\rm nt}=
\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
q_\eta=
\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
\]
So the remaining obvious question is:

> what happens if we use these physical logarithmic variables themselves as the rigid-mouth coordinates?

This stage shows that doing so diagonalizes the packet completely and turns the microscopic correction problem into an exact two-axis compiler.

---

## 1. Exact rigid-mouth physical logarithmic chart

On the rigid-mouth branch define the physical logarithmic coordinates
\[
\boxed{
U:=\ln\!\left(\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2}\right),
\qquad
V:=\ln\!\left(\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right).
}
\]
By Stage 021 these are exactly the surviving finite packet coordinates:
\[
\boxed{
q_{\rm nt}=U,
\qquad
q_\eta=V.
}
\]
So the rigid-mouth physical packet compiler is already diagonal:
\[
\boxed{
\mathbf q_{\rm rm}^{\rm phys}
:=
\begin{pmatrix}
q_{\rm nt}\\
q_\eta
\end{pmatrix}
=
M_{\rm phys}
\begin{pmatrix}
U\\
V
\end{pmatrix},
\qquad
M_{\rm phys}=I_2.
}
\]

The third direct observable on the rigid-mouth branch, `R_{\rm target}`, is then recovered exactly from the selected-branch identity
\[
R_{\rm target}\,\mathcal T^2=
\Lambda_0(1-\epsilon_\eta),
\qquad
R_{\rm target,ref}\,\mathcal T_{\rm ref}^2=
\Lambda_0(1-\epsilon_{\eta,\rm ref}),
\]
which gives
\[
\boxed{
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}\,e^{-U}.
}
\]
So the full rigid-mouth actual branch is already charted exactly by the diagonal logarithmic pair `(U,V)`.

---

## 2. Exact physical projectors and the two commuting finite legs

Because the rigid-mouth packet is diagonal in `(U,V)`, the canonical physical packet projectors are simply
\[
\boxed{
P_{\mathcal T}:=
\begin{pmatrix}1&0\\0&0\end{pmatrix},
\qquad
P_\eta:=
\begin{pmatrix}0&0\\0&1\end{pmatrix}.
}
\]
They are exact complementary projectors:
\[
P_{\mathcal T}^2=P_{\mathcal T},
\qquad
P_\eta^2=P_\eta,
\qquad
P_{\mathcal T}P_\eta=P_\eta P_{\mathcal T}=0,
\qquad
P_{\mathcal T}+P_\eta=I_2.
\]
So every rigid-mouth physical point decomposes uniquely as
\[
\boxed{
\begin{pmatrix}U\\V\end{pmatrix}
=
\begin{pmatrix}U\\0\end{pmatrix}
+
\begin{pmatrix}0\\V\end{pmatrix}.
}
\]

### 2.1 Pure transfer-shape leg

The image of `P_{\mathcal T}` is the exact finite transfer-shape leg
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2 e^U,
\qquad
\epsilon_\eta=\epsilon_{\eta,\rm ref},
\qquad
\frac{R_{\rm target}}{R_{\rm target,ref}}=e^{-U}.
}
\]
So pure transfer-shape motion changes only `\mathcal T^2` directly and compensates `R_{\rm target}` inversely.

### 2.2 Pure dressing leg

The image of `P_\eta` is the exact finite dressing leg
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2,
\qquad
\epsilon_\eta=\epsilon_{\eta,\rm ref}e^V,
\qquad
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}.
}
\]
So pure dressing motion is exactly the finite static-blind curve from Stage 020, now recognized as one coordinate axis of the physical chart.

### 2.3 Exact commutativity

Because the target ratio factorizes as
\[
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\underbrace{e^{-U}}_{\text{transfer leg}}
\underbrace{\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}}_{\text{dressing leg}},
\]
the two finite legs commute exactly.
So the rigid-mouth branch is an exact Cartesian product of

- transfer-shape motion, and
- dressing motion.

This is the physical normal-form version of the earlier packet projector calculus.

---

## 3. Exact physical-to-microscopic dependent-plane compiler

Stage 019 gives the dependent-plane packet carriers
\[
\mathbf y_{\rm nt}=\begin{pmatrix}0\\0\\q_{\rm nt}\end{pmatrix},
\qquad
\mathbf y_\eta=-q_\eta\begin{pmatrix}0\\1\\1\end{pmatrix}.
\]
Substituting the physical coordinates `q_{\rm nt}=U`, `q_\eta=V` gives the exact rigid-mouth dependent correction
\[
\boxed{
\mathbf y_{\rm rm}^{\rm dep}(U,V)
=
\begin{pmatrix}
\Delta_T\\
\Delta_{K_\eta}\\
\Delta_\mu
\end{pmatrix}
=
\begin{pmatrix}
0\\
-V\\
U-V
\end{pmatrix}.
}
\]
So the physical-to-microscopic compiler matrix is
\[
\boxed{
C_{\rm phys}^{\rm dep}
:=
\begin{pmatrix}
0&0\\
0&-1\\
1&-1
\end{pmatrix},
\qquad
\mathbf y_{\rm rm}^{\rm dep}=C_{\rm phys}^{\rm dep}
\begin{pmatrix}U\\V\end{pmatrix}.
}
\]
A left inverse is immediate:
\[
\boxed{
L_{\rm phys}^{\rm dep}
:=
\begin{pmatrix}
0&-1&1\\
0&-1&0
\end{pmatrix},
\qquad
L_{\rm phys}^{\rm dep}C_{\rm phys}^{\rm dep}=I_2.
}
\]
So the physical packet coordinates can be recovered directly from the dependent microscopic correction by
\[
\boxed{
U=\Delta_\mu-\Delta_{K_\eta},
\qquad
V=-\Delta_{K_\eta}.
}
\]

This is the cleanest rigid-mouth compiler obtained so far:

- `U` is the microscopic `\mu-K_\eta` difference,
- `V` is minus the `K_\eta` drift itself.

---

## 4. Exact microscopic images of the two physical axes

Push the physical projectors through `C_{\rm phys}^{\rm dep}`.

### 4.1 Pure transfer-shape image

Applying `P_{\mathcal T}` gives
\[
\boxed{
\mathbf y_{\mathcal T}^{\rm dep}
=
C_{\rm phys}^{\rm dep}
\begin{pmatrix}U\\0\end{pmatrix}
=
\begin{pmatrix}
0\\
0\\
U
\end{pmatrix}.
}
\]
So a pure transfer-shape defect is carried microscopically by a `\mu_W` shift only.

### 4.2 Pure dressing image

Applying `P_\eta` gives
\[
\boxed{
\mathbf y_\eta^{\rm dep}
=
C_{\rm phys}^{\rm dep}
\begin{pmatrix}0\\V\end{pmatrix}
=
-V\begin{pmatrix}0\\1\\1\end{pmatrix}.
}
\]
So a pure dressing defect is carried microscopically by the exact equal-drift `K_\eta^{(\mathrm{eff})}`–`\mu_W` ray.

This is the physical version of the Stage-019 dependent-plane theorem.
The difference is that the two microscopic pieces are now the direct images of the actual physical axes, not only of abstract quotient coordinates.

---

## 5. Exact correction compilers

Because the rigid-mouth packet is diagonal in `(U,V)`, the orbit-restoring corrections are immediate.

### 5.1 Static-only correction

To clear only the transfer-shape defect, subtract the pure transfer image:
\[
\boxed{
\Delta\mathbf y_{\rm static}
:=-\mathbf y_{\mathcal T}^{\rm dep}
=
\begin{pmatrix}0\\0\\-U\end{pmatrix}.
}
\]
After this correction,
\[
\mathbf y_{\rm rm}^{\rm dep}+
\Delta\mathbf y_{\rm static}
=
-V\begin{pmatrix}0\\1\\1\end{pmatrix},
\]
so the branch is moved exactly onto the pure dressing ray.

Thus the first static ceiling is cleared by changing only `\mu_W`.

### 5.2 Post-static dressing correction

Once the static ceiling has been cleared, the remaining orbit-restoring correction is just the opposite of the dressing image:
\[
\boxed{
\Delta\mathbf y_{\eta,\rm rest}
:=+V\begin{pmatrix}0\\1\\1\end{pmatrix}.
}
\]
So the extra step beyond the static gate is the equal positive shift in

- `K_\eta^{(\mathrm{eff})}`,
- `\mu_W`.

### 5.3 Full orbit-lock correction

Removing the whole rigid-mouth dependent correction gives
\[
\boxed{
\Delta\mathbf y_{\rm orbit}
:=-\mathbf y_{\rm rm}^{\rm dep}
=
\begin{pmatrix}
0\\
V\\
V-U
\end{pmatrix}
=
\Delta\mathbf y_{\rm static}+
\Delta\mathbf y_{\eta,\rm rest}.
}
\]
So the full orbit-restoring correction is literally the sum of

1. the static-only `\mu_W` correction,
2. the post-static equal `K_\eta`–`\mu_W` dressing correction.

This is the sharpest correction split reached so far.

---

## 6. Exact support-blindness of the physical normal form

Stage 021 already showed that the direct physical observables are support-blind:
\[
\partial_\zeta \mathcal T^2=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\qquad
\partial_{M_{\rm mix}}\mathcal T^2=0,
\qquad
\partial_{M_{\rm mix}}\epsilon_\eta=0.
\]
Therefore the physical logarithmic coordinates themselves satisfy
\[
\boxed{
\partial_\zeta U=
\partial_{M_{\rm mix}}U=
\partial_\zeta V=
\partial_{M_{\rm mix}}V=0.
}
\]
And because the microscopic compiler is linear in `(U,V)`, the dependent correction and all three correction compilers are support-blind as well:
\[
\boxed{
\partial_\zeta \mathbf y_{\rm rm}^{\rm dep}=
\partial_{M_{\rm mix}}\mathbf y_{\rm rm}^{\rm dep}=0,
}
\]
\[
\boxed{
\partial_\zeta \Delta\mathbf y_{\rm static}=
\partial_\zeta \Delta\mathbf y_{\rm orbit}=
0,
\qquad
\partial_{M_{\rm mix}} \Delta\mathbf y_{\rm static}=
\partial_{M_{\rm mix}} \Delta\mathbf y_{\rm orbit}=0.
}
\]
So coherent support enhancement cannot alter either

- the rigid-mouth physical packet, or
- the microscopic orbit-restoring correction required by that packet.

---

## 7. Cartesian orbit-lock theorem and first-order form

On the rigid-mouth actual branch,
\[
\boxed{
q_{\rm nt}=0,
\ q_\eta=0
\iff
U=0,
\ V=0
\iff
\mathcal T^2=\mathcal T_{\rm ref}^2,
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]
Because
\[
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}e^{-U},
\]
this is also equivalent to
\[
\boxed{
\mathcal T^2=\mathcal T_{\rm ref}^2,
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}
\iff
R_{\rm target}=R_{\rm target,ref},
\ \epsilon_\eta=\epsilon_{\eta,\rm ref}.
}
\]
So the rigid-mouth orbit-lock point is already a Cartesian codimension-two point in the physical logarithmic chart.

At first order the same statement becomes simply
\[
U=\delta\ln\mathcal T^2,
\qquad
V=\delta\ln\epsilon_\eta,
\]
so
\[
\boxed{
\mathbf y_{\rm rm}^{\rm dep,(1)}
=
\begin{pmatrix}
0\\
-\delta\ln\epsilon_\eta\\
\delta\ln\mathcal T^2-\delta\ln\epsilon_\eta
\end{pmatrix}.
}
\]
The static-blind line from Stages 018–020 is just the axis
\[
U=0,
\]
and its microscopic image is the pure equal-drift dressing ray.

So the earlier direct-space triangular packet was not the final normal form.
The later coherent transfer-shape compiler diagonalizes it completely.

---

## 8. Best current summary after Stage 022

The continuation from Stage 021 is now complete.

1. On the rigid-mouth actual branch, the physical logarithmic variables
   \[
   \left(\ln\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},\ \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right)
   \]
   are already the exact finite packet coordinates.
2. In that chart, the packet projector calculus is diagonal.
3. The microscopic dependent correction is exactly
   \[
   (\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V).
   \]
4. The static-only correction changes only `\mu_W`.
5. The post-static dressing correction is the equal `K_\eta`–`\mu_W` shift.
6. The whole rigid-mouth physical normal form and its microscopic correction compiler are support-blind.

So the same-charge barrier is now fully factorized on the rigid-mouth actual branch:

> first clear the pure transfer-shape axis `U`,
> then test whether the remaining pure dressing axis `V` collapses.

That is the sharpest Cartesian rigid-mouth statement reached so far.
