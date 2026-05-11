# Moving-Throat PDE — Stage 193: Exact Isotropic Grouped-`P2` Target Surface and the Scalar/Geometry Firewall

## Status

**Exact within the carried Packet-A grouped-bundle hierarchy of Stage 242, the microscopic orbit/quotient split of Stage 243, and the reduced 5PN isotropic grouped-real `P2` closure.**

This stage does **not** introduce a new constitutive law or a new physical closure.
Its job is to freeze the first new audited theorem target that the 5PN derivation made unavoidable:

1. the exact **isotropic grouped-real `P2` conservative target surface**, and
2. the exact **scalar/geometry firewall** showing that the `l=0` geometry lane cannot contaminate that target surface at linear anisotropy order.

---

## Purpose

Stage 242 compressed the reduced endgame to the finite grouped bundle packet
\[
\mathcal P_A
=
\bigl((D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4})_{A=20,21,22},\ m_{\hat 0}\bigr),
\]
and the exact branch residual packet
\[
\Delta_{\rm branch}
=
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole},\Delta_{\rm norm}).
\]
Stage 243 then separated the microscopic quotient-failure packet from pure similarity-orbit motion, so the remaining PDE problem became a clean finite-packet realization problem.

But the 5PN compression sharpened something further:

> before worrying about outgoing normalization, orbit lock, or the realized coherent placement point, the completed moving-throat PDE must first land on one exact **conservative isotropic grouped-real `P2` surface**.

That surface is the real front-end theorem target. It is where the grouped `20/21/22` bundle stops being an arbitrary three-lane packet and becomes the minimal isotropic one-pole conservative carrier needed by the higher-order chain.

The main outputs of this stage are:

1. the exact conservative isotropy surface
   \[
   a_2=b_2=a_4=b_4=0,
   \]
2. the exact one-pole surface
   \[
   \Delta_{\rm pole}=\bar\nu_4-4\bar\nu_2^{\,2}=0,
   \]
3. the exact equivalent one-parameter conservative carrier
   \[
   \widehat Y_Q^{\rm cons}(\omega)
   =
   \frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
   \]
4. and the exact reduced theorem that the `l=0` scalar/geometry lane can re-enter the grouped `l=2` carrier only at **quadratic** order in anisotropy-induced mixing.

So this stage is the correct audited start of the new post-243 extension.

---

## 1. Carry-forward Packet-A conservative compiler

Work only with the conservative grouped-lane data from Packet A:
\[
D_A^{(\mathrm{cons})}(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+O(\omega^6),
\qquad
A\in\{20,21,22\}.
\]

As in Stage 242, define the normalized conservative grouped response moments
\[
\nu_2^{(A)}:=-\frac{D_{A2}}{D_{A0}},
\qquad
\nu_4^{(A)}:=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2}.
\]

For any grouped triple \((x_{20},x_{21},x_{22})\), the exact weighted trace/anomaly map is
\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
b_x=\frac{x_{21}-x_{22}}{2},
\]
with inverse
\[
x_{20}=\bar x+4a_x,
\qquad
x_{21}=\bar x-a_x+b_x,
\qquad
x_{22}=\bar x-a_x-b_x.
\]

Apply this to the conservative response moments:
\[
(\bar\nu_2,a_2,b_2),
\qquad
(\bar\nu_4,a_4,b_4).
\]

No outgoing data are needed yet. Stage 244 is deliberately only the conservative front end.

---

## 2. Exact isotropic grouped-real `P2` conservative target surface

The first theorem target is the exact isotropy surface
\[
\boxed{
\mathcal S_{\rm iso}^{(\mathrm{cons})}
:
\quad
a_2=b_2=a_4=b_4=0.
}
\]

Because the grouped inverse map is exact, this is equivalent to the common-lane collapse
\[
\nu_2^{(20)}=\nu_2^{(21)}=\nu_2^{(22)}=\bar\nu_2,
\]
\[
\nu_4^{(20)}=\nu_4^{(21)}=\nu_4^{(22)}=\bar\nu_4.
\]

Equivalently, at the level of the conservative operator moments,
\[
D_{20,0}=D_{21,0}=D_{22,0}=:D_0,
\]
\[
D_{20,2}=D_{21,2}=D_{22,2}=:D_2,
\]
\[
D_{20,4}=D_{21,4}=D_{22,4}=:D_4.
\]

So on the exact isotropic branch the three grouped lanes collapse to one common conservative carrier
\[
D_Q^{(\mathrm{cons})}(\omega)=D_0+D_2\omega^2+D_4\omega^4+O(\omega^6),
\]
and one common normalized conservative response
\[
\widehat Y_Q^{\rm cons}(\omega)
:=
\frac{D_0}{D_Q^{(\mathrm{cons})}(\omega)}
=
1+\bar\nu_2\omega^2+\bar\nu_4\omega^4+O(\omega^6).
\]

This is the exact grouped-real `P2` conservative front-end surface the completed PDE must hit before later outgoing or orbit-side tests even matter.

---

## 3. Exact one-pole conservative surface

The second theorem target is the exact one-pole conservative identity
\[
\boxed{
\mathcal S_{\rm pole}^{(\mathrm{cons})}
:
\quad
\Delta_{\rm pole}:=\bar\nu_4-4\bar\nu_2^{\,2}=0.
}
\]

On the isotropic branch,
\[
\bar\nu_2=-\frac{D_2}{D_0},
\qquad
\bar\nu_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]
so
\[
\Delta_{\rm pole}
=
\frac{D_2^2-D_0D_4}{D_0^2}-4\frac{D_2^2}{D_0^2}
=
-\frac{3D_2^2+D_0D_4}{D_0^2}.
\]

Therefore the one-pole surface is exactly equivalent to
\[
\boxed{
D_0D_4+3D_2^2=0.
}
\]
If one rewrites the conservative moments in the underlying isotropic bundle language,
\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4),
\]
then the same surface is
\[
\boxed{
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
}
\]

So the exact grouped conservative target surface is
\[
\boxed{
\mathcal S_{\rm iso+pole}^{(\mathrm{cons})}
:
\quad
a_2=b_2=a_4=b_4=0,
\qquad
\Delta_{\rm pole}=0.
}
\]

---

## 4. Exact one-parameter conservative carrier

Once the isotropic one-pole surface is imposed, the common grouped conservative carrier is forced into a one-parameter form.

Define
\[
\boxed{
\Omega_Q^2:=-\frac{D_0}{4D_2}.
}
\]
Because
\[
\bar\nu_2=-\frac{D_2}{D_0}=\frac{1}{4\Omega_Q^2},
\]
and the one-pole surface gives
\[
\bar\nu_4=4\bar\nu_2^{\,2}=\frac{1}{4\Omega_Q^4},
\]
the unique isotropic one-pole conservative precursor through `O(ω^4)` is
\[
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
}
\]
Indeed,
\[
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}
=
1+\frac{\omega^2}{4\Omega_Q^2}+\frac{\omega^4}{4\Omega_Q^4}+O(\omega^6)
=
1+\bar\nu_2\omega^2+\bar\nu_4\omega^4+O(\omega^6).
\]

So Stage 244 freezes the exact conservative theorem target in the cleanest possible form:

> the grouped `20/21/22` conservative carrier must be isotropic and one-pole, and therefore must reduce to the
> \(\frac34+\frac14(1-\omega^2/\Omega_Q^2)^{-1}\)
> module through the audited orders.

---

## 5. Exact scalar/geometry firewall

The next question is whether the `l=0` scalar/geometry lane can contaminate this grouped `l=2` conservative target at the same order.

On the exact isotropic quadratic wall theory, the answer is **no**.

### 5.1 Exact isotropic block separation

At the quadratic level, the reference throat is `O(3)` invariant, so the wall/support operator is block-diagonal in angular momentum:
\[
\mathcal D^{(0)}(\omega)
=
\mathcal D_{l=0}(\omega)\oplus \mathcal D_{l=2}(\omega)\oplus\cdots.
\]

In particular, there is no linear `l=0 \leftrightarrow l=2` coupling on the isotropic background. So the scalar/geometry lane does not enter the grouped `l=2` conservative carrier at all on the exact isotropic branch.

### 5.2 Weak anisotropy-induced mixing appears only quadratically

Now introduce a small anisotropy parameter `χ` that turns on the first `l=0 \leftrightarrow l=2` mixing. Write the reduced block operator schematically as
\[
\mathcal D(\omega,\chi)
=
\begin{pmatrix}
\mathcal D_{0}(\omega) & \chi\,C^T(\omega)\\[4pt]
\chi\,C(\omega) & \mathcal D_{2}(\omega)I_3
\end{pmatrix},
\]
where:

- `\(\mathcal D_0\)` is the scalar/geometry block,
- `\(\mathcal D_2 I_3\)` is the isotropic grouped `l=2` block,
- `\(C\)` is the anisotropy-induced mixing vector.

Eliminating the scalar/geometry block by an exact Schur complement gives
\[
\boxed{
\mathcal D_{2,\rm eff}(\omega,\chi)
=
\mathcal D_2(\omega)I_3
-
\chi^2\,C(\omega)\,\mathcal D_0(\omega)^{-1}C(\omega)^T.
}
\]

So the entire scalar/geometry contamination of the grouped `l=2` block carries an explicit factor `χ^2`.
There is **no** `O(χ)` contamination from the scalar/geometry lane.

This is the exact reduced firewall statement:
\[
\boxed{
\partial_\chi\mathcal D_{2,\rm eff}(\omega,\chi)\big|_{\chi=0}=0.
}
\]

Therefore:

1. the scalar/geometry lane does **not** contaminate the grouped `l=2` carrier at `O(χ^0)` because the isotropic theory is block-diagonal,
2. it does **not** contaminate it at `O(χ^1)` because the first Schur-complement correction is quadratic,
3. and any re-entry of the scalar/geometry lane into the grouped `l=2` conservative moments begins only at `O(χ^2)` through anisotropy-induced mixing.

This is precisely the firewall the 5PN derivation was implicitly using and which now becomes an explicit audited theorem target.

---

## 6. Exact Stage 244 theorem target

The conservative front-end theorem target after Stage 243 is now explicit.

### 6.1 Conservative grouped-real `P2` target surface

The completed moving-throat PDE must supply Packet-A data whose conservative part lands on
\[
\boxed{
\mathcal S_{244}
:
\quad
a_2=b_2=a_4=b_4=0,
\qquad
\Delta_{\rm pole}=0.
}
\]

Equivalently, the grouped conservative carrier must satisfy
\[
\boxed{
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}
+O(\omega^6)
}
\]
through the audited orders.

### 6.2 Firewall statement

Within the exact isotropic quadratic wall closure, the `l=0` scalar/geometry lane is not a linear contaminant of this target.
Any contamination of the grouped `l=2` conservative carrier from the scalar/geometry side begins only at quadratic order in anisotropy-induced mixing:
\[
\boxed{
\delta\mathcal D_{2,\rm geom}
=
O(\chi^2).
}
\]

So Stage 244 turns a previously implicit design rule into an explicit audited theorem surface.

---

## 7. Why this is the correct next audited stage after 243

Stage 242 said the home-stretch problem depends only on Packet A and Packet B.
Stage 243 said Packet B is an exact microscopic orbit/failure split.

The next audited extension therefore should not jump immediately to outgoing normalization or the final four-condition verdict.
It should first freeze the **conservative front-end surface** that Packet A must land on before those later stages even make sense.

That is exactly what Stage 244 does.

It does three useful things for every later stage:

1. it tells us what exact conservative grouped-real `P2` state the PDE must realize,
2. it gives one exact one-parameter carrier that later outgoing stages can attach to,
3. and it removes the temptation to blame linear grouped-lane failures on the scalar/geometry lane when the exact isotropic block structure forbids that.

So this is the right starting point for the new post-243 audited sequence.

---

## 8. Script-backed status

The accompanying SymPy audit verifies:

- the exact grouped trace/anomaly inverse map on the conservative side,
- vanishing of the grouped conservative anomalies on the common-lane isotropic branch,
- the exact one-pole identity
  \(\Delta_{\rm pole}=-(3D_2^2+D_0D_4)/D_0^2\),
- exact equivalence of the one-pole surface to
  \(D_0D_4+3D_2^2=0\),
- the exact `\(\frac34+\frac14(1-\omega^2/\Omega_Q^2)^{-1}\)` series on that surface,
- and the exact Schur-complement firewall showing that scalar/geometry contamination of the grouped `l=2` carrier starts only at `O(\chi^2)`.

Supporting file:
- `moving_throat_pde_stage244_isotropic_grouped_p2_target_surface_sympy_audit.py`
