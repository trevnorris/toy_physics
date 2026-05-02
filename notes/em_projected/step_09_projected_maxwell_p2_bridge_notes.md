# Projected Maxwell → Grouped `P2` Bridge Notes

## What this extension does

This note continues the projection-first Maxwell program one step deeper into the
moving-throat PDE language.

The earlier projected-Maxwell work already established three facts:

1. the exact projected inhomogeneous Maxwell law keeps the mixed
   \(\partial_w(ZF^{w\nu})\) channel alive;
2. a **symmetric interior** projection first differs from the far-field brane law
   only at \(O(\sigma^2)\);
3. a **mouth-anchored** projection first differs already at \(O(\ell)\).

The new step is to connect those local projected-Maxwell terms directly to the
coefficients that the grouped real `P2` / full-bundle moving-throat program
actually uses:
\[
D_0,\quad D_2,\quad D_4,\quad N_0,\quad N_2,\quad N_4,
\]
and then to the derived response / prefactor objects
\[
u_2,\quad u_4,\quad P_0,\quad P_2,\quad P_4,\quad \Xi_1=\frac{P_1}{P_0}.
\]

So this is the first explicit bridge from:

- **projection-first near-throat electromagnetism**, to
- **the grouped `P_2` normalization / anisotropy language** used by the PDE and
  2.5PN/4PN/5PN notes.

---

## 1. The load-bearing project structures

The exact parent theory defines brane observables by a projection kernel \(W(w)\),
while the electromagnetic dynamics are localized by \(Z(w)\). The mixed channels
\(A_w\), \(F_{\mu w}\), and \(J^w\) remain part of the microscopic ontology outside
the strict far-field brane reduction. Projection therefore yields an exact open-system
brane electrodynamics rather than a closed copy of ordinary \(3+1\) Maxwell. fileciteturn2file0turn2file15turn2file18

The moving-throat lift then replaces the old \((a,L)\) closure by a distributed throat
shape \(R(\Omega,w,t)\), with grouped real `P2` wall/support modes appearing as literal
\(l=2\) lanes of the same geometry field. That is exactly why the grouped bundle is the
right target for the projected-Maxwell extension. fileciteturn0file17turn0file16

At the reduced-bundle level, the isotropic branch is organized by
\[
D_0=K-B_0-Z_0,\qquad
D_2=-(M+B_2+Z_2),\qquad
D_4=-(B_4+Z_4),
\]
\[
u_2=-\frac{D_2}{D_0},\qquad
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]
\[
P_0=\frac{N_0}{D_0},\qquad
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},\qquad
P_4=\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]
The exact isotropic target surface is
\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2,
\]
together with
\[
mhat_0^2\frac{N_0}{D_0}=\frac{54Gc_s^5}{5a^5c^5},
\]
and, on the constant-prefactor branch,
\[
N_2=\frac{2D_2N_0}{D_0},\qquad
N_4=\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]
These are the precise bundle targets the new projected-Maxwell piece must feed. fileciteturn3file2turn3file4turn3file15turn5file19

---

## 2. How the near-throat projection enters the bundle

The projected-Maxwell near-throat law shows that the mixed-sector throat term
\[
\langle \partial_w(ZF^{w\nu})\rangle
\]
survives exactly close to the mouth, and that mouth-anchored observation kernels
produce first corrections at \(O(\ell)\), not \(O(\sigma^2)\). fileciteturn2file0turn0file17

At bundle level, the cleanest way to represent this is:

\[
Z_n \longrightarrow Z_n + \varepsilon\,z_n,
\qquad
N_n \longrightarrow N_n + \varepsilon\,n_n,
\]
with \(\varepsilon\sim \ell\) for a mouth-anchored projection and
\(\varepsilon\sim \sigma^2\) for a symmetric interior slice.

Because only the Maxwell/mixed block is being corrected here, the conservative bundle shifts are
\[
D_0\to D_0-\varepsilon z_0,\qquad
D_2\to D_2-\varepsilon z_2,\qquad
D_4\to D_4-\varepsilon z_4,
\]
while the outgoing-transfer moments shift as
\[
N_0\to N_0+\varepsilon n_0,\qquad
N_2\to N_2+\varepsilon n_2,\qquad
N_4\to N_4+\varepsilon n_4.
\]

The SymPy derivation shows the induced first-order changes are

\[
\delta u_2=\frac{D_0 z_2-D_2 z_0}{D_0^2},
\]
\[
\delta u_4=
\frac{D_0^2 z_4-D_0(2D_2z_2+D_4z_0)+2D_2^2 z_0}{D_0^3},
\]
\[
\delta P_0=\frac{D_0 n_0+N_0 z_0}{D_0^2}.
\]

So the first useful conclusion is:

> **projection-first Maxwell feeds the moving-throat bundle in exactly the same slots
> the 5PN notes identify as load-bearing: \(Z_2,Z_4,N_0,N_2,N_4\), plus \(Z_0\) through
> denominator transport.**

That is already a good sign that this is the right missing ingredient to examine.

---

## 3. What changes the isotropic target surface — and what surprisingly cancels

The exact isotropic one-pole defect
\[
\mathcal P:=D_0(B_4+Z_4)-3(M+B_2+Z_2)^2
\]
shifts by
\[
\delta \mathcal P = D_0 z_4 - (B_4+Z_4) z_0 - 6(M+B_2+Z_2)z_2.
\]

So the first projected-Maxwell mouth corrections matter here through:

- \(z_2\) and \(z_4\) directly,
- and \(z_0\) through the static denominator.

More interesting is the **compatibility surface** from the isotropic full-bundle notes:
\[
\frac{N_0}{P_{0,\mathrm{target}}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4}.
\]
Its first projected-Maxwell shift is
\[
\delta \mathcal C
=
\frac{n_0}{P_{0,\mathrm{target}}}
-\frac{6S z_2}{T}
+\frac{3S^2 z_4}{T^2},
\qquad
S:=M+B_2+Z_2,\quad T:=B_4+Z_4.
\]

The important structural surprise is:

\[
\boxed{z_0\ \text{drops out of}\ \delta \mathcal C.}
\]

The updated SymPy proof now derives this in two exact forms and checks that
they agree:

1. solve the perturbed one-pole and normalization equations for the shifted
   wall stiffness \(K\) and subtract those two surfaces;
2. eliminate \(K\) and \(D_0\) directly to get the transported compatibility
   surface
   \[
   \mathcal C(\varepsilon)=
   \frac{N_0+\varepsilon n_0}{P_{0,\mathrm{target}}}
   -
   \frac{3(S+\varepsilon z_2)^2}{T+\varepsilon z_4}.
   \]

Because the eliminated surface itself has no \(z_0\) slot, the loss of
\(z_0\) is now a genuine algebraic consequence of the compatibility reduction.

The script now checks one stronger variant as well: impose the ratio-form
transported target
\[
P_{0,\mathrm{target}}(\varepsilon)
=
\frac{N_0+\varepsilon n_0}{D_{0,\mathrm{target}}},
\]
with an explicit transported denominator slot \(D_{0,\mathrm{target}}\). Then
the normalization solve gives
\[
K=B_0+Z_0+\varepsilon z_0 + D_{0,\mathrm{target}},
\]
so the induced compatibility transport drops not only \(z_0\) but also the
\(n_0\) numerator slot, leaving only the \(z_2,z_4\) geometry transport from
the one-pole side.

So:

- a projected-Maxwell near-throat shift of \(Z_0\) **does** change \(P_0=N_0/D_0\) directly,
- but it **does not** move the Stage-18 isotropic compatibility equation between the
  one-pole and normalization surfaces.

That is a genuinely useful sorting result: if the isotropic branch is failing the
compatibility test, the first near-throat projected-Maxwell quantities to inspect are
\[
z_2,\quad z_4,\quad n_0,
\]
not \(z_0\).

---

## 4. Constant-prefactor branch transport

On the exact constant-prefactor branch,
\[
P_2=0,\qquad P_4=0,
\]
the base conditions are
\[
N_2=\frac{2D_2N_0}{D_0},\qquad
N_4=\frac{2D_0(D_2N_2+D_4N_0)-3D_2^2N_0}{D_0^2}.
\]

The projected-Maxwell corrections then shift these conditions by explicit linear
combinations of
\[
z_0,\ z_2,\ z_4,\ n_0,\ n_2,\ n_4.
\]

So the near-throat projected-Maxwell data do not merely change the normalization
\(P_0\). They also push the system on or off the constant-prefactor branch. That matters,
because the 2.5PN/4PN bridge repeatedly singles out precisely that branch as the natural
passive/outgoing one. fileciteturn3file4turn5file19

---

## 5. Weak-axisymmetric branch: direct feed into the actual bottleneck

The grouped real weak-axisymmetric `Y_20` branch carries the exact signature
\[
x_{20}=x^{(0)}+\epsilon x^{(1)},\qquad
x_{21}=x^{(0)}+\frac{\epsilon}{2}x^{(1)},\qquad
x_{22}=x^{(0)}-\epsilon x^{(1)},
\]
so its grouped anomalies obey
\[
\boxed{b_x=3a_x,}
\]
while any rotational scalar extracted from the isotropic branch has no linear
\(P_2\) feed-down. fileciteturn5file12turn4file19

For the projected-Maxwell piece, this means:

\[
Z_{A,0}=Z_0+\epsilon\lambda_A z_0^{(1)},\qquad
N_{A,0}=N_0+\epsilon\lambda_A n_0^{(1)},
\qquad
(\lambda_{20},\lambda_{21},\lambda_{22})=\left(1,\frac12,-1\right).
\]

Then lanewise
\[
P_A=\frac{N_{A,0}}{D_{A,0}}
\]
becomes
\[
P_A
=
P_0\Bigl(1+\epsilon\lambda_A\,\Xi_1^{(\mathrm{proj})}\Bigr),
\]
with
\[
\boxed{
\Xi_1^{(\mathrm{proj})}
=
\frac{n_0^{(1)}}{N_0}
+
\frac{z_0^{(1)}}{D_0}.
}
\]

That is the direct projected-Maxwell contribution to the same weak-axisymmetric prefactor
slope that the later moving-throat notes compress to
\[
\Xi_1=\frac{P_1}{P_0}.
\]
The 5PN / actual-branch notes already identify this as the transported static bottleneck,
not a higher dynamic window. fileciteturn4file5turn0file16

So the second big conclusion is:

> **projection-first near-throat Maxwell can feed the actual weak-axisymmetric PDE bottleneck
> \(\Xi_1\) directly, with the correct grouped signature \(b=3a\).**

This is one of the strongest reasons to keep pushing this route.

---

## 6. Why this looks like a real missing piece

The project files already say that:

- mixed-sector electromagnetism, \(A_w\), \(F_{\mu w}\), and \(J^w\) are real microscopic
  channels outside the strict far-field brane reduction; projection makes the brane an open
  subsystem with leakage and unresolved-stress corrections. fileciteturn2file0turn2file3
- the grouped real `P2` lane and the full-bundle coefficients \(D_0,D_2,D_4,N_0,\dots\) are
  the right conservative / outgoing bookkeepers for the moving-throat PDE. fileciteturn5file2turn5file19turn3file4
- the current parent action gives a wall force but not yet an autonomous moving-throat PDE unless
  a wall/throat action is added, so any honest missing microscopic ingredient has to enter through
  a reduced closure like this before it can be promoted back upward. fileciteturn0file18
- the gauge sector is cleanest in a localized form \(H=Z\), which keeps the projected/reduced
  gauge parameter finite and aligned while leaving the mixed gauge invariants \(E_w\) and \(C_a\)
  fully alive. fileciteturn2file18turn2file19

Taken together, that makes the current bridge fairly compelling:

1. **Near-throat projection-first Maxwell supplies exactly the mixed-sector throat terms the
   far-field reduction suppresses too early.**
2. **Those terms feed exactly the \(Z_n\) and \(N_n\) bundle slots that the later PDE notes leave
   as the live conservative / outgoing bottleneck.**
3. **At isotropic level, the first quantities to inspect are \(z_2,z_4,n_0\).**
4. **At weak-axisymmetric level, the first quantity to inspect is the projected-Maxwell contribution
   to \(\Xi_1\).**

That is about as clean a “next theorem gate” as we can ask for without pretending the full
nonlinear throat branch is already solved.

---

## 7. The next concrete move

The most useful next derivation after this one would be:

1. pick a concrete throat-local projected ansatz for the mixed-sector data,
   enough to compute explicit
   \[
   z_0,\ z_2,\ z_4,\ n_0,\ n_2,\ n_4,
   \]
   or at least their first mouth-local slopes;
2. impose the structurally clean gauge choice \(H=Z\);
3. test whether the resulting
   \[
   \delta\mathcal C,\qquad \delta P_2,\qquad \delta P_4,\qquad \Xi_1^{(\mathrm{proj})}
   \]
   move the branch in the right direction.

If that succeeds, then projection-first Maxwell is not just a side calculation. It becomes
a credible microscopic source for the remaining PDE-side bundle data near the throat.
