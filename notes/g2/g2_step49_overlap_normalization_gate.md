# Step 49 — The overlap-integral theorem gate and first anisotropy diagnostic for the adiabatic branch

## Goal

Step 48 pulled the minimal outgoing observable surface back to the microscopic
bundle closure
\[
(D_0,D_2,D_4,N_0,N_2,N_4).
\]

The next honest move is to go one level deeper and rewrite that bundle closure in
the **actual overlap-level variables** of the isotropic moving-throat branch:

- wall baseline \(K,M\),
- BdG support moments \((B_0,B_2,B_4)\),
- conservative Maxwell/mixed moments \((Z_0,Z_2,Z_4)\),
- outgoing transfer moments \((N_0,N_2,N_4)\).

That is exactly the theorem gate Stage 6/7 was pointing to. This step makes it
explicit inside the g-2 chain.

---

## Step 49A — Exact isotropic overlap-level bundle

On the natural isotropic branch the grouped lanes collapse to common coefficients,
and the reduced bundle is

\[
\boxed{
D_0 = K - B_0 - Z_0,
}
\]
\[
\boxed{
D_2 = -(M + B_2 + Z_2),
}
\]
\[
\boxed{
D_4 = -(B_4 + Z_4).
}
\]

So the static outgoing prefactor is simply
\[
\boxed{
P_0 = \frac{N_0}{D_0}
=
\frac{N_0}{K-B_0-Z_0}.
}
\]

This is the microscopic form of the final normalization ratio already isolated by
the moving-throat notes.

---

## Step 49B — Exact overlap-integral normalization gate

The universal quadrupole target is therefore

\[
\boxed{
mhat_0^2\,\frac{N_0}{K-B_0-Z_0}
=
\frac{54Gc_s^5}{5a^5c^5}.
}
\]

That is the sharpest version of the theorem gate we have reached.

It says the full moving-throat branch does **not** need to be judged first by a huge
space of coefficients. On the isotropic branch, the entire normalization test is the
single competition between

- outgoing-transfer weight \(N_0\),
- and conservative dressed wall stiffness \(K-B_0-Z_0\).

This also gives an immediate physical reading:

- increasing \(N_0\) moves the branch upward toward the target,
- increasing \(B_0\) or \(Z_0\) softens the denominator and also moves the branch upward,
- but stability still requires \(D_0=K-B_0-Z_0>0\).

So the remaining theorem gap is now a concrete branch-balance question, not a
diffuse “more PDE” request.

---

## Step 49C — Constant-prefactor branch conditions in pure overlap variables

From Step 48, the minimal outgoing branch satisfies
\[
P_2=P_4=0.
\]

Substituting the isotropic overlap-level bundle gives

\[
\boxed{
N_2
=
-\frac{2(M+B_2+Z_2)N_0}{K-B_0-Z_0},
}
\]

and

\[
\boxed{
N_4
=
\frac{\bigl(M+B_2+Z_2\bigr)^2
-2\bigl(K-B_0-Z_0\bigr)\bigl(B_4+Z_4\bigr)}
{\bigl(K-B_0-Z_0\bigr)^2}\,N_0.
}
\]

So if the minimal constant-prefactor branch is the correct one, the higher outgoing
transfer moments are already fixed by the conservative overlap data. They are not
free branch functions.

That is a very strong reduction of the moving-throat search space.

---

## Step 49D — First weak-axisymmetric anisotropy diagnostic

Stage 7 also showed that the first natural symmetry-breaking pattern is a weak
axisymmetric quadrupole deformation with grouped-lane weights
\[
\lambda_{20}=1,\qquad
\lambda_{21}=\frac12,\qquad
\lambda_{22}=-1.
\]

If the normalization data deform as
\[
D_A = D_0 + \epsilon \lambda_A D_1 + O(\epsilon^2),
\qquad
N_A = N_0 + \epsilon \lambda_A N_1 + O(\epsilon^2),
\]
then
\[
P_A=\frac{N_A}{D_A}=P_0+\epsilon\lambda_A P_1+O(\epsilon^2),
\qquad
P_1=\frac{N_1D_0-N_0D_1}{D_0^2}.
\]

The grouped defects are therefore
\[
\boxed{
a_P=\frac{\epsilon}{4}P_1,
\qquad
b_P=\frac{3\epsilon}{4}P_1,
}
\]
so the first normalization anisotropy must satisfy
\[
\boxed{
b_P = 3 a_P.
}
\]

That gives a direct near-failure diagnostic for future PDE outputs:

- if weak grouped-lane normalization anisotropy appears and **does** satisfy
  \(b_P=3a_P\), it is consistent with a simple axisymmetric quadrupole deformation
  of the isotropic branch;
- if it **fails**, then the deformation is already more complicated than that
  minimal pattern.

---

## Main result of the step

The next honest moving-throat theorem gate for the adiabatic g-2 program is the
exact overlap-level system

\[
\boxed{
mhat_0^2\,\frac{N_0}{K-B_0-Z_0}
=
\frac{54Gc_s^5}{5a^5c^5},
}
\]
together with the minimal constant-prefactor conditions
\[
\boxed{
N_2
=
-\frac{2(M+B_2+Z_2)N_0}{K-B_0-Z_0},
}
\]
\[
\boxed{
N_4
=
\frac{\bigl(M+B_2+Z_2\bigr)^2
-2\bigl(K-B_0-Z_0\bigr)\bigl(B_4+Z_4\bigr)}
{\bigl(K-B_0-Z_0\bigr)^2}\,N_0,
}
\]
and the weak-axisymmetric diagnostic
\[
\boxed{
b_P = 3 a_P.
}
\]

So the next PDE-facing falsification test is no longer vague at all:

1. compute the actual overlap data \((B_n,Z_n,N_n)\),
2. check grouped isotropy on the natural branch,
3. test the single ratio \(mhat_0^2N_0/(K-B_0-Z_0)\),
4. and, if weak anisotropy survives, verify whether it obeys \(b_P=3a_P\).

That is the smallest honest next theorem gate after Step 49.
