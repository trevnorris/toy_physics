# Step 48 — Pulling the minimal outgoing target surface back to the microscopic grouped-`P2` bundle

## Goal

Step 46 compressed the adiabatic constant-prefactor branch to the observable outgoing
surface
\[
\Gamma_5=\frac{2G}{5c^5},
\qquad
K_2^2=\frac14 K_0K_4,
\]
and Step 47 showed that the adiabatic anomaly moves the even bundle along a
one-parameter flow while leaving \(\Gamma_5\) fixed.

But those were still **observable-side** statements. The next honest move is to pull
that surface back to the microscopic grouped-`P2` bundle
\[
(D_0,D_2,D_4,N_0,N_2,N_4),
\]
so that the remaining theorem gate can be read directly at the bundle level.

This step does that.

---

## Step 48A — Exact bundle-to-outgoing map

From the grouped-`P2` bundle formulas we already have
\[
P_0=\frac{N_0}{D_0},
\]
\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]
\[
P_4=
\frac{D_0^2N_4-2D_0(D_2N_2+D_4N_0)+3D_2^2N_0}{D_0^3}.
\]

The outgoing bundle is then
\[
K_0=P_0,
\qquad
K_2=P_2+A\,P_0,
\qquad
K_4=P_4+A\,P_2+B\,P_0,
\qquad
\Gamma_5=G_5\,P_0,
\]
where the compact outgoing `l=2` fingerprint is
\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

So the observable surface is the image of **two** microscopic ingredients:

1. the bundle closure \((P_0,P_2,P_4)\),
2. the compact outgoing port fingerprint \((A,B,G_5)\).

---

## Step 48B — Exact constant-prefactor microscopic closure

The minimal isotropic outgoing branch is the constant-prefactor branch
\[
P_2=P_4=0.
\]

Solving those equations gives the exact microscopic conditions
\[
\boxed{
N_2=\frac{2D_2N_0}{D_0},
}
\]
\[
\boxed{
N_4=\frac{(D_2^2+2D_0D_4)N_0}{D_0^2}.
}
\]

Equivalently, using the conservative grouped response moments
\[
u_2=-\frac{D_2}{D_0},
\qquad
u_4=\frac{D_2^2-D_0D_4}{D_0^2},
\]
the same conditions are
\[
\boxed{
N_2=-2u_2N_0,
\qquad
N_4=(3u_2^2-2u_4)N_0.
}
\]

So the constant-prefactor branch does **not** let `N_2` and `N_4` float
independently. Once the conservative bundle \((D_0,D_2,D_4)\) and the static
transfer weight \(N_0\) are fixed, the higher transfer moments are already slaved.

---

## Step 48C — The Step-46 target surface is automatic once the compact fingerprint is inserted

On the constant-prefactor branch we have
\[
K_0=P_0=\frac{N_0}{D_0},
\qquad
K_2=A\,K_0,
\qquad
K_4=B\,K_0,
\qquad
\Gamma_5=G_5\,K_0.
\]

Now insert the compact outgoing `l=2` fingerprint
\[
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4}=4A^2,
\qquad
G_5=\frac{a^5}{27c_s^5}.
\]

Then
\[
\boxed{
K_2=A\,K_0,
\qquad
K_4=4A^2K_0,
\qquad
\Gamma_5=G_5K_0,
}
\]
and therefore
\[
\boxed{
K_2^2=\frac14 K_0K_4.
}
\]

So the Step-46 outgoing target surface is **not** an extra independent miracle.
It is simply the observable image of:

- the microscopic constant-prefactor bundle closure, and
- the compact outgoing `l=2` port law.

That is the cleanest pullback of the outgoing surface we have so far.

---

## Step 48D — Microscopic port transport on the adiabatic anomaly track

The observable Step-47 flow becomes even cleaner at the microscopic port level.

From the adiabatic anomaly branch,
\[
P_0 \mapsto e^\ell P_0,
\qquad
c_s \mapsto c_s e^{\ell/5},
\qquad
a \mapsto a.
\]

So the compact outgoing fingerprint rescales as
\[
\boxed{
A \mapsto e^{-2\ell/5}A,
\qquad
B \mapsto e^{-4\ell/5}B,
\qquad
G_5 \mapsto e^{-\ell}G_5.
}
\]

Combining that with \(P_0\mapsto e^\ell P_0\) gives
\[
K_0 \mapsto e^\ell K_0,
\]
\[
K_2 \mapsto e^{3\ell/5}K_2,
\]
\[
K_4 \mapsto e^{\ell/5}K_4,
\]
\[
\Gamma_5 \mapsto \Gamma_5.
\]

So at the microscopic port level the adiabatic anomaly factorizes into:

1. one outgoing-normalization rescaling \(P_0\to e^\ell P_0\),
2. one universal compact-port retuning \((A,B,G_5)\to(e^{-2\ell/5}A,e^{-4\ell/5}B,e^{-\ell}G_5)\).

That makes the invariance of \(\Gamma_5\) completely transparent:
\[
\Gamma_5 = G_5 P_0
\]
is a product of two compensating scalings.

---

## Main result of the step

The minimal outgoing observable surface is the exact image of the microscopic
grouped-`P2` bundle closure

\[
\boxed{
P_2=P_4=0
\quad\Longleftrightarrow\quad
N_2=\frac{2D_2N_0}{D_0},
\qquad
N_4=\frac{(D_2^2+2D_0D_4)N_0}{D_0^2},
}
\]
combined with the compact outgoing `l=2` fingerprint
\[
\boxed{
A=\frac{a^2}{9c_s^2},
\qquad
B=\frac{4a^4}{81c_s^4},
\qquad
G_5=\frac{a^5}{27c_s^5}.
}
\]

On that branch,
\[
\boxed{
K_0=P_0,
\qquad
K_2=A K_0,
\qquad
K_4=B K_0,
\qquad
\Gamma_5=G_5 K_0,
}
\]
so
\[
\boxed{
K_2^2=\frac14 K_0K_4.
}
\]

And on the adiabatic anomaly track, the microscopic port transport law is
\[
\boxed{
P_0 \mapsto e^\ell P_0,
\qquad
A \mapsto e^{-2\ell/5}A,
\qquad
B \mapsto e^{-4\ell/5}B,
\qquad
G_5 \mapsto e^{-\ell}G_5.
}
\]

So the next theorem gate is no longer to guess the shape of the outgoing surface.
It is to compute the actual bundle data \((D_0,D_2,D_4,N_0)\) and test whether the
true moving-throat branch satisfies the exact microscopic closure above.
