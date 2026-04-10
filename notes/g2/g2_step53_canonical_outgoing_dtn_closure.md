# Step 53 — The canonical outgoing DtN branch is not just plausible, it is fixed

## Goal

Step 52 showed that the canonical outgoing grouped-`P2` branch is already a no-tuning
background and that the remaining electron sliver, if real, sits in one outgoing
branch datum. The next honest move is to stop treating the canonical branch as a
symbol and derive it directly from the explicit compact outgoing `l=2`
Dirichlet-to-Neumann model.

That is exactly what the later moving-throat notes do: they derive the outgoing
`l=2` DtN fingerprint, match it to the minimal grouped-`P2` retarded module, and
fix the outgoing-normalization factor to
\[
\chi_Q=1
\]
on the canonical compact branch. They then insert that into the reduced normalization
stack and show that, on the natural point-particle source-map branch,
\[
N_Q=1.
\]
fileciteturn20file0 fileciteturn20file5

---

## Step 53A — Exact outgoing `l=2` DtN fingerprint

The explicit outgoing spherical `l=2` DtN operator is
\[
\Lambda_2^{\rm out}(z)
=
z\frac{d}{dz}\ln h_2^{(1)}(z),
\qquad
z=\frac{a\omega}{c_s},
\]
with low-frequency expansion
\[
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\frac{z^5}{9}
+O(z^6).
\]

Normalizing by the static slot gives
\[
\widehat Y_2^{\rm out}(z)
=
-\frac{3}{\Lambda_2^{\rm out}(z)}
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\frac{z^5}{27}
+O(z^6).
\]

So the canonical outgoing odd coefficient is not free.
It is fixed directly by the explicit DtN model to
\[
\Gamma_{5,\rm can}^{\rm DtN}
=
\frac{a^5}{27c_s^5}.
\]
fileciteturn20file0

---

## Step 53B — Matching to the minimal isotropic grouped-`P2` branch

The retarded minimal grouped-`P2` module can be written as
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34
+
\frac14\,
\frac{1}{
1-\omega^2/\Omega_Q^2
-i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5
}
+O(\omega^6),
\]
with the already-fixed conservative moment match
\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{4a^5}{27c_s^5}.
\]

Expanding in the same variable \(z=a\omega/c_s\) gives
\[
\widehat Y_Q^{\rm ret}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\chi_Q\,\frac{z^5}{27}
+O(z^6).
\]

Matching the explicit DtN branch then forces
\[
\boxed{\chi_Q=1.}
\]
So the canonical outgoing grouped-`P2` branch is not just a convenient choice.
It is exactly the compact outgoing DtN branch. fileciteturn20file0turn20file18

---

## Step 53C — On the natural source-map branch this also fixes `N_Q`

The reduced normalization stack is
\[
\hat m_0^{\,2}\chi_Q N_Q=1.
\]
On the natural point-particle source-map branch
\[
\hat m_0=1+O(a^2/r^2),
\]
so in the strict point-particle limit the same relation becomes
\[
N_Q=\frac{1}{\chi_Q}.
\]

With the explicit DtN result
\[
\chi_Q=1,
\]
this collapses to
\[
\boxed{N_Q=1.}
\]

The later notes therefore fix the canonical invariant coefficients exactly:
\[
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},
\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},
\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
\]
fileciteturn20file5

---

## What this changes for the g-2 chain

This is a real strengthening of Step 52.

Before Step 53, the canonical no-tuning branch was a very plausible reduced branch.
After Step 53, it is the explicit compact outgoing DtN branch itself.

So the best honest status is now:

- the canonical outgoing background is naturally fixed;
- the canonical background predicts
  \[
  \chi_Q=1,\qquad N_Q=1;
  \]
- therefore any nonzero electron-point sliver must come from a genuine deformation
  away from that canonical outgoing branch, not from ambiguity in the canonical
  branch itself.
