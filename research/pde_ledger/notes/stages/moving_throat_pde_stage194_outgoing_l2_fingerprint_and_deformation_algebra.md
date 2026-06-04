# Moving-Throat PDE — Stage 194: Exact Outgoing `l=2` DtN Fingerprint, Fixing of `\chi_Q`, and the Isotropic Deformation Algebra

## Status

**Exact within the carried Stage 193 isotropic grouped-real `P2` conservative surface, the reduced 2.5PN outgoing compiler, and the canonical compact passive/outgoing `l=2` DtN closure already isolated in the 5PN notes.**

This stage does **not** introduce a new constitutive law.
It upgrades Stage 193 from the exact conservative isotropic one-pole front end to the first exact **retarded** isotropic branch test.

---

## Purpose

Stage 193 froze the conservative grouped-real `P2` target surface
\[
 a_2=b_2=a_4=b_4=0,
 \qquad
 \Delta_{\rm pole}=\bar\nu_4-4\bar\nu_2^{\,2}=0,
\]
and therefore the exact one-parameter conservative carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2}.
\]
That was the right front end, but it still left one reduced scalar open:

> once the grouped-real `P2` branch is isotropic and one-pole conservative, what fixes the **leading outgoing `l=2` normalization** on the actual passive/outgoing branch?

The 5PN chain narrowed that remaining freedom to the single scalar
\[
\chi_Q,
\]
which multiplies the leading odd `\omega^5` term in the isotropic retarded grouped-`P2` module. This stage freezes the next audited theorem target:

1. the exact outgoing spherical `l=2` Dirichlet-to-Neumann fingerprint,
2. the exact fixing
   \[
   \chi_Q=1
   \]
   on the canonical compact passive/outgoing branch,
3. and the exact first isotropic deformation algebra showing which DtN-side branch data can move `\chi_Q` while preserving the canonical even moments.

So Stage 194 is the first precise retarded continuation of the Stage 193 conservative surface.

---

## 1. Carry-forward retarded grouped-`P2` module

From Stage 193, the isotropic one-pole conservative carrier is already fixed through `O(\omega^4)`.
The corresponding retarded grouped-`P2` one-pole-plus-contact module is written as
\[
\boxed{
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34
+
\frac14\,
\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\,\sigma_Q^{\rm can}\,\omega^5}
+O(\omega^6).
}
\]
Here:

- `\Omega_Q` is the Stage 193 conservative pole scale,
- `\sigma_Q^{\rm can}` is the canonical compact outgoing normalization,
- `\chi_Q` is the only reduced outgoing-normalization scalar that remains open before the explicit DtN model is inserted.

Matching the even coefficients of Stage 193 to the canonical compact outgoing `l=2` branch fixes
\[
\boxed{
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}=\frac{4a^5}{27c_s^5}.
}
\]
Therefore the retarded module expands as
\[
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]
So the only leading retarded mismatch relative to the canonical compact branch is the multiplier `\chi_Q`.

---

## 2. Exact outgoing spherical `l=2` DtN fingerprint

Define the dimensionless outgoing argument
\[
\boxed{z:=\frac{a\omega}{c_s}.}
\]
Let the canonical compact outgoing partial wave be the spherical Hankel mode
\[
h_2^{(1)}(z)=j_2(z)+i\,y_2(z).
\]
Then the exact outgoing `l=2` DtN operator is
\[
\boxed{
\Lambda_2^{\rm out}(z)
=
 z\,\frac{d}{dz}\ln h_2^{(1)}(z)
=
 z\,\frac{h_2^{(1)\prime}(z)}{h_2^{(1)}(z)}.
}
\]
Its small-`z` expansion is
\[
\boxed{
\Lambda_2^{\rm out}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}
+i\,\frac{z^5}{9}
-\frac{2z^6}{27}
-i\,\frac{z^7}{27}
+O(z^8).
}
\]
So the exact static slot is
\[
\Lambda_2^{\rm out}(0)=-3.
\]

Normalize by that static slot:
\[
\boxed{
\widehat Y_2^{\rm out}(z)
:=-\frac{3}{\Lambda_2^{\rm out}(z)}.
}
\]
Then
\[
\boxed{
\widehat Y_2^{\rm out}(z)
=
1+\frac{z^2}{9}
+\frac{4z^4}{81}
+i\,\frac{z^5}{27}
-\frac{11z^6}{729}
-i\frac{z^7}{243}
+O(z^8).
}
\]
Restoring `\omega`,
\[
\boxed{
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
}
\]
This is the exact canonical compact outgoing `l=2` fingerprint.

---

## 3. Exact fixing of `\chi_Q` on the canonical compact branch

Compare the Stage 194 retarded grouped-`P2` module
\[
\widehat Y_Q^{\rm ret}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\chi_Q\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6)
\]
with the explicit outgoing DtN branch
\[
\widehat Y_2^{\rm out}(\omega)
=
1+\frac{a^2\omega^2}{9c_s^2}
+\frac{4a^4\omega^4}{81c_s^4}
+i\,\frac{a^5\omega^5}{27c_s^5}
+O(\omega^6).
\]
Matching the leading odd coefficient yields
\[
\boxed{\chi_Q=1.}
\]
So on the canonical compact passive/outgoing `l=2` DtN branch, the last reduced outgoing-normalization scalar is fixed exactly.

Equivalently, any deformed outgoing DtN branch of the form
\[
\Lambda_2^{\rm def}(z)
=
-3+\frac{z^2}{3}+\frac{z^4}{9}+i\,\xi_Q\,\frac{z^5}{9}+O(z^6)
\]
simply produces
\[
\boxed{\chi_Q=\xi_Q.}
\]
So after the canonical DtN model is inserted, the only reduced outgoing uncertainty left is deviation of the **actual** moving-throat isotropic branch from the exact compact outgoing fingerprint.

---

## 4. Exact isotropic DtN deformation algebra

To make that last statement explicit, take the first isotropic moving-throat DtN deformation family
\[
\boxed{
\Lambda_2^{\rm def}(z)
=
S\,\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6).
}
\]
Here:

- `S` is an overall mouth normalization,
- `\beta` rescales the effective outgoing argument,
- `\Sigma_0,\Sigma_2,\Sigma_4` are isotropic even throat-core self-energy data,
- `\Sigma_5` is the first extra isotropic odd `l=2` core outlet.

Expanding,
\[
\Lambda_2^{\rm def}(z)=L_0+L_2 z^2+L_4 z^4+iL_5 z^5+O(z^6),
\]
with
\[
L_0=-3S+\Sigma_0,
\qquad
L_2=\frac{S\beta^2}{3}+\Sigma_2,
\qquad
L_4=\frac{S\beta^4}{9}+\Sigma_4,
\qquad
L_5=\frac{S\beta^5}{9}+\Sigma_5.
\]
Normalize by the actual static slot:
\[
\widehat Y_2^{\rm def}(z):=\frac{L_0}{L_0+L_2 z^2+L_4 z^4+iL_5 z^5}+O(z^6).
\]
Then the exact low-frequency coefficients are
\[
\boxed{
\widehat Y_2^{\rm def}(z)
=
1-\frac{L_2}{L_0}z^2
+\left(\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}\right)z^4
-i\frac{L_5}{L_0}z^5
+O(z^6).
}
\]

### 4.1 Canonical-even matching conditions

Demand that the deformed branch preserve the canonical even fingerprint
\[
\frac{z^2}{9},
\qquad
\frac{4z^4}{81}.
\]
Then
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81}.
\]
Solving for the even isotropic core coefficients gives
\[
\boxed{
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
}
\]
So the even deformations are not free once the canonical conservative moments are held fixed.

### 4.2 Exact outgoing-normalization map

With those canonical-even constraints imposed, the retarded normalization scalar is
\[
\boxed{
\chi_Q=
\frac{3\big(S\beta^5+9\Sigma_5\big)}{3S-\Sigma_0}.
}
\]
Equivalently,
\[
\boxed{
\chi_Q-1=
\frac{3S(\beta^5-1)+\Sigma_0+27\Sigma_5}{3S-\Sigma_0}.
}
\]
So the only isotropic DtN-side branch data that can move the canonical outgoing normalization are:

1. argument deformation `\beta`,
2. static additive core shift `\Sigma_0`,
3. odd isotropic `l=2` core outlet `\Sigma_5`.

The overall scale `S` is not a separate obstruction by itself. It matters only through the normalized ratios above.

---

## 5. Carry-forward corollary on the canonical branch

If one also carries the natural point-particle source-map branch
\[
\hat m_0\to 1
\]
from the reduced 2.5PN package, then the exact fixing `\chi_Q=1` immediately returns the canonical invariant tuple
\[
\boxed{
\overline K_0=\frac{54Gc_s^5}{5a^5c^5},
\qquad
\overline K_2=\frac{6Gc_s^3}{5a^3c^5},
\qquad
\overline K_4=\frac{8Gc_s}{15ac^5},
\qquad
\overline\Gamma_5=\frac{2G}{5c^5}.
}
\]
This is not a new theorem input in Stage 194. It is the immediate reduced corollary of:

- the Stage 193 isotropic one-pole conservative surface,
- the exact outgoing `l=2` fingerprint,
- and the already-carried natural point-particle source-map branch.

---

## 6. What Stage 194 changes in the theorem problem

Stage 193 left the reduced hierarchy with one isotropic conservative carrier but one unfixed retarded normalization scalar.
Stage 194 changes that in three precise ways.

### 6.1 The canonical outgoing coefficient is no longer symbolic

The exact outgoing spherical `l=2` DtN model fixes the canonical odd coefficient directly. So `\chi_Q` is no longer a placeholder on the compact passive/outgoing branch.

### 6.2 The canonical compact passive/outgoing branch is fully normalized

Matching the exact DtN fingerprint to the retarded grouped-`P2` module gives
\[
\chi_Q=1.
\]
So on that branch the last reduced outgoing scalar is closed.

### 6.3 The remaining isotropic PDE-facing freedom is sharply localized

If the actual moving-throat branch deviates from the canonical outgoing DtN branch while preserving the canonical even moments, then the full isotropic retarded defect is carried only by
\[
(\beta,\Sigma_0,\Sigma_5).
\]
So the remaining isotropic DtN realization problem is not diffuse. It is an explicit three-parameter deformation problem.

---

## 7. Immediate next derivation step

The next clean stage is now much narrower than before.

Do **not** reopen the conservative grouped bundle.
Instead:

1. keep the Stage 193 isotropic conservative one-pole surface,
2. keep the Stage 194 exact outgoing `l=2` fingerprint,
3. compute the actual moving-throat isotropic DtN data that feed
   \[
   (\beta,\Sigma_0,\Sigma_5),
   \]
4. and test whether the realized branch lands on the canonical submanifold
   \[
   \chi_Q=1.
   \]

That is the smallest honest next theorem gate on the retarded side.
