# Moving-Throat PDE — Stage 179: Exact Higher-Odd Irrelevance Theorem, Stability of `\chi_Q`, and the Collapse of the 2.5PN Retarded Finish Line to `\Delta_Q=\chi_Q-1`

## Status

**Exact within the carried Stage-174 Packet-A branch-residual hierarchy, the Stage-176 isotropic grouped-real `P2` conservative one-pole surface, the Stage-177 exact outgoing `l=2` DtN fingerprint and isotropic deformation algebra, the Stage-178 exact source-map reduction, and the reduced 2.5PN outgoing compiler.**

This stage does **not** introduce a new constitutive law.
It answers the last natural objection left by Stage 178:

> what if the actual moving-throat quadrupole branch carries additional isotropic retarded structure beginning at `O(\omega^7)` or higher?

The answer is that such higher odd data are invisible to the reduced point-particle `2.5`PN theorem.
So the only live retarded obstruction that survives the full carried hierarchy is still the leading outgoing-normalization scalar
\[
\chi_Q.
\]

---

## Purpose

Stage 176 froze the exact isotropic grouped-real `P2` conservative one-pole carrier,
Stage 177 fixed the canonical compact passive/outgoing `l=2` fingerprint and isolated the outgoing scalar
\[
\chi_Q,
\]
and Stage 178 then reduced the last observable normalization defect to the source/outgoing factorization
\[
 m_{\hat 0}^{\,2}\chi_Q N_Q=1,
\qquad
\Delta_{\rm norm}=P_0^{\rm target}(\chi_Q^{-1}-1)
\]
once the odd observable condition is imposed.

So after Stage 178 the last reduced `2.5`PN question is already sharp:

> does anything beyond the leading `\omega^5` outgoing slot survive into the reduced theorem, or is the whole retarded finish line really measured by `\chi_Q-1` alone?

This stage answers that exactly.

The main outputs are:

1. an exact **response-side higher-odd difference identity** showing that any extra retarded tail `O(\omega^7)` or higher changes the isotropic grouped response only at `O(\omega^7)`,
2. an exact **DtN-side higher-odd difference identity** showing the same at the operator level,
3. the exact theorem that the Stage-177 outgoing-normalization compiler
   \[
   \chi_Q=
   \frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}
   \]
   is unchanged by all higher odd DtN data beginning at `O(z^7)`,
4. the exact Packet-A consequence that the source-map-reduced normalization residual of Stage 178 is likewise unchanged,
5. and the final reduced statement that the only live retarded obstruction at `2.5`PN is
   \[
   \boxed{\Delta_Q:=\chi_Q-1.}
   \]

So Stage 179 is the rigorous Stage-174–178 reformulation of the older higher-odd irrelevance result.

---

## 1. Generalized isotropic grouped-`P2` retarded module with a higher odd tail

Keep the Stage-177 isotropic retarded grouped module,
\[
\widehat Y_Q^{\rm ret,5}(\omega)
=
\frac34+
\frac14\,
\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5},
\]
and now allow an arbitrary extra isotropic retarded tail beginning at order `\omega^7`:
\[
\boxed{
\widehat Y_Q^{\rm ret,\ge 7}(\omega)
=
\frac34+
\frac14\,
\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5-i\,\mathcal T_{\ge7}(\omega)},
}
\]
with
\[
\boxed{\mathcal T_{\ge7}(\omega)=O(\omega^7).}
\]

The Stage-177 canonical outgoing scale is still
\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}=\frac{4a^5}{27c_s^5}.
\]

A useful one-coefficient representative is
\[
\mathcal T_{\ge7}(\omega)=\tau_Q\omega^7,
\]
but the theorem below is not restricted to a monomial tail.
It applies to **any** isotropic retarded correction whose first nonzero term is `O(\omega^7)`.

---

## 2. Exact response-side higher-odd difference identity

Define
\[
X_Q(\omega):=\frac{\omega^2}{\Omega_Q^2}+i\chi_Q\sigma_Q^{\rm can}\omega^5,
\qquad
H_Q(\omega):=i\,\mathcal T_{\ge7}(\omega).
\]
Then
\[
\widehat Y_Q^{\rm ret,5}(\omega)=\frac34+\frac14\frac{1}{1-X_Q(\omega)},
\]
while
\[
\widehat Y_Q^{\rm ret,\ge 7}(\omega)=\frac34+\frac14\frac{1}{1-X_Q(\omega)-H_Q(\omega)}.
\]
Subtracting gives the exact identity
\[
\boxed{
\widehat Y_Q^{\rm ret,\ge 7}(\omega)-\widehat Y_Q^{\rm ret,5}(\omega)
=
\frac{H_Q(\omega)}{4\,(1-X_Q(\omega))\,(1-X_Q(\omega)-H_Q(\omega))}.
}
\]
Since
\[
H_Q(\omega)=O(\omega^7),
\qquad
1-X_Q(\omega)=1+O(\omega^2),
\]
it follows immediately that
\[
\boxed{
\widehat Y_Q^{\rm ret,\ge 7}(\omega)-\widehat Y_Q^{\rm ret,5}(\omega)=O(\omega^7).
}
\]
So all coefficients through `O(\omega^5)` are unchanged.
In particular,
\[
\boxed{
\widehat Y_Q^{\rm ret,\ge 7}(\omega)
=
1+
\frac{\omega^2}{4\Omega_Q^2}
+
\frac{\omega^4}{4\Omega_Q^4}
+
 i\chi_Q\frac{9}{32\Omega_Q^5}\omega^5
+
O(\omega^6).
}
\]

If one keeps the first higher odd tail explicitly,
\[
\mathcal T_{\ge7}(\omega)=\tau_Q\omega^7,
\]
then
\[
\widehat Y_Q^{\rm ret,\ge 7}(\omega)
=
1+
\frac{\omega^2}{4\Omega_Q^2}
+
\frac{\omega^4}{4\Omega_Q^4}
+
 i\chi_Q\frac{9}{32\Omega_Q^5}\omega^5
+
\frac{\omega^6}{4\Omega_Q^6}
+
 i\left(
\frac{9\chi_Q}{16\Omega_Q^7}+\frac{\tau_Q}{4}
\right)\omega^7
+
O(\omega^8).
\]
So the first place the extra higher odd tail can appear is exactly `\omega^7`.
It does not modify the leading odd Burke–Thorne slot.

---

## 3. Exact DtN-side higher-odd difference identity

Work now at the isotropic DtN level.
Let the Stage-177 truncated isotropic DtN denominator be
\[
D_5(z):=L_0+L_2 z^2+L_4 z^4+iL_5 z^5,
\]
with normalized grouped response
\[
\widehat Y_2^{\rm def,5}(z):=\frac{L_0}{D_5(z)}.
\]
Allow now a completely arbitrary higher odd DtN remainder beginning at `O(z^7)`:
\[
\boxed{\mathcal L_{\ge7}(z)=O(z^7),}
\]
and define
\[
D_{\ge7}(z):=D_5(z)+\mathcal L_{\ge7}(z),
\qquad
\widehat Y_2^{\rm def,\ge7}(z):=\frac{L_0}{D_{\ge7}(z)}.
\]
Then the exact difference identity is
\[
\boxed{
\widehat Y_2^{\rm def,\ge7}(z)-\widehat Y_2^{\rm def,5}(z)
=
-\frac{L_0\,\mathcal L_{\ge7}(z)}{D_5(z)\,D_{\ge7}(z)}.
}
\]
Since `\mathcal L_{\ge7}(z)=O(z^7)` and both denominators are nonzero deformations of `L_0` near `z=0`,
\[
\boxed{
\widehat Y_2^{\rm def,\ge7}(z)-\widehat Y_2^{\rm def,5}(z)=O(z^7).
}
\]
So the normalized DtN response through `O(z^5)` is unchanged.
Equivalently,
\[
\boxed{
\widehat Y_2^{\rm def,\ge7}(z)
=
1-
\frac{L_2}{L_0}z^2
+
\left(\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}\right)z^4
-
 i\frac{L_5}{L_0}z^5
+
O(z^6),
}
\]
with no dependence on `\mathcal L_{\ge7}` before `O(z^7)`.

---

## 4. Stability of the Stage-177 isotropic deformation algebra

Parameterize the Stage-177 isotropic DtN front end by
\[
L_0=-3S+\Sigma_0,
\qquad
L_2=\frac{S\beta^2}{3}+\Sigma_2,
\qquad
L_4=\frac{S\beta^4}{9}+\Sigma_4,
\qquad
L_5=\frac{S\beta^5}{9}+\Sigma_5,
\]
and leave the entire higher odd block inside the single remainder
\[
\mathcal L_{\ge7}(z)=O(z^7).
\]
Then the exact canonical-even matching conditions are still
\[
-\frac{L_2}{L_0}=\frac19,
\qquad
\frac{L_2^2}{L_0^2}-\frac{L_4}{L_0}=\frac{4}{81},
\]
so the same Stage-177 solution follows:
\[
\boxed{
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27}.
}
\]
These formulas do not contain any higher odd coefficient.

Likewise the outgoing-normalization factor is still compiled only from the `z^5` coefficient:
\[
\boxed{
\chi_Q=-27\,\frac{L_5}{L_0}
=\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}.
}
\]
So all higher odd DtN data beginning at `O(z^7)` are invisible to the Stage-177 isotropic outgoing-normalization compiler.
They can change the response only at `z^7` and higher, but they cannot change `\chi_Q`.

---

## 5. Exact Packet-A consequence after Stage 178

Stage 178 already reduced the odd observable closure to
\[
\boxed{m_{\hat 0}^{\,2}\chi_Q N_Q=1,}
\]
and therefore the Packet-A normalization residual to
\[
\boxed{
\Delta_{\rm norm}=P_0^{\rm target}(\chi_Q^{-1}-1)
}
\]
once the odd observable condition is imposed, with
\[
P_0^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
\]
On the natural point-particle source-map branch,
\[
\boxed{m_{\hat 0}\to 1,\qquad N_Q=\chi_Q^{-1}.}
\]

Now combine this with Stage 179.
Since all higher odd retarded data beginning at `O(\omega^7)` or `O(z^7)` leave `\chi_Q` unchanged, they also leave unchanged:
\[
\boxed{N_Q,}
\qquad
\boxed{m_{\hat 0}^{\,2}\chi_Q N_Q,}
\qquad
\boxed{\Delta_{\rm norm}.}
\]
The remaining Packet-A entries are already conservative:
\[
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole}),
\]
so they are manifestly untouched by any extra higher odd retarded tail.

Therefore the entire Stage-174 branch residual packet
\[
\boxed{
\Delta_{\rm branch}
=
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole},\Delta_{\rm norm})
}
\]
is unchanged at the level relevant to the reduced point-particle `2.5`PN theorem.

---

## 6. PN meaning of the higher odd tail

The Burke–Thorne slot is the leading odd grouped-`P2` coefficient,
\[
 i\,\Gamma_5\omega^5.
\]
Any new isotropic retarded structure that first appears as
\[
 i\,\Gamma_7\omega^7
\quad\text{or higher}
\]
contains two extra powers of `\omega` relative to the `\omega^5` reaction term.
So in the standard small-velocity bookkeeping it sits **above** the universal point-particle `2.5`PN slot.

Stage 179 does not need to decide the full higher-PN interpretation of those terms.
It only needs the exact reduced statement:

> they cannot alter the `\omega^5` coefficient and therefore cannot alter the reduced point-particle `2.5`PN closure.

So the completed moving-throat PDE may well carry additional higher odd retarded data, but such data live beyond the theorem order that the present Packet-A endgame is testing.

---

## 7. Exact higher-odd irrelevance theorem

Let the isotropic grouped-real `P2` branch satisfy the carried Stage-176 conservative one-pole surface, the Stage-177 outgoing `l=2` compiler, and the Stage-178 source-map reduction.
Let the exact retarded grouped response differ from the Stage-177 one-pole model only by additional isotropic retarded structure whose first nonzero term is `O(\omega^7)` (equivalently `O(z^7)`).

Then:

1. the normalized grouped response through `O(\omega^5)` is unchanged,
2. the isotropic outgoing-normalization factor `\chi_Q` is unchanged,
3. the exact source-map-reduced Packet-A normalization residual `\Delta_{\rm norm}` is unchanged,
4. and the entire reduced point-particle `2.5`PN theorem remains sensitive only to the leading outgoing-normalization defect
   \[
   \boxed{\Delta_Q:=\chi_Q-1.}
   \]

So the only live retarded obstruction at reduced `2.5`PN order is still the single scalar `\chi_Q`.
If
\[
\boxed{\chi_Q=1,}
\]
then the reduced point-particle `2.5`PN theorem is closed regardless of all higher odd retarded data beginning at `O(\omega^7)`.

---

## 8. What Stage 179 changes in the theorem problem

Stage 178 had already reduced the finish line to the outgoing scalar `\chi_Q-1`, but one natural loophole remained: perhaps uncomputed higher odd retarded structure could re-enter the Packet-A compiler indirectly.

Stage 179 closes that loophole exactly.

### 8.1 The retarded finish line is now unique at `2.5`PN

No extra isotropic retarded datum beginning at `O(\omega^7)` can modify the reduced point-particle `2.5`PN verdict.
So the retarded finish line is not a family of hidden coefficients.
It is one number only:
\[
\chi_Q.
\]

### 8.2 The Stage-177 deformation algebra is complete at theorem order

The only isotropic DtN-side data that can move the reduced `2.5`PN verdict are still
\[
(\beta,\Sigma_0,\Sigma_5).
\]
All higher odd DtN data are outside the theorem order.

### 8.3 The Packet-A home-stretch theorem is stable

The Stage-174 minimal branch-residual packet and the Stage-178 source-map reduction already contain the whole reduced `2.5`PN theorem.
There is no hidden higher-odd loophole left inside that packet language.

---

## 9. Immediate next derivation step

The next clean continuation is now the exact conditional closure statement in the modern Packet-A language:

1. keep the Stage-176 isotropic conservative one-pole surface,
2. keep the Stage-178 exact source-map reduction,
3. keep the Stage-179 higher-odd irrelevance theorem,
4. and state the finish-line equivalence
   \[
   \Delta_{\rm branch}=0
   \iff
   \chi_Q=1
   \]
   on the natural point-particle source-map branch.

That is the natural successor to Stage 179.
