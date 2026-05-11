# Moving-Throat PDE — Stage 195: Exact Source-Map Reduction of the Canonical Outgoing Branch, the Factorization `m_{\hat 0}^{\,2}\chi_Q N_Q=1`, and the Collapse of `\Delta_{\rm norm}`

## Status

**Exact within the carried Stage 242 Packet-A branch-residual hierarchy, the Stage 244 isotropic grouped-real `P2` conservative surface, the Stage 245 exact outgoing `l=2` DtN fingerprint and deformation algebra, and the reduced natural point-particle source-map branch already isolated by the 2.5PN package.**

This stage does **not** introduce a new constitutive law.
It upgrades Stage 245 from the canonical outgoing DtN fixing `\chi_Q=1` to the exact **source-map reduction** of the last isotropic retarded normalization defect.

---

## Purpose

Stage 244 froze the exact isotropic grouped-real `P2` conservative one-pole carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
\frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
\]
and Stage 245 then inserted the exact compact passive/outgoing spherical `l=2` DtN branch and fixed the canonical outgoing scalar
\[
\chi_Q=1.
\]

But one reduced factor still sits between that retarded branch and the actual observable point-particle normalization:
\[
 m_{\hat 0}.
\]
In the Packet-A language of Stage 242, the remaining observable normalization defect is
\[
\Delta_{\rm norm}
:=
 m_{\hat 0}^{\,2}\,\bar P_0-
 \frac{54Gc_s^5}{5a^5c^5}.
\]
So the next honest theorem gate is:

> once the isotropic grouped-real `P2` branch is already conservative one-pole and the outgoing `l=2` fingerprint has fixed `\chi_Q`, how does the last observable normalization condition factor through the source map?

This stage answers that exactly.

The main outputs are:

1. the exact isotropic retarded branch ratio
   \[
   \frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}=\chi_Q N_Q,
   \]
2. the exact observable odd-closure factorization
   \[
   \boxed{m_{\hat 0}^{\,2}\chi_Q N_Q=1,}
   \]
3. the exact collapse of the Packet-A normalization residual to a **purely outgoing** scalar once the odd observable condition is imposed,
4. the natural source-map reduction
   \[
   m_{\hat 0}=1+O(a^2/r^2)
   \quad\Longrightarrow\quad
   N_Q=\frac{1}{\chi_Q}
   \]
   in the strict point-particle limit,
5. and the explicit source-map-reduced defect formulas obtained by inserting the Stage 245 DtN deformation algebra.

So Stage 246 is the natural source-map successor to Stage 245.

---

## 1. Carry-forward isotropic retarded invariant tuple

Define the exact Packet-A target scale
\[
\boxed{
P_0^{\rm target}:=\frac{54Gc_s^5}{5a^5c^5}.
}
\]
Using the Stage 245 carrier,
\[
\widehat Y_Q^{\rm ret}(\omega)
=
\frac34+
\frac14\,
\frac{1}{1-\omega^2/\Omega_Q^2-i\chi_Q\sigma_Q^{\rm can}\omega^5}
+O(\omega^6),
\]
with
\[
\Omega_Q=\frac{3c_s}{2a},
\qquad
\sigma_Q^{\rm can}=\frac{9}{8\Omega_Q^5}=\frac{4a^5}{27c_s^5},
\]
the isotropic retarded invariant tuple is
\[
\boxed{\bar P_0,}
\qquad
\boxed{\bar P_2=\frac{\bar P_0}{4\Omega_Q^2},}
\qquad
\boxed{\bar P_4=\frac{\bar P_0}{4\Omega_Q^4},}
\qquad
\boxed{\bar\Gamma_5=\chi_Q\frac{9\bar P_0}{32\Omega_Q^5}=\chi_Q\frac{a^5\bar P_0}{27c_s^5}.}
\]

Define the conservative isotropic normalization ratio
\[
\boxed{
N_Q:=\frac{\bar P_0}{P_0^{\rm target}}.
}
\]
Then the exact even ratios are
\[
\frac{\bar P_2}{\bar P_2^{\rm target}}=N_Q,
\qquad
\frac{\bar P_4}{\bar P_4^{\rm target}}=N_Q,
\]
while the odd ratio is
\[
\boxed{
\frac{\bar\Gamma_5}{\bar\Gamma_5^{\rm target}}=\chi_Q N_Q,
}
\]
where
\[
\boxed{
\bar\Gamma_5^{\rm target}:=\frac{2G}{5c^5}.
}
\]
So Stage 245's outgoing scalar `\chi_Q` multiplies the same isotropic conservative normalization ratio `N_Q` that already carries the even packet.

---

## 2. Exact factorization of the observable odd closure

The observable point-particle odd normalization condition is
\[
\boxed{
m_{\hat 0}^{\,2}\,\bar\Gamma_5=\bar\Gamma_5^{\rm target}.
}
\]
Substituting the ratio above gives the exact factorized condition
\[
\boxed{
m_{\hat 0}^{\,2}\chi_Q N_Q=1.
}
\]
Equivalently,
\[
\boxed{
N_Q=\frac{1}{m_{\hat 0}^{\,2}\chi_Q},
\qquad
m_{\hat 0}^{\,2}N_Q=\frac{1}{\chi_Q}.
}
\]

This is the exact source/outgoing/conservative factorization of the last isotropic retarded closure condition.
No further reduced ambiguity survives beyond these three scalars.

---

## 3. Exact Packet-A normalization collapse

Stage 242 defined the observable Packet-A normalization residual
\[
\boxed{
\Delta_{\rm norm}
:=
 m_{\hat 0}^{\,2}\bar P_0-P_0^{\rm target}
=
P_0^{\rm target}(m_{\hat 0}^{\,2}N_Q-1).
}
\]
Now insert the exact odd-closure factorization
\[
 m_{\hat 0}^{\,2}N_Q=\frac{1}{\chi_Q}.
\]
Then the Packet-A residual collapses **exactly** to
\[
\boxed{
\Delta_{\rm norm}
=
P_0^{\rm target}\left(\frac{1}{\chi_Q}-1\right).
}
\]
So once the observable odd normalization is imposed, the entire Packet-A normalization defect is already a purely outgoing scalar.
The source-map factor drops out identically.

Define the outgoing-normalization defect
\[
\boxed{\Delta_Q:=\chi_Q-1.}
\]
Then
\[
\boxed{
\Delta_{\rm norm}
=
-\,P_0^{\rm target}\,\frac{\Delta_Q}{1+\Delta_Q}.
}
\]
So the exact Packet-A normalization defect is equivalent to the exact outgoing-normalization defect.
No additional source-map ambiguity survives at this stage once the odd observable closure is enforced.

---

## 4. Natural source-map reduction

The carried natural orbital/worldtube STF source-map branch is
\[
\boxed{
m_{\hat 0}=1+O(a^2/r^2).
}
\]
So in the strict point-particle limit,
\[
 m_{\hat 0}\to 1.
\]
Then the exact factorization reduces to
\[
\boxed{
N_Q=\frac{1}{\chi_Q}.
}
\]
Therefore:
\[
\boxed{
N_Q=1
\iff
\chi_Q=1.
}
\]

At the level of the conservative carrier itself, the last isotropic retarded mismatch is therefore purely the outgoing-normalization scalar.
In particular,
\[
N_Q-1
=
\frac{1}{\chi_Q}-1
=
-\frac{\chi_Q-1}{\chi_Q}.
\]
For a small outgoing deviation,
\[
\boxed{
N_Q-1=-(\chi_Q-1)+O\!\big((\chi_Q-1)^2\big).
}
\]
So on the natural source-map branch the last reduced isotropic theorem gap is linearly controlled by the first outgoing-normalization defect.

---

## 5. Explicit source-map-reduced DtN deformation algebra

Stage 245 already showed that the most general first isotropic DtN deformation preserving the canonical even moments is
\[
\Lambda_2^{\rm def}(z)
=
S\Lambda_2^{\rm out}(\beta z)
+\Sigma_0+\Sigma_2 z^2+\Sigma_4 z^4+i\Sigma_5 z^5+O(z^6),
\]
with canonical-even constraints fixing
\[
\Sigma_2=-\frac{3S\beta^2-3S+\Sigma_0}{9},
\qquad
\Sigma_4=-\frac{3S\beta^4-3S+\Sigma_0}{27},
\]
and outgoing normalization
\[
\boxed{
\chi_Q=
\frac{3(S\beta^5+9\Sigma_5)}{3S-\Sigma_0}.
}
\]

### 5.1 Exact conservative ratio on the natural source-map branch

Insert this into the Stage 246 source-map reduction:
\[
N_Q=\frac{1}{\chi_Q}.
\]
Then
\[
\boxed{
N_Q
=
\frac{3S-\Sigma_0}{3(S\beta^5+9\Sigma_5)}.
}
\]
So the source-map-reduced conservative normalization ratio is already an explicit function of the same isotropic DtN deformation data that control `\chi_Q`.

### 5.2 Exact Packet-A normalization defect on the natural source-map branch

Likewise,
\[
\Delta_{\rm norm}^{\rm pt}
:=
\Delta_{\rm norm}\big|_{m_{\hat 0}\to 1}
=
P_0^{\rm target}(N_Q-1),
\]
so
\[
\boxed{
\Delta_{\rm norm}^{\rm pt}
=
P_0^{\rm target}
\left[
\frac{3S-\Sigma_0}{3(S\beta^5+9\Sigma_5)}-1
\right].
}
\]
Equivalently,
\[
\boxed{
\Delta_{\rm norm}^{\rm pt}
=
-\,P_0^{\rm target}
\frac{3S(\beta^5-1)+\Sigma_0+27\Sigma_5}
     {3(S\beta^5+9\Sigma_5)}.
}
\]
So the entire source-map-reduced Packet-A defect is already localized on the same three isotropic DtN deformation coordinates isolated in Stage 245:
\[
(\beta,\Sigma_0,\Sigma_5).
\]

### 5.3 Linearized source-map-reduced defect

Write
\[
\beta=1+\varepsilon_\beta,
\qquad
\Sigma_0=\delta\Sigma_0,
\qquad
\Sigma_5=\delta\Sigma_5,
\]
with all three small. Then
\[
\boxed{
N_Q-1
=
-5\varepsilon_\beta
-\frac{\delta\Sigma_0}{3S}
-\frac{9\,\delta\Sigma_5}{S}
+O(2).
}
\]
Hence
\[
\boxed{
\Delta_{\rm norm}^{\rm pt}
=
-\,P_0^{\rm target}
\left(
5\varepsilon_\beta
+\frac{\delta\Sigma_0}{3S}
+\frac{9\,\delta\Sigma_5}{S}
\right)
+O(2).
}
\]
So at first order the source-map-reduced normalization defect is carried only by:

1. argument deformation `\beta`,
2. static isotropic core shift `\Sigma_0`,
3. odd isotropic core outlet `\Sigma_5`.

This is the exact linearized finish-line map for the Stage 245 deformation algebra.

---

## 6. Canonical compact outgoing branch

On the canonical compact passive/outgoing branch from Stage 245,
\[
\chi_Q=1.
\]
Equivalently, in the first isotropic DtN deformation family,
\[
\beta=1,
\qquad
\Sigma_0=0,
\qquad
\Sigma_5=0.
\]
Then the Stage 246 formulas reduce exactly to
\[
\boxed{N_Q=1,}
\qquad
\boxed{\Delta_{\rm norm}=0,}
\qquad
\boxed{\Delta_{\rm norm}^{\rm pt}=0.}
\]
So once the canonical outgoing DtN fingerprint of Stage 245 is carried through the exact source-map reduction, the Packet-A normalization residual of Stage 242 closes automatically.

This is the exact bridge from:

- Stage 242's observable branch residual `\Delta_{\rm norm}`,
- through Stage 244's isotropic one-pole conservative carrier,
- through Stage 245's fixing `\chi_Q=1`,
- to the strict point-particle natural source-map branch.

---

## 7. What Stage 246 changes in the theorem problem

Stage 245 left the hierarchy with one exact outgoing scalar `\chi_Q`, but the observable point-particle branch still carried the source-map factor `m_{\hat 0}`.
Stage 246 changes that in four precise ways.

### 7.1 The observable odd closure is now fully factorized

The full odd condition is exactly
\[
 m_{\hat 0}^{\,2}\chi_Q N_Q=1.
\]
So the last observable isotropic retarded closure is one exact three-factor product.

### 7.2 The Packet-A normalization residual collapses to a purely outgoing scalar

Once the observable odd closure is imposed,
\[
\Delta_{\rm norm}=P_0^{\rm target}(\chi_Q^{-1}-1).
\]
So the source-map factor no longer appears in the normalization residual itself.

### 7.3 On the natural source-map branch the conservative carrier ratio is also purely outgoing

In the strict point-particle limit,
\[
N_Q=\chi_Q^{-1}.
\]
So the last reduced isotropic branch defect is no longer a mixed source/outgoing product.
It is purely the outgoing scalar `\chi_Q-1`.

### 7.4 The remaining isotropic PDE-facing freedom is now completely explicit

Combining Stage 245 with Stage 246 shows that the full source-map-reduced Packet-A defect is carried only by
\[
(\beta,\Sigma_0,\Sigma_5),
\]
through the exact formulas above.
So the remaining isotropic PDE realization problem is no longer diffuse.

---

## 8. Immediate next derivation step

The next clean continuation is now very sharp.

Do **not** reopen the source-map branch or the conservative grouped bundle.
Instead:

1. keep the Stage 245 exact outgoing `l=2` fingerprint,
2. keep the Stage 246 exact factorization
   \[
   m_{\hat 0}^{\,2}\chi_Q N_Q=1,
   \]
3. and prove that any extra retarded structure first entering at `O(\omega^7)` or higher is irrelevant to the 2.5PN theorem.

That is the natural next theorem gate after Stage 246.
