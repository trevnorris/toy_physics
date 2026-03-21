# 2PN axisymmetric Robin-wall PDE scaffold: current result

## 1. What this step adds

This step replaces the abstract channel-by-channel static fit with the first genuinely **explicit throat-wall PDE scaffold**.

The idea is simple:

- keep the 4D-ball DtN kernel as the bulk-side unit-test PDE anchor,
- add an **axisymmetric wall law** at the mouth,
- and ask what the smallest \(P_2\)-structured wall operator must look like in order to reproduce the already-solved raw static support data.

The answer is pleasantly sharp:

1. a finite monopole support comes from an **isotropic monopole wall stiffness**,
2. the odd dipole sector is reproduced exactly by a **first-order \(P_2\)** wall perturbation,
3. the even quadrupole sector **cannot** be reproduced by first-order \(P_2\) alone,
4. but it **is** reproduced exactly by adding one **second-order \(P_2^2\)** wall term.

So the conservative 2PN cross operator now has a concrete PDE-side wall-law interpretation.

---

## 2. Raw axisymmetric channel moments

Using the raw mouth basis from the previous step and the axisymmetric quadrupole
\[
P_2(\cos\theta)=\frac{3z^2-1}{2},
\]
the channel expectation values are

### Dipole sector
\[
\langle P_2\rangle_{1\perp}=-\frac15,\qquad
\langle P_2\rangle_{10}=\frac25,
\]
\[
\langle P_2^2\rangle_{1\perp}=\frac17,\qquad
\langle P_2^2\rangle_{10}=\frac{11}{35}.
\]

### Quadrupole sector
\[
\langle P_2\rangle_{20}=\frac27,\qquad
\langle P_2\rangle_{21}=\frac17,\qquad
\langle P_2\rangle_{22}=-\frac27,
\]
\[
\langle P_2^2\rangle_{20}=\frac37,\qquad
\langle P_2^2\rangle_{21}=\frac17,\qquad
\langle P_2^2\rangle_{22}=\frac17.
\]

These are exactly the axisymmetric splitting factors that appear in the wall-law fit.

---

## 3. Static wall-law fit

The already-solved raw support residues are
\[
Y_{0}= \frac{45}{4},\qquad
Y_{1\perp}=\frac72,\qquad
Y_{10}=4,\qquad
Y_{20}=\frac94,\qquad
Y_{21}=\frac32,\qquad
Y_{22}=\frac38,
\]
so the corresponding static impedances are
\[
Z_{0}=\frac4{45},\qquad
Z_{1\perp}=\frac27,\qquad
Z_{10}=\frac14,\qquad
Z_{20}=\frac49,\qquad
Z_{21}=\frac23,\qquad
Z_{22}=\frac83.
\]

### 3.1 Monopole
A finite monopole support requires
\[
Z_0(0)=\frac4{45}.
\]
This is the first explicit wall-side realization of the previously abstract “finite \(\Omega_0\)” requirement.

### 3.2 Odd dipole sector
The exact odd fit is already achieved by a first-order axisymmetric wall law
\[
Z_{1m}(0)=B_1+A_1\langle P_2\rangle_{1m},
\]
with
\[
B_1=\frac{23}{84},\qquad A_1=-\frac{5}{84}.
\]

This gives
\[
Z_{1\perp}=\frac27,\qquad Z_{10}=\frac14
\]
exactly.

### 3.3 Even quadrupole sector: first-order no-go
Trying the same first-order structure
\[
Z_{2m}(0)=B_2+A_2\langle P_2\rangle_{2m}
\]
fails.

Matching the \(m=0\) and \(|m|=1\) channels gives
\[
B_2=\frac89,\qquad A_2=-\frac{14}{9},
\]
which predicts
\[
Z_{22}^{\rm first\ order}=\frac43,
\]
but the solved target is
\[
Z_{22}^{\rm target}=\frac83.
\]

So a linear \(P_2\) wall is **not enough**.

### 3.4 Minimal exact quadrupole wall fit
The minimal exact extension is
\[
Z_{2m}(0)=B_2+A_2\langle P_2\rangle_{2m}+C_2\langle P_2^2\rangle_{2m},
\]
with
\[
B_2=\frac{10}{9},\qquad
A_2=-\frac{14}{3},\qquad
C_2=\frac{14}{9}.
\]

This reproduces
\[
Z_{20}=\frac49,\qquad
Z_{21}=\frac23,\qquad
Z_{22}=\frac83
\]
exactly.

So the minimal static wall law is

\[
Z_0(0)=\frac4{45},
\]
\[
Z_{1m}(0)=\frac{23}{84}-\frac{5}{84}\,\langle P_2\rangle_{1m},
\]
\[
Z_{2m}(0)=\frac{10}{9}-\frac{14}{3}\,\langle P_2\rangle_{2m}+\frac{14}{9}\,\langle P_2^2\rangle_{2m}.
\]

---

## 4. Axisymmetric source profile

The previously solved source vector is also natural in this wall language.

A two-parameter axisymmetric source
\[
S(\theta)=p_{\rm iso}+q_{\rm ax}P_2(\cos\theta)
\]
reproduces the required normalized source vector when
\[
p_{\rm iso}=\frac{11}{8},\qquad q_{\rm ax}=\frac{15}{8}.
\]

So both the support operator **and** the source sector now have an explicit axisymmetric wall-level realization.

---

## 5. Dynamic 4D-ball DtN pole equations

Keeping the 4D-ball unit-test DtN kernel
\[
\Lambda_\ell(z)= -1 + z\,\frac{J'_{\ell+1}(z)}{J_{\ell+1}(z)},
\qquad z\equiv \frac{a\Omega}{c_s},
\]
and promoting the static wall impedances directly into the pole equation
\[
\Lambda_\ell(z)+Z_{\ell m}(0)=0,
\]
gives the following lowest positive roots.

### Monopole
\[
z_0=0.591884444464394.
\]
The small-\(z\) estimate from \(\Lambda_0(z)\simeq -z^2/4\) is
\[
z_0 \simeq \frac{4}{\sqrt{45}}=0.596284793999944,
\]
so the finite monopole pole is already visible at leading order.

### Dipole
Using the isotropic odd base \(B_1=23/84\),
\[
z_1^{\rm base}=2.551215916564765,
\]
while the exact split poles are
\[
z_{1\perp}=2.561183722397930,\qquad
z_{10}=2.531063390840353.
\]

The first-order perturbative shifts about the isotropic base are
\[
z_{1\perp}^{\rm pert}=2.561219561244114,\qquad
z_{10}^{\rm pert}=2.531208627206067,
\]
which are extremely accurate. So the odd dipole sector is genuinely perturbative.

### Quadrupole
Using the isotropic even base \(B_2=10/9\),
\[
z_2^{\rm base}=4.254105628646177,
\]
while the exact split poles are
\[
z_{20}=3.901921523190568,\qquad
z_{21}=4.029116369391941,\qquad
z_{22}=4.821811915561263.
\]

First-order shifts about that base give
\[
z_{20}^{\rm pert}=3.942783453176964,\qquad
z_{21}^{\rm pert}=4.046557511666702,\qquad
z_{22}^{\rm pert}=4.980524038074339.
\]

So the perturbative ordering is correct, but the quadrupole anisotropy is large enough that the exact channel equations should be treated nonperturbatively.

---

## 6. Why this matters

This is the first point where the inner-throat program stops being merely “a target channel algebra” and becomes an actual PDE-side construction.

We now know:

1. **finite \(\Omega_0\)** comes from a finite isotropic monopole wall stiffness;
2. **dipole splitting** comes from a first-order axisymmetric \(P_2\) wall term;
3. the solved **quadrupole support hierarchy** requires a minimal second-order \(P_2^2\) wall contribution;
4. the source vector is exactly a two-parameter axisymmetric source;
5. and the full low-frequency channel structure can be attached directly to the 4D-ball DtN pole equations.

So the next derivation target is no longer vague. It is:

> derive these wall coefficients from a concrete soft-wall throat potential or stress-balance wall law, instead of fitting them at the operator level.

That is a very manageable next step, and it is exactly the wall-law gap the throat notes said had to be frozen before full throat modeling becomes coherent.
