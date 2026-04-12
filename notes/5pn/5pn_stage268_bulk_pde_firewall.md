# 5PN / Moving-Throat Continuation — Stage 268: Bulk-PDE Firewall and Boundary-Response Reading

## Purpose

This note records a clarification that matters for the rest of the 5PN program:

> we should **not** try to repair the 5PN gap by changing the parent bulk GNLS equation.

The question is whether that statement is merely intuitive or whether it is actually aligned with the present derivation chain.

The verdict is:

- **yes**, the current program supports a strong **bulk-PDE firewall** reading;
- but the strongest honest theorem-level wording is **not** “5PN divergence is purely a boundary condition” in an absolute sense;
- it is the more careful statement that **within the present reduced hierarchy**, the remaining 5PN / 2.5PN / 4PN gap has been compressed to the **moving-throat branch data**, i.e. the conservative grouped-`P2` carrier, the outgoing quadrupole normalization, and the weak-axisymmetric orbit/selector defects, **not** a retuning of the parent GNLS medium.

---

## 1. What is already frozen and therefore should not be changed

The parent `4+1` bulk theory is already fixed as a gauged GNLS matter sector plus localized Maxwell sector, with the stiff-polytropic EOS
\[
P(\rho)=K\rho^5,
\]
so the exponent `n=5` is already part of the carried exact parent setup. The exact parent fields remain
\[
\psi(\mathbf X,t),\qquad A_M,\qquad \Sigma=r-R(\Omega,w,t),
\]
and the compact moving-throat program explicitly treats this as the correct parent-theory starting point rather than a tunable late-stage ansatz.

At the lower PN levels, the hierarchy has already frozen the `n=5` calibration and the conservative 1PN–4PN carry-forward ledgers. In particular:

- lower-order optical / 1PN matching fixed the EOS exponent to `n=5`;
- the 4PN summary says the full local 4PN sector is already closed within the declared hierarchy;
- and the only remaining full-4PN gap is **the same passive/outgoing quadrupole normalization** already isolated by the 2.5PN program.

So if we were to modify the parent GNLS medium now, we would no longer be continuing the already-calibrated hierarchy. We would be starting a different hierarchy.

---

## 2. What the actual 5PN reduced gap now depends on

The 5PN notes already say the explicit Family-1 support/source side is **not** the active bottleneck anymore. The live question is whether the actual grouped-`P2` / geometry branch realizes the minimal isotropic conservative quadrupole module, and then whether the passive/outgoing `l=2` branch carries the canonical outgoing normalization.

The same notes also show that on the natural isotropic branch:

- the grouped conservative carrier is the `3/4 + 1/4` module,
  \[
  \widehat Y_Q^{\rm cons}(\omega)
  =
  \frac34+\frac14\frac{1}{1-\omega^2/\Omega_Q^2},
  \]
  so
  \[
  c_{\rm pole}=\frac14;
  \]
- the scalar/geometry `l=0` lane is orthogonal to the grouped real `l=2` bundle at linear isotropic order, so
  \[
  \epsilon_2=\epsilon_4=0
  \]
  on that branch;
- and the first nonzero geometry contamination requires explicit `l=0 \leftrightarrow l=2` mixing and enters only at
  \[
  O(\chi^2).
  \]

So the live 5PN ambiguity is already localized to the **moving-throat response branch**, not to the bulk EOS or bulk GNLS law.

---

## 3. Packet-level statement of the firewall

Stages 199–201 compressed the reduced closure problem to two finite data packets:

### Packet A — grouped branch packet
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
plus the source-map factor `mhat_0`.

### Packet B — orbit/selector packet
any equivalent form of
\[
(m_T,m_K,m_\mu),
\qquad
(R_{\rm tr},R_{\rm nt},R_\eta),
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Everything in the reduced 5PN / 2.5PN / 4PN endgame is an exact compiler output of those packets.

That means the remaining task is:

- compute the actual packet values from the moving-throat branch;
- test isotropy / `c_{\rm pole}=1/4` / outgoing normalization / weak-axisymmetric orbit lock;
- and **do not** reopen the parent GNLS action unless there is a direct contradiction with the already frozen lower-order hierarchy.

This is the precise sense in which the parent bulk PDE is “safe.”

---

## 4. The right way to phrase the user’s proposed update

The proposed reading is very close to the current theorem picture, but it needs one wording adjustment.

### Good and useful
It is correct and helpful to say:

1. the lower-order program already calibrated the background 4D medium;
2. the remaining 5PN problem is in the throat / response / boundary realization;
3. the active moving pieces are the restoring pole fraction and the first `l=0 \leftrightarrow l=2` mixing channel, not a new bulk EOS fit.

### Needs tightening
The strongest honest program-level wording is **not**
> “the 5PN divergence is purely a boundary condition”

because the compact moving-throat ledger is explicit that the remaining gap is a **realization** gap: whether the actual completed PDE realizes the reduced outgoing, mouth/core, coherent-support, and invariant structures strongly enough. So the right wording is:

> **Within the present reduced hierarchy, the unresolved 5PN data enter through the moving-throat branch / boundary-response sector, not through a retuning of the parent GNLS bulk medium.**

That is narrower, cleaner, and actually matches the current notes.

---

## 5. What to do with the black-hole / sink intuition

The “microscopic particle = subsonic / resonant interior survives” versus
“black hole = supersonic sink / resonance shredded” picture is a **reasonable physical intuition** for future branch classes, but it is **not yet** one of the program’s frozen theorem inputs.

So it should be used as:

- a motivation for why different compact objects might correspond to different throat boundary classes;

not yet as:

- a theorem we have already derived from the moving-throat PDE.

For the current 5PN derivation, the safe formal statement is simply:

- keep the exact parent bulk PDE fixed;
- classify the remaining freedom in the throat response/boundary branch.

---

## 6. Practical consequence for the next derivation move

The next derivation should therefore keep the bulk firewall explicit and continue only on the boundary/branch side:

1. keep the parent GNLS + localized Maxwell action frozen;
2. keep `n=5` and the 1PN–4PN carry-forward constants frozen;
3. treat the grouped `P2` carrier and outgoing normalization as branch data to be extracted, not bulk data to be refit;
4. continue deriving the actual grouped-`P2` / weak-axisymmetric branch formulas under that firewall.

This is the clean continuation point.

---

## 7. Recommended replacement text

A tighter version of the proposed context update is:

> **Bulk-PDE firewall for the 5PN moving-throat program.**  
> We keep the parent `4+1` GNLS + localized Maxwell bulk theory fixed, including the already-calibrated `P=K\rho^5` medium. The lower-order 1PN–4PN hierarchy has already frozen that background. Within the present reduced hierarchy, the remaining 5PN / 2.5PN / 4PN gap is localized to the moving-throat branch data: the isotropic grouped-`P2` conservative carrier, the passive/outgoing `l=2` normalization, and the weak-axisymmetric orbit/selector defects. On the natural isotropic branch this means `c_pole=1/4`, while the first geometry contamination enters only through explicit `l=0 \leftrightarrow l=2` mixing at `O(\chi^2)`. So the correct next move is to keep the bulk PDE fixed and continue extracting the throat boundary/response branch, rather than altering the parent GNLS equation.

That is the version I would carry forward.

