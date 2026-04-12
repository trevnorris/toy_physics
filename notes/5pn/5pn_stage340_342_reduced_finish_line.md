
# 5PN continuation notes — Stages 340–342: reduced finish-line compression

This session does not introduce a new microscopic branch family. It compresses the
current theorem ledger to the narrowest reduced finish line reached so far.

The guiding idea is simple:

- the Stage-200 endgame compiler still carries an 8-component reduced branch packet
  together with the 3-component orbit packet;
- the later moving-throat branch work shows that, on the **actual isotropic**
  static-geometry one-pole contact-plus-pole branch, almost all of that structure
  has already collapsed.

So the right final move is to prove that collapse cleanly and to state the exact
remaining theorem gate in the smallest possible language.

---

## Stage 340 — the Stage-200 branch packet collapses to one scalar on the actual isotropic branch

The exact Stage-200 reduced branch packet is

\[
\Delta_{\rm branch}
=
(a_2,\ b_2,\ a_4,\ b_4,\ a(P_0),\ b(P_0),\ \Delta_{\rm pole},\ \Delta_{\rm norm}).
\]

On the later actual isotropic one-pole branch, the following hold simultaneously:

1. grouped-lane isotropy kills
   \[
   a_2=b_2=a_4=b_4=a(P_0)=b(P_0)=0,
   \]
2. the static-geometry one-pole conservative branch kills
   \[
   \Delta_{\rm pole}=0,
   \]
3. the only surviving reduced branch mismatch is the outgoing-normalization scalar
   \[
   N_Q-1.
   \]

So the full 8-component reduced branch packet collapses exactly to

\[
\Delta_{\rm branch}
\longrightarrow
(0,0,0,0,0,0,0,N_Q-1).
\]

Equivalently,
\[
\Delta_{\rm branch}=0
\iff
N_Q=1
\]
once the actual isotropic one-pole branch is accepted.

That is the cleanest bridge from the Stage-200 general endgame to the later
moving-throat finish-line statement.

---

## Stage 341 — the orbit packet is exactly a three-observable invariance test

The weak-axisymmetric orbit-lock side can now be written directly in the three
physical branch observables

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta.
\]

At first order the defect triplet is

\[
\Theta_1 = \delta\ln R_{\rm tr},
\]
\[
\Xi_1 = -\,\delta\ln R_{\rm target}
        -\frac{\epsilon_\eta}{1-\epsilon_\eta}\,\delta\ln \epsilon_\eta,
\]
\[
\mathcal R_1 = \delta\ln R_{\rm target}.
\]

And the inverse map is exact:

\[
\delta\ln R_{\rm tr} = \Theta_1,
\]
\[
\delta\ln R_{\rm target} = \mathcal R_1,
\]
\[
\delta\ln \epsilon_\eta
=
-\frac{1-\epsilon_\eta}{\epsilon_\eta}\,(\mathcal R_1+\Xi_1).
\]

So the full weak-axisymmetric orbit-lock problem is now exactly equivalent to the
three-observable invariance test

\[
\Theta_1=\Xi_1=\mathcal R_1=0
\iff
\delta\ln R_{\rm tr}=0,\quad
\delta\ln R_{\rm target}=0,\quad
\delta\ln\epsilon_\eta=0.
\]

This is the direct-observable form of the older similarity-orbit / monomial-rigidity
theorem.

---

## Stage 342 — reduced completion theorem

There are now three exact pieces to combine.

### 1. Support/source theorem on the actual isotropic branch

On the actual isotropic branch,
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}^{(\rm act)}(\epsilon_{\rm blk})
=
\frac{1}{3-2\epsilon_{\rm blk}}.
\]

If a constructive explicit support/source family has ceiling `zeta_max > 1`, then
through the full admissible blocked window
\[
0\le \epsilon_{\rm blk}<\frac1{\zeta_{\max}},
\]
one has
\[
\zeta_{\rm req}^{(\rm act)} < \zeta_{\max}.
\]

Indeed, at the worst admissible edge
\[
\zeta_{\rm req,worst}
=
\frac{\zeta_{\max}}{3\zeta_{\max}-2},
\]
and
\[
\zeta_{\max}-\zeta_{\rm req,worst}
=
\frac{3\zeta_{\max}(\zeta_{\max}-1)}{3\zeta_{\max}-2}>0
\quad\text{for}\quad
\zeta_{\max}>1.
\]

So the explicit support/source side is automatic on the actual isotropic branch.
It is no longer an active reduced bottleneck.

### 2. Outgoing normalization on the natural source-map branch

On the natural source-map branch,
\[
\hat m_0\to 1,
\]
so the exact outgoing factorization reduces to
\[
N_Q=\frac{1}{\chi_Q}.
\]

Therefore
\[
N_Q=1
\iff
\chi_Q=1.
\]

So on the natural source-map branch the last outgoing defect is purely the scalar
outgoing renormalization `chi_Q - 1`.

### 3. Final reduced finish line

Inside the present reduced hierarchy, once the actual isotropic conservative precursor
is accepted, the full reduced 5PN / 2.5PN / 4PN closure problem has collapsed to

\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0,
\qquad
N_Q=1,
\]
or, on the natural source-map branch,
\[
\delta\ln R_{\rm tr}=0,
\qquad
\delta\ln R_{\rm target}=0,
\qquad
\delta\ln\epsilon_\eta=0,
\qquad
\chi_Q=1.
\]

So the reduced theorem gate is no longer:

- a large grouped branch packet,
- a separate support/source branch phase,
- and a diffuse outgoing normalization problem.

It is now exactly:

1. a **three-observable orbit-lock test**, and
2. a **one-scalar outgoing-normalization test**.

That is the narrowest honest reduced finish line reached so far.

---

## Best current verdict

Within the present reduced hierarchy:

- the explicit Family-1 support/source side is no longer an independent obstacle,
- the Stage-200 branch packet has collapsed to `N_Q - 1` on the actual isotropic branch,
- the weak-axisymmetric side is exactly the invariance of
  \[
  (R_{\rm tr},R_{\rm target},\epsilon_\eta),
  \]
- and the only remaining reduced theorem gap is whether the completed moving-throat PDE
  actually realizes

\[
\delta\ln R_{\rm tr}
=
\delta\ln R_{\rm target}
=
\delta\ln\epsilon_\eta
=
0,
\qquad
N_Q=1
\]
(or equivalently `chi_Q = 1` on the natural source-map branch).

That is as close to “finished” as the reduced program can honestly get without the
final PDE branch realization.
