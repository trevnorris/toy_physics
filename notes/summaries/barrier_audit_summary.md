
barrier_audit_summary.md — Comprehensive summary (same-charge barrier audit / mixed-bundle kill-test package)

## 0) What this package is doing

This document is the **carry-forward same-charge barrier summary** for the moving-throat / 4D toy-model program.

It is built from `barrier_audit_full.md`, but it is written in the same role as the other summary files: not as a paper draft, not as raw session notes, and not as a stage-by-stage diary. Its job is to freeze the actual theorem structure that now survives after the long same-charge audit.

The package compresses the whole chain into one readable ledger:

1. replace the original generic same-charge mixed placeholder by the **first exact one-port mixed-bundle kernel**;
2. show that the minimal static mixed sector produces only short-range product families, not a new long-range attractive law;
3. show that linear dynamic driving and the first passive/outgoing correction do **not** create conservative barrier lowering at first order;
4. derive the resonance/linewidth tradeoff and then test an explicit primitive finite-throat branch;
5. impose the exact isotropic `5`PN / `2.5`PN target surface and transport the resulting windows onto the actual branch packet;
6. run the weak-axisymmetric mechanism sieve and identify the only surviving first-order corridor;
7. compare that corridor to the stricter imported `5`PN even-gate package and isolate the pure-transfer subcorridor;
8. pull the selected-branch classifier back to the continuum placement variables and show that the first kill condition remains **static-first** rather than dynamic;
9. inject the already-known `5`PN support/source data and show that support enhancement is no longer the bottleneck;
10. factor the actual rigid-mouth branch into transfer-shape and dressing legs and isolate the post-static orbit-lock test;
11. fix the selected support slice by matching to the minimal isotropic passive/outgoing quadrupole precursor; and
12. end with the exact actual-branch compiler packet the completed moving-throat PDE still has to return.

The strongest honest carry-forward statement is:

> **The same-charge corridor is not killed by naive Maxwell shaping, by the first honest static mixed bundle, by linear outgoing-phase response, by naive resonance, by wall-only or BdG-only weak anisotropy, or by the currently known support/source data. What survives is a narrow mixed-sector / pure-transfer / selected-twin-support corridor. The first unresolved theorem gate is now the actual PDE-selected coherent orbit packet together with the canonical outgoing-normalization finish line, not broad coefficient freedom.**

So this package should be read as the barrier-audit analog of the `2.5`PN / `4`PN / `g-2` summaries: the algebraic compression is largely done, but the final physical branch verdict still sits in one sharply localized moving-throat packet.

---

## 1) Claim taxonomy and how to read this package

### 1.1 Exact

These are exact symbolic or algebraic statements once the carried one-port bundle, selected-branch definitions, and grouped-lane dictionaries are fixed.

Examples:

- the static reduced bundle determinant
  \[
  \det \mathcal K_{\rm red} = \Delta D_0;
  \]
- the exact static susceptibility law
  \[
  \delta V_{\rm mix}(r) = -\frac12 J(r)^T \mathcal K_{\rm red}^{-1} J(r);
  \]
- the static product-family theorem;
- the dynamic one-port susceptibility law and the outgoing derivative identity;
- the linear phase-lag no-go theorem;
- the exact `5`PN transport compiler from the actual grouped branch packet;
- the weak-axisymmetric grouped law
  \[
  P_A = \bar P_0(1+\epsilon \lambda_A \Xi_1);
  \]
- the selected-branch numerator/denominator classifier;
- the rigid-mouth physical normal form in
  \[
  (U,V)=\left(\ln \frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},\ \ln \frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}\right);
  \]
- the selected support ratio
  \[
  \rho_\alpha = \frac43,\qquad \zeta_{\rm req}=\frac13,\qquad \Pi_{\rm tr}=\frac43 C_{\rm mix}.
  \]

### 1.2 Exact within closure

These are exact once a stated reduced branch or hierarchy is adopted.

Examples:

- the isotropic one-port wall/BdG/Maxwell/mixed bundle;
- the explicit primitive finite-throat branch used in the pole census;
- the exact isotropic `5`PN compatibility surface;
- the weak-axisymmetric grouped reduction;
- the rigid-mouth actual branch;
- the coherent local D/N support/source branch;
- the selected twin-support branch.

So these are not heuristic guesses, but they are branch-level exact statements rather than theorems of the full unsolved two-throat PDE.

### 1.3 Reduced / controlled reduction

These are claims that rely on the same style of reduced hierarchy already used in the moving-throat PN program.

Examples:

- the one-port static and dynamic barrier audits;
- the pole/linewidth reduction;
- the grouped prefactor transport law;
- the selected-branch pullback from continuum placement data;
- the rigid-mouth packet factorization.

The package is not an unreduced theorem of the full same-charge interaction PDE.

### 1.4 Conditional actual-branch verdict

The final survival statement is still conditional.

The reduced stack now says:

- support/source enhancement is already strong enough;
- dynamic wall-window survival is not the first bottleneck;
- the remaining unresolved kill test is the actual coherent placement / orbit packet together with outgoing normalization.

So the final answer depends on the actual moving-throat branch data, not on another broad unresolved coefficient family.

---

## 2) What this package claims — and what it does not claim

### 2.1 What it claims

This package claims all of the following.

1. The first honest one-port mixed bundle generates an exact static susceptibility kernel whose denominator is the same \(\Delta D_0\) already used by the `5`PN / `2.5`PN normalization chain.
2. In the minimal static one-port setting, the mixed sector creates **no new long-range attractive family**. It produces only quadratic product kernels built from the primitive source profiles.
3. The linear dynamic one-port bundle keeps the same spatial product families and only makes their coefficients frequency dependent.
4. The first passive/outgoing correction is purely imaginary at linear order, so it produces phase lag / pumping rather than conservative barrier lowering.
5. A resonance corridor survives linearly only through a residue-to-linewidth tradeoff; naive “resonance fixes it” does not survive the audit.
6. An explicit primitive finite-throat branch can be placed on the exact isotropic `5`PN surface while still retaining a finite same-charge dynamic corridor.
7. The actual-branch transport compiler reduces the wall-window test to explicit inequalities in
   \[
   \Delta_{\rm norm},\ \hat m_0,\ \Xi_1.
   \]
8. Wall-only and pure-BdG weak-axisymmetric compensated families are killed on the concrete branch, while a constrained mixed-sector corridor survives.
9. After imposing the stricter imported `5`PN even-gate package, the surviving corridor can still be reduced further to a pure-transfer subcorridor with
   \[
   D_{01}=D_{21}=D_{41}=0,\qquad \Xi_1=\frac{N_{01}}{N_0}.
   \]
10. The selected same-charge branch is controlled by an exact numerator/denominator co-loading classifier, and the first kill condition remains **static-first**, not dynamic.
11. Pulling that classifier back to the continuum placement ratios preserves the static-first theorem.
12. Injecting the already-known `5`PN support/source data shows that support enhancement is no longer the bottleneck inside the current reduced theorem stack.
13. On the actual rigid-mouth branch, the packet factorizes completely into transfer-shape and dressing axes, both support-blind.
14. The selected support sector is fixed exactly to the symmetric lowest twin-support slice of the minimal isotropic passive/outgoing quadrupole precursor.
15. The remaining actual-branch finish line is the coherent orbit-lock / normalization packet:
   \[
   d\ln R_{\rm tr},\quad d\ln R_{\rm target},\quad d\ln \epsilon_\eta,\quad N_Q-1.
   \]

### 2.2 What it does **not** claim

This package does **not** claim:

- a full two-throat moving-throat PDE theorem for same-charge approach;
- an unconditional same-charge bound state or barrier bypass theorem;
- a new long-range attractive law produced by the mixed bundle;
- that linear outgoing response by itself can lower the barrier;
- that the actual branch must lie on the surviving reduced corridor;
- that support/source data alone decide the final verdict.

So the correct reading is:

> the reduced audit has killed the easy stories, preserved one narrow mixed-sector corridor, and pushed the real verdict all the way down to one actual branch packet.

---

## 3) Fixed inputs, notation firewall, and carried bundle objects

### 3.1 Static one-port bundle

After the stable BdG mode is integrated into the wall stiffness, the reduced one-port bundle is controlled by
\[
K_* = K-\frac{C^2}{\varpi^2},
\]
\[
\Delta = \Omega_U^2\Omega_W^2-R^2,
\]
\[
Q = G_U^2\Omega_W^2+2G_U G_W R+G_W^2\Omega_U^2,
\]
\[
P = \Omega_U^2 G_W + R G_U,
\]
\[
D_0 = K_* - \frac{Q}{\Delta}.
\]

These are the same one-port quantities that feed the isotropic normalization side:
\[
N_0=\frac{P^2}{\Delta^2},
\qquad
P_0=\frac{N_0}{D_0},
\qquad
m_{\hat 0}^{\,2}P_0=\frac{54Gc_s^5}{5a^5c^5}.
\]

So the barrier audit and the outgoing-normalization ledger are not independent.

### 3.2 Dynamic one-port bundle

The reduced dynamic bundle uses
\[
K_B(\omega)=K-M\omega^2-\frac{C^2}{\varpi^2-\omega^2},
\]
\[
A(\omega)=\Omega_U^2-\omega^2,
\qquad
W(\omega)=\Omega_W^2-\omega^2-\Pi(\omega),
\]
\[
\Delta_\Pi(\omega)=A(\omega)W(\omega)-R^2,
\]
\[
Q_\Pi(\omega)=G_U^2W(\omega)+2G_UG_WR+G_W^2A(\omega),
\]
\[
D_\Pi(\omega)=K_B(\omega)-\frac{Q_\Pi(\omega)}{\Delta_\Pi(\omega)}.
\]

The outgoing fingerprint later used in the audit is the same passive/outgoing mixed-port law already isolated in the moving-throat quadrupole program.

### 3.3 Primitive source profiles

The first primitive reduced sources are the quadrupolar/Coulomb-Hessian profile
\[
\mathcal S_Q(x)=\frac{1}{x^3},
\]
and the Yukawa/mixed profile
\[
\mathcal S_Y(x)=\frac{e^{-2\kappa x}}{x},
\qquad
\kappa=\frac{a}{\lambda}.
\]

These are the primitive spatial inputs whose product families the one-port bundle can generate.

### 3.4 Actual-branch packet from the `5`PN endgame

The actual grouped branch packet carried into the later stages is
\[
\Delta_{\rm branch}
=
(a_2,b_2,a_4,b_4,a_{P_0},b_{P_0},\Delta_{\rm pole},\Delta_{\rm norm}),
\]
with normalization slot
\[
\Delta_{\rm norm}
=
\hat m_0^{\,2}\,\bar P_0-\frac{54Gc_s^5}{5a^5c^5}.
\]

The exact mean prefactor is therefore
\[
\bar P_0
=
\frac{\Delta_{\rm norm}+\dfrac{54Gc_s^5}{5a^5c^5}}{\hat m_0^{\,2}}.
\]

### 3.5 Transported primitive-family ceilings

The first carried dynamic windows are the Stage-006 / Stage-007 primitive-family ceilings:

at the stricter `10%`-loss benchmark,
\[
P_{\rm both}^{(10)}=0.0028313316855593175,
\qquad
P_{\rm one}^{(10)}=0.0035965105896846573;
\]

at the looser `30%`-loss benchmark,
\[
P_{\rm both}^{(30)}=0.00817339430971383,
\qquad
P_{\rm one}^{(30)}=0.0116633929790174.
\]

These are still reduced branch ceilings, not full-PDE theorems, but they are the first exact windows against which the actual branch can be tested.

### 3.6 Notation firewall

The same-charge audit inherits the project firewall:

1. electric charge is still carried by \(\eta_Q, q_\star, q_{\rm eff}\), not by circulation;
2. grouped labels `20/21/22` refer to grouped real `P2` lanes, not spacetime indices;
3. the mixed channels
   \[
   A_w,\qquad J^w,\qquad F_{\mu w}
   \]
   remain part of the microscopic ontology outside the strict far-field brane reduction;
4. the same-charge scalar \(\Xi_1\) used later is not a new arbitrary parameter; inside the imported `5`PN language it is exactly the weak-axisymmetric load defect
   \[
   \Xi_{\rm load}=\frac{P_1}{P_0}.
   \]

---

## 4) Headline outputs / carry-forward formulas

This is the shortest memory ledger for future work.

### 4.1 Static-bundle admissibility

The static one-port reduced matrix satisfies
\[
\det \mathcal K_{\rm red}=\Delta D_0.
\]
So the natural admissibility conditions are
\[
\Omega_U^2>0,\qquad \Delta>0,\qquad D_0>0.
\]

### 4.2 Static same-charge response is quadratic and attractive-or-neutral

The exact static correction is
\[
\delta V_{\rm mix}(r)=-\frac12 J(r)^T \mathcal K_{\rm red}^{-1}J(r)\le 0.
\]

### 4.3 Product-family theorem

For the primitive reduced source pair \((x^{-3},e^{-2\kappa x}/x)\), the minimal static mixed bundle can generate only
\[
\frac{1}{x^6},
\qquad
\frac{e^{-2\kappa x}}{x^4},
\qquad
\frac{e^{-4\kappa x}}{x^2}.
\]
So no slower new attractive family appears.

### 4.4 Linear phase-lag no-go

With passive/outgoing mixed-port law \(\Pi(\omega)=i\Gamma_{\rm out}(\omega)\),
\[
\Re\,\delta \mathfrak V_{\rm mix}(x,\omega)=0
\qquad
\text{to first order in the outgoing port},
\]
while
\[
\overline P_{\rm abs}^{(1)}(x,\omega)
=
\frac{\omega}{2}\Gamma_{\rm out}(\omega)\,\mathcal T_J(\omega)^2
\ge 0.
\]

### 4.5 Weak-axisymmetric prefactor law

The exact weak-axisymmetric grouped signature is
\[
\lambda_{20}=1,\qquad \lambda_{21}=\frac12,\qquad \lambda_{22}=-1,
\]
and
\[
P_A=\bar P_0(1+\epsilon \lambda_A \Xi_1).
\]
So
\[
a_{P_0}=\frac{\epsilon\bar P_0\Xi_1}{4},
\qquad
b_{P_0}=\frac{3\epsilon\bar P_0\Xi_1}{4},
\qquad
b_{P_0}=3a_{P_0}.
\]

### 4.6 Actual-branch transported ceiling

On the weak-axisymmetric line, the robust all-lane ceiling collapses to
\[
\bar P_0(1+|\epsilon\Xi_1|)\le P_{\rm crit}.
\]
Equivalently,
\[
\frac{\Delta_{\rm norm}+T_{\rm quad}}{\hat m_0^{\,2}}(1+|\epsilon\Xi_1|)\le P_{\rm crit},
\qquad
T_{\rm quad}:=\frac{54Gc_s^5}{5a^5c^5}.
\]

### 4.7 Selected support ratio from the minimal isotropic precursor

Matching to the minimal isotropic passive/outgoing quadrupole module fixes
\[
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
\]

### 4.8 Selected twin-support curve

The selected branch collapses to the exact one-parameter curve
\[
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma=\frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23.
\]

### 4.9 Actual coherent placement point on that curve

On the coherent local D/N branch,
\[
\epsilon=\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
\[
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{2\epsilon}{1-\epsilon}.
\]

### 4.10 Final practical compiler packet

The smallest actual branch packet now needed by the audit is
\[
\Bigl(
\epsilon,\ \varrho_{\rm phys},\ \sigma_{\rm phys},\ \text{ranking region},
R_{\rm tr},\ R_{\rm target},\ \epsilon_\eta,
d\ln R_{\rm tr},\ d\ln R_{\rm target},\ d\ln \epsilon_\eta,\ N_Q-1
\Bigr).
\]

---

## 5) Static one-port audit: what was killed immediately

The earliest same-charge result carried into this file was already negative:

- the naive Maxwell-only story did not erase the barrier;
- pure inverse-sixth attraction was too weak by itself.

Stage 002 then replaced the generic mixed placeholder by the first exact one-port mixed-bundle kernel and forced a much sharper conclusion.

### 5.1 Exact static susceptibility law

With reduced coordinates \((q,U,W)\), the on-shell static correction is
\[
\delta V_{\rm mix}(r)=-\frac12 J(r)^T\mathcal K_{\rm red}^{-1}J(r).
\]
Because the admissible branch has \(\Delta>0\) and \(D_0>0\), the reduced matrix is positive definite and the correction is attractive-or-neutral.

So the first honest static mixed bundle does soften the potential, but only through a controlled quadratic susceptibility.

### 5.2 No new long-range family

For the first primitive source families, the exact static correction becomes a sum of only three product families:
\[
x^{-6},\qquad e^{-2\kappa x}/x^4,\qquad e^{-4\kappa x}/x^2.
\]

This kills the hope that the first mixed bundle might create a brand-new slow attractive tail. The static corridor is therefore **not** about discovering a new spatial law. It is only about coefficient renormalization of short-range families.

### 5.3 Why this already narrows the corridor hard

The same denominator \(\Delta D_0\) controls both:

- the static same-charge susceptibility, and
- the outgoing-normalization side of the `5`PN / `2.5`PN ledger.

So “large static same-charge softening” and “healthy outgoing normalization” are already coupled inside the one-port bundle.

The static corridor therefore survives only as:

- coefficient engineering of short-range families,
- perhaps helped by near-softening,
- but never as a free independent attractive law.

---

## 6) Dynamic audit: what linear time dependence killed

Stage 003 promoted the same one-port bundle to a frequency-dependent mixed kernel, and Stage 004 turned that into a pole/linewidth theorem.

### 6.1 Dynamic bundle does not change the kernel class

The dynamic one-port bundle still produces only the same three spatial product families. Time dependence merely makes their coefficients complex and frequency dependent.

So linear dynamics does **not** rescue the corridor by generating a different spatial interaction class.

### 6.2 Linear outgoing phase is not conservative softening

If the outgoing port is purely passive/outgoing at leading order, then the first correction is purely imaginary:
\[
\delta\mathfrak V_{\rm mix}(x,\omega)
=
-\frac{i}{2}\Gamma_{\rm out}(\omega)\,\mathcal T_J(\omega)^2
+O(\Gamma_{\rm out}^2).
\]
Its real part vanishes at first order, while the absorbed power is positive.

This kills the naive linear outgoing shortcut:
the first outgoing phase term feeds pumping / leakage / phase lag, not barrier lowering.

### 6.3 Naive resonance is also not enough

Stage 004 puts the dynamic corridor into simple-pole normal form and proves an exact dispersive/absorptive tradeoff.

What survives after this is only:

> a pole with sufficiently large residue-to-linewidth ratio might still amplify one of the already-existing short-range attractive families enough to matter.

So the dynamic corridor remains possible, but only as a resonant dispersive enhancement that beats its own absorptive load.

---

## 7) Primitive finite-throat branch, isotropic `5`PN compatibility, and transported windows

Stages 005–007 carried out the first explicit primitive test.

### 7.1 A concrete primitive branch exists

The audit built a concrete finite-throat primitive one-port branch, computed its quartic pole polynomial, residues, linewidths, and survival inequalities, and then calibrated that same branch onto the exact isotropic `5`PN / `2.5`PN target surface.

So by Stage 006 the same-charge corridor was no longer a schematic resonance story. It already had a genuine explicit primitive branch family.

### 7.2 Compatibility with the isotropic target does not kill the corridor

Placing that primitive branch on the exact isotropic target surface does **not** automatically destroy the wall-like dynamic corridor. Instead it leaves a finite normalization window.

That is one of the strongest intermediate positive results in the audit:

> the same finite-throat branch can satisfy the exact isotropic outgoing-normalization target and still retain nonzero dynamic headroom.

### 7.3 The actual-branch kill test becomes explicit

Stage 007 then transported the primitive ceilings to the actual branch packet.

For any ceiling \(P_{\rm crit}\), the exact actual-lane prefactors are compiled from
\[
\bar P_0,\quad a_{P_0},\quad b_{P_0},
\]
and the strongest transported sufficient condition is
\[
\max\{P_{20},P_{21},P_{22}\}\le P_{\rm crit}.
\]

On the weak-axisymmetric line this becomes
\[
\bar P_0(1+|\epsilon\Xi_1|)\le P_{\rm crit}.
\]

So after Stage 007 the same-charge corridor lives or dies, at transported level, by one explicit inequality in:
- the normalization defect \(\Delta_{\rm norm}\),
- the source-map factor \(\hat m_0\),
- and the outgoing-prefactor slope \(\Xi_1\).

That is the first real actual-branch kill test.

---

## 8) Mechanism sieve and strict even-gate compression

Stages 008–011 are where the corridor stops being a broad mixed-anisotropy idea and becomes a narrow specific one.

### 8.1 Wall-only and pure-BdG weak anisotropy are dead

The first-order conservative compensation equations kill:

- wall-only deformations generically;
- pure-BdG deformations on the concrete branch.

So neither “pure wall reshaping” nor “pure support reshaping” survives as the live weak-axisymmetric carrier.

### 8.2 Mixed-sector-only family survives

Activating only the mixed/U primitive slopes leaves a nontrivial compensated nullspace and nonzero \(\Xi_1\).

So the corridor survives this sieve only as a constrained mixed-sector family.

### 8.3 The same-charge scalar is the imported `5`PN load defect

Stage 009 proves exactly
\[
\Xi_{\rm load}=\Xi_1=\frac{P_1}{P_0}.
\]
So the barrier scalar is not a new invention. It is the same weak-axisymmetric prefactor-loading defect already isolated on the `5`PN side.

### 8.4 Strict even-gate package is stronger than the earlier compensation surface

The weaker Stage-008 first-order compensation surface and the stricter imported `5`PN even gates coincide only on the canonical branch.

On the explicit noncanonical compatibility point, imposing both forces the surviving corridor to satisfy more than just conservative-shape preservation.

### 8.5 Surviving strict corridor: the pure-transfer subcorridor

The full intersection of:

1. the Stage-008 conservative-shape preservation, and
2. the strict imported `5`PN even gates,

still leaves a two-dimensional mixed-sector corridor. On that subcorridor,
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\Xi_{\rm load}=\frac{N_{01}}{N_0}.
\]

This is the cleanest surviving mechanism in the whole file.

It means the corridor can survive all gates up to this point with the conservative one-pole bundle frozen at first order. The remaining effect is carried purely by outgoing-transfer loading.

### 8.6 Outgoing-rigidity sieve

Stage 010 then tests whether this pure-transfer corridor survives additional rigidity assumptions.

The combined interference-plus-hybridization rigidity kills the corridor, but one-rigidity survivors remain.

So the same-charge corridor remains alive only if the mixed outgoing sector is **not** over-rigid.

### 8.7 Numerator/denominator split

Stage 011 splits the pure-transfer scalar exactly as
\[
\Xi_1=2(\pi_1-\delta_1),
\qquad
\pi_1=\frac{P_{01}}{P},
\qquad
\delta_1=\frac{\Delta_{01}}{\Delta}.
\]

This produces two one-dimensional rigid survivors on the concrete branch:

- numerator-rigid \(\pi_1=0\), carrying \(\Xi_1=-2\delta_1\);
- denominator-rigid \(\delta_1=0\), carrying \(\Xi_1=2\pi_1\).

The numerator-rigid branch gives a larger static \(\Xi_1\) but split-sign dynamic response. The denominator-rigid branch gives a smaller static \(\Xi_1\) and hurts both dynamic wall-like poles.

Yet in both cases the dynamic ceilings remain looser than the transported static ceilings. So even here the corridor is still **static-first**, not dynamic-first.

---

## 9) Selected-branch classifier and the static-first theorem

Stages 012–014 convert the pure-transfer split into a universal selected-branch classifier.

### 9.1 The selected branch is an exact numerator/denominator co-loading product

The actual selected PDE branch is not literally either rigid subcorridor from Stage 011. It is an exact numerator/denominator co-loading branch with classifier
\[
\mathcal R_{ND}(\xi,\delta)
=
\frac{72\delta^2(1-\xi)}
{(9\delta+11\xi)(9\delta^2+18\delta\xi+11\xi^2)}.
\]

So the physical selected branch must be read through this classifier rather than by forcing it into one of the Stage-011 rigid axes.

### 9.2 Near softening, the selected branch is always denominator-like

The selected-branch classifier proves:

- if \(\delta\ge 8/9\), the whole stable branch is denominator-like;
- if \(0<\delta<8/9\), the branch begins numerator-like but crosses uniquely to denominator-like;
- near softening — and therefore near large normalization gain — the selected branch is always denominator-like.

So the physically relevant high-gain end of the selected branch is automatically in the denominator-like sector.

### 9.3 Dynamic windows never become the first bottleneck

Stage 013 then compiles the selected-branch classifier into actual wall-like dynamic windows and proves the key ordering theorem:

\[
B_{\rm dyn}^{(\rm both)} > B_{\rm stat}^{(\rm both)},
\qquad
B_{\rm dyn}^{(\rm nonempty)} > B_{\rm stat}^{(\rm nonempty)}.
\]

So for every selected-branch signature, the dynamic wall-window remains looser than the transported static ceiling.

This is one of the load-bearing results of the full audit:

> the first kill condition on the selected same-charge branch is static \(\Xi_1\)-budget, not the dynamic wall window.

### 9.4 Pullback to the continuum placement variables

Stage 014 then pulls the classifier all the way back to the continuum placement ratios \((\delta, M_{\rm mix}, R_{\rm target})\).

This preserves the same theorem structure:

- larger \(R_{\rm target}\) pushes the branch denominator-like;
- at fixed product scale, larger \(M_{\rm mix}\) pushes it numerator-like;
- the static-first theorem survives the physical pullback.

So the selected-branch dynamic classification is now already expressed in physical continuum placement language.

---

## 10) Known `5`PN data injection and why support/source is no longer the bottleneck

Stage 015 injects the numerically located support/source data already known from the `5`PN stack.

This is the point where the audit stops asking whether support enhancement might simply be too weak.

### 10.1 The support/source side is numerically safe

Once the known `5`PN support/source data are inserted, the support/source side is comfortably non-bottlenecked. The canonical passive/outgoing normalization side also remains exact on the natural branch.

So the same-charge corridor is no longer waiting on more support algebra.

### 10.2 The unresolved kill condition moves to the actual coherent placement packet

After this injection the first unresolved theorem gate becomes:
\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln \epsilon_\eta=0,
\]
together with
\[
N_Q=1.
\]

So the real same-charge verdict is no longer about whether the reduced support/source side can be made large enough. It is about whether the actual moving-throat branch lands at the right coherent placement / orbit-lock point.

---

## 11) Rigid-mouth packet factorization and the post-static dressing theorem

Stages 016–022 take that unresolved packet and factor it completely.

### 11.1 Direct static scalar on the branch-observable side

Stage 017 shows that the unresolved static scalar narrows to
\[
\Xi_1
=
-\delta\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,\delta\ln\epsilon_\eta.
\]

So after support/source and dynamic-window issues are removed, the static branch bottleneck is already only a two-observable combination.

### 11.2 Rigid-mouth finite packet and diagonal chart

On the actual rigid-mouth branch, the physical logarithmic variables
\[
U:=\ln\frac{\mathcal T^2}{\mathcal T_{\rm ref}^2},
\qquad
V:=\ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}},
\]
are already the exact finite packet coordinates.

The third direct observable obeys
\[
\frac{R_{\rm target}}{R_{\rm target,ref}}
=
\frac{1-\epsilon_{\eta,\rm ref}e^V}{1-\epsilon_{\eta,\rm ref}}\,e^{-U}.
\]

The canonical packet projectors are diagonal:
\[
P_{\mathcal T}=
\begin{pmatrix}1&0\\0&0\end{pmatrix},
\qquad
P_\eta=
\begin{pmatrix}0&0\\0&1\end{pmatrix}.
\]

So the packet splits exactly into:

- a pure transfer-shape leg \(U\),
- a pure dressing leg \(V\).

### 11.3 Microscopic dependent correction compiler

In this normal form the dependent microscopic correction is
\[
(\Delta_T,\Delta_{K_\eta},\Delta_\mu)=(0,-V,U-V).
\]

So:

- the static-only correction changes only \(\mu_W\);
- the post-static dressing correction is the equal \(K_\eta\)-\(\mu_W\) shift.

### 11.4 Support-blindness

Stages 020–022 prove that the whole rigid-mouth physical normal form and its microscopic correction compiler are support-blind.

This is another major compression result:
the final same-charge packet no longer depends on the explicit support lane once the rigid-mouth packet is adopted.

### 11.5 Post-static theorem

After the transfer-shape gate is cleared, the surviving amplitude is simply
\[
q_\eta = \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

Equivalently, once tracking rigidity and transfer-shape rigidity are imposed, the only remaining same-charge post-static test is:

> compute \(\epsilon_\eta\) on the actual rigid-mouth branch and check whether it is invariant.

That is the sharpest post-static barrier statement reached in the whole file.

---

## 12) Selected support slice: minimal isotropic precursor, twin-support curve, and actual coherent placement

Stages 023–025 finish the support-side classification by fixing the selected slice exactly.

### 12.1 Minimal isotropic passive/outgoing precursor fixes the support ratio

Matching to the minimal isotropic passive/outgoing quadrupole precursor gives:
\[
c_0=\frac34,
\qquad
c_1=\frac14,
\]
hence
\[
\rho_\alpha=\frac{1}{c_0}=\frac43,
\qquad
\zeta_{\rm req}=\frac{c_1}{c_0}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
\]

So the selected same-charge branch no longer scans all support sectors. It is fixed to one exact support ratio.

### 12.2 Exact support regime meaning

Because
\[
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix},
\]
the selected branch lies strictly inside the **symmetric lowest twin** window:

- mixed-only is not enough;
- the lowest symmetric twin is enough;
- non-twin asymmetry is not required.

### 12.3 Exact selected twin-support curve

The selected support slice is the exact one-parameter curve
\[
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma = \frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23.
\]

Two exact thresholds govern primitive ranking along that curve:
\[
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\]
\[
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)},
\]
with
\[
0<\varrho_{W\Lambda}<\varrho_{U\Lambda}<\frac23.
\]

So the primitive hierarchy on the selected branch splits into exactly three regions:

- low \(\varrho\): \(q_\chi>q_Z>q_\Lambda>q_W>|q_U|\);
- intermediate \(\varrho\): \(q_\chi>q_Z>q_W>q_\Lambda>|q_U|\);
- high \(\varrho\): \(q_\chi>q_Z>q_W>|q_U|>q_\Lambda\).

### 12.4 Actual coherent placement on that curve

On the coherent local D/N branch the actual selected point is set by
\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
and then
\[
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{2\epsilon}{1-\epsilon}.
\]

So the actual coherent placement state is determined as soon as the completed PDE returns \((\delta_U,\epsilon_W)\).

### 12.5 Orbit-side observables are support-blind

On the same coherent branch,
\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]
\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

The support ratio \(\zeta\) does not enter \(\epsilon\), \(R_{\rm tr}\), or \(R_{\rm target}\).
So the orbit packet is exactly support-blind.

### 12.6 Final coherent orbit-lock theorem gate

The actual coherent local D/N branch satisfies orbit lock iff
\[
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
\]
equivalently iff
\[
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln \epsilon_\eta=0.
\]

The outgoing finish line remains separate:
\[
N_Q=1,
\qquad
\text{equivalently }\chi_Q=1
\text{ on the natural source-map branch.}
\]

So the final same-charge / moving-throat endgame is now exactly:

1. extract the stationary coherent placement packet;
2. place it on the selected twin-support curve via \(\epsilon\);
3. compute
   \[
   d\ln R_{\rm tr},\quad d\ln R_{\rm target},\quad d\ln \epsilon_\eta;
   \]
4. and evaluate the outgoing finish-line defect \(N_Q-1\).

That is the real final theorem gate now.

---

## 13) Final verdict and next honest theorem gate

The same-charge corridor is still alive, but only in a sharply compressed form.

### 13.1 What is dead

The audit has already killed all of the following as primary explanations:

- naive Maxwell-only same-charge shaping;
- the hope for a new static mixed long-range attractive family;
- linear outgoing-phase barrier lowering;
- naive resonance as an unconditional rescue;
- wall-only weak-axisymmetric compensation;
- pure-BdG weak-axisymmetric compensation;
- support/source starvation as the current bottleneck.

### 13.2 What survives

What remains alive is:

1. a mixed-sector weak-axisymmetric corridor;
2. further narrowed by the strict even-gate package to a pure-transfer subcorridor;
3. further narrowed by the selected-branch classifier to a static-first denominator-vs-numerator co-loading problem;
4. further narrowed by the rigid-mouth normal form to transfer-shape plus dressing axes;
5. further narrowed by the selected support analysis to the symmetric lowest twin-support slice; and
6. finally reduced to the actual coherent orbit-lock / normalization packet.

### 13.3 Strongest honest carry-forward reading

The best present reading is:

> **Inside the current reduced stack, same-charge survival is no longer a problem of inventing new attractive laws or hoping for generic resonance. It is a problem of whether the actual coherent moving-throat branch lands on the narrow pure-transfer / selected-twin-support corridor while keeping the rigid-mouth packet locked and the canonical outgoing normalization intact.**

### 13.4 Smallest next theorem target

The next honest task is not another broad barrier model. It is to compute the exact actual-branch packet
\[
\Bigl(
\epsilon,\ \varrho_{\rm phys},\ \sigma_{\rm phys},\ \text{ranking region},
R_{\rm tr},\ R_{\rm target},\ \epsilon_\eta,
d\ln R_{\rm tr},\ d\ln R_{\rm target},\ d\ln \epsilon_\eta,\ N_Q-1
\Bigr),
\]
and then answer four concrete questions:

1. does the coherent branch sit on the selected symmetric twin-support slice as expected?
2. does the rigid-mouth packet first clear the transfer-shape gate?
3. after that, does the dressing leg collapse, i.e. is \(\epsilon_\eta\) invariant?
4. does the outgoing finish line satisfy \(N_Q=1\)?

That is now the clean finish line for the same-charge barrier audit.
