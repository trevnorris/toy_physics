g2_summary.md — Comprehensive summary (quartic electron g-2 / outgoing-normalization package)

## 0) What this package is doing

This document is the **carry-forward g-2 summary** for the moving-throat / 4D toy-model program.
It is built from the long derivation record in `g2_full_output.md`, but it is written in the same role as the other `4d*summary.md` files: not as a paper draft, and not as a raw session log, but as a **theorem ledger** that freezes what is now structurally fixed, what remains conditional, and what the next honest theorem gate actually is.

Its job is to compress the full reconstruction chain into one readable package:

1. freeze the **exact current local/staggered anomaly baseline** and the true electron-point residual;
2. replace the old sequential “charge then inertia then charge” bookkeeping by the exact moving-throat **quotient coordinates**;
3. identify the universal coherent defect as the **transfer-shape / outgoing-prefactor slope**;
4. show how the support/source, grouped real `P2`, and isotropic outgoing ledgers collapse the whole reduced problem to one scalar **normalization defect**;
5. audit explicit isotropic outlet/core realizations;
6. show what the **adiabatic no-tuning branch** actually predicts; and
7. isolate the one remaining microscopic issue if one insists on the exact electron sliver.

The strongest honest carry-forward statement is:

> **Within the present reduced moving-throat hierarchy, the canonical compact outgoing grouped-`P2` branch is naturally fixed with no tuning, and it reproduces the carried local anomaly value. The exact electron-point quartic sliver is *not* yet forced at first order. What remains is one sharply localized microscopic slippage / odd-outlet datum, not a broad fit space.**

So the present package should be read as the g-2 analog of the 2.5PN/4PN summaries: the algebraic compression is essentially finished, but the final physical branch-selection datum is still open.

---

## 1) Claim taxonomy and how to read this package

The g-2 package mixes exact algebra, reduced branch theorems, and one remaining conditional physical datum. The cleanest reading uses four classes.

### 1.1 Exact

These are exact symbolic or algebraic statements once the carried benchmark inputs, quotient definitions, and reduced branch dictionaries are fixed.

Examples:

- the exact current local/staggered benchmark law;
- the rank-`3` quotient map and nullity-`5` similarity orbit;
- the coherent-defect normal form
  \[
  \Theta_1=-C_{\rm tr}q_{\rm tr},
  \qquad
  \Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt},
  \qquad
  \mathcal R_1+\Xi_1=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta;
  \]
- the identity
  \[
  \Xi_1 = \delta\ln \mathcal T^2 = \frac{P_1}{P_0};
  \]
- the forced isotropic conservative quadrupole module
  \[
  \widehat Y_Q^{\rm cons}(\omega)=\frac34+\frac{1}{4(1-\omega^2/\Omega_Q^2)};
  \]
- the exact compact outgoing `l=2` DtN fingerprint and the canonical value `\chi_Q=1`;
- the exact finite-`f` electron target surface in isotropic DtN variables.

### 1.2 Exact within closure

These are exact once a stated reduced closure family is adopted.

Examples:

- the coherent local D/N support/source branch;
- the compensated Robin–mixed outlet family;
- the balanced two-channel throat core;
- the adiabatic-wall ground-state closure;
- the natural point-particle source-map branch;
- the minimal isotropic grouped-`P2` one-pole module.

So these statements are not “heuristics,” but they do depend on a declared branch or reduced family rather than on a completed full PDE theorem.

### 1.3 Reduced / controlled reduction

These are statements that use the same kind of reduced hierarchy already familiar from the 2.5PN–5PN moving-throat program.

Examples:

- the weak-axisymmetric grouped reduction;
- the one-port effective transfer-shape collapse;
- the selected twin-support and support-compensation reductions;
- the reduced outlet/core realization audit;
- the grouped-`P2` normalization bridge back to the point-particle g-2 anomaly law.

This package is not an unreduced theorem of the full moving-throat PDE.

### 1.4 Conditional electron-match statement

The final “exact electron value” statement is still conditional.

The package shows that:

- the canonical outgoing branch is fixed naturally;
- the last reduced mismatch is only one scalar defect;
- the no-tuning adiabatic compensated branch keeps that defect zero at first order.

So the only remaining gap is whether the **actual moving-throat branch** generates a nonzero microscopic slippage, such as a genuine off-family scalar slippage or odd mixed-port renormalization.

---

## 2) What this package claims — and what it does not claim

### 2.1 What it claims

This package claims all of the following.

1. The present local/staggered g-2 baseline can be reconstructed exactly, including its quartic term and its actual residual relative to the electron target.
2. The new moving-throat PDE language replaces old sequential bookkeeping by an exact quotient problem with:
   \[
   \dim(\text{quotient})=3,
   \qquad
   \dim(\text{similarity orbit})=5.
   \]
3. The first genuinely common quartic anomaly layer acts only on the already-existing **transport residue**, not on the entire anomaly law.
4. The coherent weak-axisymmetric defect compresses to the exact scalar
   \[
   \Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt},
   \]
   and this same scalar is exactly
   \[
   \Xi_1 = \delta\ln \mathcal T^2 = \frac{P_1}{P_0}.
   \]
5. The support/source side, the conservative grouped real `P2` side, and the outgoing normalization side all compress the isotropic point-particle finish line to one scalar defect.
6. The actual isotropic grouped-`P2` conservative branch forces the `3/4 + 1/4` module and therefore the loading ratio
   \[
   \rho_\alpha = \frac43,
   \qquad
   \zeta_{\rm req}=\frac13.
   \]
7. The explicit compact outgoing `l=2` DtN branch fixes the canonical outgoing normalization exactly:
   \[
   \chi_Q=1,
   \qquad
   N_Q=1
   \]
   on the natural point-particle source-map branch.
8. The exact electron-point quartic sliver can be rewritten as a tiny finite-`f` outgoing-normalization drift.
9. The first explicit isotropic outlet audit rules out the naive standalone channels and leaves only the compensated Robin–mixed class alive as a serious exact candidate.
10. The balanced explicit throat core collapses to a one-parameter parent family, and positive mouth sourcing selects the lower Family-1 branch without fine tuning.
11. On the adiabatic no-tuning compensated branch, the parent compensation sheet and even geometry remain frozen and the first-order outgoing defect vanishes.
12. Therefore the current no-tuning branch predicts the carried local value, not the exact electron sliver.

### 2.2 What it does **not** claim

This package does **not** claim:

- a full first-principles moving-throat PDE theorem for the exact electron-point sliver;
- a proof that the physical branch must generate a nonzero `\delta\gamma_W`, `\delta\mathfrak B_W`, or `\varepsilon_\perp`;
- a proof that the exact electron value follows automatically from the current frozen stack;
- a theorem that every explicit reduced outlet/core family considered here is the actual microscopic one.

So the correct reading is:

> the canonical outgoing background is now fixed, the anomaly bookkeeping has collapsed almost completely, and the remaining open problem is one narrow microscopic branch-selection datum.

---

## 3) Fixed inputs, notation firewall, and benchmark numbers

### 3.1 Benchmark anomaly inputs

The g-2 reconstruction carries the following benchmark data.

Define
\[
 f := \frac{\alpha_{\rm fs}}{2\pi}.
\]
The frozen local/staggered benchmark uses
\[
 \kappa = 1.177746578880.
\]
The reconstructed local value is
\[
 \boxed{g_{\rm loc} \approx 2.002319304358647956.}
\]
The target adopted in the derivation chain is
\[
 \boxed{g_e \approx 2.00231930436092.}
\]
So the residual is
\[
 \boxed{
 \Delta g = g_e-g_{\rm loc}
 \approx 2.27204390584705\times 10^{-12}.
 }
\]

### 3.2 Exact current local/staggered law

The exact rebuilt local/staggered series is
\[
\frac{g_{\rm loc}(f)}{2}
=
1+f-\frac{47}{36}f^2+c_{3,\rm total}f^3+a_{4,\rm staggered}f^4+O(f^5),
\]
with
\[
 c_{3,\rm total}
 =
 \frac{11}{6}\kappa + \frac{4-\pi}{\pi^2\kappa}
 \approx 2.23305058688410,
\]
\[
 a_{4,\rm staggered}
 =
 -\frac{55}{6}\kappa^2+\frac{4(\pi-3)}{\pi^3\kappa}
 \approx -12.6994546522869.
\]
The important interpretive point is that the current exact local law already has its own quartic coefficient. The unresolved common layer is therefore **incremental** on top of that law, not a replacement for it.

### 3.3 Quartic residual target

Define the residual quartic coefficient against the exact local baseline by
\[
 a_{4,\rm resid}:=\frac{\Delta g}{2f^4}
 \approx 0.624374101073809.
\]
Then the first common scalar tangent is fixed by
\[
 \boxed{
 \Lambda_1 = \frac{a_{4,\rm resid}}{c_{3,\rm total}}
 \approx 0.279605891931464.
 }
\]
This `\Lambda_1` is the universal quartic target carried all the way through the package.

### 3.4 Notation firewall

The g-2 package inherits the broader project firewall.

1. Electric charge is still carried by `\eta_Q, q_\star, q_{\rm eff}` and is not redefined by circulation.
2. Grouped labels `20/21/22` refer to grouped real `P2` lanes, not spacetime indices.
3. The mixed channels
   \[
   A_w,\qquad F_{\mu w},\qquad J^w
   \]
   are suppressed only in the strict far-field brane reduction. They remain part of the microscopic ontology and are crucial in the outlet audit.
4. In the g-2 chain the central defect variables are:
   \[
   q_{\rm tr},\qquad q_{\rm nt},\qquad q_\eta,
   \]
   the coherent defect
   \[
   \Xi_1,
   \]
   the outgoing branch scalar
   \[
   \chi_Q,
   \]
   and the grouped-normalization defect
   \[
   N_Q.
   \]

---

## 4) Headline outputs / memory ledger

This is the shortest high-value ledger to keep in working memory.

### 4.1 Exact quotient reduction

The microscopic grouped weak-axisymmetric drift space has dimension `8`, but the exact monomial-drift matrix has rank `3`, so the physical reduced motion lives in the quotient coordinates
\[
 q_{\rm tr},\qquad q_{\rm nt},\qquad q_\eta,
\]
while the remaining `5` directions are exact similarity-orbit motion.

### 4.2 First common scalar

The first common anomaly layer acts only on the transport residue and is controlled by the scalar tangent
\[
 \boxed{\Lambda_1\approx 0.279605891931464.}
\]

### 4.3 Exact coherent normal form

On the coherent branch the observable weak-axisymmetric defect channels are
\[
 \Theta_1=-C_{\rm tr}q_{\rm tr},
\]
\[
 \boxed{\Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt},}
\]
\[
 \mathcal R_1+\Xi_1=-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}q_\eta.
\]
So the direct quartic target is either
\[
 \Xi_1=\Lambda_1,
\]
or, on the tracking-rigid specialization,
\[
 q_{\rm nt}=\Lambda_1.
\]

### 4.4 Transfer shape and outgoing prefactor

The same coherent defect is exactly
\[
 \boxed{\Xi_1=\delta\ln \mathcal T^2=\frac{P_1}{P_0}.}
\]
So the quartic g-2 sliver is already the first normalized slope of the outgoing grouped-`P2` prefactor.

### 4.5 Microscopic demand-ratio / mixed-outgoing law

The quartic transfer-shape drift is carried entirely by mixed/outgoing microdata:
\[
 \boxed{
 \Lambda_1 = -q_\Lambda + q_Z + 2q_\chi + 2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
 }
\]
In particular, on the coherent D/N branch the support ratio `\zeta` drops out of both `R_{\rm tr}` and `R_{\rm target}`.

### 4.6 Support compensation

The coherent local D/N support lane is structurally a **fixed-target / load-compensation** mechanism. Its natural continuation lies on the coherent side
\[
 \delta\ln R_{\rm target}=0,
\]
not on a direct retargeting law.

### 4.7 Forced conservative grouped-`P2` module

If the higher conservative branch is one isotropic grouped-`P2` pole plus a static geometry completion, then the exact minimal isotropic identity forces
\[
 \boxed{
 \widehat Y_Q^{\rm cons}(\omega)
 =
 \frac34+
 \frac{1}{4}\frac{1}{1-\omega^2/\Omega_Q^2}.
 }
\]
Equivalently,
\[
 K_{\rm geom}=3K_{\rm pole},
 \qquad
 K_{\rm pole}=\frac{K_0}{4}.
\]
And the support/source side becomes automatic:
\[
 \boxed{\rho_\alpha=\frac43,\qquad \zeta_{\rm req}=\frac13.}
\]

### 4.8 Actual isotropic branch = one scalar defect

On the actual isotropic branch the geometry lane is exactly orthogonal to the grouped real `l=2` bundle at linear order, so the dynamic geometry contamination vanishes:
\[
 \epsilon_2=\epsilon_4=0.
\]
The whole reduced mismatch is therefore one scalar normalization defect,
\[
 \boxed{N_Q := \bar K_0 / \bar K_0^{\rm target}.}
\]

### 4.9 Canonical outgoing `l=2` DtN branch

The explicit compact outgoing `l=2` DtN fingerprint is
\[
 \widehat Y_2^{\rm out}(z)
 =
 1+\frac{z^2}{9}+\frac{4z^4}{81}+i\frac{z^5}{27}+O(z^6),
 \qquad
 z:=\frac{a\omega}{c_s}.
\]
Matching it to the minimal retarded grouped-`P2` module fixes
\[
 \boxed{\chi_Q=1.}
\]
On the natural source-map branch,
\[
 \hat m_0=1,
 \qquad
 \hat m_0^{\,2}\chi_Q N_Q=1,
\]
so
\[
 \boxed{N_Q=1.}
\]
Thus the canonical outgoing background is already a genuine no-tuning closure.

### 4.10 Electron-point outgoing defect

The exact electron-point quartic sliver can be encoded as a tiny finite-`f` outgoing-normalization defect:
\[
 \boxed{
 \Delta_Q^{(e)}
 =
 -\frac{\Lambda_1 f}{1+\Lambda_1 f}
 \approx -3.24631584151692\times 10^{-4}.
 }
\]
Equivalently,
\[
 \boxed{
 \chi_Q^{(e)}\approx 0.999675368415848,
 \qquad
 N_Q^{(e)}\approx 1.00032473700404,
 \qquad
 \ell := \ln(1+\Lambda_1 f)\approx 3.2468428839\times 10^{-4}.
 }
\]
So the entire electron-point deformation is only a few parts in `10^4` away from the canonical outgoing branch.

### 4.11 Exact isotropic DtN anomaly surface

The exact finite-`f` isotropic DtN branch-selection surface is
\[
 \boxed{
 3\bigl(S\beta^5+9\Sigma_5\bigr)(1+\Lambda_1 f)=3S-\Sigma_0.
 }
\]
Its tangent reduction is
\[
 \boxed{5b+\frac{a_0}{3}+9a_5=-\Lambda_1.}
\]
On the adiabatic even-frozen slice,
\[
 \beta=1,
 \qquad
 \Sigma_0=0,
 \qquad
 \Sigma_5 = \frac{S}{9}(e^{-\ell}-1),
\]
so the whole anomaly is pure odd in normalized DtN language.

### 4.12 No-tuning adiabatic prediction

On the present no-tuning adiabatic compensated branch,
\[
 \boxed{
 \Delta_Q^{(\rm nat)}=0,
 \qquad
 \chi_Q^{(\rm nat)}=1,
 \qquad
 N_Q^{(\rm nat)}=1.
 }
\]
So the branch prediction is simply
\[
 \boxed{g_{\rm pred}^{(\rm nat)} = g_{\rm loc} \approx 2.00231930435865,}
\]
with miss
\[
 \boxed{|g_e|-g_{\rm pred}^{(\rm nat)}\approx 2.27\times 10^{-12}.}
\]

---

## 5) From staggered bookkeeping to exact quotient closure

The first large structural gain of the package is that the old staged transport story is replaced by a genuine quotient geometry.

### 5.1 Exact microscopic monomials

The coherent branch is controlled by three direct microscopic monomials:
\[
 \boxed{
 \mathfrak C_{{\rm tr},*},
 \qquad
 \mathfrak C_{{\rm nt},*},
 \qquad
 \epsilon_\eta.
 }
\]
A convenient explicit choice is
\[
 \mathfrak C_{{\rm tr},*}
 =
 \left(\frac{\gamma c_{\eta U}}{K_U}\right)^{1+\delta_{U,*}}
 \left(\frac{\pi^2T_U}{L^2K_U}\right)^{1+\chi_{0,*}},
\]
\[
 \mathfrak C_{{\rm nt},*}
 =
 \frac{Z_W}{\Omega_W^2}
 \epsilon_W^{E_*}
 \delta_U^{-F_*},
\]
\[
 \epsilon_\eta = \frac{c_{\eta U}^2}{K_U K_\eta^{(\rm eff)}}.
\]
Their logarithmic drifts are exactly
\[
 \delta\ln\mathfrak C_{{\rm tr},*}=q_{\rm tr},
 \qquad
 \delta\ln\mathfrak C_{{\rm nt},*}=q_{\rm nt},
 \qquad
 \delta\ln\epsilon_\eta=q_\eta.
\]

### 5.2 Similarity orbit

The direct monomial-drift matrix has rank `3`. Therefore the exact zero-defect tangent space has dimension `5`, and the monomial-preserving weak-axisymmetric branch is the tangent space of an exact five-parameter multiplicative similarity orbit.

So the g-2 chain is not fundamentally an eight-variable microscopic drift problem. It is exactly a three-coordinate quotient problem.

### 5.3 Quartic anomaly gate in quotient variables

Because
\[
 \Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt},
\]
the quartic anomaly gate becomes
\[
 \boxed{
 A_{\rm tr}\,\delta\ln\mathfrak C_{{\rm tr},*}
 +
 \delta\ln\mathfrak C_{{\rm nt},*}
 =
 \Lambda_1.
 }
\]
On the tracking-rigid specialization this collapses to
\[
 \boxed{\delta\ln\mathfrak C_{{\rm nt},*}=\Lambda_1.}
\]
So the remaining anomaly datum is already isolated as a single quotient-level transfer-shape drift.

---

## 6) Support-demand selection and mixed/outgoing balance

### 6.1 Transfer-shape drift is support-blind at the key point

On the coherent D/N branch, the package proves that the support ratio `\zeta` does **not** enter the tracking factor `R_{\rm tr}` or the demand ratio `R_{\rm target}`. The universal transfer-shape drift is instead carried by the mixed/outgoing variables.

The sharp microscopic law is
\[
 \Lambda_1
 =
 -q_\Lambda+q_Z+2q_\chi+2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
\]
This is the first point where the anomaly is localized entirely inside the mixed/outgoing microdata rather than diffuse charge/inertia bookkeeping.

### 6.2 Support is fixed-target / load-compensation

Once the coherent local D/N support theorem is inserted, the physical support lane appears not as direct retargeting, but as a **load-compensation mechanism** at fixed demand ratio. In the package’s own language, the natural continuation lies on the coherent side
\[
 \delta\ln R_{\rm target}=0,
\]
with support increasing the available baseline while leaving the selected demand ratio unchanged.

### 6.3 Primitive hierarchy: wall-blocking dominates split-`U`

Resolving the reduced split-blocking variable back into primitive drifts gives
\[
 q_\epsilon = q_W - \beta q_U,
 \qquad
 \beta=\frac{2}{11+9\delta_{U,*}}<\frac{2}{11}.
\]
On the minimum-norm primitive coherent repair,
\[
 \boxed{q_U^{\min}=-\beta q_W^{\min}.}
\]
So the quartic repair is structurally **wall-blocking dominated**. The split-`U` drift is a smaller opposite-sign companion rather than the main carrier.

### 6.4 Support/source selector stops being the real bottleneck

The package spends a substantial number of steps scanning support-demand regimes, twin-support windows, and selector surfaces. But after the grouped real `P2` + static geometry split is enforced, those support ratios stop being the deepest issue. The decisive branch structure is already captured by the conservative/outgoing grouped-`P2` module and its normalization defect.

---

## 7) Conservative grouped-`P2` finish line and the canonical outgoing branch

### 7.1 Forced conservative module

The conservative finish line is no longer a free eight-slot residual problem. Once the branch is organized as one isotropic grouped-`P2` pole plus a static geometry completion, the exact isotropic identity
\[
 K_0 K_4 = 4 K_2^2
\]
forces
\[
 K_{\rm geom}=3K_{\rm pole},
\]
and therefore
\[
 \widehat Y_Q^{\rm cons}(\omega)
 =
 \frac34+
 \frac{1}{4(1-\omega^2/\Omega_Q^2)}.
\]
Thus the familiar `3/4 + 1/4` split is not an imported guess. It is a forced reduced consequence of the grouped real `P2` + static-geometry organization.

### 7.2 Actual isotropic branch = one normalization defect

Because the isotropic quadratic wall operator keeps the `l=0` geometry lane orthogonal to the grouped real `l=2` bundle, there is no dynamic geometry contamination at `O(\omega^2)` or `O(\omega^4)`. This means the actual isotropic passive/outgoing branch is controlled by one scalar normalization defect,
\[
 N_Q=\bar K_0/\bar K_0^{\rm target}.
\]
All other reduced support/source bookkeeping becomes automatic once that scalar is fixed.

### 7.3 Canonical compact outgoing `l=2` DtN branch

The explicit compact outgoing `l=2` DtN operator fixes the canonical branch exactly. Matching its low-frequency fingerprint to the retarded grouped-`P2` one-pole module gives
\[
 \chi_Q=1,
\]
with
\[
 \Omega_Q=\frac{3c_s}{2a},
 \qquad
 \sigma_Q^{\rm can}=\frac{4a^5}{27c_s^5}.
\]
On the natural point-particle source-map branch,
\[
 \hat m_0=1,
 \qquad
 \hat m_0^{\,2}\chi_Q N_Q=1,
\]
so the same canonical branch also fixes
\[
 N_Q=1.
\]
Thus the canonical outgoing background is now a genuine no-tuning closure, not a plausible ansatz.

---

## 8) Explicit outlet audit and parent-family realization

### 8.1 First outlet sieve

The first explicit isotropic outlet audit tests three classes.

1. **Pure Robin core.** It can reproduce the *size* of the electron outgoing defect, but it distorts the even `l=2` branch and therefore cannot be the whole story.
2. **Standalone mixed side-channel pole.** It cannot preserve the canonical even branch unless it vanishes, so the naive isolated mixed-pole model is ruled out.
3. **Hybrid Robin–mixed outlet.** This is the first serious surviving class.

### 8.2 Surviving compensated Robin–mixed outlet

On the nontrivial compensated branch,
\[
 \boxed{\rho_R=4\sigma_W,\qquad \kappa_W=\frac13.}
\]
The exact normalization law is
\[
 \chi_Q^{\rm hyb}
 =
 \frac{1-9\sigma_W\gamma_W}{1-\sigma_W}.
\]
Matching it to the electron target gives the exact anomaly family
\[
 \boxed{
 \gamma_W
 =
 \frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
 }
\]
So on the surviving outlet class the anomaly is carried entirely by a tiny odd mixed-channel detuning above its canonical value `1/9`.

### 8.3 Concrete balanced throat core

The surviving outlet can be pushed into a concrete two-channel throat core with
\[
 \rho_c=\frac{g_s^2}{K_s},
 \qquad
 \sigma_c=
 \frac{(K_sg_q-\lambda g_s)^2}{K_s^2K_q(1+r_c)},
 \qquad
 r_c=\frac{\lambda^2}{K_sK_q},
\]
\[
 \kappa_c=\frac{\kappa_0}{1+r_c},
 \qquad
 \gamma_c=\frac{\gamma_0}{1+r_c}.
\]
The canonical-even balance surface is the exact codimension-one law
\[
 \boxed{g_s^2(K_sK_q+\lambda^2)=4(K_sg_q-\lambda g_s)^2.}
\]
In parent-ratio variables
\[
 \mathfrak r := \frac{\lambda}{\sqrt{K_sK_q}},
 \qquad
 \mathfrak g := \frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
this becomes
\[
 \boxed{1+\mathfrak r^2 = 4(\mathfrak g-\mathfrak r)^2.}
\]
The same ratio also fixes the auxiliary D/N tube geometry:
\[
 \boxed{L_W=\frac{\pi a}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.}
\]
On the balance surface,
\[
 \boxed{\sigma_c = \rho_c/4.}
\]
So the surviving outlet is no longer a diffuse reduced fit. It collapses to a one-parameter parent family plus one small odd detuning.

### 8.4 Exact anomaly law on the parent family

With
\[
 x:=\Lambda_1 f,
\]
the exact electron target on the balanced hybrid outlet gives
\[
 \boxed{
 \gamma_c
 =
 \frac{\rho_c+4x}{9\rho_c(1+x)},
 \qquad
 \gamma_0
 =
 \frac{(1+\mathfrak r^2)(\rho_c+4x)}{9\rho_c(1+x)}.
 }
\]
So the anomaly is carried only by a tiny odd detuning away from the canonical values
\[
 \gamma_c=\frac19,
 \qquad
 \gamma_0=\frac{1+\mathfrak r^2}{9}.
\]

### 8.5 Positive mouth sourcing selects the lower Family-1 branch

For the carried explicit Family-1 geometry,
\[
 \mathfrak r_{F1}\approx 1.77799353547498,
\]
so the two balanced branches are
\[
 \mathfrak g_-^{F1}\approx 0.758035078944663,
 \qquad
 \mathfrak g_+^{F1}\approx 2.79795199200529.
\]
But any positive normalized mouth source satisfies
\[
 0\le \mathfrak g[\sigma]\le 1,
\]
so the upper branch is physically impossible while the lower branch is the unique admissible canonical candidate.

This lower branch also sits close to natural source laws. The self-matched derivative profile gives `\pi/4`, only about `3.61%` away in traction from the exact lower Family-1 value. So the balanced lower branch is not a fine-tuned construction.

---

## 9) Adiabatic no-tuning closure and strongest honest status

### 9.1 Lower-branch transport collapses strongly

Once the exact lower Family-1 branch is imposed, most apparent microscopic freedom co-transports. The side-tube length satisfies
\[
 \delta\ln L_W = \delta\ln a,
\]
and, after the frozen `n=5` wall-EOS relation is inserted, the actual lower-branch drift problem collapses to four irreducible drifts:
\[
 \delta\ln\mathcal Z_q,
 \qquad
 \delta\ln\rho_w,
 \qquad
 \delta\ln c_s,
 \qquad
 \delta\ln a.
\]
The old first-order off-family scalar vanishes identically on the exact lower compensated branch.

### 9.2 Adiabatic-wall ground-state closure preserves the branch

On the adiabatic-wall ground-state branch, the quartic anomaly becomes a pure core/outgoing retuning. The parent compensation ratios remain fixed:
\[
 \delta\ln\mathfrak g=0,
 \qquad
 \delta\ln\mathfrak r=0,
 \qquad
 \delta\ln r_c=0,
 \qquad
 \delta\ln L_W=0.
\]
So the anomaly is not motion of the parent compensation sheet or the even geometry. It is just a small odd-outlet retuning on top of a slightly softened loading share.

### 9.3 Pure odd normalized DtN representative

On the adiabatic even-frozen slice the exact anomaly has a pure odd normalized DtN representative:
\[
 \boxed{
 \Sigma_5(\ell)=\frac{S}{9}(e^{-\ell}-1),
 \qquad
 \beta=1,
 \qquad
 \Sigma_0=0,
 }
\]
with
\[
 \ell = \ln(1+\Lambda_1 f).
\]
In canonical normalized gauge,
\[
 \Sigma_5 = -\frac{\Lambda_1 f}{9(1+\Lambda_1 f)}.
\]
So at the normalized outgoing-DtN level the whole electron anomaly is one exact pure-odd isotropic coefficient.

### 9.4 Reciprocal grouped-normalization law

On the natural source-map branch,
\[
 \chi_Q=e^{-\ell}
 \quad\Longleftrightarrow\quad
 N_Q=e^{\ell}
 \quad\Longleftrightarrow\quad
 P_0=e^{\ell}P_0^{\rm target}.
\]
And the same rescaling can be read as motion *along* the universal target curve
\[
 P_0^{\rm target}(a,c_s)=\frac{54Gc_s^5}{5a^5c^5},
\]
because
\[
 P_0(a,c_s e^{\ell/5}) = e^{\ell} P_0(a,c_s)
\]
with `a` fixed. So the anomaly does not demand a new normalization law. It corresponds to a very small retuning along the same law.

### 9.5 But the natural compensated branch gives zero first-order defect

The package then sharpens the status once more. On the natural compensated parent family, the last reduced first-order defect collapses to one odd mixed-channel renormalization and then to one similarity-slippage scalar. On that exact parent family,
\[
 \Xi_{\rm slip}=0,
 \qquad
 \delta\mathfrak B_W=0,
 \qquad
 \delta\gamma_W=0,
 \qquad
 \Delta_Q=0
\]
at first order.

Pure grouped real `P2` anisotropy cannot rescue this because its first scalar invariant is quadratic, not linear. So linear `P2` anisotropy has no scalar feed-down into the outgoing defect.

### 9.6 Strongest honest current verdict

The current reduced derivation supports the following strongest honest conclusion.

1. The **canonical outgoing background** is naturally fixed with no tuning.
2. The **carried local anomaly value** is naturally retained.
3. The **exact electron sliver** is **not** naturally forced at first order on the no-tuning adiabatic compensated branch.
4. Matching the exact electron value would still require one real microscopic departure from that branch, such as:
   - a genuine off-family scalar slippage
     \[
     \varepsilon_\perp \neq 0,
     \]
   - a direct odd mixed-port renormalization
     \[
     \delta\gamma_W\neq 0,
     \]
   - or a beyond-first-order effect outside the present reduced no-tuning closure.

In short:
\[
 \boxed{\text{No-tuning canonical branch: yes.}}
\]
\[
 \boxed{\text{Exact electron-point sliver: not yet forced.}}
\]

---

## 10) Short list of the most important formulas to remember

If only a few formulas survive into a later session, these are the ones.

### 10.1 Exact current local baseline
\[
 \frac{g_{\rm loc}}{2}
 =
 1+f-\frac{47}{36}f^2+c_{3,\rm total}f^3+a_{4,\rm staggered}f^4+O(f^5).
\]

### 10.2 Quartic residual tangent
\[
 \boxed{\Lambda_1\approx 0.279605891931464.}
\]

### 10.3 Coherent defect scalar
\[
 \boxed{\Xi_1=A_{\rm tr}q_{\rm tr}+q_{\rm nt}=\delta\ln\mathcal T^2=\frac{P_1}{P_0}.}
\]

### 10.4 Mixed/outgoing microscopic law
\[
 \boxed{
 \Lambda_1=-q_\Lambda+q_Z+2q_\chi+2\frac{\epsilon_*}{1-\epsilon_*}q_\epsilon.
 }
\]

### 10.5 Forced conservative grouped-`P2` module
\[
 \boxed{\widehat Y_Q^{\rm cons}(\omega)=\frac34+\frac{1}{4(1-\omega^2/\Omega_Q^2)}.}
\]

### 10.6 Support/source corollary
\[
 \boxed{\rho_\alpha=\frac43,\qquad \zeta_{\rm req}=\frac13.}
\]

### 10.7 Canonical outgoing branch
\[
 \boxed{\chi_Q=1,\qquad N_Q=1.}
\]

### 10.8 Electron-point outgoing defect
\[
 \boxed{
 \Delta_Q^{(e)}=-\frac{\Lambda_1 f}{1+\Lambda_1 f}
 \approx -3.24631584151692\times10^{-4}.
 }
\]

### 10.9 Exact isotropic DtN anomaly surface
\[
 \boxed{3(S\beta^5+9\Sigma_5)(1+\Lambda_1 f)=3S-\Sigma_0.}
\]

### 10.10 Tangent isotropic branch-selection law
\[
 \boxed{5b+\frac{a_0}{3}+9a_5=-\Lambda_1.}
\]

### 10.11 Surviving compensated outlet class
\[
 \boxed{
 \rho_R=4\sigma_W,
 \qquad
 \kappa_W=\frac13,
 \qquad
 \gamma_W=\frac{\sigma_W+\Lambda_1 f}{9\sigma_W(1+\Lambda_1 f)}.
 }
\]

### 10.12 Parent-family balance law
\[
 \boxed{1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2.}
\]

### 10.13 Best current prediction versus exact electron target
\[
 \boxed{g_{\rm pred}^{(\rm nat)}=g_{\rm loc}\approx 2.00231930435865.}
\]
\[
 \boxed{|g_e|-g_{\rm pred}^{(\rm nat)}\approx 2.27\times10^{-12}.}
\]

---

## 11) The correct next theorem gate

The next honest target is **not** more baseline algebra.
That part is already compressed.

The real next theorem gate is:

\[
 \boxed{
 \text{derive the actual microscopic branch-selection law for the remaining odd/outgoing defect.}
 }
\]

Operationally, that means determining whether the completed moving-throat branch gives:

1. a genuine off-family scalar slippage `\varepsilon_\perp`,
2. a direct odd mixed-port renormalization `\delta\gamma_W` or `\delta\mathfrak B_W`,
3. or no first-order slippage at all, in which case the strongest no-tuning prediction remains the carried local value `g_{\rm loc}` rather than the exact electron-point sliver.

So the bottleneck is no longer “fit the anomaly somehow.”
It is one sharply identified microscopic branch-selection datum.
