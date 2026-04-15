# Moving-Throat PDE — Stage 225: Actual Twin-Support Placement and Coherent Orbit-Lock Compiler

## Status

**Exact within the carried selected-branch support identities, the coherent local D/N placement map, and the support-blind orbit-lock packet compiler.**

This stage does **not** yet solve the completed moving-throat PDE.
It takes the exact Stage-224 selected twin-support curve and shows how an **actual coherent branch state** is placed on that curve, how the Stage-224 ranking region is therefore fixed, and how the separate orbit-lock packet is compiled directly from the physical placement variables.

The outcome is that the endgame is no longer “find another ranking theorem.”
It is:

1. extract the stationary coherent placement state,
2. place that state on the selected twin-support curve,
3. compute the weak-axisymmetric orbit packet,
4. and test the separate outgoing finish line.

---

## Purpose

Stage 224 left the same-charge selected branch as the exact one-parameter twin-support curve

\[
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma=\frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23.
\]

That solved the **parametric** support-selection problem, but it did not yet tell us where the **actual coherent moving-throat branch** sits on that curve.

Stage 225 closes that gap.

It also packages the separate coherent orbit-lock data in the same step, so that the reduced endgame becomes a concrete observable packet rather than a loose combination of branch variables.

---

## Provenance

This stage is the moving-throat Stage-225 translation of Barrier Stage 025, sharpened by the later coherent-placement packet compiler.

Conceptually it sits directly after:

- **Stage 223**, which fixed the selected support slice through
  \[
  \Pi_{\rm tr}=\frac43\,C_{\rm mix},
  \]
- **Stage 224**, which converted that slice into the exact twin-support curve and solved the primitive carrier ranking along it,
- and the later coherent placement-map continuation, which rewrote the orbit packet directly in the physical placement variables
  \[
  (\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).
  \]

So Stage 225 is the front-end compiler that turns the selected twin-support phase diagram into an actual branch test.

---

## 0. Why this stage is needed

Before this step, the selected same-charge branch was known only as a parametric curve.
That meant the primitive hierarchy from Stage 224 was still a **phase diagram** rather than an actual prediction.

But once the actual coherent branch returns the physical placement data

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta),
\]

the selected support point is fixed immediately by the coherent support variable

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
\]

That collapses the remaining support-selection ambiguity to one realized coordinate on the Stage-224 curve.

At the same time, the coherent orbit-lock packet is carried by the support-blind observables

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta,
\]

so the two tests naturally separate:

- **selected twin-support placement**, and
- **coherent orbit lock**.

That is the exact structure this note compiles.

---

## 1. Actual selected-twin placement

On the coherent local D/N branch, the stationary placement packet is

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).
\]

The exact coherent reduction uses

\[
\boxed{
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
}
\]

The selected same-charge branch already carries the Stage-223 support demand

\[
\Pi_{\rm tr}=\frac43\,C_{\rm mix},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

Therefore the actual selected-branch coordinate on the Stage-224 twin-support curve is

\[
\varrho_{\rm phys}
:=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon)}{\pi^2}
=
\boxed{
\frac23(1-\epsilon).
}
\]

Hence the actual selected point is

\[
\boxed{
\epsilon_{*,\rm phys}=\epsilon,
\qquad
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}
=
\frac{4}{3\varrho_{\rm phys}}-2
=
\frac{2\epsilon}{1-\epsilon}.
}
\]

So the Stage-224 support curve is no longer parametric once the completed PDE returns the coherent support pair \((\delta_U,\epsilon_W)\).

---

## 2. Exact threshold rewrite in the realized variable \(\epsilon\)

Stage 224 found the two surviving primitive crossover thresholds

\[
\varrho_{W\Lambda}
=
\frac{2(1+\beta^2)}{3(2+\beta^2)},
\qquad
\varrho_{U\Lambda}
=
\frac{2(1+\beta^2)}{3(1+\beta+\beta^2)}.
\]

Using
\[
\epsilon_*=1-\frac32\varrho,
\]
these rewrite exactly as

\[
\boxed{
\epsilon_{W\Lambda}=\frac{1}{2+\beta^2},
\qquad
\epsilon_{U\Lambda}=\frac{\beta}{1+\beta+\beta^2}.
}
\]

So the realized Stage-224 ranking region is decided directly by the actual coherent support variable \(\epsilon\):

\[
\epsilon>\epsilon_{W\Lambda}
\iff
0<\varrho<\varrho_{W\Lambda},
\]

\[
\epsilon_{U\Lambda}<\epsilon<\epsilon_{W\Lambda}
\iff
\varrho_{W\Lambda}<\varrho<\varrho_{U\Lambda},
\]

\[
0<\epsilon<\epsilon_{U\Lambda}
\iff
\varrho_{U\Lambda}<\varrho<\frac23.
\]

So once the coherent branch returns \(\epsilon\), the Stage-224 primitive ranking is fixed with no additional support solve.

---

## 3. Support-lane classifier on the realized branch

The coherent support classifier is still

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\text{mixed-only enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\text{lowest symmetric twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\text{non-twin asymmetry required}.
\]

But on the selected same-charge branch one has

\[
\Pi_{\rm tr}=\frac43\,C_{\rm mix},
\]

so the realized support slice is **automatically** in the lowest symmetric twin window:

\[
\boxed{
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}.
}
\]

This is exactly the support-side meaning of the selected branch:

- the demand-side ratio from Stage 223 is fixed to
  \[
  \frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43,
  \]
- the support window is therefore the lowest symmetric twin sector,
- and the physical coherent support variable \(\zeta\) is a **separate** support quantity, not another name for the Stage-223 demand ratio.

Important notation firewall:

- the Stage-223 reduced demand-side quantity \(\zeta_{\rm req}=1/3\) belongs to the contact/pole loading problem,
- the physical coherent local D/N support ratio \(\zeta\) is a separate support variable,
- and on the lowest symmetric twin branch the physical coherent support value is \(\zeta=1\).

These must not be conflated.

---

## 4. Exact coherent placement map and reduced-state bridge

The direct coherent-branch observables are

\[
\boxed{
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
}
\]

\[
\boxed{
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
}
\]

together with \(\epsilon_\eta\).

A useful exact bridge is

\[
\Lambda_0:=\frac{27\pi^2Gc_s^5}{20a^5c^5},
\qquad
\widehat Z_W = Z_W\frac{\Lambda_0}{\Lambda}.
\]

Then

\[
R_{\rm target}
=
\Lambda_0\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{\widehat Z_W(1+\chi_0)^2},
\]

so the physical placement-map version and the older reduced-state version are exactly equivalent.

Most importantly, the support ratio \(\zeta\) does **not** enter

\[
\widehat Z_W,\qquad
\epsilon,\qquad
R_{\rm tr},\qquad
R_{\rm target}.
\]

So the coherent support ratio is absent from the extracted orbit state.

This is the first exact support-blindness theorem of the stage.

---

## 5. Exact finite orbit packet and infinitesimal compiler

The finite coherent orbit packet is the exact quotient packet in the physical placement variables.

Let the reference branch state be

\[
(\chi_{0,\rm ref},\delta_{U,\rm ref},Z_{W,\rm ref},\epsilon_{W,\rm ref},
\epsilon_{\eta,\rm ref},\Lambda_{\rm ref}),
\]

and keep the carried coherent coefficients \(E_*,F_*\).
Then

\[
\boxed{
q_{\rm tr}
=
(1+\delta_{U,*})\ln\frac{\chi_0}{\chi_{0,\rm ref}}
+
(1+\chi_{0,*})\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]

\[
\boxed{
q_{\rm nt}
=
\ln\frac{Z_W}{Z_{W,\rm ref}}
-
\ln\frac{\Lambda}{\Lambda_{\rm ref}}
+
E_*\ln\frac{\epsilon_W}{\epsilon_{W,\rm ref}}
-
F_*\ln\frac{\delta_U}{\delta_{U,\rm ref}},
}
\]

\[
\boxed{
q_\eta = \ln\frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
}
\]

So the only change from the older reduced-state language is the exact substitution

\[
\ln \widehat Z_W = \ln Z_W - \ln\Lambda + \ln\Lambda_0,
\]

with the constant \(\Lambda_0\) dropping out of the quotient packet.

### 5.1 Infinitesimal coherent packet

Let the infinitesimal coherent placement drifts be

\[
(d\ln\chi_0,\,
 d\ln\delta_U,\,
 d\ln Z_W,\,
 d\ln\epsilon_W,\,
 d\ln\epsilon_\eta,\,
 d\ln\Lambda).
\]

Then

\[
\boxed{
\Sigma_{\rm tr}
=
(1+\delta_{U,*})\,d\ln\chi_0
+
(1+\chi_{0,*})\,d\ln\delta_U,
}
\]

\[
\boxed{
\Sigma_{\rm nt}
=
(d\ln Z_W-d\ln\Lambda)
+
E_*\,d\ln\epsilon_W
-
F_*\,d\ln\delta_U,
}
\]

\[
\boxed{
\Sigma_\eta = d\ln\epsilon_\eta.
}
\]

The physical defect packet is the exact triangular normal form

\[
\Theta_1=-C_{{\rm tr},*}\Sigma_{\rm tr},
\qquad
\Xi_1=A_{{\rm tr},*}\Sigma_{\rm tr}+\Sigma_{\rm nt},
\qquad
\mathcal R_1=-\Xi_1-\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\Sigma_\eta,
\]

with the carried coherent coefficients

\[
C_{{\rm tr},*}
=
\frac{\chi_{0,*}\delta_{U,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})(1+\chi_{0,*}+\delta_{U,*})},
\]

\[
A_{{\rm tr},*}
=
\frac{2\chi_{0,*}}
{(1+\chi_{0,*})(1+\delta_{U,*})}.
\]

### 5.2 Direct-observable drift form

The strongest exact bridge in this stage is the direct-observable compiler

\[
\boxed{
\Theta_1=d\ln R_{\rm tr},
}
\]

\[
\boxed{
d\ln\epsilon
=
d\ln\epsilon_W
-
\frac{2\delta_U}{(1+\delta_U)(11+9\delta_U)}\,d\ln\delta_U,
}
\]

\[
\boxed{
d\ln R_{\rm target}
=
d\ln\Lambda-d\ln Z_W
-
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-
\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
}
\]

and therefore

\[
\boxed{
\Xi_1
=
-d\ln R_{\rm target}
-
\frac{\epsilon_{\eta,*}}{1-\epsilon_{\eta,*}}\,d\ln\epsilon_\eta
=
-d\ln\Lambda+d\ln Z_W
+
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
+
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon.
}
\]

So \(\Xi_1\) is explicitly blind to both the support lane \(\zeta\) and the selected-branch demand variable \(\epsilon_\eta\).

---

## 6. Exact two-packet split and the coherent orbit-lock theorem

The coherent branch now separates exactly into two packets.

### 6.1 Orbit-lock packet

Finite packet:

\[
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

Infinitesimal packet:

\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

This packet depends only on

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda),
\]

and is exactly blind to \(\zeta\).

### 6.2 Support / normalization packet

The support packet is

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]

\[
M_{\rm mix}
=
\frac{8Z_W(1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\]

\[
S(\zeta;\epsilon)
=
1+
\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]

\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2}.
\]

The mixed-only product law is

\[
\boxed{
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon)}{\pi^2}.
}
\]

The exact separation statement is:

- \(\zeta\) changes only the available support baseline through \(S(\zeta;\epsilon)\),
- \(\zeta\) leaves the finite orbit packet unchanged,
- \(\zeta\) leaves the infinitesimal orbit defect packet unchanged.

So support compensation cannot rescue orbit-lock failure, and orbit lock does not determine support enhancement.

### 6.3 Coherent orbit-lock theorem

The actual coherent local D/N branch satisfies orbit lock iff

\[
\boxed{
q_{\rm tr}=q_{\rm nt}=q_\eta=0,
}
\]

equivalently iff

\[
\boxed{
d\ln R_{\rm tr}=0,
\qquad
d\ln R_{\rm target}=0,
\qquad
d\ln\epsilon_\eta=0.
}
\]

The outgoing finish line remains separate:

\[
\boxed{
N_Q=1,
\qquad
\text{or equivalently }\chi_Q=1
\text{ on the natural source-map branch.}
}
\]

So the actual remaining moving-throat endgame is exactly:

1. extract the stationary coherent placement state,
2. compute the weak-axisymmetric tangent,
3. place the branch on the selected twin-support curve via \(\epsilon\),
4. evaluate
   \[
   d\ln R_{\rm tr},\qquad
   d\ln R_{\rm target},\qquad
   d\ln\epsilon_\eta,\qquad
   N_Q-1.
   \]

---

## 7. Practical Stage-225 output packet

The smallest exact front-end packet to return from the completed moving-throat branch is

\[
\boxed{
\Bigl(
\epsilon,\,
\varrho_{\rm phys},\,
\sigma_{\rm phys},\,
\text{ranking region},\,
R_{\rm tr},\,
R_{\rm target},\,
\epsilon_\eta,\,
d\ln R_{\rm tr},\,
d\ln R_{\rm target},\,
d\ln\epsilon_\eta,\,
N_Q-1
\Bigr).
}
\]

If the support-compensation lane is being tracked separately, the smallest exact support packet is

\[
\boxed{
\Bigl(
\zeta,\,
M_{\rm mix},\,
S(\zeta;\epsilon),\,
M_{\rm tr}
\Bigr).
}
\]

So Stage 225 turns the old support phase diagram and orbit-lock language into a concrete branch-evaluation checklist.

---

## 8. What this stage achieves

Stage 225 closes three specific gaps left open by Stage 224.

### 8.1 It places the actual branch on the Stage-224 support curve

The selected same-charge point is no longer a free parameter once the coherent branch returns \((\delta_U,\epsilon_W)\).
It is fixed by

\[
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{2\epsilon}{1-\epsilon}.
\]

### 8.2 It compiles the orbit packet directly from physical placement variables

The orbit-lock packet no longer needs the full microscopic kernel tuple.
It is already determined by

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda),
\]

with exact finite and infinitesimal compilers.

### 8.3 It rigorously separates orbit lock from support compensation

The support variable \(\zeta\) affects only

\[
S(\zeta;\epsilon),\qquad M_{\rm tr},
\]

and does **not** enter

\[
R_{\rm tr},\qquad R_{\rm target},\qquad \epsilon_\eta,
\qquad
(q_{\rm tr},q_{\rm nt},q_\eta),
\qquad
(\Theta_1,\Xi_1,\mathcal R_1).
\]

So the completed PDE has to satisfy the orbit-lock packet and the support packet separately on the same physical branch.

---

## 9. SymPy-backed status

The accompanying audit script verifies the exact algebra used here:

- the actual selected-branch placement identity
  \[
  \varrho_{\rm phys}=\frac23(1-\epsilon),
  \]
- the twin-support rewrite
  \[
  \sigma_{\rm phys}=\frac{2\epsilon}{1-\epsilon}
  =\frac{4}{3\varrho_{\rm phys}}-2,
  \]
- the threshold rewrite
  \[
  \epsilon_{W\Lambda}=\frac{1}{2+\beta^2},
  \qquad
  \epsilon_{U\Lambda}=\frac{\beta}{1+\beta+\beta^2},
  \]
- strict inclusion of the selected branch in the lowest symmetric twin window,
- the reduced-state bridge
  \[
  \widehat Z_W = Z_W\Lambda_0/\Lambda,
  \qquad
  R_{\rm target}
  =
  \Lambda_0\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{\widehat Z_W(1+\chi_0)^2},
  \]
- exact support-blindness of
  \[
  \epsilon,\quad R_{\rm tr},\quad R_{\rm target},
  \quad q_{\rm tr},\quad q_{\rm nt},\quad q_\eta,
  \]
  with respect to \(\zeta\),
- the infinitesimal compilers for
  \[
  d\ln\epsilon,\quad d\ln R_{\rm tr},\quad d\ln R_{\rm target},
  \]
- the direct-observable identity for \(\Xi_1\),
- and the exact two-packet split
  \[
  \partial_\zeta R_{\rm tr}
  =
  \partial_\zeta R_{\rm target}
  =
  \partial_\zeta \epsilon
  = 0,
  \qquad
  \partial_\zeta M_{\rm tr}
  =
  M_{\rm mix}\frac{1-\epsilon}{(1-\zeta\epsilon)^2}.
  \]

Supporting file:
- `moving_throat_pde_stage225_actual_twin_support_placement_and_coherent_orbit_lock_compiler_sympy_audit.py`
