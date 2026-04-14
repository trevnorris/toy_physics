# Same-Charge Barrier Audit — Stage 025
## Actual Twin-Support Placement and Coherent Orbit-Lock Compiler

## Purpose

Stage 024 closed the reduced support-selection algebra by restricting the selected
same-charge branch to the exact one-parameter twin-support curve

\[
\epsilon_* = 1-\frac{3\varrho}{2},
\qquad
\sigma = \frac{4}{3\varrho}-2,
\qquad
0<\varrho<\frac23.
\]

The remaining task is no longer another ranking theorem. It is to place the
**actual coherent moving-throat branch** on that curve and then evaluate the
separate orbit-lock packet.

This note compiles those two jobs into one exact front-end stage.

---

## 1. Actual selected-twin placement

On the coherent local D/N branch, the stationary placement packet is

\[
(\chi_0,\delta_U,Z_W,\epsilon_W,\epsilon_\eta,\Lambda,\zeta).
\]

The exact coherent reduction uses

\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right).
\]

The selected same-charge branch already carries the Stage-023 support demand

\[
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

Therefore the actual selected-branch coordinate on the Stage-024 twin-support
curve is

\[
\varrho_{\rm phys}
:=
\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=
\frac23(1-\epsilon).
\]

Hence

\[
\epsilon_{*,\rm phys}=\epsilon,
\qquad
\varrho_{\rm phys}=\frac23(1-\epsilon),
\qquad
\sigma_{\rm phys}=\frac{4}{3\varrho_{\rm phys}}-2=
\frac{2\epsilon}{1-\epsilon}.
\]

So the parametric Stage-024 curve is converted into a realized support point as
soon as the completed PDE returns the coherent packet
\((\delta_U,\epsilon_W)\).

---

## 2. Exact threshold rewrite in the realized variable \(\epsilon\)

Stage 024 found the two primitive thresholds

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
\(\epsilon_*=1-\tfrac32\varrho\),
these become

\[
\epsilon_{W\Lambda}=\frac{1}{2+\beta^2},
\qquad
\epsilon_{U\Lambda}=\frac{\beta}{1+\beta+\beta^2}.
\]

So the realized primitive ranking is decided directly by the actual coherent
support variable \(\epsilon\):

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

On the selected same-charge branch,

\[
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\]

so the realized support slice is **automatically** in the lowest symmetric twin
window:

\[
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}.
\]

Important notation point:

- the Stage-023 demand-side ratio \(\zeta_{\rm req}=1/3\) belongs to the reduced
  contact/pole loading problem,
- the physical coherent local D/N support ratio \(\zeta\) is a separate support
  variable,
- and on the lowest symmetric twin branch the physical coherent support value is
  \(\zeta=1\).

These two uses of \(\zeta\) must not be conflated.

---

## 4. Exact orbit-side observables

The coherent orbit packet is carried by

\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0},
\]

\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

together with \(\epsilon_\eta\).

The support ratio \(\zeta\) does **not** enter
\(\epsilon\),
\(R_{\rm tr}\),
or
\(R_{\rm target}\).
So the orbit packet is exactly support-blind.

---

## 5. Infinitesimal compiler

Let the actual coherent placement drifts be

\[
(d\ln\chi_0,
 d\ln\delta_U,
 d\ln Z_W,
 d\ln\epsilon_W,
 d\ln\epsilon_\eta,
 d\ln\Lambda).
\]

Then

\[
d\ln\epsilon
=
d\ln\epsilon_W
-
\frac{2\delta_U}{(1+\delta_U)(11+9\delta_U)}d\ln\delta_U.
\]

Also

\[
d\ln R_{\rm tr}
=
-
\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\Big[(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U\Big].
\]

And

\[
d\ln R_{\rm target}
=
d\ln\Lambda-d\ln Z_W
-
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-
\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon.
\]

The direct observable weak-axisymmetric defect packet is therefore

\[
\Theta_1=d\ln R_{\rm tr},
\]

\[
\Xi_1
=
-d\ln R_{\rm target}-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
=
-d\ln\Lambda+d\ln Z_W+
\frac{2\chi_0}{1+\chi_0}d\ln\chi_0+
\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]

\[
\mathcal R_1=d\ln R_{\rm target}.
\]

So \(\Xi_1\) is explicitly blind to both the support lane \(\zeta\) and the
selected-branch demand variable \(\epsilon_\eta\).

---

## 6. Exact coherent orbit-lock theorem gate

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
d\ln\epsilon_\eta=0.
\]

The outgoing finish line remains separate:

\[
N_Q=1,
\qquad
\text{or equivalently }\chi_Q=1
\text{ on the natural source-map branch.}
\]

So the actual remaining 5PN / moving-throat endgame is exactly:

1. extract the stationary coherent placement state,
2. compute the weak-axisymmetric tangent,
3. place the branch on the selected twin-support curve via \(\epsilon\),
4. evaluate
   \(d\ln R_{\rm tr}\),
   \(d\ln R_{\rm target}\),
   \(d\ln\epsilon_\eta\),
   and \(N_Q-1\).

---

## 7. Practical Stage-025 output packet

The smallest exact compiler packet to return from the actual moving-throat branch is

\[
\Bigl(
\epsilon,
\varrho_{\rm phys},
\sigma_{\rm phys},
\text{ranking region},
R_{\rm tr},
R_{\rm target},
\epsilon_\eta,
 d\ln R_{\rm tr},
 d\ln R_{\rm target},
 d\ln\epsilon_\eta,
 N_Q-1
\Bigr).
\]

That is the sharp front edge after Stage 024.
