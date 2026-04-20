
# Moving-Throat PDE — Stage 223: Selected-Branch Loading Ratio from the Minimal Isotropic Quadrupole Precursor

## Status

**Exact within the carried selected-branch normalization identities and the minimal isotropic one-pole conservative quadrupole precursor.**

This stage does **not** solve the full moving-throat PDE branch.
It closes a narrower reduced question:

> once the selected passive/outgoing quadrupole branch is normalized, what exact support loading ratio does the natural minimal isotropic conservative precursor pick?

The answer is completely sharp:

\[
\boxed{
\rho_\alpha=\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}=\frac43,
\qquad
\zeta_{\rm req}=\rho_\alpha-1=\frac13,
\qquad
\Pi_{\rm tr}=\frac43\,C_{\rm mix}.
}
\]

Equivalently,

\[
\boxed{
\varrho:=\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
=\frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}:=\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
}
\]

So the natural minimal isotropic passive/outgoing branch is **not** mixed-only and **not** non-twin.
It lands exactly on the **symmetric lowest-twin** support slice.

---

## Purpose

Stage 222 completed the rigid-mouth physical normal form and showed that the rigid-mouth actual branch is already a Cartesian product of

1. a transfer-shape axis, and
2. a dressing axis.

But the same-charge barrier / support-selection side of the program still carried one scalar ambiguity:

\[
\varrho=\frac{\pi^2\Pi_{\rm tr}}{16\Lambda}.
\]

What was still missing was the actual value selected by the passive/outgoing normalization side.
This stage fixes that value exactly by combining two already-carried ingredients:

1. the selected-branch product identities, which cancel all separate outgoing-normalization amplitudes and reduce the support test to a single loading ratio,
2. the minimal isotropic conservative quadrupole precursor, whose contact/pole split fixes that loading ratio algebraically.

So the real job of Stage 223 is:

> compress the selected-branch support demand from an open scalar selector down to one exact support slice.

---

## Provenance

This note is the moving-throat Stage-223 version of the support-selection computation carried by the selected-branch normalization track.
The barrier-audit Stage 023 records the same result in the same notation, and the matching `g2` step carries the exact contact-plus-pole extraction that fixes the loading ratio.

---

## 0. Why this stage is needed

Before this step, the selected passive/outgoing support problem was still phrased in terms of the demand parameter

\[
\varrho=\frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\]

with three possible support regimes:

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
\]

So the missing datum was no longer a whole response function.
It was only the exact ratio \(\Pi_{\rm tr}/C_{\rm mix}\) selected by the passive/outgoing branch.

This stage shows that the ratio is not free.

---

## 1. Exact selected-branch product identities

The selected-branch normalization side gives the exact product formulas

\[
\Pi_{\rm tr}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm req},
\qquad
C_{\rm mix}
=
\frac{N_Q^{(\rm target)}}{\beta_0}\,\alpha_{\rm mix}.
\]

Therefore the separate selected-mode amplitudes cancel immediately in the ratio:

\[
\boxed{
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}
=:\rho_\alpha.
}
\]

So the selected support demand depends only on the loading ratio

\[
\boxed{
\rho_\alpha:=\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}.
}
\]

In the spectral notation of the selected branch one may also write

\[
N_Q^{(\rm target)}
=
\hat m_-^{\,2}\,\beta_0\,\frac{s_-}{\lambda_-},
\]

so the same identities become

\[
\Pi_{\rm tr}
=
\hat m_-^{\,2}\frac{s_-}{\lambda_-}\alpha_{\rm req},
\qquad
C_{\rm mix}
=
\hat m_-^{\,2}\frac{s_-}{\lambda_-}\alpha_{\rm mix},
\]

and again the ratio is exactly \(\rho_\alpha\).

So after the selected-branch normalization is imposed, the support test no longer scans independent amplitudes.
It scans one loading ratio.

---

## 2. Exact contact-plus-pole inverse compiler

Write the natural minimal conservative quadrupole precursor as

\[
Y_Q^{\rm cons}(\omega)
=
c_0+\frac{c_1}{1-\omega^2/\Omega_Q^2},
\qquad
c_0+c_1=1.
\]

On the explicit support/source branch the natural contact/pole interpretation is:

- the mixed baseline contributes the static contact fraction,
- the extra support lane contributes the finite conservative pole.

So the same precursor can be written as

\[
Y_Q^{\rm cons}(\omega)
=
\frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
+
\frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
\frac{1}{1-\omega^2/\Omega_Q^2}.
\]

Introducing the loading ratio
\[
\rho_\alpha:=\frac{\alpha_{\rm req}}{\alpha_{\rm mix}},
\]
this becomes
\[
Y_Q^{\rm cons}(\omega)
=
\frac{1}{\rho_\alpha}
+
\frac{\rho_\alpha-1}{\rho_\alpha}
\frac{1}{1-\omega^2/\Omega_Q^2}.
\]

So the contact and pole fractions are exactly

\[
\boxed{
c_0=\frac{1}{\rho_\alpha},
\qquad
c_1=\frac{\rho_\alpha-1}{\rho_\alpha}.
}
\]

The inverse formulas are therefore immediate:

\[
\boxed{
\rho_\alpha=\frac{1}{c_0}=\frac{1}{1-c_1},
}
\]

and

\[
\boxed{
\zeta_{\rm req}:=\rho_\alpha-1=\frac{c_1}{c_0}.
}
\]

So the support/source loading ratio is encoded exactly in the static contact/pole split of the conservative quadrupole precursor.

### 2.1 What does not matter here

The pole location \(\Omega_Q\) controls the dynamical conservative shape of the precursor, but it does **not** enter the static loading-ratio extraction.
At this stage only the normalized contact and pole weights \((c_0,c_1)\) matter.

That is why the support-selection problem collapses so sharply.

---

## 3. Matching to the minimal isotropic quadrupole module

The minimal isotropic conservative quadrupole module is

\[
\boxed{
c_0=\frac34,
\qquad
c_1=\frac14,
\qquad
\Omega_Q=\frac{3c_s}{2a}.
}
\]

Substituting these values into the inverse compiler gives immediately

\[
\boxed{
\rho_\alpha=\frac{1}{3/4}=\frac43,
}
\]

and

\[
\boxed{
\zeta_{\rm req}=\frac{1/4}{3/4}=\frac13.
}
\]

Then the selected demand product is fixed exactly:

\[
\boxed{
\Pi_{\rm tr}
=
\rho_\alpha\,C_{\rm mix}
=
\frac43\,C_{\rm mix}.
}
\]

So the natural minimal isotropic passive/outgoing branch does **not** leave the support loading ratio open.
It fixes it once and for all.

---

## 4. Exact support-selector compiler

Now combine the selected loading ratio with the support selector

\[
\varrho:=\frac{\pi^2\Pi_{\rm tr}}{16\Lambda},
\qquad
C_{\rm mix}=\frac{8\Lambda(1-\epsilon_*)}{\pi^2}.
\]

Substituting \(\Pi_{\rm tr}=\frac43 C_{\rm mix}\) gives

\[
\varrho
=
\frac{\pi^2}{16\Lambda}\cdot\frac43\cdot\frac{8\Lambda(1-\epsilon_*)}{\pi^2}
=
\boxed{\frac{2(1-\epsilon_*)}{3}.}
\]

The required support enhancement becomes

\[
S_{\rm req}
=
\frac{\Pi_{\rm tr}}{C_{\rm mix}}
=
\boxed{\frac43.}
\]

So the selected branch is no longer scanning arbitrary demand sectors.
It is locked to one exact support ratio.

---

## 5. Regime meaning

The three support regimes were

\[
\Pi_{\rm tr}\le C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{mixed-only enough},
\]

\[
C_{\rm mix}<\Pi_{\rm tr}\le 2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{symmetric lowest twin enough},
\]

\[
\Pi_{\rm tr}>2C_{\rm mix}
\quad\Longleftrightarrow\quad
\text{non-twin asymmetry required}.
\]

But the selected branch gives

\[
\Pi_{\rm tr}=\frac43 C_{\rm mix},
\]

so one has the exact strict inequality

\[
\boxed{
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}.
}
\]

Therefore:

- mixed-only is **not** enough,
- the symmetric lowest twin **is** enough,
- non-twin asymmetry is **not** required.

So the last support ambiguity collapses from three sectors to exactly one selected support slice:

\[
\boxed{
\text{symmetric lowest twin, with }
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
}
\]

---

## 6. Exact theorem statement

The whole stage can now be read as one exact theorem.

### Selected-branch loading-ratio theorem

Assume:

1. the selected passive/outgoing normalization side satisfies
   \[
   \Pi_{\rm tr}
   =
   \frac{N_Q^{(\rm target)}}{\beta_0}\alpha_{\rm req},
   \qquad
   C_{\rm mix}
   =
   \frac{N_Q^{(\rm target)}}{\beta_0}\alpha_{\rm mix},
   \]
2. the conservative quadrupole precursor is read in the natural contact-plus-pole form
   \[
   Y_Q^{\rm cons}(\omega)
   =
   \frac{\alpha_{\rm mix}}{\alpha_{\rm req}}
   +
   \frac{\alpha_{\rm req}-\alpha_{\rm mix}}{\alpha_{\rm req}}
   \frac{1}{1-\omega^2/\Omega_Q^2},
   \]
3. the actual isotropic minimal conservative module is
   \[
   c_0=\frac34,
   \qquad
   c_1=\frac14.
   \]

Then the selected support loading is fixed exactly by

\[
\boxed{
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix},
}
\]

equivalently

\[
\boxed{
\varrho=\frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac43,
}
\]

and the selected support regime is exactly the symmetric lowest-twin slice

\[
\boxed{
C_{\rm mix}<\Pi_{\rm tr}<2C_{\rm mix}.
}
\]

So the support-selection side of the selected passive/outgoing branch is no longer an open phase diagram.
It is one exact branch value.

---

## 7. What this fixes, and what it leaves open

This stage fixes:

1. the loading ratio
   \[
   \rho_\alpha=\frac43,
   \]
2. the required support excess
   \[
   \zeta_{\rm req}=\frac13,
   \]
3. the selected support selector
   \[
   \varrho=\frac{2(1-\epsilon_*)}{3},
   \]
4. the exact regime classification:
   symmetric lowest twin is enough and non-twin asymmetry is not needed.

What remains open is not support-sector selection anymore.
The next honest question is narrower:

> once we restrict to this selected twin-support branch, how much of the remaining primitive ranking or branch-loading ambiguity survives?

So the program has moved from “which support sector?” to “what survives on the selected support curve?”

---

## 8. Best current summary after Stage 223

The selected-branch normalization side has now fixed the support ratio carried by the natural minimal isotropic passive/outgoing quadrupole branch:

\[
\boxed{
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\qquad
\Pi_{\rm tr}=\frac43 C_{\rm mix}.
}
\]

Equivalently,

\[
\boxed{
\varrho=\frac{2(1-\epsilon_*)}{3},
\qquad
S_{\rm req}=\frac43.
}
\]

So the last support ambiguity has collapsed from three sectors

- mixed-only,
- symmetric lowest twin,
- non-twin asymmetry,

to exactly one selected support slice:

\[
\text{symmetric lowest twin, with demand ratio }
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\frac43.
\]

That is the cleanest reduced loading-ratio statement available before the next primitive ranking step.

---

## 9. SymPy-backed status

The accompanying SymPy audit verifies:

- the exact selected-branch product-ratio identity
  \[
  \frac{\Pi_{\rm tr}}{C_{\rm mix}}
  =
  \frac{\alpha_{\rm req}}{\alpha_{\rm mix}}
  =\rho_\alpha,
  \]
- the contact-plus-pole compiler
  \[
  c_0=\frac{1}{\rho_\alpha},
  \qquad
  c_1=\frac{\rho_\alpha-1}{\rho_\alpha},
  \]
  together with the inverse formulas
  \[
  \rho_\alpha=\frac{1}{c_0}=\frac{1}{1-c_1},
  \qquad
  \zeta_{\rm req}=\frac{c_1}{c_0},
  \]
- the specialization of the minimal isotropic conservative module
  \[
  c_0=\frac34,
  \qquad
  c_1=\frac14,
  \]
  to
  \[
  \rho_\alpha=\frac43,
  \qquad
  \zeta_{\rm req}=\frac13,
  \]
- the exact support-selector reduction
  \[
  \varrho
  =
  \frac{\pi^2\Pi_{\rm tr}}{16\Lambda}
  =
  \frac{2(1-\epsilon_*)}{3},
  \]
- and the regime classification proving that the selected branch lies strictly in the symmetric-lowest-twin sector.

Supporting file:
- `moving_throat_pde_stage223_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_sympy_audit.py`
