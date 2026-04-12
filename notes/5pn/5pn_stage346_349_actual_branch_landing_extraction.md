# 5PN continuation notes — Stages 346–349

These stages stop trying to compress the reduced theorem any further and instead extract the **actual coherent-branch finish packet** as far as the present file stack allows.

The goal was the one left at Stage 345:

\[
(R_{\rm tr},\ R_{\rm target},\ \epsilon_\eta,\ N_Q)
\]

for the actual coherent local D/N branch, and then the exact surfaces on which the four reduced finish-line conditionals vanish.

The main outcome is sharp:

1. the current notes/scripts do fix the **actual coherent-branch formulas** for
   \[
   R_{\rm tr},\quad R_{\rm target},\quad \epsilon_\eta,
   \]
   and the natural-source-map relation
   \[
   N_Q=\chi_Q^{-1},
   \]
2. the first three finish-line conditions reduce to exact codrift surfaces in the physical branch variables,
3. the outgoing condition is the exact canonical branch condition
   \[
   \chi_Q=1
   \quad\Longleftrightarrow\quad
   N_Q=1,
   \]
   with the lower-parent compensation family giving the linear sufficient condition
   \[
   \Xi_{\rm slip}\,\delta\Pi_{\rm tan}=0,
   \]
4. but the file stack still does **not** contain a completed numerical PDE-selected branch point.  
   It contains the exact symbolic branch packet and the exact landing surfaces.

So the reduced-theorem side is finished, but the PDE-realization side is still open.

---

## Stage 346 — actual coherent local D/N branch values

On the actual coherent tracking branch,

\[
R_{\rm tr}
=
\frac{1+\chi_0/(1+\delta_U)}{1+\chi_0}
=
\frac{\chi_0+\delta_U+1}{(1+\chi_0)(1+\delta_U)},
\]

\[
\epsilon
=
\epsilon_W^{(\mathrm{split})}
=
\epsilon_W\!\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right)
=
\epsilon_W\frac{11+9\delta_U}{11(1+\delta_U)},
\]

\[
R_{\rm target}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

\[
M_{\rm mix}
=
\frac{8 Z_W (1+\chi_0)^2}{\pi^2(1-\epsilon_\eta)(1-\epsilon)},
\qquad
M_{\rm tr}=M_{\rm mix}S(\zeta;\epsilon),
\]

\[
S(\zeta;\epsilon)=1+\frac{\zeta(1-\epsilon)}{1-\zeta\epsilon},
\]

and on the natural source-map branch,

\[
N_Q=\frac{1}{\chi_Q}.
\]

The exact support-blindness check is explicit:
\[
\partial_\zeta R_{\rm tr}=0,
\qquad
\partial_\zeta \epsilon_\eta=0,
\qquad
\partial_\zeta R_{\rm target}=0.
\]

So the coherent support lane changes only the baseline through \(M_{\rm tr}=M_{\rm mix}S\); it does not move the orbit packet.

There is also an exact product law:
\[
R_{\rm target}M_{\rm mix}
=
\frac{8\Lambda(1-\epsilon)}{\pi^2}.
\]

---

## Stage 347 — the first three finish-line conditions in physical branch variables

Using logarithmic branch drifts
\[
d\ln\chi_0,\quad d\ln\delta_U,\quad d\ln\epsilon_W,\quad d\ln\epsilon_\eta,\quad d\ln Z_W,\quad d\ln\Lambda,
\]
the exact split-blocking drift is
\[
d\ln\epsilon
=
d\ln\epsilon_W
-
\frac{2\delta_U}{(1+\delta_U)(11+9\delta_U)}\,d\ln\delta_U.
\]

The three reduced finish-line conditionals are

\[
d\ln R_{\rm tr}
=
-\frac{\chi_0\delta_U}{(1+\chi_0)(1+\delta_U)(1+\chi_0+\delta_U)}
\Big[(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U\Big],
\]

\[
d\ln R_{\rm target}
=
d\ln\Lambda-d\ln Z_W
-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
\]

\[
d\ln\epsilon_\eta=d\ln\epsilon_\eta.
\]

So the exact orbit-lock / target-ratio / dressing landing surfaces are

\[
(1+\delta_U)d\ln\chi_0+(1+\chi_0)d\ln\delta_U=0,
\]

\[
d\ln\Lambda-d\ln Z_W
-\frac{\epsilon_\eta}{1-\epsilon_\eta}d\ln\epsilon_\eta
-\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
-\frac{2\epsilon}{1-\epsilon}d\ln\epsilon=0,
\]

\[
d\ln\epsilon_\eta=0.
\]

After imposing the first and third, the second becomes the actual coherent-branch codrift law
\[
d\ln\Lambda
=
d\ln Z_W
+\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
+\frac{2\epsilon}{1-\epsilon}d\ln\epsilon.
\]

---

## Stage 348 — the outgoing normalization finish surface

The exact outgoing condition on the natural source-map branch is
\[
N_Q=\chi_Q^{-1},
\]
so the fourth finish-line condition is simply
\[
N_Q=1
\quad\Longleftrightarrow\quad
\chi_Q=1.
\]

On the canonical compact passive/outgoing branch this is exactly the Stage-87/89 condition.

On the lower-parent compensation family at first order,
\[
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan},
\qquad
N_Q=\frac{1}{1+\Delta_Q},
\]
so either
\[
\Xi_{\rm slip}=0
\qquad\text{or}\qquad
\delta\Pi_{\rm tan}=0
\]
is sufficient to land on
\[
N_Q=1.
\]

So the outgoing side is not another orbit condition. It is a separate exact surface.

---

## Stage 349 — actual four-condition extractor and honest verdict

The actual coherent branch packet is therefore

\[
R_{\rm tr}
=
\frac{\chi_0+\delta_U+1}{(1+\chi_0)(1+\delta_U)},
\qquad
R_{\rm target}
=
\Lambda\,
\frac{(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\]

\[
\epsilon_\eta=\epsilon_\eta,
\qquad
N_Q=\chi_Q^{-1}.
\]

The exact four finish-line conditionals are

1. \(d\ln R_{\rm tr}=0\),
2. \(d\ln R_{\rm target}=0\),
3. \(d\ln\epsilon_\eta=0\),
4. \(N_Q-1=0\).

Their exact combined landing compiler is:

- orbit lock
  \[
  d\ln\delta_U
  =
  -\frac{1+\delta_U}{1+\chi_0}\,d\ln\chi_0,
  \]
- dressing lock
  \[
  d\ln\epsilon_\eta=0,
  \]
- target-ratio lock
  \[
  d\ln\Lambda
  =
  d\ln Z_W
  +\frac{2\chi_0}{1+\chi_0}d\ln\chi_0
  +\frac{2\epsilon}{1-\epsilon}d\ln\epsilon,
  \]
- outgoing lock
  \[
  \chi_Q=1
  \quad\text{(equivalently \(N_Q=1\)).}
  \]

This is the sharpest actual-branch extraction the present file stack supports.

### Honest end-state

What the current files **do** contain:

- the actual coherent-branch symbolic packet,
- the exact support-blind orbit packet,
- the exact outgoing normalization surface,
- and the exact combined four-condition landing surfaces.

What they **do not** yet contain:

- a completed numerical or closed-form PDE-selected point
  \[
  (R_{\rm tr}^{\rm phys},R_{\rm target}^{\rm phys},\epsilon_\eta^{\rm phys},N_Q^{\rm phys}).
  \]

So the reduced program has reached the point where **all remaining uncertainty is branch realization**.

That means the next real move is no longer another reduced compiler.
It is to solve or numerically locate the completed moving-throat branch strongly enough to read off those four actual values.
