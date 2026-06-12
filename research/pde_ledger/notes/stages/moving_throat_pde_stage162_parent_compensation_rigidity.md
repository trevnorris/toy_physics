# Moving-Throat PDE — Stage 162: Parent Compensation-Surface Rigidity and Automatic Similarity Preservation

## Goal

Stage 161 reduced the entire first-order co-evolving Family-1 normalization problem to one microscopic scalar,
\[
\Xi_{\rm slip}=\Xi_\gamma-2\Xi_L,
\qquad
\Delta_Q
=
-\frac{\sigma_*}{1-\sigma_*}\,\Xi_{\rm slip}\,\delta\Pi_{\rm tan}.
\]
So the remaining question was whether the actual tangential co-evolving mouth deformation preserves the compensated D/N similarity law
\[
\gamma_0=\frac{4L_W^2}{3\pi^2 a^2}
\]
at first order, or produces a nonzero similarity-slippage defect.

The exact parent formulas from Stage 119 sharpen that question further.
They imply two stronger statements.

1. **If the actual branch stays on the exact parent compensation family, then D/N similarity preservation is automatic:**
   \[
   \Xi_{\rm slip}=0
   \]
   identically along the family.
2. **If, in addition, the carried canonical-even gate forces \(\delta\mathfrak g=0\), then the lower compensated parent family is first-order rigid:**
   \[
   \delta\mathfrak r=0.
   \]

So the first-order reduced 2.5PN defect vanishes automatically on any co-evolving tangential deformation that remains inside the exact parent compensation family.

This is stronger than Stage 161, because it no longer treats D/N similarity preservation as an extra assumption. Inside the exact parent family, it is an identity.

---

## 1. Exact parent compensation family

Stage 119 rewrote the compensated throat-core condition in terms of the normalized parent ratios
\[
\mathfrak r:=\frac{\lambda}{\sqrt{K_sK_q}},
\qquad
\mathfrak g:=\frac{g_q\sqrt{K_s}}{g_s\sqrt{K_q}},
\]
with exact balance law
\[
1+\mathfrak r^2=4(\mathfrak g-\mathfrak r)^2,
\qquad
\mathfrak g_\pm(\mathfrak r)=\mathfrak r\pm\frac12\sqrt{1+\mathfrak r^2}.
\]
On the same parent family, the D/N mixed-tube selection law is
\[
\boxed{
\frac{L_W}{a}=\frac{\pi}{2}\sqrt{\frac{1+\mathfrak r^2}{3}}.
}
\]
Stages 115–116 give the conditional bare odd normalization formula inside the pure-scale realization posited at Stage 116:
\[
\boxed{
\gamma_0=\frac{1+\mathfrak r^2}{9}.
}
\]
Thus the even D/N scale ratio is the genuinely derived companion relation, while the odd bare normalization value is a branch-determinable pure-scale realization posit controlled by the **same single parent variable** \(\mathfrak r\), not a separate derived parent theorem.

---

## 2. Exact automatic similarity preservation on the parent family

Differentiate the two parent-family formulas, treating the \(\gamma_0\) formula as the exact conditional consequence of the posited pure-scale realization:
\[
\ln\gamma_0
=
\ln(1+\mathfrak r^2)-\ln 9,
\]
\[
\ln\left(\frac{L_W}{a}\right)
=
\ln\frac{\pi}{2}-\frac12\ln 3+\frac12\ln(1+\mathfrak r^2).
\]
Therefore
\[
\delta\ln\gamma_0
=
\frac{2\mathfrak r}{1+\mathfrak r^2}\,\delta\mathfrak r,
\qquad
\delta\ln\left(\frac{L_W}{a}\right)
=
\frac{\mathfrak r}{1+\mathfrak r^2}\,\delta\mathfrak r.
\]
Subtracting,
\[
\boxed{
\delta\ln\gamma_0
-
2\,\delta\ln\left(\frac{L_W}{a}\right)
=0.
}
\]
So the Stage 161 similarity-slippage scalar vanishes identically on the exact parent family:
\[
\boxed{
\Xi_{\rm slip}=0.
}
\]
This result does **not** require the extra Stage 159 even-preservation gate.
It is true for any infinitesimal motion that stays on the exact parent compensation family.

Equivalently, inside that posited pure-scale realization, the bare mixed side-channel follows the exact formula-level deformation of the canonical compact outgoing branch all along that family.

---

## 3. Lower compensated branch rigidity under the canonical-even gate

Now add the carried canonical-even result from the co-evolving Family-1 outlet analysis:
\[
\delta\mathfrak g=0.
\]
On the lower compensated parent branch,
\[
\mathfrak g_-(\mathfrak r)
=
\mathfrak r-\frac12\sqrt{1+\mathfrak r^2},
\]
so
\[
\delta\mathfrak g
=
\left(
1-\frac{\mathfrak r}{2\sqrt{1+\mathfrak r^2}}
\right)
\delta\mathfrak r.
\]
The prefactor is exactly
\[
\boxed{
1-\frac{\mathfrak r}{2\sqrt{1+\mathfrak r^2}}
=
\frac{4+3\mathfrak r^2}{2\sqrt{1+\mathfrak r^2}\bigl(2\sqrt{1+\mathfrak r^2}+\mathfrak r\bigr)}
>0
}
\]
for every real \(\mathfrak r\). Therefore the lower compensated branch is first-order rigid:
\[
\boxed{
\delta\mathfrak g=0
\quad\Longrightarrow\quad
\delta\mathfrak r=0.
}
\]
So on the actual lower branch selected by the positive-source Family-1 throat core, the even-preservation gate does more than kill the mouth-coupling drift. It pins the normalized static/mixed hybridization itself.

On the Family-1 value
\[
\mathfrak r_{F1}\approx 1.77799353547498,
\]
one finds
\[
\left.\frac{d\mathfrak g_-}{d\mathfrak r}\right|_{\mathfrak r_{F1}}
\approx 0.564199521046343,
\qquad
\delta\mathfrak r
\approx
1.77242263188285\,\delta\mathfrak g.
\]
So the rigidity is numerically substantial, not marginal.

---

## 4. Collapse of all first-order D/N similarity defects

Once \(\delta\mathfrak r=0\), all first-order parent-family D/N quantities freeze:
\[
\delta r_c
=
2\mathfrak r_*\,\delta\mathfrak r
=0,
\]
\[
\delta\ln\left(\frac{L_W}{a}\right)=0,
\qquad
\delta\ln\gamma_0=0,
\qquad
\delta\kappa_0=0,
\qquad
\delta\gamma_0=0.
\]
Therefore the bare similarity defects of Stage 161 vanish individually,
\[
\boxed{
\delta\varepsilon_\kappa=0,
\qquad
\delta\varepsilon_\gamma=0,
\qquad
\delta\mathfrak B_W=0,
}
\]
and so does the renormalized odd outlet defect,
\[
\boxed{
\delta\gamma_W=0.
}
\]
Substituting back into the Stage 161 defect law,
\[
\boxed{
\Delta_Q=0,
\qquad
N_Q-1=0
}
\]
at first order.

So inside the exact parent compensation family, the co-evolving tangential mouth deformation is harmless at linear order.
The first-order reduced 2.5PN obstruction disappears automatically.

---

## 5. Best current theorem statement after Stage 162

Inside the explicit co-evolving Family-1 closure, together with the exact parent compensation family of Stage 119:

1. the D/N mixed-tube scale ratio is an exact function of the single parent variable \(\mathfrak r\), and the bare odd normalization follows the conditional pure-scale realization formula carried by Stage 116,
2. therefore the D/N similarity-slippage scalar vanishes identically along the full parent compensation family,
   \[
   \Xi_{\rm slip}=0,
   \]
3. and on the lower compensated branch the carried canonical-even gate \(\delta\mathfrak g=0\) forces
   \[
   \delta\mathfrak r=0,
   \]
   so the co-evolving tangential deformation is actually tangent to a fixed-\(\mathfrak r\) slice of that family,
4. hence
   \[
   \delta\mathfrak B_W=0,
   \qquad
   \delta\gamma_W=0,
   \qquad
   \Delta_Q=0,
   \qquad
   N_Q-1=0
   \]
   at first order.

This is the strongest positive result reached so far.

The remaining PDE-facing question is now even narrower:

> does the true moving-throat co-evolving core stay on the exact parent compensation family to first order?

If it does, the first-order reduced 2.5PN quadrupole-normalization defect vanishes automatically.
