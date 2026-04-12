# 5PN continuation notes — Stages 258–261 exact one-pole regime bridge back to the selected moving-throat grouped-`P2` branch

This session continued directly from the exact `\Lambda_{\rm EM}` refresh and the caution not to re-freeze the canonical `c_{\rm pole}=1/4` split too early.

The next clean move was to propagate the exact support-regime theorems back into the **generic isotropic one-pole conservative grouped-`P2` carrier** and then into the geometry-lane contamination variables `\epsilon_2,\epsilon_4`, while keeping the blocking variable `\epsilon_{\rm blk}` explicit instead of silently setting it to zero again.

The payoff is sharper than before:

1. the entire mixed / twin / non-twin support classification is exact in one-pole language,
2. the same split becomes an exact phase portrait in `\bigl(\epsilon_2,\epsilon_4\bigr)`,
3. blocking changes the branch demand with a sign flip at the universal threshold `c_{\rm pole}=1/2`,
4. and on the exact `\Lambda_{\rm EM}`-refreshed Family-1 branch the admissible non-twin corridor actually **widens** with blocking even though the actual isotropic point stays twin-safe throughout.

So the continuation point is now much cleaner: if the true moving-throat grouped-`P2` branch ever leaves the universal twin-safe strip, we know the exact `c_{\rm pole}` and `\bigl(\epsilon_2,\epsilon_4\bigr)` corridor in which explicit Family-1 rescue is still possible.

---

## Stage 258 — exact one-pole selected-branch regime split

Take the generic isotropic one-pole conservative carrier
\[
\widehat Y_Q^{\rm cons}(\omega)
=
 c_{\rm geom}
+
\frac{c_{\rm pole}}{1-\omega^2/\Omega_Q^2},
\qquad
c_{\rm geom}+c_{\rm pole}=1.
\]
Then the exact selected-branch demand variables are
\[
\rho_\alpha = \frac{1}{1-c_{\rm pole}},
\qquad
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}},
\]
and, exactly as in Stage 254,
\[
\frac{\Pi_{\rm tr}}{C_{\rm mix}} = \rho_\alpha.
\]
So the support-regime split is independent of blocking in one-pole language:
\[
\Pi_{\rm tr}\le C_{\rm mix}
\iff
\rho_\alpha\le1
\iff
c_{\rm pole}\le0,
\]
\[
C_{\rm mix}<\Pi_{\rm tr}\le2C_{\rm mix}
\iff
1<\rho_\alpha\le2
\iff
0<c_{\rm pole}\le\frac12,
\]
\[
\Pi_{\rm tr}>2C_{\rm mix}
\iff
\rho_\alpha>2
\iff
c_{\rm pole}>\frac12.
\]
So mixed-only sufficiency is already impossible on every physical positive-pole branch.
The universal symmetric-lowest-twin strip is exactly
\[
0<c_{\rm pole}\le\frac12.
\]
The actual isotropic branch remains
\[
c_{\rm pole}=\frac14,
\qquad
\rho_\alpha=\frac43,
\qquad
\zeta_{\rm req}=\frac13,
\]
so it lies strictly inside the twin-safe regime.

On the exact `\Lambda_{\rm EM}`-refreshed Family-1 branch at `\epsilon_{\rm blk}=0`, the hard ceiling is still
\[
c_{\rm pole}<c_{\rm pole,max}^{(F1)}(0)
\approx 0.7116102605226109,
\]
so there remains a finite non-twin corridor
\[
\frac12<c_{\rm pole}<c_{\rm pole,max}^{(F1)}(0).
\]

---

## Stage 259 — exact geometry-lane regime surface in `\bigl(\epsilon_2,\epsilon_4\bigr)`

Using the exact Stage-75 obstruction formula
\[
c_{\rm pole}
=
\frac{1+\epsilon_4}{4(1+\epsilon_2)^2},
\]
the selected-branch demand becomes
\[
\rho_\alpha
=
\frac{4(1+\epsilon_2)^2}{4\epsilon_2^2+8\epsilon_2-\epsilon_4+3},
\qquad
\frac{\Pi_{\rm tr}}{C_{\rm mix}}=\rho_\alpha.
\]
So the exact regime surfaces are now explicit in geometry-contamination language:

### Mixed-only boundary
\[
c_{\rm pole}=0
\iff
1+\epsilon_4=0.
\]

### Universal lowest-twin boundary
\[
c_{\rm pole}=\frac12
\iff
1+\epsilon_4 = 2(1+\epsilon_2)^2.
\]
Equivalently,
\[
c_{\rm pole}-\frac12
=
\frac{(1+\epsilon_4)-2(1+\epsilon_2)^2}{4(1+\epsilon_2)^2},
\]
so the sign of
\[
(1+\epsilon_4)-2(1+\epsilon_2)^2
\]
exactly determines whether the branch is twin-safe or non-twin.

At `\epsilon_{\rm blk}=0`, the exact `\Lambda_{\rm EM}`-refreshed Family-1 upper ceiling becomes
\[
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(0)
(1+\epsilon_2)^2
\approx
2.8464410420904435
(1+\epsilon_2)^2.
\]
So the exact unblocked admissible non-twin strip is
\[
2(1+\epsilon_2)^2
<
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(0)(1+\epsilon_2)^2.
\]
The actual isotropic grouped-`P2` point
\[
\epsilon_2=\epsilon_4=0
\]
has exact twin margin
\[
2(1+0)^2-(1+0)=1,
\]
so it sits safely below the non-twin surface.

---

## Stage 260 — exact blocking monotonicity and asymmetry demand

Keeping `\epsilon_{\rm blk}` explicit, the exact support demand is
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk})
=
\frac{c_{\rm pole}}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}}.
\]
Its exact derivatives are
\[
\frac{\partial \zeta_{\rm req}}{\partial c_{\rm pole}}
=
\frac{1-\epsilon_{\rm blk}}{\bigl[1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}\bigr]^2}>0,
\]
\[
\frac{\partial \zeta_{\rm req}}{\partial \epsilon_{\rm blk}}
=
\frac{c_{\rm pole}(1-2c_{\rm pole})}{\bigl[1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}\bigr]^2}.
\]
So the branch ordering in `c_{\rm pole}` is exact, and blocking has a sign flip at the universal threshold `c_{\rm pole}=1/2`:
\[
c_{\rm pole}<\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}>0,
\]
\[
c_{\rm pole}=\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}=0,
\]
\[
c_{\rm pole}>\frac12
\Rightarrow
\partial_{\epsilon_{\rm blk}}\zeta_{\rm req}<0.
\]
So blocking hurts the twin-safe side but helps the non-twin side.

The exact excess beyond the symmetric-twin threshold is
\[
\Delta_\zeta:=\zeta_{\rm req}-1
=
\frac{(1-\epsilon_{\rm blk})(2c_{\rm pole}-1)}{1-\epsilon_{\rm blk}-(1-2\epsilon_{\rm blk})c_{\rm pole}}.
\]
Therefore:
- `\Delta_\zeta<0` in the twin-safe strip,
- `\Delta_\zeta=0` on the universal boundary `c_{\rm pole}=1/2`,
- `\Delta_\zeta>0` on the non-twin side.

The actual isotropic branch stays
\[
\zeta_{\rm req}^{(act)}=\frac{1}{3-2\epsilon_{\rm blk}},
\]
so blocking increases its demand but never drives it out of the twin-safe regime on the admissible Family-1 branch.

---

## Stage 261 — exact blocked Family-1 non-twin corridor

Let `\zeta_{\max}^{(F1)}` be the exact refreshed Family-1 support ceiling in support-ratio language. Solving
\[
\zeta_{\rm req}(c_{\rm pole};\epsilon_{\rm blk}) = \zeta_{\max}^{(F1)}
\]
for `c_{\rm pole}` gives the exact blocked hard ceiling
\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})
=
\frac{\zeta_{\max}^{(F1)}(1-\epsilon_{\rm blk})}{1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}}.
\]
Because `\zeta_{\max}^{(F1)}>1`, its derivative is strictly positive:
\[
\frac{d c_{\rm pole,max}^{(F1)}}{d\epsilon_{\rm blk}}
=
\frac{\zeta_{\max}^{(F1)}(\zeta_{\max}^{(F1)}-1)}{\bigl[1+(1-2\epsilon_{\rm blk})\zeta_{\max}^{(F1)}\bigr]^2}>0.
\]
So the hard Family-1 pole ceiling grows with blocking.

Exact endpoint values:
\[
c_{\rm pole,max}^{(F1)}(0)
=
\frac{\zeta_{\max}^{(F1)}}{1+\zeta_{\max}^{(F1)}},
\qquad
\epsilon_{\rm blk}^{crit} = \frac{1}{\zeta_{\max}^{(F1)}},
\qquad
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}^{crit})=1.
\]
So the admissible blocked non-twin corridor is exactly
\[
\frac12<c_{\rm pole}<c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}),
\]
and its width
\[
c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})-\frac12
\]
increases monotonically from
\[
0.2116102605226109\quad (\epsilon_{\rm blk}=0)
\]
to
\[
0.5\quad (\epsilon_{\rm blk}=\epsilon_{\rm blk}^{crit}).
\]
In geometry-contamination language the exact blocked non-twin strip becomes
\[
2(1+\epsilon_2)^2
<
1+\epsilon_4
<
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk})(1+\epsilon_2)^2.
\]
The upper coefficient therefore rises from
\[
4c_{\rm pole,max}^{(F1)}(0)
\approx 2.8464410420904435
\]
to
\[
4c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk}^{crit})=4.
\]

A final exact invariance is that the maximal asymmetry demand on the blocked Family-1 ceiling stays fixed:
\[
\Delta_{\zeta,\max}
=
\zeta_{\rm req}(c_{\rm pole,max}^{(F1)}(\epsilon_{\rm blk});\epsilon_{\rm blk})-1
=
\zeta_{\max}^{(F1)}-1.
\]
So blocking widens the admissible `c_{\rm pole}` / `\bigl(\epsilon_2,\epsilon_4\bigr)` corridor, but it does **not** raise the maximal support asymmetry demanded at the Family-1 hard ceiling, because that ceiling is defined at fixed `\zeta_{\max}^{(F1)}`.

---

## Net result after Stage 261

The exact one-pole / geometry-lane branch picture is now much sharper.

What is exact after this pass:

1. the selected-branch support regime split is exact in `c_{\rm pole}` language,
2. the same split is an exact phase portrait in `\bigl(\epsilon_2,\epsilon_4\bigr)`,
3. blocking hurts the twin-safe side but helps the non-twin side, with the sign flip fixed at `c_{\rm pole}=1/2`,
4. the exact `\Lambda_{\rm EM}`-refreshed Family-1 non-twin corridor widens with blocking,
5. and the maximal asymmetry demand on the Family-1 ceiling stays fixed at `\zeta_{\max}^{(F1)}-1`.

So the live theorem gap is now even narrower than it was at Stage 257.
It is no longer about whether the explicit support/source side can tolerate one-pole deformation in the abstract.
It is:

> if the actual moving-throat grouped-`P2` branch ever leaves the universal twin-safe strip, does its exact `c_{\rm pole}` / `\bigl(\epsilon_2,\epsilon_4\bigr)` placement remain inside the widened blocked Family-1 corridor where explicit support rescue is still possible?
