# 4PN generic-frame ordinary-translation notes

This note freezes the next local 4PN chart step after the Hamiltonian-chart canonical-slice audit.

## 1. What was translated

The previous audit froze a sparse exact **generic-frame Hamiltonian** representative of the
local comparable-mass residual,

after subtracting the natural Hamiltonian self/static seed. The next question was whether
that representative survives quartic translation back to the **ordinary** chart without
reintroducing the earlier ordinary-chart confusion.

The answer is yes, but with one important qualification.

## 2. Exact quartic residual compiler theorem

Once the lower-order ordinary ledger through 3PN is frozen, the quartic perturbative
Legendre compiler has the general form

\[
H_4 = F[L_1,L_2,L_3](v_0) - L_4(v_0),
\qquad v_0=M^{-1}p,
\]
where the feedback functional \(F\) depends only on the lower-order data.

Therefore, if we split the ordinary 4PN block as
\[
L_4 = L_4^{\rm seed} + \Delta L_4,
\]
and choose the seed in the **Hamiltonian-aligned** convention, then the residual obeys
\[
\boxed{
\Delta H_4 = -\Delta L_4.
}
\]
So the quartic generic-frame **residual** compiler is an exact sign flip, just as the 3PN
fixed-chart generic-frame compiler was for the 3PN residual.

## 3. Canonical ordinary residual

If
\[
\Delta H_{4,\rm loc}^{\rm can}
=
\frac{Gpq}{r}Q_{\rm can}
+
\frac{G^2pq}{r^2}T_{\rm can}
+
\frac{G^3pq}{r^3}S_{\rm can}
+
\frac{G^4}{r^4}U_{\rm can}
+
\frac{G^5}{r^5}W_{\rm can},
\]
then the translated canonical ordinary residual is simply
\[
\boxed{
\Delta L_{4,\rm loc}^{\rm can}
=
-
\Delta H_{4,\rm loc}^{\rm can}.
}
\]
So the sparse canonical Hamiltonian representative immediately gives a sparse canonical
ordinary representative of the **local comparable-mass residual**.

## 4. Exact COM verification

The translated canonical ordinary residual was checked directly against the exact fixed-chart
ordinary local 4PN target. The comparison must be made relative to the **ordinary seed aligned
with the chosen Hamiltonian self/static seed**, not relative to the earlier natural one-body
ordinary seed.

With that aligned seed, the COM reduction matches exactly in every local block:

- \(G/r\): all 5 slots match,
- \(G^2/r^2\): all 4 slots match,
- \(G^3/r^3\): all 3 slots match,
- \(G^4/r^4\): both slots match,
- \(G^5/r^5\): the static slot matches.

The ordinary residual nu-degree ceilings are unchanged by the translation:
\[
Q:(4,4,4,4,3),\quad
T:(3,3,3,3),\quad
S:(3,3,3),\quad
U:(2,2),\quad
W:(2).
\]
Equivalently: the ordinary residual itself is **not** the source of the earlier ordinary-chart
nu-tail problem.

## 5. The seed-alignment correction

The exact fixed-chart local ordinary target decomposes as
\[
\boxed{
L_{4,\rm target}^{\rm loc}
=
L_{4,\rm seed}^{\rm natural}
+
\delta L_{4,\rm seed}
+
\Delta L_{4,\rm loc}^{\rm can},
}
\]
where
\[
\delta L_{4,\rm seed}

equiv
L_{4,\rm seed}^{\rm aligned}-L_{4,\rm seed}^{\rm natural}.
\]
This \(\delta L_{4,\rm seed}\) is exactly the chart/seed object that was previously appearing as
an ordinary-chart obstruction.

Its decisive structural property is that it injects new \(\nu^4\) content into the middle/upper
ordinary blocks:
\[
G^2/r^2:
\text{degree }4,
\qquad
G^3/r^3:
\text{degree }4,
\qquad
G^4/r^4:
\text{degree }4.
\]
That is precisely the pattern seen earlier in the COM ordinary translation audit.

So the earlier ordinary-chart issue is not a failure of the Hamiltonian generic-frame lift.
It is the separate problem of realizing the aligned ordinary seed (or equivalently
\(\delta L_{4,\rm seed}\)) in the generic frame.

## 6. What is now settled

This step settles the local comparable-mass chart question as far as the Hamiltonian-derived
residual is concerned:

1. the sparse canonical Hamiltonian representative does translate cleanly back to the ordinary chart;
2. the translated ordinary residual is exact and COM-correct;
3. the residual itself carries no new ordinary-side span obstruction;
4. the only remaining local chart problem is the explicit generic-frame realization of the
   aligned ordinary seed / seed-alignment correction.

## 7. What remains next

So the local-first 4PN program is now split even more sharply:

\[
\boxed{
\text{(A) Hamiltonian-derived generic-frame residual: solved},
}
\]
\[
\boxed{
\text{(B) ordinary seed-alignment correction: next local chart problem},
}
\]
\[
\boxed{
\text{(C) nonlocal tail / hereditary bridge: still separate}.
}
\]

That is the clean carry-forward status after this translation step.
