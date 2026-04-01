# 4PN generic-frame aligned-seed lift notes

This note freezes the next local 4PN chart step after the ordinary residual translation audit.

## 1. What the actual local problem was

The previous audit had already settled the comparable-mass local residual:

- the Hamiltonian-chart generic-frame lift exists,
- the translated ordinary residual exists,
- and that ordinary residual matches the exact fixed-chart ordinary target **once the target is measured relative to the Hamiltonian-aligned ordinary seed**.

So the remaining local issue was never the comparable-mass residual itself. It was the
seed-alignment correction

after subtracting the natural one-body ordinary seed,

a.k.a.

delta_seed = L4,seed^(aligned) - L4,seed^(natural).

At COM level, the local interaction part of this correction is nonzero in the

- Q block (
  G/r
  ),
- T block (
  G^2/r^2
  ),
- S block (
  G^3/r^3
  ),
- U block (
  G^4/r^4
  ),

and vanishes in the top static W block (
G^5/r^5
).

The decisive degree pattern is:

- Q correction: degree 4 in \(\nu\),
- T correction: degree 4,
- S correction: degree 4,
- U correction: degree 4,
- W correction: zero.

The pure-kinetic COM correction is separate and was not part of the local interaction lift.

## 2. Why the old ordinary comparable-mass scaffold failed

The old local ordinary comparable-mass interaction scaffold was the same constant-coefficient
family already used in the earlier local span audit:

- \(Q\): mass degree 0,
- \(T\): mass degree 1,
- \(S\): mass degree 2,
- \(U\): mass degree 3,
- \(W\): mass degree 4 static.

That scaffold lifts the comparable-mass local residual just fine in the Hamiltonian chart,
but it does **not** lift the ordinary seed-alignment correction.

The exact blockwise result is:

- \(Q_\delta\) is already liftable in the old \(Q\) family,
- \(W_\delta=0\),
- but \(T_\delta\), \(S_\delta\), and \(U_\delta\) are **not** in the image of the old constant-coefficient ordinary interaction scaffold.

So the earlier ordinary-chart obstruction was real, but it was also narrowly located: it sits in
the seed sector, not in the comparable-mass residual sector.

## 3. Minimal structured seed-sector enlargement

The exact lift works once the **seed sector** is enlarged in a structured way.

Write the original local generic-frame block families as

- \(\mathcal Q\),
- \(\mathcal T\),
- \(\mathcal S\),
- \(\mathcal U\),
- \(\mathcal W\),

using the same invariant variables \((a,b,c,d,e,p,q)\) as in the earlier local scaffold audits.

Then the seed-alignment correction lies in

\[
Q_\delta \in \mathcal Q,
\]
\[
T_\delta \in \mathcal T \oplus (pq)\,\mathcal T,
\]
\[
S_\delta \in \mathcal S \oplus (pq)\,\mathcal S,
\]
\[
U_\delta \in \mathcal U \oplus (p^2q^2)\,\mathcal U,
\]
\[
W_\delta = 0.
\]

This is the key structural result of the step.

The blockwise coefficient-space ranks in those enlarged seed spaces are full:

- \(Q\): rank \(20/20\),
- \(T\): rank \(16/16\),
- \(S\): rank \(12/12\),
- \(U\): rank \(8/8\).

So the ordinary seed-alignment correction is not an inexplicable extra obstruction. It is an
exactly liftable generic-frame object once the seed sector, rather than the residual sector,
is dressed in this specific way.

## 4. Exact canonical representative of the seed-alignment correction

Using the exact Gauss–Jordan lift with all null coordinates set to zero gives one sparse canonical
representative:

- \(Q_\delta\): 13 nonzero directions,
- \(T_\delta\): 14 nonzero directions,
- \(S_\delta\): 12 nonzero directions,
- \(U_\delta\): 8 nonzero directions,
- \(W_\delta=0\).

The compact block formulas are:

### \(Q_\delta\)
\[
Q_\delta = -\frac{1}{32}\Big(
214 a^3 c - 698 a^2 b^2 - 476 a^2 b c - 122 a^2 b d^2 - 290 a^2 b d e - 181 a^2 b e^2
+ 16 a^2 c^2 - 8 a^2 c d^2 - 18 a^2 d^2 e^2 - 54 a^2 d e^3 - 27 a^2 e^4
- 476 a b^2 c - 181 a b^2 d^2 - 290 a b^2 d e - 122 a b^2 e^2
+ 30 a d^3 e^3 + 15 a d^2 e^4
+ 214 b^3 c + 16 b^2 c^2 - 8 b^2 c e^2 - 27 b^2 d^4 - 54 b^2 d^3 e - 18 b^2 d^2 e^2
+ 15 b d^4 e^2 + 30 b d^3 e^3
\Big).
\]

### \(T_\delta\)
\[
T_\delta = -\frac{1}{96}\Big(
2814 a^2 b p^2 q + 1065 a^2 b p - 4725 a^2 b q + 2256 a^2 c q
- 3495 a^2 d^2 p^2 q + 408 a^2 d^2 p + 8031 a^2 d e p^2 q - 1782 a^2 d e p
+ 2814 a b^2 p q^2 - 4725 a b^2 p + 1065 a b^2 q
- 80 a d^3 e q + 2736 a d^2 e^2 p^2 q - 592 a d^2 e^2 p + 80 a d^2 e^2 q
+ 2256 b^2 c p + 8031 b^2 d e p q^2 - 1782 b^2 d e q
- 3495 b^2 e^2 p q^2 + 408 b^2 e^2 q
+ 2736 b d^2 e^2 p q^2 + 80 b d^2 e^2 p - 592 b d^2 e^2 q - 80 b d e^3 p
+ 432 d^3 e^3 p^2 q + 432 d^3 e^3 p q^2 + 576 d^3 e^3 p + 576 d^3 e^3 q
\Big).
\]

### \(S_\delta\)
\[
S_\delta = \frac{1}{192}\Big(
3672 a^2 p^3 q + 4464 a^2 p^2 q^2 - (10660+3\pi^2)a^2 p^2 + (1220-3\pi^2)a^2 p q
- 11304 a d^2 p^3 q - 7464 a d^2 p^2 q^2 - (2460-9\pi^2)a d^2 p^2 + (9\pi^2-2292)a d^2 p q
+ 4464 b^2 p^2 q^2 + 3672 b^2 p q^3 + (1220-3\pi^2)b^2 p q - (10660+3\pi^2)b^2 q^2
- 7464 b e^2 p^2 q^2 - 11304 b e^2 p q^3 + (9\pi^2-2292)b e^2 p q - (2460-9\pi^2)b e^2 q^2
+ 960 d^3 e q^2
- 11848 d^2 e^2 p^3 q - 19136 d^2 e^2 p^2 q^2 - 8512 d^2 e^2 p^2
- 11848 d^2 e^2 p q^3 - 8512 d^2 e^2 q^2
+ 960 d e^3 p^2
\Big).
\]

### \(U_\delta\)
\[
U_\delta = -\frac{pq}{96}\Big(
372 a p^3 q^2 + 108 a p^2 q^3 + (30\pi^2-1799)a p + (9\pi^2-1200)a q
+ 108 b p^3 q^2 + 372 b p^2 q^3 + (9\pi^2-1200)b p + (30\pi^2-1799)b q
+ 7740 d^2 p^3 q^2 + 2340 d^2 p^2 q^3 - (6953+96\pi^2)d^2 p - (2688+27\pi^2)d^2 q
+ 2340 e^2 p^3 q^2 + 7740 e^2 p^2 q^3 - (2688+27\pi^2)e^2 p - (6953+96\pi^2)e^2 q
\Big).
\]

and
\[
W_\delta = 0.
\]

## 5. Natural seed, aligned seed, and one full local generic-frame candidate

The natural generic-frame local ordinary seed is the direct symmetric promotion of the exact
one-body 4PN gate:

\[
Q_{\rm nat} = \frac{75}{128}(a^4+b^4),
\]
\[
T_{\rm nat} = \frac{59}{16}(q a^3 + p b^3),
\]
\[
S_{\rm nat} = \frac{203}{32}(q^2 a^2 + p^2 b^2),
\]
\[
U_{\rm nat} = \frac{31}{32}(q^3 a + p^3 b),
\]
\[
W_{\rm nat} = \frac{1}{16}(q^4+p^4).
\]

Then the aligned local seed is simply

\[
L_{4,\rm seed}^{\rm aligned,loc}
=
L_{4,\rm seed}^{\rm natural,loc}
+
\delta L_{4,\rm seed}^{\rm loc},
\]
with \(\delta L_{4,\rm seed}^{\rm loc}\) built from the blocks above.

Finally, combining that aligned seed with the already solved canonical ordinary residual
\(\Delta L_{4,\rm loc}^{\rm can}\) from the previous step gives one exact local generic-frame
ordinary candidate:

\[
\boxed{
L_{4,\rm loc}^{\rm gen}
=
L_{4,\rm seed}^{\rm natural,gen}
+
\delta L_{4,\rm seed}^{\rm gen}
+
\Delta L_{4,\rm loc}^{\rm can}.
}
\]

Its COM reduction reproduces the full fixed-chart local ordinary 4PN target block by block.

## 6. What is now settled

This step settles the local ordinary generic-frame **existence** issue.

More precisely:

1. the comparable-mass residual lift was already solved previously;
2. the seed-alignment correction is now solved as a structured generic-frame lift;
3. the natural seed plus the lifted seed-alignment correction plus the solved canonical residual
   gives one exact local ordinary generic-frame candidate;
4. therefore the local-first 4PN program no longer has a local existence bottleneck.

## 7. What remains open

What remains is sharper than before:

- the **uniqueness / interpretation** of this aligned-seed generic-frame sector,
- and, still completely separate, the **nonlocal tail / hereditary quadrupole bridge**.

So the 4PN local-first program now has both:

- an exact Hamiltonian-chart generic-frame representative of the local comparable-mass residual,
- and an exact generic-frame ordinary realization of the aligned local seed correction.

That is the carry-forward status after this step.
