
# 4PN local Hamiltonian-to-ordinary translation notes

This note freezes the next exact local 4PN result after the target-import audit.

## 1. What was done

Starting from the fixed-chart **local** 4PN ADM Hamiltonian target and the exact quartic perturbative Legendre compiler, I translated the reduced COM local Hamiltonian target into the reduced COM ordinary-Lagrangian target.

To make the theorem gate exact rather than heuristic, I used the full 21-slot COM basis:

- free kinetic: 6 slots,
- \(G/r\): 5 slots,
- \(G^2/r^2\): 4 slots,
- \(G^3/r^3\): 3 slots,
- \(G^4/r^4\): 2 slots,
- \(G^5/r^5\): 1 slot.

The carried lower-order ordinary blocks \(L_1,L_2,L_3\) were the exact reduced COM 1PN, 2PN, and 3PN blocks already frozen in the hierarchy.

## 2. Exact quartic COM map

The 21-slot quartic COM map is diagonal.

If the ordinary 4PN coefficients are denoted \(L_1,\dots,L_{21}\), then each Hamiltonian slot has the form
\[
h_i = F_i(\nu) - L_i,
\]
with no cross-slot mixing. Equivalently, the Jacobian of the full 21-slot map is exactly
\[
\frac{\partial h_i}{\partial L_j} = -\delta_{ij}.
\]

This is the exact local 4PN analog of the diagonal reduced COM compiler phenomenon already seen at 3PN.

## 3. Exact ordinary local 4PN target

After importing the fixed local Hamiltonian target, the exact ordinary local target is:

### Free block
\[
L_1^{\rm target}
=
\frac{7+71\nu-1135\nu^2+4082\nu^3-4497\nu^4}{256},
\qquad
L_2^{\rm target}=L_3^{\rm target}=L_4^{\rm target}=L_5^{\rm target}=L_6^{\rm target}=0.
\]

So the free block remains a one-slot block after translation.

### \(G/r\) block
\[
L_7^{\rm target}
=
\frac{150+512\nu-7292\nu^2+9285\nu^3+10070\nu^4}{256},
\]
\[
L_8^{\rm target}
=
\frac{\nu(2414\nu^3-1723\nu^2+346\nu-16)}{64},
\]
\[
L_9^{\rm target}
=
\frac{3\nu^2(270\nu^2-191\nu+30)}{128},
\]
\[
L_{10}^{\rm target}
=
-\frac{5\nu^3(34\nu-13)}{64},
\]
\[
L_{11}^{\rm target}
=
\frac{35\nu^3(2\nu-1)}{256}.
\]

### \(G^2/r^2\) block
\[
L_{12}^{\rm target}
=
-\frac{7504\nu^4+22415\nu^3+3297\nu^2+340\nu-944}{256},
\]
\[
L_{13}^{\rm target}
=
-\frac{\nu(17648\nu^3+17503\nu^2-12188\nu+1872)}{256},
\]
\[
L_{14}^{\rm target}
=
-\frac{\nu(21888\nu^3-5105\nu^2+6915\nu-2916)}{768},
\]
\[
L_{15}^{\rm target}
=
\frac{\nu(5760\nu^3-2673\nu^2+11510\nu-2952)}{1280}.
\]

### \(G^3/r^3\) block
The exact \(G^3/r^3\) leading coefficient is
\[
L_{16}^{\rm target}
=
\frac{
30412800\nu^4+128822400\nu^3-3987675\pi^2\nu^2+258968192\nu^2-76341312\nu-1294650\pi^2\nu+23385600
}{3686400},
\]
\[
L_{17}^{\rm target}
=
\frac{\nu(6144000\nu^3+15650400\nu^2+281025\pi^2\nu+5897600\nu-7914912+166050\pi^2)}{153600},
\]
\[
L_{18}^{\rm target}
=
\frac{\nu(5836800\nu^3+4899200\nu^2-13348224\nu+579825\pi^2\nu-11250\pi^2+3289536)}{245760}.
\]

### \(G^4/r^4\) block
\[
L_{19}^{\rm target}
=
-\frac{
614400\nu^4+1382400\nu^3-3916025\pi^2\nu^2+40000704\nu^2-3160350\pi^2\nu+22640704\nu-1190400
}{1228800},
\]
\[
L_{20}^{\rm target}
=
-\frac{\nu(27648000\nu^3+89856000\nu^2+21708480\nu+7681425\pi^2\nu-5340450\pi^2+114494656)}{3686400}.
\]

### \(G^5/r^5\) block
\[
L_{21}^{\rm target}
=
-\frac{-6430720\nu^2+555225\pi^2\nu^2-16243104\nu+1403325\pi^2\nu-14400}{230400}.
\]

## 4. One-body gate

The ordinary target satisfies the strict one-body gate exactly:
\[
L_1\to \frac{7}{256},\quad
L_7\to \frac{75}{128},\quad
L_{12}\to \frac{59}{16},\quad
L_{16}\to \frac{203}{32},\quad
L_{19}\to \frac{31}{32},\quad
L_{21}\to \frac{1}{16},
\qquad (\nu\to 0),
\]
and every subleading slot vanishes in that limit.

## 5. Natural local self/static seed

The natural local self/static seed in the same reduced COM ordinary basis is
\[
L_1^{\rm seed}=\frac{7}{256}(X_A^9+X_B^9),
\qquad
L_7^{\rm seed}=\frac{75}{128}(X_A^8+X_B^8),
\qquad
L_{12}^{\rm seed}=\frac{59}{16}(X_A^7+X_B^7),
\]
\[
L_{16}^{\rm seed}=\frac{203}{32}(X_A^6+X_B^6),
\qquad
L_{19}^{\rm seed}=\frac{31}{32}(X_A^5+X_B^5),
\qquad
L_{21}^{\rm seed}=\frac{1}{16}(X_A^4+X_B^4),
\]
with \(X_A=(1+\Delta)/2\), \(X_B=(1-\Delta)/2\), and \(X_A X_B=\nu\).

The exact residual beyond that seed is nonzero in the free block and in every interaction block except the one-body limit, so the local comparable-mass 4PN content is genuinely nontrivial.

## 6. Structural comparison with the constant-coefficient local scaffold

From the earlier local scaffold audit, the constant-coefficient generic-frame local interaction scaffold has the following **COM \(\nu\)-degree ceilings**:

\[
G/r:4,\qquad
G^2/r^2:3,\qquad
G^3/r^3:3,\qquad
G^4/r^4:2,\qquad
G^5/r^5:2.
\]

The imported **Hamiltonian** target has exactly the same blockwise ceilings:
\[
G/r:4,\qquad
G^2/r^2:3,\qquad
G^3/r^3:3,\qquad
G^4/r^4:2,\qquad
G^5/r^5:2.
\]

But the translated **ordinary** target has
\[
G/r:4,\qquad
G^2/r^2:4,\qquad
G^3/r^3:4,\qquad
G^4/r^4:4,\qquad
G^5/r^5:2.
\]

So the quartic compiler itself creates extra \(\nu^4\) tails in the ordinary chart in the middle and upper local blocks:
- \(G^2/r^2\),
- \(G^3/r^3\),
- \(G^4/r^4\),

while leaving \(G^5/r^5\) at the same degree ceiling.

## 7. Main consequence

This is the key exact local 4PN lesson from the translation audit:

\[
\boxed{
\text{The fixed local 4PN generic-frame lift should be performed in the Hamiltonian chart first.}
}
\]

Reason: the Hamiltonian target matches the constant-coefficient generic-frame scaffold blockwise, while the translated ordinary target contains new \(\nu^4\) tails in the \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks that are not present in the fixed Hamiltonian target and are not reachable from the present constant-coefficient local ordinary scaffold.

So the next clean local step is **not** an ordinary-chart lift done before the generic-frame solve. It is the **Hamiltonian-chart generic-frame lift**, followed only afterward by translation back to the ordinary chart if desired.
