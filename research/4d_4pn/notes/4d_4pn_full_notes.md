# 4PN kickoff notes

## Step 1 — exact one-body 4PN gate

Using the exact isotropic Schwarzschild worldline Lagrangian,

\[
L_{\rm exact}/m
=
-c^2\sqrt{\left(\frac{1-u/2}{1+u/2}\right)^2-(1+u/2)^4\frac{v^2}{c^2}},
\qquad u\equiv U/c^2,
\]

its 4PN \((c^{-8})\) coefficient is

\[
\frac{L_{4\mathrm{PN},1\text{-body}}}{m}
=
\frac{7}{256}\frac{v^{10}}{c^8}
+\frac{75}{128}\frac{Uv^8}{c^8}
+\frac{59}{16}\frac{U^2v^6}{c^8}
+\frac{203}{32}\frac{U^3v^4}{c^8}
+\frac{31}{32}\frac{U^4v^2}{c^8}
+\frac{1}{16}\frac{U^5}{c^8}.
\]

## Step 1a — direct continuation of the 3PN one-body package

Continue the 3PN denominator-style self sector by

\[
D(u)=1-4u+8u^2+d_3u^3+d_4u^4,
\]

with carried lower-order values

\[
d_3=-\frac{45}{4},
\qquad
\mu_{\rho,3}=\frac14,
\qquad
s_{24}=-\frac1{16}.
\]

Using the same direct continuation as at 3PN,

\[
L_{\rm cand}
=
-c^2(1-u)\sqrt{1-\frac{v^2/c^2}{D(u)}}
-
\frac{U^2}{2c^2}
+
\frac{U^3}{4c^4}
-
\frac{\mu_{\rho,3}U^4}{2c^6}
+
\frac{\mu_{\rho,4}U^5}{2c^8}
+
s_{24}\frac{U^2v^4}{c^6}
+
s_{34}\frac{U^3v^4}{c^8}
+
s_{26}\frac{U^2v^6}{c^8},
\]

one finds:

- a quartic denominator plus a quartic static gate is **not enough**;
- two genuinely new self-sector slots survive at 4PN: \(U^3v^4/c^8\) and \(U^2v^6/c^8\).

The minimal direct-continuation repair ledger is

\[
\boxed{
\mu_{\rho,4}=\frac18,
\qquad
d_4=\frac{205}{16},
\qquad
s_{34}=-\frac{15}{32},
\qquad
s_{26}=-\frac1{16}.
}
\]

So the 4PN one-body gate already shows that 4PN is not a one-parameter extension of 3PN.

## Step 2 — exact quartic Legendre compiler

Let

\[
L=L_0+\varepsilon L_1+\varepsilon^2L_2+\varepsilon^3L_3+\varepsilon^4L_4,
\]

with quadratic \(L_0\) and constant Newtonian mass matrix \(M\). Define

\[
v_0=M^{-1}p,
\qquad
A_0=(\partial_vL_1)|_{v_0},
\qquad
B_0=(\partial_vL_2)|_{v_0},
\qquad
D_0=(\partial_vL_3)|_{v_0},
\]

\[
C_0=(\partial_v^2L_1)|_{v_0},
\qquad
E_0=(\partial_v^2L_2)|_{v_0},
\qquad
T_0=(\partial_v^3L_1)|_{v_0}.
\]

Then the exact quartic perturbative Legendre coefficients are

\[
H_1=-L_1(v_0),
\]

\[
H_2=-L_2(v_0)+\frac12A_0^TM^{-1}A_0,
\]

\[
H_3=-L_3(v_0)+A_0^TM^{-1}B_0-\frac12A_0^TM^{-1}C_0M^{-1}A_0,
\]

\[
\boxed{
H_4=
-L_4(v_0)
+A_0^TM^{-1}D_0
+\frac12B_0^TM^{-1}B_0
-B_0^TM^{-1}C_0M^{-1}A_0
-\frac12A_0^TM^{-1}E_0M^{-1}A_0
+\frac12A_0^TM^{-1}C_0M^{-1}C_0M^{-1}A_0
+\frac16T_0[M^{-1}A_0,M^{-1}A_0,M^{-1}A_0].
}
\]

This is the compiler that has to be frozen before any trustworthy local 4PN comparable-mass solve.

## Immediate consequence for the 4PN roadmap

The correct next move is now clear:

1. keep the full 4PN theorem split as
   \[
   L_{4\mathrm{PN}}=L_{4\mathrm{PN}}^{\rm local}+L_{4\mathrm{PN}}^{\rm tail},
   \]
2. use the one-body gate above to freeze the admissible local self-sector,
3. use the quartic compiler above as the unique ordinary/Hamiltonian bridge,
4. only then build the local generic-frame 4PN scaffold and import a fixed-chart local 4PN target,
5. and keep the tail sector symbolic until the quadrupole normalization issue is closed.
# 4PN local scaffold notes

This note records the first local comparable-mass 4PN structural audit after the one-body gate and quartic compiler were frozen.

## 1. Natural local 4PN self/static seed

Promoting the exact one-body 4PN Schwarzschild target symmetrically gives the natural local seed

\[
L_{4,\mathrm{seed}}=
\frac{7}{256}(m_A v_A^{10}+m_B v_B^{10})
+\frac{75Gm_A m_B}{128r}(v_A^8+v_B^8)
+\frac{59G^2m_A m_B}{16r^2}(m_B v_A^6+m_A v_B^6)
\]
\[
\qquad
+\frac{203G^3m_A m_B}{32r^3}(m_B^2 v_A^4+m_A^2 v_B^4)
+\frac{31G^4m_A m_B}{32r^4}(m_B^3 v_A^2+m_A^3 v_B^2)
+\frac{G^5m_A m_B}{16r^5}(m_B^4+m_A^4),
\]
understanding the whole block as the coefficient of \(c^{-8}\).

This seed isolates the exact one-body/self/static content before any genuine comparable-mass local 4PN residual is solved.

## 2. Raw local 4PN generic-frame interaction scaffold

Using the same exchange-symmetric invariant language as the 3PN generic-frame scaffold,

- \(a=v_A^2\),
- \(b=v_B^2\),
- \(c=v_A\!\cdot\!v_B\),
- \(d=v_A\!\cdot\!n\),
- \(e=v_B\!\cdot\!n\),
- \(p=m_A\),
- \(q=m_B\),

and imposing strict test-mass vanishing on the branch \(m_A\to0\), \(v_B\to0\), the raw local 4PN residual basis sizes are

- \(G/r\) degree-8 block: **52** directions,
- \(G^2/r^2\) degree-6 block: **46** directions,
- \(G^3/r^3\) degree-4 block: **29** directions,
- \(G^4/r^4\) degree-2 block: **10** directions,
- \(G^5/r^5\) static block: **2** directions.

So the raw local 4PN interaction scaffold has **139** exchange-symmetric directions before any contact/gauge reduction.

## 3. COM projection and interaction-slot counts

Reducing to the COM frame gives the local 4PN interaction-slot structure

- \(G/r\): 5 COM slots,
- \(G^2/r^2\): 4 COM slots,
- \(G^3/r^3\): 3 COM slots,
- \(G^4/r^4\): 2 COM slots,
- \(G^5/r^5\): 1 COM slot,

for a total of **15 COM interaction slots**.

These are the interaction analog of the full 21-slot local 4PN COM ordinary-Lagrangian basis, whose remaining 6 slots are pure kinetic.

## 4. Blockwise COM image ranks

The raw generic-frame scaffold is already blockwise COM-complete:

- \(52\to5\) with nullity **47**,
- \(46\to4\) with nullity **42**,
- \(29\to3\) with nullity **26**,
- \(10\to2\) with nullity **8**,
- \(2\to1\) with nullity **1**.

So the total COM interaction rank is **15**, with total COM-nullity **124**.

This is the main structural result of the first local 4PN scaffold audit:

> There is no early kinematic shortage at local 4PN. The raw constant-coefficient exchange-symmetric generic-frame interaction scaffold is already large enough to span the full COM interaction image blockwise.

The real difficulty is therefore not blockwise undercompleteness but the very large COM-null sector, which must later be quotiented by the true fixed-chart Hamiltonian conditions rather than by COM data alone.

## 5. Minimal pivot-spanning subset

A convenient minimal spanning subset for the blockwise COM image is:

### \(G/r\) block
- \(a^2 b^2\)
- \(a^2 b d^2 + a b^2 e^2\)
- \(a^2 d^2 e^2 + b^2 d^2 e^2\)
- \(a d^2 e^4 + b d^4 e^2\)
- \(d^4 e^4\)

### \(G^2/r^2\) block
- \(a^2 b p + a b^2 q\)
- \(a^2 d^2 p + b^2 e^2 q\)
- \(a d^2 e^2 p + b d^2 e^2 q\)
- \(d^3 e^3 p + d^3 e^3 q\)

### \(G^3/r^3\) block
- \(a^2 p^2 + b^2 q^2\)
- \(a d^2 p^2 + b e^2 q^2\)
- \(d^2 e^2 p^2 + d^2 e^2 q^2\)

### \(G^4/r^4\) block
- \(a p^2 q + b p q^2\)
- \(d^2 p^2 q + e^2 p q^2\)

### \(G^5/r^5\) block
- \(p^2 q^2\)

These pivots already span the full 15-slot COM interaction image and are a clean starting point for the eventual fixed-chart local 4PN target import.

## 6. Immediate consequence for the roadmap

The next mathematically clean local step is now sharply defined:

1. import a fixed-chart local 4PN target,
2. apply the exact quartic Hamiltonian compiler,
3. solve the 15-slot local interaction residual inside a generic-frame basis while quotiening the 124-dimensional COM-null sector by the fixed-chart Hamiltonian conditions.

The tail/hereditary sector remains separate and should still be treated as the quadrupole-bridge problem rather than folded into the local solve.
# 4PN local target-import notes

## What is now frozen

The fixed-chart **local** 4PN target is now imported in the reduced COM Hamiltonian chart.
The source choice is:

- **Primary local target:** Jaranowski–Schafer, *Derivation of local-in-time fourth post-Newtonian ADM Hamiltonian for spinless compact binaries*, arXiv:1508.01016, Eq. (8.41e).
- **Local/nonlocal split:** Damour–Jaranowski–Schafer, *Nonlocal-in-time action for the fourth post-Newtonian conservative dynamics of two-body systems*, arXiv:1401.4548.
- **Cross-check families:** Marchand–Bernard–Blanchet–Faye (harmonic/Fokker 4PN) and Foffa–Sturani (EFT 4PN regularized local Lagrangian).

This fixes the correct source logic for our local-first program:


action/Hamiltonian at 4PN = local instantaneous sector + separate nonlocal tail sector.

So the present import audit deliberately freezes only the local instantaneous Hamiltonian target.

## Independent one-body check

The strict one-body 4PN Hamiltonian gate can be derived directly from our already frozen one-body Schwarzschild Lagrangian by the exact quartic Legendre compiler. The resulting reduced one-body 4PN Hamiltonian coefficients are

\[
\frac{7}{256}p^{10},
\qquad
\frac{45}{128}\frac{p^8}{r},
\qquad
\frac{13}{8}\frac{p^6}{r^2},
\qquad
\frac{105}{32}\frac{p^4}{r^3},
\qquad
\frac{105}{32}\frac{p^2}{r^4},
\qquad
-\frac{1}{16}\frac{1}{r^5}.
\]

The imported ADM local target matches this gate exactly in the strict test-mass limit \(\nu\to0\).

## Local 4PN COM slot structure

The imported local COM Hamiltonian naturally organizes into

- 1 pure free-kinetic slot,
- 5 slots in the \(G/r\) block,
- 4 slots in the \(G^2/r^2\) block,
- 3 slots in the \(G^3/r^3\) block,
- 2 slots in the \(G^4/r^4\) block,
- 1 slot in the \(G^5/r^5\) block.

So the local COM Hamiltonian has **16 slots total**, of which **15 are interaction slots**.

This matches the structural count already seen in the local generic-frame scaffold audit.

## Exact theorem gate answered

The first exact local target-import question was:

> Do the upper local blocks \(G^4/r^4\) or \(G^5/r^5\) contain \(\nu^3\) or \(\nu^4\) tails?

Answer:

- the free kinetic slot goes up to \(\nu^4\),
- the \(G/r\) block goes up to \(\nu^4\),
- the \(G^2/r^2\) block goes up to \(\nu^3\),
- the \(G^3/r^3\) block goes up to \(\nu^3\),
- the **\(G^4/r^4\) block stops at \(\nu^2\)**,
- the **\(G^5/r^5\) block stops at \(\nu^2\)**.

So the imported fixed-chart local target does **not** force any \(\nu^3\) or \(\nu^4\) tails in the upper local interaction blocks.

## Immediate consequence for the roadmap

This removes one possible obstruction.

The next exact local step is no longer source selection. It is:

1. translate the frozen local Hamiltonian target back through the quartic compiler,
2. compare the resulting ordinary-Lagrangian target against the local generic-frame scaffold,
3. and determine whether any seed/contact/gauge refinement is still needed before the actual local 4PN comparable-mass solve.

The tail/hereditary bridge remains separate and should still be treated as the quadrupole-bridge problem, not folded into the local solve.

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
# 4PN Hamiltonian-chart generic-frame lift notes

This note freezes the next exact local 4PN step after the Hamiltonian-to-ordinary translation audit.

## 1. What was actually tested

The earlier local-scaffold audit showed that the raw constant-coefficient exchange-symmetric generic-frame local 4PN interaction scaffold is already **slot-complete** after COM projection:

- \(G/r\): 5 COM interaction slots,
- \(G^2/r^2\): 4 COM interaction slots,
- \(G^3/r^3\): 3 COM interaction slots,
- \(G^4/r^4\): 2 COM interaction slots,
- \(G^5/r^5\): 1 COM interaction slot.

But the later Hamiltonian-to-ordinary audit showed that the **ordinary** target develops extra \(\nu^4\) tails in the \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks after quartic translation, whereas the fixed local **Hamiltonian** target matches the scaffold degree ceilings exactly.

So the next exact question was:

> Does the constant-coefficient generic-frame scaffold span the full fixed-chart **Hamiltonian** local 4PN comparable-mass residual, or is there still a seed/obstruction problem hiding there?

## 2. Natural Hamiltonian seed

Using the exact one-body 4PN Hamiltonian gate
\[
\frac{7}{256}p^{10},
\qquad
\frac{45}{128}\frac{p^8}{r},
\qquad
\frac{13}{8}\frac{p^6}{r^2},
\qquad
\frac{105}{32}\frac{p^4}{r^3},
\qquad
\frac{105}{32}\frac{p^2}{r^4},
\qquad
-\frac{1}{16}\frac{1}{r^5},
\]
the natural reduced COM Hamiltonian self/static seed is

\[
K^{\rm seed}=\frac{7}{256}(X_A^9+X_B^9),
\]
\[
Q_1^{\rm seed}=\frac{45}{128}(X_A^8+X_B^8),
\qquad
T_1^{\rm seed}=\frac{13}{8}(X_A^7+X_B^7),
\qquad
S_1^{\rm seed}=\frac{105}{32}(X_A^6+X_B^6),
\]
\[
U_1^{\rm seed}=\frac{105}{32}(X_A^5+X_B^5),
\qquad
W_1^{\rm seed}=-\frac{1}{16}(X_A^4+X_B^4),
\]
with \(X_A=(1+\Delta)/2\), \(X_B=(1-\Delta)/2\), \(X_A X_B=\nu\), and all subleading slots initially zero.

A first exact outcome is that this seed already captures the **entire free 4PN COM Hamiltonian slot**:
\[
\Delta K = 0.
\]

So the generic-frame local lift only has to solve the **interaction** residual.

## 3. Polynomial coefficient-space formulation

To test the lift exactly, each COM block was decomposed into explicit \(\nu\)-coefficient space:

- \(G/r\): \(5\) slots, each with \(\nu,\nu^2,\nu^3,\nu^4\) coefficients \(\Rightarrow 20\) coefficients,
- \(G^2/r^2\): \(4\) slots, each with \(\nu,\nu^2,\nu^3\) coefficients \(\Rightarrow 12\),
- \(G^3/r^3\): \(3\) slots, each with \(\nu,\nu^2,\nu^3\) coefficients \(\Rightarrow 9\),
- \(G^4/r^4\): \(2\) slots, each with \(\nu,\nu^2\) coefficients \(\Rightarrow 4\),
- \(G^5/r^5\): \(1\) slot, with \(\nu,\nu^2\) coefficients \(\Rightarrow 2\).

This is stricter than just counting COM interaction slots. It asks whether the scaffold spans the full **polynomial data** of the imported local Hamiltonian target.

## 4. Exact blockwise ranks

The constant-coefficient generic-frame Hamiltonian scaffold achieves full rank in every block:

\[
G/r:\quad 52 \to 20 \quad\text{with nullity }32,
\]
\[
G^2/r^2:\quad 46 \to 12 \quad\text{with nullity }34,
\]
\[
G^3/r^3:\quad 29 \to 9 \quad\text{with nullity }20,
\]
\[
G^4/r^4:\quad 10 \to 4 \quad\text{with nullity }6,
\]
\[
G^5/r^5:\quad 2 \to 2 \quad\text{with nullity }0.
\]

So the total interaction coefficient-space dimensions are

- rows: \(20+12+9+4+2=47\),
- columns: \(52+46+29+10+2=139\),
- total rank: \(47\),
- total nullity: \(92\).

The important point is not the large nullity — that is expected in a generic-frame scaffold — but the fact that the **rank is maximal in every block**.

## 5. Main theorem outcome

There is **no Hamiltonian-chart obstruction** analogous to the ordinary-chart \(\nu^4\)-tail problem.

More concretely:

1. the free 4PN COM Hamiltonian slot is already absorbed by the natural Hamiltonian seed;
2. every local interaction residual slot of the imported fixed-chart COM Hamiltonian target lies in the image of the constant-coefficient exchange-symmetric generic-frame scaffold;
3. no extra seed refinement is needed before the Hamiltonian-chart generic-frame local lift;
4. one exact generic-frame representative exists blockwise and reduces back to the imported COM Hamiltonian residual exactly.

So the correct local-first program is now sharply confirmed:

\[
\boxed{
\text{fixed local 4PN Hamiltonian target}
\;\longrightarrow\;
\text{generic-frame Hamiltonian lift}
\;\longrightarrow\;
\text{ordinary translation only afterward}.
}
\]

## 6. Contrast with the ordinary-chart translation audit

This is the clean structural comparison:

- In the **Hamiltonian** chart, the imported local target has blockwise degree ceilings
  \[
  G/r:4,\qquad G^2/r^2:3,\qquad G^3/r^3:3,\qquad G^4/r^4:2,\qquad G^5/r^5:2,
  \]
  exactly matching the constant-coefficient scaffold.

- In the translated **ordinary** chart, the target becomes
  \[
  G/r:4,\qquad G^2/r^2:4,\qquad G^3/r^3:4,\qquad G^4/r^4:4,\qquad G^5/r^5:2,
  \]
  which is why the ordinary-chart lift looked obstructed.

So the previous “Hamiltonian first” recommendation is no longer just a heuristic. It is now backed by an explicit exact rank test.

## 7. What remains open after this step

The local Hamiltonian-chart lift is now **structurally open but not obstructed**. The remaining tasks are:

1. fix a canonical generic-frame representative inside the \(92\)-direction null family,
2. import or reconstruct the full fixed-chart **generic-frame** local 4PN Hamiltonian target if we want uniqueness in a chosen chart rather than mere existence,
3. then translate that representative back through the quartic compiler to the ordinary chart,
4. while keeping the tail/hereditary quadrupole bridge separate.

So the local 4PN program is now in a better position than the previous ordinary-chart result suggested: the real local obstruction is no longer span, but only chart-fixing / canonical-choice bookkeeping.

# 4PN Hamiltonian-chart canonical-slice notes

This note freezes the next exact local 4PN step after the Hamiltonian-chart generic-frame lift audit.

## 1. What was actually fixed

The previous lift audit proved **existence**: every local interaction block of the fixed-chart COM 4PN Hamiltonian residual lies in the image of the constant-coefficient exchange-symmetric generic-frame scaffold, with total nullity \(92\).

The unresolved local question was therefore no longer span. It was:

> What is a clean canonical generic-frame Hamiltonian representative inside that 92-direction COM-blind family?

The present audit answers that.

## 2. The COM-null ideals at 4PN

The generic-frame invariants satisfy the same COM identities as at 3PN:

\[
\mathcal C_1 = p a + q c,\qquad
\mathcal C_2 = p c + q b,\qquad
\mathcal C_3 = p d + q e,
\]
\[
\mathcal C_4 = ab-c^2,\qquad
\mathcal C_5 = ae-cd,\qquad
\mathcal C_6 = bd-ce.
\]

The 4PN local Hamiltonian-chart nullities are

\[
Q:32,\qquad T:34,\qquad S:20,\qquad U:6,\qquad W:0.
\]

The exact ideal-membership result is:

- all \(32\) \(Q\)-block null vectors lie in the **full** COM ideal
  \[
  \langle \mathcal C_1,\mathcal C_2,\mathcal C_3,\mathcal C_4,\mathcal C_5,\mathcal C_6\rangle;
  \]
- all \(34\) \(T\)-block null vectors lie already in the **linear-momentum** ideal
  \[
  \langle \mathcal C_1,\mathcal C_2,\mathcal C_3\rangle;
  \]
- all \(20\) \(S\)-block null vectors lie already in the same linear-momentum ideal;
- all \(6\) \(U\)-block null vectors also lie in that linear-momentum ideal;
- the \(W\)-block has **zero nullity**.

So the whole 92-direction ambiguity is purely COM-blind algebraic freedom, and the top static \(G^5/r^5\) block is already fully fixed once the COM target is fixed.

## 3. Canonical quotient slice

The canonical slice used in the audit is the exact Gauss–Jordan lift with **all null coordinates set to zero**.

This produces a sparse explicit generic-frame Hamiltonian representative with

\[
Q:11,\qquad T:12,\qquad S:9,\qquad U:4,\qquad W:2
\]
nonzero scaffold directions.

## 4. Explicit canonical representative

Write the local comparable-mass Hamiltonian residual as

\[
\Delta H_{4,\mathrm{loc}}^{\mathrm{can}}
=
\frac{Gpq}{r}Q_{\mathrm{can}}
+
\frac{G^2pq}{r^2}T_{\mathrm{can}}
+
\frac{G^3pq}{r^3}S_{\mathrm{can}}
+
\frac{G^4}{r^4}U_{\mathrm{can}}
+
\frac{G^5}{r^5}W_{\mathrm{can}}.
\]

Then one exact canonical representative is:

### \(G/r\) block
\[
\begin{aligned}
Q_{\mathrm{can}}={}&
-\frac{11}{64}a^2b^2
+\frac{5}{256}(a^2bc+ab^2c)
-\frac{3}{32}(a^2bd^2+ab^2e^2)
+\frac{1}{64}(a^2bde+ab^2de)
\\
&-\frac{27}{64}(a^2c^2+b^2c^2)
-\frac{9}{64}(a^2d^2e^2+b^2d^2e^2)
+\frac{3}{128}(a^2de^3+b^2d^3e)
+\frac{3}{64}(a^2e^4+b^2d^4)
\\
&-\frac{5}{32}(ad^2e^4+bd^4e^2)
+\frac{5}{64}(ad^3e^3+bd^3e^3)
-\frac{35}{256}(d^5e^3+d^3e^5).
\end{aligned}
\]

### \(G^2/r^2\) block
\[
\begin{aligned}
T_{\mathrm{can}}={}&
-\frac{87}{128}(a^2bp+ab^2q)
-\frac{2227}{256}(a^2bq+ab^2p)
+\frac{63}{64}(a^2cq+b^2cp)
\\
&+\frac{49}{16}(a^2d^2p+b^2e^2q)
-\frac{435}{64}(a^2dep+b^2deq)
+\frac{2435}{256}(a^2e^2p+b^2d^2q)
\\
&-\frac{183}{16}(ad^2e^2p+bd^2e^2q)
-\frac{8305}{768}(ad^2e^2q+bd^2e^2p)
+\frac{889}{192}(ad^3eq+bde^3p)
\\
&-\frac{5343}{1280}(d^3e^3p+d^3e^3q)
+\frac{325}{128}(d^4e^2q+d^2e^4p)
-\frac{369}{160}(d^5eq+de^5p).
\end{aligned}
\]

### \(G^3/r^3\) block
\[
\begin{aligned}
S_{\mathrm{can}}={}&
\left(-\frac{3307423}{28800}+\frac{40483\pi^2}{16384}\right)(a^2p^2+b^2q^2)
+\left(-\frac{211189}{19200}+\frac{2749\pi^2}{8192}\right)(a^2pq+b^2pq)
\\
&+\left(-\frac{184897}{900}+\frac{34985\pi^2}{8192}\right)abpq
+\left(\frac{139241}{1200}-\frac{12507\pi^2}{2048}\right)(ad^2p^2+be^2q^2)
\\
&+\left(\frac{63347}{1600}-\frac{1059\pi^2}{1024}\right)(ad^2pq+be^2pq)
+\left(-\frac{716971}{9600}+\frac{10389\pi^2}{2048}\right)(adep^2+bdeq^2)
\\
&+\left(-\frac{16727}{384}-\frac{35655\pi^2}{16384}\right)(d^2e^2p^2+d^2e^2q^2)
+\left(-\frac{51193}{960}-\frac{36405\pi^2}{8192}\right)d^2e^2pq
\\
&+\left(\frac{23533}{1280}-\frac{375\pi^2}{8192}\right)(d^3eq^2+de^3p^2).
\end{aligned}
\]

### \(G^4/r^4\) block
\[
\begin{aligned}
U_{\mathrm{can}}={}&
\left(\frac{930047}{9600}-\frac{551243\pi^2}{49152}\right)(ap^2q+bpq^2)
+\left(\frac{500761}{19200}-\frac{21837\pi^2}{8192}\right)(apq^2+bp^2q)
\\
&+\left(\frac{274387}{1600}-\frac{62047\pi^2}{49152}\right)(d^2p^2q+e^2pq^2)
+\left(\frac{3401779}{57600}-\frac{28691\pi^2}{24576}\right)(d^2pq^2+e^2p^2q).
\end{aligned}
\]

### \(G^5/r^5\) block
\[
W_{\mathrm{can}}=
\left(-\frac{169799}{2400}+\frac{6237\pi^2}{1024}\right)(p^3q+pq^3)
+\left(-\frac{609427}{3600}+\frac{44825\pi^2}{3072}\right)p^2q^2.
\]

## 5. Exact verification

The audit verifies directly that this canonical representative reduces back to the imported fixed-chart local COM Hamiltonian residual in every slot:

- \(Q\): all 5 slots match exactly,
- \(T\): all 4 slots match exactly,
- \(S\): all 3 slots match exactly,
- \(U\): both slots match exactly,
- \(W\): the static slot matches exactly.

So this is not just one heuristic sparse guess. It is an exact canonical quotient-slice representative.

## 6. What this does and does not settle

This step **does** settle the local algebraic existence problem at the generic-frame Hamiltonian level:

- the 92-direction ambiguity is purely COM-blind;
- the relevant null ideals are now identified;
- and one exact sparse canonical representative is frozen.

This step does **not** yet prove fixed-chart uniqueness against the true full generic-frame ADM Hamiltonian target. For that we would still need either:

1. the fixed-chart generic-frame local 4PN Hamiltonian target itself, or
2. an exact generic-frame Hamiltonian compiler uniqueness theorem strong enough to eliminate every remaining chart ambiguity directly.

So the local-first program is now:

\[
\boxed{
\text{canonical Hamiltonian generic-frame representative}
\;\longrightarrow\;
\text{generic-frame quartic translation to the ordinary chart}
\;\longrightarrow\;
\text{tail bridge kept separate}.
}
\]

That is the next clean move.
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
# Local 4PN master/referee assembly — result ledger

## What this step accomplished

This step closes the **local** 4PN chain end to end.

The new referee audit rebuilds the full fixed-chart local 4PN ordinary target from the imported local ADM Hamiltonian target using the exact quartic Hamiltonian-to-ordinary compiler, then verifies that one exact generic-frame ordinary representative of the **entire local sector** exists after combining:

1. the natural generic-frame local seed,
2. the generic-frame lift of the aligned-seed correction,
3. the already solved canonical comparable-mass local residual.

So the local-first 4PN program is now no longer blocked by local existence.

---

## Exact ingredients used

The master audit reuses the exact stage results already frozen in the current 4PN chain:

- `4pn_local_hamiltonian_to_ordinary_audit.py`
  - exact 21-slot quartic COM map,
  - exact local Hamiltonian target,
  - exact ordinary target,
  - strict one-body gate,
  - natural local self/static seed.

- `4pn_generic_frame_ordinary_translation_audit.py`
  - exact sign-flip theorem for the canonical comparable-mass residual,
  - exact aligned ordinary seed,
  - exact ordinary residual slots.

- `4pn_hamiltonian_chart_generic_frame_lift_audit.py`
  - COM reduction map,
  - block-slot extractor,
  - exact image-matrix machinery for the generic-frame lift.

---

## The three decisive checks

### 1. The free block is already fixed by the aligned seed

The aligned ordinary seed reproduces the entire free block of the local ordinary target exactly:

- slot 1 matches the full translated free coefficient,
- slots 2–6 remain zero.

So the free local 4PN ordinary sector does **not** need a separate generic-frame rescue.

### 2. The canonical local residual translates exactly by sign flip

The Hamiltonian-derived generic-frame comparable-mass local residual satisfies

\[
\Delta L_{4,\mathrm{loc}}^{\mathrm{can}} = -\Delta H_{4,\mathrm{loc}}^{\mathrm{can}}.
\]

This removes the earlier ordinary-chart worry about the true comparable-mass residual: the only remaining ordinary obstruction comes from **seed alignment**, not from the comparable-mass local lift itself.

### 3. The aligned-seed correction is exactly generic-frame liftable

The seed-alignment correction

\[
\delta L_{4,\mathrm{seed}} = L_{4,\mathrm{seed}}^{(\mathrm{aligned})} - L_{4,\mathrm{seed}}^{(\mathrm{natural})}
\]

is liftable in the minimal structured seed spaces

\[
Q,
\qquad
T \oplus (pq)T,
\qquad
S \oplus (pq)S,
\qquad
U \oplus (p^2 q^2)U,
\qquad
W=0.
\]

That is the exact generic-frame repair of the old ordinary-chart obstruction.

---

## Final theorem statement for the local sector

Define

\[
L_{4,\mathrm{loc}}^{(\mathrm{gen})}
=
L_{4,\mathrm{seed}}^{(\mathrm{natural,gen})}
+
\delta L_{4,\mathrm{seed}}^{(\mathrm{gen})}
+
\Delta L_{4,\mathrm{loc}}^{(\mathrm{can,gen})}.
\]

Then the master audit verifies that the COM reduction of this generic-frame candidate reproduces the exact fixed-chart local ordinary 4PN target slot by slot in all interaction blocks:

- `Q` block: 5 slots,
- `T` block: 4 slots,
- `S` block: 3 slots,
- `U` block: 2 slots,
- `W` block: 1 slot.

Equivalently, the **entire local 4PN ordinary interaction sector is now assembled exactly**.

---

## Canonical structured-lift sparsity

The final master audit also records one canonical sparse realization of the aligned-seed lift with the following structured basis sizes / nonzero-direction counts:

- `Q`: basis size 52, nonzero directions 13,
- `T`: basis size 92, nonzero directions 14,
- `S`: basis size 58, nonzero directions 12,
- `U`: basis size 20, nonzero directions 8,
- `W`: no lift required.

These are not yet uniqueness theorems for the aligned-seed sector, but they show that the local 4PN lift is concretely realizable rather than merely formal.

---

## What remains open after this step

This step does **not** close the full conservative 4PN theorem.

What remains open is:

1. the sharper uniqueness / interpretation question for the aligned-seed sector,
2. the **tail / hereditary bridge**, which stays separate from the local theorem by construction.

So the next mathematically clean step is to turn back to the tail side and build the quadrupole bridge on top of the now-closed local scaffold.
# 4PN tail / hereditary bridge audit — result ledger

## What this step accomplished

This step freezes the **nonlocal conservative 4PN tail bridge** after the local 4PN theorem.

The main result is that the tail side is much narrower than the local side was:

1. the exact GR conservative 4PN tail is one nonlocal STF-quadrupole functional,
2. its local logarithmic shadow occupies only the **degree-2 `G^4/r^4` block**,
3. the source representation is the same canonical STF quadrupole already isolated by the 2.5PN program,
4. and the remaining toy-model gap can be parameterized by **one scalar transport coefficient**.

So after the local 4PN theorem, the hereditary problem is **not** another large generic-frame lift. It is a sharply focused quadrupole-transport theorem.

---

## Exact GR tail functional now frozen

From the standard 4PN conservative literature, the nonlocal time-symmetric tail Hamiltonian has coefficient

\[
\alpha_{\rm tail}^{\rm GR}=\frac{G^2 M}{5c^8}.
\]

Equivalently,

\[
H_{4\rm PN}^{\rm tail,sym}(t)
=
-\frac{1}{5}\frac{G^2M}{c^8}
I_{ij}^{(3)}(t)
\,\mathrm{Pf}_{2s/c}\!\int_{-\infty}^{+\infty}\frac{dv}{|v|}
I_{ij}^{(3)}(t+v).
\]

The local logarithmic shadow is controlled by

\[
F_{\rm tail}(t)=\frac{2}{5}\frac{G^2M}{c^8}\bigl(I_{ij}^{(3)}(t)\bigr)^2.
\]

---

## Universal Newtonian quadrupole shadow

Using the Newtonian order-reduced STF quadrupole

\[
I_{ij}=\mu\,\mathrm{STF}(x_i x_j),
\]

the audit verified the exact order-reduced third derivative

\[
I_{ij}^{(3)}
=
-\frac{2GM\mu}{r^3}
\left(4x_{\langle i}v_{j\rangle}-3\frac{\dot r}{r}x_{\langle i}x_{j\rangle}\right),
\]

and therefore

\[
\bigl(I_{ij}^{(3)}\bigr)^2
=
\frac{8}{3}\frac{G^2M^2\mu^2}{r^4}
\left(12v^2-11\dot r^2\right).
\]

This gives the exact local logarithmic tail shadow

\[
\frac{F_{\rm tail}}{\mu}
=
\frac{16}{15}\,\nu\,\frac{U^4}{c^8}
\left(12v^2-11\dot r^2\right),
\qquad
U\equiv\frac{GM}{r}.
\]

So the hereditary shadow sits **only** in the local 4PN `U`-block (`G^4/r^4` with degree-2 kinematics).

That is the exact reason the tail problem is structurally much smaller than the local lift problem.

---

## Exact bridge to the 2.5PN quadrupole coefficient

The 2.5PN audit already froze the canonical Burke–Thorne coefficient target as

\[
\gamma_{\rm GR}=\frac{2G}{5c^5}.
\]

The present tail audit observes the exact arithmetic relation

\[
\alpha_{\rm tail}^{\rm GR}
=
\frac{GM}{2c^3}\,\gamma_{\rm GR}.
\]

So in GR the 4PN conservative tail coefficient is exactly **half a monopole-scattering factor** `(GM/c^3)` times the 2.5PN quadrupole coefficient.

This is the cleanest bridge between the already narrowed 2.5PN quadrupole route and the 4PN hereditary sector.

---

## Minimal toy-model tail bridge ansatz

Let

\[
\gamma_{\rm toy}
\equiv
\hat m_0^2\Gamma_5
\]

be the canonically normalized STF quadrupole coefficient already isolated by the 2.5PN notes.

Then the remaining hereditary gap can be parameterized by one scalar transport factor `Theta_tail`:

\[
\boxed{
\alpha_{\rm tail}^{\rm toy}
=
\Theta_{\rm tail}
\frac{GM}{2c_s^3}
\gamma_{\rm toy}.
}
\]

On the `c_s=c` branch,

- if `gamma_toy = gamma_GR`, and
- if `Theta_tail = 1`,

then the exact GR 4PN conservative tail coefficient is recovered automatically.

So the tail bridge has become maximally narrow:

- the **representation/source problem is already solved** by the 2.5PN quadrupole audit,
- the **local 4PN basis problem is already solved** by the local referee chain,
- and the remaining issue is only whether the moving-throat model produces the correct **monopole-scattering transport coefficient**.

---

## What remains open after this step

This step does **not** derive the full hereditary kernel from the moving-throat PDE.

What remains open is the single scalar theorem:

\[
\Theta_{\rm tail}\stackrel{?}{=}1
\]

on the same canonical STF quadrupole branch already isolated by the 2.5PN program.

Equivalently:

- no new tensor channel is needed,
- no new generic-frame existence solve is needed,
- the remaining 4PN hereditary problem is a **quadrupole tail-transport normalization theorem**.

That is the clean next target before any claim of a full conservative 4PN theorem.
# 4PN tail / hereditary bridge — result ledger

## What this step accomplished

This step freezes the **nonlocal 4PN bridge coefficient**.

The local 4PN sector was already assembled exactly in the previous referee step. What remained open was the separate hereditary/tail side. The new audit shows that the 4PN tail coefficient is **not** a new independent datum. It is fixed by the *same* quadrupole normalization that the 2.5PN program had already isolated as the last remaining universal gap.

Equivalently, once the orbital/worldtube STF quadrupole is canonically normalized, the conservative 4PN tail coefficient follows automatically.

---

## GR reference structure

The standard conservative 4PN dynamics splits into

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]

On the action side, the conservative tail term may be written in the time-symmetric form

\[
S_{\mathrm{tail}}
=
\frac{G^2 M}{5c^8}
\operatorname{Pf}_{2s/c}
\iint \frac{dt\,dt'}{|t-t'|}
I_{ij}^{(3)}(t)
I_{ij}^{(3)}(t').
\]

On the Hamiltonian side the corresponding nonlocal term is

\[
H_{\mathrm{tail}}(t)
=
-
\frac{G^2 M}{5c^8}
I_{ij}^{(3)}(t)
\operatorname{Pf}_{2s/c}
\int_{-\infty}^{+\infty}
\frac{dv}{|v|}
I_{ij}^{(3)}(t+v).
\]

So the universal GR tail coefficient is

\[
C_{\mathrm{tail}}^{\mathrm{GR}}=
\frac{G^2 M}{5c^8}.
\]

---

## Canonical bridge to the 2.5PN Burke–Thorne coefficient

The canonical Burke–Thorne local odd quadrupole coefficient is

\[
\gamma_{\mathrm{GR}}=
\frac{2G}{5c^5}.
\]

The new audit verifies the exact identity

\[
\boxed{
C_{\mathrm{tail}}^{\mathrm{GR}}
=
\frac{GM}{2c^3}\,\gamma_{\mathrm{GR}}.
}
\]

This is the cleanest structural statement of the 4PN bridge.

It says that the conservative hereditary coefficient is exactly the 2.5PN quadrupole coefficient multiplied by the standard monopole-scattering factor

\[
\frac{GM}{2c^3}.
\]

So if the 2.5PN quadrupole coefficient is correct, then the 4PN tail coefficient is fixed automatically.

---

## Harmonic-frequency structure of the tail kernel

Using a regulated finite-part proxy for the kernel

\[
2\int_0^{\infty} d\tau\,e^{-\varepsilon\tau}
\frac{\cos(\omega\tau)-1}{\tau}
=
-\ln\!\left(1+\frac{\omega^2}{\varepsilon^2}\right),
\]

the nonlocal tail action reduces on a monochromatic mode to the expected logarithmic form

\[
K_{\mathrm{tail}}(\omega)
\sim
-
\frac{G^2 M}{5c^8}
\omega^6
\ln|\omega|.
\]

Relative to the local Burke–Thorne kernel

\[
K_{\mathrm{BT}}(\omega)\sim \gamma_{\mathrm{GR}}\omega^5,
\]

the ratio is

\[
\frac{K_{\mathrm{tail}}}{K_{\mathrm{BT}}}
\sim
\frac{GM\,\omega}{c^3}\ln|\omega|,
\]

which is the expected 1.5PN tail promotion above the leading 2.5PN quadrupole reaction channel.

---

## Toy-model quadrupole branch version

The 2.5PN quadrupole program already reduced the remaining universal normalization problem to the canonical invariant quantity

\[
\gamma_{\mathrm{quad}}^{\mathrm{eff}}
=
\mathcal N_Q\frac{a^5}{27c_s^5}
=
\overline\Gamma_5
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

The new bridge result then gives the toy-model 4PN tail coefficient as

\[
\boxed{
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}
\gamma_{\mathrm{quad}}^{\mathrm{eff}}.
}
\]

Equivalently,

\[
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}\,\overline\Gamma_5
=
\frac{9GM}{2c^3}
\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

So the entire hereditary 4PN bridge is controlled by the same normalized quadrupole branch data as the 2.5PN odd channel.

---

## Exact matching target inherited from the 2.5PN normalization problem

Imposing Burke–Thorne normalization gives

\[
\gamma_{\mathrm{quad}}^{\mathrm{eff}}=
\frac{2G}{5c^5},
\]

which is equivalent to

\[
\mathcal N_Q^{\mathrm{target}}
=
\frac{54Gc_s^5}{5a^5c^5}.
\]

Substituting that into the 4PN bridge coefficient yields exactly

\[
C_{\mathrm{tail}}^{\mathrm{toy}}
=
\frac{GM}{2c^3}
\frac{2G}{5c^5}
=
\frac{G^2 M}{5c^8}
=
C_{\mathrm{tail}}^{\mathrm{GR}}.
\]

So the GR 4PN hereditary coefficient follows automatically once the 2.5PN quadrupole normalization gap is closed.

---

## Main theorem statement from this step

Let the canonically normalized orbital/worldtube STF quadrupole branch satisfy the isotropic passive/outgoing low-frequency relations already isolated in the 2.5PN audit. Then:

1. the unique compatible conservative hereditary 4PN coefficient is
   \[
   C_{\mathrm{tail}}=
   \frac{GM}{2c^3}\,\gamma_{\mathrm{quad}}^{\mathrm{eff}},
   \]
2. therefore the 4PN tail sector introduces **no new independent quadrupole normalization datum**,
3. and full GR matching of the hereditary coefficient is equivalent to the already known 2.5PN normalization target
   \[
   \gamma_{\mathrm{quad}}^{\mathrm{eff}}=
   \frac{2G}{5c^5}.
   \]

So the remaining full-4PN gap is still the same narrow passive/outgoing quadrupole normalization problem already identified by the 2.5PN program.

---

## What remains open after this step

This step does **not** close the full conservative 4PN theorem unconditionally.

What remains open is:

1. the final derivation of the passive/outgoing quadrupole normalization from the completed moving-throat PDE,
2. the explicit insertion of that normalized tail module into a full end-to-end 4PN referee/master script.

But the bottleneck is now much sharper:

- the **local** 4PN sector is already closed,
- the **tail** coefficient is now reduced to the same invariant normalization problem as 2.5PN,
- so the next full-4PN theorem step is an assembly step, not a new open-ended search.
# 4PN conditional referee master — result ledger

## What this step accomplished

This step freezes the sharpest **full conservative 4PN theorem statement** currently available inside the declared closure hierarchy.

It does not try to solve the full moving-throat PDE. Instead it replays the three decisive already-solved stages and then proves the interface theorem tying them together:

1. the **2.5PN master audit** already reduces the universal dissipative bridge to the canonically normalized STF quadrupole branch, with one narrow normalization gap left open;
2. the **local 4PN referee master** already closes the entire local instantaneous 4PN sector exactly;
3. the **4PN hereditary bridge audit** shows that the hereditary coefficient is not a new free datum, but is fixed by the same quadrupole normalization.

So the full conservative 4PN problem is now conditionally reduced to the **same** passive/outgoing quadrupole normalization problem already isolated on the 2.5PN side.

---

## Stage-audit status now frozen

The master replay confirms:

- `2_5pn_master_session_sympy_audit.py` passes,
- `4pn_local_referee_master_sympy_audit.py` passes,
- `4pn_tail_hereditary_bridge_audit.py` passes.

That means the following three legs are simultaneously stable:

- the narrowed **quadrupole normalization gap** from 2.5PN,
- the exact **local 4PN** generic-frame ordinary representative,
- the exact **4PN hereditary bridge coefficient relation**.

---

## Exact interface theorem

Let the canonically normalized STF quadrupole odd coefficient be

\[
\gamma_{\rm quad}^{\rm eff}.
\]

Then the unique compatible conservative 4PN hereditary coefficient is

\[
\boxed{
C_{\rm tail}=rac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
}
\]

So if the 2.5PN normalization target holds,

\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
\]

then automatically

\[
\boxed{
C_{\rm tail}=\frac{G^2M}{5c^8}.
}
\]

This is the exact GR conservative 4PN tail coefficient.

---

## Equivalent branch forms

The same result may be written in either of the two branch languages already isolated in the 2.5PN program.

### Geometric normalization form

\[
\gamma_{\rm quad}^{\rm eff}=\mathcal N_Q\frac{a^5}{27c_s^5},
\qquad
\mathcal N_Q^{\rm target}=\frac{54Gc_s^5}{5a^5c^5}.
\]

### Canonical invariant-pair form

\[
\gamma_{\rm quad}^{\rm eff}
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

On the Burke–Thorne target branch,

\[
\overline K_0^{\rm target}=\frac{64G}{45c^5}\Omega_Q^5,
\qquad
\overline K_2^{\rm target}=\frac{16G}{45c^5}\Omega_Q^3,
\]

and substitution gives

\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
\qquad
C_{\rm tail}=\frac{G^2M}{5c^8}.
\]

So **no new independent normalization datum opens at 4PN**.

---

## Strongest honest theorem statement available now

Within the declared closure hierarchy:

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}},
\]

with

- `L_4PN^local` already assembled exactly by the local referee chain,
- `L_4PN^tail` fixed by the same STF quadrupole normalization that controls the 2.5PN Burke–Thorne channel.

Therefore the remaining gap between the present hierarchy and a fully unconditional conservative 4PN theorem is **not** a new 4PN-specific normalization constant. It is exactly the same narrow passive/outgoing quadrupole normalization problem already isolated by the 2.5PN program.

---

## Best next move after this step

The next sharp target is no longer local 4PN algebra. That part is done.

The next real theorem gate is the one already exposed by the 2.5PN notes:

- derive the passive/outgoing quadrupole normalization from the moving-throat side,
- or derive enough of the higher conservative grouped real `P2` bundle to fix the canonical pair `(\overline K_0,\overline K_2)` on the natural branch.

Once that is closed, both

- the normalized 2.5PN odd quadrupole coefficient, and
- the normalized 4PN conservative hereditary coefficient

close simultaneously.
# Full conservative 4PN status — conditional theorem summary

## 0. What this summary is doing

This document freezes the strongest honest **full conservative 4PN theorem statement** currently available inside the declared closure hierarchy of the 4D toy-model program.

It is not claiming a solved moving-throat PDE. It is doing something narrower and more useful:

1. recording exactly what is already closed at 4PN,
2. stating the precise interface between the solved **local** sector and the **hereditary/tail** sector,
3. and isolating the single remaining condition that separates the present result from a fully unconditional conservative 4PN theorem.

The decisive outcome is:

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}},
\]

with

- the **local** sector already assembled exactly in generic-frame ordinary form,
- the **tail** coefficient reduced to the same STF quadrupole normalization problem already isolated by the 2.5PN program,
- and therefore **no new independent 4PN normalization datum** opening beyond that 2.5PN quadrupole gap.

---

## 1. Claim taxonomy

This summary uses the same four-way claim taxonomy already adopted by the 1PN, 2PN, and 3PN writeups.

### 1.1 Exact

“Exact” means a symbolic identity once the declared lower-order ledger, target chart, and basis choices are fixed.

At 4PN this includes:

- the strict Schwarzschild one-body gate through \(c^{-8}\),
- the exact quartic perturbative Legendre compiler,
- the exact reduced COM local Hamiltonian/ordinary translation,
- the exact Hamiltonian-chart generic-frame local span result,
- the exact canonical-slice representative of the Hamiltonian-chart local residual,
- the exact generic-frame lift of the aligned ordinary seed correction,
- the exact local referee reconstruction of the full local ordinary target,
- and the exact coefficient bridge
  \[
  C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
  \]

### 1.2 Controlled reduction

“Controlled reduction” means the same already-declared reductions used below 4PN:

- projection / Poisson-hook / worldtube reduction from the 4D parent theory,
- strict one-body and COM reductions,
- fixed-chart local Hamiltonian import,
- and the grouped real quadrupole reduction inherited from the 2.5PN narrowing.

So nothing here is being claimed as an unreduced theorem of the full moving-throat PDE.

### 1.3 Protocol closure

“Protocol closure” means the result lives inside the same declared closure hierarchy already used successfully at 1PN, 2PN, and 3PN.

At 4PN this includes:

- carrying the lower-order self/static ledger forward,
- using the Hamiltonian-aligned seed convention for the ordinary 4PN local chart,
- and treating the hereditary coefficient as determined by the same passive/outgoing STF quadrupole normalization problem already isolated by 2.5PN.

### 1.4 Full assembly

“Full assembly” means the local theorem and the hereditary bridge have both been inserted into one consistent theorem envelope:

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]

What remains conditional is not another generic 4PN coefficient search. It is one sharply identified quadrupole-normalization theorem.

---

## 2. Lower-order inputs carried forward

The present 4PN statement is built on a lower-order backbone that is already frozen.

### 2.1 Conservative 1PN and 2PN status

The conservative 1PN and 2PN sectors are already solved **within a declared closure hierarchy**, not as assumption-free theorems of a fully solved moving-throat PDE. The carried conservative 1PN ledger includes

\[
\kappa_\rho=1,
\qquad
n=5,
\qquad
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3,
\]

with exact conservative EIH equality in the solved scope. The conservative 2PN sector is likewise closed against the generic-frame ADM target in the same hierarchy, and already sharpens the particle/throat ontology into a carried dipole wake, a real \(P_0\oplus P_2\) mouth/support layer, and a separate geometry-closure channel. See `4d_1pn_full_summary.md` and `4d_2pn_summary.md`.

### 2.2 Conservative 3PN status

The conservative 3PN two-body ledger is already fully assigned in the fixed ADM chart. Its exact split is

\[
\Delta \hat L_3^{\rm GR}
=
\underbrace{\Delta l_1 v^8}_{\text{compiler image of free COM Hamiltonian}}
+
\underbrace{L_{P_2}^{\rm mid}}_{\text{exact grouped real }P_2\text{ closure}}
+
\underbrace{\Delta l_{15}^{(g)}U^4}_{\text{unique geometry completion}}.
\]

So the 3PN bottleneck is no longer algebraic coefficient matching; the open issue is the deeper moving-throat derivation of the grouped-\(P_2\) and geometry lanes.

### 2.3 2.5PN narrowing that controls 4PN

The 2.5PN audit already narrowed the universal point-particle route to the **orbital/worldtube STF quadrupole**. The scalar and dipole branches were demoted above universal point-particle 2.5PN on the strict small-body branch, while the surviving retarded STF quadrupole branch was characterized as

\[
\mathcal M_{ab}^R(\omega)
=
\bigl[m_0+m_2\omega^2+m_4\omega^4+i\Gamma_5\omega^5+\cdots\bigr]\delta_{ab},
\qquad
m_0\neq 0,
\quad
\Gamma_5>0.
\]

The remaining gap was already narrowed to the **final passive/outgoing normalization** of that branch.

So 4PN should be read as the next conservative use of a quadrupole route that 2.5PN had already isolated.

---

## 3. Local 4PN theorem chain

## 3.1 Exact one-body 4PN gate

The strict isotropic Schwarzschild one-body target at 4PN is

\[
\frac{L_{4\mathrm{PN},1\text{-body}}}{m}
=
\frac{7}{256}\frac{v^{10}}{c^8}
+\frac{75}{128}\frac{Uv^8}{c^8}
+\frac{59}{16}\frac{U^2v^6}{c^8}
+\frac{203}{32}\frac{U^3v^4}{c^8}
+\frac{31}{32}\frac{U^4v^2}{c^8}
+\frac{1}{16}\frac{U^5}{c^8}.
\]

Continuing the 3PN self-sector packaging shows that 4PN is **not** a one-parameter continuation. The minimal direct-continuation one-body repair ledger is

\[
\boxed{
\mu_{\rho,4}=\frac18,
\qquad
 d_4=\frac{205}{16},
\qquad
 s_{34}=-\frac{15}{32},
\qquad
 s_{26}=-\frac1{16}.
}
\]

So 4PN opens a quartic static gate, a quartic denominator datum, and two genuinely new self-sector slots.

## 3.2 Exact quartic Legendre compiler

For

\[
L=L_0+\varepsilon L_1+\varepsilon^2L_2+\varepsilon^3L_3+\varepsilon^4L_4,
\]

with quadratic \(L_0\), constant Newtonian mass matrix \(M\), and

\[
v_0=M^{-1}p,
\quad
A_0=(\partial_vL_1)|_{v_0},
\quad
B_0=(\partial_vL_2)|_{v_0},
\quad
D_0=(\partial_vL_3)|_{v_0},
\]

\[
C_0=(\partial_v^2L_1)|_{v_0},
\quad
E_0=(\partial_v^2L_2)|_{v_0},
\quad
T_0=(\partial_v^3L_1)|_{v_0},
\]

the exact quartic compiler is

\[
H_4=
-L_4(v_0)
+A_0^TM^{-1}D_0
+\frac12B_0^TM^{-1}B_0
-B_0^TM^{-1}C_0M^{-1}A_0
-\frac12A_0^TM^{-1}E_0M^{-1}A_0
+\frac12A_0^TM^{-1}C_0M^{-1}C_0M^{-1}A_0
+\frac16T_0[M^{-1}A_0,M^{-1}A_0,M^{-1}A_0].
\]

This is the exact local Hamiltonian/ordinary bridge used throughout the rest of the 4PN chain.

## 3.3 Natural local seed and raw scaffold

Promoting the exact one-body target symmetrically gives the natural local self/static seed

\[
L_{4,\mathrm{seed}}=
\frac{7}{256}(m_A v_A^{10}+m_B v_B^{10})
+\frac{75Gm_A m_B}{128r}(v_A^8+v_B^8)
+\frac{59G^2m_A m_B}{16r^2}(m_B v_A^6+m_A v_B^6)
\]
\[
\qquad
+\frac{203G^3m_A m_B}{32r^3}(m_B^2 v_A^4+m_A^2 v_B^4)
+\frac{31G^4m_A m_B}{32r^4}(m_B^3 v_A^2+m_A^3 v_B^2)
+\frac{G^5m_A m_B}{16r^5}(m_B^4+m_A^4),
\]

understanding the whole expression as the coefficient of \(c^{-8}\).

The raw exchange-symmetric generic-frame local interaction scaffold then has block sizes

- \(52\) in the \(G/r\) degree-8 block,
- \(46\) in the \(G^2/r^2\) degree-6 block,
- \(29\) in the \(G^3/r^3\) degree-4 block,
- \(10\) in the \(G^4/r^4\) degree-2 block,
- \(2\) in the \(G^5/r^5\) static block,

for a total of \(139\) raw interaction directions.

After COM projection these blocks span exactly the full 15 interaction slots of the local reduced 4PN sector. So the local issue was never early span shortage; it was chart fixing and null-family control.

## 3.4 Why the lift had to be Hamiltonian-first

Importing the fixed local 4PN ADM target shows that the local instantaneous Hamiltonian has one free slot plus 15 interaction slots. Its strict test-mass limit matches the one-body 4PN Hamiltonian gate exactly.

The key structural lesson is:

- in the **Hamiltonian** chart the local interaction blocks have \(\nu\)-degree ceilings
  \[
  (4,3,3,2,2)
  \]
  across the \(G/r\), \(G^2/r^2\), \(G^3/r^3\), \(G^4/r^4\), \(G^5/r^5\) blocks,
- and these exactly match the constant-coefficient generic-frame scaffold.

By contrast, after quartic translation the **ordinary** target becomes

\[
(4,4,4,4,2),
\]

so the ordinary chart develops extra \(\nu^4\) tails in the middle and upper interaction blocks.

That proves the generic-frame lift must be done **in the Hamiltonian chart first**, with translation back to the ordinary chart only afterward.

## 3.5 Hamiltonian-chart generic-frame lift and canonical slice

In the Hamiltonian chart, the blockwise coefficient-space ranks are maximal:

\[
G/r:20/20,
\qquad
G^2/r^2:12/12,
\qquad
G^3/r^3:9/9,
\qquad
G^4/r^4:4/4,
\qquad
G^5/r^5:2/2.
\]

So there is **no Hamiltonian-chart span obstruction**. The full interaction coefficient-space rank is 47, with a 92-direction COM-blind null family.

That null family is purely algebraic. Its blockwise nullities are

\[
Q:32,
\qquad
T:34,
\qquad
S:20,
\qquad
U:6,
\qquad
W:0,
\]

with the \(Q\)-nulls lying in the full COM ideal and the \(T/S/U\)-nulls already lying in the linear-momentum ideal. The top \(W\)-block has no nullity at all.

Fixing the canonical quotient slice by setting all null coordinates to zero yields one sparse exact generic-frame Hamiltonian representative of the local comparable-mass residual.

## 3.6 Ordinary translation and the seed-alignment correction

Once the lower-order ledger is frozen and the 4PN seed is chosen in a Hamiltonian-aligned convention, the quartic residual compiler is an exact sign flip:

\[
\boxed{
\Delta H_4=-\Delta L_4.
}
\]

So the canonical Hamiltonian residual translates directly to a canonical ordinary residual.

The earlier ordinary-chart obstruction turned out not to come from this residual. It came from the **seed-alignment correction**

\[
\delta L_{4,\rm seed}
=
L_{4,\rm seed}^{(\rm aligned)}-L_{4,\rm seed}^{(\rm natural)}.
\]

This correction injects the extra \(\nu^4\) content into the ordinary \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks.

The exact lift theorem is that this seed-alignment correction is generic-frame liftable in the structured seed spaces

\[
Q,
\qquad
T\oplus(pq)T,
\qquad
S\oplus(pq)S,
\qquad
U\oplus(p^2q^2)U,
\qquad
W=0.
\]

So the true local ordinary obstruction was narrow and completely seed-sector in origin.

## 3.7 Final local theorem

Define

\[
L_{4,\mathrm{loc}}^{(\mathrm{gen})}
=
L_{4,\mathrm{seed}}^{(\mathrm{natural,gen})}
+
\delta L_{4,\mathrm{seed}}^{(\mathrm{gen})}
+
\Delta L_{4,\mathrm{loc}}^{(\mathrm{can,gen})}.
\]

Then the local referee master verifies that its COM reduction reproduces the full fixed-chart local ordinary 4PN target block by block.

So the local-first 4PN program is no longer blocked by local existence. The entire local instantaneous 4PN ordinary sector now has an exact generic-frame representative.

What remains locally is a sharper **uniqueness/interpretation** question for the aligned-seed sector, not another existence problem.

---

## 4. Tail / hereditary 4PN bridge

## 4.1 Exact GR tail structure

The conservative 4PN dynamics in GR splits into a local instantaneous sector plus a time-symmetric nonlocal tail sector. The nonlocal tail Hamiltonian is

\[
H_{4\rm PN}^{\rm tail,sym}(t)
=
-\frac{1}{5}\frac{G^2M}{c^8}
I_{ij}^{(3)}(t)
\,\mathrm{Pf}_{2s/c}\!\int_{-\infty}^{+\infty}\frac{dv}{|v|}
I_{ij}^{(3)}(t+v),
\]

with universal tail coefficient

\[
C_{\rm tail}^{\rm GR}=\frac{G^2M}{5c^8}.
\]

Its local logarithmic shadow is controlled by

\[
F_{\rm tail}(t)=\frac{2}{5}\frac{G^2M}{c^8}\bigl(I_{ij}^{(3)}(t)\bigr)^2.
\]

## 4.2 Universal Newtonian quadrupole shadow

Using the Newtonian order-reduced STF quadrupole

\[
I_{ij}=\mu\,\mathrm{STF}(x_i x_j),
\]

the exact third derivative is

\[
I_{ij}^{(3)}
=
-\frac{2GM\mu}{r^3}
\left(4x_{\langle i}v_{j\rangle}-3\frac{\dot r}{r}x_{\langle i}x_{j\rangle}\right),
\]

so

\[
\bigl(I_{ij}^{(3)}\bigr)^2
=
\frac{8}{3}\frac{G^2M^2\mu^2}{r^4}
\left(12v^2-11\dot r^2\right).
\]

Hence the local logarithmic tail shadow is

\[
\boxed{
\frac{F_{\rm tail}}{\mu}
=
\frac{16}{15}\,\nu\,\frac{U^4}{c^8}
\left(12v^2-11\dot r^2\right),
\qquad
U\equiv\frac{GM}{r}.
}
\]

This lives only in the degree-2 \(G^4/r^4\) block. That is why the hereditary side is much smaller than the local lift problem.

## 4.3 Exact bridge to the 2.5PN coefficient

Let the canonically normalized STF quadrupole odd coefficient be \(\gamma_{\rm quad}^{\rm eff}\). The exact coefficient bridge is

\[
\boxed{
C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
}
\]

In GR,

\[
\gamma_{\rm GR}=\frac{2G}{5c^5},
\qquad
C_{\rm tail}^{\rm GR}=\frac{G^2M}{5c^8},
\]

and indeed

\[
C_{\rm tail}^{\rm GR}=\frac{GM}{2c^3}\,\gamma_{\rm GR}.
\]

So the hereditary 4PN coefficient is exactly the 2.5PN quadrupole coefficient multiplied by the standard monopole-scattering factor.

## 4.4 Toy-model quadrupole branch forms

The 2.5PN program had already reduced the remaining quadrupole normalization problem to the canonically normalized branch quantity

\[
\gamma_{\rm quad}^{\rm eff}
=
\mathcal N_Q\frac{a^5}{27c_s^5}
=
\overline\Gamma_5
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]

Therefore the toy-model 4PN hereditary coefficient is

\[
\boxed{
C_{\rm tail}^{\rm toy}
=
\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}
=
\frac{GM}{2c^3}\,\overline\Gamma_5.
}
\]

On the Burke–Thorne target branch,

\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
\]

which implies automatically

\[
C_{\rm tail}^{\rm toy}=\frac{G^2M}{5c^8}=C_{\rm tail}^{\rm GR}.
\]

So **no new independent normalization datum opens at 4PN**.

---

## 5. Strongest honest current theorem statement

Within the declared closure hierarchy, the strongest current 4PN statement is:

### Theorem (conditional conservative 4PN statement)

Let

\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]

Then:

1. **Local sector.**
   The entire local instantaneous 4PN ordinary sector already has an exact generic-frame representative.

2. **Hereditary coefficient.**
   The unique compatible hereditary coefficient is
   \[
   C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
   \]

3. **No new 4PN normalization datum.**
   The 4PN tail sector introduces no new quadrupole-normalization constant beyond the same passive/outgoing STF quadrupole normalization already isolated by the 2.5PN audit.

4. **GR recovery condition.**
   If
   \[
   \gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
   \]
   then automatically
   \[
   C_{\rm tail}=\frac{G^2M}{5c^8},
   \]
   and the full conservative 4PN local+tail coefficient structure matches the standard GR target.

So the remaining gap between the present hierarchy and a fully unconditional conservative 4PN theorem is **exactly the same narrow passive/outgoing quadrupole normalization problem already isolated by the 2.5PN program**.

---

## 6. What is solved, what is conditional, what remains open

### 6.1 Already solved inside the present hierarchy

- the exact one-body 4PN gate,
- the exact quartic Legendre compiler,
- the full local 4PN ordinary generic-frame existence theorem,
- the canonical local Hamiltonian residual representative,
- the generic-frame lift of the aligned-seed correction,
- the exact local referee reconstruction,
- the universal local logarithmic tail shadow,
- and the exact coefficient bridge
  \[
  C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
  \]

### 6.2 Conditional part

The only conditional step left for a full conservative 4PN theorem is the same one already isolated on the 2.5PN side:

\[
\gamma_{\rm quad}^{\rm eff}\stackrel{?}{=}\frac{2G}{5c^5}
\]

on the natural passive/outgoing STF quadrupole branch.

### 6.3 Still open beyond the current hierarchy

- a first-principles moving-throat PDE derivation of the passive/outgoing quadrupole normalization,
- the deeper microscopic derivation of the grouped real conservative quadrupole data,
- the final uniqueness/interpretation theorem for the aligned-seed local sector,
- and any dissipative / radiative completion beyond the already narrowed quadrupole route.

---

## 7. The correct next theorem gate

The next clean target is not more local 4PN algebra. That part is already done.

The next actual theorem gate is:

\[
\boxed{
\text{derive the passive/outgoing STF quadrupole normalization on the moving-throat branch.}
}
\]

Equivalently, the task is to determine the canonical pair

\[
(\overline K_0,\overline K_2)
\]

on the natural branch. Once that pair is fixed,

- the normalized 2.5PN odd quadrupole coefficient closes,
- the normalized 4PN conservative hereditary coefficient closes,
- and the present conditional 4PN theorem becomes unconditional inside the same hierarchy.

---

## 8. Core package artifacts

The main 4PN package artifacts now are:

- `4pn_onebody_audit.py`
- `4pn_quartic_legendre_audit.py`
- `4pn_local_scaffold_audit.py`
- `4pn_local_target_import_audit.py`
- `4pn_local_hamiltonian_to_ordinary_audit.py`
- `4pn_hamiltonian_chart_generic_frame_lift_audit.py`
- `4pn_hamiltonian_chart_canonical_slice_audit.py`
- `4pn_generic_frame_ordinary_translation_audit.py`
- `4pn_generic_frame_aligned_seed_lift_audit.py`
- `4pn_local_referee_master_sympy_audit.py`
- `4pn_tail_bridge_audit.py`
- `4pn_tail_hereditary_bridge_audit.py`
- `4pn_conditional_referee_master_audit.py`

This summary is the human-readable theorem ledger for that package.
