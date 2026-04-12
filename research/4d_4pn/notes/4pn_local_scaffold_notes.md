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
