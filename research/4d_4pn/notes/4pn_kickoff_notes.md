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
