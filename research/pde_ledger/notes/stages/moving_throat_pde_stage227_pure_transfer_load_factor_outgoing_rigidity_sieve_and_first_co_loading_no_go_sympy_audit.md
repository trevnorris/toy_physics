# Moving-Throat PDE — Stage 227: Pure-Transfer Load Factor, Outgoing-Rigidity Sieve, and the First Co-Loading No-Go

## Status

**Exact within the explicit finite-throat one-port weak-axisymmetric logarithmic-slope closure built on the Stage-223 compatibility branch and the Stage-226 pure-transfer subcorridor.**

This stage does **not** solve the full moving-throat PDE.
It takes the strict same-charge survivor from Stage 226,
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\Xi_{\rm load}=\frac{N_{01}}{N_0},
\]
and asks the next natural microscopic question:

> if the conservative one-pole bundle is frozen at first order, what one-port outgoing load factor is actually carrying the same-charge signal, and which rigidity assumptions kill that mechanism?

---

## Purpose

Stage 226 already reduced the live same-charge corridor to a very specific mechanism:
\[
D_{01}=D_{21}=D_{41}=0,
\qquad
\Xi_1=\Xi_{\rm load}=\frac{N_{01}}{N_0}.
\]

So the next honest step is no longer another generic mixed-anisotropy scan.
The surviving effect is already telling us what it is:

> it is a **pure outgoing-load** problem.

The main outputs of this stage are:

1. the exact pure-transfer theorem
   \[
   \Xi_1=\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta};
   \]
2. the exact one-port load-factor decomposition
   \[
   \Lambda:=\frac{P}{\Delta}
   =
   \frac{G_W}{\Omega_W^2}\,\frac{1+I}{1-H},
   \]
   with
   \[
   I=\frac{RG_U}{\Omega_U^2G_W},
   \qquad
   H=\frac{R^2}{\Omega_U^2\Omega_W^2};
   \]
3. the exact pure-transfer identity
   \[
   \Xi_1=2\,\delta\ln\Lambda;
   \]
4. the outgoing-rigidity sieve:
   \[
   i=0,\qquad h=0,\qquad m=0,\qquad i=h=0;
   \]
5. the first exact same-charge **co-loading no-go**
   \[
   i=h=0
   \quad\Longrightarrow\quad
   \text{only the trivial pure-transfer drift survives;}
   \]
6. and the transported dynamic ceilings on the remaining one-dimensional rigid survivors.

So after Stage 227, the question is no longer

> does some pure-transfer effect survive?

It is now

> which outgoing-load channels are allowed to move together, and which rigidity assumptions kill the same-charge corridor outright?

---

## 1. Frozen input carried forward

### 1.1 Concrete finite-throat one-port sample branch

Keep the same explicit finite-throat branch from Stages 223–226:

- wall and brane-like internal coordinate on the lowest N/N zero mode,
- trapped support and mixed coordinate on the lowest D/N half-wave,
- exact overlap constant
  \[
  \kappa=\frac{2\sqrt2}{\pi}.
  \]

The primitive parameters remain
\[
(\lambda_B,\lambda_U,\lambda_W,\lambda_R,\Omega_U,\Omega_W,\varpi,M)
=
\left(\frac12,\frac3{10},\frac25,\frac14,1,\frac75,2,1\right).
\]

The one-port mixed primitives are
\[
G_U=\lambda_U,
\qquad
G_W=\kappa\lambda_W,
\qquad
R=\kappa\lambda_R,
\]
\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\qquad
P=\Omega_U^2G_W+RG_U,
\qquad
N_0=\frac{P^2}{\Delta^2}.
\]

### 1.2 Stage-226 pure-transfer corridor

Stage 226 showed that on this noncanonical compatibility branch, the full intersection of

1. Stage-225 conservative-shape preservation, and
2. the stricter imported `5`PN even-gate package,

is equivalent to
\[
D_{01}=D_{21}=D_{41}=0.
\]

So on that two-dimensional mixed-sector subcorridor the same-charge scalar becomes
\[
\boxed{
\Xi_1=\frac{N_{01}}{N_0}.
}
\]

That is the starting point of this stage.

---

## 2. Exact pure-transfer load theorem

Differentiate
\[
N_0=\frac{P^2}{\Delta^2}.
\]
Then on any weak first-order deformation,
\[
\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta}.
\]
Because the pure-transfer corridor has
\[
D_{01}=0,
\]
the surviving same-charge scalar is exactly
\[
\boxed{
\Xi_1
=
\frac{N_{01}}{N_0}
=
2\frac{P_{01}}{P}
-
2\frac{\Delta_{01}}{\Delta}.
}
\]

So the surviving same-charge mechanism is literally the logarithmic slope of the one-port load factor
\[
\Lambda:=\frac{P}{\Delta}.
\]

Equivalently,
\[
\boxed{
\Xi_1=2\,\delta\ln\Lambda.
}
\]

This is already a major sharpening of the Stage-226 verdict.
The mechanism is no longer “mixed-sector enhancement” in general.
It is one exact load-factor slope.

---

## 3. Exact one-port factorization

Write
\[
P=\Omega_U^2G_W+RG_U,
\qquad
\Delta=\Omega_U^2\Omega_W^2-R^2.
\]
Then
\[
P=\Omega_U^2G_W\left(1+\frac{RG_U}{\Omega_U^2G_W}\right),
\qquad
\Delta=\Omega_U^2\Omega_W^2\left(1-\frac{R^2}{\Omega_U^2\Omega_W^2}\right).
\]
So the load factor becomes
\[
\boxed{
\Lambda
=
\frac{P}{\Delta}
=
\frac{G_W}{\Omega_W^2}\,
\frac{1+I}{1-H},
}
\]
with exact one-port invariants
\[
\boxed{
I=\frac{RG_U}{\Omega_U^2G_W},
\qquad
H=\frac{R^2}{\Omega_U^2\Omega_W^2}.
}
\]

Define the microscopic logarithmic drifts
\[
m:=\delta\ln\!\left(\frac{G_W}{\Omega_W^2}\right),
\qquad
i:=\delta\ln I,
\qquad
h:=\delta\ln H.
\]
In the primitive mixed slope variables
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}),
\]
these are
\[
m=x_{\lambda_W}-2x_{\Omega_W},
\]
\[
i=x_{\lambda_R}+x_{\lambda_U}-x_{\lambda_W}-2x_{\Omega_U},
\]
\[
h=2x_{\lambda_R}-2x_{\Omega_U}-2x_{\Omega_W}.
\]

Then on the Stage-226 pure-transfer corridor,
\[
\boxed{
\Xi_1
=
2\left[
m+\frac{I}{1+I}\,i+\frac{H}{1-H}\,h
\right].
}
\]

So the outgoing-load scalar has three transparent pieces:

1. a direct mixed-leg / port-weight piece \(m\),
2. an interference-ratio piece \(i\),
3. a hybridization-ratio piece \(h\).

### 3.1 Exact sample-branch coefficients

On the concrete sample branch,
\[
I=\frac{3}{16},
\qquad
H=\frac{25}{98\pi^2}.
\]
So the exact Stage-227 load law is
\[
\boxed{
\Xi_1
=
2m+\frac{6}{19}i+\frac{50}{98\pi^2-25}\,h
}
\]
on the pure-transfer corridor.

This is the first exact microscopic decomposition of the surviving same-charge scalar.

---

## 4. The outgoing-rigidity sieve

Now impose the first natural rigidity filters.

### 4.1 Combined interference and hybridization rigidity

If both
\[
i=0,
\qquad
h=0,
\]
are imposed on the pure-transfer corridor, the exact reduced `2 x 2` rigidity matrix on the Stage-226 basis has nonzero determinant:
\[
\det[(i,h)|_{\rm pure\ transfer}]
=
-\frac{19(-25+98\pi^2)(200+147\pi^2)(441\pi^2+4400)}
{6(8670000+14894275\pi^2+2117682\pi^4)}
\neq 0.
\]

So the combined rigidity system has only the trivial solution:
\[
\boxed{
D_{01}=D_{21}=D_{41}=0,\quad i=0,\quad h=0
\ \Longrightarrow\
x_{\rm mixed}=0.
}
\]

This is the first exact same-charge co-loading no-go of the audit.

It says:

> if both the interference ratio and the hybridization ratio are frozen, the pure-transfer same-charge mechanism dies on this concrete branch.

### 4.2 Single-rigidity survivors

The situation is different if only one rigidity is imposed.

The audit finds:

- \(i=0\) leaves a one-dimensional survivor,
- \(h=0\) leaves a one-dimensional survivor,
- \(m=0\) also leaves a one-dimensional survivor.

So the corridor is narrower, but it does not die under a **single** outgoing-rigidity filter.

This is useful physically.
It means the same-charge branch still has room to live if one of the outgoing subchannels is rigid, but not if both \(i\) and \(h\) are rigid simultaneously.

---

## 5. Concrete unit directions and same-charge gain scales

Using ambient Euclidean normalization in the mixed primitive space
\[
(x_{\lambda_U},x_{\lambda_W},x_{\lambda_R},x_{\Omega_U},x_{\Omega_W}),
\]
the audit gives the following unit directions.

### 5.1 `i = 0` survivor
\[
v_i\approx
(0.45280825,\,-0.29424612,\,-0.82815170,\,-0.04054866,\,0.14458380),
\]
with
\[
|\Xi_1(v_i)|\approx 1.26576248.
\]

### 5.2 `h = 0` survivor
\[
v_h\approx
(0.66561963,\,-0.38941932,\,0.46712837,\,0.03609301,\,0.43103536),
\]
with
\[
|\Xi_1(v_h)|\approx 2.04509123.
\]

### 5.3 `m = 0` survivor
\[
v_m\approx
(0.13386239,\,-0.10586713,\,-0.98242900,\,-0.05389175,\,-0.05293356),
\]
with
\[
|\Xi_1(v_m)|\approx 0.29342952.
\]

The `m = 0` result is especially interesting.
It means the direct mixed-leg factor can be frozen while the same-charge signal is carried entirely by a correlated interference/hybridization deformation.

So the pure-transfer corridor is not synonymous with a raw mixed-leg effect.

---

## 6. Transported dynamic ceilings on the rigid survivors

Interpret the ambient microscopic mixed-sector drift amplitude as
\[
\|x_{\rm mixed}\|_2=t.
\]
If the operator norm of `\Xi_1` on a corridor is `\sigma`, then the transported Stage-224 ceiling law becomes
\[
|\epsilon|t \le \frac{\text{budget}}{\sigma}.
\]

### 6.1 Reference pure-transfer corridor
\[
\sigma_{\rm transfer}\approx 2.31561904.
\]
So the stricter `10%`-loss ceilings are
\[
|\epsilon|t \lesssim 0.15889070
\qquad
(\text{both wall-like poles survive}),
\]
\[
|\epsilon|t \lesssim 0.31854077
\qquad
(\text{nonempty wall-like corridor}).
\]

### 6.2 `i = 0` survivor
\[
\sigma_i\approx 1.26576248.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 0.29067881,
\qquad
|\epsilon|t \lesssim 0.58274682.
\]

### 6.3 `h = 0` survivor
\[
\sigma_h\approx 2.04509123.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 0.17990900,
\qquad
|\epsilon|t \lesssim 0.36067783.
\]

### 6.4 `m = 0` survivor
\[
\sigma_m\approx 0.29342952.
\]
So the stricter `10%`-loss ceilings become
\[
|\epsilon|t \lesssim 1.25389678,
\qquad
|\epsilon|t \lesssim 2.51378617.
\]

The `m = 0` branch leaves the most headroom simply because the same-charge scalar is much smaller per unit ambient microscopic drift there.

That does **not** make it automatically the physically best mechanism, but it does make it the least constrained by the transported dynamic ceiling.

---

## 7. What Stage 227 changes

Before this stage, the strongest statement was only

> the strict even-gate survivor is a pure-transfer mixed corridor.

After this stage, the picture is much sharper.

1. The surviving same-charge scalar is exactly one **outgoing load-factor slope**:
   \[
   \Xi_1=2\,\delta\ln\Lambda.
   \]
2. On the concrete branch, that load factor splits into three pieces:
   \[
   m,\qquad i,\qquad h.
   \]
3. Freezing both \(i\) and \(h\) kills the corridor outright.
4. Freezing only one of \(i,h,m\) still leaves a one-dimensional same-charge survivor.
5. Even with \(m=0\), the same-charge effect can survive through correlated interference/hybridization motion.

So the best current summary is:

> the idea survives Stage 227, but only as a very structured outgoing co-loading effect. Pure transfer is real, but it is not generically “just the mixed leg.” And the first exact no-go now says that simultaneous interference and hybridization rigidity kills the mechanism on this concrete branch.

That is a real narrowing, and it is the right kind of narrowing.

---

## 8. Script-backed status

The accompanying SymPy audit verifies:

- the exact Stage-226 pure-transfer corridor on the explicit finite-throat compatibility branch;
- the exact pure-transfer theorem
  \[
  \Xi_1=\frac{N_{01}}{N_0}=2\frac{P_{01}}{P}-2\frac{\Delta_{01}}{\Delta};
  \]
- the exact one-port factorization
  \[
  \Lambda=\frac{G_W}{\Omega_W^2}\frac{1+I}{1-H};
  \]
- the exact microscopic slope law
  \[
  \Xi_1=2\left[m+\frac{I}{1+I}i+\frac{H}{1-H}h\right]
  \]
  and its explicit sample-branch specialization
  \[
  \Xi_1=2m+\frac{6}{19}i+\frac{50}{98\pi^2-25}h;
  \]
- the nonzero pure-transfer rigidity determinant for the combined `i=h=0` filter;
- the one-dimensional unit survivors for `i=0`, `h=0`, and `m=0`;
- the corresponding same-charge gain scales \(\sigma_i,\sigma_h,\sigma_m\);
- and the transported stricter `10%`-loss dynamic ceilings on the reference pure-transfer corridor and on each rigid survivor.

Supporting file:
- `moving_throat_pde_stage227_pure_transfer_load_factor_outgoing_rigidity_sieve_and_first_co_loading_no_go_sympy_audit.py`
