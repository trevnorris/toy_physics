# Stage V2-20 — Weak-form and numerical branch-extraction preparation

## Status

This stage is a **pre-solver theorem gate**, not a nonlinear branch solution.

The previous stages collapsed the isotropic grouped-\(P_2\) finish line to one scalar normalization test,

\[
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
=
\frac{54Gc_s^5}{5a^5c^5},
\]

together with conservative one-pole, constant-prefactor, tail-transport, stability, and freeze-order gates.

V2-20 answers the next practical question:

> When a future moving-throat solver outputs a branch, exactly what weak-form data must be extracted, and how do we turn those data into \(K,M,B_n,Z_n,N_n,\widehat m_0,\Theta_{\rm tail}\) without refitting the target?

## Architectural patch retained from V2-04

The throat is now treated as an **open finite-radius conduit**.

The exit data are:

\[
R(L)>0,
\]

not

\[
R(L)=0.
\]

The AC support modes may see an effective Neumann condition through an open-exit impedance mismatch, but this is not a hard geometric cap. The zero/DC superfluid current is allowed to exit and is recorded through the open-system leakage bookkeeping rather than through the finite-support AC mode equation.

In this stage, the far end is represented by a Robin/impedance condition,

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0.
\]

The Neumann support limit is

\[
Y_L(\omega)\to0.
\]

For a low-impedance open expansion, this is the support-coordinate version of the organ-pipe reflection result from V2-04.

## 1. Weak form for one wall harmonic

For a real harmonic sector \((l,m)\), write

\[
V_l(w)=K_\eta(w)+l(l+1)T_\Omega(w).
\]

The densitized wall equation is

\[
\mu_\eta q_{tt}
-\partial_w(T_w q_w)
+V_l q
=
S_l.
\]

Using \(q(w,t)=q(w,\omega)e^{-i\omega t}\), this becomes

\[
-\omega^2\mu_\eta q
-\partial_w(T_w q_w)
+V_lq
=
S_l.
\]

For a test function \(\varphi\) satisfying the mouth constraint, the weak form is

\[
\int_0^L
\left[
T_w q_w\varphi_w
+
(V_l-\omega^2\mu_\eta)q\varphi
\right]dw
+
Y_L(\omega)q(L)\varphi(L)
=
\int_0^L S_l\varphi\,dw.
\]

The SymPy script verifies the integration-by-parts identity

\[
Tq_w\varphi_w-\partial_w(Tq_w\varphi)
=
-\partial_w(Tq_w)\varphi.
\]

It also verifies the open-exit sign convention: the endpoint variation is

\[
\left[T_w q_w(L)+Y_Lq(L)\right]\varphi(L),
\]

so free exit variations impose

\[
T_w q_w(L)+Y_Lq(L)=0.
\]

## 2. Galerkin extraction matrices

Choose a frozen branch basis \(\{\chi_i(w)\}\) satisfying the mouth condition. In a numerical solver this basis may be finite elements, splines, spectral functions, or computed eigenfunctions. The extracted matrices are

\[
M_{ij}^{(l)}
=
\int_0^L\mu_\eta(w)\chi_i(w)\chi_j(w)\,dw,
\]

\[
K_{ij}^{(l)}
=
\int_0^L
\left[
T_w(w)\chi_i'(w)\chi_j'(w)
+
V_l(w)\chi_i(w)\chi_j(w)
\right]dw
+
Y_L(0)\chi_i(L)\chi_j(L).
\]

The script checks a two-function Dirichlet-mouth open-exit toy basis,

\[
\chi_1=w/L,\qquad \chi_2=(w/L)^2,
\]

and derives

\[
M
=
\begin{pmatrix}
L\mu/3 & L\mu/4\\
L\mu/4 & L\mu/5
\end{pmatrix},
\qquad
\det M=\frac{L^2\mu^2}{240}>0.
\]

The stiffness matrix is symmetric and is positive in the tested positive-coefficient sample.

For the AC support limit \(Y_L\to0\), the same script checks the D/N support seed,

\[
q(w)=\sin\left(\frac{\pi w}{2L}\right),
\]

with

\[
q(0)=0,\qquad q_w(L)=0.
\]

This keeps the support ladder while avoiding a capped geometry.

## 3. Wall/BdG/Maxwell extraction slots

For each grouped lane

\[
A\in\{20,21,22\},
\]

the solver must output or permit extraction of the following frozen primitives.

### 3.1 Wall data

For the selected wall/support coordinate,

\[
M_A,
\qquad
K_A.
\]

These are obtained from the weak-form matrices and the selected normalized wall eigenvector or response coordinate.

### 3.2 Stable BdG support data

For positive-energy stable BdG modes,

\[
c_{A\alpha},
\qquad
\varpi_{A\alpha}>0.
\]

The extracted moments are

\[
B_{A0}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^2},
\]

\[
B_{A2}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^4},
\]

\[
B_{A4}=\sum_\alpha\frac{c_{A\alpha}^2}{\varpi_{A\alpha}^6}.
\]

The script verifies the low-frequency Schur-complement expansion

\[
K-M\omega^2-\frac{g^2}{\varpi^2-\omega^2}
=
\left(K-\frac{g^2}{\varpi^2}\right)
-
\left(M+\frac{g^2}{\varpi^4}\right)\omega^2
-
\frac{g^2}{\varpi^6}\omega^4
+\cdots.
\]

It also verifies the two-mode Hankel/Stieltjes positivity identity

\[
B_0B_4-B_2^2
=
\frac{
w_1w_2(\lambda_1-\lambda_2)^2
}{
\lambda_1^3\lambda_2^3
}
\ge0.
\]

### 3.3 Conservative Maxwell/mixed data

For each mixed port \(r\), the numerical branch must output

\[
\Omega_{U,Ar},\quad
\Omega_{W,Ar},\quad
R_{Ar},\quad
g_{U,Ar},\quad
g_{W,Ar}.
\]

Define

\[
\Delta_{Ar}=\Omega_{U,Ar}^2\Omega_{W,Ar}^2-R_{Ar}^2,
\]

\[
S_{Ar}=\Omega_{U,Ar}^2+\Omega_{W,Ar}^2,
\]

\[
Q_{Ar}
=
g_{U,Ar}^2\Omega_{W,Ar}^2
+
2g_{U,Ar}g_{W,Ar}R_{Ar}
+
g_{W,Ar}^2\Omega_{U,Ar}^2,
\]

\[
H_{Ar}=g_{U,Ar}^2+g_{W,Ar}^2.
\]

Then

\[
Z_{A0}^{(r)}=\frac{Q_{Ar}}{\Delta_{Ar}},
\]

\[
Z_{A2}^{(r)}
=
\frac{Q_{Ar}S_{Ar}-H_{Ar}\Delta_{Ar}}{\Delta_{Ar}^2},
\]

\[
Z_{A4}^{(r)}
=
\frac{
Q_{Ar}(S_{Ar}^2-\Delta_{Ar})-S_{Ar}H_{Ar}\Delta_{Ar}
}{
\Delta_{Ar}^3
}.
\]

The full conservative mixed data are

\[
Z_{An}=\sum_r Z_{An}^{(r)}.
\]

The internal stability gate is

\[
\Delta_{Ar}>0
\]

for every active mixed port.

### 3.4 Outgoing transfer data

For each outgoing mixed port,

\[
P_{Ar}
=
\Omega_{U,Ar}^2g_{W,Ar}
+
R_{Ar}g_{U,Ar}.
\]

The outgoing-transfer moments are

\[
N_{A0}^{(r)}
=
\frac{P_{Ar}^2}{\Delta_{Ar}^2},
\]

\[
N_{A2}^{(r)}
=
\frac{
2P_{Ar}(P_{Ar}S_{Ar}-\Delta_{Ar}g_{W,Ar})
}{
\Delta_{Ar}^3
},
\]

\[
N_{A4}^{(r)}
=
\frac{
\Delta_{Ar}^2g_{W,Ar}^2
-2\Delta_{Ar}P_{Ar}^2
-4\Delta_{Ar}P_{Ar}S_{Ar}g_{W,Ar}
+3P_{Ar}^2S_{Ar}^2
}{
\Delta_{Ar}^4
}.
\]

The full transfer data are

\[
N_{An}=\sum_rN_{An}^{(r)}.
\]

The dark-port failure is

\[
P_{Ar}=0
\]

for all active ports, which gives

\[
N_{A0}=0.
\]

## 4. Total grouped operator and response extraction

The conservative lane operator is

\[
D_A(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+O(\omega^6),
\]

with

\[
D_{A0}=K_A-B_{A0}-Z_{A0},
\]

\[
D_{A2}=-(M_A+B_{A2}+Z_{A2}),
\]

\[
D_{A4}=-(B_{A4}+Z_{A4}).
\]

Then

\[
u_{2,A}=-\frac{D_{A2}}{D_{A0}},
\]

\[
u_{4,A}=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2}.
\]

The outgoing prefactor is

\[
P_{0,A}=\frac{N_{A0}}{D_{A0}},
\]

\[
P_{2,A}=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},
\]

\[
P_{4,A}
=
\frac{
D_{A0}^2N_{A4}
-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})
+3D_{A2}^2N_{A0}
}{D_{A0}^3}.
\]

The script verifies these formulas algebraically and confirms that the constant-prefactor branch is

\[
N_2=\frac{2D_2N_0}{D_0},
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}.
\]

## 5. Grouped projectors and anisotropy extraction

The grouped metric is

\[
G_{\rm grp}=\mathrm{diag}(1,2,2).
\]

For any grouped triple \(x=(x_{20},x_{21},x_{22})\),

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}5,
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}2.
\]

The inverse map is

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

The script verifies projector completeness and idempotence.

For a pure weak-axisymmetric perturbation,

\[
\lambda=(1,\tfrac12,-1),
\]

the script verifies

\[
\bar x=x_0,
\qquad
b_x=3a_x.
\]

So future branch output must be classified as:

- isotropic if \(a=b=0\),
- pure weak-axisymmetric if \(b=3a\),
- more general anisotropy otherwise.

## 6. Final target residual packet

The numerical branch extraction should report these residuals:

\[
R_{\rm pole}
=
D_0(B_4+Z_4)-3(M+B_2+Z_2)^2,
\]

\[
R_{\rm norm}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}\frac{N_0}{D_0}
-
\frac{54Gc_s^5}{5a^5c^5},
\]

\[
R_{P2}=P_2,
\qquad
R_{P4}=P_4,
\]

\[
R_{\rm tail}
=
\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

The script verifies that

\[
R_{\rm norm}=0
\]

is equivalent to

\[
\gamma_{\rm quad}^{\rm eff}
=
\frac{2G}{5c^5}.
\]

## 7. Mandatory stability gates

A branch cannot be accepted unless the following are true before target evaluation:

\[
\mu_\eta>0,
\qquad
T_w>0,
\]

\[
\Delta_{Ar}>0
\quad
\text{for every active mixed port},
\]

\[
D_0>0,
\]

\[
B_0B_4-B_2^2\ge0
\]

for the stable BdG moment sequence, and

\[
B_4+Z_4>0
\]

on a stable nondegenerate one-pole branch.

The branch also fails for the quadrupole route if the active outgoing ports are all dark:

\[
N_0=0.
\]

## 8. Freeze manifest

The script emits a freeze manifest requiring:

1. open finite-radius exit geometry;
2. no hard cap;
3. mouth boundary class;
4. exit impedance boundary class;
5. DC leakage policy;
6. weak-form convention;
7. densitized measure convention;
8. extraction slot list;
9. stability gates;
10. target residual packet;
11. no-refit rule.

The manifest hash produced by the script is

```text
067c626bb61456fad945f5b3f7fa4d10c19e38f9083bd835896e1d052261e390
```

This hash is not a physical invariant. It is a bookkeeping certificate for the V2-20 extraction schema. A future branch run should include its own branch-packet hash before target residuals are evaluated.

## 9. Solver protocol produced by this stage

A future numerical PDE branch test should run in this order:

1. **Freeze the packet**: parent action, gauge convention, open-exit boundary class, source projection, wall action, support basis, port list, and extraction formulas.
2. **Solve the stationary branch** with \(R(L)>0\).
3. **Linearize** the GNLS/BdG, Maxwell/mixed, and wall/interface equations.
4. **Select positive-energy stable modes** and record any excluded negative-Krein modes.
5. **Extract wall matrices** \(K_A,M_A\).
6. **Extract BdG moments** \(B_{A0},B_{A2},B_{A4}\).
7. **Extract mixed-sector moments** \(Z_{A0},Z_{A2},Z_{A4}\).
8. **Extract outgoing transfer moments** \(N_{A0},N_{A2},N_{A4}\).
9. **Project grouped data** into \(\bar x,a_x,b_x\).
10. **Compute residuals** \(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4},R_{\rm tail}\).
11. **Report pass/fail** without changing the frozen branch packet.

## 10. SymPy result

The script reports:

```text
checks_total: 28
checks_passed: 28
checks_failed: 0
```

The algebraic scaffold passes.

The important limitation is unchanged:

> V2-20 prepares the branch extraction theorem gate. It does not prove that the actual moving-throat PDE branch lands on the target surface.
