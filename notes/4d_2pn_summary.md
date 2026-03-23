4d_2pn.tex — Comprehensive summary (full conservative 2PN assembly paper)

## 0) What this paper is doing

This paper is the **full conservative two-body 2PN assembly** for the 4D toy-model program.

It starts from three things that are already frozen before the paper begins:

1. the exact 4D parent reduction backbone from `4d.tex`:
   - exact projection map,
   - exact projected continuity with leakage,
   - exact longitudinal identity,
   - controlled quasi-static Poisson hook,
   - controlled coherent-defect / worldtube particle reduction;

2. the full conservative Newtonian + 1PN two-body assembly from `4d_1pn_full.tex`, including exact EIH match;

3. the corrected Paper 7 / Paper 8 ontology and notation firewall:
   - electric charge sign is the topological puncture orientation
     \[
     \eta_Q=\pm 1,
     \]
   - microscopic charge branch
     \[
     q_* = \eta_Q e_*,
     \qquad e_*>0,
     \]
   - canonically normalized brane charge
     \[
     q_{\rm eff} = \frac{q_*}{\sqrt{Z_{\rm int}}},
     \qquad
     e_{\rm eff} = \frac{e_*}{\sqrt{Z_{\rm int}}},
     \]
   - circulation is **not** the definition of electric charge;
   - historical gravity-side bare `q=1` is now written 
     \[
     \kappa_\rho=1.
     \]

The paper’s main target is the complete conservative near-zone two-body ledger through order
\[
c^{-4},
\]
i.e. the full conservative **2PN** sector.

That means closing all of the following at once:

- free-particle sextic kinematics \(v^6\),
- one-body / self terms,
- comparable-mass quartic-velocity cross terms at order \(G/r\),
- quadratic-velocity terms at order \(G^2/r^2\),
- and the static cubic term at order \(G^3/r^3\).

The paper’s strongest statement is:

> **Within a stated closure hierarchy, the derived two-body Lagrangian Legendre-transforms exactly to the standard generic-frame ADM Hamiltonian through 2PN.**

This is not claimed as an assumption-free theorem of the fully solved moving-throat PDE. It is a **full conservative 2PN derivation within a declared closure hierarchy**.

---

## 1) Claim taxonomy and how to read the paper

The paper uses the same four claim classes as the 1PN full paper.

### 1.1 Exact

These are identities or exact algebraic steps once the relevant inputs are declared.

Examples:
- carried exact projection formulas,
- carried exact longitudinal identity,
- exact perturbative Legendre-transform identity,
- exact residual solve once the basis is fixed,
- exact final equality to ADM after full assembly.

### 1.2 Controlled reduction

These use an explicit regime assumption or reduced description.

Examples:
- quasi-static Poisson hook,
- coherent-defect / small-body worldtube reduction,
- low-frequency DtN packaging,
- Bernoulli exactification of the self sector,
- quadratic local mass scaling,
- Family-1 wall reductions.

### 1.3 Protocol closure

These are fixed only **within a declared response hierarchy** rather than by a fully solved bulk PDE.

Examples:
- the one-body invariant denominator ansatz,
- the carried adiabatic breathing closure from 1PN,
- the specific low-frequency throat-response packaging used in the appendices.

### 1.4 Full assembly

Once all declared ingredients are inserted, the final equality to the standard target is exact algebra.

That is the status of:
\[
L_{\rm cons}^{\rm derived}
= L_0 + c^{-2}L_1 + c^{-4}L_2 + O(c^{-6}),
\qquad
H_2 = H_{2\mathrm{PN}}^{\rm ADM}.
\]

---

## 2) What the paper claims — and what it does not claim

### 2.1 What it claims

The paper claims all of the following.

1. The 4D projection / Poisson-hook / worldtube backbone carries forward unchanged.
2. The missing one-body 2PN response slot is fixed by a minimal DtN-invariant denominator that preserves the frozen 1PN slot.
3. Minimal Bernoulli exactification plus quadratic local mass scaling gives the self/static 2PN seed.
4. The remaining comparable-mass discrepancy is a finite algebraic residual that solves uniquely in a compact invariant basis.
5. The assembled conservative two-body 2PN Lagrangian Legendre-transforms exactly to the generic-frame ADM target.
6. The new 2PN cross content already admits a constructive throat-response interpretation:
   - carried odd dipole wake,
   - new even \(P_0\oplus P_2\) support sector,
   - separate geometry-closure channel,
   - and an explicit monopole breathing response module.

### 2.2 What it does **not** claim

The paper explicitly does **not** claim:

- a fully solved moving-throat bulk PDE theorem;
- an assumption-free derivation of every new 2PN response datum;
- any dissipative, radiative, or non-conservative result;
- spin couplings, radiation reaction, or higher-PN completion;
- a completed dynamic derivation of all appendix-side pole scales.

So the correct reading is:

> the conservative 2PN two-body sector is closed **within the declared hierarchy**, while the fully dynamical throat-response problem remains open.

---

## 3) Charge ontology / notation firewall

This paper is downstream of the corrected 4D / EM papers.

Important notational firewall:

1. Electric charge is **not** part of the main conservative 2PN proof.
2. When charge is mentioned, it uses the corrected notation
   \[
   \eta_Q = \pm 1,
   \qquad
   q_* = \eta_Q e_*,
   \qquad
   q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}}.
   \]
3. Quantized circulation belongs to the magnetic/vortical sector, not to electric charge.
4. The historical gravity-side scalar coefficient once written as bare `q=1` is now
   \[
   \kappa_\rho=1.
   \]
5. The comparable-mass residual coefficients
   \[
   q_1,\dots,q_7
   \]
   in this paper are **not** electric charges; they are just residual-basis coefficients.
6. The appendix-side vector
   \[
   J = \left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right)
   \]
   is a **throat-response source-weight vector**, not an electromagnetic current.
7. There are two different scalar gap objects later in the appendices:
   - pure-potential geometry-closure deficit
     \[
     \Delta_{\rm geom}=\frac{281}{80},
     \]
   - monopole wall add-on
     \[
     \Delta K_{00}=\frac{109}{280}.
     \]
   These are **not the same object**.

---

## 4) Headline outputs / carry-forward constants

This is the fastest “memory ledger” for future work.

### 4.1 Carried lower-order backbone

Frozen from earlier papers:
\[
\kappa_\rho = 1,
\qquad
n=5,
\qquad
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3.
\]

Frozen parity-even wake coefficients:
\[
C_\parallel = -\frac72,
\qquad
C_L = -\frac12.
\]

Frozen Newtonian and 1PN two-body blocks:
\[
L_0 = \frac12 m_A v_A^2 + \frac12 m_B v_B^2 + \frac{Gm_A m_B}{r},
\]
\[
L_1 =
\frac18(m_Av_A^4+m_Bv_B^4)
+\frac{Gm_A m_B}{r}
\left[
\frac32(v_A^2+v_B^2)-\frac72 v_{AB}-\frac12 v_{An}v_{Bn}
\right]
-\frac{G^2m_A m_B(m_A+m_B)}{2r^2}.
\]

### 4.2 One-body 2PN denominator closure

Exact isotropic Schwarzschild target through 2PN requires
\[
D_{\rm target}(u)=1-4u+8u^2+O(u^3),
\qquad
u\equiv U/c^2.
\]

The paper closes this with the minimal invariant denominator
\[
D_{\rm eff}(u)
=
C_s(u)\left[1+\mu(G(u)-1)^2\right],
\qquad
\mu = \frac{32768}{3249},
\]
so that
\[
D_{\rm eff}(u)=1-4u+8u^2+O(u^3).
\]

Important sub-results:
\[
\alpha=0,
\qquad
g_1=\frac{57}{64},
\qquad
\mu = \frac{8}{g_1^2} = \frac{32768}{3249}.
\]

### 4.3 Static 2PN gate

Quadratic local mass scaling plus the pure-static isotropic test-mass gate fix
\[
\lambda_\rho = \frac12.
\]

Hence the two-body static 2PN block is
\[
L_{2,\rm static}
=
\frac{G^3m_A m_B(m_A^2+m_B^2)}{4r^3}.
\]

### 4.4 Frozen self/static seed

The self/static input to the comparable-mass lift is
\[
L_{2,\mathrm{self+static}}=
\frac{m_A v_A^6+m_B v_B^6}{16}
+\frac{7Gm_A m_B}{8r}(v_A^4+v_B^4)
+\frac{2G^2m_A m_B}{r^2}(m_Bv_A^2+m_Av_B^2)
+\frac{G^3m_A m_B(m_A^2+m_B^2)}{4r^3}.
\]

### 4.5 Unique comparable-mass residual solve

The exact residual solve gives
\[
q_1=-\frac74,
\quad q_2=-\frac14,
\quad q_3=\frac{11}{8},
\quad q_4=\frac14,
\quad q_5=-\frac58,
\quad q_6=\frac32,
\quad q_7=\frac38,
\]
\[
t_1=0,
\quad t_2=\frac{11}{8},
\quad t_3=-\frac{15}{4},
\quad t_4=0,
\quad t_5=0,
\quad t_6=\frac{15}{8},
\]
\[
s_1=\frac54.
\]

Noteworthy zeroes:
\[
t_1=t_4=t_5=0.
\]

### 4.6 Final new 2PN block

The full coefficient of \(c^{-4}\) inserted into the carried Newtonian + 1PN ledger is
\[
\begin{aligned}
L_{2,\rm full}=
&\frac{m_Av_A^6+m_Bv_B^6}{16}
+\frac{Gm_A m_B}{r}
\Big[
\frac78(v_A^4+v_B^4)
-\frac74 v_{AB}(v_A^2+v_B^2)
-\frac14 v_{An}v_{Bn}(v_A^2+v_B^2)
\\
&\hspace{6em}
+\frac{11}{8}v_A^2v_B^2
+\frac14 v_{AB}^2
-\frac58(v_A^2v_{Bn}^2+v_B^2v_{An}^2)
+\frac32 v_{AB}v_{An}v_{Bn}
+\frac38 v_{An}^2v_{Bn}^2
\Big]
\\
&+\frac{G^2m_A m_B}{r^2}
\Big[
\Big(2m_B+\frac{11}{8}m_A\Big)v_A^2
+\Big(2m_A+\frac{11}{8}m_B\Big)v_B^2
-\frac{15}{4}(m_A+m_B)v_{AB}
+\frac{15}{8}(m_Av_{An}^2+m_Bv_{Bn}^2)
\Big]
\\
&+\frac{G^3m_A m_B}{4r^3}(m_A^2+5m_A m_B+m_B^2).
\end{aligned}
\]

### 4.7 Main theorem

After full assembly,
\[
H_2 = H^{\rm ADM}_{2\rm PN},
\qquad
\Delta H_2 \equiv H_2 - H^{\rm ADM}_{2\rm PN} = 0.
\]

### 4.8 Zero-frequency throat-operator data already fixed in the appendices

Odd channels:
\[
R_{1\perp}=\frac72,
\qquad
R_{10}=4.
\]

Even support channels:
\[
R_0 = R_{20} = R_{21} = R_{22} = 1.
\]

Scalar source vector:
\[
J=
\left(
\frac{4}{\sqrt5},\frac54,0,0,0,0
\right).
\]

Pure-potential geometry-closure deficit:
\[
\Delta_{\rm geom}=\frac{281}{80}.
\]

Monopole wall add-on:
\[
\Delta K_{00}=\frac{109}{280}.
\]

### 4.9 Open dynamic observables

The dynamic pole scales that remain genuinely open are
\[
\Omega_{1\perp},
\quad
\Omega_{10},
\quad
\Omega_0,
\quad
\Omega_{20},
\quad
\Omega_{21},
\quad
\Omega_{22},
\quad
\Omega_g.
\]

Optional near-spherical simplification sometimes mentioned:
\[
\Omega_{20}=\Omega_{21}=\Omega_{22},
\]
but this is **not** required by the main proof.

---

## 5) Minimal dictionary and invariant basis used throughout

Two-body generic-frame variables:
\[
\mathbf r = \mathbf x_A - \mathbf x_B,
\qquad
r=|\mathbf r|,
\qquad
\mathbf n = \mathbf r/r.
\]

Velocity invariants:
\[
v_A^2 = \mathbf v_A\cdot\mathbf v_A,
\qquad
v_B^2 = \mathbf v_B\cdot\mathbf v_B,
\qquad
v_{AB}=\mathbf v_A\cdot\mathbf v_B,
\]
\[
v_{An}=\mathbf v_A\cdot\mathbf n,
\qquad
v_{Bn}=\mathbf v_B\cdot\mathbf n.
\]

Local pair potentials:
\[
U_A = \frac{Gm_B}{r},
\qquad
U_B = \frac{Gm_A}{r}.
\]

One-body limit:
\[
U = \frac{GM}{r},
\qquad
u \equiv U/c^2.
\]

PN expansion convention:
\[
L_{\rm cons}=L_0+\frac{1}{c^2}L_1+\frac{1}{c^4}L_2+O(c^{-6}).
\]

Important normalization warning:
- throughout this paper, \(L_2\) means the **coefficient of \(c^{-4}\)**;
- sometimes formulas are displayed with the explicit overall \(c^{-4}\) suppressed.

---

## 6) Main-text proof chain

### 6.1 Carried backbone from the earlier 4D and 1PN papers

The 2PN paper does **not** rederive the lower-order program. It imports:

1. the exact projection map
   \[
   \mathcal P_W[Q](t,\mathbf x)=\int W(w)Q(t,\mathbf x,w)\,dw,
   \qquad \int W(w)dw=1;
   \]

2. exact projected continuity with leakage
   \[
   \partial_t \rho_{\rm br} + \nabla_3\cdot \mathbf J_{\rm br} = S_{\rm leak};
   \]

3. exact longitudinal identity after Helmholtz split
   \[
   \rho_{\rm br}\,\nabla_3^2\phi_{\rm br}
   =
   S_{\rm leak} - \partial_t\rho_{\rm br}
   -(\nabla_3\rho_{\rm br})\cdot(\nabla_3\phi_{\rm br}+\mathbf v_T);
   \]

4. controlled quasi-static Poisson hook
   \[
   \nabla_3^2\phi_{\rm br} \simeq \frac{S_{\rm eff}}{\rho_0};
   \]

5. controlled coherent-defect / worldtube particle reduction;

6. frozen exact 1PN EIH match.

Operationally, the 2PN paper is only allowed to add the new order-\(c^{-4}\) block. It is **not** allowed to retune anything in the carried Newtonian + 1PN ledger.

### 6.2 One-body 2PN closure and test-mass target

The exact isotropic Schwarzschild test-mass expansion is
\[
\frac{L_{\rm Schw}}{m}
=
-c^2 + \frac12 v^2 + U
+\frac{1}{c^2}\left(\frac18 v^4+\frac32Uv^2-\frac12U^2\right)
+\frac{1}{c^4}\left(\frac1{16}v^6+\frac78Uv^4+2U^2v^2+\frac14U^3\right)
+O(c^{-6}).
\]

The carried raw optical denominator gives the right free \(v^6\) and \(Uv^4\) structure, but it leaves the wrong \(U^2v^2\) coefficient. The exact target requires that coefficient to be
\[
2.
\]

The reduced one-body kinetic packaging is written as
\[
L_{\rm red,1body}
=
-mc^2(1-u)
\sqrt{1-\frac{v^2/c^2}{D(u)}},
\qquad
D(u)=1-4u+d_2u^2+O(u^3).
\]

Expanding gives
\[
\frac{L_{\rm red,1body}}{m}
=
-c^2+\frac12v^2+U
+\frac{1}{c^2}\left(\frac18v^4+\frac32Uv^2\right)
+\frac{1}{c^4}\left(\frac1{16}v^6+\frac78Uv^4+\left(6-\frac{d_2}{2}\right)U^2v^2\right)+O(c^{-6}).
\]

Matching the target forces
\[
d_2=8.
\]

The DtN-side invariant packaging uses
\[
Z_{00}(\omega;u)=Z_2(u)\omega^2+Z_4(u)\omega^4+O(\omega^6),
\]
\[
C_s(u)=\frac{Z_4/Z_2^3}{(Z_4/Z_2^3)_0},
\qquad
G(u)=\frac{Z_4/Z_2^2}{(Z_4/Z_2^2)_0}.
\]

On the cylinder / Neumann-bottom unit-test branch,
\[
Z_2=-\frac{L}{c_s^2},
\qquad
Z_4=-\frac{L^3}{3c_s^4},
\]
so
\[
C_s(u)=\frac{c_s^2(u)}{c_{s0}^2},
\qquad
G(u)=\frac{L(u)}{L_0}.
\]

Under the carried freeze,
\[
C_s(u)=1-4u,
\qquad
G(u)=1+\frac{57}{64}u+O(u^2).
\]

The minimal invariant ansatz is
\[
D_{\rm eff}(u)=C_s(u)\left[1+\alpha(G-1)+\beta(G-1)^2\right].
\]

Preserving the frozen 1PN slot forces
\[
\alpha=0,
\]
and then exact isotropic Schwarzschild matching fixes
\[
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right]
=1-4u+8u^2+O(u^3).
\]

This is the paper’s **one-body 2PN closure**.

### 6.3 Minimal self/static 2PN sector

The carried Bernoulli continuation gives
\[
C_s(u)=1-4u.
\]

Minimal Bernoulli exactification of only the optical denominator gives
\[
L_{\rm self,min}
=
-mc^2(1-u)\sqrt{1-\frac{v^2/c^2}{1-4u}}.
\]

Its 2PN expansion is
\[
\frac{L_{\rm self,min}}{m}
=
-c^2+\frac12v^2+U
+\frac{1}{c^2}\left(\frac18v^4+\frac32Uv^2\right)
+\frac{1}{c^4}\left(\frac1{16}v^6+\frac78Uv^4+6U^2v^2\right)
+O(c^{-6}).
\]

Important raw self coefficients:
\[
\frac1{16},\qquad \frac78,\qquad 6.
\]

The raw self sector already gets \(v^6\) and \(Uv^4\) right, but overshoots the target \(U^2v^2\) coefficient by
\[
6-2=4.
\]
That is exactly the slot repaired by the one-body denominator closure above.

Static sector: quadratic local mass scaling is written as
\[
m_A^{\rm eff}=m_A\left(1-\frac{U_A}{c^2}+\lambda_\rho\frac{U_A^2}{c^4}\right),
\qquad
m_B^{\rm eff}=m_B\left(1-\frac{U_B}{c^2}+\lambda_\rho\frac{U_B^2}{c^4}\right).
\]

This yields
\[
L_{\rm static}
=
\frac{Gm_A m_B}{r}
-\frac{G^2m_A m_B(m_A+m_B)}{2c^2r^2}
+\frac{\lambda_\rho G^3m_A m_B(m_A^2+m_B^2)}{2c^4r^3}
+O(c^{-6}).
\]

The exact isotropic one-body target requires the \(U^3\) coefficient to be \(+1/4\), so
\[
\lambda_\rho = \frac12.
\]

Then the one-body self/static seed becomes exactly the Schwarzschild target,
\[
\frac{L_{\mathrm{self+static},1body}}{m}
=
-c^2+\frac12v^2+U
+\frac{1}{c^2}\left(\frac18v^4+\frac32Uv^2-\frac12U^2\right)
+\frac{1}{c^4}\left(\frac1{16}v^6+\frac78Uv^4+2U^2v^2+\frac14U^3\right)
+O(c^{-6}).
\]

Symmetric two-body promotion then gives the frozen self/static seed
\[
L_{2,\rm self+static}
=
\frac{m_Av_A^6+m_Bv_B^6}{16}
+\frac{7Gm_A m_B}{8r}(v_A^4+v_B^4)
+\frac{2G^2m_A m_B}{r^2}(m_Bv_A^2+m_Av_B^2)
+\frac{G^3m_A m_B(m_A^2+m_B^2)}{4r^3}.
\]

### 6.4 Comparable-mass ADM lift

The key algebraic simplification is the perturbative Legendre transform:
\[
L=L_0+\epsilon L_1+\epsilon^2L_2,
\qquad \epsilon=c^{-2},
\]
with quadratic \(L_0\), then
\[
H_1=-L_1(v_0),
\qquad
H_2=-L_2(v_0)+\frac12 A_0^TM^{-1}A_0,
\]
where
\[
v_0=p/m,
\qquad
A_0=\left.\frac{\partial L_1}{\partial v}\right|_{v=v_0}.
\]

Because the full lower-order block is already frozen, any genuinely new 2PN Lagrangian term enters the Hamiltonian only through the explicit
\[
-L_{2,\rm new}(v_0)
\]
piece.

So once the one-body/self/static seed is fixed, the remaining comparable-mass mismatch against the generic-frame ADM target is just a finite algebraic residual.

The compact invariant basis is organized into three groups.

#### Quartic \(G/r\) group
\[
\frac{Gm_A m_B}{r}\times
\Big\{
 v_{AB}(v_A^2+v_B^2),
 v_{An}v_{Bn}(v_A^2+v_B^2),
 v_A^2v_B^2,
 v_{AB}^2,
 v_A^2v_{Bn}^2+v_B^2v_{An}^2,
 v_{AB}v_{An}v_{Bn},
 v_{An}^2v_{Bn}^2
\Big\}.
\]

#### Quadratic \(G^2/r^2\) group
\[
\frac{G^2m_A m_B}{r^2}\times
\Big\{
 m_Bv_A^2+m_Av_B^2,
 m_Av_A^2+m_Bv_B^2,
 (m_A+m_B)v_{AB},
 (m_A+m_B)v_{An}v_{Bn},
 m_Bv_{An}^2+m_Av_{Bn}^2,
 m_Av_{An}^2+m_Bv_{Bn}^2
\Big\}.
\]

#### Static cross mass polynomial
\[
\frac{G^3}{r^3}m_A^2m_B^2.
\]

The solve is unique and yields the coefficients recorded above:
\[
\{q_i,t_i,s_1\}.
\]

### 6.5 Full two-body 2PN assembly and exact ADM match

Once the lower-order ledger, one-body denominator closure, self/static seed, and comparable-mass residual solution are inserted, the final reduced two-body Lagrangian is exactly the standard generic-frame conservative ADM theory through 2PN.

This is the main theorem of the paper.

### 6.6 Fixed versus symbolic quantities

The paper makes a very sharp distinction between:

1. quantities already fixed by the derivation chain,
2. quantities fixed only within a declared response protocol,
3. appendix-side zero-frequency throat data fixed by the solved target,
4. and quantities that must remain symbolic because they are genuine dynamic throat observables or higher-order/finite-size effects.

This distinction is one of the main deliverables of the paper.

### 6.7 Verification and reproducibility

The paper’s verification archive has three layers:

1. inherited 1PN WL referee suite,
2. compact 2PN WL harness,
3. symbolic 2PN derivation ledger (`4d_2pn_full_notes.md`).

The paper is explicit that the archive is not a black-box derivation engine. It is a referee / reproducibility ledger for the analytic derivation.

---

## 7) What is already fixed vs what remains open

### 7.1 Fixed by the paper / carried as frozen data

These should be treated as fixed in future work.

#### Lower-order carry-forward data
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
\[
C_\parallel=-\frac72,
\qquad
C_L=-\frac12.
\]

#### 2PN one-body closure
\[
\alpha=0,
\qquad
\mu=\frac{32768}{3249},
\qquad
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G-1)^2\right].
\]

#### Static 2PN gate
\[
\lambda_\rho=\frac12.
\]

#### Comparable-mass cross coefficients
\[
q_i,
\quad t_i,
\quad s_1
\]
with the exact values already listed.

#### Full conservative two-body 2PN assembly
\[
H_2=H_{2\rm PN}^{\rm ADM}.
\]

#### Zero-frequency throat operator data
\[
R_{1\perp}=\frac72,
\qquad R_{10}=4,
\qquad R_0=R_{20}=R_{21}=R_{22}=1,
\]
\[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
\qquad
\Delta_{\rm geom}=\frac{281}{80},
\qquad
\Delta K_{00}=\frac{109}{280}.
\]

### 7.2 Fixed only within the declared response hierarchy

These are fixed **within the paper’s chosen conservative hierarchy**, but not yet as assumption-free PDE theorems.

- the specific invariant denominator ansatz;
- the carried adiabatic breathing closure;
- the breathing slope
  \[
  g_1=\frac{57}{64};
  \]
- the Family-1 wall-profile branch used for the final coupled monopole breathing normalization.

### 7.3 Still genuinely open

These remain open after the paper.

1. The dynamic pole scales
   \[
   \Omega_{1\perp},\Omega_{10},\Omega_0,\Omega_{20},\Omega_{21},\Omega_{22},\Omega_g.
   \]
2. A fully dynamical first-principles moving-throat PDE derivation of those pole data.
3. A more systematic finite-size and higher-multipole worldtube extension.
4. Dissipative, radiative, and higher-PN sectors.
5. Spin couplings and other particle-identity questions beyond the conservative scalar/tensor response problem.

---

## 8) Appendix digest — what each appendix adds

The appendices are not decorative. They contain the constructive throat-response program that sits underneath the main theorem.

### 8.1 Notation and conventions appendix

This appendix freezes:
- coordinates,
- pair variables,
- PN bookkeeping,
- claim-class labels.

Operationally important conventions:
- generic frame is used throughout;
- no center-of-mass specialization in the residual solve;
- \(L_2\) means the coefficient of \(c^{-4}\);
- appendix-side operators are often presented at zero frequency first and then dynamically completed.

### 8.2 Carry-forward formulas from the earlier 4D and 1PN papers

This appendix collects the exact / carried formulas used without rederivation.

#### Exact projection formulas
\[
\mathcal P_W[Q](t,\mathbf x)=\int W(w)Q(t,\mathbf x,w)\,dw,
\qquad
\rho_{\rm br}=\mathcal P_W[\rho],
\qquad
\mathbf J_{\rm br}=\int W(w)\mathbf j_{xyz}\,dw.
\]

Projected continuity with leakage:
\[
\partial_t\rho_{\rm br}+\nabla_3\cdot\mathbf J_{\rm br}=S_{\rm leak}.
\]

Longitudinal identity:
\[
\rho_{\rm br}\,\nabla_3^2\phi_{\rm br}
=
S_{\rm leak}-\partial_t\rho_{\rm br}
-(\nabla_3\rho_{\rm br})\cdot(\nabla_3\phi_{\rm br}+\mathbf v_T).
\]

#### Controlled Poisson hook / Newtonian point-particle limit
\[
\nabla_3^2\phi_{\rm br}\simeq \frac{S_{\rm eff}}{\rho_0},
\qquad
\Phi_N = \lambda_\Phi\,\phi_{\rm br}.
\]

Worldtube point-particle reduction:
\[
M\ddot{\mathbf X}_{\rm cm}
=
-M\nabla_3\Phi_N + O(Q\nabla^3\Phi_N),
\]
with quadrupole suppression
\[
\frac{|\mathbf F_{\rm quad}|}{|\mathbf F_{\rm mono}|}=O(a^2/r^2).
\]

#### Scalar/optical 1PN worldline formulas
\[
L_{\rm sc}
=
-M_0\left(1+\kappa_\rho\frac{\Phi_N}{c^2}\right)c^2
\sqrt{1-\frac{v^2/c^2}{1+(n-1)\Phi_N/c^2}}.
\]

Frozen outputs:
\[
\kappa_\rho=1,
\qquad
n=5,
\qquad
\text{quartic free coefficient }=\frac18,
\qquad
\text{self coefficient }=+\frac32.
\]

#### Static and response ledgers at 1PN
\[
L_{\rm stat}^{(AB)}=-\frac{G^2m_A m_B(m_A+m_B)}{2c^2r_{AB}^2}.
\]

Self ledger split:
\[
\beta_{\rm 1PN}=\kappa_\rho+\kappa_{\rm add}+\kappa_{\rm PV}.
\]

Frozen values:
\[
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3.
\]

Breathing closure data carried forward:
\[
E_w:E_f:E_{\rm PV}=11:2:5,
\qquad
\frac{d\ln a_*}{d\ln\rho}=-\frac{57}{64}.
\]

#### Frozen wake formulas
The isotropic wake reconstruction carries forward
\[
\alpha^2=\frac34,
\qquad
a_H=0,
\qquad
a_L=\pm\frac{\sqrt3}{2},
\qquad
K_{\rm vec}=\frac{2}{\pi^2},
\]
which imply
\[
C_\parallel=-\frac72,
\qquad
C_L=-\frac12.
\]

This appendix is basically the paper’s “do not reopen this” ledger.

### 8.3 One-body DtN closure details appendix

This appendix gives the detailed low-frequency one-body algebra.

#### Cylinder / Neumann-bottom unit-test branch
\[
Z_2=-\frac{L}{c_s^2},
\qquad
Z_4=-\frac{L^3}{3c_s^4}.
\]
Hence
\[
C_s(u)=\frac{c_s^2(u)}{c_{s0}^2},
\qquad
G(u)=\frac{L(u)}{L_0}.
\]

Current conservative freeze:
\[
\frac{c_s^2}{c^2}=1-4u,
\qquad
a(u)=\frac{L(u)}{L_0}=1+\frac{57}{64}u+O(u^2).
\]

The appendix records the fuller series
\[
a(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3),
\]
so
\[
\frac{Z_2(u)}{Z_2(0)}=1+\frac{313}{64}u+\frac{2862917}{131072}u^2+O(u^3),
\]
\[
\frac{Z_4(u)}{Z_4(0)}=1+\frac{683}{64}u+\frac{10301487}{131072}u^2+O(u^3).
\]

Then
\[
C_s(u)=1-4u,
\qquad
G(u)=1+\frac{57}{64}u+\frac{298821}{131072}u^2+O(u^3).
\]

#### Minimal denominator ansatz and uniqueness
The general quadratic invariant ansatz is
\[
D_{\rm eff}(u)=C_s(u)\left[1+\alpha(G-1)+\beta(G-1)^2\right].
\]
Preserving the frozen 1PN coefficient forces
\[
\alpha=0.
\]
Then matching the exact isotropic target fixes
\[
\mu g_1^2=8,
\qquad
g_1=\frac{57}{64},
\qquad
\mu=\frac{32768}{3249} \approx 10.0855647892.
\]

So the final one-body denominator is
\[
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right].
\]

#### Relation to earlier raw resonance-proxy fit
The raw resonance proxy is
\[
D_{\rm raw}(u)
=
\frac{[Z_2/Z_4](u)}{[Z_2/Z_4](0)}
=
\frac{1-4u}{G(u)^2}
=1-\frac{185}{32}u+\frac{324075}{65536}u^2+O(u^3).
\]

The multiplicative port factor needed to convert this into the exact target is
\[
P_{\rm port}(u)=G(u)^2\left[1+\mu(G(u)-1)^2\right]
=1+\frac{57}{32}u+\frac{875093}{65536}u^2+O(u^3).
\]

So the old port fit factorizes cleanly as
\[
D_{\rm eff}(u)=D_{\rm raw}(u)\,P_{\rm port}(u).
\]

This is one of the paper’s strongest “old intuition was not arbitrary” results.

### 8.4 Bernoulli exactification, quadratic local mass scaling, and the static hierarchy appendix

This appendix spells out the self/static formulas in full detail.

#### Exact \(n=5\) continuation of the barotropic sector
\[
P(\rho)=K_{\rm EOS}\rho^5,
\qquad
U_{\rm int}(\rho)=\frac{K_{\rm EOS}}{4}\rho^5,
\qquad
h(\rho)=\frac{5K_{\rm EOS}}{4}\rho^4.
\]

Background normalization is fixed by
\[
c_s^2(\rho_0)=c^2.
\]

Exact Bernoulli continuation then gives
\[
h(\rho)=h_0+\Phi,
\qquad
\rho_{\rm Bern}(\Phi)=\rho_0\left(1+\frac{4\Phi}{c^2}\right)^{1/4},
\qquad
\frac{c_s^2(\Phi)}{c^2}=1+\frac{4\Phi}{c^2}.
\]

With \(\Phi=-U\), this becomes
\[
\frac{c_s^2}{c^2}=1-4u.
\]

The optical index is
\[
N(U)=\left(1-\frac{4U}{c^2}\right)^{-1/2}
=1+2\frac{U}{c^2}+6\frac{U^2}{c^4}+O(c^{-6}).
\]

This is the source of the raw self coefficients \(1/16,7/8,6\).

#### Minimal vs aggressive self exactification
Minimal exactification:
\[
L_{\rm self,min}
=-mc^2\left(1+\frac{\Phi}{c^2}\right)
\sqrt{1-\frac{v^2}{c^2(1+4\Phi/c^2)}}.
\]

Aggressive exactification would Bernoulli-dress the mass prefactor too, but that generates a forbidden \(+3/2\,U^2/c^2\) static term already at 1PN. So the paper keeps only the optical denominator exactified.

#### Static hierarchy
Quadratic local mass scaling yields not only the two-body \(G^3/r^3\) term but a compact higher-body cubic-static hierarchy.

The appendix records selected 3-body and 4-body terms explicitly. These are not needed for the main two-body theorem, but they are the first concrete higher-body predictions of the same local-scaling rule.

#### Worldtube finite-size gate
The appendix gives the sharp regime statement:
\[
\frac{F^{(2)}_{\rm finite-size}}{F^{(2)}_{\rm universal}}
=\left(\frac{ac^2}{GM}\right)^2.
\]

So the universal conservative 2PN force is valid only when the body remains deep in the controlled small-worldtube regime.

### 8.5 Perturbative Legendre transform and ADM residual-solve details appendix

This appendix is the fully explicit algebraic heart of the ADM lift.

Main ingredients:
- perturbative Legendre identity,
- explicit residual basis,
- exact solve for \(q_i,t_i,s_1\),
- exact zero residual against the generic-frame ADM target.

The main lesson is that once the lower-order ledger and the self/static seed are frozen, the remaining comparable-mass problem is not mysterious; it is a finite linear residual solve.

This appendix is the symbolic basis for the statements in `4d_2pn_full_notes.md`.

### 8.6 Tensor-wake decomposition appendix

This appendix gives the first constructive rewrite of the solved added comparable-mass cross block.

#### Universal kinetic dressing of the frozen 1PN wake
The carried 1PN vector-wake kernel is
\[
\mathcal W^{(1)}_{AB}=C_\parallel v_{AB}+C_Lv_{An}v_{Bn},
\qquad C_\parallel=-\frac72,
\quad C_L=-\frac12.
\]

The first two quartic 2PN coefficients are exactly
\[
-\frac74 v_{AB}(v_A^2+v_B^2)-\frac14 v_{An}v_{Bn}(v_A^2+v_B^2)
=
\frac12(v_A^2+v_B^2)\,\mathcal W^{(1)}_{AB}.
\]

So those two coefficients are **not new fits**; they are a universal kinetic dressing of the frozen 1PN wake with
\[
\sigma=\frac12.
\]

#### Minimal rank-2 tensor basis for the quartic residual
Define the pair-projector invariants:
\[
T_A=v_A^2-v_{An}^2,
\qquad
L_A=v_{An}^2,
\qquad
T_B=v_B^2-v_{Bn}^2,
\qquad
L_B=v_{Bn}^2,
\]
\[
S_{AB}=(v_{AB}-v_{An}v_{Bn})^2-\frac12T_AT_B,
\]
\[
M_{AB}=2(v_{AB}-v_{An}v_{Bn})v_{An}v_{Bn}.
\]

Then the residual quartic sector is exactly
\[
\frac{Gm_A m_B}{r}
\Big[
\frac32T_AT_B
+\frac14S_{AB}
+1\cdot M_{AB}
+\frac34(T_AL_B+L_AT_B)
+\frac94L_AL_B
\Big].
\]

So the minimal five-channel tensor kernel has coefficients
\[
k_{TT}=\frac32,
\qquad
k_S=\frac14,
\qquad
k_M=1,
\qquad
k_{TL}=\frac34,
\qquad
k_{LL}=\frac94.
\]

The scalar \((T,L)\) subsector matrix is
\[
K_{\rm scalar}=
\begin{pmatrix}
\frac32 & \frac34\\[3pt]
\frac34 & \frac94
\end{pmatrix},
\]
which is positive:
\[
\det K_{\rm scalar}=\frac{45}{16}>0,
\qquad
\operatorname{tr}K_{\rm scalar}=\frac{15}{4}>0,
\]
with eigenvalues
\[
\lambda_\pm = \frac{15\pm 3\sqrt5}{8}.
\]

#### Local-potential dressing of the \(G^2/r^2\) block
Using
\[
U_A=\frac{Gm_B}{r},
\qquad U_B=\frac{Gm_A}{r},
\]
the quadratic-velocity block is reproduced by:
- purely parallel local-potential dressing of the frozen 1PN wake;
- plus diagonal transverse/longitudinal tensor-potential channels.

Exact coefficients:
\[
\tau_\parallel=\frac{15}{14},
\qquad
\beta_T=\frac{11}{8},
\qquad
\beta_L=\frac{13}{4}.
\]

#### Full constructive decomposition of the added cross block
The appendix writes
\[
L_{2,\rm cross}
=
L_{2,\rm vec,kin}
+
L_{2,\rm tensor,quartic}
+
L_{2,\rm quad,pot}
+
L_{2,\rm static,cross},
\]
where the four pieces are exactly reconstructed from the coefficients above.

#### First 3-body predictions
The same local-potential module already predicts nontrivial pair-\(AB\) / background-\(C\) three-body coefficients:
- \(-15/4\) on the pair-velocity cross terms,
- \(11/8\) on the transverse quadratic terms,
- \(15/8\) on the longitudinal quadratic terms.

So the constructive module already makes nontrivial many-body predictions beyond the matched two-body sector.

### 8.7 \(P_0\oplus P_2\) mouth-port rebuild and support-minus-closure EFT appendix

This appendix sharpens the tensor-wake decomposition into a body-local mouth-port theory and then into a signed auxiliary-channel EFT.

#### Carried odd dipole ports from 1PN
In the pair frame \(\mathbf n=\hat z\), write
\[
\mathbf v = \mathbf u + d\,\mathbf n,
\qquad \mathbf u\cdot\mathbf n = 0.
\]

The frozen 1PN wake is already an odd dipole kernel:
\[
L^{\rm wake}_{1\rm PN}
=-\Big[\Pi_A^{(1\pm)}\cdot\Pi_B^{(1\pm)} + \Pi_A^{(10)}\Pi_B^{(10)}\Big],
\]
with odd ports
\[
\Pi^{(1\pm)} = \sqrt{\frac72}\,\mathbf u,
\qquad
\Pi^{(10)} = 2d.
\]

#### Scalar diagonalization into \(P_0\) and \(P_{20}\)
The scalar variables are
\[
T=u^2,
\qquad L=d^2,
\]
with kernel
\[
K_{TL}=
\begin{pmatrix}
\frac32 & \frac34\\[3pt]
\frac34 & \frac94
\end{pmatrix}.
\]

Define
\[
P_0=T+L=v^2,
\qquad
P_{20}=-T+2L=3(v\cdot n)^2-v^2.
\]

Then the kernel diagonalizes exactly to
\[
K_{P_0P_{20}}=
\begin{pmatrix}
\frac54 & 0\\[3pt]
0 & \frac14
\end{pmatrix}.
\]

Canonical scalar ports:
\[
\Pi^{(0)}=\frac{\sqrt5}{2}v^2,
\qquad
\Pi^{(20)}=\frac12\big(3(v\cdot n)^2-v^2\big).
\]

Real quadrupole partners:
\[
\Pi^{(21c)}=\sqrt2(v\cdot n)u_x,
\qquad
\Pi^{(21s)}=\sqrt2(v\cdot n)u_y,
\]
\[
\Pi^{(22c)}=\frac{u_x^2-u_y^2}{2\sqrt2},
\qquad
\Pi^{(22s)}=\frac{2u_xu_y}{2\sqrt2}.
\]

Then the genuinely new quartic residual is exactly a positive Gram overlap:
\[
Q^{\rm new}_{2\rm PN}
=
\Pi_A^{(0)}\Pi_B^{(0)}
+\Pi_A^{(20)}\Pi_B^{(20)}
+\Pi_A^{(21c)}\Pi_B^{(21c)}
+\Pi_A^{(21s)}\Pi_B^{(21s)}
+\Pi_A^{(22c)}\Pi_B^{(22c)}
+\Pi_A^{(22s)}\Pi_B^{(22s)}.
\]

So the new even 2PN support sector is exactly **one monopole-like scalar channel plus the full real quadrupole multiplet**.

#### Source shifts and support-minus-closure factorization
Using
\[
U_A=\frac{Gm_B}{r},
\qquad U_B=\frac{Gm_A}{r},
\]
the solved \(U\times \Pi\) block excites only the scalar axisymmetric ports:
\[
-\frac{15}{4}(U_A+U_B)(\mathbf v_A\cdot\mathbf v_B)
+
U_A\left(\frac{4}{\sqrt5}\Pi_B^{(0)}+\frac54\Pi_B^{(20)}\right)
+
U_B\left(\frac{4}{\sqrt5}\Pi_A^{(0)}+\frac54\Pi_A^{(20)}\right).
\]

Define shifted support bundles
\[
\mathcal S_A=
\Big(
\Pi_A^{(0)}+\frac{4}{\sqrt5}U_A,
\Pi_A^{(20)}+\frac54U_A,
\Pi_A^{(21c)},
\Pi_A^{(21s)},
\Pi_A^{(22c)},
\Pi_A^{(22s)}
\Big),
\]
and similarly for \(\mathcal S_B\).

Support alone would generate a static coefficient
\[
\left(\frac{4}{\sqrt5}\right)^2 + \left(\frac54\right)^2 = \frac{381}{80},
\]
while the solved target static coefficient is only
\[
\frac54 = \frac{100}{80}.
\]

So the exact closure deficit is
\[
\Delta_{\rm geom} = \frac{381}{80}-\frac{100}{80} = \frac{281}{80}.
\]

The negative closure channel is purely potential:
\[
\mathcal C_A = \sqrt{\frac{281}{80}}U_A,
\qquad
\mathcal C_B = \sqrt{\frac{281}{80}}U_B.
\]

Hence the full even block factorizes exactly as
\[
L^{\rm even}_{2\rm PN} = \mathcal S_A\cdot\mathcal S_B - \mathcal C_A\mathcal C_B.
\]

This is the appendix where the added 2PN cross sector first becomes a genuine **support-minus-closure throat-response theory**.

#### Signed auxiliary-channel EFT interpretation
Positive support channel auxiliary:
\[
L_q = -\frac12 q^2 + q(S_A+S_B) \quad \Rightarrow \quad +S_AS_B.
\]

Negative closure channel auxiliary:
\[
L_g = +\frac12 g^2 + g(C_A+C_B) \quad \Rightarrow \quad -C_AC_B.
\]

Odd dipole auxiliary:
\[
L_d = -\frac12 d^2 + d(D_A-D_B) \quad \Rightarrow \quad -D_AD_B.
\]

So the solved 2PN sector is exactly the on-shell result of a small signed auxiliary-channel EFT.

### 8.8 Minimal irreducible-channel throat operator and dynamic DtN scaffolding appendix

This appendix turns the constructive decompositions into the smallest low-frequency throat operator that reproduces the solved added conservative cross target.

#### Minimal zero-frequency channel content
Odd channels:
\[
1\perp:\ \mathbf u,
\qquad
10:\ d.
\]

Even channels:
\[
0,\ 20,\ 21c,\ 21s,\ 22c,\ 22s.
\]

Separate geometry-closure channel:
\[
g.
\]

So the minimal irreducible static channel set is
\[
\{1\perp,10,0,20,21c,21s,22c,22s,g\}.
\]

Odd sector dressings extracted from the solved operator:
\[
\eta_\perp=\frac{15}{14},
\qquad
\eta_\parallel=\frac{15}{16}.
\]

The full added cross operator is summarized as
\[
L^{\rm add}_{2\rm PN,cross}
=
\frac{Gm_A m_B}{r}
\left[
L^{\rm add}_{\rm odd}
+
(\Pi_A+JU_A)\cdot(\Pi_B+JU_B)
-
\Delta_{\rm geom}U_AU_B
\right].
\]

#### Immediate 3-body lift
Promoting the body-local pair potentials to local 3-body ones immediately induces background-\(C\) lifts. The appendix writes the explicit velocity and static 3-body terms. The main message is:

> once the zero-frequency operator is known, it already carries nontrivial higher-body content.

#### Isotropic 4D-ball DtN unit test
The regular interior Helmholtz problem in a 4D ball gives diagonal harmonic channels with positive \(\ell=1,2\) towers, but it fails the solved static operator in exactly the ways expected.

So the isotropic cavity is a **useful no-go benchmark**, not the final answer.

#### Minimal axisymmetric pole-completed DtN scaffolding
The appendix then builds a one-pole-per-channel dynamic continuation. The support-minus-closure matrix is promoted to frequency-dependent admittances \(Y_a(\omega)\), and the scalar source vector becomes
\[
J_{\rm eff}(\omega)=\big(J_0Y_0(\omega),J_{20}Y_{20}(\omega),0,0,0,0\big).
\]

The effective pure-potential coefficient becomes
\[
S_{\rm eff}(\omega)=J_0^2Y_0(\omega)+J_{20}^2Y_{20}(\omega)-\Delta_{\rm geom}Y_g(\omega),
\]
with
\[
S_{\rm eff}(0)=\frac54.
\]

Its low-frequency expansion is
\[
S_{\rm eff}(\omega)
=
\frac54
+\omega^2\left(
\frac{16}{5\Omega_0^2}
+\frac{25}{16\Omega_{20}^2}
-\frac{281}{80\Omega_g^2}
\right)
+O(\omega^4).
\]

So the appendix isolates the first nontrivial dynamic scalar data very sharply.

This appendix ends with the key status statement:

> the entire zero-frequency throat operator is fixed; what remains open is the dynamic pole structure.

### 8.9 Family-1 soft-wall and wall-stress program appendix

This appendix gives the constructive Family-1 chain that turns the solved even support/source sector into an explicit local wall law on a flared soft wall.

#### Solved support/source sector already sits in the Family-1 flare basis
The full passed static support data can be written in the same low-order flare basis
\[
\{1,\mu^2,\mu^4\},
\qquad \mu=\cos\theta.
\]

The support operator is packaged via
\[
Z_{\ell m}
=
\langle z_{\rm base}(\mu)\rangle_{\ell m}
+
(\lambda_\ell-2)\langle z_{\rm curv}(\mu)\rangle_{\ell m},
\qquad
\lambda_\ell=\ell(\ell+1),
\]
with
\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
z_{\rm curv}(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

The scalar source sector collapses into the same basis:
\[
S(\mu)=\frac{11}{8}+\frac{15}{8}P_2(\mu)=\frac{7}{16}+\frac{45}{16}\mu^2.
\]

#### Exact reduced wall Hessian
The static support data admit an exact reduced wall-Hessian on the truncated \(\ell=0,1,2\) basis. Solving the reduced form gives
\[
K_{\rm mono}=\frac{4}{45},
\qquad
T_0=\frac{23}{135},
\]
\[
A_0=\frac{9095}{3528},
\qquad
A_1=-\frac{25559}{21168},
\qquad
B_0=-\frac{109}{56},
\qquad
B_1=\frac{1765}{3024}.
\]

So the passed static support data are not just a Robin table; they are the Hessian of a concrete reduced wall law.

A sharp sub-result: the minimal linear-gradient interpretation points to a flare width of order the throat radius,
\[
\frac{z_m}{a_0}\approx 1.0113880097.
\]

#### Strict isotropic steep-wall no-go
A strict isotropic steep-wall pullback of the flare geometry fails. The best-fit residual is not small:
\[
\sum_i \Delta_i^2 \approx 0.5536733194,
\qquad
\max_i|\Delta_i|\approx 0.4390744081.
\]

So isotropic penetration plus geometry pullback is too small a model.

#### Exact traction-balance completion
The smallest exact repair is one base pressure profile plus one tangential wall-stress / curvature-compliance profile.

The exact coefficients are
\[
b_0=\frac{17}{56},
\quad
b_2=-\frac{5}{56},
\quad
t_0=\frac{593}{672},
\quad
t_2=-\frac{1553}{672},
\quad
t_4=\frac78,
\]
so
\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4.
\]

The induced pressure profile is
\[
p(\mu)=\frac{571}{672}-\frac{5141}{672}\mu^2+7\mu^4.
\]

This local wall model kills all five \(\ell=1,2\) support residuals exactly.

#### Minimal anisotropic PDE-side completion and universal monopole identity
Promoting only the tangential wall moment to the smallest axisymmetric profile
\[
B_{\rm tan}(\mu)=B_0+B_2\mu^2
\]
repairs the strict isotropic no-go.

On the natural Family-1-like branch,
\[
A\approx -0.2819232193,
\qquad
B_0\approx 0.6489147092,
\qquad
B_2\approx -1.0854346039,
\]
\[
q\approx 2.3709157172,
\qquad
r\approx 2.7583434741.
\]

The structural payoff is the universal monopole identity
\[
K_{00}=K_{1\perp}+\frac12K_{10}-\frac1{10}K_{20}-\frac15K_{21}-\frac15K_{22}.
\]

Using the solved support targets gives the raw boundary-layer monopole
\[
K_{00}^{\rm BL}=-\frac{757}{2520},
\]
so the separate monopole add-on required to hit the target \(K_{00}=4/45\) is
\[
\Delta K_{00}^{\rm mono} = \frac{109}{280}.
\]

So the support-sector repair does **not** eliminate the separate monopole channel; it predicts its exact size.

#### Support/source locking and canonical scalar source vector
The exact source profile obeys
\[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu).
\]
So once the base support profile is fixed, the source profile is fixed automatically.

Using canonical normalization, the scalar source vector becomes
\[
J_0=\frac{2}{\sqrt5}\left(p_{\rm iso}+\frac{q_{\rm ax}}{3}\right),
\qquad
J_{20}=\frac23 q_{\rm ax},
\]
which at the passed dipole support values gives exactly
\[
J=\left(\frac{4}{\sqrt5},\frac54\right).
\]

This is the same canonical scalar source vector that appears in the mouth-port and irreducible-channel appendices.

The full local wall operator is therefore
\[
\mathcal O_{\rm loc}
=
z_{\rm base}(\mu)+\frac12\{-\Delta_S-2,\,t(\mu)\},
\]
with channel values
\[
K_{00}^{\rm raw}=-\frac{757}{2520},
\qquad
K_{1\perp}=\frac27,
\qquad
K_{10}=\frac14,
\qquad
K_{20}=\frac49,
\qquad
K_{21}=\frac23,
\qquad
K_{22}=\frac83.
\]

Then the full static even-wall completion is simply
\[
\mathcal O_{\rm full}=\mathcal O_{\rm loc}+\frac{109}{280}\,\mathbb P_{00}.
\]

This appendix is where the support/source data stop being just a solved operator list and become a concrete Family-1 wall program.

### 8.10 Geometry-breathing and dynamic monopole closures appendix

This appendix isolates and dynamizes the scalar monopole channel left over after the local wall law closes the \(\ell=1,2\) support/source sector.

Important warning repeated from the paper:
- \(\Delta_{\rm geom}=281/80\) is the pure-potential closure deficit in the even cross EFT;
- \(\Delta K_{00}=109/280\) is the separate monopole wall add-on.

This appendix is about the second object.

#### Static geometry-Hessian closure in \((a,L)\)
Use the cylinder-like throat geometry
\[
V(a,L)=\frac{4\pi}{3}a^3L,
\qquad
A(a,L)=4\pi a^2L+\frac{8\pi}{3}a^3,
\]
and the minimal curvature-completed geometry energy
\[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L)+\kappa_b\frac{a^2}{L}.
\]

The uniform monopole port couples through volume work, and integrating out the reduced geometry variables gives the exact static compressibility
\[
\Delta K_{00}^{\rm geom}=\frac{g^TH_0^{-1}g}{V_0^2}.
\]

This reinterprets the old “monopole auxiliary” as a genuine reduced geometry compressibility.

Without the curvature completion (\(\beta=0\)), the baseline Hessian has
\[
\det H_{\rm base}<0,
\]
so literal \(P_{\rm vac}V+\sigma A\) geometry is not a passive 2-DOF Hessian and cannot by itself close the monopole channel.

Representative EM worked point:
\[
\Lambda_{\rm EM}=\frac{\sqrt2\pi}{x_{01}}\approx 1.8474865771,
\qquad
\rho=\frac{1}{10},
\qquad
\beta=12,
\qquad
\Sigma_*\approx 0.2076143291835488854.
\]

#### Exact dynamic monopole susceptibility — two-pole breathing branch
Dimensionless geometry variables:
\[
\mathbf Q_{\rm geom}=\left(\frac{\delta a}{a_0},\frac{\delta L}{a_0}\right).
\]

Affine inertia metric (uniform-slice version):
\[
\widehat M=
\begin{pmatrix}
\frac35 & 0\\[3pt]
0 & \frac1{12}
\end{pmatrix}.
\]

Scaled frequency variable:
\[
s=\omega^2\frac{\rho_{\rm eff}V_0a_0^2}{\Sigma}.
\]

Exact reduced geometry response:
\[
Y_{\rm geom}(s)=\frac{1}{\Sigma}\,\bar g^T(\widehat H-s\widehat M)^{-1}\bar g.
\]

Raw monopole wall channel:
\[
K_{00}^{\rm raw}(s)=-\frac{757}{2520}+Y_{\rm geom}(s).
\]

Because \(\widehat H\) and \(\widehat M\) are positive on the worked branch, this is a two-pole Stieltjes function:
\[
Y_{\rm geom}(s)=\frac{R_-}{1-s/\lambda_-}+\frac{R_+}{1-s/\lambda_+},
\qquad
R_\pm>0,
\qquad
R_-+R_+=\frac{109}{280}.
\]

At the EM affine worked point,
\[
\lambda_-\approx 5.23115613,
\qquad
\lambda_+\approx 230.99360624,
\]
\[
R_-\approx 0.00327376153,
\qquad
R_+\approx 0.38601195275.
\]

Physical pole frequencies:
\[
\Omega_\pm^2 = \frac{\Sigma_*}{\rho_{\rm eff}V_0a_0^2}\,\lambda_\pm.
\]

#### Controlled one-pole reduction / historical breathing auxiliary
The low-frequency Padé reduction is
\[
Y_{\rm geom}(s)\approx \frac{\Delta_0}{1-s/\lambda_{\rm eff}},
\qquad
\lambda_{\rm eff}=\frac{\Delta_0}{\Delta_2},
\qquad
\Delta_0=\frac{109}{280}.
\]

At the EM affine worked point,
\[
\lambda_{\rm eff}\approx 169.48205088,
\]
and on the natural low-frequency band
\[
0\le s\le 0.1\lambda_-,
\]
the maximum relative error is only
\[
\max|\delta Y/Y|\approx 8.87\times 10^{-5}.
\]

So the historical single breathing auxiliary is a controlled low-frequency reduction of the exact two-pole geometry response.

#### Family-1 bulk inertia and the \(n=5\) Thomas–Fermi branch
The appendix then derives the effective bulk inertia from the Family-1 parent PDE on the filled-to-endcap Thomas–Fermi branch.

Key outputs:
\[
\rho_{\rm eff}^{\rm TF}(n)=c_0(n)\rho_c,
\]
with exact closed form
\[
c_0(n)=\frac{\sqrt\pi\,\Gamma\left(1+\frac{1}{n-1}\right)}{2\,\Gamma\left(\frac32+\frac{1}{n-1}\right)}.
\]

For the frozen \(n=5\) branch,
\[
\rho_{\rm eff}^{\rm TF}(5)=
\frac{\sqrt\pi\,\Gamma(1/4)}{6\Gamma(3/4)}
\left(\frac{m_\psi\Omega_{\rm in}^2L_0^2}{10K_{\rm EOS}}\right)^{1/4}.
\]

The EOS-sensitive axial inertia moment becomes
\[
\widehat M_{\rm TF}(5)=
\begin{pmatrix}
\frac35 & 0\\[3pt]
0 & \frac1{14}
\end{pmatrix}.
\]

With this replacement, the exact reduced response is
\[
Y_{\rm geom}^{\rm TF}(s)=\frac{1}{\Sigma}\,\bar g^T(\widehat H-s\widehat M_{\rm TF}(5))^{-1}\bar g.
\]

At the same EM worked point,
\[
\lambda_-\approx 5.92556258,
\qquad
\lambda_+\approx 237.91117494,
\]
\[
R_-\approx 0.00262800,
\qquad
R_+\approx 0.38665771,
\qquad
R_-+R_+=\frac{109}{280},
\]
with effective Padé pole
\[
\lambda_{\rm eff}^{\rm TF}\approx 188.17695898,
\]
and maximum relative error still tiny:
\[
\max\text{ rel.err.}\approx 7.10\times10^{-5}.
\]

#### Wall-completed dynamic monopole closure
The separated-order Family-1 sidewall and endcap reductions are then tightened by the direct coupled full-profile wall completion.

The final coupled full-profile inertia data are
\[
R_{\rm mass}^{\rm full}\approx 0.886313972989725,
\qquad
\widehat M_{aa}^{\rm full}\approx 0.563114968953987,
\qquad
\widehat M_{LL}^{\rm full}\approx 0.065829228119349.
\]

Using these values in the same geometry Hessian gives
\[
\lambda_-\approx 6.405572392138922,
\qquad
\lambda_+\approx 254.444968136936126,
\]
\[
R_-\approx 0.002552474771738,
\qquad
R_+\approx 0.386733239513976,
\qquad
R_-+R_+=\frac{109}{280}.
\]

So the final monopole wall channel is
\[
K_{00}(s)
=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+}.
\]

At zero frequency,
\[
K_{00}(0)=\frac{4}{45}
\]
exactly.

Relative to the sharp-wall TF baseline, the pole ratios are
\[
\frac{\Omega_-^2}{\Omega_{-,\rm TF}^2}\approx 1.2196655544,
\qquad
\frac{\Omega_+^2}{\Omega_{+,\rm TF}^2}\approx 1.2066780945.
\]

The one-pole Padé reduction remains excellent:
\[
\lambda_{\rm eff}\approx 202.923516367519028,
\qquad
\max\text{ rel.err.}\approx 6.89\times 10^{-5}
\quad \text{on }0\le s\le 0.1\lambda_-.
\]

So the scalar breathing branch is now near-final in the conservative low-frequency sense.

### 8.11 Final coupled Family-1 response module appendix

This appendix is the constructive endpoint of the throat-side program in the paper.

#### Static law and zero-frequency even/odd packaging remain unchanged
The direct coupled wall completion does **not** change the solved static support/source law:
\[
z_{\rm base}(\mu)=\frac{17}{56}-\frac{5}{56}\mu^2,
\qquad
t(\mu)=\frac{593}{672}-\frac{1553}{672}\mu^2+\frac78\mu^4,
\]
with exact source profile
\[
S(\mu)=10-\frac{63}{2}z_{\rm base}(\mu)=\frac{7}{16}+\frac{45}{16}\mu^2.
\]

The zero-frequency odd/even packaging also remains unchanged:
\[
L^{\rm add}_{\rm odd}
=
\frac12(v_A^2+v_B^2)L^{\rm wake}_{1\rm PN}
-\frac{15}{4}(U_A+U_B)(\mathbf v_A\cdot\mathbf v_B),
\]
\[
L^{\rm add}_{\rm even}
=
\Pi_A\cdot\Pi_B
+U_AJ\cdot\Pi_B
+U_BJ\cdot\Pi_A
+(J\cdot J-\Delta_{\rm geom})U_AU_B.
\]

With
\[
J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
\qquad
\Delta_{\rm geom}=\frac{281}{80},
\]
this gives
\[
J\cdot J = \frac{381}{80},
\qquad
\frac{381}{80}-\frac{281}{80}=\frac54,
\]
which is exactly the solved static cross coefficient.

So the zero-frequency residual against the exact added comparable-mass 2PN target still vanishes identically.

#### Direct coupled full-profile wall completion
The balanced thin-layer-consistent reference branch uses
\[
\alpha_r=10,
\quad \varepsilon_r=0.05,
\quad p_r=2,
\qquad
\chi_{\rm cap}=4,
\quad \varepsilon_z=0.05,
\quad \alpha_{\rm cap}=0.8,
\quad p_z=2,
\]
with full coupled profile
\[
\widetilde\rho(x,y)
=
\Bigg[
1-y^2
-\alpha_r S\Big(\frac{x-1}{\varepsilon_r}\Big)^{p_r}
-\alpha_{\rm cap} S\Big(\frac{y-1}{2\varepsilon_z}\Big)^{p_z}
\Bigg]_+^{1/4}.
\]

The exact coupled mass and inertia integrals then reproduce the full-profile values quoted above.

#### Separated-order approximation was already very accurate
Earlier separated-order values were already close:
\[
R_{\rm mass}^{\rm sep}\approx 0.8846236634,
\qquad
\widehat M_{aa}^{\rm sep}\approx 0.5623810783,
\qquad
\widehat M_{LL}^{\rm sep}\approx 0.0671965962.
\]

Relative corrections are tiny for overall mass and radial inertia, and still modest for axial inertia. So the coupled full-profile step is a genuine tightening, not a rescue of an unstable approximation.

#### Final practical summary of the throat-response program
The final coupled Family-1 response module is best summarized channel-theoretically:
\[
\boxed{
\text{carried odd dipole wake}
\;\oplus\;
\text{exact }P_0\oplus P_2\text{ even support/source sector}
\;\oplus\;
\text{coupled two-pole monopole breathing channel }K_{00}(s)
}
\]

The first two pieces reproduce the solved added comparable-mass conservative 2PN cross block exactly at zero frequency. The third gives the best current dynamic Family-1 normalization of the scalar monopole branch.

This is the appendix that justifies the paper’s statement that the **conservative 2PN Family-1 throat-response module is effectively in hand**.

---

## 9) Verification / WL / reproducibility ledger

The paper’s reproducibility archive is split across explicit WL harnesses and symbolic notes.

### 9.1 Inherited 1PN referee suite

Artifacts:
- `4d_gravity_1pn_master_harness.wl`
- `vector_wake_rebuild.wl`
- `4d_full_1pn_derivation_with_vectorwake.wl`

What they freeze:
- \(\kappa_\rho=1\),
- \(n=5\),
- \(\kappa_{\rm add}=1/2\),
- \(\kappa_{\rm PV}=3/2\),
- \(\beta_{\rm 1PN}=3\),
- parity-even wake coefficients,
- exact conservative EIH equality.

Reported summaries:
- `PASS count = 60`, `FAIL count = 0` for the broad master harness;
- wake module recovers
  \[
  a_H=0,
  \qquad
  a_L=\pm\frac{\sqrt3}{2},
  \qquad
  K_{\rm vec}=\frac{2}{\pi^2},
  \qquad
  C_\parallel=-\frac72,
  \qquad
  C_L=-\frac12;
  \]
- full 1PN assembly script reports
  `Passes: 37`, `Fails: 0`, `Skips: 0`.

### 9.2 Compact 2PN executable harness

Artifact:
- `4d_gravity_2pn_master_harness.wl`

What it checks:
- carried 1PN regressions;
- exact \(n=5\) Bernoulli continuation;
- minimal self exactification;
- quadratic local mass scaling hierarchy;
- generic cubic adiabatic-elimination template;
- self/static test-mass reduction;
- isotropic Schwarzschild target coefficients;
- universal finite-size gate.

Reported summary:
- `Passes: 42`, `Fails: 0`, `Skips: 0`.

### 9.3 Symbolic 2PN derivation ledger

Artifact:
- `4d_2pn_full_notes.md`

What it records:
- perturbative Legendre identity;
- compact residual basis;
- unique solved \(q_i,t_i,s_1\);
- exact zero-residual full ADM match;
- exact zero-frequency throat-module reconstruction;
- appendix-side constructive data.

Important interpretation:
- this is a symbolic derivation ledger, not yet a single all-in-one public PASS/FAIL harness.

### 9.4 What the archive verifies and what it does not

It verifies the paper **as stated**:
- conservative near-zone 2PN derivation,
- within the declared closure hierarchy,
- with explicit carried lower-order regressions,
- and explicit symbolic residual bookkeeping.

It does **not** verify:
- a fully solved moving-throat PDE,
- dissipation / radiation / higher-PN sectors,
- spin couplings,
- or a finished dynamic derivation of all appendix-side pole data.

---

## 10) Short list of the most important formulas to remember

If you only need the highest-value equations for future work, keep these.

### 10.1 Main lower-order carry-forward data
\[
\kappa_\rho=1,
\quad
n=5,
\quad
\kappa_{\rm add}=\frac12,
\quad
\kappa_{\rm PV}=\frac32,
\quad
\beta_{\rm 1PN}=3,
\quad
C_\parallel=-\frac72,
\quad
C_L=-\frac12.
\]

### 10.2 One-body 2PN closure
\[
D_{\rm eff}(u)=C_s(u)\left[1+\frac{32768}{3249}(G(u)-1)^2\right]
=1-4u+8u^2+O(u^3).
\]

### 10.3 Static 2PN gate
\[
\lambda_\rho=\frac12.
\]

### 10.4 Full 2PN block
\[
L_{\rm cons}=L_0+\frac{1}{c^2}L_1+\frac{1}{c^4}L_{2,\rm full}+O(c^{-6}).
\]

### 10.5 Exact match statement
\[
H_2=H_{2\rm PN}^{\rm ADM}.
\]

### 10.6 Zero-frequency throat operator data
\[
R_{1\perp}=\frac72,
\quad R_{10}=4,
\quad R_0=R_{20}=R_{21}=R_{22}=1,
\quad J=\left(\frac{4}{\sqrt5},\frac54,0,0,0,0\right),
\quad \Delta_{\rm geom}=\frac{281}{80}.
\]

### 10.7 Monopole wall closure
\[
\Delta K_{00}=\frac{109}{280},
\qquad
K_{00}(s)
=
-\frac{757}{2520}
+\frac{R_-}{1-s/\lambda_-}
+\frac{R_+}{1-s/\lambda_+}.
\]

### 10.8 Final constructive throat-response summary
\[
\boxed{
\text{carried odd dipole wake}
\;\oplus\;
\text{exact }P_0\oplus P_2\text{ even support/source sector}
\;\oplus\;
\text{coupled two-pole monopole breathing channel}
}
\]

That is the cleanest single-line summary of where the 2PN program stands.

---

## 11) Bottom-line interpretation

This paper does three things at once.

1. It upgrades the program from full conservative 1PN to full conservative 2PN.
2. It keeps the proof line clean and short by separating it from the larger constructive throat-response program.
3. It turns the new 2PN comparable-mass structure into a **small throat-response theory** rather than just a solved polynomial.

The main conservative two-body 2PN coefficient ledger is now closed within the declared hierarchy.

What remains open is **not** the static conservative answer. What remains open is the **dynamic throat-response completion**:
- derive the open pole scales,
- derive the full dynamic DtN operator from the inner-throat PDE,
- and then move beyond conservative near-zone physics.

That is why this paper is a major structural milestone in the project.
