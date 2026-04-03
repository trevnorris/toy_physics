4d_4pn_summary.md — Comprehensive summary (full conservative 4PN assembly paper)

## 0) What this paper is doing

This paper is the **full conservative 4PN assembly** for the 4D toy-model program, distilled from the long derivation record in `4d_4pn_full_notes.md` and formatted as a carry-forward session summary in the same spirit as `4d_3pn_summary.md`.

It starts from five things that are already frozen before the 4PN work begins:

1. the exact 4D parent reduction backbone from `4d.tex`;
2. the full conservative Newtonian + 1PN assembly from `4d_1pn_full.tex`;
3. the full conservative 2PN assembly from `4d_2pn.tex`;
4. the full conservative 3PN assembly from `4d_3pn.tex`;
5. the 2.5PN narrowing that reduced the universal higher-order bridge to the **orbital/worldtube STF quadrupole** branch.

The target is the complete conservative near-zone ledger through order
\[
c^{-8},
\]
i.e. the full conservative **4PN** sector.

That means closing, in one chain:

- the strict one-body 4PN Schwarzschild gate,
- the quartic ordinary/Hamiltonian compiler,
- the full **local instantaneous** 4PN generic-frame ordinary sector,
- the separate **hereditary/tail** coefficient,
- and the interface theorem showing whether 4PN opens any new normalization datum beyond the already-known 2.5PN quadrupole gap.

The strongest honest carry-forward statement is:

> **Within the same declared closure hierarchy used at 1PN, 2PN, and 3PN, the entire local instantaneous 4PN sector is already assembled exactly, and the only remaining condition for a full conservative 4PN theorem is the same passive/outgoing STF quadrupole normalization already isolated by the 2.5PN program.**

Equivalently,
\[
L_{\rm cons}=L_0+c^{-2}L_1+c^{-4}L_2+c^{-6}L_3+c^{-8}L_4+O(c^{-10}),
\]
and at 4PN the theorem split is
\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}},
\]
with the local piece already closed and the tail piece fixed by the same quadrupole normalization that controls the 2.5PN Burke–Thorne channel.

This is **not** claimed as an assumption-free theorem of a fully solved moving-throat PDE. It is a **full conservative 4PN derivation within a declared closure hierarchy**, with one sharply identified conditional interface left open.

---

## 1) Claim taxonomy and how to read the paper

The 4PN summary uses the same four claim classes as the 1PN, 2PN, and 3PN full papers.

### 1.1 Exact

These are symbolic identities once the lower-order ledger, target chart, and basis choices are fixed.

At 4PN this includes:

- the strict Schwarzschild one-body gate through \(c^{-8}\),
- the exact quartic perturbative Legendre compiler,
- the exact reduced COM local Hamiltonian/ordinary translation,
- the exact Hamiltonian-chart generic-frame local span result,
- the exact canonical-slice representative of the Hamiltonian-chart local residual,
- the exact generic-frame lift of the aligned-seed ordinary correction,
- the exact local referee reconstruction of the full local ordinary target,
- the exact universal Newtonian quadrupole shadow of the GR tail,
- and the exact coefficient bridge
  \[
  C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
  \]

### 1.2 Controlled reduction

These statements use the same already-declared reductions carried from the lower-order program:

- projection / Poisson-hook / worldtube reduction from the 4D parent theory,
- strict one-body and COM reductions,
- fixed-chart local Hamiltonian import,
- and the grouped real quadrupole reduction inherited from the 2.5PN narrowing.

So nothing here is being claimed as an unreduced theorem of the full moving-throat PDE.

### 1.3 Protocol closure

These are fixed only inside the same declared closure hierarchy already used successfully below 4PN.

At 4PN this includes:

- carrying forward the lower-order self/static ledger,
- using the Hamiltonian-aligned seed convention for the local ordinary 4PN chart,
- and treating the hereditary coefficient as controlled by the same passive/outgoing STF quadrupole normalization already isolated by 2.5PN.

### 1.4 Full assembly

Once the carried lower-order ledger, the fixed local ADM target, the quartic compiler, the Hamiltonian-first lift, and the hereditary bridge are inserted, the final equality is exact algebra inside the declared hierarchy:
\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]

What remains conditional is **not** another generic 4PN coefficient search. It is one already-narrow quadrupole-normalization theorem.

---

## 2) What the paper claims — and what it does not claim

### 2.1 What it claims

The paper claims all of the following.

1. The strict one-body 4PN Schwarzschild gate is fixed exactly.
2. The quartic perturbative Legendre-transform identity closes the exact ordinary/Hamiltonian bridge through 4PN.
3. The raw exchange-symmetric generic-frame local scaffold is already COM-complete blockwise, so local 4PN never had an early span shortage.
4. The fixed local 4PN lift must be done **Hamiltonian-first**, because the Hamiltonian degree ceilings match the constant-coefficient scaffold while the translated ordinary target develops extra \(\nu^4\) tails.
5. In the Hamiltonian chart there is **no span obstruction**: the generic-frame local interaction residual is fully reachable, block by block.
6. The canonical comparable-mass local residual translates back to the ordinary chart by exact sign flip once the lower-order ledger and aligned seed convention are fixed:
   \[
   \Delta H_4=-\Delta L_4.
   \]
7. The entire local instantaneous 4PN ordinary sector has an exact generic-frame representative.
8. The hereditary 4PN coefficient is not a new independent datum; it is fixed by the same canonically normalized STF quadrupole branch already isolated by the 2.5PN program:
   \[
   C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
   \]
9. Therefore the remaining 4PN gap is exactly the same passive/outgoing quadrupole normalization problem that already survived 2.5PN.

### 2.2 What it does **not** claim

The paper explicitly does **not** claim:

- a fully solved moving-throat bulk PDE theorem;
- an assumption-free microscopic derivation of the passive/outgoing quadrupole normalization;
- a final uniqueness/interpretation theorem for every aligned-seed local representative;
- any new dissipative or radiative theorem beyond the already narrowed quadrupole route;
- spin couplings, strong-field completion, or higher-PN completion beyond the present conservative 4PN scope.

So the correct reading is:

> the conservative 4PN two-body sector is closed **except for one already-known quadrupole normalization condition**, and that remaining condition is inherited from the 2.5PN program rather than opened freshly by 4PN.

---

## 3) Lower-order inputs, notation firewall, and bookkeeping conventions

### 3.1 Carried lower-order conservative backbone

Frozen from earlier papers:
\[
\kappa_\rho=1,
\qquad
n=5,
\qquad
\kappa_{\rm add}=\frac12,
\qquad
\kappa_{\rm PV}=\frac32,
\qquad
\beta_{\rm 1PN}=3.
\]

The exact 4D parent reduction backbone also carries forward unchanged:

- exact projection map,
- exact projected continuity with leakage,
- exact longitudinal identity,
- quasi-static Poisson hook,
- and coherent-defect / worldtube particle reduction.

From 2PN onward the defect is already more than a scalar monopole. The carried ontology includes:

- an odd dipole wake sector,
- a real \(P_0\oplus P_2\) mouth/support layer,
- and a separate geometry-closure channel.

From 3PN the conservative grouped real \(P_2\) sector is already the live higher-order conservative payload, while the remaining universal dissipative route had already been narrowed at 2.5PN to the orbital/worldtube STF quadrupole branch.

### 3.2 Why the 2.5PN narrowing controls 4PN

The 2.5PN audit already demoted the scalar and dipole dangers above universal point-particle 2.5PN on the strict small-body branch, and isolated the surviving retarded quadrupole branch as
\[
\mathcal M_{ab}^R(\omega)
=
\bigl[m_0+m_2\omega^2+m_4\omega^4+i\Gamma_5\omega^5+\cdots\bigr]\delta_{ab},
\qquad
m_0\neq 0,
\quad
\Gamma_5>0.
\]

So 4PN should be read as the next **conservative** use of a quadrupole route that 2.5PN had already isolated, not as the opening of a new unrelated bridge.

### 3.3 PN bookkeeping and variables

The conservative expansion convention is
\[
L_{\rm cons}=L_0+c^{-2}L_1+c^{-4}L_2+c^{-6}L_3+c^{-8}L_4+O(c^{-10}).
\]
So throughout this summary, “4PN” means the coefficient of \(c^{-8}\).

Generic-frame two-body invariants:
\[
\mathbf r=\mathbf x_A-\mathbf x_B,
\qquad
r=|\mathbf r|,
\qquad
\mathbf n=\mathbf r/r,
\]
\[
v_A^2=\mathbf v_A^2,
\qquad
v_B^2=\mathbf v_B^2,
\qquad
v_{AB}=\mathbf v_A\!\cdot\!\mathbf v_B,
\]
\[
v_{An}=\mathbf v_A\!\cdot\!\mathbf n,
\qquad
v_{Bn}=\mathbf v_B\!\cdot\!\mathbf n.
\]

Local pair potentials:
\[
U_A=\frac{Gm_B}{r},
\qquad
U_B=\frac{Gm_A}{r}.
\]
In the strict one-body gate:
\[
U=\frac{GM}{r}.
\]

Reduced COM notation:

- ordinary-Lagrangian coefficients: \(l_i\),
- Hamiltonian coefficients: \(h_i\),
- residuals beyond the chosen seed: \(\Delta l_i,\Delta h_i\).

For the raw local generic-frame scaffold it is convenient to use the invariant shorthand
\[
a=v_A^2,
\qquad
b=v_B^2,
\qquad
c=\mathbf v_A\!\cdot\!\mathbf v_B,
\qquad
d=\mathbf v_A\!\cdot\!\mathbf n,
\qquad
e=\mathbf v_B\!\cdot\!\mathbf n,
\]
\[
p=m_A,
\qquad
q=m_B.
\]

### 3.4 COM versus generic-frame conventions

The settled 4PN convention is:

- use COM variables whenever they simplify the exact target or the reduced compiler map,
- but treat the **Hamiltonian chart** as the authoritative chart for the 4PN local generic-frame lift.

The key lesson is that the translated ordinary local target develops extra \(\nu^4\) tails in the middle and upper interaction blocks, while the fixed local Hamiltonian target matches the constant-coefficient scaffold degree ceilings exactly. So the correct 4PN sequence is:
\[
\text{fixed local Hamiltonian target}
\;\longrightarrow\;
\text{generic-frame Hamiltonian lift}
\;\longrightarrow\;
\text{ordinary translation only afterward}.
\]

### 3.5 Notation firewall

The corrected 4D / EM / PN notation firewall remains in force.

1. Electric charge is **not** part of the conservative 4PN proof.
2. When charge appears at all, it uses
   \[
   \eta_Q=\pm1,
   \qquad
   q_*=\eta_Q e_*,
   \qquad
   q_{\rm eff}=\frac{q_*}{\sqrt{Z_{\rm int}}}.
   \]
3. Quantized circulation belongs to the magnetic/vortical sector, not to electric charge.
4. The historical gravity-side scalar coefficient once written as bare \(q=1\) remains
   \[
   \kappa_\rho=1.
   \]
5. The grouped quadrupole labels \(20,21,22\) always refer to grouped real \(P_2\) channels, not tensor indices.
6. In the raw 4PN generic-frame scaffold, the symbols \(a,b,c,d,e,p,q\) are invariant shorthands, whereas in the quadrupole-branch formulas below the throat-radius scale will be written explicitly as \(a_{\rm th}\) to avoid conflict.

---

## 4) Headline outputs and carry-forward constants

### 4.1 Exact one-body 4PN Lagrangian target

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

### 4.2 Minimal one-body repair ledger

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
So 4PN opens

- a quartic static gate,
- a quartic denominator datum,
- and two genuinely new self-sector slots.

### 4.3 Exact one-body 4PN Hamiltonian gate

In reduced COM Hamiltonian variables, the strict one-body 4PN Hamiltonian gate is
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
This is the exact one-body Hamiltonian seed data used by the local Hamiltonian-first lift.

### 4.4 Exact quartic perturbative Legendre compiler

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
This is the exact ordinary/Hamiltonian bridge used throughout the local 4PN chain.

### 4.5 Natural local 4PN self/static seed

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

### 4.6 Raw local scaffold sizes and COM slot counts

The raw exchange-symmetric generic-frame local interaction scaffold has block sizes
\[
G/r:52,
\qquad
G^2/r^2:46,
\qquad
G^3/r^3:29,
\qquad
G^4/r^4:10,
\qquad
G^5/r^5:2,
\]
for a total of \(139\) raw interaction directions.

After COM projection these span exactly the full 15 interaction slots of the reduced local 4PN sector:
\[
G/r:5,
\qquad
G^2/r^2:4,
\qquad
G^3/r^3:3,
\qquad
G^4/r^4:2,
\qquad
G^5/r^5:1.
\]
So local 4PN never had an early kinematic shortage.

### 4.7 Hamiltonian-chart degree ceilings and exact rank theorem

The decisive local structural comparison is:

- in the **Hamiltonian** chart the imported local interaction blocks have \(\nu\)-degree ceilings
  \[
  (4,3,3,2,2)
  \]
  across the \(G/r, G^2/r^2, G^3/r^3, G^4/r^4, G^5/r^5\) blocks;
- after quartic translation the **ordinary** target becomes
  \[
  (4,4,4,4,2),
  \]
  so the translated ordinary target develops extra \(\nu^4\) tails in the middle and upper blocks.

That is why the lift has to be Hamiltonian-first.

In the Hamiltonian chart the blockwise coefficient-space ranks are maximal:
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
So the full interaction coefficient-space dimensions are
\[
\text{rank}=47,
\qquad
\text{nullity}=92.
\]
There is **no Hamiltonian-chart span obstruction**.

### 4.8 Canonical quotient slice and aligned-seed liftability

The 92-dimensional Hamiltonian-chart null family is purely COM-blind algebraic freedom. Its blockwise nullities are
\[
Q:32,
\qquad
T:34,
\qquad
S:20,
\qquad
U:6,
\qquad
W:0.
\]

Fixing the canonical quotient slice by setting all null coordinates to zero gives one sparse exact generic-frame Hamiltonian representative with
\[
Q:11,
\qquad
T:12,
\qquad
S:9,
\qquad
U:4,
\qquad
W:2
\]
nonzero scaffold directions.

After translation back to the ordinary chart, the remaining obstruction was found to be purely a **seed-alignment** issue:
\[
\delta L_{4,\rm seed}
=
L_{4,\rm seed}^{(\rm aligned)}-L_{4,\rm seed}^{(\rm natural)}.
\]
This correction is exactly generic-frame liftable in the structured seed spaces
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
A canonical sparse realization of the aligned-seed lift has nonzero-direction counts
\[
Q:13,
\qquad
T:14,
\qquad
S:12,
\qquad
U:8,
\qquad
W:0.
\]

### 4.9 Final exact local theorem

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
Then the local referee master verifies that the COM reduction of this generic-frame candidate reproduces the full fixed-chart local ordinary 4PN target block by block.

So the **entire local instantaneous 4PN ordinary interaction sector is now assembled exactly**.

### 4.10 Exact GR tail functional and local logarithmic shadow

The standard conservative 4PN GR tail coefficient is
\[
C_{\rm tail}^{\rm GR}=\frac{G^2M}{5c^8}.
\]
The time-symmetric nonlocal Hamiltonian is
\[
H_{4\rm PN}^{\rm tail,sym}(t)
=
-\frac{1}{5}\frac{G^2M}{c^8}
I_{ij}^{(3)}(t)
\,\mathrm{Pf}_{2s/c}\!\int_{-\infty}^{+\infty}\frac{dv}{|v|}
I_{ij}^{(3)}(t+v),
\]
and the local logarithmic shadow is controlled by
\[
F_{\rm tail}(t)=\frac{2}{5}\frac{G^2M}{c^8}\bigl(I_{ij}^{(3)}(t)\bigr)^2.
\]
Using the Newtonian STF quadrupole
\[
I_{ij}=\mu\,\mathrm{STF}(x_i x_j),
\]
one obtains the exact order-reduced third derivative
\[
I_{ij}^{(3)}
=
-\frac{2GM\mu}{r^3}
\left(4x_{\langle i}v_{j\rangle}-3\frac{\dot r}{r}x_{\langle i}x_{j\rangle}\right),
\]
so that
\[
\bigl(I_{ij}^{(3)}\bigr)^2
=
\frac{8}{3}\frac{G^2M^2\mu^2}{r^4}
\left(12v^2-11\dot r^2\right),
\]
and therefore
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
This occupies only the degree-2 \(G^4/r^4\) block.

### 4.11 Exact bridge to the 2.5PN coefficient and toy-model branch forms

Let the canonically normalized STF quadrupole odd coefficient be \(\gamma_{\rm quad}^{\rm eff}\). The exact bridge is
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

The toy-model quadrupole branch forms inherited from the 2.5PN program are
\[
\gamma_{\rm quad}^{\rm eff}
=
\mathcal N_Q\frac{a_{\rm th}^5}{27c_s^5}
=
\overline\Gamma_5
=
9\frac{\overline K_2^{5/2}}{\overline K_0^{3/2}}.
\]
Hence
\[
C_{\rm tail}^{\rm toy}
=
\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}
=
\frac{GM}{2c^3}\,\overline\Gamma_5.
\]
On the Burke–Thorne target branch,
\[
\gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5}
\quad\Longrightarrow\quad
C_{\rm tail}^{\rm toy}=\frac{G^2M}{5c^8}=C_{\rm tail}^{\rm GR}.
\]
Equivalently,
\[
\mathcal N_Q^{\rm target}=\frac{54Gc_s^5}{5a_{\rm th}^5c^5},
\]
and in the invariant-pair language the Burke–Thorne target branch is
\[
\overline K_0^{\rm target}=\frac{64G}{45c^5}\Omega_Q^5,
\qquad
\overline K_2^{\rm target}=\frac{16G}{45c^5}\Omega_Q^3.
\]

### 4.12 No new independent 4PN normalization datum

Because the hereditary coefficient is fixed by
\[
C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff},
\]
**no new independent normalization datum opens at 4PN** beyond the same passive/outgoing STF quadrupole normalization already isolated by the 2.5PN program.

### 4.13 Final theorem split

The final 4PN split is therefore
\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}},
\]
with

- the local sector already closed exactly in generic-frame ordinary form,
- the tail coefficient fixed by the same canonically normalized STF quadrupole branch as 2.5PN,
- and no new 4PN-specific normalization constant left to determine.

---

## 5) Main local 4PN proof chain

### 5.1 Strict one-body 4PN gate

The chain begins with the exact Schwarzschild one-body gate and the observation that the 3PN self-sector packaging does not extend trivially to 4PN. The minimal one-body repair ledger already shows that 4PN opens a quartic static gate, a quartic denominator correction, and two genuinely new self-sector slots.

### 5.2 Exact quartic compiler

The quartic perturbative Legendre-transform identity is then frozen before any comparable-mass local solve is attempted. This is what makes the later Hamiltonian-first program exact rather than heuristic.

### 5.3 Raw local scaffold and COM completeness

The exchange-symmetric generic-frame local scaffold is constructed block by block. Its raw size is large (139 directions), but after COM projection it spans exactly the 15 interaction slots of the reduced local 4PN sector. So the right question is no longer “is there enough room?” but “which chart and which quotient slice should be used?”

### 5.4 Fixed local target import and the Hamiltonian-first lesson

Importing the fixed local ADM Hamiltonian target shows that the Hamiltonian degree ceilings match the constant-coefficient scaffold exactly, while the translated ordinary target develops extra \(\nu^4\) tails in the middle and upper interaction blocks. This proves that the generic-frame local lift must be performed in the Hamiltonian chart first.

### 5.5 Hamiltonian-chart generic-frame lift

In the Hamiltonian chart, the blockwise ranks are maximal, so there is no span obstruction. The only ambiguity is a 92-dimensional COM-blind null family. That ambiguity is algebraic rather than physical.

### 5.6 Canonical quotient slice

Setting all null coordinates to zero gives one exact canonical generic-frame Hamiltonian representative of the local comparable-mass residual. This solves the existence problem in the Hamiltonian chart.

### 5.7 Ordinary translation and seed alignment

The canonical comparable-mass residual translates to the ordinary chart by exact sign flip once the lower-order ledger and aligned seed are fixed:
\[
\Delta H_4=-\Delta L_4.
\]
The earlier ordinary-chart obstruction is then exposed as purely a **seed-alignment** issue, not a failure of the comparable-mass local lift itself.

### 5.8 Generic-frame lift of the aligned-seed correction

The aligned-seed correction is shown to be generic-frame liftable in the structured seed sectors
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
So the old ordinary obstruction is repaired exactly inside a narrow, intelligible seed sector.

### 5.9 Local referee reconstruction

Combining

- the natural generic-frame seed,
- the generic-frame lift of the aligned-seed correction,
- and the canonical comparable-mass local residual,

produces one exact generic-frame ordinary representative of the **entire local instantaneous 4PN sector**. The local referee suite verifies this block by block in the COM chart.

---

## 6) Tail / hereditary 4PN bridge

### 6.1 Exact GR tail structure

The conservative 4PN GR dynamics splits into a local instantaneous sector plus a time-symmetric nonlocal tail functional with universal coefficient
\[
C_{\rm tail}^{\rm GR}=\frac{G^2M}{5c^8}.
\]

### 6.2 Universal Newtonian quadrupole shadow

The local logarithmic shadow of the tail sits only in the degree-2 \(G^4/r^4\) block:
\[
\frac{F_{\rm tail}}{\mu}
=
\frac{16}{15}\,\nu\,\frac{U^4}{c^8}
\left(12v^2-11\dot r^2\right).
\]
So the hereditary side is structurally much smaller than the local lift problem.

### 6.3 Exact coefficient bridge to the 2.5PN quadrupole route

The hereditary coefficient is related exactly to the canonically normalized STF quadrupole odd coefficient by
\[
C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
\]
This is the cleanest bridge between the already narrowed 2.5PN route and the conservative 4PN tail sector.

### 6.4 No new 4PN-specific normalization gap

Because the tail coefficient is entirely controlled by \(\gamma_{\rm quad}^{\rm eff}\), 4PN does **not** open a new independent normalization constant. Once the passive/outgoing STF quadrupole normalization is fixed on the natural branch, the 2.5PN and 4PN coefficients close together.

---

## 7) Final conditional conservative 4PN theorem

### 7.1 Theorem (conditional conservative 4PN statement)

Let
\[
L_{4\mathrm{PN}}^{\mathrm{cons}}
=
L_{4\mathrm{PN}}^{\mathrm{local}}
+
L_{4\mathrm{PN}}^{\mathrm{tail}}.
\]
Then, inside the declared closure hierarchy:

1. **Local sector.** The entire local instantaneous 4PN ordinary sector already has an exact generic-frame representative.
2. **Hereditary coefficient.** The unique compatible tail coefficient is
   \[
   C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
   \]
3. **No new 4PN normalization datum.** The 4PN tail sector introduces no new quadrupole-normalization constant beyond the same passive/outgoing STF quadrupole normalization already isolated by the 2.5PN audit.
4. **GR recovery condition.** If
   \[
   \gamma_{\rm quad}^{\rm eff}=\frac{2G}{5c^5},
   \]
   then automatically
   \[
   C_{\rm tail}=\frac{G^2M}{5c^8},
   \]
   and the full conservative 4PN local+tail coefficient structure matches the standard GR target.

So the remaining gap between the present hierarchy and a fully unconditional conservative 4PN theorem is **exactly the same narrow passive/outgoing quadrupole normalization problem already isolated by the 2.5PN program**.

### 7.2 What is solved, what is conditional, and what remains open

Already solved inside the present hierarchy:

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

Conditional part:
\[
\gamma_{\rm quad}^{\rm eff}\stackrel{?}{=}\frac{2G}{5c^5}
\]
on the natural passive/outgoing STF quadrupole branch.

Still open beyond the present hierarchy:

- a first-principles moving-throat PDE derivation of the passive/outgoing quadrupole normalization;
- the deeper microscopic derivation of the grouped real conservative quadrupole data;
- the final uniqueness/interpretation theorem for the aligned-seed local sector;
- and any dissipative/radiative completion beyond the already narrowed quadrupole route.

---

## 8) Verification / reproducibility ledger

The core 4PN package artifacts are:

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

These scripts collectively verify:

- the strict one-body gate,
- the quartic ordinary/Hamiltonian compiler,
- the local Hamiltonian import and Hamiltonian-first structural lesson,
- the Hamiltonian-chart generic-frame lift and canonical slice,
- the ordinary residual translation and aligned-seed correction,
- the local referee reconstruction of the full ordinary target,
- the universal quadrupole tail shadow,
- and the exact 4PN hereditary bridge to the 2.5PN STF quadrupole coefficient.

What the suite does **not** verify is the missing moving-throat PDE-side normalization of the passive/outgoing quadrupole branch.

---

## 9) Mental checklist for future work

If I pick this program up again later, the non-negotiable 4PN carry-forward points are:

1. **4PN is split, not monolithic.** Always treat
   \[
   L_{4\mathrm{PN}}^{\mathrm{cons}}=L_{4\mathrm{PN}}^{\mathrm{local}}+L_{4\mathrm{PN}}^{\mathrm{tail}}
   \]
   as the basic theorem envelope.
2. **The local 4PN problem is already solved for existence.** Do not reopen local span questions unless I am specifically studying uniqueness or interpretation of the aligned-seed sector.
3. **The local lift must be Hamiltonian-first.** The ordinary-chart \(\nu^4\)-tail obstruction is a chart artifact tied to quartic translation, not a failure of the Hamiltonian local scaffold.
4. **The true ordinary obstruction was seed alignment, not comparable-mass local liftability.** The aligned-seed correction is generic-frame liftable in a narrow structured seed sector.
5. **The hereditary problem is tiny compared with the local one.** Its local logarithmic shadow lives only in the degree-2 \(G^4/r^4\) block.
6. **The 4PN tail coefficient is not a new free constant.** It is fixed by the same canonically normalized STF quadrupole branch as 2.5PN:
   \[
   C_{\rm tail}=\frac{GM}{2c^3}\,\gamma_{\rm quad}^{\rm eff}.
   \]
7. **No new independent normalization datum opens at 4PN.** The only remaining theorem gate is the same passive/outgoing quadrupole normalization already exposed by the 2.5PN audit.
8. **The next real theorem target is not more local algebra.** It is to derive the passive/outgoing STF quadrupole normalization — equivalently the canonical pair \((\overline K_0,\overline K_2)\) on the natural branch — from the moving-throat side.
9. **Once that quadrupole normalization closes, 2.5PN and 4PN close together.** That is the whole strategic reason the present 4PN theorem is worth stating in this conditional form.
