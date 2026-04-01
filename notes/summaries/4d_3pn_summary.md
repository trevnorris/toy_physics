4d_3pn.tex — Comprehensive summary (full conservative 3PN assembly paper)

## 0) What this paper is doing

This paper is the **full conservative two-body 3PN assembly** for the 4D toy-model program.

It starts from four things that are already frozen before the paper begins:

1. the exact 4D parent reduction backbone from `4d.tex`:
   - exact projection map,
   - exact projected continuity with leakage,
   - exact longitudinal identity,
   - controlled quasi-static Poisson hook,
   - controlled coherent-defect / worldtube particle reduction;

2. the full conservative Newtonian + 1PN assembly from `4d_1pn_full.tex`, including exact EIH equality in the solved scope;

3. the full conservative 2PN assembly from `4d_2pn.tex`, including exact ADM equality in the solved scope;

4. the 2.5PN narrowing that reduced the live higher-order universal point-particle route to the **orbital/worldtube STF quadrupole branch**, packaged conservatively as the grouped real
   \[
   P_{20}\oplus P_{21}\oplus P_{22}
   \]
   bundle.

The paper’s main target is the complete conservative near-zone two-body ledger through order
\[
c^{-6},
\]
i.e. the full conservative **3PN** sector.

That means closing all of the following at once:

- free-particle octic kinematics \(v^8\),
- one-body / self terms through \(U^4\),
- the full comparable-mass middle interaction block,
- the static quartic term at order \(G^4/r^4\),
- and the generic-frame uniqueness problem in the fixed ADM chart.

The strongest carry-forward statement is:

> **Within the same declared closure hierarchy used at 1PN and 2PN, the conservative 3PN two-body ledger is now fully assigned in the fixed ADM chart.**

The final 3PN split is
\[
\Delta \hat L_3^{\rm GR}
=
\underbrace{\Delta l_1 v^8}_{\text{compiler image of the free COM Hamiltonian}}
+
\underbrace{L_{P_2}^{\rm mid}}_{\text{exact grouped real }P_2\text{ closure}}
+
\underbrace{\Delta l_{15}^{(g)}U^4}_{\text{unique geometry completion}}.
\]

This is **not** claimed as an assumption-free theorem of the fully solved moving-throat PDE. It is a **full conservative 3PN derivation within a declared closure hierarchy**.

---

## 1) Claim taxonomy and how to read the paper

The 3PN paper uses the same four claim classes as the 1PN and 2PN full papers.

### 1.1 Exact

These are identities or exact algebraic statements once the lower-order ledger, basis choice, and target chart are fixed.

Examples:

- the exact isotropic Schwarzschild one-body target through \(c^{-6}\),
- the exact cubic perturbative Legendre-transform identity,
- the exact diagonal COM ordinary/Hamiltonian compiler,
- the exact generic-frame Hamiltonian compiler in the fixed ADM chart,
- the exact invertibility of the richer grouped real \(P_2\) middle-block map,
- the exact sigma-collapse that removes the apparent static \(P_0/g\) ambiguity,
- and the exact pure-kinetic collapse to the universal free relativistic two-body COM Hamiltonian.

### 1.2 Controlled reduction

These use an explicit regime assumption or reduced description.

Examples:

- the carried 4D-to-brane projection and Poisson-hook reduction,
- the coherent-defect / worldtube particle reduction,
- the reduced isotropic COM 15-slot basis,
- the grouped real \(P_2\) packaging of the higher-order conservative quadrupole sector,
- and the Hamiltonian-first COM reduction rule used to fix the generic-frame representative.

### 1.3 Protocol closure

These are fixed only **within a declared response hierarchy** rather than by a fully solved bulk PDE.

Examples:

- carrying forward the lower-order adiabatic/response closures,
- the one-body 3PN repair ledger as the next target inside the same self-sector packaging,
- the richer local grouped-\(P_2\) constitutive family \((T_A,S_A,V_A)\),
- and the interpretation of the leftover static slot as a geometry-side completion.

### 1.4 Full assembly

Once the carried lower-order ledger, the exact 3PN targets, the fixed ADM chart, and the declared closure hierarchy are inserted, the final equality to the standard target is exact algebra.

That is the status of
\[
L_{\rm cons}=L_0+c^{-2}L_1+c^{-4}L_2+c^{-6}L_3+O(c^{-8}),
\qquad
\Delta H_3=-\Delta L_3,
\]
together with the final 3PN split quoted above.

---

## 2) What the paper claims — and what it does not claim

### 2.1 What it claims

The paper claims all of the following.

1. The exact 4D parent reduction backbone carries forward unchanged.
2. The strict one-body 3PN Schwarzschild gate fixes a minimal three-parameter repair ledger:
   \[
   \mu_{\rho,3}=\frac14,
   \qquad
   d_3=-\frac{45}{4},
   \qquad
   s_{24}=-\frac1{16}.
   \]
3. The exact isotropic COM target is imported, solved, and reduced to a finite 15-slot comparable-mass residual vector.
4. In the fixed ADM chart the exact generic-frame cubic compiler becomes
   \[
   \Delta H_3=-\Delta L_3,
   \]
   so the generic-frame ordinary representative is unique once the Hamiltonian target is imposed.
5. The richer grouped real \(P_2\) local family
   \[
   (T_{20},T_{21},T_{22},S_{20},S_{21},S_{22},V_{20},V_{21},V_{22})
   \]
   closes the full 9-slot 3PN middle block exactly.
6. The remaining static slot is uniquely geometry-side.
7. The apparent pure-kinetic leftover is only the ordinary-Lagrangian counterimage of the universal free relativistic two-body COM Hamiltonian.
8. Therefore the conservative 3PN bottleneck is no longer algebraic. What remains is the deeper moving-throat derivation of the grouped-\(P_2\) constitutive law and the geometry completion.

### 2.2 What it does **not** claim

The paper explicitly does **not** claim:

- a fully solved moving-throat bulk PDE theorem;
- an assumption-free microscopic derivation of every grouped-\(P_2\) or geometry coefficient;
- any dissipative, radiative, or 2.5PN outgoing-normalization result;
- spin couplings, radiation reaction, or strong-field completion;
- or a completed dynamic derivation of the passive/outgoing quadrupole normalization on the true moving-throat branch.

So the correct reading is:

> the conservative 3PN two-body sector is closed **within the declared hierarchy**, while the fully dynamical throat-response problem remains open.

---

## 3) Lower-order inputs, notation firewall, and bookkeeping conventions

### 3.1 Carried Newtonian / 1PN / 2PN backbone

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

From 2PN the defect is already known to be more than a scalar monopole. The carried ontology includes:

- an odd dipole wake sector,
- a real \(P_0\oplus P_2\) mouth/support layer,
- and a separate geometry-closure channel.

### 3.2 Why the 2.5PN notes narrow the next conservative target to the grouped real \(P_2\) bundle

The 2.5PN program already narrowed the live universal point-particle route to the orbital/worldtube STF quadrupole branch and repackaged the conservative higher-order question as the grouped real
\[
P_{20}\oplus P_{21}\oplus P_{22}
\]
bundle. That means the real higher-order conservative payload is no longer “derive more PN somehow,” but rather: determine the grouped low-frequency even coefficients, test isotropy, and then feed that data into the remaining narrow 2.5PN outgoing-normalization bridge.

### 3.3 PN bookkeeping and invariant variables

The conservative expansion convention is
\[
L_{\rm cons}=L_0+c^{-2}L_1+c^{-4}L_2+c^{-6}L_3+O(c^{-8}).
\]
So throughout this summary, “3PN” means the coefficient of \(c^{-6}\).

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
- residuals beyond the frozen self/static seed: \(\Delta l_i,\Delta h_i\).

Useful 3PN grouped notation:
\[
d\equiv \dot r,
\qquad
u_\perp^2\equiv v^2-d^2,
\qquad
\nu\equiv\frac{m_A m_B}{(m_A+m_B)^2}.
\]

### 3.4 COM versus generic-frame conventions

The settled convention is:

- use the COM basis whenever it simplifies the exact target or the diagonal compiler map,
- but treat the **generic-frame Hamiltonian compiler** as the authoritative uniqueness statement.

The important lesson is that a naive ordinary-Lagrangian COM reduction of the imported generic-frame 3PN target does **not** reproduce the exact COM ordinary target. The correct route is:

1. import the generic-frame target,
2. apply the exact cubic Legendre lift first,
3. then reduce to COM.

That Hamiltonian-first route recovers the exact standard 3PN COM Hamiltonian.

### 3.5 Notation firewall

The corrected 4D / EM / PN notation firewall remains in force.

1. Electric charge is **not** part of the 3PN conservative proof.
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
5. The grouped labels \(20,21,22\) always refer to grouped real \(P_2\) channels, not tensor indices.
6. In the static-completion section it is useful to write
   \[
   p\equiv m_A,
   \qquad
   q\equiv m_B.
   \]

---

## 4) Headline outputs and carry-forward constants

### 4.1 Exact one-body 3PN target

The exact isotropic Schwarzschild one-body 3PN target is
\[
\boxed{
\frac{L_{3\mathrm{PN},1\mathrm{body}}}{m}
=
\frac5{128}v^8
+\frac{11}{16}Uv^6
+\frac{47}{16}U^2v^4
+\frac{13}{8}U^3v^2
-\frac18U^4.
}
\]
This is the coefficient of \(c^{-6}\).

### 4.2 Minimal one-body repair ledger

Extending the carried 2PN denominator-style self sector to 3PN shows that a single cubic denominator repair is **not** enough. The exact one-body repair ledger is
\[
\boxed{
\mu_{\rho,3}=\frac14,
\qquad
 d_3=-\frac{45}{4},
\qquad
 s_{24}=-\frac1{16}.
}
\]
Interpretation:

- \(\mu_{\rho,3}=1/4\) is the cubic static gate,
- \(d_3=-45/4\) is the cubic denominator datum,
- \(s_{24}=-1/16\) is the genuinely new self slot needed in the \(U^2v^4\) channel.

### 4.3 Natural 3PN self/static seed

Promoting the exact one-body target symmetrically gives the natural two-body self/static seed
\[
\boxed{
L_{3,\mathrm{seed}}
=
\frac{5}{128}(m_Av_A^8+m_Bv_B^8)
+\frac{11Gm_A m_B}{16r}(v_A^6+v_B^6)
}
\]
\[
\boxed{
\qquad
+\frac{47G^2m_A m_B}{16r^2}(m_Bv_A^4+m_Av_B^4)
+\frac{13G^3m_A m_B}{8r^3}(m_B^2v_A^2+m_A^2v_B^2)
-\frac{G^4m_A m_B}{8r^4}(m_B^3+m_A^3).
}
\]
Again, this whole block is the coefficient of \(c^{-6}\).

### 4.4 Exact COM compiler data

In the standard isotropic 15-slot COM basis, the cubic Legendre map is diagonal. The most important carry-forward formulas are
\[
h_1=-l_1+\frac{3}{16}\nu-\frac{21}{16}\nu^2+\frac94\nu^3,
\qquad
h_2=-l_2,
\qquad
h_3=-l_3,
\qquad
h_4=-l_4,
\qquad
h_5=-l_5,
\]
\[
h_6=-l_6+\frac14+\frac78\nu-\frac{35}{8}\nu^2-\frac{21}{4}\nu^3,
\]
and similarly slot by slot for the remaining coefficients. The key structural point is that, once the lower-order ledger is frozen, the COM ordinary/Hamiltonian compiler is exact, explicit, and invertible.

### 4.5 Exact COM GR target and residual data

Subtracting the frozen self/static seed from the imported exact GR COM target gives the exact comparable-mass residual vector
\[
\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16},
\qquad
\Delta l_2=\Delta l_3=\Delta l_4=\Delta l_5=0,
\]
\[
\Delta l_6=\frac{\nu(38-116\nu-57\nu^2)}{16},
\qquad
\Delta l_7=\frac{\nu^2(20-69\nu)}{16},
\qquad
\Delta l_8=\frac{3\nu^2(3-11\nu)}{16},
\qquad
\Delta l_9=\frac{5\nu^3}{16},
\]
\[
\Delta l_{10}=\frac{\nu(129-98\nu+52\nu^2)}{16},
\qquad
\Delta l_{11}=\frac{\nu(-3+52\nu+124\nu^2)}{16},
\qquad
\Delta l_{12}=\frac{\nu(-5+11\nu+48\nu^2)}{12},
\]
\[
\Delta l_{13}=-\frac{\nu(244+3\pi^2+1272\nu+96\nu^2)}{192},
\qquad
\Delta l_{14}=\frac{\nu(452+3\pi^2-384\nu-224\nu^2)}{64},
\]
\[
\Delta l_{15}=\frac{\nu(-908+63\pi^2)}{96}.
\]
This 15-slot vector is the exact 3PN target data against which every later grouped-\(P_2\), geometry, and compiler statement is checked.

### 4.6 Generic-frame uniqueness theorem

The generic-frame interaction scaffold has block sizes
\[
24,\qquad 17,\qquad 8,\qquad 1
\]
for the \(G/r\), \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks, for a total of 50 interaction directions.

But once the lower-order 1PN/2PN ledger is frozen and the Hamiltonian target is imposed, the exact cubic compiler reduces to the slotwise sign flip
\[
\boxed{\Delta H_3=-\Delta L_3.}
\]
So in the fixed ADM chart the generic-frame compiler is simply
\[
\boxed{C_{\rm gen}=-I_{50}.}
\]

### 4.7 Exact grouped real \(P_2\) middle-block closure

The grouped real conservative 3PN front end is controlled by the three \(O(\omega^2)\) data
\[
u_2^{(20)},\qquad u_2^{(21)},\qquad u_2^{(22)},
\]
or equivalently by
\[
\bar u_2=\frac{u_2^{(20)}+2u_2^{(21)}+2u_2^{(22)}}{5},
\qquad
 a_2=\frac{2u_2^{(20)}-u_2^{(21)}-u_2^{(22)}}{10},
\qquad
 b_2=\frac{u_2^{(21)}-u_2^{(22)}}{2}.
\]

The exact 3PN grouped closure uses the richer local family
\[
(T_{20},T_{21},T_{22},S_{20},S_{21},S_{22},V_{20},V_{21},V_{22}),
\]
with definitions given in Section 6 below. The exact middle-block matrix has constant determinant
\[
\boxed{\det M_{\rm mid}=-\frac4{27}.}
\]
So the richer grouped compiler is exactly invertible on the full nine-slot middle block \((l_6,\dots,l_{14})\).

### 4.8 Unique geometry completion of the static slot

The richer grouped compiler predicts a specific static companion
\[
\boxed{
\Delta l_{15}^{(P_2\,\rm pred)}
=
\frac{\nu(293-308\nu-102\nu^2)}{24}.
}
\]
The remaining static gap is
\[
\boxed{
\Delta l_{15}^{(0/g)}
=
\Delta l_{15}^{\rm GR}-\Delta l_{15}^{(P_2\,\rm pred)}
=
\frac{\nu(408\nu^2+1232\nu-2080+63\pi^2)}{96}.
}
\]
The apparent \(P_0/g\) family collapses algebraically because
\[
\boxed{\nu(p^3+q^3)=(1-3\nu)(p^2q+pq^2).}
\]
So the full generic-frame static completion is uniquely
\[
\boxed{
\Delta L_{3,\rm ct}^{\rm scalar}
=
\frac{G^4pq}{96r^4}
\bigl(408\nu^2+1232\nu-2080+63\pi^2\bigr)(p^2q+pq^2).
}
\]
Its Hamiltonian image is simply
\[
\Delta H_{3,\rm ct}^{\rm scalar}=-\Delta L_{3,\rm ct}^{\rm scalar}.
\]

### 4.9 Pure-kinetic collapse to the free COM Hamiltonian

The generic-frame free relativistic two-body Lagrangian contributes
\[
\boxed{L_{3,\rm free}^{\rm gen}=\frac{5}{128}(m_Av_A^8+m_Bv_B^8).}
\]
The exact free COM Hamiltonian has
\[
\boxed{
h_1^{\rm free}=\frac{-5+35\nu-70\nu^2+35\nu^3}{128},
\qquad
h_2=h_3=h_4=h_5=0.
}
\]
So the apparent leftover COM pure-kinetic term is fixed automatically:
\[
\boxed{\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16}.}
\]
It is not a new throat-response datum.

### 4.10 Final conservative 3PN split

The exact comparable-mass 3PN residual in COM form is
\[
\boxed{
\Delta \hat L_3^{\rm GR}
=
\underbrace{\Delta l_1v^8}_{\text{compiler image of the free COM Hamiltonian}}
+
\underbrace{L_{P_2}^{\rm mid}}_{\text{exact grouped real }P_2\text{ closure}}
+
\underbrace{\Delta l_{15}^{(g)}U^4}_{\text{unique geometry completion}}.
}
\]

### 4.11 Main theorem

Within the same declared closure hierarchy as the 1PN and 2PN papers,
\[
\mu_{\rho,3}=\frac14,
\qquad
 d_3=-\frac{45}{4},
\qquad
 s_{24}=-\frac1{16},
\]
\[
\Delta H_3=-\Delta L_3,
\qquad
\det M_{\rm mid}=-\frac4{27},
\]
and the exact comparable-mass 3PN residual is fully assigned by the free-Hamiltonian compiler image, the grouped real \(P_2\) middle-block closure, and the unique geometry completion.

### 4.12 What remains open for the 2.5PN bridge

What remains open is no longer algebraic 3PN bookkeeping. The open conservative-to-radiative bridge problem is now: derive the grouped real \(P_2\) constitutive law and the unique geometry completion microscopically from the moving-throat dynamics, then insert that data into the already narrowed 2.5PN outgoing quadrupole normalization problem.

---

## 5) Main proof chain

### 5.1 Carried lower-order conservative backbone

The 3PN paper begins from a lower-order hierarchy that is already frozen: the exact 4D reduction backbone, the full conservative 1PN assembly, the full conservative 2PN ADM match, and the 2.5PN narrowing to the grouped real quadrupole bundle. So 3PN is not solving a fresh EFT; it is extending a already-constrained conservative ladder.

### 5.2 Exact one-body 3PN Schwarzschild gate

The strict one-body target is the exact isotropic Schwarzschild expansion through \(c^{-6}\). This freezes the only acceptable test-mass answer for the 3PN lift.

### 5.3 Minimal one-body repair ledger

Trying to continue the 2PN denominator strategy with only one new cubic denominator parameter fails. The \(U^3v^2\) slot can be repaired, but the \(U^2v^4\) slot remains wrong. This is why 3PN opens the new one-body ledger
\[
(\mu_{\rho,3},d_3,s_{24}).
\]

### 5.4 Natural self/static seed

The exact one-body target is then promoted symmetrically into the natural 3PN self/static two-body seed. This isolates the universal one-body/self/static content before any comparable-mass residual solve begins.

### 5.5 Exact cubic Legendre-transform identity

The 3PN lift from ordinary Lagrangian data to Hamiltonian data is governed by the exact cubic perturbative Legendre identity. Once the 1PN and 2PN ledger is frozen, the genuinely new 3PN piece enters the Hamiltonian in a rigid algebraic way.

### 5.6 Reduced COM 15-slot basis and diagonal compiler

In the isotropic COM basis the cubic Legendre map diagonalizes. This makes the imported GR COM target solvable slot by slot and turns the 3PN problem into a finite algebraic residual rather than a blind variational search.

### 5.7 Imported exact GR COM target and residual extraction

Importing the exact GR 3PN COM Hamiltonian target and subtracting the frozen self/static seed isolates the exact residual vector \((\Delta l_1,\dots,\Delta l_{15})\). This is the first point where the comparable-mass 3PN problem becomes finite and explicit.

### 5.8 Generic-frame residual scaffold

On the generic-frame side, the interaction basis organizes into the 24/17/8/1 scaffold for the \(G/r\), \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks. So the generic-frame problem is also finite from the start.

### 5.9 COM projection, middle-block obstruction, and seed repair

Projecting that scaffold to COM reveals a small \(\nu^3\) obstruction in the middle blocks if one keeps only the naive natural seed split. A minimal \(\nu\)-dressed seed repair removes this mismatch without changing the strict one-body limit. The obstruction is therefore bookkeeping/placement, not a failure of the full scaffold.

### 5.10 COM-null ideal and contact/gauge orbit

After seed repair, the remaining unfixed COM-blind directions form a clean COM-null ideal. The actual ordinary contact/gauge orbit is smaller: 11-dimensional inside a larger 27-dimensional COM-blind family. So COM blindness and ordinary gauge freedom are not identical.

### 5.11 Naive ordinary COM reduction fails

A naive COM reduction of the imported generic-frame ordinary 3PN target does **not** reproduce the exact COM ordinary target. This is the same subtlety known from standard 3PN literature. So naive ordinary COM reduction cannot be the final generic-frame quotient principle.

### 5.12 Hamiltonian-first lift fixes the issue

If the generic-frame ordinary target is lifted through the exact cubic Legendre map **before** COM reduction, the reduced Hamiltonian reproduces the exact standard COM 3PN Hamiltonian. This makes the Hamiltonian-first route the authoritative one.

### 5.13 Fixed-chart uniqueness theorem

Once the lower-order ledger is frozen, the generic-frame Hamiltonian compiler reduces the whole 50-slot residual to
\[
\Delta H_3=-\Delta L_3.
\]
So in the fixed ADM chart the generic-frame representative is unique. At that point the algebraic lift is no longer the bottleneck; the remaining work moves into the grouped real \(P_2\) module, the static scalar/geometry lane, and the pure-kinetic interpretation.

---

## 6) The grouped real \(P_2\) conservative module

### 6.1 Why grouped real \(P_{20}\oplus P_{21}\oplus P_{22}\) is the right 3PN front end

The 2.5PN notes already reduced the surviving universal route to the orbital/worldtube STF quadrupole branch. In grouped real language, that means the right conservative front end is the grouped real
\[
P_{20}\oplus P_{21}\oplus P_{22}
\]
bundle, not an unfocused collection of generic higher-order coefficients.

### 6.2 Grouped conservative data and isotropy variables

At 3PN the grouped conservative data are
\[
u_2^{(20)},\qquad u_2^{(21)},\qquad u_2^{(22)}.
\]
Equivalently,
\[
\bar u_2=\frac{u_2^{(20)}+2u_2^{(21)}+2u_2^{(22)}}{5},
\qquad
 a_2=\frac{2u_2^{(20)}-u_2^{(21)}-u_2^{(22)}}{10},
\qquad
 b_2=\frac{u_2^{(21)}-u_2^{(22)}}{2}.
\]
The exact isotropy gate for the minimal branch is therefore
\[
\boxed{a_2=0,\qquad b_2=0.}
\]

### 6.3 Minimal local grouped-\(P_2\) front end and why it fails

The direct time-local \(O(\omega^2)\) grouped front end is
\[
L_{P_2,\omega^2}^{\rm front}
=
\frac{c_2^{(20)}}{3r^2}d^2(3u_\perp^2-U)^2
+
\frac{c_2^{(21)}}{r^2}u_\perp^2(u_\perp^2-d^2-U)^2
+
\frac{c_2^{(22)}}{r^2}d^2u_\perp^4.
\]
It has uniform virial weight five, so the unique local isotropic one-step demotion is multiplication by one power of \(U^{-1}\sim r\). After that demotion the image lands only in the 9 interaction slots \((l_6,\dots,l_{14})\), has rank three, and obeys six exact relations:
\[
2l_6+2l_7+l_8=0,
\qquad
l_9=l_6+l_7,
\qquad
l_{10}=-2l_6,
\]
\[
l_{11}+l_{12}=2l_6,
\qquad
l_{13}=l_6,
\qquad
l_{14}=-\frac16 l_{11}.
\]
The exact 3PN GR target violates all six. So the minimal local grouped lift is too rigid.

### 6.4 Richer local grouped-\(P_2\) family

The grouped ontology itself is large enough; only the minimal constitutive choice was too small. Define
\[
C_{20}^2=\frac16(3d^2-v^2-2U)^2,
\qquad
C_{21}^2=2d^2u_\perp^2,
\qquad
C_{22}^2=\frac12u_\perp^4,
\]
\[
T_{20}=\frac13U d^2(3u_\perp^2-U)^2,
\qquad
T_{21}=Uu_\perp^2(u_\perp^2-d^2-U)^2,
\qquad
T_{22}=U d^2u_\perp^4,
\]
\[
S_A=U^2 C_A^2,
\qquad
V_A=\frac{v^2}{U}S_A,
\qquad A\in\{20,21,22\}.
\]
Then the exact 9-family middle-block compiler built from
\[
(T_{20},T_{21},T_{22},S_{20},S_{21},S_{22},V_{20},V_{21},V_{22})
\]
has
\[
\boxed{\det M_{\rm mid}=-\frac4{27}.}
\]

### 6.5 Exact GR grouped coefficients

Writing the richer grouped lift as
\[
L_{P_2}^{\rm mid}
=
\sum_{A\in\{20,21,22\}}\bigl(\lambda_A T_A+\sigma_A S_A+\tau_A V_A\bigr),
\]
the exact fit to the solved GR middle block is
\[
\lambda_{20}^{\rm GR}=\frac{\nu\,(-276\nu^2-3154\nu+3\pi^2+2672)}{32},
\]
\[
\lambda_{21}^{\rm GR}=\frac{\nu\,(2568\nu^2+2236\nu-3044-3\pi^2)}{192},
\]
\[
\lambda_{22}^{\rm GR}=\frac{\nu\,(7650\nu^2+32858\nu-30136-33\pi^2)}{96},
\]
\[
\sigma_{20}^{\rm GR}=\frac{\nu\,(-102\nu^2-308\nu+293)}{16},
\]
\[
\sigma_{21}^{\rm GR}=\frac{\nu\,(-4188\nu^2-23778\nu+21\pi^2+22006)}{192},
\]
\[
\sigma_{22}^{\rm GR}=\frac{\nu\,(3906\nu^2+2478\nu-2791-3\pi^2)}{48},
\]
\[
\tau_{20}^{\rm GR}=\frac{3\nu\,(-154\nu^2-87\nu+38)}{32},
\]
\[
\tau_{21}^{\rm GR}=\frac{\nu\,(-6672\nu^2-5824\nu+3\pi^2+4412)}{384},
\]
\[
\tau_{22}^{\rm GR}=\frac{\nu\,(-2790\nu^2-3367\nu+3\pi^2+3386)}{96}.
\]
With this coefficient set, the richer grouped lift reproduces the exact nine-slot GR middle block
\[
(l_6,\dots,l_{14})=(\Delta l_6,\dots,\Delta l_{14})
\]
with no middle-block residual left.

### 6.6 Structural corollaries

Two corollaries matter.

First, every element in \((T_A,S_A,V_A)\) carries at least one power of \(U\), so the grouped module still forces
\[
l_1=l_2=l_3=l_4=l_5=0
\]
identically.

Second, the grouped compiler predicts a specific static companion:
\[
\boxed{l_{15}^{(P_2\,\rm pred)}=l_{10}+l_{11}+l_{12}+2(l_6+l_7+l_8+l_9).}
\]
That is the formula that produces the explicit static gap quoted in Section 4.8.

---

## 7) Static scalar/geometry completion

### 7.1 Residual static slot after grouped closure

After exact grouped middle-block closure, the 3PN COM residual splits as
\[
\Delta \hat L_3^{\rm GR}
=
\Delta l_1v^8
+
L_{P_2}^{\rm mid}[\lambda_A^{\rm GR},\sigma_A^{\rm GR},\tau_A^{\rm GR}]
+
\Delta l_{15}^{(0/g)}U^4.
\]
So the grouped real \(P_2\) module has done all it can do; the remaining conservative scalar job is the one-dimensional static remainder.

### 7.2 U-block no-go inside the old constant-coefficient generic scaffold

In the old generic-frame interaction scaffold, the only static basis polynomial is
\[
U_1=p^2q+pq^2
\]
inside the common prefactor \(G^4pq/r^4\). Its COM image is always proportional to
\[
\nu U^4,
\]
linear in \(\nu\). But the exact static remainder contains \(\nu^2\) and \(\nu^3\) after the overall \(\nu\) is factored out. So the old constant-coefficient static block is too small to carry the actual 3PN scalar gap.

### 7.3 Apparent \(P_0\) versus geometry ambiguity

Once one allows body-local scalar pieces as well as pair/geometry pieces, the static completion initially looks like a two-family problem: a body-local \(P_0\)-type family and a pair/geometry family. Algebraically, however, this ambiguity is fake.

### 7.4 Sigma-collapse and unique geometry completion

The exact mass identity
\[
\nu(p^3+q^3)=(1-3\nu)(p^2q+pq^2)
\]
collapses the apparent two-family static ambiguity to one unique pair/geometry representative. The final static completion is therefore exactly the geometry-side counterterm quoted earlier:
\[
\Delta L_{3,\rm ct}^{\rm scalar}
=
\frac{G^4pq}{96r^4}
\bigl(408\nu^2+1232\nu-2080+63\pi^2\bigr)(p^2q+pq^2).
\]

### 7.5 Hamiltonian image of the static completion

Because the fixed-chart generic-frame compiler is simply
\[
\Delta H_3=-\Delta L_3,
\]
the Hamiltonian image of the static completion is immediately
\[
\Delta H_{3,\rm ct}^{\rm scalar}=-\Delta L_{3,\rm ct}^{\rm scalar}.
\]
There is no further algebraic freedom to resolve.

---

## 8) Pure-kinetic sector and the free two-body Hamiltonian

### 8.1 Why the apparent pure-kinetic residual is not a new throat datum

After grouped middle-block closure and unique static completion, the only apparently separate COM slot is
\[
\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16},
\qquad
\Delta l_2=\Delta l_3=\Delta l_4=\Delta l_5=0.
\]
The exact kinetic audit shows that this does **not** require a new 3PN kinetic throat-response datum.

### 8.2 Generic-frame free kinematics

At generic-frame level, the exact free ordinary relativistic two-body Lagrangian is simply the sum of the two free-body square roots. Through 3PN it gives
\[
L_{3,\rm free}^{\rm gen}=\frac{5}{128}(m_Av_A^8+m_Bv_B^8).
\]
There is no mixed comparable-mass pure-kinetic interaction monomial.

### 8.3 Exact free COM Hamiltonian target

The exact free relativistic two-body COM Hamiltonian is
\[
H_{\rm free}^{\rm COM}
=
\sqrt{m_A^2c^4+P^2c^2}+
\sqrt{m_B^2c^4+P^2c^2},
\]
and its reduced 3PN coefficient is
\[
h_1^{\rm free}=\frac{-5+35\nu-70\nu^2+35\nu^3}{128},
\qquad h_2=h_3=h_4=h_5=0.
\]

### 8.4 Collapse of the ordinary pure-kinetic residual

Applying the exact COM compiler gives
\[
l_1=F_1(\nu)-h_1,
\qquad
F_1(\nu)=\frac{3}{16}\nu-\frac{21}{16}\nu^2+\frac94\nu^3.
\]
Subtracting the true free COM Hamiltonian target from the naive free ordinary seed gives
\[
\Delta l_1=h_1^{\rm seed}-h_1^{\rm free}=\frac{3\nu(3\nu-1)(4\nu-1)}{16}.
\]
So the whole pure-kinetic leftover is exactly the ordinary-Lagrangian compensation required to reproduce the universal free relativistic two-body COM Hamiltonian.

---

## 9) Final conservative 3PN theorem

### 9.1 Final exact split in COM form

Within the declared hierarchy, the exact comparable-mass 3PN coefficient in the isotropic COM chart is
\[
\boxed{
\Delta \hat L_3^{\rm GR}
=
\underbrace{\Delta l_1v^8}_{\text{compiler image of the free COM Hamiltonian}}
+
\underbrace{L_{P_2}^{\rm mid}}_{\text{exact grouped real }P_2\text{ middle-block closure}}
+
\underbrace{\Delta l_{15}^{(g)}U^4}_{\text{unique geometry completion}}.
}
\]
The two scalar coefficients worth keeping visible are
\[
\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16},
\qquad
\Delta l_{15}^{(g)}=\frac{\nu(408\nu^2+1232\nu-2080+63\pi^2)}{96}.
\]

### 9.2 Final exact split in generic-frame form

Once the lower-order ledger is frozen, the fixed ADM chart obeys
\[
\Delta H_3=-\Delta L_3.
\]
So the earlier COM-blind null directions and the smaller ordinary contact/gauge orbit are intermediate bookkeeping structures only. They do not survive as true algebraic freedom after the full generic-frame Hamiltonian target is imposed.

### 9.3 Equality to the standard 3PN ADM target

The Hamiltonian-first lift recovers the standard 3PN ADM-type target in exactly the same sense that the full conservative 1PN paper closes to EIH and the full conservative 2PN paper closes to ADM: exact final algebra **within the declared hierarchy**.

### 9.4 Fixed versus open quantities

What should now be treated as fixed carry-forward 3PN data is:
\[
\mu_{\rho,3}=\frac14,
\qquad
 d_3=-\frac{45}{4},
\qquad
 s_{24}=-\frac1{16},
\]
\[
\det M_{\rm mid}=-\frac4{27},
\qquad
\Delta l_{15}^{(g)}=\frac{\nu(408\nu^2+1232\nu-2080+63\pi^2)}{96},
\qquad
\Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16}.
\]
What remains open is no longer algebraic coefficient matching. It is the deeper moving-throat derivation of the grouped real \(P_2\) constitutive map and the geometry-side completion, followed by their insertion into the already narrowed 2.5PN quadrupole bridge.

---

## 10) Interface to the 2.5PN program

### 10.1 Why 3PN is the first conservative grouped-\(P_2\) gate

The 2.5PN notes already reduced the live higher-order universal point-particle route to the grouped real
\[
P_{20}\oplus P_{21}\oplus P_{22}
\]
quadrupole bundle. The present 3PN result is the first conservative theorem-level test of exactly that same branch.

### 10.2 The \(a_2=b_2=0\) isotropy test and grouped low-frequency data

At 3PN the conservative grouped data are still only the three \(O(\omega^2)\) numbers
\[
u_2^{(20)},\qquad u_2^{(21)},\qquad u_2^{(22)}.
\]
So the first grouped conservative branch test remains
\[
\boxed{a_2=0,\qquad b_2=0.}
\]
If it passes, the 3PN conservative data collapse to a single isotropic trace datum \(\bar u_2\). If it fails, the minimal isotropic quadrupole branch is already dead before any outgoing normalization work is attempted.

### 10.3 What the present 3PN result fixes for the minimal outgoing \(l=2\) branch

The main positive interface result is that the grouped real \(P_2\) ontology survives its first honest conservative stress test. The richer local grouped-\(P_2\) compiler closes the full 9-slot middle block exactly, so the throat side does not need a qualitatively new non-\(P_2\) conservative carrier for that sector. The remaining static slot is uniquely geometry-side, and the pure-kinetic slot is only a compiler image of the free COM Hamiltonian, so neither competes with the grouped quadrupole module for the middle-block burden.

### 10.4 What the full moving-throat PDE still has to supply

The remaining 2.5PN theorem gap is therefore not “derive 3PN somehow.” It is narrower. The moving-throat PDE still has to supply:

1. the actual low-frequency constitutive origin of the grouped real \(P_2\) coefficients;
2. the isotropy verdict \(a_2=b_2=0\) on the true branch;
3. and then the final passive/outgoing quadrupole normalization on that branch.

So the 3PN theorem should be read as removing the conservative algebraic uncertainty **beneath** the 2.5PN outgoing-normalization gap, not as closing that gap itself.

---

## 11) Verification / reproducibility ledger

### 11.1 Role of the 3PN referee suite

The 3PN verification archive plays the same role as the earlier referee suites: it is a symbolic reproducibility ledger, not a black-box derivation engine. It verifies exact identities, compiler maps, residual extractions, grouped-\(P_2\) closure, sigma-collapse, and the final theorem package **inside the declared closure hierarchy**.

### 11.2 Master script

The main verification artifact is

- `3pn_referee_master_sympy_audit.py`

It is a standalone SymPy referee script that vendors the source text of every stage audit, replays them in isolated namespaces, and stops immediately on the first failed identity. It uses only the Python standard library plus SymPy.

### 11.3 Embedded stage audits carried inside the master script

The master script vendors and replays the following stage audits:

1. `3pn_onebody_audit.py`
2. `3pn_grouped_p2_audit.py`
3. `3pn_comparable_mass_audit.py`
4. `3pn_com_linear_map_audit.py`
5. `3pn_com_gr_target_audit.py`
6. `3pn_generic_frame_com_projection_audit.py`
7. `3pn_seed_repair_and_com_null_ideal_audit.py`
8. `3pn_contact_generator_and_gauge_orbit_audit.py`
9. `3pn_generic_frame_target_import_audit.py`
10. `3pn_hamiltonian_level_lift_audit.py`
11. `3pn_generic_frame_hamiltonian_compiler_audit.py`
12. `3pn_grouped_p2_target_pack_audit.py`
13. `3pn_grouped_p2_middle_block_test_audit.py`
14. `3pn_grouped_p2_richer_lift_audit.py`
15. `3pn_static_p0_geometry_counterterm_audit.py`
16. `3pn_sigma_collapse_and_unique_geometry_completion_audit.py`
17. `3pn_pure_kinetic_collapse_audit.py`
18. `3pn_conservative_theorem_audit.py`

### 11.4 What the suite verifies

The suite verifies all of the following.

- the exact one-body 3PN Schwarzschild gate,
- the minimal one-body repair ledger,
- the cubic Legendre-transform identity,
- the exact COM ordinary/Hamiltonian compiler,
- the exact imported GR COM target and residual extraction,
- the generic-frame projection / seed repair / COM-null structure,
- the Hamiltonian-first fixed-chart uniqueness theorem,
- the grouped real \(P_2\) target pack,
- the minimal-lift no-go,
- the richer grouped middle-block closure with \(\det M_{\rm mid}=-4/27\),
- the static sigma-collapse and unique geometry completion,
- the pure-kinetic collapse to the free COM Hamiltonian,
- and the final conservative 3PN theorem ledger.

### 11.5 What the suite does **not** verify

The suite does **not** verify a fully solved moving-throat PDE, and it does not replace the future microscopic derivation of the grouped real \(P_2\) constitutive law or the geometry-side completion. The fixed-versus-open distinction is the same one used by the earlier conservative papers: exact algebra once the hierarchy is declared, but no false claim that the full moving-throat dynamics is already solved.

---

## 12) Mental checklist for future work

If a future session needs the shortest durable memory of this paper, the non-negotiable carry-forward points are:

1. **The one-body 3PN gate is fixed** by
   \[
   \mu_{\rho,3}=\frac14,
   \qquad d_3=-\frac{45}{4},
   \qquad s_{24}=-\frac1{16}.
   \]
2. **The exact COM residual is known** and the generic-frame fixed-chart compiler is simply
   \[
   \Delta H_3=-\Delta L_3.
   \]
3. **Minimal grouped \(P_2\) fails, richer grouped \(P_2\) succeeds.** The necessary conservative family is
   \[
   (T_A,S_A,V_A),
   \qquad
   \det M_{\rm mid}=-\frac4{27}.
   \]
4. **The grouped module closes the whole 9-slot middle block, but not the kinetic or static leftovers.**
5. **The static leftover is uniquely geometry-side,** with
   \[
   \Delta l_{15}^{(g)}=\frac{\nu(408\nu^2+1232\nu-2080+63\pi^2)}{96}.
   \]
6. **The pure-kinetic leftover is not new throat physics.** It is only the ordinary-Lagrangian counterimage of the universal free COM Hamiltonian:
   \[
   \Delta l_1=\frac{3\nu(3\nu-1)(4\nu-1)}{16}.
   \]
7. **The remaining open problem is now dynamical, not algebraic:** derive the grouped real \(P_2\) constitutive law and the geometry completion microscopically, then close the narrow 2.5PN passive/outgoing normalization gap on the true moving-throat branch.
