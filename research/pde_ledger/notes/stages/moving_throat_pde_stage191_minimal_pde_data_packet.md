# Moving-Throat PDE — Stage 191: Minimal PDE Data Packet, Exact Branch/Orbit Residuals, and the Home-Stretch Theorem

## Status

**Exact within the carried grouped real `P2` / coherent-branch / outgoing-compiler hierarchy already fixed in Stages 239–241.**

This stage does **not** add a new microscopic closure law.
It packages the existing reduced theorem chain into the smallest exact finite data packet that the completed moving-throat PDE still has to supply.

---

## Purpose

Stage 239 completed the first-order **branch-observable** front end in the direct observables
\[
R_{\rm tr},
\qquad
\mathfrak N_*:=\mathcal T^2 R_{\rm tr}^{B_*},
\qquad
\epsilon_\eta,
\]
and showed that their logarithmic drifts are an exact first-order compiler for
\[
(\Theta_1,\Xi_1,\mathcal R_1).
\]

Stage 240 carried the same branch-observable data into the exact **transfer-shape / grouped-response / outgoing-prefactor** compiler,
so that the isotropic grouped bundle is described by the low-frequency coefficients
\[
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
\qquad
A\in\{20,21,22\},
\]
together with the source-map factor `\(m_{\hat 0}\)`.

Stage 241 then separated the **direct transfer-shape defect** from the **selected-branch dressing residual**, proved the support-blindness theorem for the direct transfer shape, and promoted the main weak-axisymmetric no-go filters. After that stage, the remaining gap was no longer “what microscopic channel is missing?” but rather:

> what is the smallest exact finite packet of branch data the completed PDE must return so that the final reduced verdict is automatic?

This stage answers that question.

The main outputs are:

1. the exact **Packet A** compiler from grouped-lane low-frequency data to the final branch residual packet
   \[
   \Delta_{\rm branch},
   \]
2. the exact **Packet B** interconversion between the three equivalent orbit-lock packets
   \[
   (m_T,m_K,m_\mu),
   \qquad
   (\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta),
   \qquad
   (q_{\rm tr},q_{\rm nt},q_\eta),
   \]
3. the exact reduced closure criterion
   \[
   \Delta_{\rm branch}=0,
   \qquad
   \Delta_{\rm orbit}=0,
   \]
4. and the sharp statement that **everything else in the reduced 2.5PN / 4PN / 5PN endgame is downstream algebra of those two packets.**

So this stage is the finite-packet endpoint of the moving-throat PDE derivation program in its present hierarchy.

---

## 1. Carry-forward grouped-lane bundle data

### 1.1 Packet A: the grouped low-frequency PDE bundle

For each grouped real `P2` lane
\[
A\in\{20,21,22\},
\]
carry forward the exact conservative wall/worldtube operator and outgoing-transfer expansions
\[
D_A^{(\mathrm{cons})}(\omega)
=
D_{A0}+D_{A2}\omega^2+D_{A4}\omega^4+O(\omega^6),
\]
\[
N_A(\omega)
=
N_{A0}+N_{A2}\omega^2+N_{A4}\omega^4+O(\omega^6).
\]
These coefficients already absorb everything computed in the previous stages:

- geometry-only wall data,
- BdG support self-energy,
- conservative localized-Maxwell / mixed-sector self-energy,
- and the static part of the outgoing transfer block.

The exact finite grouped bundle packet is therefore
\[
\boxed{
\mathcal P_A
:=
\Bigl(
(D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4})_{A=20,21,22},
\ m_{\hat 0}
\Bigr).
}
\]

This is the smallest branch-side packet that still retains all low-frequency information needed by the reduced endgame compiler.

### 1.2 Exact normalized grouped response and prefactor moments

From `\(\mathcal P_A\)` define the exact grouped response moments
\[
\boxed{
\nu_2^{(A)}:=-\frac{D_{A2}}{D_{A0}},}
\qquad
\boxed{
\nu_4^{(A)}:=\frac{D_{A2}^2-D_{A0}D_{A4}}{D_{A0}^2}.}
\]
To avoid collision with the bulk velocity symbol `v`, this stage writes these grouped response moments as `\(\nu_2^{(A)},\nu_4^{(A)}\)` rather than `u_2^{(A)},u_4^{(A)}`. They are the same objects used in the later endgame notes.

Likewise define the outgoing prefactor moments
\[
\boxed{P_0^{(A)}:=\frac{N_{A0}}{D_{A0}},}
\]
\[
\boxed{P_2^{(A)}:=\frac{D_{A0}N_{A2}-2D_{A2}N_{A0}}{D_{A0}^2},}
\]
\[
\boxed{
P_4^{(A)}:=
\frac{D_{A0}^2N_{A4}-2D_{A0}(D_{A2}N_{A2}+D_{A4}N_{A0})+3D_{A2}^2N_{A0}}
{D_{A0}^3}.
}
\]

These formulas are exact compiler identities. They do not introduce any new closure choice.

### 1.3 Exact grouped trace/anomaly decomposition

For any grouped triple
\[
x=(x_{20},x_{21},x_{22}),
\]
define the weighted grouped trace/anomaly variables
\[
\boxed{
\bar x:=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\qquad
 a_x:=\frac{2x_{20}-x_{21}-x_{22}}{10},
\qquad
 b_x:=\frac{x_{21}-x_{22}}{2}.
}
\]
The exact inverse map is
\[
\boxed{x_{20}=\bar x+4a_x,}
\qquad
\boxed{x_{21}=\bar x-a_x+b_x,}
\qquad
\boxed{x_{22}=\bar x-a_x-b_x.}
\]

Apply this decomposition to the grouped response and prefactor moments:
\[
(\bar\nu_2,a_2,b_2),
\qquad
(\bar\nu_4,a_4,b_4),
\qquad
(\bar P_0,a_{P_0},b_{P_0}).
\]

Then the exact grouped isotropy gate is simply
\[
a_2=b_2=a_4=b_4=a_{P_0}=b_{P_0}=0.
\]

### 1.4 Exact one-pole and normalization defects

The one-pole conservative defect is
\[
\boxed{\Delta_{\rm pole}:=\bar\nu_4-4\bar\nu_2^{\,2}.}
\]
It measures failure of the isotropic grouped response to land on the one-pole conservative identity.

The normalization defect is
\[
\boxed{
\Delta_{\rm norm}
:=
 m_{\hat 0}^{\,2}\bar P_0-
 \frac{54Gc_s^5}{5a^5c^5}.
}
\]
It is exactly the invariant normalization gap already isolated by the outgoing quadrupole bridge.

### 1.5 The exact final branch residual packet

Collecting the previous quantities gives the exact branch-side verdict vector
\[
\boxed{
\Delta_{\rm branch}
:=
\bigl(
 a_2,\ b_2,\ a_4,\ b_4,\ a_{P_0},\ b_{P_0},\ \Delta_{\rm pole},\ \Delta_{\rm norm}
\bigr).
}
\]

This is the exact finite packet whose zero-set encodes the grouped-lane isotropy, one-pole, and outgoing-normalization demands.

### 1.6 Why `P_2` and `P_4` are compiled but not retained in `\Delta_{\rm branch}`

Packet A contains `\(N_{A2},N_{A4}\)` because they are part of the exact low-frequency bundle and because they compile the even outgoing prefactor moments `\(P_2^{(A)},P_4^{(A)}\)`.

But the **minimal reduced verdict packet** does not need extra residual entries built from `P_2` and `P_4`.
Within the present hierarchy, the leading odd quadrupole coefficient depends only on the isotropic static prefactor `\(\bar P_0\)`, while the grouped conservative side is already monitored by
\[
(a_2,b_2,a_4,b_4,\Delta_{\rm pole}).
\]
So `P_2` and `P_4` remain exact derived diagnostics, but they are not additional independent verdict coordinates in the minimal home-stretch packet.

---

## 2. Carry-forward orbit-lock data

### 2.1 Packet B: three exact representations of the same finite orbit verdict

The orbit side is already known to be exactly representable by any one of three equivalent packets.

#### (i) Residual mismatch ratios
\[
\boxed{(m_T,m_K,m_\mu).}
\]

#### (ii) Invariant ratios
To avoid collision with the exact Stage 239 tracking factor `\(R_{\rm tr}\)`, this note writes the finite orbit/invariant ratios as
\[
\boxed{(\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta).}
\]
The later endgame notes denote the same packet by `(R_tr,R_nt,R_eta)`.

#### (iii) Quotient coordinates
\[
\boxed{(q_{\rm tr},q_{\rm nt},q_\eta).}
\]

So the exact orbit-side packet is
\[
\boxed{\mathcal P_B:=\text{any one of the three equivalent packets above}.}
\]

### 2.2 Exact interconversion formulas

Let `\(\chi_{0,*}\)` and `\(F_*\)` be the carried coherent-branch coefficients.
Then the exact packet interconversion laws are
\[
\boxed{\mathfrak R_{\rm tr}=m_T^{1+\chi_{0,*}},}
\qquad
\boxed{\mathfrak R_{\rm nt}=\frac{m_\mu}{m_K m_T^{F_*}},}
\qquad
\boxed{\mathfrak R_\eta=\frac{1}{m_K}.}
\]
Taking logarithms gives
\[
\boxed{q_{\rm tr}=\ln \mathfrak R_{\rm tr},}
\qquad
\boxed{q_{\rm nt}=\ln \mathfrak R_{\rm nt},}
\qquad
\boxed{q_\eta=\ln \mathfrak R_\eta.}
\]
Conversely,
\[
\boxed{m_T=e^{q_{\rm tr}/(1+\chi_{0,*})},}
\qquad
\boxed{m_K=e^{-q_\eta},}
\qquad
\boxed{m_\mu=e^{q_{\rm nt}-q_\eta+F_* q_{\rm tr}/(1+\chi_{0,*})}.}
\]

### 2.3 Exact orbit residual packet

The sharpest orbit-lock verdict is therefore the quotient-coordinate packet
\[
\boxed{
\Delta_{\rm orbit}:=
(q_{\rm tr},q_{\rm nt},q_\eta).
}
\]
Its zero-set is exactly orbit lock:
\[
\boxed{
\Delta_{\rm orbit}=0
\iff
m_T=m_K=m_\mu=1
\iff
\mathfrak R_{\rm tr}=\mathfrak R_{\rm nt}=\mathfrak R_\eta=1.
}
\]

### 2.4 Relation to the Stage 239 quotient directions

At the infinitesimal level, `\(\Delta_{\rm orbit}\)` is just the finite upgrade of the same three quotient directions whose tangent compiler was fixed in Stages 238–239.
So Packet B is not a new object. It is the exact finite orbit-lock version of the already-frozen three-coordinate quotient structure.

---

## 3. Exact home-stretch theorem

### 3.1 Statement

\[
\boxed{\textbf{Theorem (Stage 242 home-stretch theorem).}}
\]

Within the carried grouped real `P_2` / coherent-branch / outgoing-prefactor hierarchy of Stages 239–241,
all reduced endgame diagnostics needed for the final GR-compatible test are exact compiler outputs of the two finite packets
\[
\mathcal P_A
\qquad\text{and}\qquad
\mathcal P_B.
\]
In particular,
\[
\boxed{
\text{reduced closure is complete}
\iff
\Delta_{\rm branch}=0
\quad\text{and}\quad
\Delta_{\rm orbit}=0.
}
\]

### 3.2 Meaning

This theorem says that the moving-throat PDE no longer has to “show 5PN somehow,” or “produce the anomaly correction somehow,” or “demonstrate 2.5PN/4PN closure somehow.”
It only has to supply the actual branch values needed to evaluate two exact finite residual vectors:
\[
\Delta_{\rm branch},
\qquad
\Delta_{\rm orbit}.
\]
Everything else is compiler algebra already fixed by the previous stages.

### 3.3 Decision tree

The reduced verdict is now immediate.

1. **Compile Packet A** into `\(\Delta_{\rm branch}\)`.
2. **Compile Packet B** into `\(\Delta_{\rm orbit}\)`.
3. If
   \[
   \Delta_{\rm branch}\neq 0,
   \]
   then the branch fails the reduced GR-compatible grouped test.
4. If
   \[
   \Delta_{\rm branch}=0
   \quad\text{but}\quad
   \Delta_{\rm orbit}\neq 0,
   \]
   then the branch is isotropic and normalized on the grouped side but is still off the exact similarity orbit.
5. Only if
   \[
   \Delta_{\rm branch}=0
   \quad\text{and}\quad
   \Delta_{\rm orbit}=0
   \]
   is the reduced closure complete inside the present hierarchy.

---

## 4. What is now settled, and what the completed PDE still has to compute

### 4.1 Settled inside the present hierarchy

The following are no longer active algebraic bottlenecks.

1. The first-order branch-observable front end is fixed exactly by Stage 239.
2. The transfer-shape / outgoing-prefactor compiler is fixed exactly by Stage 240.
3. The direct-defect / dressing split and the main weak-axisymmetric no-go filters are fixed exactly by Stage 241.
4. The final finite residual packets `\(\Delta_{\rm branch}\)` and `\(\Delta_{\rm orbit}\)` are now explicit.

So the remaining gap is no longer compression, notation, or algebra.

### 4.2 What the completed moving-throat PDE still has to return

The completed PDE still has to compute the actual branch values of:

- the grouped bundle coefficients
  \[
  (D_{A0},D_{A2},D_{A4},N_{A0},N_{A2},N_{A4}),
  \qquad A\in\{20,21,22\},
  \]
- the source-map factor `\(m_{\hat 0}\)`,
- and one exact form of the orbit-lock packet
  \[
  (m_T,m_K,m_\mu),
  \quad\text{or}\quad
  (\mathfrak R_{\rm tr},\mathfrak R_{\rm nt},\mathfrak R_\eta),
  \quad\text{or}\quad
  (q_{\rm tr},q_{\rm nt},q_\eta).
  \]

That is the entire remaining home-stretch theorem gap in the current hierarchy.

---

## 5. Script-backed status

The accompanying SymPy audit verifies:

- the exact single-lane response and prefactor compilers from `\((D_0,D_2,D_4,N_0,N_2,N_4)\)`,
- the exact grouped trace/anomaly inverse formulas and weighted-projector identities,
- the exact compilation of `\(\Delta_{\rm branch}\)` from Packet A,
- isotropic collapse of the grouped anisotropy coordinates,
- vanishing of `\(\Delta_{\rm branch}\)` on the isotropic one-pole normalized test branch,
- the exact orbit-packet interconversion formulas,
- and the zero-set equivalence of the three orbit-lock packet representations.

Supporting files:

- `moving_throat_pde_stage242_minimal_pde_data_packet_sympy_audit.py`
- `moving_throat_pde_stage242_minimal_pde_data_packet_sympy_audit_output.txt`
