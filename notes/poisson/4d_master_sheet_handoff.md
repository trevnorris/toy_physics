# Paper 7 Hard‑Mode 4D — Master Sheet + Handoff (Stiff Water v1)

**Purpose:** single, self‑contained “carry forward” sheet: frozen modeling decisions, full equation stack, what is already verified/derived, and the concrete next steps needed to finish the 4D unification (Poisson/EM/leakage) and move toward a paper.

---

## 0) Arena, fields, and state

**Bulk coordinates:**
- Spatial: \(\mathbf X=(x,y,z,w)\in\mathbb R^4\)
- Time: \(t\)
- Operators: \(\nabla_4=(\partial_x,\partial_y,\partial_z,\partial_w)\), \(\nabla_4^2=\partial_x^2+\partial_y^2+\partial_z^2+\partial_w^2\)
- Brane: operationally around \(w=0\)

**Master state (what evolves):**
\[
\boxed{\ \big(\psi(\mathbf X,t),\ A_M(\mathbf X,t),\ a(t),\ L(t)\big)\ }
\]
where \(\psi\) is the superfluid/order parameter, \(A_M\) is a genuine 4+1D gauge field (needed for an EM sector), and \((a,L)\) are reduced geometry DOFs for the throat.

---

## 1) Frozen “Stiff Water” configuration (decisions)

This is the *committed* baseline used for derivations and simulations.

### 1.1 Confinement potential family
**Frozen:** Family‑1 “modulated brane trap” (geometry‑as‑potential).

**Intuition:** far from the defect the brane is stiff (strong \(w\) confinement). Inside the throat corridor, confinement is relaxed.

### 1.2 Geometry energy (closing force)
**Frozen:** surface term (contractile)\; \(E_{\rm geom}\supset \sigma\,\mathcal A\) (optionally add small curvature regularizer later only if numerically required).

### 1.3 Closure law for geometry
**Frozen:** dynamical closure (“breathing mode”): \(a(t),L(t)\) respond to changing stresses/energies.

### 1.4 System constraint / closure
**Frozen:** open system with **fixed chemical potential** \(\mu\) (bulk as reservoir).

### 1.5 Brane projection weight
**Frozen:** weighted projection with brane ground‑state weight
\[
\mathcal P_W[f](x,y,z,t)=\int_{-\infty}^{\infty}W(w)f(x,y,z,w,t)\,dw,\qquad W=|\chi_0(w)|^2.
\]

### 1.6 Measurement region
**Frozen:** brane 2‑sphere \(\Gamma=S^2\) at \(R_3=r_{\rm port}\), \(w=0\), outward normal \(\hat n=\hat R_3\).

### 1.7 Flux definition
**Frozen:** **multi‑output**: mouth flux **and** bulk leakage flux.

### 1.8 Drive protocol
**Frozen:** wave/vorticity injection ("spark") — drive the throat region by injecting rotational/wave energy (not by manually moving walls).

### 1.9 Refill physics knob
**Frozen:** **comoving refill** (incoming bulk support travels with the defect → reduced drag/impedance).

---

## 2) Final freezes (Section 9 resolution) — items 10.1–10.5

### 2.1 Charge neutrality
**Frozen:** jellium / background subtraction
\[
J^0(\mathbf X)=q\big(|\psi(\mathbf X)|^2-\rho_0\big).
\]

### 2.2 Gauge localization
**Frozen:** dielectric / kinetic prefactor waveguide
\[
\mathcal L_{\rm EM}= -\frac{1}{4\mu_0}\,Z(w)\,F_{MN}F^{MN},\qquad Z(w)=\exp\!\big(-w^2/\lambda_{\rm conf}^2\big).
\]

### 2.3 Exact potential shapes
**Frozen:**
- Far‑field brane trap: harmonic\; \(V_w(w)=\tfrac12 m\omega_w^2 w^2\).
- Throat envelope: super‑Gaussian order 4\; \(E(r)=\exp(-r^4/a^4)\).

### 2.4 Boundary conditions
**Frozen:** sponge layers (CAP)
\[
V_{\rm sim}(\mathbf X)=V_{\rm conf}(\mathbf X)-i\,\Gamma_{\rm sponge}(\mathbf X).
\]

### 2.5 Shape field
**Frozen for single‑particle studies:** no shape field yet. Geometry is only \(a(t),L(t)\). Shape‑field \(\xi(\mathbf x)\) reserved for N‑body/interaction studies.

---

## 3) Governing equations (the “master bundle”)

### 3.1 EOS (non‑negotiable)
\[
P(\rho)=K\rho^5,\qquad c_s^2(\rho)=\frac{dP}{d\rho}=5K\rho^4,
\]
\[
U(\rho)=\frac{K}{4}\rho^5,\qquad h(\rho)=\frac{dU}{d\rho}=\frac{5K}{4}\rho^4.
\]

### 3.2 Matter PDE: gauged 4D GNLS
Define covariant derivatives (minimal coupling):
\[
D_t=\partial_t+\frac{i q}{\hbar}A_0,\qquad D_i=\partial_i-\frac{i q}{\hbar}A_i,\quad i\in\{x,y,z,w\}.
\]

**GNLS:**
\[
\boxed{\ i\hbar D_t\psi=\left[-\frac{\hbar^2}{2m}D_iD_i+V_{\rm conf}(\mathbf X;a,L)+h(|\psi|^2)\right]\psi\ }
\]
with \(h(|\psi|^2)=\tfrac{5K}{4}|\psi|^8\).

### 3.3 Current and continuity
**Gauge‑covariant mass/number current:**
\[
J^i_{\rm mass}=\frac{\hbar}{m}\,\Im\big(\psi^* D^i\psi\big),\qquad i\in\{x,y,z,w\}.
\]
**Charge current:** \(J^i_{\rm ch}=q\,J^i_{\rm mass}\).

**Continuity:**
\[
\boxed{\ \partial_t\rho+\partial_i(\rho v^i)=0\ }\qquad (\rho\equiv|\psi|^2).
\]

### 3.4 Madelung variables and Euler equation
Madelung split: \(\psi=\sqrt\rho\,e^{i\theta}\).

Gauge‑covariant velocity:
\[
\boxed{\ v_i=\frac{\hbar}{m}\left(\partial_i\theta-\frac{q}{\hbar}A_i\right)\ }
\]

Quantum pressure potential:
\[
Q(\rho)\equiv -\frac{\hbar^2}{2m}\frac{\nabla_4^2\sqrt\rho}{\sqrt\rho}.
\]

Euler‑like form (bulk 4D):
\[
\boxed{\ m(\partial_t+v_j\partial_j)v_i\ =\ q\,(E_i+v_j B_{ij})\ -\ \partial_i\big(V_{\rm conf}+h(\rho)+Q\big)\ }
\]
where \(B_{ij}=\partial_iA_j-\partial_jA_i\) and \(E_i= -\partial_t A_i-\partial_i A_0\) (sign conventions as usual).

**Automatic identity (vorticity ↔ gauge):**
\[
\Omega_{ij}\equiv\partial_i v_j-\partial_j v_i\ =\ -\frac{q}{m}B_{ij}.
\]
This is the cleanest statement of the “extra sector”: vorticity is not an arbitrary fluid degree; in the gauged theory it is tied to \(F_{ij}\).

### 3.5 Gauge PDE: 4+1D Maxwell with localization
Field strength: \(F_{MN}=\partial_MA_N-\partial_NA_M\).

Action density (plus gauge fixing):
\[
\mathcal L_{\rm EM}= -\frac{1}{4\mu_0}Z(w)F_{MN}F^{MN}\ -\ \frac{1}{2\xi\mu_0}(\partial_MA^M)^2.
\]

EOM:
\[
\boxed{\ \partial_M\big(Z(w)F^{MN}\big)+\frac{1}{\xi}\partial^N(\partial\cdot A)=\mu_0\,J^N_{\rm ch}+\mu_0\,J^N_{\rm ext}\ }
\]

### 3.6 Confinement potential (frozen functional form)
Use the harmonic far‑field trap and a throat gate/envelope. In the “stiff water” baseline:
\[
V_{\rm conf}(r,w;a,L)=V_w(w)\,[1-E(r;a)]\ +\ V_{\rm wall/cap}(r,w;a,L)\qquad (E=\exp(-r^4/a^4)).
\]
(Implementation detail: \(V_{\rm wall/cap}\) is a smooth barrier producing walls and endcaps; it must remain differentiable in \(a,L\) so \(\partial_{a,L}V_{\rm conf}\) exist.)

### 3.7 Geometry energy and wall law
Geometry energy (baseline):
\[
E_{\rm geom}(a,L)=P_{\rm vac}V(a,L)+\sigma A(a,L).
\]
With 4D “tube” measures:
\[
V(a,L)=\frac{4\pi}{3}a^3L,\qquad A(a,L)=4\pi a^2L+2\cdot\frac{4\pi}{3}a^3.
\]

Dynamic closure:
\[
M_a\ddot a+C_a\dot a=-\partial_a H_{\rm tot},\qquad
M_L\ddot L+C_L\dot L=-\partial_L H_{\rm tot}.
\]

Total energy (with EM):
\[
H_{\rm tot}=H_{\rm fluid}[\psi;a,L]+H_{\rm EM}[A;a,L]+E_{\rm geom}(a,L)+H_{\rm drive/diss}.
\]

Non‑negotiable force ledger identity:
\[
F_{a,L}(t)\equiv-\partial_{a,L}H_{\rm tot}=F^{(\psi)}_{a,L}+F^{(\rm EM)}_{a,L}+F^{(\rm geom)}_{a,L}+\cdots
\]
with
\[
F^{(\psi)}_{a,L}= -\int d^4X\,\rho(\mathbf X,t)\,\partial_{a,L}V_{\rm conf}(\mathbf X;a,L).
\]

### 3.8 Brane observables (projection + ports)
Projection (frozen): \(\rho_{\rm brane}=\mathcal P_W[\rho]\), \(\mathbf J_{\rm brane}=\mathcal P_W[\mathbf J]\).

Effort variable (frozen):
\[
 u\equiv\delta h(\rho_{\rm brane})\approx \left.\frac{dh}{d\rho}\right|_{\rho_{\rm brane,0}}\delta\rho_{\rm brane}=5K\rho_{\rm brane,0}^3\,\delta\rho_{\rm brane}.
\]

Measurement surface: \(\Gamma=S^2\) at \(R_3=r_{\rm port}\), \(w=0\).

Port basis: \(P_{\ell m}=Y_{\ell m}(\theta,\phi)\), at least \(\ell\le 2\).

Port amplitudes:
\[
 u_i(t)=\int_{\Gamma}\overline{P_i}\,u\,d\mu,\qquad
 j_i(t)=\int_{\Gamma}\overline{P_i}\,(\mathbf J_{\rm brane}\cdot\hat n)\,d\mu.
\]

**Multi‑flux outputs:**
- Mouth flux: \(j_{\rm mouth}(\mathbf s,t)=\mathbf J_{\rm brane}\cdot\hat n\) on \(\Gamma\).
- Leakage monitor (bulk exchange):
\[
J_w=\frac{\hbar}{2im}(\psi^*\partial_w\psi-\psi\,\partial_w\psi^*),\qquad
j_{\rm leak}(t)=\int_{\Omega_{xyz}}\big(J_w\big)\,d^3x\Big|_{w=\pm W_{\rm cut}}\ \text{(outward oriented)}.
\]

### 3.9 Open‑system closure + boundaries
Chemical‑potential reservoir coupling must be implemented explicitly (targeting \(\mu\) in a far‑field zone) *and* the domain must be non‑reflecting.

Matter absorbing layer (frozen):
\[
V\to V-i\Gamma_{\rm sponge}(\mathbf X).
\]
(For EM, use a compatible damping / sponge / PML‑style layer as needed.)

---

## 4) What we learned / what the math “naturally opened”

1. **4D GNLS → (continuity + Euler) is clean and matches the stiff EOS.** The symbolic decomposition yields the expected enthalpy term \(\nabla h\) and the expected quantum pressure term \(\nabla Q\).

2. **In the gauged model, vorticity becomes an EM‑controlled sector** via \(\Omega_{ij}=-(q/m)B_{ij}\). This is a natural, not forced, bridge to an EM interpretation.

3. **Poisson‑type structure is *available* but not automatic.** The longitudinal sector (after Helmholtz splitting + low‑frequency/static limit) yields Poisson candidates for a scalar potential; however, leakage terms, quantum pressure, and the transverse/gauge sector are genuine additional contributions.

4. **“Leakage” is a first‑class diagnostic, not a nuisance.** Because the brane is a projection of a 4D system, brane‑observed continuity generically contains boundary/projection terms sourced by \(J_w\). That is exactly the mathematical place where “refill faster than drain” ideas can live.

5. **Geometry selection is not guaranteed by \(P_{\rm vac}V+\sigma A\) alone.** Any preferred aspect ratio must emerge from field contributions (cavity stress, flow stress) or later shape‑field upgrades.

---

## 5) Remaining undecided / to freeze later (make explicit)

These do *not* block continuing, but must be fixed before a paper draft is final:

- Numerical values / nondimensionalization: \(K,m,\hbar,q,\mu,\omega_w,\sigma,P_{\rm vac},\lambda_{\rm conf}\).
- Gauge choice: \(\xi\) and practical gauge‑fixing strategy in numerics.
- Exact reservoir coupling form for “fixed \(\mu\)” (local \(\mu_{\rm loc}\) definition, coupling region, strength \(\gamma_\mu\)).
- Exact wave/vorticity injection operator (how \(\psi_{\rm drive}\) is implemented without violating constraints).
- EM boundary treatment (sponge/PML details) consistent with \(Z(w)\).
- When to upgrade from (a,L) DOFs to a shape field \(R(w,t)\) (trigger conditions).

---

## 6) Next steps (paper‑solidification plan)

### A) “Derivation‑complete” deliverables
1. **Write the unified action/Hamiltonian explicitly** (matter + EM + geom + localization + gauge‑fix + drive/diss) and derive all EOMs from it.
2. **Produce the reduced brane equations in a controlled limit** via a \(w\)-mode expansion (next section) and show exactly which terms survive:
   - longitudinal Poisson candidate,
   - transverse/vorticity (vector Poisson / magnetostatics),
   - leakage source terms.
3. **Energy/force ledger checks**: numerically verify \(-\partial_{a,L}H\) equals measured stresses (flow + EM).

### B) Minimal “nontrivial” \(w\)-ansatz (recommended)
Use a two‑mode HO truncation in \(w\):
\[
\psi(x,y,z,w,t)\approx \psi_0(x,y,z,t)\,\chi_0(w)\ +\ \epsilon\,\psi_1(x,y,z,t)\,\chi_1(w).
\]
This is the smallest ansatz that:
- respects the far‑field harmonic trap,
- allows genuine \(w\)-structure (mixing / leakage),
- and makes the projection \(W=|\chi_0|^2\) operationally meaningful.

### C) Minimal simulation program (to test “extra sectors”)
1. **No‑gauge baseline (fluid only):** confirm brane longitudinal sector reproduces the 1/r potential in the static limit *and* quantify leakage term size.
2. **Gauge on + localization:** check brane fields are 3D‑like (Coulomb scaling) and measure \(\Omega\leftrightarrow B\) numerically.
3. **Wave drive (spark):** inject rotational energy; test if a(t) stabilizes against \(\sigma\)-closure, and identify whether support is dominated by EM stress or flow stress.

### D) “Paper readiness” checklist
We are close to being able to write Paper 7, but the paper needs:
- One frozen model spec (this sheet) + notation.
- A short derivation section that is self‑contained (from the action).
- A controlled‑limit section (when Poisson emerges; when EM emerges; what extra terms appear).
- A minimal numerical demo set with convergence checks:
  - stability of a single throat,
  - brane field scaling,
  - leakage diagnostics,
  - (optionally) two‑throat superposition/interaction sanity check.

---

## 7) Repro commands / scripts (reference)

Symbolic derivation harness (SymPy):
- derive4d script: see `derive4d_brane_reduction_sympy.py`

Key runtime flags used this session included:
- `--report` (summary)
- `--no_gauge` (strip gauge sector)
- `--wrapup` (collect core identities)
- `--simplify` / `--full_simplify` (algebra cleanup)

(Exact CLI list should live at the top of the script under `--help`.)

