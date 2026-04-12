# 5PN stages 207–208 — explicit overlap extractor and isotropic packet compiler

## What this turn accomplished

This turn pushed the Stage-205/206 normalized bridge one level deeper into the actual moving-throat overlap model.

The main new result is that the explicit finite-throat prototype already compresses to one **finite isotropic branch state**, and the two surviving endgame packets do **not** require the same pieces of that state.

That is the first honest extractor-level answer to the question “what data does the moving-throat branch really have to provide?”

---

## Stage 207 — explicit finite-throat overlap extractor

### Prototype used

I fixed the explicit finite-throat radial/axial prototype on
\[
s\in[0,L],
\]
with
\[
\chi_\eta(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{L},
\qquad
\phi_{\rm DN}(s)=\sqrt{\frac{2}{L}}\sin\!\frac{\pi s}{2L},
\]
and chose the simplest support/mixed profiles
\[
u(s)=\chi_\eta(s),
\qquad
w(s)=\phi_{\rm DN}(s).
\]

The exact overlaps are
\[
I_{\eta u}=1,
\qquad
I_{\eta\phi}=I_{\eta w}=I_{uw}=\frac{8}{3\pi}.
\]

### Overlap-renormalized couplings

If the raw microscopic inputs are
\[
(K,M,\lambda_B^{\rm raw},\varpi,c_{\eta U}^{\rm raw},\lambda_W^{\rm raw},\gamma^{\rm raw},K_U,\mu_U,K_W,\mu_W,T_U,\sigma),
\]
then the overlap extractor gives
\[
C=\frac{8}{3\pi}\lambda_B^{\rm raw},
\qquad
c_{\eta U}^{\rm eff}=c_{\eta U}^{\rm raw},
\qquad
\lambda_W^{\rm eff}=\frac{8}{3\pi}\lambda_W^{\rm raw},
\qquad
\gamma^{\rm eff}=\gamma^{\rm raw}.
\]

So in this explicit prototype the geometry renormalizes only the wall/support and wall/mixed amplitudes that actually feel the nontrivial half-wave overlap.

### Extracted isotropic branch state

The prototype then compresses to the finite isotropic state
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U),
\]
with
\[
\Omega_U=\sqrt{K_U/\mu_U},
\qquad
\Omega_W=\sqrt{K_W/\mu_W},
\]
\[
\chi_0=\frac{\gamma^{\rm raw} c_{\eta U}^{\rm raw}}{K_U},
\qquad
\epsilon_\eta=\frac{(c_{\eta U}^{\rm raw})^2}{K K_U},
\]
\[
Z_W=\frac{(\lambda_W^{\rm eff})^2}{K K_W},
\qquad
\epsilon_W=\frac{(\gamma^{\rm raw})^2(\lambda_W^{\rm eff})^2\sigma}{K_U K_W},
\qquad
\delta_U=\frac{\pi^2T_U}{L^2K_U}.
\]

This is exactly the coherent local-kernel variable set, now obtained from an explicit overlap model rather than inserted abstractly.

### Selected-branch transfer identity

Writing
\[
\epsilon
=
\epsilon_W\left(1-\frac{2}{11}\frac{\delta_U}{1+\delta_U}\right),
\]
the direct transfer shape is
\[
\mathcal T^2
=
\frac{Z_W(1+\chi_0)^2}{\Omega_W^2(1-\epsilon)^2}.
\]

The selected-branch demand ratio becomes
\[
R_{\rm target}
=
\frac{\Lambda(1-\epsilon_\eta)(1-\epsilon)^2}{Z_W(1+\chi_0)^2},
\qquad
\Lambda=
\frac{27\pi^2Gc_s^5K_W}{20a^5c^5\mu_W},
\]
so the exact product identity is
\[
R_{\rm target}\,\mathcal T^2
=
\Lambda_0(1-\epsilon_\eta),
\qquad
\Lambda_0=
\frac{27\pi^2Gc_s^5}{20a^5c^5}.
\]

That is the first explicit overlap-level version of the selected-branch identity.

---

## Stage 208 — isotropic Packet-A / Packet-B compiler from the extracted state

Stage 208 takes the Stage-207 state and compiles the endgame packets directly.

### Packet A on the isotropic branch

The conservative isotropic operator moments are
\[
D_0=K-B_0-Z_0,
\qquad
D_2=-(M+B_2+Z_2),
\qquad
D_4=-(B_4+Z_4),
\]
with the usual normalized branch defects
\[
\Delta_{\rm pole}
=
\bar u_4-4\bar u_2^2
=
-\frac{3D_2^2+D_0D_4}{D_0^2},
\]
\[
P_0=\frac{N_0}{D_0},
\qquad
\Delta_{\rm norm}
=
m_{\hat 0}^2P_0-rac{54Gc_s^5}{5a^5c^5}.
\]

Because the grouped lanes are already collapsed on the isotropic branch, the full Stage-199 branch packet reduces exactly to the scalar pair
\[
(\Delta_{\rm pole},\Delta_{\rm norm}).
\]

### Packet B from the same extracted state

Using the reference-branch exponents \((\chi_{0,*},\delta_{U,*},E_*,F_*)\), the exact direct monomials are
\[
\mathfrak C_{\rm tr}=\chi_0^{1+\delta_{U,*}}\,\delta_U^{1+\chi_{0,*}},
\]
\[
\mathfrak C_{\rm nt}=\frac{Z_W}{\Omega_W^2}\,\epsilon_W^{E_*}\,\delta_U^{-F_*},
\qquad
\mathfrak C_\eta=\epsilon_\eta.
\]

So the orbit packet is compiled by
\[
(\mathfrak C_{\rm tr},\mathfrak C_{\rm nt},\mathfrak C_\eta)
\to
(R_{\rm tr},R_{\rm nt},R_\eta)
\to
(q_{\rm tr},q_{\rm nt},q_\eta).
\]

### The main new structural theorem

The two packets separate cleanly.

**Packet A depends only on**
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W)
\quad\text{plus}\quad m_{\hat 0},
\]
and is blind to
\[
(\epsilon_W,\delta_U).
\]

**Packet B depends only on**
\[
(\chi_0,\delta_U,\epsilon_\eta,\epsilon_W,Z_W,\Omega_W),
\]
and is blind to the explicit support pair
\[
(C,\varpi),
\]
as well as the wall inertia pair
\[
(M,\Omega_U)
\]
once the extracted invariants are formed.

So the corrected combined isotropic endgame state is
\[
(K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U)
\]
plus the source factor
\[
m_{\hat 0}.
\]

That is the cleanest current statement of what the actual moving-throat branch still has to provide.

---

## What changed in the roadmap

Before this turn, the next step was loosely “extract isotropic branch data from the moving-throat overlap model.”

After this turn, the next theorem gate is much sharper:

1. extract the exact 11-scalar isotropic branch state
   \[
   (K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,\epsilon_W,Z_W,\delta_U),
   \]
2. extract or compute the source factor \(m_{\hat 0}\),
3. compile
   \[
   \Delta_{\rm pole},\qquad \Delta_{\rm norm},\qquad (q_{\rm tr},q_{\rm nt},q_\eta),
   \]
4. and test whether they all vanish on the actual moving-throat branch.

That is the direct continuation point.
