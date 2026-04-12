# 5PN continuation — normalized Packet-A / Packet-B bridge from the explicit overlap model

## What this session did

The previous session turned the Stage 199–201 endgame into working code for the **branch side** of the 5PN problem:

- Stage 202: exact Packet-A compiler,
- Stage 203: minimal isotropic single-mode overlap bridge,
- Stage 204: linear grouped outlet obstructions.

The main gap after that was structural rather than algebraic:

> Packet A and Packet B were still living in different microscopic languages.

Packet A was phrased in the explicit overlap variables
\((K,M,C,\varpi,G_U,G_W,R,\Omega_U,\Omega_W)\),
while the later weak-axisymmetric / similarity-orbit notes were phrased in the normalized coherent variables
\((\chi_0,\epsilon_\eta,Z_W,\delta_U,\dots)\).

This session closes that gap.

I built two new exact bridges:

1. a **normalized coherent isotropic Packet-A bridge**,
2. a **normalized coherent Packet-B / orbit bridge**.

So the explicit finite-throat overlap model and the later normalized monomial/orbit theorem now talk to each other directly.

---

## 1. Stage 205 — normalized coherent isotropic Packet-A bridge

### 1.1 Normalized coherent dictionary

Starting from the minimal isotropic overlap model, define the normalized coherent variables

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\qquad
Z_W = \frac{M G_W^2}{K\Omega_W^2}.
\]

Then the microscopic couplings are reconstructed exactly as

\[
G_U = \Omega_U\sqrt{\frac{\epsilon_\eta K}{M}},
\qquad
G_W = \Omega_W\sqrt{\frac{Z_W K}{M}},
\qquad
R = \chi_0\,\Omega_U\Omega_W\sqrt{\frac{Z_W}{\epsilon_\eta}}.
\]

So the isotropic overlap packet can now be written in the same normalized coherent variables used later in the 5PN notes.

---

### 1.2 The exact mixed-sector control parameter

The Stage-205 script shows that the conservative/outgoing mixed block enters through the single exact combination

\[
\epsilon_{\rm mix}
:=
\frac{\chi_0^2 Z_W}{\epsilon_\eta}.
\]

Then the one-port isotropic denominator becomes

\[
\Delta = \Omega_U^2\Omega_W^2\bigl(1-\epsilon_{\rm mix}\bigr).
\]

This is the cleanest exact reduction of the mixed block I have so far.

---

### 1.3 Exact isotropic conservative/outgoing coefficients in normalized form

The script gives the first exact normalized formulas for the isotropic Packet-A front end.

The BdG support moments remain

\[
B_0=\frac{C^2}{\varpi^2},
\qquad
B_2=\frac{C^2}{\varpi^4},
\qquad
B_4=\frac{C^2}{\varpi^6}.
\]

The mixed conservative static slot becomes

\[
Z_0
=
\frac{K}{M}
\frac{\epsilon_\eta + Z_W(1+2\chi_0)}{1-\epsilon_{\rm mix}}.
\]

The outgoing-transfer static slot becomes

\[
N_0
=
\frac{K}{M\Omega_W^2}
\frac{Z_W(1+\chi_0)^2}{(1-\epsilon_{\rm mix})^2}.
\]

The higher conservative mixed slots are packaged as

\[
Z_2 = \frac{K}{M}\,\Sigma_2,
\qquad
Z_4 = \frac{K}{M}\,\Sigma_4,
\]

with exact closed expressions \(\Sigma_2,\Sigma_4\) produced in the script.

So the isotropic Packet-A bundle is now explicitly controlled by

- the **support pair** \((C,\varpi)\),
- the **wall pair** \((K,M)\),
- the **transport scales** \((\Omega_U,\Omega_W)\),
- the **normalized mixed ratios** \((\chi_0,\epsilon_\eta,Z_W)\).

---

### 1.4 Exact support-blind theorem for the outgoing-transfer bundle

One of the strongest exact results of this session is:

\[
\frac{\partial N_0}{\partial C}
=
\frac{\partial N_0}{\partial \varpi}
=
\frac{\partial N_2}{\partial C}
=
\frac{\partial N_2}{\partial \varpi}
=
\frac{\partial N_4}{\partial C}
=
\frac{\partial N_4}{\partial \varpi}
=0.
\]

So on the isotropic branch, the entire outgoing-transfer packet
\((N_0,N_2,N_4)\)
is **exactly support-blind** in the explicit BdG pair \((C,\varpi)\).

That does **not** mean the full Packet-A verdict is support-blind, because

\[
D_0 = K - B_0 - Z_0
\]

still depends on \(B_0=C^2/\varpi^2\), and similarly for \(D_2,D_4\).

So the correct separation is:

- the **outgoing transfer side** is support-blind,
- the **conservative wall operator** is not.

This is precisely the kind of structural split we needed to make the roadmap honest.

---

### 1.5 Exact normalized compatibility surface

Combining the 5PN one-pole gate with the 2.5PN/4PN normalization target still gives the same exact scalar compatibility law,

\[
\frac{N_0}{P_{0,\rm target}}
=
\frac{3(M+B_2+Z_2)^2}{B_4+Z_4},
\qquad
P_{0,\rm target}=
\frac{54Gc_s^5}{5a^5c^5m_{\hat 0}^{\,2}}.
\]

But now it is written in normalized coherent language:

\[
\frac{K}{M\Omega_W^2}
\frac{Z_W(1+\chi_0)^2}{(1-\epsilon_{\rm mix})^2}
\frac{1}{P_{0,\rm target}}
=
\frac{3\left(M+\frac{C^2}{\varpi^4}+\frac{K}{M}\Sigma_2\right)^2}
{\frac{C^2}{\varpi^6}+\frac{K}{M}\Sigma_4}.
\]

That is the first exact normalized form of the isotropic Packet-A theorem gate.

So the remaining isotropic Packet-A problem is no longer “solve all overlap data somehow.”
It is:

> compute the scalar set
> \((K,M,C,\varpi,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W,m_{\hat0})\)
> on the physical branch.

Then the entire isotropic Packet-A verdict is immediate.

---

## 2. Stage 206 — normalized coherent Packet-B / orbit bridge

### 2.1 Exact normalized monomials

The later weak-axisymmetric/similarity-orbit notes were encoded in the direct monomials

\[
\mathfrak C_{{\rm tr},*},
\qquad
\mathfrak C_{{\rm nt},*},
\qquad
\epsilon_\eta.
\]

This session rewrites them directly in the normalized coherent variables.

With

\[
\chi_0 = \frac{R G_U}{\Omega_U^2 G_W},
\qquad
\epsilon_W = \frac{R^2\sigma}{\Omega_U^2\Omega_W^2},
\qquad
\epsilon_\eta = \frac{M G_U^2}{K\Omega_U^2},
\]

the exact monomials are

\[
\mathfrak C_{\rm tr}
=
\chi_0^{1+\delta_{U,*}}\,
\delta_U^{1+\chi_{0,*}},
\]

\[
\mathfrak C_{\rm nt}
=
\frac{M G_W^2}{K\Omega_W^4}
\epsilon_W^{E_*}
\delta_U^{-F_*},
\]

\[
\mathfrak C_\eta = \epsilon_\eta.
\]

So Packet B is now coded directly in the same microscopic language as Packet A.

---

### 2.2 Exact invariant packet and quotient coordinates

Relative to a reference branch,

\[
R_{\rm tr} = \frac{\mathfrak C_{\rm tr}}{\mathfrak C_{{\rm tr},\rm ref}},
\qquad
R_{\rm nt} = \frac{\mathfrak C_{\rm nt}}{\mathfrak C_{{\rm nt},\rm ref}},
\qquad
R_\eta = \frac{\epsilon_\eta}{\epsilon_{\eta,\rm ref}}.
\]

Then the exact quotient-coordinate packet is just

\[
q_{\rm tr} = \ln R_{\rm tr},
\qquad
q_{\rm nt} = \ln R_{\rm nt},
\qquad
q_\eta = \ln R_\eta.
\]

The new script verifies that this packet roundtrips exactly through the existing Stage-200/201 common compiler.

So Packet B is now operational in normalized coherent variables rather than only in abstract quotient notation.

---

### 2.3 Exact normalized monomial-drift matrix

The script then re-derives the direct monomial-drift matrix from first principles in the normalized variable set

\[
(d\ln G_W,
 d\ln G_U,
 d\ln R,
 d\ln K,
 d\ln M,
 d\ln\Omega_U,
 d\ln\Omega_W,
 d\ln\delta_U).
\]

The exact result is

\[
M_{\rm norm}=
\begin{pmatrix}
-(1+\delta_{U,*}) & 1+\delta_{U,*} & 1+\delta_{U,*} & 0 & 0 & -2(1+\delta_{U,*}) & 0 & 1+\chi_{0,*} \\
2 & 0 & 2E_* & -1 & 1 & -2E_* & -(4+2E_*) & -F_* \\
0 & 2 & 0 & -1 & 1 & -2 & 0 & 0
\end{pmatrix}.
\]

The script verifies this exactly against the Stage-13 matrix.

So the Stage-13 normalized similarity-orbit theorem is now implemented as code rather than only as notes.

---

### 2.4 Exact triangular zero-defect solve

Setting the three monomial drifts to zero and solving the system gives the same triangular co-drift laws as the notes.

The script recovers exactly:

\[
d\ln\delta_U
=
-\frac{1+\delta_{U,*}}{1+\chi_{0,*}}
\bigl(d\ln R + d\ln G_U - d\ln G_W - 2d\ln\Omega_U\bigr),
\]

\[
d\ln M
=
d\ln K - 2d\ln G_U + 2d\ln\Omega_U,
\]

\[
d\ln\Omega_W
=
\frac{d\ln G_W - d\ln G_U + (1-E_*)d\ln\Omega_U + E_* d\ln R - \tfrac{F_*}{2}d\ln\delta_U}{E_*+2}.
\]

So the zero-defect tangency problem is now executable:

> once the actual branch gives drifts of
> \((G_W,G_U,R,K,\Omega_U)\),
> the required co-drifts of
> \((\delta_U,M,\Omega_W)\)
> are fixed exactly.

That is the cleanest current reduced theorem gate on the Packet-B side.

---

## 3. Combined reading of Stages 205 and 206

The main gain of this session is that Packet A and Packet B now share one microscopic vocabulary.

### Shared microscopic sector

Both packets depend on the same normalized coherent core:

\[
(K,M,G_U,G_W,R,\Omega_U,\Omega_W).
\]

Equivalently, after normalization,

\[
(K,M,\Omega_U,\Omega_W,\chi_0,\epsilon_\eta,Z_W).
\]

### Packet-A-only data

Packet A still needs the explicit support pair

\[
(C,\varpi),
\]

because the conservative wall operator depends on

\[
B_0,B_2,B_4.
\]

### Packet-B-only data

Packet B needs the tracking variable

\[
\delta_U,
\]

plus the reference-branch exponents and calibration constants

\[
(\chi_{0,*},\delta_{U,*},E_*,F_*).
\]

So the full reduced 5PN data request has now sharpened further.

The completed moving-throat branch does **not** need to hand us unrelated large symbolic objects.
It needs to hand us:

1. the shared normalized coherent core,
2. the support pair \((C,\varpi)\),
3. the tracking variable \(\delta_U\),
4. and the source factor \(m_{\hat0}\).

From that, the present code can now produce:

- the isotropic Packet-A residual,
- the orbit Packet-B residual,
- and the exact zero-defect tangency conditions.

---

## 4. Best continuation point after this session

The next clean move is now very explicit.

### Step 1 — extract isotropic branch scalars
From the moving-throat branch, compute

\[
K,
\ M,
\ C,
\ \varpi,
\ \Omega_U,
\ \Omega_W,
\ \chi_0,
\ \epsilon_\eta,
\ Z_W,
\ \delta_U,
\ m_{\hat0}.
\]

### Step 2 — run Packet A
Feed those into Stage 205 to get

- \(B_n,Z_n,N_n,D_n\),
- \(u_2,u_4,P_0\),
- \(\Delta_{\rm pole},\Delta_{\rm norm}\).

### Step 3 — run Packet B
Feed the same normalized core plus \(\delta_U\) into Stage 206 to get

- \((R_{\rm tr},R_{\rm nt},R_\eta)\),
- \((q_{\rm tr},q_{\rm nt},q_\eta)\),
- and the exact zero-defect tangency test.

### Step 4 — only then reopen anisotropy
Once the isotropic packet and orbit packet are both computed, the weak-axisymmetric Stage-204 obstruction test becomes the right next refinement.

So after this session, the continuation point is no longer vague.
It is:

> compute the actual isotropic branch scalars in the shared normalized coherent language.

That is now the shortest honest route back into the full 5PN endgame.
