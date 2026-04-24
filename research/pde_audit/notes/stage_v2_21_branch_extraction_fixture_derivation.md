# Stage V2-21 — Branch-Extraction Fixture and Observable Packet

## 0. Purpose

This stage turns the Volume 2 weak-form scaffold into an executable branch-extraction fixture.

The fixture does **not** solve the nonlinear moving-throat PDE. It does the next necessary thing: it accepts a frozen branch packet, extracts the grouped real `P2` observables, enforces open-throat and stability gates, and evaluates the target residuals without changing branch data after the residuals are known.

The implementation file is:

```text
stage_v2_21_branch_extraction_fixture.py
```

The fixture writes a machine-readable extraction packet:

```text
stage_v2_21_branch_extraction_packet.json
```

and includes a sample manifest:

```text
stage_v2_21_sample_branch_manifest.json
```

---

## 1. Input manifest

The branch manifest has the form

```json
{
  "schema": "stage_v2_21_branch_extraction_fixture/v1",
  "branches": [
    {
      "name": "branch name",
      "geometry": {
        "R_exit": 0.35,
        "boundary_class": "open_impedance",
        "Y_L_limit": 0.0
      },
      "target": {
        "constants": {
          "G": 1.0,
          "c_s": 1.0,
          "c": 1.0,
          "a": 1.0,
          "mhat0": 1.0,
          "S_port": 1.0,
          "theta_tail": 1.0
        }
      },
      "lanes": {
        "20": { "...": "lane data" },
        "21": { "...": "lane data" },
        "22": { "...": "lane data" }
      }
    }
  ]
}
```

Each grouped lane may be supplied in either of two ways.

### 1.1 Primitive reduced-mode input

```json
{
  "K": 1.6,
  "M": 0.9,
  "bdg_modes": [
    {"coupling": 0.42, "varpi": 2.7}
  ],
  "mixed_ports": [
    {"Omega_U": 3.0, "Omega_W": 4.0, "R": 0.65, "g_U": 0.28, "g_W": 0.52}
  ]
}
```

### 1.2 Direct coefficient input

```json
{
  "direct_coefficients": {
    "K": 1.0,
    "M": 1.0,
    "B0": 0.0,
    "B2": 0.0,
    "B4": 3.0,
    "Z0": 0.0,
    "Z2": 0.0,
    "Z4": 0.0,
    "N0": 10.8,
    "N2": -21.6,
    "N4": -54.0
  }
}
```

The direct-coefficient path is included only as an algebraic self-test and for importing future pre-extracted numerical coefficients. The primitive path is the intended route for real branch data.

---

## 2. Open-throat gate

The fixture enforces the open-organ-pipe patch by requiring

\[
R_{\rm exit}>0,
\]

and

\[
\texttt{boundary\_class}=\texttt{open\_impedance}.
\]

It rejects or fails hard-cap geometry:

\[
R(L)=0
\]

as a branch-realization condition for the support problem. In this fixture, the support-mode Neumann-like condition is represented by the open impedance limit

\[
T_w q_w(L,\omega)+Y_L(\omega)q(L,\omega)=0,
\qquad
Y_L\to0
\]

for AC support modes, while DC leakage is allowed through the finite exit.

---

## 3. Lane-level extraction formulas

For each grouped lane

\[
A\in\{20,21,22\},
\]

the fixture computes BdG moments, Maxwell/mixed moments, outgoing-transfer moments, conservative operator moments, response moments, and prefactor moments.

### 3.1 BdG support moments

For stable positive-energy support coordinates with couplings \(c_\alpha\) and frequencies \(\varpi_\alpha>0\),

\[
B_0=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^2},
\]

\[
B_2=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^4},
\]

\[
B_4=\sum_\alpha \frac{c_\alpha^2}{\varpi_\alpha^6}.
\]

These are the Schur-complement support moments from the stable BdG closure.

### 3.2 Conservative Maxwell/mixed moments

For a port with

\[
\Omega_U,
\qquad
\Omega_W,
\qquad
R,
\qquad
g_U,
\qquad
g_W,
\]

define

\[
\Delta=\Omega_U^2\Omega_W^2-R^2,
\]

\[
S=\Omega_U^2+\Omega_W^2,
\]

\[
Q=g_U^2\Omega_W^2+2g_Ug_WR+g_W^2\Omega_U^2,
\]

\[
H=g_U^2+g_W^2.
\]

Then

\[
Z_0=\frac{Q}{\Delta},
\]

\[
Z_2=\frac{QS-H\Delta}{\Delta^2},
\]

\[
Z_4=\frac{Q(S^2-\Delta)-SH\Delta}{\Delta^3}.
\]

For multiple mixed ports, the fixture sums these contributions over ports.

### 3.3 Outgoing-transfer moments

Define the outgoing-port numerator

\[
P=\Omega_U^2g_W+Rg_U.
\]

Then

\[
N_0=\frac{P^2}{\Delta^2},
\]

\[
N_2=\frac{2P(PS-\Delta g_W)}{\Delta^3},
\]

\[
N_4=
\frac{
\Delta^2 g_W^2
-2\Delta P^2
-4\Delta PSg_W
+3P^2S^2
}{\Delta^4}.
\]

Again, for multiple ports, the fixture sums over ports.

### 3.4 Conservative wall operator

The lane operator moments are

\[
D_0=K-B_0-Z_0,
\]

\[
D_2=-(M+B_2+Z_2),
\]

\[
D_4=-(B_4+Z_4).
\]

### 3.5 Normalized response

The normalized conservative response is

\[
Y(\omega)=\frac{D_0}{D_0+D_2\omega^2+D_4\omega^4+O(\omega^6)}.
\]

The fixture verifies symbolically that

\[
u_2=-\frac{D_2}{D_0},
\]

\[
u_4=\frac{D_2^2-D_0D_4}{D_0^2}.
\]

### 3.6 Outgoing prefactor

The outgoing prefactor is

\[
\mathrm{Pref}(\omega)
=
\frac{D_0\left(N_0+N_2\omega^2+N_4\omega^4\right)}{
\left(D_0+D_2\omega^2+D_4\omega^4\right)^2
}.
\]

The fixture verifies symbolically that

\[
P_0=\frac{N_0}{D_0},
\]

\[
P_2=\frac{D_0N_2-2D_2N_0}{D_0^2},
\]

\[
P_4=
\frac{
D_0^2N_4
-2D_0(D_2N_2+D_4N_0)
+3D_2^2N_0
}{D_0^3}.
\]

---

## 4. Grouped real `P2` decomposition

For any grouped triple

\[
x=(x_{20},x_{21},x_{22}),
\]

the fixture computes

\[
\bar x=\frac{x_{20}+2x_{21}+2x_{22}}{5},
\]

\[
a_x=\frac{2x_{20}-x_{21}-x_{22}}{10},
\]

\[
b_x=\frac{x_{21}-x_{22}}{2}.
\]

The inverse map is verified symbolically:

\[
x_{20}=\bar x+4a_x,
\]

\[
x_{21}=\bar x-a_x+b_x,
\]

\[
x_{22}=\bar x-a_x-b_x.
\]

The anisotropy norm is

\[
A_x^2=4a_x^2+\frac45b_x^2.
\]

The axisymmetric weak-splitting diagnostic is also reported:

\[
b_x-3a_x.
\]

---

## 5. Target residual packet

The fixture evaluates the isotropic trace residuals

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
\]

\[
R_{P4}=P_4,
\]

\[
R_{\rm tail}=\Theta_{\rm tail}\left(\frac{c}{c_s}\right)^3-1.
\]

It also reports the equivalent quadrupole coefficient residual

\[
\gamma_{\rm eff}-\gamma_{\rm GR}
=
\widehat m_0^{\,2}\mathcal S_{\rm port}P_0\frac{a^5}{27c_s^5}
-
\frac{2G}{5c^5}.
\]

---

## 6. Stability gates

The fixture currently enforces the reduced stability gates:

\[
D_0>0,
\]

\[
B_4+Z_4>0,
\]

\[
N_0\ne0,
\]

\[
M>0,
\]

and for every primitive mixed port,

\[
\Delta=\Omega_U^2\Omega_W^2-R^2>0.
\]

These are reduced checks. A future PDE run should also supply the full positivity certificate for the BdG signature and full wall/mixed Hamiltonian.

---

## 7. Built-in fixture branches

The script includes two built-in branches.

### 7.1 Calibrated direct-coefficient branch

This is an algebraic self-test, not a physical primitive realization.

It uses normalized constants

\[
G=c=c_s=a=\widehat m_0=\mathcal S_{\rm port}=\Theta_{\rm tail}=1,
\]

so

\[
P_0^{\rm target}=\frac{54}{5}=10.8.
\]

The direct coefficients are chosen so that

\[
D_0=1,
\qquad
M+B_2+Z_2=1,
\qquad
B_4+Z_4=3,
\]

which gives

\[
D_0(B_4+Z_4)=3(M+B_2+Z_2)^2.
\]

Then

\[
N_0=10.8,
\]

and the constant-prefactor branch fixes

\[
N_2=\frac{2D_2N_0}{D_0}=-21.6,
\]

\[
N_4=\frac{N_0(D_2^2+2D_0D_4)}{D_0^2}=-54.
\]

The output shows this branch passing all gates up to floating-point roundoff.

### 7.2 Primitive open-throat demo branch

This branch uses actual reduced primitive data: two BdG support modes and one Maxwell/mixed port. It is stable and open-throat, but it is not tuned to the target surface.

The output is intentionally a falsifying packet:

- open gate passes,
- stability gate passes,
- grouped anisotropy vanishes because the lanes are isotropic,
- but the one-pole, normalization, and constant-prefactor residuals do not vanish.

This demonstrates the role of the fixture: it does not rescue a branch. It reports the residuals.

---

## 8. Output from this run

The script reported:

```text
symbolic_checks: 13/13 passed
```

For the calibrated branch:

```text
target_packet_pass: True
D0_bar: 1
N0_bar: 10.8
P0_bar: 10.8
P0_target: 10.8
R_pole: 0
R_norm: 0
R_P2: 0
R_P4: 1.4210854715202e-14
R_tail: 0
```

For the primitive open-throat demo branch:

```text
target_packet_pass: False
D0_bar: 1.546870256628791
N0_bar: 0.001146719334105218
P0_bar: 0.0007413157821033734
P0_target: 10.8
R_pole: -2.459877868885326
R_norm: -10.7992586842179
R_P2: 0.0009676829082644174
R_P4: 0.0008900366908790024
R_tail: 0
```

---

## 9. Interpretation

V2-21 changes the workflow from handwritten coefficient checks to a reproducible branch-extraction protocol.

A future PDE or numerical branch run should now do this:

1. freeze the branch manifest before target evaluation,
2. extract primitive wall, BdG, Maxwell/mixed, and outgoing-port data,
3. write them into the manifest schema,
4. run the fixture,
5. read the residual packet.

The important discipline is that the fixture has no rescue logic. It does not change support cardinality, boundary class, gauge convention, or port normalization after the residuals are known.

---

## 10. Carry-forward status

The fixture is ready for mock data and reduced numerical branch data.

The next technical stage should be either:

1. **V2-22A:** build a profile-to-coefficient adapter, where actual sampled profiles \(\chi_\eta(s),\phi_\alpha(s),u_r(s),w_r(s)\) are integrated numerically to produce \(B_n,Z_n,N_n\); or
2. **V2-22B:** build a weak-axisymmetric manifest adapter that adds first-order lane slopes and automatically reports \(\Xi_1=P_1/P_0\), \(K_1\), and the hidden-even residual.

The first option is the more direct continuation toward an actual PDE branch. The second is the more direct continuation toward the same-charge / weak-axisymmetric prefactor packet.
