# Parent Throat Action — Actual Branch Export Checklist

## Purpose

This short note records the smallest parent-complete data packet that an actual moving-throat branch should export if we want to test it against the isotropic and weak-axisymmetric targets without post-target retuning.

In the local audit bundle this note is now also the frozen **V2-style reduced
branch packet** for the parent-background Galerkin demo:

- `branch_id = v2_local_parent_background_galerkin_demo`,
- `pre_target_freeze = true`,
- `target_blind = true`,
- `no_post_residual_refit = true`,
- `boundary_class = open_impedance_demo`,
- `branch_freeze_hash = fdf4aae0223c1faf`.

That metadata is intentionally narrower than the full PDE schema in the compact
program. It says this script is a reduced open-branch negative control, not yet
the completed moving-throat PDE branch.

---

## A. Isotropic branch packet

Export the following from one frozen isotropic grouped-`P2` branch.

### Geometry / wall

- stationary throat profile `R_0(w)`,
- grouped wall profile `\beta_2(w)`,
- parent wall coefficient functions `\mu_\eta,T_w,T_\Omega,K_\eta`,
- source normalization `\widehat m_0`.

### Support

- stable grouped-`P2` BdG couplings and pole data needed to build
  `B_0,B_2,B_4`.

### Conservative Maxwell/mixed

- port data needed to build
  `Z_0,Z_2,Z_4`.

### Outgoing transfer

- port data needed to build
  `N_0,N_2,N_4`.

### Derived wall integrals

- `M_\Sigma = \int \mu_\eta\beta_2^2`,
- `K_\Sigma = \int [T_w(\beta_2')^2 + (K_\eta+6T_\Omega)\beta_2^2]`.

### Algebraic tests

1. `D_0 = K_\Sigma-B_0-Z_0`, `D_2 = -(M_\Sigma+B_2+Z_2)`, `D_4 = -(B_4+Z_4)`.
2. One-pole surface:
   `D_0(B_4+Z_4) = 3(M_\Sigma+B_2+Z_2)^2`.
3. Static normalization:
   `\widehat m_0^2 N_0 / D_0 = 54Gc_s^5/(5a^5c^5)`.
4. Constant-prefactor branch if relevant:
   `P_2=0`, `P_4=0`.

---

## B. Weak-axisymmetric packet

Export the following from one frozen weak-axisymmetric branch around that isotropic branch.

### Wall slopes from parent action

- `\delta\mu_\eta, \delta T_w, \delta T_\Omega, \delta K_\eta`,
- derived slope integrals
  `\delta M_\Sigma`, `\delta K_\Sigma`.

### Support slopes

- `B_{01},B_{21},B_{41}`.

### Conservative mixed slopes

- `Z_{01},Z_{21},Z_{41}`.

### Outgoing slope

- `N_{01}`.

### Algebraic tests

1. `D_{01}=\delta K_\Sigma-B_{01}-Z_{01}`,
   `D_{21}=-(\delta M_\Sigma+B_{21}+Z_{21})`,
   `D_{41}=-(B_{41}+Z_{41})`.
2. Even gates:
   `K_1 = D_{21}+D_{01}/9`,
   `H_{\rm even}=D_{41}-(2/3)D_{21}-D_{01}/27`.
3. Prefactor slope:
   `\Xi_1 = N_{01}/N_0 - D_{01}/D_0`.
4. If one imposes the canonical even-compensated branch, verify the exact wall-slope solve:

   `\delta K_\Sigma = B_{01}+Z_{01}+27(B_{41}+Z_{41})`,

   `\delta M_\Sigma = -(B_{21}+Z_{21}) + 3(B_{41}+Z_{41})`.

---

## C. Why this packet is the right next target

This packet is small enough to be exported from an actual reduced branch and large enough to decide whether the parent-promoted wall block really helps.

If the packet fails, the failure will be localized:

- wall branch too stiff / too soft,
- support moments wrong,
- conservative mixed dressing wrong,
- outgoing slope wrong,
- or weak-axisymmetric compensation impossible on the realized branch.

That is exactly the kind of target-blind diagnostic the program now needs.

In the PDE-program vocabulary, the isotropic outputs are now reported directly
as the residual packet
\[
(R_{\rm pole},R_{\rm norm},R_{P2},R_{P4}),
\]
with
\[
R_{\rm pole}=D_0(B_4+Z_4)-3(M_\Sigma+B_2+Z_2)^2,
\]
\[
R_{\rm norm}=\widehat m_0^{\,2}P_0-\frac{54Gc_s^5}{5a^5c^5},
\]
and \(R_{P2}=P_2\), \(R_{P4}=P_4\).

---

## D. Concrete audit packet

The accompanying SymPy script no longer builds this section by choosing moments
from the target equations and then asserting those same targets back.
It now separates two logically different checks.

### D.1 Branch-derived wall block

Choose an actual isotropic wall branch

\[
R_0(w)=e^{-w^2/2},
\qquad
\beta_2(w)=e^{-w^2/2},
\]

with parent coefficients built from genuinely distinct \(w\)-profiles

\[
\mu_\Sigma(R,w)=1+\Bigl(1+\frac{w^2}{4}\Bigr)R,
\]

\[
T_{w,\Sigma}(R,w)=1+\Bigl(1+\frac{w^2}{6}\Bigr)R,
\qquad
T_{\Omega,\Sigma}(R,w)=\frac{1+\left(1+\frac{w^2}{8}\right)R}{6},
\]

and a quadratic potential
\[
U_\Sigma(R,w)=\frac12\,m(w)R^2
\]
whose coefficient \(m(w)\) is chosen so that the wall background equation is
satisfied by the displayed \(R_0(w)\). This gives the exported wall profiles

\[
\mu_\eta(w)=1+\Bigl(1+\frac{w^2}{4}\Bigr)e^{-w^2/2},
\]

\[
T_w(w)=1+\Bigl(1+\frac{w^2}{6}\Bigr)e^{-w^2/2},
\qquad
T_\Omega(w)=\frac{1+\left(1+\frac{w^2}{8}\right)e^{-w^2/2}}{6},
\]

and the closure-derived stiffness profile
\[
K_\eta(w)=w^2-1+e^{-w^2/2}\Bigl(\frac{w^2}{2}+\frac{w^4}{12}\Bigr).
\]

This produces nontrivial exported wall data

\[
M_\Sigma=\sqrt{\pi}\,\frac{36+13\sqrt6}{36},
\qquad
K_\Sigma=\sqrt{\pi}\,\frac{24+13\sqrt6}{24}.
\]

So the wall block is now an actual parent-derived profile export rather than a
constant packet, and the three directly chosen parent coefficients no longer
share the same \(1+R_0\) skeleton.

### D.2 Honest untuned support/mixed/transfer packet

Separately, the script now derives the support, conservative mixed, and
outgoing coefficients from an explicit **Galerkin resolvent of the actual
parent-background wall operator**, not from hand-written moment kernels and no
longer from the bare harmonic-oscillator inverse.

The exported profiles are

\[
\phi_B(w)=\beta_2(w),
\qquad
\phi_Z(w)=R_0(w)\beta_2(w),
\qquad
\phi_N(w)=\Bigl(1+\frac{w^2}{2}\Bigr)\beta_2(w),
\]

and they are tested against the actual grouped-`P2` wall operator on the chosen
branch
\[
\mathcal L_{\rm parent}[f]
=
-\partial_w\!\bigl(T_w\,\partial_w f\bigr)
+
\bigl(K_\eta+6T_\Omega\bigr)f.
\]

Let \(\psi_0,\psi_2,\psi_4\) be the normalized even Hermite functions used as
a Galerkin basis. The script checks their orthonormality directly and forms the
exact parent-background mass and stiffness matrices
\[
M^{\rm parent}_{mn}=\int dw\,\mu_\eta\,\psi_m\psi_n,
\qquad
K^{\rm parent}_{mn}
=
\int dw\,\Bigl[T_w\,\psi_m'\psi_n'+(K_\eta+6T_\Omega)\psi_m\psi_n\Bigr].
\]

On the live branch these matrices are symmetric and nondegenerate, with
\[
\det M_{\rm parent}=\frac{1475267\sqrt6+3638214}{1889568},
\qquad
\det K_{\rm parent}=\frac{3329998587+1455364309\sqrt6}{34012224}.
\]

Then:

- the support source has coefficients
  \[
  c_B=(\pi^{1/4},0,0),
  \]
- the mixed source has coefficients
  \[
  c_Z=\Bigl(\frac{\sqrt6\,\pi^{1/4}}{3},-\frac{\sqrt3\,\pi^{1/4}}{9},\frac{\pi^{1/4}}{18}\Bigr),
  \]
- the outgoing source has coefficients
  \[
  c_N=\Bigl(\frac{5\pi^{1/4}}{4},\frac{\sqrt2\,\pi^{1/4}}{4},0\Bigr).
  \]

The support and outgoing sources lie exactly in that three-mode basis, while
the mixed source \(\phi_Z=R_0\beta_2=e^{-w^2}\) is still exported through its
\((n=0,n=2,n=4)\) basis projection.

For the mixed profile, the script also records the explicit projection error
\[
\|\phi_Z-\Pi_{0,2,4}\phi_Z\|^2
=
\sqrt{\pi}\,\frac{162\sqrt2-229}{324},
\]
while the support and outgoing projection residuals vanish. So the truncation
is still measured rather than hidden.

The low-frequency response coefficients are now obtained by solving the
parent-background resolvent order by order:
\[
K_{\rm parent}a_0=c,\qquad
K_{\rm parent}a_2+M_{\rm parent}a_0=0,\qquad
K_{\rm parent}a_4+M_{\rm parent}a_2=0.
\]
The script verifies the corresponding Galerkin equations directly against
\(\mathcal L_{\rm parent}\) and then exports
\[
(B_0,B_2,B_4)\approx(0.7816221907762512,\,-0.6449280928722579,\,0.5321477935090640),
\]
\[
(Z_0,Z_2,Z_4)\approx(0.5083075562298627,\,-0.4143228212361644,\,0.3408547831911427),
\]
\[
(N_0,N_2,N_4)\approx(1.311583928278745,\,-1.062315116465984,\,0.8724932052126446).
\]

The script now also builds an **operator-adapted spectral basis** inside the
truncated Galerkin space by solving the generalized eigenproblem
\[
K_{\rm parent}v=\lambda\,M_{\rm parent}v.
\]
Numerically, this gives the parent-background eigenvalue ladders
\[
\lambda^{(3)}\approx(1.211928842993,\ 4.903138976811,\ 8.889028852429),
\]
\[
\lambda^{(4)}\approx(1.211926229000,\ 4.896292471678,\ 8.854823697872,\ 12.868952366471),
\]
\[
\lambda^{(5)}\approx(1.211919555609,\ 4.896014419797,\ 8.839287758009,\ 12.846837390175,\ 16.848273681346),
\]
with the script checking that the adapted modes are \(M_{\rm parent}\)-
orthonormal and diagonalize \(K_{\rm parent}\) to numerical tolerance.

The basis-growth audit is now therefore done in two layers:

1. exact projection-residual control in the Hermite trial basis,
2. operator-adapted coefficient export in the generalized eigenbasis.

At the projection level, extending the basis from \((0,2,4)\) to
\((0,2,4,6)\) drops the mixed-source residual exactly by
\[
\|\phi_Z-\Pi_{0,2,4,6}\phi_Z\|^2-\|\phi_Z-\Pi_{0,2,4}\phi_Z\|^2
=
-\frac{5\sqrt{\pi}}{17496},
\]
so the residual improves from
\[
\sqrt{\pi}\,\frac{162\sqrt2-229}{324}
\]
to
\[
\sqrt{\pi}\,\frac{8748\sqrt2-12371}{17496}.
\]

The corresponding low-frequency coefficients move only slightly:
\[
(B_0,B_2,B_4)_{4\text{-mode}}
\approx(0.7816273402373896,\,-0.6449331985854891,\,0.5321531025367578),
\]
\[
(Z_0,Z_2,Z_4)_{4\text{-mode}}
\approx(0.5083076399936368,\,-0.4143227047792866,\,0.3408544970748318),
\]
\[
(N_0,N_2,N_4)_{4\text{-mode}}
\approx(1.311690146078847,\,-1.062371646707426,\,0.8725297055800806),
\]
with every printed coefficient shifting by less than \(2\times 10^{-4}\) from
the 3-mode export.

The operator-adapted 5-mode audit then pushes one level further. The mixed
projection residual drops again, now by
\[
\|\phi_Z-\Pi_{0,2,4,6,8}\phi_Z\|^2-\|\phi_Z-\Pi_{0,2,4,6}\phi_Z\|^2
=
-\frac{35\sqrt{\pi}}{1259712},
\]
so the residual becomes
\[
\sqrt{\pi}\,\frac{629856\sqrt2-890747}{1259712}.
\]

The adapted 5-mode coefficients are
\[
(B_0,B_2,B_4)_{5\text{-mode}}
\approx(0.7816325278051803,\,-0.6449410093680708,\,0.5321624788442048),
\]
\[
(Z_0,Z_2,Z_4)_{5\text{-mode}}
\approx(0.5083084120357607,\,-0.4143247482933136,\,0.3408577693739790),
\]
\[
(N_0,N_2,N_4)_{5\text{-mode}}
\approx(1.3117016125078909,\,-1.0623879656363835,\,0.8725480934398114).
\]

The 4-mode to 5-mode drift stays below \(3\times 10^{-5}\) for every printed
coefficient, so the current parent-background packet is still a finite Galerkin
truncation, but it is no longer just a single untested basis choice and it is
no longer tied only to the raw Hermite trial basis.

More importantly, the untuned isotropic target residues themselves are stable
under this basis growth. On the 4-mode packet the script finds
\[
\text{one-pole}\approx -13.134593938872376,
\qquad
\text{static}\approx -10.33719584868593,
\]
\[
P_2\approx 0.37009844569768474,
\qquad
P_4\approx 0.8889149882257381,
\]
while on the adapted 5-mode packet it finds
\[
\text{one-pole}\approx -13.134428427110006,
\qquad
\text{static}\approx -10.337190829817821,
\]
\[
P_2\approx 0.37009832625683897,
\qquad
P_4\approx 0.8889265110475489.
\]
So every printed target residue drifts by less than \(2\times 10^{-4}\) across
the 3-mode \(\to\) 4-mode \(\to\) adapted 5-mode audit, while remaining far
from zero.

These quantities are still **not** tuned to the isotropic targets.

On that untuned packet, the script finds that the four isotropic target tests

\[
D_0(B_4+Z_4)-3(M_\Sigma+B_2+Z_2)^2,
\qquad
\widehat m_0^2P_0-\frac{54Gc_s^5}{5a^5c^5},
\qquad
P_2,
\qquad
P_4
\]

are all nonzero.

In the current live run the untuned parent-background packet gives
\[
\text{one-pole residue}\approx -13.13467188936131,
\]
\[
\text{static residue}\approx -10.33723418002391,
\]
\[
P_2\approx 0.3700576386561055,
\qquad
P_4\approx 0.8888349036886998.
\]

That is the important changed result of the strengthened audit:

> a concrete parent-derived wall block can be exported cleanly, but an honest
> untuned branch packet does **not** automatically satisfy the one-pole,
> normalization, or constant-prefactor targets.

So the current step-19 script is now evidence for three narrower claims:

1. the parent-throat wall block can be realized as real profile integrals,
2. a concrete parent-background Galerkin resolvent export can be built from
   explicit wall profiles and a measured mixed-source truncation,
3. the full target-satisfying branch export is still an open realization problem.
