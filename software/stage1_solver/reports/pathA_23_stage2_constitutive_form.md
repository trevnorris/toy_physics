# pathA_23 Stage 2 constitutive form

## VERDICT

`FAIL_UNSPECIFIED_SUBSTRUCTURE` - with the full first-gradient in-plane strain content present, the gate can detect either
zero or nonzero symmetric shear, but the independently motivated substructure record does not determine which limit the medium
realizes. The symmetric shear modulus
\[
\mu_{\rm br}
\]
is therefore not derived as zero and not derived as nonzero.

Conditional-status flag: `NOT_CONDITIONAL_UNDERDETERMINED_FAILURE`. This is not a conditional success from a postulated shear
law. It is a Stage-2 failure to derive the needed substructure well enough to choose \(\mu_{\rm br}=0\) or
\(\mu_{\rm br}>0\).

Stage 2 stops here. No Stage-3 no-leak closure, Stage-4 spectrum pass, Stage-5 Maxwell/gauge result, charge result, or cone
result is claimed.

## Substructure Characterization

The tested question is not "does the GP/NLS mean-field have shear?" It does not: a scalar fluid has
\(\mu_{\rm shear}=0\). The Track-B question is whether the cohesive substructure beneath that mean-field carries a finite
in-plane symmetric shear modulus.

The redo uses this independent criterion:

| physical limit | independent criterion | resulting \(\mu_{\rm br}\) |
| --- | --- | --- |
| simple fluid / soap-film limit | surface tension or cohesion exists, but no persistent in-plane neighbor memory; tangential shear relaxes as flow | \(0\) |
| elastic membrane / coherent-network limit | finite correlation length plus persistent neighbor memory over the probe window, with a specified affine/network free energy | \(>0\) |
| current record | cohesion, healing/correlation length, tautness, and high-frequency elasticity are motivated, but persistent in-plane neighbor locking, a network free energy, and the relevant frequency regime are not independently specified | undetermined |

The earlier failure was tautological because it selected \(W=W(J)\) before the test. This redo does not select \(J\)-only.
The algebra keeps the deviatoric invariant and then asks the substructure classifier to decide \(\mu_{\rm br}\).

Outcome: the record does not justify the simple-fluid limit either. It also does not justify the coherent-network limit as a
derivation. Treating "elastic membrane" as already proving persistent neighbor memory would be a postulate, not a derivation.

The gyrostat / Cosserat / MacCullagh route remains refused as an independent derivation. Spinning gyrostats, directors, or a
micro-rotation field would have to be supplied by the medium's substructure before asking for the MacCullagh answer; adding
them here would be reverse-engineering. Symmetric Cauchy shear is a legitimate separate invariant, but its coefficient is still
not fixed by the present substructure data.

## Energy And Able-To-Fail Gate

Use
\[
\theta=D_a u^a,\qquad
e_{ab}=\frac12(D_a u_b+D_b u_a),\qquad
e_{\langle ab\rangle}=e_{ab}-\frac13\theta\delta_{ab}.
\]

The full first-gradient symmetric-strain energy tested is
\[
\mathcal U_\parallel
=\frac12K_{\rm br}\theta^2
+\mu_{\rm br}e_{\langle ab\rangle}e^{\langle ab\rangle}.
\]

The scripts verify that:

| branch/control | classifier input | algebraic result |
| --- | --- | --- |
| actual record | cohesion/healing motivated; neighbor memory and network free energy unspecified | `FAIL_UNSPECIFIED_SUBSTRUCTURE` |
| able-to-fail no-shear control | simple fluid / soap-film facts: no persistent neighbor memory | \(\mu_{\rm br}=0\), transverse block zero |
| able-to-fail shear control | coherent elastic network facts: persistent neighbor memory plus specified network energy | \(\mu_{\rm br}>0\), transverse stiffness detected |

For a Fourier mode,
\[
K_{ab}(k)=\mu_{\rm br}k^2\delta_{ab}
+\left(K_{\rm br}+\frac13\mu_{\rm br}\right)k_ak_b,
\]
with characteristic polynomial
\[
(\lambda-\mu_{\rm br}k^2)^2
\left(\lambda-\left(K_{\rm br}+\frac43\mu_{\rm br}\right)k^2\right).
\]

So the gate is genuinely able to fail both directions:

- If \(\mu_{\rm br}=0\), the law earns `FAIL_NO_TRANSVERSE_STIFFNESS` with spectrum
  \(0,0,K_{\rm br}k^2\).
- If \(\mu_{\rm br}>0\), the law earns a shear-carrying Cauchy candidate with two transverse eigenvalues
  \(\mu_{\rm br}k^2\), but also the longitudinal signature
  \((K_{\rm br}+4\mu_{\rm br}/3)k^2\), to be flagged as `FAIL_CAUCHY_STRAY_LONGITUDINAL` in Stage 4.
- The current record lands neither branch because \(\mu_{\rm br}\) is not derived.

## Stress, Couple Stress, Boundary Work

For the symmetric Cauchy branch, the force stress is
\[
\sigma_{ab}
=\frac{\partial\mathcal U_\parallel}{\partial(D_a u_b)}
=K_{\rm br}\theta\delta_{ab}+2\mu_{\rm br}e_{\langle ab\rangle}.
\]

This stress is symmetric. The first-gradient symmetric-strain branch has no independent spin or curvature dependence, so
\[
M_{cab}=0.
\]

For a brane region \(M\) with outward in-brane unit normal \(\nu_a\),
\[
\delta U_\parallel
=-\int_M d^3\xi\,D_a\sigma_{ab}\,\delta u_b
+\int_{\partial M}d^2S\,\nu_a\sigma_{ab}\,\delta u_b.
\]
The boundary traction object is therefore
\[
t_b^{\rm brane}=\nu_a\sigma_{ab},
\]
and the finite-thickness force density is
\[
f_b^{(4)}=B_\ell(w)D_a\sigma_{ab}.
\]

Units are restored:

| quantity | dimension |
| --- | --- |
| \(K_{\rm br},\mu_{\rm br},\sigma_{ab},t_b^{\rm brane}\) | \(E/L^3=M L^{-1}T^{-2}\) |
| \(D_a\sigma_{ab}\) | \(E/L^4=M L^{-2}T^{-2}\) |
| \(B_\ell D_a\sigma_{ab}\) | \(E/L^5=M L^{-3}T^{-2}\) |
| \(\int dt\,d^2S\,t_b\delta u_b\) | action |

These are hand-off objects only if a later user gate postulates or independently derives a Cauchy/network branch. They are not
selected by this Stage-2 verdict.

## Angular Momentum

The local angular-momentum closure condition is
\[
\sigma_{ab}-\sigma_{ba}+D_cM_{cab}+\dot s_{ab}=0.
\]

For the symmetric-strain branch,
\[
\sigma_{ab}=\sigma_{ba},\qquad M_{cab}=0,
\]
and no intrinsic spin \(s_{ab}\) is needed. Thus a Cauchy symmetric-shear branch closes angular momentum without a Cosserat
sector.

The MacCullagh contrast remains different:
\[
\mathcal U_R=\frac12\mu_R r_ar_a,\qquad
r_a=\epsilon_{abc}D_bu_c,
\]
\[
\sigma^R_{ab}=\mu_R(D_a u_b-D_b u_a).
\]
This first-gradient stress is antisymmetric. Without an independently derived spin/couple sector it triggers
`FAIL_ANTISYMMETRIC_STRESS_NO_SPIN_CLOSURE`. The curl-only potential also keeps the known C5 obstruction: the potential is
invariant under \(u\to u+\nabla\chi\), but \( \frac12\rho_\parallel \dot u^2 \) is not invariant for time-dependent \(\chi\)
unless a scalar-potential or constraint sector is supplied.

## Hamiltonian And Modes

For the symmetric-strain branch,
\[
H_k=\frac12\rho_\parallel |\dot u_k|^2
+\frac12u^*_{a,k}K_{ab}(k)u_{b,k}.
\]

The stiffness eigenvalues are
\[
\mu_{\rm br}k^2,\qquad
\mu_{\rm br}k^2,\qquad
\left(K_{\rm br}+\frac43\mu_{\rm br}\right)k^2.
\]
Therefore \(H_k\) is nonnegative for all Fourier modes when
\[
\rho_\parallel>0,\qquad \mu_{\rm br}\ge0,\qquad K_{\rm br}+\frac43\mu_{\rm br}\ge0.
\]
A simple sufficient condition is \(K_{\rm br}\ge0,\mu_{\rm br}\ge0\). The scripts include negative controls showing that
\(\mu_{\rm br}<0\) gives a transverse ghost and \(K_{\rm br}+4\mu_{\rm br}/3<0\) gives a longitudinal ghost.

Cosserat/micropolar bookkeeping branch, not selected:
\[
\mathcal U_{\rm Cosserat}
=\frac12K_c\theta^2+\mu_c e_{\langle ab\rangle}e^{\langle ab\rangle}
+\frac12\kappa_c(2\varpi_a-r_a)^2
+\frac12A_cD_a\varpi_bD^a\varpi^b+\cdots .
\]
For a transverse pair at \(k=k\hat z\), the checked dynamic block is
\[
\det\begin{pmatrix}
(\mu_c+\kappa_c)k^2-\rho_\parallel\omega^2 & -2\kappa_ck\\
-2\kappa_ck & 4\kappa_c+A_ck^2-I_c\omega^2
\end{pmatrix}=0.
\]
At \(k=0\), it contains an acoustic mode plus an optic micro-rotation mode
\[
\omega_0^2=\frac{4\kappa_c}{I_c}.
\]
Using this branch would require deriving the micro-rotation variables and classifying the gap scale. Hiding the extra mode by
choosing a large gap is a C6 tuning unless independently motivated.

## Provenance Table

| ingredient | Stage-2 status | provenance |
| --- | --- | --- |
| GP/NLS mean-field fluidity | mean-field has zero shear | independently known, but it tests the wrong object if used alone |
| cohesive substructure / healing length / tautness | motivated physical picture | not enough to derive a shear modulus |
| persistent in-plane neighbor memory | required to select elastic-network limit | unspecified |
| affine/network free energy | required to derive \(\mu_{\rm br}>0\) | unspecified |
| simple-fluid / soap-film absence of neighbor memory | would derive \(\mu_{\rm br}=0\) | not established for the substructure |
| \(K_{\rm br}\) compression modulus | coefficient in the test family | phenomenological |
| \(\mu_{\rm br}\) symmetric shear modulus | not derived zero or nonzero | `FAIL_UNSPECIFIED_SUBSTRUCTURE` |
| MacCullagh \(\mu_R r^2\) | contrast only | reverse-engineered unless gyrostat/spin substructure is independently supplied |
| Cosserat \(\varpi_a,\kappa_c,A_c,I_c\) | contrast only | extra variables/gap not derived; conditional/tuning risk |
| \(B_\ell(w)\) finite-thickness profile | carried from Stage 0 | postulated regulator/profile |
| brane-to-bulk coupling | carried from Stage 0 | `SYMMETRY_ALLOWED_POSTULATED` |

## Stage-3 Hand-Off

Because Stage 2 lands `FAIL_UNSPECIFIED_SUBSTRUCTURE`, there is no selected constitutive law to roll into Stage 3. If the user
later explicitly accepts a postulated or newly derived Cauchy/network branch, the parametric objects to balance against the
bulk are
\[
\sigma_{ab}=K_{\rm br}\theta\delta_{ab}+2\mu_{\rm br}e_{\langle ab\rangle},
\qquad
M_{cab}=0,
\]
\[
t_b^{\rm brane}=\nu_a\sigma_{ab},
\qquad
f_b^{(4)}=B_\ell D_a\sigma_{ab}.
\]

The Stage-1 carry-forward remains unchanged: the no-leak token is a deferral, not a pass. Stage 3 would have to balance the
brane traction/couple objects against
\[
T_{na}=T_{wa}+(T_{ww}\delta_{ab}-T_{ab})D_bu_w
\]
and the mouth traction \(t_A^a\). The no-leak conditions are special, not generic, and can fail near a draining defect. This
report gives no form any unearned no-leak credit.

## Script Pointers

Primary Mathematica:

- `software/stage1_solver/tools/pathA_23_stage2_constitutive_form.wl`
- Output: `software/stage1_solver/_scratch/pathA_23_stage2_constitutive_form_mathematica.json`
- Run: `timeout 600 math -script software/stage1_solver/tools/pathA_23_stage2_constitutive_form.wl`
- Result: 32/32 checks, exit 0, token `FAIL_UNSPECIFIED_SUBSTRUCTURE`, \(\mu_{\rm br}\) decision `UNDETERMINED`.

Independent SymPy cross-check:

- `software/stage1_solver/tools/pathA_23_stage2_constitutive_form_sympy.py`
- Output: `software/stage1_solver/_scratch/pathA_23_stage2_constitutive_form_sympy.json`
- Run: `timeout 600 python3 software/stage1_solver/tools/pathA_23_stage2_constitutive_form_sympy.py`
- Result: 32/32 checks, exit 0, token `FAIL_UNSPECIFIED_SUBSTRUCTURE`, \(\mu_{\rm br}\) decision `UNDETERMINED`.

Load-bearing checks:

- the deviatoric invariant \(e_{\langle ab\rangle}e^{\langle ab\rangle}\) is present alongside compression;
- the substructure classifier can land `ZERO`, `NONZERO`, or `UNDETERMINED`;
- the actual record lands `UNDETERMINED` because neighbor memory/network free energy/frequency regime are unspecified;
- the Fourier stiffness extractor detects transverse stiffness when \(\mu_{\rm br}>0\) and zero transverse stiffness when
  \(\mu_{\rm br}=0\);
- stress, boundary traction, angular-momentum closure, and Hamiltonian positivity/ghost controls are derived for the symmetric
  branch;
- MacCullagh and Cosserat branches are retained as contrasts with the antisymmetric-stress and micro-rotation-gap obstructions.

Bookkeeping checks:

- dimensions for energy density, stress, boundary work, brane force density, finite-thickness force density, Cosserat spin
  inertia, and Cosserat curvature modulus.

What the scripts do not verify:

- a microscopic formula or numerical value for \(\mu_{\rm br}\);
- an independently specified lattice, neighbor graph, gyrostat, director, or Cosserat micro-rotation;
- Stage-3 leak closure against the bulk stress/throat profile;
- Stage-4 off-brane \(u_w\) status or full spectrum;
- Stage-5 Maxwell/gauge structure, charge, or cone locking.
