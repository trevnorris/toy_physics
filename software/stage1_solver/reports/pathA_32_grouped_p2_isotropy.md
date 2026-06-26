# PathA-32 grouped-P2 isotropy result

Computed verdict: `ISOTROPY_CALIBRATED` (`ISOTROPY_CALIBRATED`).

The two engines computed the real l=2 basis, angular Gram matrix, Laplacian eigenvalues, grouped response coefficients, raw-D defects, normalized u-defects, and counterfactual verdict flips. The final rung is calibrated because the wall profile and radial/support scalars remain frozen calibration inputs rather than derived from the Gate-1 R0 support equation.

## Key computed checks

- Gram matrix equals I5: `True`.
- Computed -Delta_S2 eigenvalues: `{'20': '6', '21c': '6', '21s': '6', '22c': '6', '22s': '6'}`.
- K2 angular coefficient equals computed lambda_m: `True`.
- Isotropic raw-D defects: `{'a_D0': '0', 'b_D0': '0', 'a_D2': '0', 'b_D2': '0', 'a_D4': '0', 'b_D4': '0'}`.
- Normalized u-defects: `{'a2': '0', 'b2': '0', 'a4': '0', 'b4': '0'}`.
- Stability guard: `True`; denominator guard: `True`.

## Probe outcomes

- Pure-prefactor anisotropy: `FAIL_ANISOTROPIC_BRANCH`; raw-D moves, normalized defects stay zero = `True`.
- Sector-selective anisotropy: `FAIL_ANISOTROPIC_BRANCH`; raw-D and normalized u move = `True`.
- m-dependent profile: `FAIL_ANISOTROPIC_BRANCH`.
- Degenerate beta2=0: `FAIL_STABILITY`.
- Wrong eigenvalue coefficient: `FAIL_NOT_COVARIANT`.
- Singular denominator guard probe: `FAIL_SINGULAR_RESPONSE`.
- Tautology hash probe: `FAIL_TAUTOLOGICAL`.
- Static response probe: `FAIL_STATIC_RESPONSE`.

## Fixed probe self-ablations

- `singular_denominator`: with mutation `FAIL_SINGULAR_RESPONSE`, without mutation `ISOTROPY_CALIBRATED`, fail suppressed = `True`.
- `wrong_eigenvalue`: with mutation `FAIL_NOT_COVARIANT`, without mutation `ISOTROPY_CALIBRATED`, fail suppressed = `True`.
- `degenerate_beta_zero`: with mutation `FAIL_STABILITY`, without mutation `ISOTROPY_CALIBRATED`, fail suppressed = `True`.
- `tautology_hash_collision`: with mutation `FAIL_TAUTOLOGICAL`, without mutation `ISOTROPY_CALIBRATED`, fail suppressed = `True`.
- `static_drop_inertia`: with mutation `FAIL_STATIC_RESPONSE`, without mutation `ISOTROPY_CALIBRATED`, fail suppressed = `True`.

## Able-to-fail aggregate

- Computed probe gate flags: `{'pure_prefactor_anisotropy': True, 'sector_selective_anisotropy': True, 'm_dependent_profile': True, 'degenerate_beta_zero': True, 'wrong_eigenvalue': True, 'singular_denominator': True, 'tautology_hash_collision': True, 'static_drop_inertia': True}`.
- Neutering any one probe flips aggregate false: `True`.

## Engine agreement

- Status: `pass`.
- Max symbolic delta: `0.0` with tolerance `1e-10`.
- Max numeric delta: `2.220446049250313e-16` with tolerance `1e-08`.
- Per-lane `D_A,n` max numeric delta: `0.0`.
- Per-lane `D_A,n` deltas: `{'20.D0': 0.0, '20.D2': 0.0, '20.D4': 0.0, '21.D0': 0.0, '21.D2': 0.0, '21.D4': 0.0, '22.D0': 0.0, '22.D2': 0.0, '22.D4': 0.0}`.

## Input partition

- Derived inputs: `['explicit real l=2 harmonics', '5x5 angular Gram from S2 integrals', 'per-harmonic -Delta_S2 eigenvalues', 'K2 angular coefficient selected from computed lambda_m', 'ungrouped channel angular self-overlaps', 'raw-D and normalized defect algebra from assembled lanes', 'counterfactual probe verdicts']`.
- Calibration inputs: `['R0(w) linearized isotropic reference', 'beta2(w) radial profile', 'mu_eta(w)', 'T_w(w)', 'K_eta(w)', 'T_Omega(w)', 'Mtilde radial mass scalar', 'Ktilde radial stiffness scalar excluding angular T_Omega', 'TomegaTilde angular radial scalar', 'B0tilde,B2tilde,B4tilde support scalars', 'Z0tilde,Z2tilde,Z4tilde mixed/Maxwell scalars', 'physical calibration window for positivity and denominator guards', 'Gate-1 D/N boundary provenance']`.

Deferred: the 54/5 quadrupole normalization, outgoing odd N coefficients, and solved nonlinear branch data remain Gate 4/5/6 work.
