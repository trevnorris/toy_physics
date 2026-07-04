THROAT_ELECTRIC_LOCALIZED_COULOMB

# pathA_38 Throat-Body Electric Localization

Computed headline: `THROAT_ELECTRIC_LOCALIZED_COULOMB`.

This run derives the radial branch from the transverse wall spectrum. The zero mode is first solved as an eigenfunction of `O f=m^2 K_parallel f`; only then is its radial factor solved.
The computed norm is `N0=8/(3*ell)` from the integral density `2/(ell**2*cosh(w/ell)**4)`.
The solved static/dynamic exponents are `p_static=2` and `p_dynamic=2`.

The compact source projections are `q_h(+)=2*QE*tanh(b/ell)/b` and `q_h(-)=-2*QE*tanh(b/ell)/b` from independent integrals.
The pure even gravity overlap is computed as `0`; the orthogonalized no-monopole source gives `0`.
The projected Green kernel is `3*QE**2*ell*tanh(b/ell)**2/(8*pi*R*b**2)`.
The YAML records the static zero, shape, and relative bound-mode kernels; the nonzero massive bound terms are Yukawa.

The sign matrix is `U++=3*QE**2*ell*tanh(b/ell)**2/(8*pi*R*b**2)` and `U+-=-3*QE**2*ell*tanh(b/ell)**2/(8*pi*R*b**2)`.

Able-to-fail classifier self-test:
- main branch -> `THROAT_ELECTRIC_LOCALIZED_COULOMB`
- delocalized ablation -> `FAIL_DELOCALIZED_BULK_1_OVER_R3` with `p=3`
- ghost ablation -> `FAIL_GHOST_INSTABILITY` with `G0_sign=-1`
- no-monopole ablation -> `FAIL_NO_MONOPOLE`
- pinned-branon ablation -> `FAIL_PINNED_BRANON`
- Yukawa ablation -> `FAIL_YUKAWA`

Provenance split:
- derived symbolic: transverse spectrum, normalizability, radial exponents, source projections, interaction signs
- computed ablations: delocalized, ghost, no-monopole, pinned-branon, Yukawa, gravity-mixing
- dimensional firewall: symbolic units plus negative controls
- calibrated/deferred: `Q_E`, full nonlinear throat source compactness, nonlinear operator parity mixing

Engine agreement status: `ENGINE_AGREE`.

Run commands:
- `timeout 600 python3 software/stage1_solver/tools/pathA_38_throat_body_electric_sympy.py` -> exit `0`
- `timeout 600 math -script software/stage1_solver/tools/pathA_38_throat_body_electric.wl` -> exit `0`
