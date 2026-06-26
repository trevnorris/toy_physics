# PathA-33 quadrupole normalization result

Computed verdict: `QUAD_CALIBRATED` (`QUAD_CALIBRATED`).

The earned part is the outgoing l=2 fingerprint shape, the squared-denominator prefactor algebra, the a^-5 target scaling, the Gamma5/chi_Q equivalence, and the restored-units dimensional closure. The magnitude remains calibrated: G and the PN 2/5 input keep the assembled 54/5 on the `external_bridge_input` rung.

## Derived fingerprint

- SymPy: u2=`a**2/(9*c_s**2)`, u4=`4*a**4/(81*c_s**4)`, v5=`a**5/(27*c_s**5)`.
- Mathematica: u2=`a^2/(9*cs^2)`, u4=`(4*a^4)/(81*cs^4)`, v5=`a^5/(27*cs^5)`.
- Derived chi_Q: `1`; incoming gives `-1`.

## Prefactor algebra

- P0=`N0/D0`.
- P2=`(D0*N2 - 2*D2*N0)/D0**2`.
- P4=`(D0**2*N4 - 2*D0*D2*N2 - 2*D0*D4*N0 + 3*D2**2*N0)/D0**3`.
- N/D self-check: plain P2=`(D0*N2 - D2*N0)/D0**2` versus correct P2=`(D0*N2 - 2*D2*N0)/D0**2`; the missing term is `-D2*N0` versus `-2*D2*N0`.

## Dimensional result

- Mu-free gate: `[P0_raw]` = `T^2`, `[(c_s/a)^2]` = `T^-2`, `[P0_phys]` = `1`; dimensional_ok = `True`.
- `mu_hat0` diagnostic: `non-able-to-fail (mu_hat0 free carrier)`; `[mu_hat0]` = `L^-1 T^-1 M^-1/2`, LHS/RHS = `L^-2 T^-2 M^-1` / `L^-2 T^-2 M^-1`, diagnostic pass = `True`.
- Section 3d drop-normalization probe verdict: `FAIL_DIMENSIONAL`.
- Section 3d' corrupt-port-dimension probe verdict: `FAIL_DIMENSIONAL`.

## Provenance

- Decomposition: `54/5 = 2*27/5`; earned factor class = `derived_in_gate`, calibrated factor class = `external_bridge_input`.
- Assembled 54/5 class: `external_bridge_input`.

## Probe verdicts

- 3a incoming: `FAIL_FINGERPRINT`; standing: `FAIL_FINGERPRINT`.
- 3b imposed dissipation: `FAIL_NOT_OUTGOING`.
- 3c wrong scaling: `FAIL_SCALING`.
- 3d dimensional break: `FAIL_DIMENSIONAL`.
- 3d' corrupt port dimension: `FAIL_DIMENSIONAL`.
- 3e equivalence break: `FAIL_EQUIVALENCE`.
- 3f partition mislabels: `FAIL_PROVENANCE_PARTITION`, `FAIL_PROVENANCE_PARTITION`.
- 3g wrong prefactor object: `FAIL_PREFACTOR_ALGEBRA`.

## Engine agreement

- Status: `pass`.
- Max symbolic delta: `0.0` (tol `1e-10`).
- Max numeric delta: `3.642919299551295e-16` (tol `1e-08`).

Deferred: the numerical branch data `(D_n,N_n)`, port scalars, actual branch a-scaling, cross-l reconciliation, and derivation of G/the magnitude remain outside Gate 4.
