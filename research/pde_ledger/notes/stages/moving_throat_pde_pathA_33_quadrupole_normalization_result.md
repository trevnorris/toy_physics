# PathA-33 quadrupole normalization feed note

- Verdict: `QUAD_CALIBRATED` (`QUAD_CALIBRATED`).
- DtN fingerprint: u2=`a**2/(9*c_s**2)`, u4=`4*a**4/(81*c_s**4)`, v5=`a**5/(27*c_s**5)`, chi_Q=`1`.
- Prefactor object: `P(omega)=D0*N(omega)/Dcons(omega)^2`; plain `N/D` fails the P2 factor-of-two self-check.
- Dimensions: mu-free gate `[P0_raw]=T^2`, `[(c_s/a)^2]=T^-2`, `[P0_phys]=1`, dimensional_ok=`True`; drop-normalization probe=`FAIL_DIMENSIONAL`, corrupt-port-dimension probe=`FAIL_DIMENSIONAL`.
- Diagnostic only: `[mu_hat0]=L^-1 T^-1 M^-1/2`, non-able-to-fail (mu_hat0 free carrier).
- Provenance: fingerprint 27 is `derived_in_gate`; PN 2/5, G, and the assembled 54/5 magnitude are `external_bridge_input`.
- Engine agreement: `pass`, max symbolic delta `0.0`, max numeric delta `3.642919299551295e-16`.
- Deferred to Gate 6: numerical branch data `(D_n,N_n)`, port scalars, actual branch scaling, and on-solution branch selection.
