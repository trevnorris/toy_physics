# Saved omega=1 current balance: reporting scope

Source/JSON assessment only, 2026-09-29. No scientific payload was restored and
no new observable was calculated. This note accompanies the proposed execution-4
result report; it is not an expected result for that premise instrument.

The four final-baseline outgoing/incident ratios quoted in
`S11c_d_saved_balance_numeric_first_assessment.md` are **four incident columns
of the baseline solve**, not one ratio for each physical case. Separate saved
responses cover the other cases. The quoted baseline values straddle one and
lie within the declared absolute current reporting resolution of `1e-6`.
The appropriate statement is: **no loss resolved at that resolution in the
saved finite calculation at its stated settings**. This is neither proof of
zero physical loss nor a certified upper bound.

Two limitations must be attributed to their actual calculations:

- The **retained current series** omits parent pure-second-order contributions.
  Such field/current terms can interfere with the nonzero incident baseline
  at precisely the order relevant to a leading loss. The saved coefficient-
  polynomial contrast-halving comparisons do not supply these missing terms
  or a Born/full power anchor.
- The **full finite original-contrast solve** is a different saved object.
  Do not automatically attribute the retained series' omission to that solve.
  Its current ratios still refer to a finite discretization, finite matching
  domain and positive Abel regulator, with the reported selected numerical
  sensitivities. Those comparisons are not absolute physical error bounds.

In particular, the observed change between Abel regulators `0.2` and `0.1`
does not determine the absolute regulator error or a removable additive
absorption floor. The regulator weights half-line Fourier operands; it is
distinct from physical face permeability. The homogeneous uniform controls
do not duplicate that regulated profile experiment.

No fixed-point premise result, all-grade matched-end localization, physical
loss-side power pairing, omega=3 conclusion or calibrated analog-light band
follows from these saved ratios alone. Report any new premise result separately.

Source routes: `S11c_d_finite_scattering_domain_report.md`, its checkpoint,
`S11c_d_remaining_case_response_report.md`, and the source addresses/hash
provenance already recorded in `S11c_d_saved_balance_numeric_first_assessment.md`.
