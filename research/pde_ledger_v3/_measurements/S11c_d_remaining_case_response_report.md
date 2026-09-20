# Four-case finite and continuum responses

All four approved cases now have complete finite and independent-grade
responses, each with 645 unknowns and four incident columns. The three new
cases use their actual 70/80/70 rows, 140/185/162 native terms, 25/35/25 source
amplitudes and case-specific end/current maps. The baseline is copied from its
accepted result. No quadrature, mode or baseline response was repeated.

All systems have full rank. Maximum finite scaled equation residual is
8.04e-16; maximum continuum residual is 3.80e-14 including the baseline.
The largest new independent coefficient difference is 1.78e-12. Actual
baseline, mixed inverse/forcing terms, open/closed amplitudes, current
denominators and phase/normalization maps are retained.

Same-density anchoring comparisons use identical end coordinates:

| Density rule | Finite amplitude change | Finite current-ratio change | Evaluated continuum amplitude change |
| --- | ---: | ---: | ---: |
| RHO4_CONSTANT | 7.44e-7 | 3.36e-7 | 6.87e-7 |
| RHOBR_CONSTANT | 5.99e-8 | 1.19e-7 | 6.02e-8 |

These are below the declared amplitude 1e-4/current 1e-6 absolute resolution.
No resolved tiny reflection, loss or anchoring effect is inferred. Full common
grid field differences are 4.95e-4 and 1.55e-6 in the inherited coefficient
frame; these are retained separately from channel amplitudes. Cross-density
S matrices were not subtracted. New finite total-current ratios lie between
0.9999995945 and 1.0000001045.

Formal coefficient-polynomial remainder maxima decrease by about four on
halving both parameters. At (eta,sigma)=(0.01,0.001), the new LAB/RHOBR and
MAT/RHOBR field-coefficient differences are 3.22e-5; MAT/RHO4 is 8.82e-3.
These are Taylor truncation diagnostics, not physical observable error bounds
or accuracy for omitted parent pure-second-order terms.

The original 687.49-second run completed all numerical packets and all three
individual emission/replay tails. Raw concatenation then failed because each
part starts its own codec reference table. Finish01 completed all four saved
numerical validations in 32.62 seconds, then found that baseline checks remain
in the original accepted producer. Finish02 joined those checks to all six
baseline artifacts, reused every completed numerical validation, and globally
re-encoded the original decoded payloads. It completed in 62.54 seconds with
exit zero, empty stderr and checks/stdout byte identity. Neither repair repeated
a solve, integration or individual emission.

Final acceptance verifies 95 current/frozen sources, 158 inputs, 275 unchanged
copies, all 269 original artifacts and 188 final artifacts. All 2,024 decoded
payloads, four local source-line indices and 182,924 metadata paths are retained;
1,004 export keys are unique. Publication fa684a3a is annex-verified:

- Output: scripts/out/S11c_d_remaining_case_response.out, 16,310,842 bytes.
- MD5E key: MD5E-s16310842--14c206142892fae92f3f6c5ca247e8f9.out.
- SHA256: c30b166b9dfa31123cdd753a430b2d14c0ce978e2eaffda605151ac74786e1c3.

The numerical settings remain 129 coefficients per field, interval/source 64,
momentum 4, profile 14, regulator 0.1, momentum orders 16/4/4 and source/profile
256/512. Positive regulator and approximate modal boundaries remain explicit.
Case current bookkeeping, relevant practical controls, bounded case searches
and final all-case engine/exports remain work. The baseline 32-point search is
closed without a resolved candidate, not a certified empty spectrum. No broad
quadrature campaign or further angular doubling is queued.
