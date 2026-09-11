# S10 coefficient and sign controls — 2026-09-11

The remaining two supplied action controls now have Lean proofs in
[Scalar.lean](S10Controls/Scalar.lean), under the `S10Controls` library.

The density is defined as

```text
L_c(J) = (rho/2)|J_time|² - (c mu/2) S_curl(J).
```

`c=s>0`, `s != 1` gives `XCOEF_SCALE`; `c=-1` gives `XFORM_SIGNFLIP`.
The definition is exactly the baseline density with stiffness coefficient
`c mu`. This identity transfers the already proved density derivative,
integrated relative action, compact-test variational principle, local PDE,
and phase average. The controls enter at the action; no matrix is edited in
place or used as a substitute for the action derivation.

## Checked conclusions

For `rho>0`, `mu>0`, `k != 0`, and `c != 0`, the nonzero squared-frequency
candidate is

```text
z = (c mu/rho)|k|².
```

The spectral proof uses an arbitrary real `z`, not an already nonnegative
`omega²`. The polynomial quadratic action's amplitude derivative is checked,
and its specialization at real `omega²` equals the plane-wave density and
the operator obtained from integrated stationarity.

| Sector | Full amplitude space | N2 | N3 |
|---|---|---:|---:|
| `z=0` | span{k} | 1 | 0 |
| Nonzero candidate | k-perpendicular space | D-1 | D-1 |

Every nonzero real squared-frequency root with a nonzero amplitude is in this
table. The counts hold for every wavevector direction. Nonzero amplitudes on
the negative branch exist for D>=2; at D=1 its candidate value has zero nullity.

- Positive `c` gives a positive real cone frequency and integrated stationary
  plane waves with the baseline transverse count. If `c != 1`, its squared
  frequency differs from the baseline value.
- Negative `c`, including the sign control `c=-1`, gives a negative real
  squared-frequency root with the same transverse nullity. There is no nonzero
  stationary real cosine wave at a nonzero real frequency for this action.

The second result preserves the distinction between a negative spectral branch
and a propagating real-frequency wave. It does not construct an exponentially
growing spacetime solution, a retarded resolvent, or an interface scattering
operator. Those are separate mathematical constructions.

The coefficient-control root obeys the exact scaling ratio
`r(lambda k)/r(k)=lambda²` whenever `r(k) != 0`. The new
[anisotropic scaling proof](S10Anisotropic/Scaling.lean) establishes the same
ratio for its direction-dependent extra branch. Existing baseline scaling
covers the ordinary, full-gradient, and divergence-only frequency formulas.
No scaling ratio at the zero root is asserted.

## Verification and remaining S10 work

The full project build passes. The controls root now audits the original 32
stiffness-form declarations plus 16 scalar-control declarations; the anisotropic
root audits 48 declarations after adding its two scaling results. Together with
S9 and the baseline, these four libraries audit 151 declarations. The later
[Q6/Q7 audit](Q6_Q7_RESULT.md) adds 99, for 250 in the full build. All use only
`propext`, `Classical.choice`, and `Quot.sound`, with no proof admissions.
See [CONTROLS_VERIFICATION.txt](CONTROLS_VERIFICATION.txt).

The six S10 action families now have formal variational and mode-classification
coverage. The later Q6/Q7 audit proves expression-tree dimensions for the actions
and roots and the explicit Levi-Civita comparison. The [matrix extension](MATRIX_RESULT.md)
adds dimensions for matrices, minors, bases and residuals. Their actual CAS
expression/emission bridge and full comparator/export integration remain open.
The supplied physical action and dimension are premises throughout.

The S11c-d construction and its frozen upstream exports are outside this
increment. None of its operators, interface assumptions, continuation rules,
resolvents, or scattering claims is altered by these scalar-control proofs.
