RETURN_RESIDUAL_PREDICTION

# pathA_29 Brane-Bulk Return v3

Computed headline: `RETURN_RESIDUAL_PREDICTION`.

The executable family is a flat finite slab with our brane at `w=0` and an adjacent return/absorber at `w=d`.
The geometry is postulated; the return response is derived from projected continuity and the solved Helmholtz transport phase.
`T0(0)=1/(epsilon0 + 1)`, so the residual is bounded but lower order: `p_res(ell0)=1`, `p_res(ell1)=3`.

Check B was run only on the admissible DC-sink completions:
- `destructuring_absorbing`: solved compact-cell spectrum and a branch-specific 3D radial equation, derived `p=2`.
- `bloch_stack`: solved q=0 Bloch spectrum and a separate 3D radial equation, derived `p=2`.

Counterfactual guard: `multiplied the solved static zero-mode Green function by r**(-4), changing 1/r to 1/r**5` was rejected with residual `5/(pi*d*r**7)`.

The radiation/Sommerfeld boundary is recorded as `ac_check_a_only` and is not used as a Check-B branch.
The signed local source accounting gives `Z=-M0*epsilon0/(epsilon0 + 1)`; under v3 this is accounting, while `Z<0` is the drain admissibility premise.

The static-dynamic consistency control used separate traces:
- dynamic trace `a98499d0b28117d8504dcd5a31a5fec9b59f5771bb3cee98e523fac838e17e03`
- static trace `d56c35ec0aa20ef526ee9c2362964fef9cd98db8a7ef45e1645ed6f2b9a3d41b`

The mandatory no-go control constructs the anti-localizing half-line warp `mu(w)=exp(2*k_warp*w)`.
Its zero mode is non-normalizable, the continuum Green integral gives `p=3`, and the same classifier returns `RETURN_NOGO`.

pde_ledger feed: open-item #9 is not closed; the deliverable is the falsifiable residual radiation prediction tied to the drain strength. The gravity-range item passes inside the localizing flat-slab family because both DC-sink completions give `p=2`.

Downstream: the full nonlinear brane-bulk return closure remains track-3 work.
