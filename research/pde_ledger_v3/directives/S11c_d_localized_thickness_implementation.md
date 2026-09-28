# Fixed-input thickness forcing and outgoing action: implementation for review

Status: prepared method/build proposal, **no scientific result or launch**.
The worker is `../_measurements/S11c_d_localized_thickness_response.py`.
This implements the [thickness-first plan](S11c_d_localized_thickness_response_plan.md)
for the two existing grade `(0,1)` cases, blocks 16 and 17. Physical inputs,
saved redundant residue-column frame, five coupled fields and units stay fixed.
The six zero end fields found in the saved inspection are not a waiver of
density or mixed-grade work. Only the two thickness cases enter this worker.

## Source and proposed transform

Restore the saved native forcing/domain bundles, reference context/units and
accepted outgoing candidate. Require exact seed/residue, zero-lift, forcing
sign, physical input and source-domain joins. In these cases the actual source
is `F01 = -L01 exp(i*k0*z) A0`; the lifted reference action is zero. No native
producer or previous end-field construction is called.

Use the inherited convention

    Fhat(k) = integral exp(-i*k*z) F(z) dz
    F(z) = integral exp(+i*k*z) Fhat(k) dk / (2*pi).

After removing the known plane character, the local source is a polynomial
`P(tanh(z/L))`. Divide the actual polynomial by `1-T^2`, and save the quotient,
remainder and both endpoint values before checking that the remainder and
endpoints are zero. No undecided zero flag counts as zero. This exposes a
finite sum of `sech(x)^2 tanh(x)^n` terms with `x=z/L`.

Write `B0(q)=pi*q/sinh(pi*q/2)` with its removable value `B0(0)=2` and `p0=1`.
The transform of `sech(x)^2 tanh(x)^n` is `B0(q)*pn(q)`, where

    p_(n+1) = (n*p_(n-1) - i*q*p_n)/(n+2).

For n=0 the previous term is zero. The recurrence comes from the derivative
of `sech(x)^2 tanh(x)^n` and integration by parts; the endpoint term vanishes.
For the base identity, the substitution `u=(1+tanh(x))/2` gives
`2*Beta(1-i*q/2,1+i*q/2)`, hence the stated B0. The scaled transform is
`L*B0(L*(k-k0))*pn(L*(k-k0))`. The worker handles only the finite polynomial
profiles it actually finds (degree bound six after removing `1-T^2`), saves
every rule application, and independently compares each used basis order
against 60-digit mpmath quadrature at dimensionless q=0,1,2. These are small
dictionary checks, not new physical-frequency or profile-response sweeps.

The saved nonlocal terms have exactly two supported forms. The worker checks
their actual bound-variable roles, order, limits and remaining free symbols:

1. `c integral dk_out dzs exp(i*k_out*z-i*k_out*zs) M(k_out) H(zs)`.
   Fourier transformation in z gives `2*pi*c*M(k)*Hhat(k)`. H includes the
   inherited plane character and the same localized polynomial dictionary.
2. `c integral dk_out dp dzs exp(i*k_out*z-i*p*zs) M(k_out,p)
   exp(i*k0*zs) hhat(k_out-p)`.
   The z and zs characters give two factors of `2*pi`; the input momentum
   is the inherited k0, and the output is
   `(2*pi)^2*c*M(k,k0)*hhat(k-k0)`. The existing hhat has its original L factor.

The original ordered integral and transformed return are saved together.
These rules propose a distributional reading tested on the explicit smooth
localized profile, rather than an assertion of absolute Fubini convergence
for the original oscillatory integrals. In the second rule the plane delta
can multiply the native multiplier only after its smoothness at real p=k0
is justified. No profile Abel regulator is present in these thickness terms.
The review must assess this reading of the inherited source convention.

## Real-axis domain and coupled action

Each transformed term has the common factor `B0(L*(k-k0))` times an amplitude.
The implementation inspects the actual amplitude before whole-matrix
simplification. It permits rational dependence on k and one source radical
`sqrt(-a-b*k^2)=i*r`, with source-derived positive a,b and
`r=sqrt(a+b*k^2)>0`. Negative half powers receive the same branch substitution.
Momentum-dependent denominator factors must reduce to rational polynomials
in r. The saved real or imaginary coefficient list of each numerator and
denominator must have one strict sign and a nonzero coefficient. Unsupported
factors or undecided signs stop with saved evidence, without a different
method or automatic retry.

Such factors are smooth on the real line. Since r has a positive lower bound,
the same-sign polynomial argument also bounds reciprocal denominators and
their derivatives by powers of momentum. Multiplication by B0 and its
removable zero-transfer continuation then makes this fixed forcing transform
smooth and rapidly decreasing. The packet includes the actual native source
expressions so that this conclusion can be checked against the implementation,
including multipliers collapsed at p=k0. It is not a global operator-domain
certification task.

Multiply the full saved 5-by-5 inverse by the complete Fhat. Preserve the
saved symmetric exclusion intervals and signed-current delta coefficients:

    V(z) = PV integral exp(i*k*z) R(k) Fhat(k) dk/(2*pi)
           + sum_j exp(i*kj*z) Dj Fhat(kj)/(2*pi).

The pole terms use the actual smooth-source value at both inherited poles,
including the removable transfer value when kj=k0. No scalar replacement,
new inversion, LU solve, roots, modes, current normalization or new radiation
condition is introduced. The proposed claim is an action of the saved
momentum-space distribution on these two specific sources. Whether the
previous separated-point candidate supports that extension is an explicit
method-review question. No pointwise coincident-position kernel value,
general diagonal distribution extension or complex-frequency retarded
equivalence is asserted.

## Selected controls and persistence

At k=0,+1/L,-1/L, evaluate the actual saved symbol/inverse and new forcing
at 80 digits, and save `P*(R*Fhat)-Fhat` and its relative entry norm. This is
a new-source contraction, not a repeated solve. Save the pole null residual
on each new source. Use declared diagnostic thresholds 1e-10 for regular
source equations/profile dictionary and 1e-8 for inherited pole contractions;
these are algebraic implementation checks, not observable accuracy claims.

Omit one actually addressed thickness native term from the transformed source,
hold the baseline fixed and require a changed forcing and spectral response
at the same small probe set. Save the term, both matrices and actual movement.
Do not substitute the old density mutation for this control. Contract actual
saved row/field/column units through all five intermediate rows.

Every restored payload and new operation has complete journal inputs and
returns. The whole case also has a containing operation receipt, so a timeout
between smaller operations retains its full source operands. Guards examine
profile/domain returns only after persistence. The main failure handler saves
the incomplete operation stack, completed index, artifacts and source hashes.

Normal budget: one shared guarded worker around the existing normalization
supervisor, 900s outer/840s native, 2GiB, zero swap, one CPU, nice15, 32 tasks,
one native thread. Cost is unmeasured; do not infer a duration exception from
earlier work. No automatic retry, overlap, unguarded fallback or polling.
Before any future launch, arm the existing local hook for session
`01a0e01b-ef84-7192-817f-584cda5d339b` and pin the actual helper/input/review
disposition hashes in a fresh gate.

Success produces two source-specific Fourier/PV response candidates and
their controls, with integrals unevaluated and saved-output acceptance still
pending. An unsupported source/domain rule stops before that claim. Full
response/Green/FORM/A11/A12, density/mixed grades, physical current
normalization and radiating coverage remain open. Practical toy-model scope
still governs; optional wording is not grounds for a repeated review cycle.
