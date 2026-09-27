# Retained-order bookkeeping P1–P4: fidelity record

Local verification and author correspondence validation passed in recorded run7.
Independent Claude/Grok reviews returned CLEAR with no blocking findings.
P1–P4 is complete at its bounded scope; work is paused at the user's request.
The result concerns supplied finite polynomials; it supplies neither a physical
scattering solution nor a parent-theory Taylor error bound.

## Object and normalization

The unchanged reviewed F1–F4 definition is
`pair B x y = sum_i conj(x_i) * (B y)_i`, and `flux B a = Re(pair B a a)`.
Every matrix entry and both coherent cross terms remain. Matrices need not be
Hermitian and coefficients need not be positive. Real-valued `flux` is defined
by taking a real part; identifying it with an actual physical current requires
the application's reality and normalization premises.

Shared physics §3c supplies the amplitude rectangle, with epsilon factored out:

```
a(eta,sigma) = a00 + eta*a10 + sigma*a01 + eta*sigma*a11
Delta       = eta*a10 + sigma*a01 + eta*sigma*a11
eta=t, sigma=r*t:
a0=a00, a1=a10+r*a01, a2=r*a11.
```

Here r is an arbitrary supplied real number. Its physical identification with
Wbar0/L_W, the admissible physical values and units are application premises;
Lean does not choose or certify a physical width. Delta includes the zero-jet
contrast term. Specializing to one path can lose coefficient information and
does not identify the two original grades. Real scalar action is represented
by explicit complex casts; t and r remain real in the theorems.

For supplied `a(t)=a0+t*a1+t^2*a2` and `B(t)=B0+t*B1+t^2*B2`, the theorem proves
the actual contraction equals `sum_d t^d*q_d`, through degree six. The seven
definitions contain all 27 ordered contractions
`q_d = sum_(i+j+k=d) pair B_j a_i a_k`, with i,j,k in {0,1,2}.
In particular,

```
q2 = B0[a0,a2] + B1[a0,a1] + B2[a0,a0]
   + B0[a1,a1] + B1[a1,a0] + B0[a2,a0].
```

`retained_with_remainder` groups degrees three to six into the exact expression
`t^3*(q3+t*q4+t^2*q5+t^3*q6)`. `real_flux_expansion` takes its real part.
This remainder belongs to the supplied polynomial model. It does not bound
omitted amplitude/current coefficients of a parent theory. In particular, the
ellipsis in the shared physical current expansion is not silently set to zero
in a theorem about the parent current: the formal `current` is explicitly a
degree-two polynomial.

## Denominators and observables

Over any field, nonzero j0 gives the unique coefficients
`c0=n0/j0`, `c1=(n1-j1*c0)/j0`, `c2=(n2-j1*c1-j2*c0)/j0` satisfying the three
coefficient equations. The actual polynomial residual is

```
(j0+t*j1+t^2*j2)*(c0+t*c1+t^2*c2) - (n0+t*n1+t^2*n2)
 = t^3*(j1*c2+j2*c1) + t^4*j2*c2.
```

No analytic expansion, interval of nonvanishing denominator or estimate for a
physical quotient follows automatically. At j0=0,n0 nonzero the leading
equation is impossible. At j0=n0=0 there is no general uniqueness claim: the
all-zero leading equation admits distinct constants, as explicitly witnessed.
Higher-order singular quotients are outside this contract.

For real epsilon, amplitude scaling multiplies flux by epsilon squared. The
algebraic cancellation theorem uses nonzero epsilon; a *defined physical
fraction* additionally requires a nonzero incident denominator. That second
condition appears in `scaled_denominator_nonzero`. Totalized field division in
Lean is not a substitute for an application's denominator check.

`subtracted_total` retains both cross terms between baseline and induced
amplitudes. Its iff corollary identifies exactly when their real sum vanishes,
so total-minus-baseline then equals the induced quadratic form. With a0=0,
q0=q1=0 and q2=B0[a1,a1], independent of a2,B1,B2. With nonzero baseline an
additional second-order amplitude can alter q2: the scalar witness B0=1,a0=1,
a1=0,a2=3 gives q2=6. This witness varies the retained a2; it does not construct
an omitted pure eta² or sigma² parent-theory term. Neither that witness nor the
baseline-free specialization computes the actual physical baseline or licenses
a strong-edge extrapolation.

## Compact connection to the native calculation

The fidelity anchors are `directives/S11c_d_SHARED_PHYSICS.md` §3c–§3d,
`_measurements/S11c_d_continuum_currents.py` and
`_measurements/S11c_d_continuum_response.py`. The bounded instrument extracts
only the original current `multiply`, `quadratic`, `subtract`,
`amplitude_bookkeeping`, `quotient`, and response `adjoint`, `lambda_series`
function ASTs. No production module imports, driver, saved scientific operand
or S11c report generation ran.

Native `current.multiply` is unrestricted convolution. It must not be replaced
with the engine's componentwise-cutoff product when testing higher grades.
This describes the downstream `current.quadratic` operation. The native current
coefficients supplied by upstream `gram` already use the engine's rectangular
cutoff; Lean accepts B0,B1,B2 as supplied and does not derive that upstream current.
Native `lambda_series` combines supplied grade (p,q) into degree p+q with r^q.
Exact small fixtures test component reconstruction, the path map, all seven
degrees, full two-channel complex off-diagonal currents, actual evaluation and
the higher-degree remainder against a separate triple-loop contraction.
The scalar example a=(1,2,3),B=(1,4,5) has coefficients
`(1,8,31,72,107,96,45)`. Its second total-minus-baseline coefficient is 26,
whereas the induced coefficient is 4.
Lean applies the path before contraction; the native fixture contracts before
applying the path. Their agreement is tested on the fixture, not separately
formalized as a general commutation theorem. The native two-channel current
fixtures are Hermitian; the Lean identities permit arbitrary complex matrices.

The native quotient uses each incident column's diagonal after the full
quadratic contraction. It is not a matrix inverse and does not drop
off-diagonal current entries inside that contraction. A two-column fixture
gives c0=(1,1),c1=(1,1),c2=(1,3), with zero recurrence residuals; a separate
fixture retains the complex value i/2. The native helper has no zero guard.
The Lean complex-field positive uses real numerals; the genuinely non-real
quotient witness is the native i/2 fixture.
An isolated input with zero leading denominator produced nonfinite output and
is recorded solely as an invalid-domain witness. No saved physical denominator
or physical output defect was established.

The native report records Python 3.10.12, NumPy 1.26.4, original source and
selected AST hashes, thirty identities/examples, eleven wrong-formula controls
and one invalid-domain witness. Fixtures use binary-exact small numbers and
exact array equality. These are targeted translation tests outside Lean's
kernel, not a proof of arbitrary NumPy floating-point computation.

## Verification and repair disposition

Run7 built all seven local objects fresh in an isolated directory, including
the unchanged F1–F4 Current import. It passed 33 standard-axiom audits, sixteen
paired mathematical rejections and twenty positives. Every rejected control
has exactly one unsolved `False` in `contract_control` and no other diagnostic.
Three q2 omission pairs repeat a positive, leaving eighteen distinct positive
statements. These controls use canonical witnesses/arithmetic; they are not
independent derivations or canonical source replacements.
In particular, `path_ratio` directly checks a rectangle value. The proved
`rectangle_path` identity and native `wrong_mixed_ratio_power` control establish
the disclosed path connection; the Lean pair is not a mutated path polynomial.

Earlier attempts exposed an absent real scalar-action instance, an overly
broad import's internal interpreter memory limit, pointwise conjugation/cast
normalization, and a propositional simplification omission. Repairs narrowed
imports and changed representation/proof code without changing the mathematical
domain, coefficients or claims. Pointwise star identities, explicit real casts
and direct cancellation with j0 nonzero avoid broad normalization. The final
repair added proved logical rewrites only to the zero-model paired control;
all canonical sources and statements remained unchanged. No compiler, linter,
resource or unreduced-goal failure counted as a mathematical rejection.
The failed attempts and pre-repair hash validations remain in scratch records.

The guarded run used one CPU, 2 GiB total memory, no swap and 32 task slots;
observed peak was 268177408 bytes, with no recorded max/OOM/swap events. Lean
used one worker, a 2048 MiB internal allocator limit, strict warnings and a
180-second process timeout with whole-group cleanup. Native run1 evidence was
hash-revalidated rather than rerun. Fifteen clean pinned package repositories,
six direct Mathlib source/object pairs, all local source/log/object hashes,
567 historical files and 236 preserved objects matched. The transitive external
cache remains the pinned baseline, not a fresh source rebuild of all Mathlib.

Author validation: `_measurements/S11_lean_bookkeeping_validation.json`.
Verification summary: `SCATTERING_BOOKKEEPING_VERIFICATION.txt`.
Claude and Grok independently cleared fixed packet
`5c0a6c321424ae09c70e61442f73d1982c9aba49ac7d33d0503b610cb2baa9de`.
Both read the source and recorded evidence; neither reran proofs or native
checks, or independently verified the external filesystem behind the hashes.
Author read-only closure checks revalidated those bindings and preservation.
Optional findings and their dispositions are in
`SCATTERING_BOOKKEEPING_FIDELITY_REVIEW.md`. Only live documentation changed;
the approved snapshot, canonical proofs, instruments and historical evidence
remain unchanged. See `_measurements/S11_lean_bookkeeping_closure.json`.
No commit is authorized. Work is paused; no subsequent increment has started.
