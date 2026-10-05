# Tail/Abel and inverse-stability fidelity boundary

Fidelity record for [ANALYTIC_ERROR_COVERAGE.md](ANALYTIC_ERROR_COVERAGE.md).
The fresh recorded suite passed, including all mathematical controls. Claude
and Grok independently cleared statement fidelity; bounded T1–T4 is complete.
See [ANALYTIC_ERROR_FIDELITY_REVIEW.md](ANALYTIC_ERROR_FIDELITY_REVIEW.md). Focused builds were not accepted as recorded reuse evidence.

## Mathematical object

`Tail` uses actual Bochner integrals of complex amplitudes on the real line.
Its top-level truncation theorem requires integrable g, an almost-everywhere
strongly measurable bounded b, and a measurable kept set s. These imply
integrability of b*g. The primitive omitted-integral norm inequality alone is
not an integrability certificate; the actual truncation theorem supplies that
premise. `tailMass` is the integral of norm g over the complement, and
`firstMoment` is the integral of abs(x)*norm(g). No real-part projection occurs.

`Abel` defines the actual regulator as `exp(-a*abs x)`. Its combined theorem
requires a>=0 and a finite first absolute moment, and proves

```
norm (integral b*g - integral_on s (exp(-a*abs x) * b*g))
  <= B * (tailMass + a * firstMoment).
```

The proof bounds the difference of the full integrals and the tail of the
regulated integral separately. It does not take a pointwise limit in momentum
space or discard a step contribution. Zero a and zero amplitudes are included.
The set `{x | abs x <= R}` is measurable for every real R (including empty
negative-R sets); the physical cutoff uses R>=0. Measures are general, so
Dirac measures provide exact nonzero controls. The physical measure is Lebesgue.

`Stability` starts from a continuous linear equivalence A between normed spaces
and a bounded perturbation E. Completeness of the domain suffices for the
Neumann construction; a continuous linear equivalence then transports the
Banach property to the codomain. `perturbEquiv` uses Mathlib's `Units.oneSub`
to construct a genuine inverse. `perturbEquiv_eq` identifies its forward map
with A+E. Intermediate estimates for a supplied equivalent B are joined to
that existence theorem by `exists_controlled_inverse`; invertibility of A+E
is not silently assumed. The norms kappa and epsilon are upper bounds, and
their strict product margin gives the positive denominator `1-kappa*epsilon`.

The inverse-difference identity, source perturbation identity and fixed bounded
observation map are proved on the full spaces. They do not certify a particular
real-frequency outgoing realization or an approximate channel extraction map.

## Compact source connection

`S11_lean_analytic_error_source_check.py` extracts only the original Fourier
normalization and Abel half-line construction from `EdgeReduction.__init__`.
It does not import the production module, construct the closed operator or
read/replace pinned exports. It checks the original operands exactly with
SymPy, including the native Gaussian computation of Fourier mass 2*pi:

```
negative half = 1/(a-i*s), positive half = 1/(a+i*s)
constant      = 2*a/(a^2+s^2)
step even     = a/(a^2+s^2), step odd = -i*s/(a^2+s^2).
```

With s=L_W*(k-k'), the measure factor L_W/(2*pi) gives the Poisson kernel of
width a/L_W and the step kernel (P_h-i Q_h)/2. The supplied development input
has L_W=10. Wrong width multiplication, PV sign, missing step half and discarded
constant are distinct nonzero symbolic controls. The half-line operands and Fourier mass are executed from selected native
source. The checker itself supplies the Jacobian/substitution and physical
measure mapping; correspondence of that mapping is verified by inspection of
`EdgeReduction.hat` and `prescribe` (native lines 1660–1680, 1730–1756). Both
reviewers also inspected these methods. In particular, the wrong-width control
tests this explicit map, not execution of those two full native methods. This
is a compact tested/inspected translation, not a kernel-certified transform
implementation.
The raw native `delta_mass` field remains an unevaluated Meijer-G expression
in the report; it is not accepted as a proved Abel mass integral. The checked
Fourier mass comes from the separate original Gaussian computation above.

For application to the reduced operator, partition contributions by treatment:

| Contribution | Required estimate/identification | Current status |
|---|---|---|
| Local differential terms | Graph/Sobolev domain and coefficient-error bounds | Not supplied by this theorem |
| Ordinary profile and regular nonlocal terms | Integrable folded amplitudes, retained/complement domains, uniform tail constants; justified integral order | Conditional T1 estimate only |
| Constant/step distributional terms | Fixed origin, native transform/measure, Fourier pairing with actual trial/test spaces, uniform first moments | Native scalar operands identified; application premises open |
| Solve and channels | Common outgoing domain, operator-norm error, inverse margin, channel extraction/current bounds and their errors | Conditional T2 estimate only |

The categories are an application-obligation inventory, not a claim that every
native term already satisfies an estimate. There is no certified numerical
epsilon, parent-theory bound or class-wide scattering convergence result here.
The actual selected profile/parameter source hashes are recorded in the source
check. Shared-source drift must be reconciled; never overwrite the other
session's work to match an older hash.

## Controls and verification

The finite suite has twelve paired false statements and sixteen positive
executions representing fifteen distinct statements: `inverse_difference_positive`
and `operator_error_omission_positive` share the correct bound 1.
They exercise actual nonzero complex/real integral witnesses, omitted tails
or moments, regulator sign/absolute position/physical width, lost inverse
conditioning, omitted source/operator error, observation amplification and the
singular critical margin. Every accepted false statement must reduce to False
in its named theorem; syntax, import, timeout and environment errors are not
mathematical rejections. Passing zero-regulator, zero-amplitude, integrability
and strict-margin cases expose nonempty admissible domains.

These are paired statement mutations with concrete counterexamples, not
automated source replacements in the canonical general theorems. In particular,
the inverse controls use exact scalar arithmetic for A=1, E=-1/2, source=1
and source error=1/4; the critical case uses E=-1. They detect the selected
lost factors or omitted terms, without claiming sensitivity to every possible
edit of a general theorem. The integral controls use actual Dirac integrals;
they do not establish physical Lebesgue moments or a solver's operator bounds.

Focused proof repairs so far concern measurability of the explicitly defined
regulator, finite-measure integral simplification, coercion reduction, lemma
argument order and unnecessary section instances. Claims and normalizations
are as stated above. No old proof or object has been rebuilt for this increment.
The recorded suite built four modules plus the 41-root axiom audit freshly,
and ran every control with one worker. All twelve mutants had exactly one
mathematical `False` diagnostic; all sixteen positives passed. There were no
timeouts or instrument errors in this recorded run. External package pins, clean
tracked sources and direct imported Mathlib source/object hashes are checked;
local transitive source and input/output object hashes guard any later reuse.
The existing pinned transitive Mathlib cache remains the dependency trust
baseline; the full library was not rebuilt. All nineteen preceding D4B objects
and the recorded historical proof bytes remain unchanged. Current lakefile and
README additions belong to T1–T4, not a revision of the D4B evidence.

| Contract | Principal declarations | Remaining application premises |
|---|---|---|
| T1 | `truncation_error_bound`, `abel_pairing_bound`, `abel_truncation_bound` | Uniform folded-amplitude tails and first moments on the actual trial/test spaces |
| T2 | `perturbEquiv_eq`, `exists_controlled_inverse`, `solution_error_bound`, `observation_error_bound` | Outgoing/graph domain, common operator norm, inverse margin and channel errors |
| T3 | `abelFactor`, `physicalWidth`, `approved_width`; compact native source check | Legitimate distributional pairing and composition for the complete reduced operator |
| T4 | `Controls.lean`, recorded paired controls and audit root | Both independent reviews CLEAR; physical application premises remain open |

Stop after the contract and reviews clear. Stronger Fourier-space rates,
limiting absorption, operator-specific global estimates and the other two
analytic questions remain separate work.

For a real Fourier transfer, its unit-modulus phase may be included in b
without enlarging B (`real_phase_norm`). Complex transfer requires a separate
compatible bound and is not automatically covered by that unit-modulus claim.
Optional wrappers and stronger rates were not added under the stopping rule.
