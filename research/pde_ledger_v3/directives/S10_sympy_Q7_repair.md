# Make the SymPy engine's Q7 emit the density its own action used

One attempt. If the object cannot be emitted as specified, say so and stop — a recorded "this could not be
done, here is why" is the wanted outcome in that case, not a workaround.

**File:** `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py`
Do not commit. Do not touch the Mathematica engine, its committed output, `steps/`, `paper/`, or
`REBUILD_HANDOFF.md`.

## The defect, measured

`:320` builds the action from **the package's own** stiffness form:

```python
stiffness = stiffness_density(package.stiffness, gradient)
```

`:1538` then hardcodes the form for Q7, regardless of package:

```python
curl_stiffness = sp.expand(stiffness_density("curl", gradient_symbols))
...
curl_difference = sp.expand(curl_stiffness - curl_dot)
emit ... Q7_STIFFNESS, Q7_CURL_DOT, Q7_DIFFERENCE
```

Consequences, both measured by a review leg on committed output:

- `PY_S10_*_D3_Q7_STIFFNESS` is **byte-identical across all five packages**, while the Mathematica side's
  `Q7_PACKAGE_STIFFNESS_DENSITY` tracks each package's form.
- `*_q7_difference` compares the curl density against the curl norm **for every package**, so the emitted
  residual carries no information about the package. The cross-engine row reads `py=0` against a nonzero
  WL residual.

The second point is the important one: as written, Q7 on this engine is a check that cannot fail. It asks
whether the curl density equals the curl norm and answers that question five times.

## What Q7 is for

The Mathematica side emits, per package: the density its action used, the ordinary curl norm, and their
residual. Q7 is the `D=3` bridge that tests whether a package's general-`D` antisymmetric-derivative form
reduces to the ordinary curl-squared — which is what guards the coefficient carrying the wave speed.

A package whose form is not the curl has no reason to match the curl norm. That is the content the
comparison is supposed to expose, and this engine currently cannot express it.

## What to build

**The object:** the stiffness density **this package's action uses**, with the independent gradient
symbols `g_ij` substituted in — the same density, from the same source, that `:320` passes into the
action. Do not re-derive it, and do not select it by a second lookup that could drift from `:320`.

Emit all three objects, as the shared spec requires: the package density, the ordinary curl norm, and
their residual. Emit both operands and the residual — a residual alone carries no information.

Keep the independent-symbol construction. `g_ij` are independent symbols standing for `∂_i u_j`; they are
not a `k × a` amplitude curl.

## Acceptance — report what it reads, do not iterate toward a value

1. **Regenerate this engine's committed output** and paste the `Q7_` lines for every package at `D=3`.
   Report what the residual reads per package. Do not state, target, or iterate toward any particular
   value — whatever it reads is the measurement.
2. **Provenance, by mutation:** change the package's declared stiffness form and show that the emitted
   density moves with it. An emitted density that does not respond to the form its action was built from
   is the defect this repair exists to remove.
3. **Non-vacuity:** show that the emitted residual is not identically zero by construction — that it is
   capable of being nonzero for some package. A check that cannot fail is what is being replaced.
4. The full unit suite, with a one-line cause for any failure.
5. Confirm the Mathematica engine and its committed output are untouched.

## Constraints

- Only this engine's source and its own regenerated output may change.
- Do not launch Mathematica or `wolframscript` — this machine has a two-seat licence and the WL output is
  already committed and unchanged.
- No script over 600 seconds; a timeout means reformulate, never raise the limit.
- The cross-engine tag names differ between engines for these objects (`Q7_STIFFNESS` versus
  `Q7_PACKAGE_STIFFNESS_DENSITY`). That is a known shared-spec naming gap and is **out of scope** — do not
  rename tags to force a cross-engine pairing.

## Report back — under 20 lines

1. What changed, with line numbers.
2. The `Q7_` emissions per package at `D=3`, literal.
3. The mutation result from acceptance 2.
4. Whether the residual can be nonzero, with what you did to establish it.
5. Anything you measured that contradicts this directive — including if the object cannot be emitted as
   specified, in which case stop and say why.

Stop there. Do not commit.
