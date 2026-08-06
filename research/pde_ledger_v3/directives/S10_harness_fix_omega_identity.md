# One bug: a declared symbol identity is missing, and it makes S9 report a false DISAGREE

Fix this one thing. Do not refactor anything else, and do not commit.

## The bug

`reduction/checks_S9.yaml` declares `factored_determinant`, comparing
`WL_S9_DETERMINANT` against `PY_S9_MAIN_DET_M_FACTORED`. The rebuilt harness reports `DISAGREE`. The old
harness reported `AGREE`, and the old harness was right.

The two emitted operands:

```
WL:  omega**2 * rhoBr  * (-kx**2*muR  - ky**2*muR  - kz**2*muR  + omega**2*rhoBr )**2
PY:  omega2   * rho_br * ( kx**2*mu_R + ky**2*mu_R + kz**2*mu_R - omega2*rho_br  )**2
```

The inner factor differs by an overall sign and is squared, so that is not the difference. The difference
is that this engine carries the squared frequency as a single atomic symbol `omega2`, while the other
writes `omega**2`.

Verified, with the declared spelling map applied to the WL side:

```
spelling transliteration only            -> residual != 0        (the DISAGREE)
spelling transliteration + omega2=omega**2 -> residual == 0
control: same substitution, PY operand doubled -> residual != 0  (so the check is not vacuous)
```

## Why it regressed

The previous harness held `omega2 -> omega**2` and `omegaSquared -> omega**2` in a hardcoded
`_SYMBOL_ALIASES` table inside the comparator. The rebuild instruction said an identity between
differently-named quantities "belongs in the registry as a declared relation, or it does not hold", and
the entries were removed without determining which. It holds. It needs declaring.

## Scope, measured — do not over-generalise

```
S9    wl: omega=18   omega2=0    omegaSquared=0
      py: omega=9    omega2=27   omegaSquared=0
S10   wl: omegaSquared=28   omega2=0   omega=0
      py: omegaSquared=82   omega2=0   omega=0
```

Only S9 has a cross-engine mismatch on this symbol. In S10 both engines already spell it `omegaSquared`,
and `checks_S10.yaml:137` already declares `omegaSquared: [0, -2, 0]` as a primitive dimension. So S10
needs no identity for cross-engine comparison. Do not add one there to be safe; an unnecessary
substitution is a way to manufacture agreement.

## What to build

1. **A declared identity section in the config**, at the same level as `symbol_naming`, in which an engine
   may declare that an atomic symbol it emits stands for an expression. Each entry carries the engine, the
   symbol, the expression, and a reason. Applied only to the declaring engine.

2. **Apply it in the cross-engine path** where the naming map is already applied
   (`engine_output_checks.py:1097`), after shape selection and action normalisation, and record what was
   applied per row exactly as `naming_applied` is recorded.

3. **Report it on the verdict line**, beside `naming=[...]`, so a reader can see that a comparison assumed
   an identity. An identity that is applied but not visible in the report is the defect this replaces.

4. **Declare the S9 entry** and nothing else: this engine's `omega2` stands for `omega**2`, because the S9
   spectrum is solved as a polynomial in the squared frequency and that engine carries it as one variable.

5. **An identity must not be able to reconcile operands that genuinely differ.** Add a test that a
   declared identity applied to a deliberately altered operand still reports `DISAGREE`.

## Constraints

- Do not change `symbolic_equal`'s equality semantics.
- Do not add an identity to `checks_S10.yaml`.
- Do not restore `_SYMBOL_ALIASES` or any hardcoded table in the comparator.
- A physics disagreement must not exit non-zero.
- No script over 600 seconds. Do not launch Mathematica; `/tmp/s9_wl.txt` exists and the suite needs it.

## Report back — under 15 lines

1. The config entry you declared and the code path that applies it.
2. Literal stdout of both configs' cross-engine summary lines, before and after.
3. The result of the non-vacuity test in item 5.
4. Anything you measured that contradicts the above.

Stop there. Do not commit.
