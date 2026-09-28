**Verdict: NEEDS REVISION**

The method holds up: the transform rules, the smoothness argument and the source-specific PV/delta action are all sound for these two sources. Two implementation defects need fixing before the single no-retry launch. I only read the packet; nothing was executed, so the SymPy behaviour in S1 is my reasoning, not a test result.

## Substantive findings

**S1. Phase removal will stop both native rules at their first term.** (`S11c_d_localized_thickness_response.py`: `remove_phase` l.257–259; phases at l.309, l.325, l.328)

- `remove_phase` expands the source expression with `sp.expand_mul`, which by default also expands exponents inside `exp`. The phase it removes, `exp(-phase)`, is not expanded.
- SymPy only spreads a plain rational factor over a sum automatically (e.g. `2*(x+y)`). It does not do this for `I*(k0-kout)*zs`, `I*(k0-p)*zs` or `-I*ell*(kout-p)*xi`.
- So in rule 1, rule 2's multiplier and rule 2's profile, the combined exponent is equal to zero but not written as zero. The leftover `exp` still contains `zs` or `kout,p`.
- The guards then stop the run: `'unrecognized localized source factor'`, `'plane source leaves a multiplier only'` or `'localized fixed profile transform'`. The failure is safe, but it would use up the one bounded launch.
- The local path works because `I*k0*z` contains no sum.
- **Repair:** in `remove_phase`, use `sp.exp(-sp.expand(phase))`, or expand the exponents of the product before `powsimp`.

**S2. The term-omission control can't fail.** (`construct_case` l.484–497)

- `changed` is `forcing` minus the selected term's own transform, so `forcing - changed` is just that term's transform.
- That is nonzero whether or not the term was ever added into `amplitudes`. The "responsive" check therefore only shows the term's transform is nonzero, not that it reached the assembled source.
- **Repair:** build the omitted source through the same assembly loop, skipping that term's address. Require that the full rebuild equals the baseline and that the omitted rebuild differs from it. Everything else stays the same.

**S3. The unit contraction isn't saved before its check** (l.499–505). `unit_checks` goes to `require` without being saved, which contradicts the directive's "guards only after persistence" rule. **Repair:** save it to JSON before the `require`.

## Answers to the review questions

1. **Transform rules: correct.**
   - I re-derived the base transform `πq/sinh(πq/2)` with its removable value 2, and the recurrence `I_{n+1}=(n I_{n-1}-iq I_n)/(n+2)` from `(sech²tanhⁿ)' = n sech²tanh^{n-1} − (n+2) sech²tanh^{n+1}`. Both are correct.
   - Rule 1 gives `2π c M(k) Ĥ(k)`. It is actually absolutely convergent: the inner `zs` integral decays and the outer one is controlled by B0.
   - Rule 2 gives `(2π)² c M(k,k0) ĥ(k−k0)`. The saved profile integrand already carries the factor L=10, and the worker correctly uses scale 1 for that profile.
   - Signs are correct: the source is the negated total, `direct == -total`, and the native terms are negated.
   - Assembly addresses match the producer: `row` / `frameColumn` in `S11c_d_two_asymptote_forcing_continue.py` l.509–516.
   - Bound-variable roles and whole-real limits are checked, and unexpected structures stop the run.
2. **Branch and denominator checks: sufficient.**
   - `sqrt(-100k²-4) → +i·r` matches SymPy's principal branch, which is also used numerically, and r ≥ 2 on the real axis.
   - If all coefficients of a polynomial in r share one strict sign, it can't vanish for r > 0. Its reciprocal and derivatives grow at most polynomially, and `dr/dk = bk/r` is bounded.
   - Any explicit k in a denominator, or a factor mixing momentum and r whose sign can't be decided, stops the run.
   - Rule 2 is certified in p for all real p, including p = k0. After collapse, the kout dependence is certified again through `certificate(amplitude)`.
   - Multiplying by B0 gives exponential decay, so the complete forcing transform is smooth and rapidly decreasing. I found no missing condition.
3. **PV/delta action: justified, and needs no correction.**
   - For a smooth, integrable source transform the momentum-space action is actually better defined than the separated-point kernel. It needs no diagonal extension and no evaluation of the kernel where source and field positions coincide.
   - The only requirements are ones the saved candidate already provides: simple real poles only at the two saved momenta (its residue identities and denominator certificates), no other real singularities, and bounded tails (its 50 tail checks).
   - The PV then exists for a smooth transform, and the delta terms correctly use the value at each pole. That includes the removable value at k = k0: here k0 = −√555/30 is itself one of the poles, and B0 at zero is correctly 2.
4. **Checks:**
   - The source joins, zero-lift, sign and no-Abel checks are adequate.
   - The dictionary quadrature is genuinely independent of the recurrence.
   - The equation probes and pole-null checks only confirm the saved inverse and the delta coefficients' null space. They say nothing about whether the forcing transform is right, but they are honestly labelled as contractions.
   - The omission control is circular (S2).
   - Failure persistence is otherwise sound: every `checked`, pole and mutation step saves first and the failure handler saves the whole stack. S3 is the only exception.
   - Runtime looks plausible for 840 s. There are no adjugate derivatives or solves, only small exact Poly/EX arithmetic and 80-digit substitutions into the saved inverse.

## Optional suggestions (not grounds for another review cycle)

- **Probe points:** the probe momenta 0 and ±1/L sit 7.9–8.9 widths out in the B0 tail, so the forcing is about 1e-4 of its peak there. Probes at k0 ± 1/L would be more informative.
- **Checkpoint status:** assert the status string of the pinned `outgoingCheckpoint`, not just its hash.
- **Probe order:** run `profile_probe` before the case construction.
- **End-to-end check:** add one mpmath quadrature of a local forcing entry against its assembled transform.
- **Summary field:** report `used_orders` per case instead of the running total across both cases.

**Smallest necessary repair:** the one-line exponent expansion in `remove_phase`, the omission control rebuilt through the assembly loop, and saving `unit_checks` before its `require`. No change to the method is needed.