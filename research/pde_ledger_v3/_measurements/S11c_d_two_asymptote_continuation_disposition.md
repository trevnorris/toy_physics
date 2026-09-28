# Two-asymptote continuation: local denominator-check repair

The user said **“Holy crap that's a lot. Well, let's get going”** after the
scope check and remaining-work rundown. This resumes the identified local
repair and unfinished end-lift/forcing stage, under the existing safeguards.
It does not authorize a new external review or a different physics method.

The [completed diagnostic](S11c_d_two_asymptote_entry_diagnostic_report.md)
established the actual failure: the selected first-order denominator
coefficient has `is_zero = None`, while its printed real part is a sum of
strictly negative terms. The original guard's cheap property test was
insufficient for that representation. The first production failure and the
completed diagnostic remain unchanged.

## Bounded correction

`S11c_d_two_asymptote_continue.py` retains the existing `is_zero is False`
test where it succeeds. Only when the property is undecided, it derives the
actual real and imaginary components, checks exact reconstruction of the
original coefficient, and examines their expanded terms' exact signs. A
component whose nonzero terms are all strictly positive or all strictly
negative proves the coefficient nonzero. The actual expressions, flags,
terms and reconstruction residuals are saved before enforcing the result.
If neither component establishes the claim, the stage stops with evidence.
No coefficient, expected sign or numeric threshold is hard-coded.

The quotient and residue formulas, branch choice, Laurent degree bound,
end-field construction, native action, formal carrier controls and limited
forcing-domain stop remain unchanged. Fourteen existing scientific helper
ASTs and the full construction after the regular-coefficient loop boundary
are checked against the original source. This is a local instrument repair,
not a general certification project or a new physical assumption.

## Exact saved-return boundary

The continuation restores all **61** completed parent journal returns and
all **12** completed diagnostic returns without executing their functions.
It also copies the completed non-journal unit/native/source/premise bundles
and source-join summaries byte-for-byte. The old construction prefix is
replaced by state restoration, so it does not rerun the native census,
grade/unit joins, end-pencil mappings or completed residual calculations.

The first new work checks that the saved diagnostic tuple matches the
current `block17-regular-0-0` operands, then certifies its already saved
denominator coefficient. Its branch/numerator/denominator series are reused.
Branch series and polynomial-series returns are also reused when their
complete exact arguments match; this changes no algebra. New intermediates
are journaled before later guards. A nested failure preserves both the
actual failing input and its enclosing incomplete operation stack.

The original failure, 279 parent artifacts, 59 diagnostic artifacts and
101 original input records are pinned. A new `continuation` directory is
used; no earlier root is restarted or overwritten.

## Authority, limits and stop

Literal parent build-review verdicts remain Claude **NEEDS REVISION**, Grok
**CLEAR FOR THIS BOUNDED INGREDIENT**, with the prior local disposition.
This repair has a local disposition, **not a fresh independent CLEAR**.
No reviewer rerun, peer sharing or new submission is needed for the resolved
property-test limitation, and none is launched. The earlier independent
method review is retained; it is not replaced by static checks.

The fresh gate pins the actual worker, complete input manifest, this
disposition, existing review adjudication and approval record. The shared
`scripts/s11c_guarded_run.py` wraps the existing normalization supervisor:
900 seconds outside / 840 native seconds, 2 GiB memory, zero swap, one CPU,
nice 15, 32 tasks and one native thread. No previous duration exception is
inherited. Global lock, overlap refusal, native limits and fail-closed
containment remain. The silent local hook must arm before science and target
session `01a0e01b-ef84-7192-817f-584cda5d339b`. No automatic retry or polling.

The retained stop is `END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED`, if all actual
construction/control checks pass. The nonlocal plane-jet extension and
whole-line forcing pairing remain explicitly unresolved; formal carrier
responsiveness is not a nonzero integrated response. No response, full
Green/FORM, A11/A12 or radiating witness is accepted by this stage alone.
There are no new inputs, roots, modes, LU solves, producer reconstructions,
profile integrations, current normalization or centre mechanics. Practical
toy-model scope, all earlier evidence, Lean work, shared guard and protected
builder suffix remain intact.
