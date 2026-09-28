# Endpoint-binding correction and saved-return forcing continuation

The user instructed **“Continue your next step”** after the
[completion report](S11c_d_two_asymptote_continuation_report.md) identified
the endpoint-binding error and the small repair followed by a saved-return
continuation. This authorizes one corrected ordinary bounded run. It does
not authorize an automatic retry, a new physics method, new external review
or inheritance of the old outgoing-constructor duration exception.

The failed worker remains unchanged. Its source SHA256 is
`cc68e3fc4ec16814b0433a9b52e2c58b42ce1fddb39a68d99628dd34c31758be`.
The [checkpoint](S11c_d_two_asymptote_continuation_checkpoint.json) preserves
all 536 complete returns, 2,366 files and the incomplete native-binding input.

The new `S11c_d_two_asymptote_forcing_continue.py` restores all 536 journal
returns by copying their input/value bytes and loading their saved returns;
it calls none of the prior operation functions. It also copies completed
non-journal scientific outputs, including two regular matrices, all eight
end-field bundles, residual summaries, source/unit/native bundles and the
already constructed auxiliary partition. End fields and source derivatives
are read from saved returns, not reconstructed. The original complete
Laurent/end-field loop is absent from this worker.

The sole mathematical-expression handling correction is inside `alpha_safe`:

```python
return visit(expression.xreplace(endpoint_values), {})
```

Previously alpha-renaming changed a bound variable inside a known endpoint
Limit before `bind` attempted an exact lookup. This correction substitutes
the same accepted endpoint map before renaming. It evaluates no Limit or
Integral, changes no endpoint value, and keeps the unknown-Limit guard.
The source map and actual native Limit operands are saved in a new journal
operation and concise JSON before native binding. The unfinished binding
call copies and restores its exact previously saved input; it requires
structural identity with the restored native source before proceeding.

All remaining forcing, decomposition, carrier reduction, actual native-term
omission, separate end-field control and domain formulas are unchanged.
Raw plane actions still save before reconstruction guards. Formal carrier
responsiveness is not a nonzero integrated response. The expected stop is
`END_LIFT_SAVED_FORCING_DOMAIN_UNRESOLVED` only if actual construction and
controls pass. Any new failure saves its inputs, full completed returns and
posthashes, then stops without retry.

The prior independent method review remains in force through its local
disposition: Claude's literal verdict was **NEEDS REVISION**, Grok's was
**CLEAR FOR THIS BOUNDED INGREDIENT**. This local implementation repair is
not a fresh independent CLEAR and does not trigger a new submission.

Use the unchanged shared guard around the existing normalization supervisor:
900s outer / 840s native, 2 GiB, zero swap, one CPU, nice15, 32 tasks and one
native thread. The new `forcing-continuation` result directory must have the
silent hook armed first for session `01a0e01b-ef84-7192-817f-584cda5d339b`.
No overlap, fallback, retry, new roots/modes/LU, producer reconstruction,
profile integration, physical-input change or current normalization.

Forcing pairing and the nonlocal plane-jet extension remain unverified.
Full response, Green/FORM, A11/A12 and a radiating witness remain unaccepted;
the fixed slice is evanescent. Preserve prior failures, accepted reference
and outgoing artifacts, incident history, Lean work, shared guard and the
protected builder suffix. Keep this repair at the practical toy-model level.
