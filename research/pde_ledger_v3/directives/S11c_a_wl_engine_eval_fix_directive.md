# S11c-a Wolfram engine — evaluation-error fix directive

## Authority and boundary

Edit `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl` **in place**. Its
product is the flushed stdout tag stream; that is its only write. `CLAUDE.md` binds. The physics authority
remains `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` (`2926c71c`); the engine's build
directive is `research/pde_ledger_v3/directives/S11c_a_wl_build_directive.md` (`acbdd236`). This is a
**mechanical evaluation-error fix**, ⛔ not a physics or structural change: preserve every emitted
`WL_S11CA_<QUANTITY>` tag, the blindness (import nothing), the tag grammar, and the derivations. Add no
expected value (`CLAUDE.md` rule 5).

## The measured defect

Running the finished engine (`wolframscript -file …`) exits 0 but prints these evaluation messages to
stderr (the only other stderr line, `WolframScript.conf`, is a benign wolframscript config warning to
ignore):

```
Part::pkspec1: The expression index cannot be used as a part specification.
AssociateTo::argrx: AssociateTo called with 3 arguments; 2 arguments are expected.
```

Both are suppressed by Mathematica after a few occurrences, so the counts in stderr understate them.
⛔ These are real bugs: a Mathematica message means a subexpression did **not** evaluate as written, so a
payload can be silently wrong or a control-operand map silently incomplete.

## What to fix — root cause only, ⛔ never suppression

⛔⛔ **Do not silence either message with `Off[…]`, `Quiet[…]`, `Check[…]`, or by deleting the offending
computation.** Fix the construct so it evaluates correctly and the message no longer arises.

### Bug 1 — `AssociateTo::argrx` (a rule-list passed as separate arguments)

Several `AssociateTo` calls pass **two rules as two separate arguments** instead of one list of rules.
`AssociateTo[assoc, key -> val]` takes exactly two arguments; two associations at once must be
`AssociateTo[assoc, {k1 -> v1, k2 -> v2}]` (a list), or two separate `AssociateTo` calls. Confirmed sites
include the form-control virtual/evolution pair and the uniform-limit virtual/evolution pair, e.g.:

```
AssociateTo[formBaseOperand, virtualKey -> objectRecord[virtualBase, dimensionZero],
  evolutionKey -> objectRecord[evolutionBase, dimensionEvolution]];      (* 3 args — WRONG *)
AssociateTo[uniformS11caOperand, virtualKey -> objectRecord[…], evolutionKey -> objectRecord[…]];  (* WRONG *)
```

Find **every** such multi-rule `AssociateTo` in the file and rewrite each so both rules are inserted (list
form, or two calls). The visible effect of the bug is that `WL_S11CA_CONTROL_FORM_*` and
`WL_S11CA_UNIFORM_LIMIT_*` are **missing** the virtual-work (`VIRTUAL_WORK…`) and evolution
(`EVOLUTION_MASS_BALANCE…`) keyed entries for the affected cases; after the fix those keyed entries must be
present in the emitted operand maps.

### Bug 2 — `Part::pkspec1` (a part index that is not a valid part spec)

Somewhere an `expr[[index]]` (or similar) evaluates with `index` — or the head being indexed — not a valid
part specification (e.g. an unbound iterator symbol, or a scalar/association where a list was assumed). It
appears **first** in the run, so it may sit in an early/primary derivation, not only in the controls —
localise it (Mathematica's `Stack[]`/`Trace`/binary-searching the tasks) and fix the actual construct so
the part access is well-defined. ⛔ Do not guard it away.

## Confirm

Re-run the finished engine and verify: (1) stderr contains **only** the benign `WolframScript.conf` line —
no `Part::pkspec1`, no `AssociateTo::argrx`, no other message; (2) all the same `WL_S11CA_<QUANTITY>` tags
still emit (the primary T-0…T-i set, the §2/§3 supplied tags, and the REP_INVARIANCE / CONTROL_INDEPENDENCE
/ CONTROL_FORM / UNIFORM_LIMIT / HOMOGENEITY / DIMENSIONS / LOCAL_TAG_NAMES controls); (3) the CONTROL_FORM
and UNIFORM_LIMIT operand maps now carry their virtual-work and evolution keyed entries. Report the before/
after stderr and the fixed sites. Demonstration runs go under scratch, ⛔ never `mathematica/out/`; one
kernel at a time (two-seat licence); kill a demo at 600 s no-output or 6 GB RSS. Run with
`--sandbox danger-full-access`.

## Conflicts

If a message traces to a genuine physics ambiguity rather than a syntax slip, ⛔ do not invent a value or
suppress it — report it under the spec's §8 builder report with what is undefined and why.
