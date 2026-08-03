# DIRECTIVE — S11b-B blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bB_interface_assembly_mathematica_audit.wl`

Run with `math -script <path>`. Iterate until it completes with no errors and no unevaluated output. Then
stop and exit — ⛔ do not write a report, a summary, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It is an **independent** check on a SymPy audit that does not yet exist. If it agrees with that audit
because it copied from something, the check is worthless and the step is not verified.
⇒ **The read-bar list is §0b of the shared specification below**, byte-identical to the other engine's.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads.**
- Prefix every output tag with `WL_`.
- Strip `ConditionalExpression[0, …]` when checking that something vanishes; test poles with
  `1/expr == 0`.
- Keep total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically rather
  than raising a limit.
- ⚠ **On a non-uniform background do not silently assume plane waves are eigenmodes.** Where you must
  restrict to a tractable form, ⛔ say so in the emitted value rather than in a comment.

---

