# DIRECTIVE — S11b-A blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bA_interface_response_mathematica_audit.wl`

Run it with `math -script <path>`. Iterate until it completes with no errors and no unevaluated output.
Then stop and exit — ⛔ do not write a report, a summary document, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It exists to be an **independent** check on a SymPy audit of the same physics. If it agrees with that
audit because it copied from something, the check is worthless and the step is not verified.

⇒ **The read-bar list is §0b of the shared specification below**, and it is byte-identical to the one the
other engine receives. ⛔ Follow it exactly.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads** — the shape the 44
  Mathematica scripts already in this corpus use, and the source of their reliability.
- Prefix every output tag with `WL_`, so `S11BA_PROJECTION_FINITE` prints as `WL_S11BA_PROJECTION_FINITE`.
- Strip `ConditionalExpression[0, …]` wrappers when checking that something vanishes.
- Test for a pole with `1/expr == 0`, ⛔ not `expr == Infinity`.
- Keep total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically rather
  than raising a limit.

---

