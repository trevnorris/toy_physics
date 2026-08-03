# DIRECTIVE — S11b-A blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bA_interface_response_mathematica_audit.wl`

Run it with `math -script <path>`. Iterate until it completes with no errors and no unevaluated output.
Then stop and exit — ⛔ do not write a report, a summary document, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It exists to be an **independent** check on a SymPy audit that does not yet exist. If it agrees with that
audit because it copied from something, the check is worthless and the step is not verified.

⛔ **DO NOT READ, open, grep, `cat`, `git show`, or otherwise inspect:**
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/` — any file
- any `.py` under `/var/projects/toy_physics/research/pde_ledger_v3/`
- `/var/projects/toy_physics/research/pde_audit/` — any file. ⚠ It contains prior work on adjacent
  physics whose current validity is **not established**; the other engine is barred from it too, and an
  asymmetry here would make one engine a transcriber and the other a deriver.
- any file whose name contains `PREREGISTERED`

⭐ **You do not need any of them.** Every input is in the specification below. If you believe an input is
missing, emit that tag as `NOT_ESTABLISHED` and name what is missing. ⛔ Do not go looking.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads** — the shape the 44
  Mathematica scripts already in this corpus use, and the source of their reliability.
- Prefix every output tag with `WL_`, so `S11BA_PROJECTION_FINITE` prints as `WL_S11BA_PROJECTION_FINITE`.
- Strip `ConditionalExpression[0, …]` wrappers when checking that something vanishes.
- Test for a pole with `1/expr == 0`, ⛔ not `expr == Infinity`.
- Keep total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically rather
  than raising a limit.

---

