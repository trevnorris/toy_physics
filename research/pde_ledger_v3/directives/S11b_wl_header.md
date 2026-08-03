# DIRECTIVE — S11b blind Mathematica audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11b_brane_bulk_interface_mathematica_audit.wl`

Run it with `math -script <path>`. Iterate until it runs to completion with no errors and no unevaluated
output. Then stop and exit — ⛔ do not write a report, a summary document, or a second script.

## ⛔⛔ THIS SCRIPT IS BLIND. THAT IS ITS ENTIRE PURPOSE.

It exists to be an **independent** check on a SymPy audit that does not yet exist. If it agrees with that
audit because it copied from something, the check is worthless and the step is not verified.

⛔ **DO NOT READ, open, grep, `cat`, `git show`, or otherwise inspect:**
- `/var/projects/toy_physics/research/pde_ledger_v3/reduction/` — any file, especially `quantities.yaml`
  and `relations.yaml`
- any `.py` file anywhere under `/var/projects/toy_physics/research/pde_ledger_v3/`
- `/var/projects/toy_physics/research/pde_audit/` — any file
- any file whose name contains `PREREGISTERED`

⭐ **You do not need any of them.** Every input is in the specification below. If you believe an input is
missing, emit the relevant tag as `NOT_ESTABLISHED` and name what is missing. ⛔ Do not go looking.

## Script conventions

- **Standalone, print-only, no arguments, no exports, no external file reads.** This matches the 44
  Mathematica scripts already in this corpus, and their reliability comes from exactly that shape.
- Prefix every output tag with `WL_` — so `S11B_PROJECTION_IDENTITY_FINITE` is printed as
  `WL_S11B_PROJECTION_IDENTITY_FINITE`.
- Strip `ConditionalExpression[0, …]` wrappers when checking that something vanishes.
- To test whether a quantity has a pole, test `1/expr == 0` rather than `expr == Infinity`.
- Keep the total runtime under **10 minutes**. If a computation will not finish, reduce it symbolically
  rather than raising a limit.

---

