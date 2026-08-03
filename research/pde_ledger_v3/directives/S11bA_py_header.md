# DIRECTIVE — S11b-A SymPy audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bA_interface_response_sympy_audit.py`

Run it. Iterate until it exits 0. Then stop and exit — ⛔ do not write a report or a summary document.

## ⛔⛔ THIS SCRIPT IS ONE OF TWO INDEPENDENT ENGINES.

A blind Mathematica audit of the same physics exists. **The disagreement between the two engines is the
test**, so an agreement reached by reading a shared source certifies nothing.

⇒ **The read-bar list is §0b of the shared specification below**, and it is byte-identical to the one the
other engine receives. ⛔ Follow it exactly.

⚠⚠ **Note that §0b bars `reduction/`, and that is deliberate for this sub-step.** Normally the SymPy engine
is the only one that sees the registry. ⭐ **Not here:** S11b-A derives no brane quantity and earns no
registry row — the assembly that would is deferred. Giving one engine access to recorded relations lets it
silently identify two symbols the other is treating as independent, and the two then "agree" on an
identification only one of them made. With the registry out of reach, both engines face **identical**
information, which is the condition their agreement is supposed to certify.

⚠ **A6's dimensions must be DERIVED from this specification's equations.** ⛔ Not read from any table.

## Script conventions

- Output one tag per line, `TAG: value`, with no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.
- ⛔ Do not write a test file, a fixture, or a gate — this sub-step produces one script and its printed
  tags.

---

