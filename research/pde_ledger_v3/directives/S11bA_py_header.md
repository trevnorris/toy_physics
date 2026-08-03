# DIRECTIVE — S11b-A SymPy audit

**Deliverable (absolute path):**
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11bA_interface_response_sympy_audit.py`

Run it. Iterate until it exits 0. Then stop and exit — ⛔ do not write a report or a summary document.

## ⛔ WHAT YOU MUST NOT READ

- `/var/projects/toy_physics/research/pde_ledger_v3/mathematica/S11bA_*` — a blind Mathematica audit of
  this same physics exists. It has been moved out of the tree; ⛔ do not retrieve it from git either
  (`git show`, `git cat-file`, `git checkout`). **The disagreement between the two engines is the test.**
- `/var/projects/toy_physics/research/pde_audit/` — any file. ⚠ Prior work on adjacent physics whose
  current validity is **not established**. The other engine is barred from it too; an asymmetry here would
  make one engine a transcriber and the other a deriver, and their agreement would then mean nothing.
- any file whose name contains `PREREGISTERED`.

⛔ `git status` and `git diff` are fine. `git show` and `git cat-file` are not.

## ⛔⛔ NO REGISTRY WORK IN THIS SUB-STEP

⚠ Normally the SymPy engine is the only one that sees `reduction/`, and that is deliberate. ⭐ **Not here.**
This sub-step derives no brane quantity and earns no registry row; the assembly that would is deferred to
S11b-B.

⇒ ⛔ **Do not read, import, or modify `reduction/`** — not `quantities.yaml`, not `relations.yaml`, not
`registry_read.py`, not the gates. ⚠ **The reason is a measured one:** giving one engine access to
recorded relations lets it silently identify two symbols the other engine is treating as independent, and
the two then "agree" on an identification only one of them made. With the registry out of reach, both
engines face **identical** information, which is the condition their agreement is supposed to certify.

⚠ **A6's dimensions must be DERIVED from this specification's equations.** ⛔ Do not read them from any
table, registry, or gate.

## Script conventions

- Output one tag per line, `TAG: value`, with no `WL_` prefix.
- Keep total runtime under **10 minutes**; runners get `timeout 600`.
- ⛔ Do not write a test file, a fixture, or a gate — this sub-step produces one script and its printed
  tags.

---

