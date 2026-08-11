# Build directive — CUT the witness exchange, and simplify

⭐ **Target:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md`.
⛔ **Edit ONLY that file.** ⛔ Do not commit. ⛔ Do not run either engine.

⚠⚠ **USER DECISION: `§Q8c`, the witness exchange, is CUT.** ⭐ Everything that exists **only** to serve it
goes with it. ⛔ This is a removal round — ⛔ nothing is replaced, and ⛔ no new mechanism is introduced.

⚠ **Why**, so you can judge scope correctly: the exchange required the two engines to share a coordinate
vocabulary, and the blind build exists so that they do not. ⭐ It was the single largest source of defects
in this repair.

---

## ⭐⭐ C1 · REMOVE the witness exchange and everything that serves only it

⭐ Remove:
- ⭐ **`§Q8c`** in full — the two-pass protocol and every `WITNESS<w>_*` object.
- ⭐ **`§4`'s permitted-input subsection** — it exists only so a witness point could be received.
  ⚠ ⛔ **Do NOT touch `§4`'s quoted block itself**, which is shared verbatim with the other steps.
- ⭐ **`§8`'s `WITNESS<w>` tag-grammar lines**, its `<w>` index rule, and `WITNESS_ORDERING` wherever `§8`
  lists what is emitted once per package-and-dimension.
- ⭐ **Corollary 4's witness exception** — ⭐ it returns to permitting exactly package, dimension and
  quantity.
- ⭐ The **transport-scope** clause and its `POINT_TRANSPORT_ELIGIBLE_LOCI` / `POINT_TRANSPORT_EXCLUDED_LOCI`
  tags, and the `§9` row that records their boundary.
- ⭐ The **per-object reading of `WITNESS<w>_POINT_COVERAGE`** — that token disappears with `§Q8c`.

## ⭐ C2 · REMOVE the shared payload-key vocabulary, and RESTORE `§Q6r`'s own map

⚠ `§8`'s vocabulary subsection was declared **only** so a transported point could be keyed by it.
⭐ Remove it, and ⭐ **restore `§Q6r`'s own closed coefficient map verbatim**, exactly as it read before —
⛔ so `§Q6r` again carries the names it uses and depends on nothing outside itself.
⛔ Do not leave `§Q6r` referencing a section that no longer declares anything.

## ⭐ C3 · KEEP these — ⛔ they are not transport machinery

⛔ **Do not remove, weaken or re-open:**
- ⭐ **The definition of a point** — *an assignment to exactly the solve variables its locus names* — at
  both `§5`'s `_REAL_WITNESS` and `§Q8b`'s `STRATUM<s>_POINT`. ⚠ It is a **definition**, ⛔ not transport.
- ⭐ **All of `§Q8b`'s component-scoped machinery**: the component-scoped `Q3`/`Q4` objects, the count
  record and its status/certificate/change-locus families, the coverage token, the `POINT_EVIDENCE_` infix.
- ⭐ **The component comparison**: the difference reduced modulo **each engine's own** defining equations,
  both carried as operands.
- ⭐ `§5`'s locus protocol in full, its typed branch statuses and undecided handling, and `_CANONICAL_LOCUS`
  as it stands.

## ⭐⭐ C4 · SAY PLAINLY, IN `§9`, WHAT IS NO LONGER ESTABLISHED

⚠ Without the exchange, where one engine admits a component and the other cannot decide, the two emit
typed, differing objects and ⭐ **that difference is a finding for the orchestrator to adjudicate** — ⛔ but
**no computation resolves it**, and the build establishes nothing about the counterpart's component from
this engine.
⭐ **Record that in `§9`**, which is where this file states what it does not establish.
⛔ Do not describe the removed mechanism, ⛔ do not say it was considered, ⛔ and do not state any measured
case.

---

## ⛔ Constraints

- ⛔⛔ **THE SPEC SAYS WHAT TO COMPUTE. ⛔ NEVER what anything equals or was measured.**
- ⛔⛔ **`§4`'s quoted block must stay byte-identical** to `directives/S10_SHARED_PHYSICS.md:111-113`.
  ⭐ Verify and report.
- ⛔ Do not renumber sections. ⛔ Do not reflow untouched prose.
- ⚠⚠ **A REMOVAL FALSIFIES DOWNSTREAM SENTENCES THE SAME WAY AN EDIT DOES.** ⭐⭐ **Grep the whole file for
  every reference to anything you deleted** — witness, two-pass, transport, the vocabulary, the coverage
  token — ⛔ and for sentences quoting a **consequence** of their existence. ⚠ **Seven of eight rounds left
  exactly one stale sentence; a removal round is the easiest place to leave several.**

## Deliverables

1. The edited file.
2. ⭐ Per `C1`–`C4`: what you removed or added, with line numbers.
3. ⭐ The `§4` byte-identity result.
4. ⭐ Your stale-reference sweep: every place that mentioned a removed object, and what you did with it.
5. ⛔ Anything you could not remove cleanly, or anything in `C3` that turned out to depend on `C1` ⇒
   **STOP AND REPORT.** ⭐ Reporting is success.
