# Decision list v2 — repairing `C17` and `C18` in S11's shared spec

**Orchestrator-written, 2026-08-10, folding both legs.** ⚠ v1 (`767828f4`) was **BLOCKED by both** — ⛔ it is
a record. ⭐ The **repaired spec** gets its own two legs; ⛔ this list is not re-legged (rule 7).

⭐⭐ **ORDERING, settled by both legs independently, each by opening the outputs:** ⭐ **repair the spec
BEFORE either engine is built.** At `XFORM_EXTRA, D = 2` both engines emit the **same defining equation**
and diverge only in CAS-specific branch and admissibility handling ⇒ ⭐ the measured divergence is a finding
about the **PROSE**, ⛔ not about the physics. ⚠ Building first re-pays for a known non-physics split.

---

## The decisions

**T1 · A stratum's `Q3`/`Q4` are computed ON THE COMPONENT, in its free parameters — and every count-valued
tag carries its STATUS on that component:** constant (with the certificate), **varies** (with the sub-locus
where it changes), or undecided.
⛔⛔ **A bare generic integer is not admissible.** ⚠ Both legs measured that `S1` alone **relocated** `C17`:
the component-scoped count is the generic one, and a sub-locus where the count moves gets **no tag at all**.

**T2 · A point-local evaluation is EVIDENCE, ⛔ never the component's answer.**
⭐ An engine that cannot build `T1`'s object emits its point-scoped result **and an explicit
incomplete-coverage token**; ⛔ no component-level count may be claimed from it.
⚠ Without this, both engines take the labelled fallback and `C17` survives **with a label on it**.

**T3 · ⛔⛔ WITHDRAWN, 2026-08-11 (user decision). ⛔ Do not build this.**
⭐⭐ **Why, and it is a finding rather than a retreat:** the exchange requires the two engines to share a
**coordinate vocabulary**, ⛔ and the blind build exists precisely so they do not. ⚠ Nine rounds
established that every route to that correspondence violates a rule of the same file ⇒
`DEFECT_REGISTER.md#C18` carries the full reason. ⭐ `§9` now records that no computation resolves the gap.

~~**T3 · ⭐⭐ WITNESS EXCHANGE — where one engine has a point and the other has none.**~~
⭐ Each engine evaluates the **other's** emitted point against **its own** `_EQUATIONS` and **its own** `M`,
and emits the residual.
⇒ ⭐⭐ **The divergence becomes a computation instead of an incomparability.** ⚠ This is the only decision
that closes the measured `XFORM_EXTRA, D = 2` case; ⛔ `T4`+`T5` alone leave it as silence.

**T4 · ⛔ A BRANCH IS NEVER SILENTLY DROPPED.** ⭐ An excluded branch is emitted with the test's operands; an
undecidable one gets an explicit **undecided** token — ⛔ never a bare `False`, which makes an omission
indistinguishable from a refutation.

**T5 · Undecided is an explicit COVERAGE FINDING that BOUNDS what may be claimed.**
⛔ Not agreement, ⛔ not disagreement, ⛔ **and not silence.** ⚠ Both legs: *"not compared"* is not an
acceptable terminal state for a region where the **mode count moves**.

**T6 · The canonical locus object is required where the system is POLYNOMIAL, and explicitly
`NOT_APPLICABLE` otherwise** — with the real status typed: proved empty · proved non-empty **with an exact
witness** · undecided.
⚠⚠ **Measured by both legs, independently:** the equations this step actually produces **carry radicals**.
SymPy rejects a Gröbner basis on them; Wolfram returns an object still containing the radical. ⛔ And a
**complex** basis is not the **real** locus. ⇒ ⛔ v1's `S3` was unbuildable.

**T7 · ⛔ THE COMPARATOR CONTRACT IS FROZEN BEFORE IT SEES EITHER OUTPUT, and it must:**
- ⛔⛔ **reject a native boolean as a residual operand.** ⚠ **Measured in the live comparator:** the parsers
  turn the token `False` into a **Python bool**, `residual_is_zero` is `value == 0`, and `False` vs `0`
  scores as **AGREEMENT** (so does `True` vs `1`). ⭐ Lowercase `false` correctly mismatches.
  ⚠ **Not realized to date** — S10's join carries **0** boolean payloads in 562 keys. ⛔ **But S11's locus
  protocol requires three booleans per locus in both engines**, and S11's old engines already show **400**
  boolean pairs in a 3042-key join ⇒ ⭐ this is on the critical path.
- ⭐ treat undecided per `T5`.

---

## ⛔ What this list does not decide

- ⛔ Any **admissibility algorithm** — ⚠ `S11:232-237` measures that the two CASes have unequal real-domain
  capability; ⛔ requiring one symmetrically asks an engine for what it cannot do.
- ⛔ **Recursive stratification to canonicalise components** — ⚠ cylindrical-decomposition-grade, ⛔ no two
  independently built engines implement it identically.
- ⭐ `C19`, `C20`, the export chain ⇒ `S11_export_chain_decisions_v2.md`.
