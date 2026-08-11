# Decision list v2 — the export chain

**Orchestrator-written, 2026-08-10, folding the ten findings that blocked v1.**
⚠ v1 (`36589024`) is a **record**; ⛔ do not read its `E1`–`E7` as decisions.

⛔ **This list is not itself a build directive.** It settles the chain questions the **S11 PY decision
list** must answer; ⭐ that list gets the two legs, before any builder launches (rule 7).

⭐⭐ **What v1 got wrong, in one line:** it made cross-step collisions impossible by construction — ⛔ and
the collision **is** the measurement. Two steps deriving one object must be able to **meet**.

---

## The decisions

**F1 · Storage keys stay FLAT. `D5` is unchanged; ⛔ there is no producer prefix.**
⭐ A step writes the object's name. A later step re-deriving that object writes **the same key**.

**F2 · Before writing a key that exists in the imported `LEDGER`, the writer compares the OBJECT.**
⭐ This is `DEFECT_REGISTER.md:675`, which already prescribed it.
- ⭐ **Same object** ⇒ a **re-derivation**: emit both operands and the residual, then guard.
- ⛔ **Different object** ⇒ a **finding that fails loudly**, naming both. ⛔ Never a silent overwrite.
⚠ Deciding "same object" is the load-bearing part and belongs to the S11 PY list.

**F3 · A re-derived row carries its own evidence, in the row.**
⛔ `corroborated_steps` alone is an agreement claim with no operands in the file that carries it, and
`S11:527-529` **specifies `Q6r` to propagate it** ⇒ a consumer forwards a claim it cannot check.
⭐ A consumer reading **only** the merged export must be able to recompute the residual.

**F4 · S10's export is REGENERATED under `F3`.** ⚠ ⭐ Every value must be **byte-identical** apart from the
fields `F3` adds — ⛔ if any value moves, it was not a regeneration.

⛔⛔ **F4's "⛔ NOT re-keyed, tag names do not move" is SUPERSEDED by `F8`, 2026-08-11.** ⚠ It rested on a
measurement taken before S11's namespace was censused; that census found S10 keys claimed by two different
objects. ⛔ Do not cite F4 as authority against a rename.

**F5 · `C19` is a REAL deviation, ⛔ not a traversal gap.** `S10_SHARED_PHYSICS.md:197` orders
*"Emit `M_A`, `M_B`, and `M_A − M_B`"*; both engines emit `Q2_MATRIX_A/B`.
⭐ S10's step record must **disclose** it. ⛔ The rename itself is its own gated workstream — ⚠ injectivity
across the worklist first — and ⭐ **S11's build does not depend on it**, because S11's spec fixes its own
quantity names.

**F6 · An export is published only if every declared `MAIN` cell completed** — ⭐ or it carries
machine-readable per-cell completeness that every consumer checks. ⛔ One or the other, chosen before the
writer is built (`S11:895`).

**F7 · A `Q6r` lookup that cannot resolve is reported as unresolved, generically.**
⛔ No placeholder, ⛔ no exception, ⛔ no membership stated anywhere a builder can read.

---

## ⛔⛔ F8 — **BLOCKED RECORD, 2026-08-11. ⛔ NOT A DECISION. ⛔ DO NOT BUILD FROM IT.**

⚠⚠ **Two legs blocked it**, and `F1` and `F4`'s supersession above go with it — ⛔ neither survives as
written. ⭐ **What replaced all of it: `directives/S11_sympy_export_decisions_codex.md`**, which supersedes
`F1` outright — ⭐ **a key is an opaque row locator, ⛔ not an object name** — and moves identity into a
typed per-row claim. ⇒ ⛔ **There is no rename**, so ⛔ every control below is moot.

⭐ **Why it is kept:** it records the third and last attempt to make a **key string** carry object identity.
⛔ The measurement that ended the approach: an **already-wrong** frozen migration map passed **all six**
controls below, and ⛔ nothing bound a destination name to the object it named. ⇒ ⭐ `_legs/round3_*`.

---

### (blocked) A key names the object — and for an object derived from an action, that includes which action

⚠⚠ **The defect, and it is S10's name, ⛔ not S11's.** `root_ordering_d3` names a **slot in a procedure**
that every step has, ⛔ not an object. S10's value for it is the spectrum of a brane that resists shear
only; S11's is the spectrum of a brane that resists shear **and compression** — the two differ in exactly
the term the step exists to test. ⇒ ⛔ One key, two different objects.

⭐⭐ **F1 is unamended: keys stay flat, and a collision is still the measurement.** ⛔ F8 adds no machinery
to prevent collisions. ⭐ It makes the name say enough that **only real meetings collide.**

**F8a · The partition is MEASURED, ⛔ never declared.** ⭐ Mutate the action; ⛔ if the object does not move,
it is **universal** — one object for every step — and keeps a **bare** name, and its collision across steps
is the cross-step check working. ⭐ If it moves, it is **action-dependent** and its key must say which
action. ⚠ This is the form ablation every review leg already runs (rule 14), and the shared spec already
draws this line in physics terms at `S11_SHARED_PHYSICS.md:604-608`.

**F8b · The qualifier identifies the ACTION, ⛔ not the step, ⛔ not the package label, ⛔ not the question
the step was asking.** ⚠ A step's subject (`brane_mode_…`, `stray_longitudinal_…`) fails `D11`: it says why
an object was requested, ⛔ not what it is — and two steps may analyse the same material while one step
analyses eight.

⚠⚠ **A support-level slug is NOT an identity.** Measured on §7: `MAIN`, `XCOEF_BSCALE` and `XCOEF_BSIGN`
all carry curl-plus-divergence support and differ only in a coefficient or a sign; `MAIN` and `XKIN_ANISO`
share `W` and differ in `T`. ⇒ ⛔ a slug naming which primitives appear **aliases genuinely different
materials**.

**F8c · Every export carries the CANONICAL DESCRIPTOR of the action it was derived from, as its own row,
computed from the LIVE action object.** ⛔ Not from the package label, ⛔ not from the step name, and ⛔ never
from any computed result — ⚠ an identity derived from the answer lets a wrong answer pick a different key
and escape comparison.
⭐ The key carries a short readable slug; ⭐ the descriptor is what makes the slug auditable:
- ⛔ two steps whose descriptors are **equal** while their slugs differ ⇒ a finding — a real meeting was
  missed;
- ⛔ two steps whose slugs are **equal** while their descriptors differ ⇒ a finding — different materials
  collided.

⭐⭐ **The ledger key namespace is written by ONE engine.** The Wolfram engine imports nothing and writes no
ledger, so ⛔ **no cross-CAS canonical serialization is required, and none may be built** — that is the
shared-vocabulary problem that `C18`'s `T3` was cut for. ⭐ The rule must be stable **across steps**, which
is a far weaker requirement.

**F8d · Dependency scope is the FULL action pair `(T, W)`, for every action-dependent object.**
⚠ ⭐ **Known cost, recorded rather than engineered around:** an object depending only on `W` will not meet
its counterpart under a different `T`. ⛔ Do not build a per-object dependency scope to recover it — ⭐ revisit
only if a step exports more than one package.

### ⛔⛔ The rename is a RE-POINTING, and re-pointing is the failure where EVERY repo check passes

⚠ Measured precedent: re-pointing one table line made the light cone an `ω²` and **all** repo checks passed.
⇒ ⭐ The controls are the deliverable, ⛔ not the rename:

1. ⭐ **A migration map frozen BEFORE any edit**, pinned against the committed artifact — ⛔ **total** over
   every S10-owned row and ⛔ **injective**.
2. ⭐ **Every row preserved through the map** — payload and classification fields object-identical, except
   the evidence fields `F3` authorises.
3. ⭐ **Every reference rewritten through the same map** — a `dimension_key` must resolve to the same record
   it resolved to before.
4. ⭐ **The destination name bound to non-payload evidence**: recompute the descriptor from the live action
   and check it against what the key claims.
5. ⛔⛔ **A DELIBERATE DERANGEMENT MUST FAIL THE GATE.** Swap two destination keys — ⭐ including rows whose
   payloads are equal, zero, or symmetry-invariant — and require a failure. ⚠ Payload equality alone
   **cannot** police a rename; that is why control 4 exists.
6. ⭐ **The meeting rules exercised as fixtures**: same action + same quantity across steps ⇒ same key and an
   `F2` residual; different action + same quantity ⇒ different keys; an under-specified collision ⇒ a loud
   `F2` failure.

---

## ⛔ What this list does not decide

- ⭐ The **"same object" predicate** of `F2`, and the row shape of `F3` ⇒ the S11 PY decision list.
- The `C19` rename worklist ⇒ its own gate, per `F5`.
- `C17`, `C18`, S10's requirements registers.
