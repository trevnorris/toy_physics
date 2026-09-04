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
⚠⚠ **⛔ AMENDED BY `F9`:** a producer prefix **is** written when the override would be wrong. ⛔ Read `F9`.

**F2 · Before writing a key that exists in the imported `LEDGER`, the writer compares the OBJECT.**
⭐ This is `DEFECT_REGISTER.md:675`, which already prescribed it.
- ⭐ **Same object** ⇒ a **re-derivation**: emit both operands and the residual, then guard.
- ⛔ **Different object** ⇒ a **finding that fails loudly**, naming both. ⛔ Never a silent overwrite.
⚠ Deciding "same object" is the load-bearing part and belongs to the S11 PY list.
⚠⚠ **⛔ THE SECOND BRANCH IS SUPERSEDED BY `F9c`** — a differing object is the **ordinary** case and ⛔ must
not stop the run. ⭐ The first branch survives as `F9b`. ⭐ `F9` settles the predicate. ⛔ Read `F9`.

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
written. ⇒ ⛔ **There is no rename**, so ⛔ every control below is moot.

⛔⛔ **STALE POINTER CORRECTED 2026-08-12.** ⚠ This block used to say `F8` was replaced by
`directives/S11_sympy_export_decisions_codex.md`. ⛔ **That file is design 4 of the four blocked designs**
(`DEFECT_REGISTER.md#c20`) — ⚠ it moved identity into an **authored per-row claim** and nothing recomputes
it. ⭐ **What actually replaced all of it is `F9`**, below.

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

## ⭐⭐ F9 — WHEN A STEP'S EXPORT WRITES A KEY THE IMPORTED LEDGER ALREADY HAS

**Orchestrator-written 2026-08-11, on the user's decision.** ⭐ It closes `F2`'s open part — the
*"same object"* predicate — and amends `F1`: ⛔ a producer prefix is no longer universally forbidden; it is
what a step writes **when the override would be wrong**.

⭐⭐ **The user's instruction, and it is the whole of F9:** *"Only add `s11_` if it would override something
it shouldn't. The point was to override the times when it should."*

### ⭐⭐ THE RULE

**F9a · The key is ABSENT from the imported `LEDGER`** ⇒ ⭐ write the bare key.

**F9b · PRESENT, and the two objects are PROVED EQUAL** ⇒ a **re-derivation**. ⭐ Write the bare key.
⚠⚠ **This is the override-when-it-should case.**

**F9c · PRESENT, and equality is NOT PROVED** — ⭐ proved different, ⛔ **or the test returned no decision**
⇒ ⭐ write the **prefixed** key, and ⛔ leave the imported row exactly as it stands. ⚠ The report
distinguishes proved-different from undecided: ⭐ `§5`'s locus protocol already fixes that reading —
⛔ **do not re-word it.**

⭐⭐ **WHICH ASSUMPTIONS DECIDE IT: only the ones the two live objects themselves carry.** ⛔ Not a premise
this step asserts on top of them.
⚠⚠ **Measured, on the real colliding row:** the same pair is **UNDECIDED** under the objects' own
assumptions and **PROVED EQUAL** under this step's joint premise set — ⛔ and the second reading
**overwrites a predecessor whose step never asserted that premise.** ⇒ ⭐ the objects' own assumptions are
the conservative reading and the only one both steps hold.

⛔⛔ **F9c SUPERSEDES `F2`'s second branch.** ⚠ `F2` calls a differing object *"a finding that fails
loudly"*; ⭐ under `F9` it is the **ordinary case** — two steps analysing two different materials — so ⛔ it
must **not** stop the run. ⭐ `F2`'s first branch is unchanged and is now `F9b`.

⭐⭐ **THE COMPARISON MUST BE TOTAL OVER THE ROW SHAPES THE IMPORT ACTUALLY CARRIES**, ⭐ elementwise into
containers, ⭐ three-valued. ⚠⚠ **Measured: most of the imported namespace is NOT an expression** — tuples,
booleans, strings, relations and predicates ⇒ ⛔ a routine built around a subtraction residual raises on
them, falls back to structural equality, ⛔ and then `F9b` never fires — which deletes the
override-when-it-should case entirely.

⚠ ⭐ **`F9b` CARRIES NO RESIDUAL.** ⛔ The branch is entered *because* equality was proved, so a residual is
structurally zero and can never fire ⇒ ⭐ **the informative object is the OPERAND PAIR**, which a consumer
recomputes for itself. ⛔ Delete the residual; ⛔ do not repair it.

⭐⭐ **That is the whole rule.** ⛔ There is no fourth relation, ⛔ no substitution to derive, ⛔ no action to
locate.

⚠ ⭐ **F9 IS PY-ONLY.** The WL engine imports nothing and writes no ledger ⇒ ⛔ there is no joint property
here, and ⛔ this must not enter `S11_SHARED_PHYSICS.md`.

### ⛔⛔ (blocked record) ROUND 1 PROPOSED A THIRD OUTCOME — ⛔ DO NOT RE-PROPOSE IT

⚠ **A SUPERSESSION branch:** *not equal, but this step's object carries the imported one as a special case
⇒ write the bare key anyway.* ⭐ **Two legs killed it, and the measurement is one line:**
⛔⛔ **`B_comp > 0` IS A DECLARED PREMISE** (`S11_SHARED_PHYSICS.md:69`, `:119`, `:1027`, `:1126`)
⇒ ⭐⭐ **S11's model does NOT contain S10's as a special case OF ITSELF.** ⛔ The specialisation lies outside
the premise set.
⇒ ⭐ Whether a later result supersedes an earlier one is a claim for the **step record**, ⛔ never a ledger
key. ⚠ Evidence: `/tmp/f9_leg_{codex,grok}/`, ⛔ outside the tree.

### ⛔ What F9 does NOT decide

- ⛔ Whether S10 is regenerated ⇒ `F4`, still open.
- ⛔ The **model-level vs step-level** split — ⭐ the real distinction. ⭐⭐ **MADE via the bind test in
  `export_ledger_bind_closure_design.md` (`c04e071f`)** on S11c-c1's evidence (F10 below made a first, category-
  based cut and is now superseded by that design): the LEDGER carries the **bind-closure**, step-level
  diagnostics are emit-only.
- ⛔ The **tag stream** — `§8`'s grammar is untouched.
- ⚠ Whether the prefix generalises from `s11_` to `<step>_`.
- ⛔⛔ **WHAT "BARE" MEANS TO A LATER READER.** ⚠ Measured: ⛔ no presently instantiated broken consumer, ⭐ and
  a latent one — after `F9c` the bare key still answers with the **predecessor**. ⇒ ⭐ **a routed-key
  contract is owed BEFORE any consumer of an `F9c` object is built.**

---

## ⭐ F10 · THE LEDGER CARRIES THE MODEL-LEVEL REGISTER; STEP-LEVEL DIAGNOSTICS ARE EMIT-ONLY

⛔⛔ **SUPERSEDED 2026-09-03 by `export_ledger_bind_closure_design.md` (committed `c04e071f`, two-leg-gated).**
F10's correct half — the EMIT-vs-EXPORT channel split and the `.out`-is-the-comparator-channel fact — survives
in the design's §D1/§D3. Its wrong halves are struck: **category** membership (the design uses the **bind
test** — export iff a later step binds the row's **F9 write-key**, not by label); **"carry the whole inherited
LEDGER"** (the design carries only the bind-closure, via a **generate/fold** topology over a frozen base);
**"D3 catches the closure gap"** (false — the `D3` round-trip does **not** catch under-export; the design
supplies a manifest + edge-wise closure guard); and the **export-every-primary no-manifest fallback** (a missing
manifest is now a **hard error**). ⛔ Do not build from the F10 text below — read the design. It is retained as
a record.

⭐⭐ **This makes the model-level-vs-step-level split `F9`'s "does not decide" (above) left open "for a third
step's evidence."** The third step's evidence is S11c-c1's: the full-primary export measured **145 MB**, of
which **84.7 MB was two dissipation-AUDIT rows**, exceeding GitHub's 100 MB plain-git file cap — and an
`*_exports.py` must stay plain-git (⛔ never annexed; it is a hash-chained input the next step imports). The
model-level register the next step consumes (the DtN operator + closed face response) was ~5 MB; the rest was
step-level evidence.

⭐ **Two channels, and the split is between them:**
- **EMIT** (the stdout tag stream → the annexed `out/*.out`): **every** primary object **and** every control
  operand and residual, unchanged and complete. This is what the cross-engine comparator and the review legs
  read, and the `.out` is git-annex+GIN backed, so it carries **no** size cap.
- **EXPORT** (the plain-git `*_exports.py` LEDGER → imported only by the next chain step): the **model-level
  register** — the chain-**consumed** model-definition objects (operators, kernels, closed responses, and the
  knobs/symbols a later step binds) **as NAMED by the emitting step's own spec's consume-set**, plus their
  free-symbol closure, merged onto the carried-forward inherited LEDGER. Derived **step-level diagnostics**
  (dissipation/energy audits, Hermitian/reactive splits, regime/parity/term-origin views, loci/noninvertibility
  conditions) are **emit-only** — in the `.out`, ⛔ not in the LEDGER.

⚠ **This REFINES `D1`, it does not overturn it.** `D1`'s anti-under-export default stands — the safe direction
is still to export, and the builder still makes **no per-object judgement**: it exports exactly the objects the
**spec** names in its consume-set (plus free-symbol closure plus the carried-forward LEDGER), ⛔ never a set the
builder chooses. The change from `D1` is only that the boundary is the **spec-declared consume-set**, ⛔ not
"every emitted primary object." `D1`'s real fear — a builder silently dropping something a consumer needs — is
answered by making the keep-set **spec-driven**, not builder-discretionary.

⛔ **Guards (an under-export here breaks a downstream step, so these bind):**
- The emitting step's spec **MUST declare its consume-set** (the objects the next step imports). A step with no
  declared consume-set falls back to `D1` (export every primary object). ⛔ F10 is not a licence to narrow a
  spec that did not declare one.
- The engine's `D3` in-run round-trip catches a **free-symbol-closure gap** (a kept object referencing a symbol
  neither inherited nor freshly exported).
- A diagnostic is emit-only **only** because the spec's consume-set does not name it. If a **later** step turns
  out to need one, it is promoted into that step's import then — the object is still in the annexed `.out`,
  recoverable. ⛔ Emit-only is never "computed then discarded."

⚠ **SAFE because the only LEDGER consumer is the next chain step.** The cross-engine comparator and the review
legs read the `.out` tag streams, ⛔ not the LEDGER (`S11c_b_cross_engine_comparator.py`: *"S11c_b_exports.py is
not used as a tag stream"*); `F9` is PY-ONLY and the WL engine writes no ledger — so nothing that reads the
`.out` loses anything under F10.

⭐ **Scope.** Applies going forward to any step whose spec declares a consume-set; ⛔ prior committed exports are
unchanged. **First applied: S11c-c1** (ratified by the user 2026-09-03; gated by the S11c-c1 build-directive
review — Grok + Codex — where Codex's must-fix is what elevated this from a build-directive "refinement" to this
list-level rule). ⛔ F10 does **not** authorize dropping a knob, a structural relation, or any object a later
step binds.

---

## ⭐⭐ OWED TO THE BUILD REVIEW — ⛔ NOT to this list, and ⛔ NOT to the engine's own guards

⚠⚠ **Eight review legs on this section produced these, and ⛔ every one of them is a property of a SCRIPT
THAT DOES NOT EXIST YET.** ⭐ Reviewing them as prose measures the document instead of the artifact.
⇒ ⭐⭐ **They are the checklist for the legs that review the built PY engine, where a FORM ablation can
actually be run** (rule 14). ⛔ Do not fold them back into this list.

1. ⛔⛔ **The published payload must MOVE under a FORM change to the action** — for every row, ⛔ not only
   the comparison. ⚠ Measured: a hand-typed payload survives an ablation that only deletes the comparison.
   ⭐ The obligation is `S11_SHARED_PHYSICS.md` **§4**'s *reached by computation* and its deletion test —
   ⚠ ⛔ **not** corollary 5, whose letter admits a payload copied out of a live imported row.
2. ⛔ **The class tag is policed on an `F9b`.** ⚠ S10's committed guard asserts it
   (`S10_brane_mode_spectrum_sympy_audit.py:2085`); ⛔ an object-only rule drops it, and `class` is what
   makes the accumulated ledger *the true list of knobs*.
3. ⛔ **An imported row is never silently dropped.** ⚠ Measured: **50** live `dimension_key` references into
   **34** target rows; ⛔ dropping one target makes `Q6r`'s cross-step check vacuous and `F7` reports it as
   an ordinary unresolved lookup.
4. ⛔⛔ **A `dimension_key` may not be re-pointed with the payload left intact** — ⚠ every residual reads
   zero. ⭐ This is the re-pointing failure where **every repo check passes**; ⚠ **12** such references
   already cross a `D` boundary.
5. ⛔ **Two of the step's own objects may not collide on one key** — ⚠ S10's writer raises on it
   (`:2072-2074`); ⛔ a rule quantified over *published rows* loses it.
6. ⛔⛔ **AN EMPTY INTERSECTION BETWEEN THE ENGINE'S DERIVED KEYS AND THE IMPORTED `LEDGER` IS A FINDING,
   ⛔ NOT A PASS.** ⚠ The engine emits both key sets as operands; ⭐ a builder whose key spelling never meets
   the import takes the absent branch on every row, ⛔ `F9b` never runs, and the run looks identical to one
   that met the import and agreed. ⇒ ⭐ **`F1` is what the spelling must satisfy** — ⛔ check it against `F1`,
   ⛔ never against a transform this file does not supply.
7. ⛔ **A consumer checking the comparison must use its OWN routine** — ⚠ reusing the writer's makes the
   check audit its own input. ⭐ The comparison's own obligations are in the rule above, ⛔ not here: ⚠ they
   are what the builder must DO, and a builder never reads this section.

⛔ **A leak prohibition was written here and is CUT.** ⚠ It said no key's outcome may appear *"anywhere a
builder can read"* — ⛔ unenforceable, and `DEFECT_REGISTER.md#c20` already states the value, the
membership and the count. ⭐ A denylist means the architecture is wrong: ⛔ blindness is enforced by
**absence**, by bounding what the builder is handed, ⛔ never by a sentence forbidding it.

---

## ⛔ What this list does not decide

- ⭐ The row shape of `F3` ⇒ the S11 PY decision list.
- The `C19` rename worklist ⇒ its own gate, per `F5`.
- `C17`, `C18`, S10's requirements registers.
