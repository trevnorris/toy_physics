> ⛔⛔ **DEFERRED 2026-08-05 alongside `S10_spec_rewrite_directive_v2.md`. ⛔ NOT APPLIED.**
> ⭐ This is the fold of **two independent review legs** on that directive — ⭐ keep it, it is the expensive
> part. ⚠ A partial application was started and **discarded**; the reviewed baseline is commit `200dcf63`.
> ⚠ ⭐ **Three items in here are NOT contract-specific and were pulled forward into S10** — the missing
> distinctness premises (`F3`), the supplied normalisation that kills a two-route check (`F2`), and the
> honest scoping of the extraction residual (`F4`). ⇒ `S10_spec_Q7_repair_directive.md`.

# Fix `S10_spec_rewrite_directive_v2.md` — ⭐ a DECISION LIST, ⛔ not a rewrite

**File:** `/var/projects/toy_physics/research/pde_ledger_v3/directives/S10_spec_rewrite_directive_v2.md`
⛔ **Edit in place. ⛔ Do not commit. ⛔ Do not modify any other file.**

⭐ You wrote this file. ⭐ **Two independent legs reviewed it and both said "do not hand to the builder."**
⚠ Both agreed the core idea is sound and worth keeping — ⭐ the two-residual split, the CAS-computed
Levi-Civita reference, the pre-simplification weight guard, the `LOCAL_` resolution, and a closed
schema-checked format that **killed all six of the attacks you named yourself.**
⛔ **Do not redesign. ⭐ Close the items below.**

---

## ⛔⛔ F1 — TWO ATTACKS GOT THROUGH. ⭐ Both re-open a defect this directive exists to prevent.

**F1a — the REFERENT GAP.** A leg left `q7` **byte-identical** to the oracle and redefined its
*referents*: set the package weight to `1` and folded the sign and the rational factor into the density.
⛔ **Accepted, no rejection.** ⚠ That is precisely the earlier "drop the sign and rational factor" defect,
re-entering by a new door.
⭐ **Root cause:** only `AUTHOR_Q7_ORACLE` carries `BEGIN`/`END` markers, so acceptance step 2's
*"apply the same AST comparison to the action, Route A, and dimension definitions"* has **nothing to
compare against**, and the weight and density columns are never named there.
⭐ **Fix:** add `AUTHOR_PACKAGE_TABLE`, `AUTHOR_ROUTE_A`, `AUTHOR_DIMENSION` oracle blocks in the same
canonical form, and make step 2 compare **all of them** byte-wise.

**F1b — the FORBIDDEN-SHAPE MATCHER IS LITERAL.** The banned quotient/reconstruction shape was
reintroduced with the `multiply` operands **swapped** and passed, because the pattern is matched as a
literal array. ⭐ **Fix:** apply the contract's canonical operand ordering to **both** the pattern and the
candidate before matching, or match modulo commutativity. ⚠ Canonical ordering is currently required only
for payload serialisation, ⛔ not for the authoring-time matcher.

## ⛔⛔ F2 — a supplied constant SILENTLY KILLS a two-route check

The directive orders a bare rational normalisation into the spec for Route A's matrix.
⛔⛔ **With it supplied, the two routes' difference is identically the zero matrix** — ⇒ the cross-route
residual becomes a **guaranteed pass**, and it is guaranteed by a number we handed over.
⚠ **And its correctness depends on an Euler–Lagrange sign convention this directive never writes.**
⭐ **Fix:** ⛔ supply **no** bare rational. ⭐ **Write the Euler–Lagrange sign convention as an equation**
and let the normalisation **follow from it by computation.** ⇒ the two routes stay genuinely two.

## ⛔⛔ F3 — a missing premise KILLS TWO CONTROLS

⚠ **Measured:** the retained premise families omit the **distinctness** conditions for the anisotropy scale
and the coefficient scale. ⛔ With the degenerate value admitted, the anisotropic package's action is
**identical to the main package's**, and the same collapse hits the coefficient-scale package.
⇒ ⛔ **two of the six controls are dead.**
⭐ **Fix:** add both distinctness premises to the retained set. ⚠ The existing instruction *"do not add a
premise to force a solver to decide"* is being read as licence to omit them — ⭐ **say explicitly that a
premise which keeps a CONTROL DISTINCT is required, and is not the same thing as a premise that forces a
decision.**

## ⛔ F4 — ⭐ SCOPE THE EXTRACTION RESIDUAL HONESTLY. ⛔ Do not let it overclaim.

⚠ **Measured, and it is the most important correction here:** for **every row of your own supplied table**
the three reports satisfy the decomposition **identically**. ⇒ they are independent as **AST owners**, ⛔ but
**algebraically dependent**. ⇒ ⛔ the extraction residual is **zero on any conforming build**, so it can
catch a **transcription slip** and ⛔ **never a wrong shared formula** — which is the failure mode the
exercise exists to defend against.
⭐ **The reference residual is the only diagnostic that can fail on physics**, because the Levi-Civita side
is **CAS-computed** while the density side is **typed**.
⭐ **Fix:** say this in the directive **and** require the rewritten spec to say it. ⛔ Do not delete the
extraction residual — ⭐ it still catches transcription and it is cheap — ⛔ but ⛔ **never describe it as
policing the physics.**
⚠ ⭐ **Also record the limit that cannot be fixed here:** the contract can enforce independence in the
**graph**; ⛔ it **cannot** force it at engine **runtime** — if both engines internally compute the density
as the quotient, the extraction residual is zero regardless. ⭐ State it as a known limitation.

## ⛔ F5 — the frozen oracle is UNCONSTRUCTIBLE, and it contradicts its own model

- **F5a** The Levi-Civita vector node juxtaposes the symbol and the gradient reference with **no
  `multiply` node**, binds only two of three indices leaving one **unbound with no binder or vector
  constructor**, and refers to an **indexed family** through a leaf grammar that has no indexed-reference
  form. ⇒ ⛔ it cannot be built. ⭐ Add an indexed-leaf form and an explicit binder, or lift the
  "may not edit" clause for this node.
- **F5b** The oracle references ids that are **per-package**, **per-stratum**, or **pattern variables**,
  while the model demands one globally unique id with exactly one owner and rejects dangling references.
  ⇒ ⛔ **acceptance steps 1 and 2 cannot both pass.** ⭐ Add a template/parameter binding mechanism and an
  explicit exemption for pattern variables, ⭐ both stated **in** the model.
- **F5c** Several node keys used by the oracle are **outside** the model's own exhaustive node list.
  ⭐ Extend the list to name them.
- **F5d** A three-element list of names is syntactically indistinguishable from a composite expression
  under the closed grammar. ⭐ Disambiguate it.

## ⛔ F6 — the builder GRADES ITS OWN ARTIFACT

⚠ The acceptance validator is authored by the same builder whose output it judges, and its schema lives in
the builder's report rather than in a locked artifact. ⇒ ⛔ **a check whose trigger lives inside the thing
it checks is that thing agreeing with itself.**
⭐ **Fix:** either embed the full schema and allowlists **in this directive** as oracle blocks, or require a
**fixed external check script at a stated path that the builder may not edit.** ⭐ Prefer the second.

## ⛔ F7 — close the underspecified schemas

`classifiers`, `quantity_registry`, and `report_contract` are demanded as part of a **closed** model while
their own shapes are left to the builder. ⭐ Specify each, ⛔ or drop it from the closed model. ⚠ Leaving
them open is what turns a schema into builder thrash.

## ⛔ F8 — the oracle carries RATIONALE and HISTORY into a document that forbids both

The frozen block carries purpose/rationale fields and a record of a previously forbidden shape into the
target spec, whose own rules forbid exactly that.
⭐ **Fix:** move both into the **authoring-time validator**, where they belong.

## ⭐⭐ F9 — SERIALISATION: use YAML, ⛔ not JSON

⛔ **Standing project rule: no JSON for LLM-consumed artifacts — YAML or markdown; JSON is
machine-to-machine only.** ⚠ The rewritten spec is read by the two engines' **builders**, ⇒ it is
LLM-consumed by definition.
⭐ **Your guarantee does not depend on JSON** — it comes from a **closed schema plus unique definition
ownership**, and YAML gives both. ⭐ Convert the contract and every oracle block to YAML, keep canonical
ordering and byte-comparison, ⛔ and keep every equation human-legible inside it.

## ⚠ F10 — a classifier that manufactures cross-engine disagreement

⚠ **Measured:** one engine's CAS lacks real cylindrical algebraic decomposition, so the sign classifier
returns the undecided status where the other returns a determined one — ⇒ ⛔ **a systematic cross-engine
status difference with NO physical cause**, which would be read as a finding.
⭐ **Fix:** the status must record **which CAS capability was exercised**, and ⭐ a status difference
attributable to solver capability must be **excluded from the disagreement count** and reported in its own
category. ⛔ Do not force either engine to guess.

## ⚠ F11 — non-blocking, ⛔ but state them as accepted risks in the directive

- The legacy-name inventory is **thousands** of names to classify with no payload inspection — ⭐ a scale
  risk against a short handoff. ⭐ Say how it is staged.
- A number of committed tags carry no scope token and need one invented.
- One engine places a scope token inside the quantity field where the grammar puts it earlier — ⭐ name it
  as an engine obligation.

---

## ⛔ WHAT NOT TO CHANGE

⛔ The settled `§Q7` physics. ⛔ The two-residual split. ⛔ The CAS-computed reference. ⛔ The
pre-simplification weight guard. ⛔ The `LOCAL_` resolution. ⛔ The twelve harness config rows you named —
⭐ a leg verified them as exactly right.

## Report back — ⛔ under 40 lines

1. One line per `F1`–`F11`: fixed / partially / not, with line numbers.
2. ⭐ **Re-run the two attacks that succeeded and paste the rejection category for each.**
3. ⭐ Anything in this list you believe is **wrong**. ⭐ This is wanted.
4. ⛔ Do not report any engine value or the direction of any engine comparison.
