# Script rebuild — state and what's next

Read `CLAUDE.md` first (how we work), then this (where we are).
Process detail: `docs/development_pipeline.md`. Defects: `DEFECT_REGISTER.md`.

Last updated 2026-08-06.

---

## The scope, corrected — three steps need rebuilding, not one

The measured defect was in **three independently-built steps**: engines emitting physics conclusions as
typed sentences with no CAS object behind them, missed by eight review legs. S9 and S10 have been rebuilt.
**S11, S11b-A and S11b-B were built under the same pattern and are not trustworthy until rebuilt.**

| step | state |
|---|---|
| S9 | rebuilt |
| S10 | rebuilt; harness closed 2026-08-06; record still to write |
| S11 `stray_longitudinal` | built under the old pattern → **rebuild** |
| S11b-A `interface_response` | built under the old pattern → **rebuild** |
| S11b-B `interface_assembly` | built under the old pattern → **rebuild** |
| S11b-C non-uniform coupling | **never built** |

The debt is bounded and it is exactly these three. S12 onward does not exist yet, so every remaining step
in every remaining sector gets built under the new pattern the first time. This tax is paid once.

**And the rebuilds now have an instrument.** S9 and S10 were rebuilt largely blind — the harness could not
produce a report on S10 at all, so cross-engine comparison never ran. S11 onward is rebuilt against a
comparator that compares the action, scores coverage against a declaration, and cannot report success on
an empty layer.

### The pattern each rebuilt step must follow

1. Both engines emit **CAS objects** as tagged payloads. No typed conclusions; a residual asserted zero
   always prints `0` and carries no information — emit both operands and the residual, then guard.
2. A **`reduction/checks_S<n>.yaml`** declaring what must be compared: cross-engine pairs, the control
   matrix, dimension sources per engine and package, and action families where the step has an action.
   `reduction/` currently holds only `checks_S9.yaml` and `checks_S10.yaml`. This is the mechanical
   tie-in to the harness.
3. The harness runs and the **step record cites what it measured**, not what the author concluded.
4. **Requirements-register entries** fall out of the same read ⇒ `SUBSTRATE_REQUIREMENTS.md`.

S10 is the worked example of all four. Its record, written next, is the template the three rebuilds copy.

---

## Where S10 stands

| | state |
|---|---|
| **The harness** | **Done**, `bad20207`, six review legs. Runs on both configs. |
| **The action comparison** | **Done.** 13 Lagrangian + 13 Euler–Lagrange rows AGREE, verified twice independently. |
| **SymPy's Q7** | **Repaired**, `bad20207`. Emits the density its action used; nonzero for exactly the two non-curl packages. |
| **The `g` symbol map** | **Open.** No `q7_stiffness` row is compared on any package. Both legs endorse declaring it — see below. |
| **Q7 spec repair** | **Open.** `directives/S10_spec_Q7_repair_directive.md` (`fcda865a`), three review rounds done. |
| **The step record** | **Not written.** Must carry the disclosures below. |
| **TeX card** | Not written. |

---

## What's next, in order

1. **Declare the gradient-symbol map** `g{i}{j} ↔ g{i}x{j}`. Both legs judge it justified. Source it from
   the two engine definitions **with line numbers** — WL `gradient[[coordinateIndex, fieldIndex]]` with
   `gradientSymbols[[row,column]] = g{row}x{column}` (`.wl:147-151,1304-1311`); PY `gradient[i][j] =
   diff(fields[j], xvec[i])` mapped to `g{i+1}{j+1}` (`:309`, `:1539-1543`).
   ⛔ **Do not source it from the payloads.** Every Q7 object on both engines is invariant under
   transposing the index pair, so the payloads cannot discriminate the correct map from its transpose.
   State that limitation inside the declaration.
2. **Q7 spec repair** — build from `fcda865a`, then two legs. The shared spec still tells both engines what
   Q7's residual is expected to be, which is a leak into the one artifact both engines read.
3. **Write the S10 step record** — orchestrator's job — **together with pass 1 of the requirements
   register** (S9 and S10). They are the same reading job; doing them separately means reading the same
   records twice. This record is the template for the three rebuilds.
4. **TeX card** — Codex, with its own two legs.
5. **Then the rebuilds:** S11, S11b-A, S11b-B, each to the four-point pattern above. Then **S11b-C**, which
   has never been built and is the differentiator from MacCullagh.

---

## Disclosures the S10 step record must carry

Measured this session unless noted. These live only here and in `bad20207`'s message; they are expensive
to reconstruct.

- **Cross-engine agreement did not improve.** `DISAGREE 32→28`, `NAMING_MISMATCH 17→21`, coverage numerator
  `541→537`. Four rows moved from "compared and wrong" to "not compared". Do not write the repair up as
  improving agreement.
- **Q7 compares FORM only.** It drops the stiffness sign and coefficient, so `XFORM_SIGNFLIP`'s payload is
  byte-identical to `MAIN`'s while its action's stiffness term carries the opposite sign (recovered ratios
  `−mu_R/2` vs `+mu_R/2`; `XCOEF_SCALE` `−mu_R·s/2`). A zero residual is a claim about the **form**, never
  about the **term**. Both engines share the convention.
- **No `q7_stiffness` row is compared on any package**, pending item 1.
- **S9's dimension layer reported `checked=1219 homogeneous=1219`.** Rebuilt, that is **312** actual
  comparisons. S9 is a signed-off step; its cited figure was ~4× weaker than recorded.
- **Cross-engine agreement rests on declared spellings.** `lambda_scale→lambdaScale` moves 30 rows,
  `D↔braneDimension` 13, `s↔coefficientScale` 14 (8 from DISAGREE). XCOEF_SCALE's agreement rests entirely
  on the third; present it as a declaration.
- **The 26 action agreements are real** — verified twice independently, once by rebuilding the comparator
  with a stricter key, once by recovering each package's density from the emitted Lagrangian. Single-
  difference controls (order, variable, function, argument, coefficient, form) each still reach DISAGREE.
- **`XFORM_SIGNFLIP` and `XFORM_ANISO` are not form controls** despite the prefix — they vary a coefficient
  under the `curl` form. Only `FULLGRAD` and `DIVONLY` change the stiffness functional.
- **Prior art.** The light-sector algebra is MacCullagh's rotational aether (1839); curl-only stiffness
  giving `D−1` transverse modes plus a zero-restoring longitudinal is textbook (Whittaker Ch. V). **Cite
  it.** Prior art is an oracle, never a premise.
- **S10 linearises about `v₀ = 0`, which *is* MacCullagh's regime.** So agreement with him is agreement
  inside his own domain of validity and is **not** evidence our medium differs from a 19th-century aether.
  This is the most over-claimable statement in the sector.
- **What is genuinely ours**, for the Departure field: the `D`-sweep and that the bulk never enters (he had
  no bulk); the in-plane/out-of-plane split; the longitudinal kept deliberately as the charge anchor. And
  for each of the three things that killed the elastic-aether programme — longitudinal modes, preferred
  frames, matter–aether coupling — what this model does instead, and whether it escapes or relocates the
  problem. S11 establishes that all three reduce to the S11b interface law.

---

## Deferred deliberately — kept, not lost

- **`NAMING_MISMATCH` cannot fail on a bare-symbol payload**, so the 13 route-name rows
  (`M_B` / `quadraticFormRoute`) report the same verdict whether the engines took one route or two.
- **The bijection's own verification line is untested** — deleting it leaves the suite green.
- **`OpaqueDerivative` is a `sp.Symbol` subclass**, so same-text instances share cached state and
  `_map_tree` re-stamps in place. Unreachable with today's declarations; **must be fixed before any field
  or coordinate spelling is ever declared.**
- **`DEFECT_REGISTER.md#f7`** — the comparison kernel matches truthiness between booleans and numbers. A
  tripwire reports both selected and tree-level exposure instead of the fix.
- The contract rewrite, §8's registry, the sign classifier, Q6d.

---

## Known limits of the instrument

- The canonical derivative key uses `sorted(set(arguments))`, so a repeated coordinate in a dependence
  list is not distinguished. Multiplicity is not physical, so this is recorded rather than open.
- A stiffness form the gradient substitution cannot cover has no completeness guard — it dies in the
  dimension walker after printing a malformed tag.
- Exit code carries no signal on either config; both exit 2 in the healthy state. Judge by counters.
- Runtime margin is thin: modest form ablations pushed both engines past 600 s, always in the same package.

## Pinned

- **`s10-as-built`** (`9309da70`) — the spec and both engines as they ran before this rebuild. ⚠ SymPy's
  engine and output have since changed (`bad20207`); the record must cite the **new** blobs for anything
  Q7-related.
- **`wip-2026-08-05-unreviewed`** (`92461853`) — two builds committed without review and reverted.
  Reachable, prior art only, nothing in it is known-good.
