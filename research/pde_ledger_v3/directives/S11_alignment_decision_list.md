# S11 engine alignment — decision list

⭐ Both S11 engines ran to completion and their outputs are on disk. ⛔ **Both were built against the
PRE-REPAIR spec.** `directives/S11_SHARED_PHYSICS.md` was then repaired and closed at `f49a1684`
(4 rounds, 8 legs). ⇒ each engine owes **one** aligned round before anything downstream is written.

⛔ **Do not write `reduction/checks_S11.yaml` first.** ⭐ It would freeze a tag vocabulary that this round
changes.

---

## What is measured about the current outputs

⭐ Surveyed from the committed outputs, ⛔ not inferred:

| | WL | PY |
|---|---|---|
| tag lines | 3750 | 3539 |
| distinct suffixes after stripping prefix/package/D-cell | 229 | 219 |
| **shared** suffixes | **190** | **190** |
| `STRATUM1_POINT_RESIDUAL` tags | **0** | **10** |
| distinct `PREMISE_*` suffixes | **17** | **22** |
| `PREMISE_INVENTORY` | absent | absent |

⭐ Packages, identical in both: `MAIN`, `XCOEF_BSCALE`, `XCOEF_BSIGN`, `XFORM_CURLONLY`, `XFORM_DIVONLY`,
`XFORM_EXTRA`, `XFORM_TRACELESS`. ⭐ Dimension cells `D2 … D5`, identical in both.

⇒ ⭐⭐ **The shared spec worked.** The pre-spec survey measured the two engines sharing **exactly one** tag
suffix; they now share 190. ⇒ the eventual pair list is **mechanical**.

---

## A1 · The five items the repaired spec owes each engine

⭐ From `REBUILD_HANDOFF.md`, ⭐ **verify each against the spec yourself**; ⛔ do not inherit this list:

1. ⭐ Per-premise tags collapse into a **single engine-local premise inventory**. ⚠ Measured: 17 and 22
   distinct premise suffixes today, no inventory tag in either engine.
2. ⭐ `c_s0` joins the `KW_SIGN` and `KW_ZERO_LOCUS` admissibility tests.
3. ⭐ `Q10`'s unconditional pinned failure object is adopted. ⚠ The handoff records WL as already emitting
   it and PY as not; ⛔ **the orchestrator could not confirm either from the outputs** — ⭐ establish the
   real state before changing anything.
4. ⛔ **PY emits 10 `STRATUM1_POINT_RESIDUAL` tags for an object the repaired spec deletes; WL emits 0.**
   ⇒ ⚠ this is a **specification artifact**, ⛔ not a physics disagreement.
5. ⭐ Fingerprints for the headline objects at shared probe points.

## ⛔⛔ A2 · The failure this round must not repeat

⚠ **Measured in this step, 2026-08-06:** a fix directive told **engine 1** to drop `POINT_RESIDUAL` on
corollary-3 grounds without noticing that spec §Q8b **names** it. ⇒ WL emitted 0 and PY emitted 5, and the
resulting cross-engine disagreement was **manufactured by the directive**, not by the physics.

⇒ ⭐⭐ **A change to one engine's reading of a shared clause is a SPEC question first.** ⛔ Never repair one
engine against a clause without stating what the other engine must do with the same clause.
⭐ Every item above must be expressed so that **both** engines read the same requirement.

## ⛔⛔ A3 · The defect class that cost engine 1 nine rounds

⚠ **All nine of engine 1's defects were one class:** a tag that **declares what the run used**, assembled
from a literal beside the computation instead of read out of it — the stripped factor, the simplifier, the
sort key, the bulk premises, `[s]`'s dimension, the supplied action premises, and finally the tag whose only
job is to inventory the others.

⭐⭐ **NONE was visible by reading. EVERY one was visible by mutation.** ⚠ Across three rounds the ablating
leg found **8 of 9**.

⇒ ⛔ **Item 1 creates exactly such a tag.** A premise inventory is a declaration of what the run assumed.
⭐ It must be **read out of the engine's own assumption set**, ⛔ never assembled from a literal list.
⭐ **Audit it by mutation before the first review round**, ⛔ not after.
⚠ A premise stating an **absence** (`v₀ = 0`, no dissipation, frozen wall width) cannot drive a
construction ⇒ ⭐ mark that explicitly as the honest outcome; ⛔ do not manufacture a consumer for it.

## ⚠ A4 · Asymmetries to resolve deliberately, ⛔ not silently

⭐ Measured, ⭐ decide each and say which it is — engine idiom, spec gap, or physics:

- **`DIM_REGISTRY_DECLARED` / `DERIVED` / `RESIDUAL` / `SOURCE_LOCUS` are PY-only.** ⇒ S11's registry
  comparison is **single-engine**. ⚠ It also subtracts a `D = 3` specialisation rather than the registry
  law — ⭐ a separate change owns that; ⛔ do not fix it here, ⭐ but do not entrench it either.
- **WL-only `*_CONDITIONS` tags** — a `Reduce`/`Solve` idiom with no PY counterpart.
- **`L` (PY) vs `LAGRANGIAN` (WL)**, **`EULER_LAGRANGE` (PY) vs `EULER_LAGRANGE_SYSTEM` (WL)** — ⭐ same
  object, different spelling. ⛔ A naming difference is not a disagreement.

## ⚠ A5 · The probe set is the artifact that can go quietly wrong

⛔ A probe point where a headline object **degenerates** makes the comparison vacuous while reporting green.
⭐ Requirements: several **distinct exact rational** points; values satisfying the spec's positivity
premises so each point is admissible; per-`D` tables for `D ∈ {2,3,4,5}`; and ⭐ a check that **flags any
headline object evaluating to zero at every probe point.**
⛔ **Never floats** — that is the 1989 caution and it is the named failure mode of this comparison method.

## ⛔ A6 · Out of scope

⛔ Do not write `reduction/checks_S11.yaml`. ⛔ Do not touch `reduction/` at all — it is under concurrent
change. ⛔ Do not modify S9's or S10's engines, records, or committed outputs. ⛔ Do not commit.
⛔ **Do not declare `N6_BASIS` comparable.** ⚠ A nullspace basis is not canonical; S10 carries 11 rows
reporting DISAGREE on representation alone, and the spec already marks `N6` display-only.

## ⚠ A7 · A cross-step comparison exists and is the ORCHESTRATOR's to run

⭐ `XFORM_CURLONLY` is the same physical system S10 computes as its `MAIN` action. ⇒ once both engines are
aligned, four independently built engines will have computed one system — which this ledger has never had.
⛔⛔ **This must not reach a builder.** ⚠ Pointing a builder at committed S10 rows supplies a target it can
converge on; the cross-step language was **deliberately removed** from the spec for that reason.

---

## Rules for whoever writes the directive

- ⭐ Name objects. ⛔ Every time you are about to write *how*, write *what* instead.
- ⛔ No expected value, no residual, no acceptance number anywhere a builder can read.
- ⭐ Each engine gets its **own** build from the shared spec plus a short engine-specific directive; ⛔ no
  builder reads the other engine's change.
- ⚠ ⭐ **Verify this list against the spec and the outputs before using it.** ⛔ If an item is wrong, say so
  and stop — ⭐ that is worth more than the deliverable. ⚠ An earlier decision list from this orchestrator
  carried a false measurement and was correctly refused.
