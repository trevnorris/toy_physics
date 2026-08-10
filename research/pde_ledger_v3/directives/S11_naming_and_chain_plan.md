# Plan — naming, the chain, and the order S11 gets built in

**Status: PROPOSED, 2026-08-10. ⛔ Nothing here is applied.** Written for review before execution.

Three defects surfaced while trying to launch S11's PY engine. They interact, and the order in which they
are fixed changes what has to be rebuilt. ⭐ This proposes one route and states three alternatives.

---

## What was measured

⭐ Every number below came from a command against a committed artifact.

| measurement | result |
|---|---|
| `S10_SHARED_PHYSICS.md:195-199` names the route matrices | **`M_A`**, **`M_B`** |
| what **both** S10 engines emit | ⛔ `Q2_MATRIX_A`, `Q2_MATRIX_B` |
| `Q2_MATRIX` in `S10_py_directive.md` / `S10_wl_directive.md` | ⛔ **0** and **0** — handed to neither builder |
| S10 py tag suffixes · wl · shared | 4233 · 2983 · ⛔ **562** |
| S10 `MAIN` `D3` quantities · S11 spec quantities · shared strings | 304 · 80 · ⛔ **1** |
| that one shared string | `ROOT_ORDERING` — ⚠ and the two steps mean **different objects** by it |
| `S10_exports.LEDGER['root_ordering_d3']` | `class DERIVED`, `(0, mu_R*(k1²+k2²+k3²)/rho_br)` — S10's root spectrum |
| `ROOT_ORDERING` in `S10_SHARED_PHYSICS.md` | ⛔ **absent** — S10's engine coined the key |
| entries in `S10_exports.LEDGER` carrying a `dim` field | ⛔ **0** (`S9_REWRITE_PLAN#D5`'s sketch was never built) |

⇒ `DEFECT_REGISTER.md#C19` (names are the engines', not the spec's) and `#C20` (one accidental cross-step
key collision).

---

## ⭐⭐ The thing that decides the plan

⚠⚠ **Unifying the vocabulary across steps would turn one collision into many.**

S11's `MAIN` is a **different action** from S10's `MAIN` — `(μ_R/2)S_curl + (B_comp/2)S_div` against
curl-only. Both steps run a similar question list, so both compute a `DET_M`, a `ROOT_ORDERING`, a
`COEFFICIENT_ORDERING`. ⇒ ⭐ **the two vocabularies are currently disjoint by accident, and that accident is
the only thing preventing widespread key collisions.** Make both engines emit the spec's names and the
collisions arrive at scale.

⇒ ⭐⭐ **These are two different namespaces and must be decided separately:**

| namespace | scope | what it is for |
|---|---|---|
| **emitted tag name** | py ↔ wl, **within** one step | the cross-engine comparator's join key |
| **`LEDGER` export key** | **across** steps | the chain's carry-forward and overwrite key |

⛔ Conflating them is what produced the S11 PY decision list's critical defect: a rule that made the export
key a lowercased tag name, which is correct **within** a step and collides **across** steps.

---

## ⭐ The proposed route — ⚠ REVISED after two review legs, 2026-08-10

⛔⛔ **The first draft of this section recommended deferring S10's retrofit. That is REVERSED.** ⚠ Both
legs, and the adjudication below, moved it. ⭐ The superseded reasoning is in `git show 20f57adf`.

**1 · Producer-scope every newly computed or recomputed `LEDGER` result.** ⭐ Imported rows this step does
not touch are preserved unchanged. ⇒ ⭐⭐ **collisions become impossible by construction**, and the count
that the two legs disputed stops being load-bearing.

**2 · ⛔ A re-derivation is NOT an overwrite. Emit both operands and the residual.** ⚠⚠ `S9_REWRITE_PLAN#4`'s
*"overwrite what it derives"* **destroys operand A and leaves agreement ASSERTED** — ⛔ that is a direct
violation of `CLAUDE.md` rule 2, in the chain's own carry-forward mechanism. ⭐ Cross-step corroboration
becomes an explicit `(previous, current, residual)` triple.

**3 · The spec's name is the emitted name, in both engines, at every step.** ⚠ `D12`'s mechanism stands
(⭐ emitted strings only; internal `rhoBr`/`muR` untouched) but ⛔ **`D12`'s DIRECTION is wrong** — the `.wl`
is not the deviant one. ⚠ Injectivity across the worklist **first**.

**4 · Retrofit S10 BEFORE S11's engines are built — ⭐ and the reason is DEPENDENCY, ⛔ not collisions.**
⭐ S11's PY engine **binds to S10's export keys**. Retrofitting those keys afterwards invalidates every
binding and forces S11 to be rebuilt. ⇒ ⭐ do it once, first.

**5 · Then S11 PY → WL → comparator**, its decision list rewritten against 1–4 and its legs' five findings,
and re-legged.

**6 · S10's record must disclose** that its emitted names were the engines', not the spec's — ⭐ needed
whether or not the retrofit lands, because today a reader cannot get from the record's `M_A` to the
comparator's `Q2_MATRIX_A` without knowing an undocumented convention.

---

## ⭐⭐ What the legs measured, and where they DISAGREED

⚠⚠ **The two legs contradicted each other on the load-bearing question.** ⭐ That disagreement, not either
verdict, is what settled the plan.

| | Grok | Codex |
|---|---|---|
| proxy for S11's emissions | the **old S11 engine's committed output** | **S11's spec text** |
| current colliding keys | **8** | **1** |
| verdict | ⛔ Alternative A alone is wrong | ⭐ Alternative A is correct |

⭐ **Adjudicated by measurement — Grok's count is right and Codex's is an artifact.** S11's repaired spec
names `ROOT<r>_N2_RANK` / `_N2_NULLITY` / `_N7_BASIS_COUNT` (`:354`, `:359`); S10 already emits the concrete
`ROOT1_N2_RANK`, `ROOT1_N3_STACKED_RANK`. Codex's extraction pattern **misses placeholder-indexed
families**, so it never saw them. ⇒ the collisions are **systematic**, ⛔ not accidental, and they arrive
with the **new** engine because the names come from the spec.
⚠ ⭐ **A third `C19` instance fell out of it:** `S10_SHARED_PHYSICS.md` names those `N`-families **zero**
times — S10's engines coined those too.
⭐ Codex's **methodological** caution is still right: the full post-retrofit count is **undecidable** without
the rename map, since `S11:940-942` lets an engine name any object the spec does not.

⇒ ⭐⭐ **Under item 1 the count is moot either way, which is why item 1 leads.**

## ⛔ Where this plan's first draft was WRONG — ⚠ both errors favoured its own proposal

1. ⛔ **Alternative A was straw-manned.** *"A large mechanical rename of 4233 + 2983 tags"* quotes **output
   instances**, ⛔ not rename work: the engines carry **106** `emit(` and **137** `emit[` call sites, a
   subset of which change. ⚠ Output cardinality presented as work cardinality.
2. ⛔ **Alternative D's cost was overstated.** *"Deletes the chain's only demonstrated corroboration"* is
   false — ⭐ producer-scoped keys **retain both operands** and permit an explicit comparison. ⚠ What
   disappears is destructive overwrite-as-comparison, which item 2 shows is a **defect**.

⚠ ⭐ **Both errors ran the same way — inflating the alternatives against the route I had chosen.**

---

## ⛔ The alternatives, as reviewed

**A · Retrofit S10 first.** ⭐ **Adopted, as item 4** — ⛔ but on dependency grounds, ⛔ not the collision
argument the first draft used, and ⛔ only when paired with items 1–2.

**B · Forward-only** (S11 emits spec names, S10 keeps its convention).
⛔ *Rejected:* S11 binds to S10's keys, so the retrofit debt becomes a forced S11 rebuild.

**C · Minimal** (`C20` detector only). ⛔ *Rejected:* leaves `C19` live and undisclosed.
⚠ ⭐ **The detector is still required** even under item 1 — naming discipline cannot prevent a future
accident.

**D · Step-scope everything.** ⭐ **Substantially adopted as item 1**; ⛔ the first draft's objection to it
was wrong (see above).

**E · Codex's synthesis** — retrofit first · preserve untouched rows · producer-scope new results ·
lifecycle replacement as explicit old/new operands **plus a residual** · build S11 only against the
regenerated export. ⭐⭐ **This is the route, and items 1–5 above are it.**

⭐ **Not amending `S10_SHARED_PHYSICS.md` was the right call, and both legs agree** — ⛔ amending it would
make engine convergence its own naming authority and corrupt a correct specification.

---

## ⛔ What this plan does not decide

- ⛔ The producer-scope **key syntax**, and which rows are exempt as genuinely-carried objects.
  ⚠ That is the retrofit's decision list, and it is where a wrong rule silently corrupts the chain.
- ⛔ `C17` and `C18`, which stay open and are repaired by nothing here.
