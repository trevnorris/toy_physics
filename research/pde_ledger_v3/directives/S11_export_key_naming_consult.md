# Consult — what a storage key must say, so that only real meetings collide

⭐ **This is a targeted design consult, ⛔ not a review.** One question. Answer it with reasons and failure
modes; ⛔ do not audit anything else.

## The situation, measured

A chain of steps each import the previous step's flat `LEDGER` and export a merged one
(`research/pde_ledger_v3/scripts/S10_exports.py`, 617 rows, `{display, value, value_kind, class, step}`).

Two steps analyse **different actions**:

```
S10 MAIN:  L = (ρ_br/2)·Σ_j (∂_t u_j)²  −  (μ_R/2)·S_curl[∂u]
S11 MAIN:  L = (ρ_br/2)·Σ_j (∂_t u_j)²  −  (μ_R/2)·S_curl[∂u]  −  (B_comp/2)·S_div[∂u]
```

S10 named its exports after the **kind** of object with no mention of which system — `root_ordering_d3`,
`q1_lagrangian_expanded_d3`, `q7_stiffness_d3`. S11 computes those same kinds of object from its own
action, so a flat key is claimed by two different objects.

⭐ **The partition is measurable, not a matter of taste.** Mutating the stiffness form and re-reading the
committed output of `scripts/out/S10_brane_mode_spectrum_sympy_audit.out` at `D = 3`:

```
MAIN quantities comparable across both form controls   220
  unchanged under both mutations                       104     one object for every step
  moved under a form mutation                          116     one per material
```

The shared spec already draws exactly this line in physics terms — `directives/S11_SHARED_PHYSICS.md:604-608`
distinguishes objects "fixed by §2 and the same for every package" from those where "mutating the stiffness
functional must move them".

## The settled constraints — ⛔ an answer that violates one of these is not an answer

- **`F1`** (`directives/S11_export_chain_decisions_v2.md`): storage keys are **flat**; ⛔ **no producer
  prefix.** A later step re-deriving **one object** writes **the same key**, so that the writer compares
  them. ⭐ **The collision is the measurement** — this is the reason the chain exists.
- **`F2`**: before writing a key that exists in the import, compare the **object**. Same ⇒ re-derivation,
  emit both operands and the residual. ⛔ Different ⇒ a finding that fails loudly.
- **`D11`** (`S9_REWRITE_PLAN.md:268-281`): ⭐ **a name belongs to the object, ⛔ not to the engine that
  computed it.** Both engines emit the same name; comparison is a lookup, ⛔ never a hand-written pair
  table. ⚠ The alternative cost 3,121 lines of hand-written name pairs.
- ⭐ **The design goal: only REAL meetings should collide.** A collision between two genuinely identical
  objects is the check working. A collision between two different objects under an under-specified name is
  the defect. ⛔ Machinery that prevents all collisions destroys the check (this has already been proposed
  twice and blocked twice).

## The question

**For an object that MOVES when the action is mutated, what should its storage key say, so that it cannot
be confused with another step's object of the same kind — while an object that does NOT move keeps a bare
name and still meets its counterpart across steps?**

Two candidate schemes are on the table. ⭐ Judge them, and ⭐ propose a better one if there is one:

1. **Name the material** — the physical system the object belongs to, e.g. a qualifier derived from the
   stiffness functional (`curl_only_…` vs a name for curl-plus-compression).
2. **Name the step's subject** — the physics question the step asks, e.g. `brane_mode_…` vs
   `stray_longitudinal_…`.

## What the answer must address

- ⭐ Which scheme keeps `D11` true — that the name belongs to the **object**, so that **two independently
  built engines analysing the same material produce the same key without consulting each other**.
- ⚠ A qualifier naming the material must be **derivable from the action itself**. Say how, or say why it
  cannot be — ⛔ and if it cannot, that is a finding, not a detail.
- ⚠ Two different steps may analyse **the same material**. Under each scheme, do their objects meet? ⭐ They
  should — that is a real meeting and the check exists to catch it.
- ⚠ One step analyses **eight packages**, only one of which is exported today. Does the scheme still work if
  a later step exports more than one?
- ⛔ What each scheme does when a builder must name an object the spec does not name — the spec permits this
  explicitly (`S11_SHARED_PHYSICS.md:1084`), so the namespace **cannot be closed in advance**.
- ⭐ The failure mode of each: what wrong physics survives, and what real comparison is lost.

## Also answer

⚠ The rename applies to an **already committed** artifact: 429 of S10's 617 rows are its own tag-derived
objects. Re-pointing a name is the failure mode where **every check in the repository passes and the
physics silently moves**. ⭐ Name the controls that make a rename safe. ⛔ Do not assume the ones I have in
mind; state yours, and say what each would catch.

⛔ Read-only on the working tree. ⛔ Any script under `/tmp`, wrapped in `timeout 600`; report the absolute
paths of anything you run and its literal stdout.
