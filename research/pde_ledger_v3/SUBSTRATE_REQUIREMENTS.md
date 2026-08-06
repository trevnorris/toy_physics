# Substrate requirements — what the force sectors oblige the unbuilt steps to deliver

**Status: seeded, not populated.** Two entries, both found by accident. The population passes are
scheduled below.

## Why this exists

v3 is requirements-first: each sector states what it needs, and the knit happens last. The light sector
(S9, S10, S11, S11b) is built on **S1–S8, none of which have been run** — there are no step records for
them. That is the design, not an oversight: the substrate gets built once every sector has said what it
requires, so it can be checked against all of them at once rather than rebuilt per sector.

The bet only pays if the requirements are **captured**. Measured 2026-08-05: exactly two obligations are
recorded, both from S11, both about one quantity, both as inline anchors in a 1200-line plan. Everything
else a built sector assumes about the substrate is implicit prose spread across four step records.

`DEFECT_REGISTER.md` records what is *wrong*. `ANSATZ_LEDGER.md` records what is *postulated and not
derived*. Neither records **what a later step must deliver for an earlier result to stand.** That is this
file.

## The failure this prevents

At Phase 5 the sectors are knitted and the substrate is built against their combined requirements. If a
requirement was never written down, one of two things happens: the substrate is built without it and the
sector's result silently loses its footing, or it is rediscovered late and the substrate is rebuilt. Both
were survivable with one sector. With five they are not.

The S7 entry below is the shape of the problem. Standing at S7, the honest answer is "the slab width is
not selected" — a real result, bankable, and **insufficient**. Only S11 knows that the question must be
asked again with the width made position- and time-dependent. Nobody at S7 would think to ask.

## Entry schema

| field | meaning |
|---|---|
| `id` | stable identifier, cited from both ends |
| `source` | the step that needs it |
| `target` | the step that must deliver it |
| `requirement` | **the object required**, named — not a derivation path for it |
| `on failure` | what breaks in the source step if the target cannot deliver |
| `status` | OPEN · DELIVERED · RETIRED (with the commit that changed it) |

Name the object; do not prescribe how the target should obtain it. A requirement that specifies a
derivation route manufactures arguments about the route — see `CLAUDE.md` rule 3.

Distinguish carefully:
- a **requirement** is an obligation on a *future* step;
- a **postulate** is a value or form assumed *here* and not derived — that belongs in `ANSATZ_LEDGER.md`;
- a **defect** is something already wrong — that belongs in `DEFECT_REGISTER.md`.

A postulate with a named retirement condition generates a requirement. `B_comp` below is exactly that.

---

## Entries

### R-S6-01 — `B_comp` must be retired or re-affirmed

- **source** S11 (`steps/S11_stray_longitudinal.md`) · **target** S6 · **status** OPEN
- **requirement** — the brane's compression modulus `B_comp` (`Q.brane.B_comp`), as a derived quantity or
  an explicitly re-affirmed postulate.
- **on failure** — S11's knob count stays an upper bound and the longitudinal mode keeps a postulated
  constant underneath it. Entered as a postulated knob on the user's explicit call (2026-08-02) precisely
  so the retirement would be visible when it happens.
- **note** — recorded at `V3_STEP_PLAN.md` `{#s6-b-comp-callback}`, written into S6 on purpose: a note
  living only in S11's record is a note nobody reads on arrival at S6.

### R-S7-01 — the slab-width flat direction, under a non-uniform width

- **source** S11 · **target** S7 · **status** OPEN
- **requirement** — whether the slab-width flat direction survives when the width is made **position- and
  time-dependent**. A uniform-width statement does not answer it.
- **on failure** — if the flatness is a genuine flat direction, the wall offers no resistance to
  thickening, so `B_wall = 0` and by the series law `B_comp = 0`; S10's longitudinal zero is never lifted
  and **S11's propagating mode does not exist.** S11's own answer is that gradients lift it — a wave
  modulates the width, tilting and stretching the interfaces at cost proportional to `σ_wall|∇W|²`, so
  flat at `k=0` and stiff as `k²`. S7 must confirm or refute that.
- **note** — `V3_STEP_PLAN.md` `{#s7-b-comp-callback}`. The charge anchor rests on this.

---

## Population passes

**Pass 1 — S9 and S10.** Scheduled next. These two were rebuilt carefully and their records are the most
trustworthy, so they are the right place to calibrate what an entry looks like and how fine-grained to go.

**Pass 2 — S11, S11bA, S11bB, then the remaining sectors as they close.** Once the shape is settled.

Method for a pass, per step record: read it for every statement of the form *this assumes*, *this rests
on*, *this is postulated at*, *this is deferred to*, *X enters at S**n***, and for every named retirement
condition. Each becomes an entry keyed by the step that must deliver. Where the record asserts something
about a step that has no record yet, that is a requirement by definition.

Two things to watch for, both seen already:
- A requirement can point **sideways**, not only backwards — S11's questions reduce to the S11b interface
  law, which is a sibling step, not a substrate one.
- The same object can be required by several sectors. Charge and magnetism ride the same brane–bulk
  coupling as light, so an entry may gain sources rather than being duplicated.

**Honest sizing.** Pass 1 is small — four records exist. The work grows with each sector that closes, and
the whole point is that it grows *incrementally* rather than being reconstructed at Phase 5, where an
implicit assumption is exactly what gets missed.

## When this is done

Phase 5 knits the sectors and builds the substrate. At that point this file is the checklist the
substrate is tested against — the thing that answers *"is the brane–bulk we built sufficient to carry all
five forces"* by enumeration rather than by recollection.
