# Cross-step consistency by dataflow — the export chain

**Status: PROPOSAL, 2026-08-07.** ⛔ Not built. ⭐ Written for review before anything is deleted.

---

## The requirement, unchanged since v3 began

> Each step must be consistent with every other step, and each step internally dimensionally consistent.

⭐ That is the whole reason v2 was abandoned for v3.

## What is actually delivered today — measured

| ask | state |
|---|---|
| internal dimensional consistency | ⭐ **delivered.** `dimensional_homogeneity_gate.py` reports `HOMOGENEOUS=5`, and as of `b8ade918` the declared dimensions are corroborated against five engines' own symbolic derivations. |
| cross-step consistency | ⛔ **~20%, and the 20% is broken.** |

⛔ Of the **five** relations in `relations.yaml`, exactly **one** (`R4`) is wired into an executable check.
`R1`, `R2.xi_h`, `R2.h0` and `R5` have **no `registry_residual` row at all** — verified by repository-wide
search, twice, by two independent parties. ⛔ And `R4`'s single row binds a **squared** speed (S9) and an
**ω²** (S10) to a **linear-speed** quantity, so its residual is nonzero by construction and has never
closed.

⇒ ⛔⛔ **The central requirement of the v3 rebuild is essentially unbuilt.**

## Why it went this way — the root cause

⭐ The registry exists to **reconcile two independent copies** of a quantity: a step derives `μ_R`, the
registry declares `μ_R`, and a residual asks whether they agree.

⛔ **But no quantity has two independent sources.** Each is derived in exactly one place. So the residual
layer reconciles a derivation against a **hand-copied restatement of itself** — and every defect this
project has fought (a declaration drifting from a derivation, a wrong `D` coefficient invisible to the
gate, a witness needing a witness) is *created by* maintaining that second copy.

⚠ Measured cost of the second copy, this session alone: five build rounds, ten review legs, four
recurrences of the same defect one level up, and no physics.

---

## The architecture

### ⭐ Python engines: a chain

Each step's `.py` engine **imports** the previous step's export record and **exports** its own. Cross-step
consistency is then true **by construction** — there is no second copy to disagree with.

### ⭐ Wolfram engines: siloed

Each step's `.wl` engine **imports nothing**. It derives everything from the action and its stated
premises, exactly as it does today.

⚠ **This is already the case.** Measured registry references: `S9-wl` **0**, `S10-wl` **0**, `S9-py` **0**,
`S10-py` **5**, `S11-py` **2**. ⇒ the proposal formalises what the code already does and replaces the
registry with a direct record.

### ⭐⭐ Why the asymmetry is the point

⛔ A pure chain has a hole: a wrong export propagates downstream silently, because nothing has an
independent opinion.
⭐ The siloed WL side **never imported it**, so it disagrees. ⇒ **consistency by dataflow on one side,
independent re-derivation on the other, and the disagreement is the measurement.**

⚠ Note where the catch happens: a wrong value is caught **at the step that derived it**, where both
engines compute it independently — ⛔ not at every later step that consumes it. ⭐ That is sufficient, and
the record should say so rather than implying continuous re-verification.

## The export record

One record per exported object, **computed by the engine**, ⛔ never hand-written:

```
name · value · dimension (as the engine derived it) · derived|declared · source locus
```

⚠ `dimension` is **solved by the engine**, ⛔ not declared anywhere. S9, S10 and S11 already emit exactly
this — `DIM_HOMOGENEITY_ACTION`, `DIM_COEFFICIENTS`, the symbolic dimension solves. ⇒ the machinery exists;
it is not currently used as the interface.

## The three checks

| check | catches | status |
|---|---|---|
| **cross-engine agreement** — the two engines' exports match | an error in one engine | ⭐ exists — `engine_output_checks.py` |
| **per-engine homogeneity** — each engine's own expressions are dimensionally consistent under its own derived dimensions | a dimensional error, at the step that made it | ⭐ exists, inside the engines |
| **import ⊄ export** — nothing a step marks `derived` appears in what it imported | ⛔ a tautology: importing the answer, then "deriving" it | ⛔ new, ~20 lines |

⭐⭐ **The third check is what protects blind derivation**, and it is mechanical: compare two lists.
⭐ A step may import what it **consumes**; ⛔ importing what it **produces** is the defect.
⭐ WL needs no such rule — it imports nothing.

## What is declared per step

Two short lists, ⛔ and nothing else:

1. **py's import list** — which objects this step consumes.
2. **the cross-engine comparison list** — which tags pair between the two engines.

⚠ The second is the only real recurring cost. ⭐ Measured on S11: once a shared spec exists, the engines
share **190 of ~225** tag suffixes, so the pairing is mechanical — ⛔ not S10's 3,121-line hand-written
list, which was written before the spec existed.

## What is deleted

- `quantities.yaml` and `relations.yaml` **as the cross-step medium**
- `registry_residual` bindings
- `registry_dimension_witness*`
- `engine_dimension_pin.py`'s reconciliation role
- `dimensional_homogeneity_gate.py`

### ⚠ Why the homogeneity gate goes too

⭐ It was proposed for retention and that was inconsistent. ⛔ It grows linearly — ~26 lines per relation,
~25 per quantity, so 500–1500 hand-maintained lines across 20+ steps — and it checks a **restatement**
using **declared** dimensions while the engine already checks its **derivation** using **derived** ones.

⛔ Its discriminating power is measurably low. All three measured 2026-08-07:

```
wrong D-coefficient in all three brane laws  ->  HOMOGENEOUS=5  PASS
every law deleted                            ->  HOMOGENEOUS=5  PASS
common-mode shift, registry + all 5 engines  ->  HOMOGENEOUS=5  PASS
```

⇒ it constrains **differences** between declarations, and the quantity of interest cancels out of every
difference. ⭐ That is why the engine pin had to be built at all.

⭐ **Phase 5 is the one case that genuinely wants a global statement** — when sectors are knitted and light's
constants are checked against gravity's. ⛔ That is **one event**, not per-step upkeep, and it should be
built then against export records.

## What survives

⭐ Two CAS engines per step. ⭐ The cross-engine comparison. ⭐ The in-engine dimensional solves.
⭐ The ablation controls — ⚠ they test **physics**, not bookkeeping, and are the only thing that has ever
caught a typed conclusion.

## Per-step cost

| | now | proposed |
|---|---|---|
| two engines | required | required |
| `checks_S<n>.yaml` | 194–3,121 lines | comparison list only |
| `quantities.yaml` entries | ~25 lines each | — |
| `relations.yaml` entries | ~26 lines each | — |
| harness battery extension | required | required |
| import list | — | a few lines |

## Migration

⭐ S9, S10, S11 already emit the objects an export needs. ⚠ The open question is whether they emit **enough**
without engine changes — ⛔ that is a read, and it should be answered before anything is deleted.

⛔ Nothing is deleted until the replacement runs on one step and reproduces what the current checks report.

---

## ⚠ Known risks — ⭐ stated, ⛔ not resolved

1. ⛔ **A wrong export propagates.** Mitigation is the dual-engine check at the step of origin. ⇒ if both
   engines make the same error, nothing catches it — ⚠ unchanged from today, and it is why the ablation
   controls stay.
2. ⚠ **An object only one engine computes has no contention.** Today: `DIM_REGISTRY_*` is py-only,
   `*_CONDITIONS` is wl-only. ⭐ Those simply are not compared; ⛔ the record must not imply they are.
3. ⚠ **Representation divergence.** A nullspace basis is not canonical; S10 carries 11 rows reporting
   DISAGREE on representation alone. ⭐ The comparison list must exclude such objects deliberately.
4. ⚠ **The chain is only as good as the export being machine-generated.** ⛔ The moment an export record is
   hand-edited, the second-copy defect returns.
