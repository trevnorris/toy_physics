# Cross-step consistency by direct import — v2 of this proposal

**Status: PROPOSAL, 2026-08-08.** ⛔ Nothing built, nothing deleted.
⚠ **v1 of this document was rejected by two legs.** ⭐ Their findings are answered inline; ⛔ the parts they
were right about are not re-argued.

---

## The requirement, unchanged since v3 began

> Each step consistent with every other step · each step internally dimensionally consistent ·
> **and every value that is not derived tracked as a free knob, with what would retire it.**

⚠ **The third clause was in the original ask and v3 never built it.** Measured 2026-08-08:
`ANSATZ_LEDGER.md` exists only under `research/pde_ledger/redteam_adversarial/` — the **old** ledger.
`quantities.yaml` carries `kind: parameter` ×6, which is a **type**, ⛔ not a provenance.
`relations.yaml` carries `provenance_status: DERIVED-EXECUTED` ×5, which describes **relations**, ⛔ not
knobs. ⇒ **v3 has no free-knob register.** The registry was meant to be one and is not.

## What is delivered today

| ask | state |
|---|---|
| internal dimensional consistency | ⭐ delivered |
| cross-step consistency | ⛔ one of five relations wired, and that one **miswired** — `R4` binds a squared speed (S9) and an ω² (S10) to a linear-speed quantity, so its residual is nonzero by construction. `R1`, `R2.xi_h`, `R2.h0`, `R5`: **no row at all.** ⭐ Confirmed independently by two legs. |
| free-knob register | ⛔ **does not exist** |

## ⭐⭐ The diagnosis, corrected

The registry was doing **two jobs** and failing both:

| job | this proposal |
|---|---|
| **data medium** carrying values between steps | ⭐ **delete it — use direct Python import.** No file, no schema, no intermediary. |
| **free-knob register** — what is *not* derived | ⭐ **build it.** Small, genuinely declarative, and the one thing that *should* be hand-maintained, because it records a judgement rather than a copy of a computation. |

⚠ v1 of this proposal replaced the medium with an **export record**. ⛔ Both legs killed that: no engine emits
the required fields (`named c_gamma/c_L: 0`, `class field: 0`, `source_locus: 0` in five of six outputs,
and S11-py's 18 are **copied from the registry**), and the adapter needed to bridge it would be *"the
unacknowledged third list that recreates the copied-declaration problem."* ⭐ Direct import needs no record.

---

## The architecture

### ⭐ Python: chain by direct import

```python
# S9_light_requires_shear_sympy_audit.py
def derive():
    ...
    return {...}                 # live SymPy objects

if __name__ == "__main__":
    emit_all(derive())           # tag output unchanged

# S10_brane_mode_spectrum_sympy_audit.py
from S9_light_requires_shear_sympy_audit import derive as derive_S9
```

⭐ Cross-step consistency is then true **by construction** — there is no second copy to disagree with, and
no serialisation, canonicalisation or CAS-equality problem, because the objects never leave Python.

### ⭐ Wolfram: siloed

Each `.wl` engine imports nothing and derives from the action and its premises. ⚠ **Already true** —
measured registry references: `S9-wl` 0, `S10-wl` 0, `S9-py` 0, `S10-py` 5, `S11-py` 2.

### ⭐ The free-knob register

One entry per quantity **no step derives**: the symbol, its dimension, why it is free, and **what would
retire it**. ⚠ `c_s0` is the worked example — a premise in *both* engines (`bulkSpeedPremise = cs0 > 0`;
`sp.Symbol("c_s0", positive=True)`), derived nowhere.

⭐⭐ **This corrects a mis-framing in v1.** A leg reported premise-only objects as *"unguarded by the
asymmetry."* ⛔ They were never meant to be guarded by cross-engine contention — ⭐ they are **free knobs**,
and the check on them is that they are **declared**, not that two engines agree.

---

## The checks

| check | catches | status |
|---|---|---|
| **cross-engine agreement** | an error in one engine, at the step that made it | ⭐ exists |
| **per-engine homogeneity**, from each engine's own derived dimensions | a dimensional error, at the step that made it | ⭐ exists |
| **derived-object provenance** — an object claimed `derived` traces to **this step's own action / EL path** | ⛔ importing what determines the answer, then "deriving" it | ⛔ new |
| **knob completeness** — every free symbol not derived in-step and not imported is on the register | ⛔ a silently hardcoded value | ⛔ new |

⚠ **The third check replaces v1's `import ⊄ export`, which both legs killed as a spelling check.** Name
disjointness does not block *functional determination*: import `μ_R` and `ρ_br`, "derive" `√(μ_R/ρ_br)`
with no spectral work, and the names never collide. ⭐ Provenance to the step's own action is the property
that actually distinguishes them; ⚠ **it is the harder thing to build and should be scoped before it is
promised.**

## What happens to `relations.yaml`

⭐ **It survives, with its job narrowed**, and this answers a leg's objection to deleting the gate.

- ⛔ `R4`, `R5` **leave it.** Those quantities *are* derived; under direct import they are dataflow, and a
  residual against a restatement is the second copy this proposal removes.
- ⭐ `R1`, `R2.xi_h`, `R2.h0` **stay.** No v3 step derives any medium quantity — so these are **constraints
  among free knobs**, saying those knobs are not independent. ⇒ that is a real and checkable statement,
  ⛔ and it is not dataflow.

⇒ ⭐ `relations.yaml` becomes the **knob-constraint document**, and the homogeneity gate keeps a narrowed,
honest claim over it. ⚠ It grows only when a new *constraint among undetermined quantities* appears —
⛔ not once per step.

---

## Costs — ⭐ stated at the size the legs measured, ⛔ not the size v1 claimed

1. ⚠ **One refactor per Python engine** — wrap the derivation in a function, keep `__main__` unchanged.
   ⭐ Done once per engine. ⛔ It is an engine change and needs the usual gates.
2. ⚠ **Upstream re-runs downstream.** S11 → S10 → S9 is roughly **37 minutes** measured
   (S9-py fast, S10-py 165 s, S11-py 1855 s). ⭐ Tolerable; cacheable later if it bites.
3. ⛔⛔ **The cross-engine comparison list does NOT get short.** ⚠ A leg measured why: S10's 3,121 lines are
   **structural** — 13 cells × ~73 basenames = 690 rows — ⛔ not naming debt. S11 shares 190 of ~225
   suffixes, which proves **naming overlap, ⛔ not comparable semantics**: 68 types remain unmatched, and
   shared tags still need cardinality, domain, selector and equivalence rules. ⭐ Templating (`tag_templates`,
   which S10 already uses) is what keeps it bounded, ⛔ not shared naming.
4. ⚠ **Recurring costs v1 omitted**, from the legs: export selection, completeness checking, downstream
   invalidation, migration adapters, and ablation tests for the chain itself.
5. ⛔ **The provenance check is unscoped.** ⭐ It is the load-bearing new check and the only one without a
   worked design.

## What is NOT claimed

- ⛔ **Not** that two engines existing guarantees contention. ⭐ It exists only for objects **both engines
  derive from the action**. Everything else is a declared knob, and the register is its control.
- ⛔ **Not** that a common-mode error is caught. Both engines share premises and the action form; ⚠ that is
  the main shared channel and no layer here closes it. ⭐ Independent re-derivation stays mandatory.
- ⛔ **Not** that the migration is free. ⭐ Nothing is deleted until the replacement runs on one step and
  reproduces what the current checks report.

## Migration order

1. ⭐ Build the **free-knob register** first. ⛔ It deletes nothing, it is the missing original ask, and it
   is the cheapest thing here.
2. ⭐ Fix `R4`'s binding to the squared speed both engines already emit; ⭐ wire `R5` the same way.
   ⇒ **that is the cross-step check, available today, without any of the above.**
3. ⭐ Only then pilot direct import on **one** pair of steps, with nothing deleted.
