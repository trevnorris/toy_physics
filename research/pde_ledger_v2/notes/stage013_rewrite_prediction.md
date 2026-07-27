# stage013 — PREDICTION, written BEFORE the `.wl` emission and before any comparison

HEAD `ae7e4a7d`. Stage `ledger_stage013_breathing_harmonic_mk_projection`. Basis **`(L,M,T)`** (`.py:176`).

## What the committed `.out` renders today (`:83-90`)

| what | how it is printed | live or hardcoded? |
|---|---|---|
| `L0=(1,0,0)`, `beta=(-1,0,0)`, `muEta=(-1,1,0)`, `Tw=(1,1,-2)`, `rAL=(0,0,0)` | one hardcoded string, `.wl:446` | **hardcoded print** |
| `K_eta = Tw*beta^2 = (-1,1,-2)` | one hardcoded string, `.wl:447` | **hardcoded print** |
| `[M_AB] entries = <\|aa->M, aL->M, LL->M\|>` | `KeyValueMap` over `dim["MDims"]`, `.wl:448` | **LIVE** |
| `[K_AB] entries = <\|aa->M T^-2, …\|>` | `KeyValueMap` over `dim["KDims"]`, `.wl:449` | **LIVE** |
| `[M_AB] shared = {0,1,0}`, `[K_AB] shared = {0,1,-2}`, `[K/M] = {0,0,-2}` | computed | **LIVE** |

## PREDICTIONS

**P1 — emit from LIVE VALUES, never from the print strings.** The two hardcoded `Print` lines are the
catalogued vacuous pattern (task #22): comparing them against the `.py`'s identical literals proves
nothing. **But the underlying symbol dimensions are independently declared in each engine** and feed
real computation, so emitting them *from the live variables* is a genuine check on transcription drift
between the two engines' inputs — the same status as stage012's posted axioms. The distinction is
**print-string vs live-value**, not literal vs derived.
⚠ If the `.wl` holds no live binding for a sourced dim, it is `UNREACHABLE` — do not synthesise one.

**P2 — expect ≈9 live records minimum**: `M_AB` entries ×3 (`aa`, `aL`, `LL`), `K_AB` entries ×3,
`M_AB` shared, `K_AB` shared, `K/M` ratio. Plus up to 6 more (5 sourced + `K_eta`) if live bindings
exist. So 9–15.

**P3 — values, from the committed `.out`:** `M_AB` shared `{0,1,0}` = `M`; `K_AB` shared `{0,1,-2}` =
`M T⁻²`; ratio `{0,0,-2}` = `T⁻²`; every `M_AB` entry `M`; every `K_AB` entry `M T⁻²`.
None should move.

**P4 — PASS tally unchanged; `.out` reproduces byte-identically after `$NNNNN` normalisation.**

**P5 — ⭐ `K_eta` enters the canonical table as `M L⁻¹T⁻²`** (`= Tw·beta²` = `(1,1,-2)+2·(-1,0,0)`).
This is the quantity whose three-stage spread was investigated: **013 line-density `M L⁻¹T⁻²`** vs
016 volume-density `M L⁻³T⁻²` vs 023 reduced scalar `M T⁻²`, all integrating to `M T⁻²` under their own
measure. It matches `parameter_register.md:179`. **Prediction: when 016 and 023 are later converted,
the table flags this group `NEEDS_ADJUDICATION` — and that is the system working correctly, not a
defect.** Nothing should be renamed to make it disappear.

**P6 — the transposition ablation will be weak-ish.** Many stage013 quantities are pure-`M` or
`M T⁻²` with a zero `L` exponent, so an L↔M swap moves them but M↔T less so. Expect fewer detectors
than stage012's 16/18. Report the count; do not assume.

## What would falsify the approach
- A rendered value changes ⇒ not behaviour-preserving.
- An emitted record traces to a hardcoded `Print` string rather than a live binding ⇒ vacuous.
- Zero detectors on any transposition ⇒ stage013 has no working gate (the stage021 condition).
