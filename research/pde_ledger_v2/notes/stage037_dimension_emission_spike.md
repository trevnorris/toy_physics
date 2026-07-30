# stage037 — dimension-emission feasibility spike (2026-07-27)

**Verdict: `ROUTE_EXISTS` for stage037.** Run deliberately out of order as a **measurement, not a
conversion**. Nothing tracked was modified; stage037 remains unconverted.

⚠ **This note exists because the spike's working artifacts lived in gitignored
`_scratch/spike037/` and `_scratch/OUT_spike037.txt` — ⛔ both since deleted and not retained, so this
note is now the whole record.** Without it, a fresh clone or a
context-compacted reader could not audit the word "prototyped". It records the **route, the
measured numbers, and the honest limits** — but see *Reproducing this*: it is **not** a
clone-runnable reproduction, and saying otherwise would repeat the overstatement this note was
written to correct.

## Why it was run

`notes/rewrite_reference_table.md` §5.6 recorded stage037 as one of *"the five stages where step (a) is
genuinely impossible"*, and `manifests/DIMENSION_REWRITE.md` §3b noted that **D2 reopened** it — the
print-only rule that blocked it was lifted, but **D3 (engine independence) was not relaxed**.

So the real question was never cost. It was: **can the `.wl` derive exponent vectors from its own
symbolic content, or would every available route reduce to re-transcribing constants that originate in
the `.py`?** The second case is a check that **cannot fail** — worse than no check — and it would have
applied to a whole class of stages.

Seven remaining stages are not more-of-the-same — **035, 036, 037, 044, 008, 038, 042** — and in the
prior remaining order they sat at positions **14–16 and 22–25**, i.e. an infeasibility would have
surfaced only after most of the cheap stages were spent.

## The route

stage037's `.wl` holds **no dimension exponent-vector objects at all** — its firewall is symbolic
monomial substitution over `massScale` / `lengthScale` / `timeScale`. (Its only two integer lists,
`{3,1}` and `{-1,1}`, are unrelated to units.) The route does not look for tuples; it *manufactures*
them from live expressions:

```wolfram
dimensionAxes = <| "M" -> massScale, "L" -> lengthScale, "T" -> timeScale |>;

dimensionScale[expression_] := Factor[
  Together[(expression /. unitRules[False])/expression]
];

dimensionExponents[expression_] :=
  Exponent[dimensionScale[expression], #] & /@ Values[dimensionAxes];
```

Substitute the scaled units into the **original** expression, divide by the original, cancel to a pure
scale monomial, then read the exponents off it. The `axes=` label is built from
`StringRiffle[Keys[dimensionAxes], ","]` — **derived from the live axis map, not typed** (see below).

## What was MEASURED

| | |
|---|---|
| quantities emitted | **21 of 21**, zero waivers |
| names | `q_T, mu_R, rho_br, A_E, c_gamma_squared, c_E, D_V, A_V, kernel_B00, U_B, F_B0, F_B1, F_B2, F_Br, U_A2, r_BA, delta_BA, r_cone, Delta_U, s_1, s_2` |
| prototype size | **+33 `.wl` lines** incl. a perturbation hook; ~27 without |
| runtime | **~5 s**, exit 0 |
| stage PASS tally | **16/16 preserved** |
| `.out` re-run | **required** — the committed `.out` currently goes straight from `PASS UNITS_RESTORED` to the next section |

**Able-to-fail, demonstrated — ⚠ AGAINST THE THEN-CURRENT COMPARATOR.** A seeded `.py`-only error
was caught by the **real** comparator (loaded via `importlib`, not a reimplementation):

```
COVERAGE|stage=stage037|py=21|wl=21|compared=21|py_only=0|wl_only=0
MISMATCH U_B: py={L=3, T=-2, M=1}; wl={M=1, L=2, T=-2}
RESULT|stage=stage037|status=FAIL|mismatches=1
```

⛔ **Do not re-run that probe against the current comparator and expect this output.** It predates
`35710cee`: the `COVERAGE|` spelling is pre-rename (now `ARTIFACT_NAME_SET|`), and — more
importantly — the probe hand-built a sidecar, which **sidecar freshness now rejects before any
mismatch comparison happens**. The result stands as a record of what was measured on the day; it is
**not** a reproduction against today's tooling. Reproducing it now requires a hash-valid sidecar.

**Fabrication guard fired correctly** — perturbing the emitted expression moved the record:
`U_B -> R U_B` gave `{1, 3, -2}` against a clean `{1, 2, -2}`.

## Independence (D3) — three separate signs

1. **Different algorithm.** Wolfram rescales-and-cancels; the `.py` walks the expression tree with an
   atom-dimension dictionary (`ledger_stage037_…_sympy_audit.py:607`, `:635`).
2. **Different axis order, chosen independently.** Wolfram `M,L,T`; Python `(L,T,M)`
   (`…stage037…py:604`). The comparator pairs by **label**, not position — so this is a live
   demonstration that axis-labelled pairing does real work.
3. **A seeded one-sided error is caught** (above).

## ⚠ The honest limit

**8 of the 21 records are effectively primitive declarations** — `q_T`, `mu_R`, `rho_br`, `A_E`,
`c_gamma_squared`, `c_E`, `s_1`, `s_2`. Comparing those detects **transcription divergence between the
engines, not wrong upstream physics**: a value that is faithfully and independently encoded by both
engines but wrong in the shared source report passes. The other 13 are algebraically derived from
those declarations and from live expressions.

⇒ **Independent corroboration that the physics-vs-model leg (§4-c) is mandatory**, reached from a
different direction than the stage018 naming defect that first forced it.

## Cost — measured vs estimated, kept apart

- **Measured (037 only):** ~27–33 `.wl` lines, ~5 s, 16/16 PASS preserved, `.out` re-run required.
- **Estimated:** **0.5–1 engineer-day per stage including review.** Review effort was **not** measured
  for any stage. Treat as a planning figure.
- **Review classification (a judgement, and a firm one):** this is **new mathematics, not print-only**.
  A reviewer must check scale-factor cancellation, Laurent exponent extraction, axis labelling, and the
  nonzero/homogeneity assumptions.

## Transfer to 035 and 036 — ANALYSED, NOT PROTOTYPED

⛔ **Neither was executed. Do not cite these as measured.**

- **035 — looks directly transferable and cheaper.** `unitScalingObject[]` (`…stage035….wl:411`)
  already computes a nine-key association of live scaling monomials ≈ **20 component quantities**.
  Extract exponents from `restored` rather than `expectedScaling` and it needs no new dimension
  evaluator. Obstacle: flattening array-valued keys and aligning semantic names — rendering work, not
  a mathematical barrier.
- **036 — looks transferable but dearer.** `dimensionResiduals[]` discards identities and returns only
  zero residuals, so D2 restructuring must build a named association over **13 logical groups / 29
  component quantities**, flatten it, and get the new data flow reviewed. ⛔ The old *"printing them
  means re-transcribing constants … invisible by construction"* verdict is **avoidable**: 036 has
  native primitive scaling rules (`…stage036….wl:232`) and live kernel/interaction/force expressions
  (`:158`) to derive from.
- **042 and 044 are NOT covered by this spike.** 044 is frozen pending 044-v2; 042 has a separate
  problem (stiffness basis, and a disputed guard-invocation count — `manifests/DIMENSION_REWRITE.md`
  §9 says the guard runs **once**, the reference table says twice; **re-measure before relying on
  either**).

## ⭐ Bonus finding — a better emission pattern than the one in use

The spike derives `axes=` from `Keys[dimensionAxes]`. **By inspection, all five converted `.wl` files
hardcode their axes strings** (`004:220`, `011:407`, `012:568`, `013:447`, `018:387-393`). The stage018
adversarial leg demonstrated the stale-label risk **for stage018 only** — it did not audit the other
four. If a `.wl`'s internal axis order ever changed, a typed label would not follow.

⇒ **Adopt the derived-label form from stage016 onward.** The five already converted are bounded,
recorded debt — not worth retrofitting mid-flight.

## Reproducing this — ⚠ NOT YET CLONE-RUNNABLE, and that is a real gap

⛔ **This note preserves the measured result and the route; it is NOT a complete reproduction.**
The three helper definitions above are the mathematical core, but `emitDimensionRecords[]`, its call
site, and the seeded-error probe existed **only** in gitignored
`_scratch/spike037/` (`ledger_stage037_route_b_boost_structural_relation_spike.wl`,
`compare_transcription_probe.py`, `parse_wolfram_output.py`, `survey_integer_lists.wl`) — ⛔ **that tree has
since been deleted and is not retained**, so those four files no longer exist anywhere. The rest of
the prototype was a copy of the committed stage037 `.wl`.

⇒ **Before citing this spike as reproducible**, track the complete emission block and an
updated, hash-valid probe, and re-run both against the current comparator. Until then this note
supports the claim *"a route was demonstrated on 2026-07-27"* and **not** *"anyone can rerun it
today"*.
