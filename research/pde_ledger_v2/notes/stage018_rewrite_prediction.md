# stage018 — PREDICTION, written BEFORE the `.wl` emission and before any comparison

HEAD `a6a28d0c`, tree clean. Stage `ledger_stage018_dtn_hankel_fingerprint`.
Basis **`(L,M,T)`** — `.py:163` (`Dim` docstring "Exact exponent vector in (L, M, T) order"),
`.wl:95` (`pairs = {{"L", d[[1]]}, {"M", d[[2]]}, {"T", d[[3]]}}` — a positional→label binding,
not prose). ⭐ **The two engines' axis orders are independently evidenced, and they agree.**

## What the committed `.out` renders today

| `.out` | what | how printed | live or hardcoded? |
|---|---|---|---|
| `:71` | `dimension order = (L,M,T); [a]=L, [c_s]=L*T^-1.` | one string, `.wl:386` | ⛔ **HARDCODED** — the catalogued task-#22 defect named in `DIMENSION_REWRITE.md` §5 |
| `:72` | `computed dimensions = <\|u2 -> T^2, u4 -> T^4, v5 -> T^5\|>` | `dimText[u2Dim]` etc., `.wl:387` | **LIVE** |
| `:77` | `corrupted u2 expression = a^2/(9*cs^3); dimension = L^-1 T^3` | `dimText[corruptedU2Dim]`, `.wl:392` | **LIVE** |

## PREDICTIONS

**P1 — every value is REACHABLE; stage018 should need no D2 relaxation and no waiver.**
All ten dimension objects are **top-level globals** in the `.wl`, not `Module` locals:
`zeroDim`/`dimL`/`dimSpeed`/`dimT2`/`dimT4`/`dimT5` at `:84-89`, and `u2Dim`/`u4Dim`/`v5Dim`/
`corruptedU2Dim` at `:229-233`. So emission is **print-only** — the cheapest case in the corpus, and
the reason 018 is A1. Predict **zero `UNREACHABLE`, zero waivers**, matching stage013.

**P2 — expect 10 records.** 4 genuinely computed by `dimOf` walking an expression
(`u2`, `u4`, `v5`, `corrupted_u2`) + 6 declared literals (`a`, `c_s`, the three expected `T^n`
targets, `zero`).

⚠ **That ratio is worse than stage013's** (which was 10 computed of 15). **6 of 10 records here will
be literal-vs-literal**, so the comparator can only catch a transcription divergence between the two
engines' declarations — a *shared* wrong declaration passes. The physics-vs-`model_map.md` leg
(§4 (f2)) carries proportionally more of the load at 018 than at any converted stage so far.

**P3 — the predicted values** (`(L,M,T)`, from the committed `.out` and `.py:191-196`):

| name | exponents | labelled | provenance |
|---|---|---|---|
| `a` | `{1, 0, 0}` | L | declared literal both engines |
| `c_s` | `{1, 0, -1}` | L T⁻¹ | declared literal both engines |
| `u2` | `{0, 0, 2}` | T² | **computed** — `a²/(9c_s²)` ⇒ `2(1,0,0) − 2(1,0,−1)` |
| `u4` | `{0, 0, 4}` | T⁴ | **computed** — `4a⁴/(81c_s⁴)` |
| `v5` | `{0, 0, 5}` | T⁵ | **computed** — `a⁵/(27c_s⁵)` |
| `corrupted_u2` | `{-1, 0, 3}` | L⁻¹ T³ | **computed** — `a²/(9c_s³)`, the able-to-fail tooth |
| expected `T2`/`T4`/`T5` | `{0,0,2}`/`{0,0,4}`/`{0,0,5}` | T²/T⁴/T⁵ | declared literal both engines |
| `zero` | `{0, 0, 0}` | 1 | declared literal both engines |

None should move.

**P4 — ⭐ EVERY exponent's M slot is ZERO. stage018 has no mass content at all.**
This is a structural fact about the stage (a DtN fingerprint in `a` and `c_s` only), not an omission.
Consequence — **predicted transposition detectors, computed from the table above** (a quantity detects
a swap iff its two swapped exponents differ):

| transposition | detects | predicted count |
|---|---|---|
| **L↔M** | `a`, `c_s`, `corrupted_u2` | **3 of 10** |
| **M↔T** | all but `zero` and `a` | **8 of 10** |
| **L↔T** | all but `zero` | **9 of 10** |

⚠ **L↔M at 3/10 is the weakest detector rate of any converted stage** (004 was 17/20, 012 16/18).
It is still non-zero, so 018 does not hit the stage021 no-gate condition — but record the measured
numbers and compare them against these three predictions. **A measured L↔M count of 0 would mean the
emission is not reading the length-bearing declarations, and is a defect, not a property of the stage.**

**P5 — ⛔ `.wl:386` is a hardcoded literal print and must NOT become a `DIM|` record.**
`[a]=L, [c_s]=L*T^-1` is typed into a string. `DIMENSION_REWRITE.md` §5 names this exact line as one
of two pre-existing instances and says *do not add a third*. The `a` and `c_s` records must read the
live `dimL` / `dimSpeed` (or `dimRules`) bindings via `ToString[InputForm[…]]`, exactly as stage013's
`K_eta` does at `.wl:458`. Do not delete `:386` either — it is a catalogued finding, not this task's
to repair.

**P6 — ⭐ the sidecar makes three printed claims FALSE, and this is the first stage where that bites.**
stage018 asserts **zero file I/O** in three places:
- `.py:4` docstring — "Standalone, print-only, no arguments, no file I/O."
- `.py:658` — "EXPORTS: … **no file artifact is written**."
- `.py:686` — "Engine: … **zero file I/O**."

`emit_dimension_sidecar` writes `scripts/<stem>.dimensions.txt`. After the rewrite all three are
untrue. None is an `expect_*` marker, so **no PASS breaks and nothing detects it** — which is exactly
why it has to be caught here by reading. **The `.py` prose must be corrected to state the sidecar**;
stage013 set the precedent (its docstring now reads "with audit results on stdout and labelled
dimensions in a deterministic sidecar"). The `.wl`'s matching claim stays TRUE — it prints to stdout.

⭐ **This generalises: `016`, `021`, `023`, `027` carry the same claim** (`016:4`, `016:874`;
`021:4`, `021:725`, `021:742`; `023:1016`, `023:1035`; `027:4`, `027:833`). Settle the wording once
here and reuse it. **stage011/012/013 did not have the claim, so this never came up before.**

**P7 — stdout byte-identity WILL break, and that is authorized.** P6's prose fix changes stdout.
`018` is inside the 001–028 range where `scripts/output/*.txt` exists, so gate 2 applies — but **D1
makes it a diagnostic, not a blocking gate**. This is the first live exercise of D1. Re-baseline
`scripts/output/ledger_stage018_*.txt`; do not contort the code to preserve the old bytes.

**P8 — ⛔ the forbidden-token guard is a live landmine for the `.py` rewrite.**
`.py:666` checks `set(globals())` against `forbidden_tokens()` = `{mu_hat0, mtilde0, N0, D0,
P0_phys, build_dimensions}` (`.py:624-632`, assembled by string concatenation to dodge a source
grep — see `feedback-grep-acceptance-dodgeable`). The `.wl` runs the same check over
`Names["Global`*"]` (`.wl:316-323`, `:449`). **The rewrite must not introduce a global named
`build_dimensions`.** The membership test is exact, not substring — so the sibling scaffolding this
stage will be modelled on is *safe as spelled*: stages 011/012/013 name their helper
`build_dimensional_block` (`011:330`, `012:527`, `013:380`), which does not collide. The risk is a
rewrite that *shortens* it to `build_dimensions` while copying — an easy, silent, guard-tripping
edit. The four module imports
(`Dimension`, `DimensionBasis`, `dim_residual`, `emit_dimension_sidecar`) are all safe.
`.py:661` additionally requires live symbolic names ⊆ `{z, a, c_s}`.

**P9 — naming, under D4.** ⛔⛔ **P9 WAS FALSIFIED. The claim below is WRONG and is kept, struck
through, as the record of a wrong pre-registration — do not act on it.**

> ~~`c_s` here is the **same physical quantity** as stage012's `clean_walk.cs_dim` (`L T⁻¹`, sound
> speed): standardising the name creates a real cross-stage comparison, which D4 sanctions.~~

**What is actually true.** They are **different quantities that happen to share `L T⁻¹`**:
- stage018's is **`c_s0` = √(5Kρ0⁴/m)** at the **bulk asymptotic** density. It has its own register
  row, `parameter_register.md:129`, which names *"LIVE units carrier at II-G4a (018)"* as its site and
  states outright *"distinct from the frozen-wall `c_S` (011–017)"*.
- stage012's is **`c_S` = 5Kρ\*⁴/m** at the **wall** density `rho_star`
  (`ledger_stage012_…_sympy_audit.py:539`; its own prose at `:673` — *"c_S^2 = 5*K*rho_star^4/m is R1
  at rho_star"*).

**stage018 denies the identity in its own artifacts**: the stage note at `:142-143` (*"⚠ Distinct
from the frozen-wall Helmholtz speed `c_S` … evaluated at the wall density ρ\* vs the bulk ρ0. `c_S`
is not a live symbol here"*), and a **live assertion** that the frozen-wall speed symbol is not
present (`.wl:452-453`, `.py:662-663`).

⛔ **So merging them is precisely the operation D4 forbids** — making *different* quantities share a
name. And it is undetectable by construction: both really are `L T⁻¹`, so **no dimensional check can
ever catch it**, while `generate_canonical_dimension_table.py:293` groups on the **scope-stripped**
name and would have published the pair as an `AGREE` group — a *manufactured* cross-stage
confirmation, in the very artifact Part VII's dimensional firewall consumes.

⭐ **Two independent parties made this same error.** This note asserted it, and the build reached it
separately (it never read this file — verified: its only references to this path are `git status`
lines). "Same dimension + same word *sound speed*" is an **attractive** wrong merge. That is the case
for the physics leg being mandatory: engine-vs-engine agreement is structurally blind to it.

✅ **The `a` half of P9 held.** `a` is the DtN mouth/pin radius (`register:132`, `CONV`, `= ħ/(m_GNLS
c_s0)`); it shares `L` with stage011 `LengthDim` / stage013 `symbol_dims.L0` / stage004 `L` but is
**not** the same quantity, so it was correctly left unmerged.

**Revised prediction: ZERO new groupings, and one name to fix.**
⚠ Grouping is a review flag, never independent verification (§10).

⛔ **PROCESS FINDING — this file is a contamination channel.** It sits untracked in the working tree,
so every agent can read it, and the physics leg did (it cited `:105-106` as "the stated basis for the
merge"). A pre-registration left in the tree **functions as a supplied answer** — the same failure as
stating an expected result in a directive, arriving through a side file instead. Keep the next
stage's prediction **outside the repo** until its build and reviews are done, then commit it.

## What would falsify the approach
- A rendered value changes ⇒ not behaviour-preserving.
- Any emitted record traces to `.wl:386` rather than a live binding ⇒ vacuous, and a third instance
  of the catalogued defect.
- **Zero L↔M detectors** ⇒ the length declarations are not really in the comparison.
- A `DIM|` count below 10 without a stated, locus-cited `UNREACHABLE` reason ⇒ silent under-coverage,
  the failure the comparator's non-empty floor exists to catch.
