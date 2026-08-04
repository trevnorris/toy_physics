# BUILD — the registry row generator (cross-engine gated)

⛔⛔ **Build ONE tool. ⛔ Do not modify either audit script. ⛔ Do not hand-edit the registry.**

## Why this exists

Registry rows are currently transcribed by hand from an audit script's printed tags. ⚠ That is a manual
step in a project where transcription is a known failure mode, and it means a row can enter the ledger
carrying **one engine's** value.

⭐⭐ **The rule this tool enforces: a value enters the registry only if BOTH independently-written engines
report it and they AGREE.** Disagreement produces a **blocked row with both values**, ⛔ never a silent
choice.

⚠ This would have caught a real defect: one engine's face-response `μ`-coefficient differed from the
other's by a factor of `ρ_br⁰/ρ_m`. Under this tool that row is refused, not entered.

## Deliverable

`/var/projects/toy_physics/research/pde_ledger_v3/reduction/generate_rows.py`

Run as a standalone script. ⚠ Keep runtime well under 10 minutes; it is pure text processing.

## Inputs

- Mathematica tag output: `/home/trevnorris/s11bB_build/wl_run5.txt` (tags prefixed `WL_`)
- SymPy tag output: `/home/trevnorris/s11bB_build/py_run3.txt` (tags prefixed `S11BB_`)
- The schema the rows must satisfy: `research/pde_ledger_v3/reduction/registry_schema.yaml`
- The existing registry, **read-only, for shape and to avoid duplicate ids**:
  `research/pde_ledger_v3/reduction/quantities.yaml`, `relations.yaml`

⛔ **The tool must take these as arguments or module-level constants, ⛔ not hard-code one step's paths as
its only mode.** ⚠ It will be reused for later steps.

## What it must do

**1. Map the two tag sets to each other via an EXPLICIT correspondence table** that lives in the script.
⛔ **No fuzzy matching, no normalisation heuristics, no "strip the prefix and hope".** ⚠ The two engines use
different symbol names (`la0`/`lv0`/`lx0` vs `Lambda_A0`/`Lambda_V0`/`Lambda_X0`, `rhom` vs `rho_m`,
`w0` vs `W0`, and others you must read out of the two outputs).
⭐ **Any tag the table does not cover is REPORTED AS UNMAPPED, not guessed.**

**2. Compare, per mapped quantity:**
- the **dimension exponent triple** — this is directly comparable and is the primary target;
- the **route kind** (`independent` vs `definitional`).

**3. Emit two artifacts:**
- **`reduction/_generated/quantities_S11bB.yaml`** — proposed rows, schema-conforming, for every quantity
  where **both engines agree**. Each row must carry a `source_locus` with the **real** path and line range
  in the script that produced it, and must record honestly whether the stage uses a shared dimensions
  module.
- **`reduction/_generated/BLOCKED.md`** — every quantity that is **not** emitted, with the reason:
  `DISAGREEMENT` (give both values), `WL_ONLY`, `PY_ONLY`, or `UNMAPPED`.

⛔ **Do NOT write into `quantities.yaml` or `relations.yaml`.** Insertion into the live registry is a
separate reviewed step. This tool only proposes.

**4. Relations.** ⚠ Relation statements are prose and algebra in different notations, so ⛔ **do not attempt
to auto-compare them.** ⭐ Emit them to **`reduction/_generated/RELATION_CANDIDATES.md`** as candidates
with both engines' text side by side, clearly marked as **requiring review before insertion**. ⛔ Do not
generate relation YAML rows.

**5. Exit non-zero if anything is blocked**, so the condition is visible rather than buried in a report.

## ⛔⛔ CONSTRAINTS

1. ⛔ **The audit scripts must not import this tool, and this tool must not be imported by them.** ⚠ An
   audit script that can reach the registry breaks the blindness this step was built on.
2. ⛔ **No heuristic that resolves a disagreement.** ⛔ No "prefer the more detailed value", no "prefer the
   one with more digits", no tie-breaks. ⚠ In the one real disagreement so far, the **more** detailed
   engine was the wrong one.
3. ⛔ **Do not invent a value, a dimension, a qid, or a provenance line.** If a field cannot be filled from
   the inputs, block the row and say which field was missing.
4. ⛔ **Do not modify the schema** to make rows fit. If a row cannot satisfy the schema, block it and report
   which constraint it violates.
5. ⭐ Write the tool so a later step can run it by pointing at two different output files.

## Output

The tool, plus a run of it on the two inputs, plus a report **under 25 lines**:
- how many quantities were emitted, and how many blocked by each reason;
- the full list of blocked ones with their reasons;
- confirmation that neither audit script was modified and that the live registry files were not written.

⛔ Then stop and exit.
