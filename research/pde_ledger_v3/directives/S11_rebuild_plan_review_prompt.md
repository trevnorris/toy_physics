# Independent review — the S11 rebuild plan, before any repair or build starts

## Artifact under review

The block headed **"S11 — THE CURRENT STATE, MEASURED 2026-08-09"** in
`/var/projects/toy_physics/research/pde_ledger_v3/REBUILD_HANDOFF.md`.

⚠ It runs from that heading down to (and excluding) the heading
`"## S11 — where it actually is (⚠ PROVENANCE, superseded by the block above)"`. ⛔ **Everything after that
second heading is explicitly marked stale provenance and is NOT under review** — but ⭐ do read it if you
want to check whether the new block wrongly discards something still true.

It is **orchestrator-written**. It is the plan two engine rebuilds will be driven from. **Nothing has been
repaired or built yet.** You are one of two independent legs; the other is not visible to you.

## ⭐⭐ REQUIRED READING ORDER — ⛔ do not read the plan first

For a document, blindness comes from **reading order**. Reading the plan first anchors you to its framing,
which is the thing under test.

1. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` — **914 lines**, the spec the plan proposes to
   repair. Read it properly; it is the primary source for most of what follows.
2. `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — the closed, reviewed record for the
   **previous** step. ⭐ The plan's two central claims are derived from what this record measured.
3. **Write down, before step 4:** what would *you* fix in the S11 spec before building engines from it, and
   what would you require of the build order? ⭐ Keep that list.
4. **Only now** read the plan block.

## What to check

### 1. ⭐⭐ ARE THE TWO CLAIMED SPEC DEFECTS REAL? — this is the most important question

The plan asserts two defects in a spec that is **already closed after 4 rounds and 8 legs**. If either is
wrong, a repair round will damage a correct spec.

**(a) "S11 has no inertia control."** The plan states that all seven packages vary only the stiffness
functional `W`, that the kinetic term is identical across all of them, and that an exact-token search of
all 914 lines finds no `aniso` / `isotropic inertia` / `s_rho`.
⭐ **Check this yourself against §2 and §7.** Run the searches. Report what you find.
⭐ **Then the harder question: does it MATTER?** The plan's argument is by analogy to S10, where the
transverse count needed two structural premises and the stiffness controls probed one. ⚠ **Does that
analogy hold for what S11 actually computes?** S11's question list is §6. If S11's load-bearing outputs are
not of a kind that an inertia change could move, the defect is real but inert, and the plan overstates it.
⛔ Say so if that is what you find.

**(b) "C16 applies verbatim: §Q8a/§Q8b enumerate strata from minors of `M_r`, never from the stacked
`[M_r; kᵀ]` that governs `nu_T`."**
⭐ Check `:321-340`, `:534`, `:549` and anything else bearing on it. ⭐ Is the stacked matrix genuinely
absent from the stratum construction? ⚠ **And is the consequence stated correctly** — that a locus where
the transverse count *alone* moves cannot be found by construction?

### 2. ⭐⭐ WHAT IS MISSING FROM THE PLAN?

⭐ **Compare your step-3 list against it.** What would you repair or require that the plan does not? This is
the plan's most important failure mode and it is **invisible from the plan alone**.
⚠ In particular: are there **other** defects in this 914-line spec that a rebuild would bake into both
engines? ⛔ The spec being closed is not evidence it is correct.

### 3. Does the plan leak an answer a builder could converge on?

⭐ The plan deliberately says *"the object to add is a control that breaks inertial isotropy on one axis
while holding `W` fixed"* and *"do not state its expected effect."* ⚠ **Is that line drawn correctly?**
⛔ A prohibition leaks as surely as an assertion. ⚠ Also check the plan's treatment of
`W_XFORM_CURLONLY` = S10's `MAIN` action: it says withhold that from both builders. Is withholding right,
and is it achievable given what the builders will otherwise read?

### 4. Is the measured state table accurate?

⭐ Check every factual claim: file sizes, the `registry_read` import at
`scripts/S11_stray_longitudinal_sympy_audit.py:21`, the absence of any `S10_exports` import, the spec line
count, the WL engine's stated tag count and runtime, the five "open items" and their cited spec lines
(`:549`, `:647`, `:691-704`, `:887-901`). ⛔ **Open them. A path that resolves is not a source.**

### 5. Is the build order right?

⭐ Six gates: spec repair → PY engine → WL engine → run both → comparator → record/card/registers.
⚠ Is any gate missing, in the wrong order, or under-specified? ⭐ Does anything in it contradict
`S9_REWRITE_PLAN.md` or what S9/S10 actually did? ⚠ Does the plan correctly assign **who reviews what**
(orchestrator-written → Codex + Grok; Codex-written → fresh Claude agent + Grok)?

### 6. Does the plan commission anything harmful?

⚠ A repair round that touches a **correct** part of a spec breeds new defects in the material it changes.
⛔ Is any instruction in the plan broader than the defect it names? ⭐ Name what would be damaged.

### 7. Unsatisfiable or fabrication-forcing?

⚠ Is any requirement satisfiable only by inventing something? ⭐ The plan flags one existing case
(`PREMISE_INVENTORY` is exempt from the spec's live-read rule at `:887-901`, because several supplied
premises assert an **absence** and have no CAS object to read from). ⛔ Are there others it missed?

## Physics filter

Report a finding only if it catches a way the **physics could be wrong**, or a way the **plan would cause a
wrong engine to be built**. ⛔ Not style, not formatting, not "a builder might misread this" absent a
concrete reading that produces a wrong build.

## Method

- ⭐⭐ **Quote both sides for every finding**: the plan's text, and the spec/record/source text it fails
  against. A finding without both quotations is not usable.
- ⭐ For any claim about a file, a line, or a count: **state what you opened or ran, and what you found.**
- ⛔ Do **not** edit the plan, the spec, or anything else. Read-only.
- ⛔ Do **not** start the repair or write any engine code.
- ⭐ End with: **which items on your step-3 list the plan handles, and which it misses.**

⚠ **A leg that returns "nothing survives the filter" is weak evidence.** If that is genuinely your
conclusion, say so plainly and state what you checked that could have failed and did not.
