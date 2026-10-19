# S9b oracle check: the flow-only condition against the river/acoustic-flow prior art (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-19. **Author of the check:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names
the object the check computes, so that content is reviewed until clear (spec row). The builder starts only after
that. The check is a Codex-written instrument, reviewed by a fresh Claude agent and Grok until clear (G1, E1), before
its output is used.

**Why it exists.** The user chose "Queue it after record"
(`directives/_measurements/S9b_record_r2_review_disposition.md` at `4760e6e2`, lines 19–21). The accepted S9b record
(`ff3bac63`) names these oracles and applies none. M3: prior art is an oracle, checked against our computed result,
never a premise.

## Deliverable
- A committed check script under `scripts/`, its output transcript under `scripts/out/`, and a short run record under
  `directives/_measurements/`. The builder chooses the names. A name may name an object, never its value or a
  comparison outcome.

## Sources (read-only)
- **S9b's objects.** The accepted engine transcripts at `503736ff`:
  - `scripts/out/S9b_light_bending_sympy_audit.out`;
  - `mathematica/out/S9b_light_bending_mathematica_audit.out`.

  Both are git-annex content: run `datalad get` first. The tags are each engine's `R_FLOW_ONLY_B_*` condition and
  implied-exchange objects: PY's squared velocity, general solution with its integration constant, velocity
  components and `IMPLIED_JN`; and WL's `GeneralSquaredVelocity`, `SquaredVelocity`, constant domains, cases and
  `Exchange`. The record `steps/S9b_light_bending.md` at `ff3bac63`, §3, transcribes them, but the transcripts are
  the source.
- **The oracles.** Those whose equations the scoping survey opened (`directives/O2_steady_brane_balance_scoping.md`,
  §5, at `36f59ff1`):
  - Hamilton and Lisle (2008), "The river model of black holes", AJP 76, 519–532, manuscript gr-qc/0411060, §§I–II,
    equations (2) and (4);
  - Visser (1998), "Acoustic black holes", CQG 15, 1767–1791, manuscript gr-qc/9712010, §§2 and 7–8.

  §5 also gives each oracle's conditions for comparison.
- **The repository's survey** `docs/medium_requirements_and_prior_art.md` (repository root), lines 246 and 263–266.

## Model point
- **Regime.** S9b's exterior domain, `r > b_min`, on the flow-only restriction: `δ ≡ 0`, `ξ_w ≡ 0`, with `V` and
  `ρ_br` live. The restriction is a substitution the build defined, not a premise added here.
- **Freezes.** `δ ≡ 0` and `ξ_w ≡ 0` define the restriction. Each engine retains the zero integration-constant
  branch by a supplied counting, in its own form. PY labels the excluded branch
  `EXCLUDED_BY_SUPPLIED_V_OVER_C0_O_EPSILON_HALF`. The record transcribes WL's as the supplied `V/c0=O(r^(-1/2))`
  counting. Either counting is a supplied input, not a derived result. No other quantity is held fixed.
- **Transfer.** Only through the printed identification map below.

## What the check computes
For each oracle above and each engine, the check prints the following.
1. **The oracle's object.** Hamilton–Lisle's river velocity profile and its mass-function form. Visser's flowing-fluid
   background: its velocity, its continuity and density conditions, and its acoustic line element as the paper
   states them. Each is quoted verbatim, with section and equation number.
2. **The identification map.** The correspondence between S9b's symbols and the oracle's, for each condition §5
   names:
   - velocity, and the frame it is measured in;
   - free-fall versus lab time;
   - units, and `c₀` against the oracle's light or sound speed;
   - spatial dimension;
   - mass calibration: S9b's orbital `GM` against the oracle's mass parameter;
   - for Visser, also continuity, exchange and density, against S9b's implied `j_n` and live `ρ_br`.

   Each entry either cites the source that fixes it or prints `UNMATCHED`.
3. **The residual.** The difference between each engine's printed object and the oracle's object, under each
   identification the sources leave open. It is computed at S9b's retained order. Where the oracle's object is exact
   or covers more than S9b's restriction, the check prints the order used and what is left out. The comparison
   covers the general solution with its integration constant live, as well as the retained branch.

## What must be true of the check
1. Every S9b operand is read from the transcript, not retyped, so the check fails if its transcript is missing. Each
   engine's object is used as that engine printed it. An unevaluated object is evaluated without changing its
   integrand, limits or domain.
2. The oracle's objects are supplied inputs. They are quoted verbatim with their citation, and the check cannot
   falsify them.
3. No identification is chosen by its residual. Every identification the sources leave open is printed with its own
   residual.
4. The check prints operands, map and residuals, then guards. It states no conclusion.
5. The run record names the command, the transcript and its sha256.

## Out of scope
- The register's second-route question for these oracles: why the prior result was rejected in its own time, and
  whether that objection is answered here. It is a literature question and stays open with the register.
- Unruh (1981): §5 records its full text as unverified. Cahill (2003) and Sbitnev: named in the survey (line 246),
  not opened in §5, and the survey flags Cahill's status (263–264). Each stays open, with that reason.
- Any edit to the S9b record or the register. The check's outcome goes to the user first.
