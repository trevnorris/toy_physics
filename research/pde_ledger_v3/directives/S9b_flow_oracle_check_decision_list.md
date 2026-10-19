# S9b oracle check: the flow-only condition against the river/acoustic-flow prior art (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-19. **Author of the check:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). It also names
the object the check computes, so that content is reviewed until clear (spec row). The builder starts only after
that. The check is a Codex-written instrument, reviewed by a fresh Claude agent and Grok until clear (G1, E1), before
its output is used.

It was folded once on 2026-10-19: three findings, all confirmed
(`directives/_measurements/S9b_flow_oracle_check_decision_list_review_disposition.md`). The fold changes what the check
computes, so the revised content gets its own review before the builder starts. That review's second pass found
two defects in comparison 3, repaired on 2026-10-19
(`directives/_measurements/S9b_flow_oracle_check_decision_list_review_p2_disposition.md`).

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

  Both are git-annex content: run `datalad get` first. The record `steps/S9b_light_bending.md` at `ff3bac63`
  transcribes the objects below (§2 lines 68–73, §3), but the transcripts are the source.
  - **Flow-only velocity and exchange.** PY: `PY_S9B_R_FLOW_ONLY_B_IMPLIED_JN_*`, its `squared`, `general_solution`
    with `s9b_velocity_integration_constant`, its signed velocity components and their exchange equalities, with
    their domains. WL: `WL_S9B_R_FLOW_ONLY_B_IMPLIED_JN_*`, its `GeneralSquaredVelocity`, `SquaredVelocity`,
    constant domains, and `Cases/*/ImpliedExchange` with each case's attached domains and boundary qualifications.
    `WL_LOCAL_S9B_MASS_BALANCE` is the unsubstituted balance, not the flow-only operand.
  - **The supplied optical input.** The leading-eikonal dispersion `(ω − Vⁱkᵢ)² = cγ² gⁱʲkᵢkⱼ` (record lines
    68–73). It is a supplied input of S9b, not an S9b result. Each engine prints it in its own symbols:
    - PY: `PY_LOCAL_S9B_VECTOR_DISPERSION`, a polynomial in `s9b_chi`, `s9b_A`, `s9b_U` and the wave-vector
      components;
    - WL: the `Dispersion` field of `WL_LOCAL_S9B_SUPPLIED_INPUTS`, in `cSquared`, `aMetric` and `vLocal`.

    That tag's `Metric` field is the spatial metric, not the dispersion. `WL_LOCAL_S9B_FERMAT_OBJECT` is constructed
    by the engine, not supplied. Where the transcript does not bind a dispersion symbol to `V`, `cγ` or `g`, that
    binding is an open identification.
- **The oracles.** Those whose equations the scoping survey opened (`directives/O2_steady_brane_balance_scoping.md`,
  §5, at `36f59ff1`):
  - Hamilton and Lisle (2008), "The river model of black holes", AJP 76, 519–532, manuscript gr-qc/0411060, §§I–II:
    the river line element (1), the velocity profile (2) and its mass-function form (4);
  - Visser (1998), "Acoustic black holes", CQG 15, 1767–1791, manuscript gr-qc/9712010:
    - the acoustic metric (4);
    - §2's continuity (5) and Euler (6) equations;
    - §7's background and its line element (53);
    - §8's velocity (54) and line element (55);
    - §10's eikonal equation (86) and ray form (87). (85) there is the metric with its overall factor removed.

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
Three comparisons, for each oracle that states the object and for each engine. Each prints, in turn:
- the oracle's operand, quoted verbatim with section and equation number;
- S9b's operand;
- the identification-map entries it uses;
- its residual, or `UNMATCHED` where one side has no operand.

1. **Flow profile.**
   - *S9b operand:* the flow-only squared velocity and each signed component. It is taken with the general solution's
     integration constant live, and again on the retained branch.
   - *Oracle operands:* Hamilton–Lisle's river velocity profile, (2), and its mass-function form, (4). Visser's §7 and
     §8 background velocities (§7's background and (54)), each with the conditions the paper attaches to it.
   - *Residual:* the difference of the profiles, under each identification the sources leave open.
2. **Continuity and exchange.**
   - *S9b operand:* each engine's printed implied exchange, with `ρ_br` live.
   - *Oracle operands:* Visser's background continuity equation (5) and the density conditions of its §7 and §8
     backgrounds, quoted as stated, including §7's statement about which profiles continuity admits. Hamilton–Lisle
     states no continuity operand, so it is `UNMATCHED`.
   - *Residual:* S9b's exchange less the oracle's, as a function of live `ρ_br`, with the set of `ρ_br` on which it
     vanishes. No oracle density law is substituted into S9b's operand.
3. **Optical characteristics.**
   - *S9b operand:* each engine's printed dispersion (Sources, above), on the flow-only restriction.
     `WL_LOCAL_S9B_RESTRICTION_FLOW_ONLY` sets only the two profiles to zero and does not rewrite the dispersion. This
     comparison checks a supplied input against the oracle, not an S9b result, and its output says so.
   - *Oracle operands:* the null cones of Hamilton–Lisle's (1) and of Visser's (4), (53) and (55), each with its
     overall factor printed. Visser's eikonal form (86) and ray form (87) are quoted with them.
   - *Residual:* the difference of the characteristic cones at leading eikonal order, under each open
     identification of the dispersion's symbols, with every overall conformal factor printed explicitly. S9b prints
     no line element, so a full line-element residual is `UNMATCHED`.

**The identification map.** One entry for each condition §5 names, each citing the source that fixes it or printing
`UNMATCHED`:
- velocity, and the frame it is measured in;
- free-fall versus lab time;
- units, and `c₀` and `cγ` against the oracle's light or sound speed;
- spatial dimension;
- mass calibration: S9b's orbital `GM` against the oracle's mass parameter;
- for Visser, also:
  - background continuity and force accounting (Euler's equation);
  - barotropic, inviscid and irrotational flow;
  - the equation of state and density conditions of each background;
  - scalar acoustics against S9b's transverse optical branch;
  - no shear-bearing sheet and no off-sheet exchange.

A profile residual is not read as a force balance or as the acoustic comparison.

**Order.** Every residual is at S9b's retained order. Where an oracle's object is exact, or covers more than S9b's
restriction, the check prints the order used and what is left out.

**What the output can support.** Comparison 1 can support at most a conditional correspondence of profiles. No
comparison establishes equivalence of a wave equation, a force law or a material law.

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
