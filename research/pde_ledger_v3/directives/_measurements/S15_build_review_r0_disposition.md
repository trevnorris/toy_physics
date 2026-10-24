# S15 build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_build_review_r0_lookups.md`, beside this file (generator
`_scratch/s15/gen/s15_build_r0_lookups.sh`). "M-n" cites a section of that file. Every computed result quoted here is
a leg's saved stdout; the legs' scripts were not re-run.

## The round

- **Artifact.** The Codex build of S15 to the decision list cleared at `91248b67`, uncommitted on `a9acfb1d`:
  `scripts/S15_two_drain_region.py`, `mathematica/S15_two_drain_region.wl`, their transcripts (about 24.6 MB and
  9.8 MB) and `directives/_measurements/S15_build_run_record.md`. All five matched the baseline after the round. The
  working tree was unchanged apart from the lookups file (M0).
- **The runs.** The final guarded runs, `_scratch/s11c/s15/sympy_20261024_03` and `wolfram_20261024_02`, exited 0 with
  limits verified and empty stderr (M2). The first build run stopped before any CAS ran, because the brief named a
  guard log path that does not resolve under `_scratch/s11c/`. The orchestrator added the link and resumed the same
  builder session. The contract did not change.
- **Prompt and legs.** `_scratch/s15/build_review/s15_build_review_prompt_r0.md`, sha256 `5ebb9c2b…` (M0). Both legs
  received it identically, in parallel:
  - a fresh Claude agent (Opus): `/var/projects/_review_scratch/s15-build-r0-claude/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-build-r0-grok/REVIEW.md`, exit 0.
- **Verdicts.** Claude: NEEDS REVISION, one finding. Grok: CLEAR (M1).

## Finding

**H1. Each transcript states the region's balance in a sentence that no computation produces.** (Claude 1.)
**CONFIRMED.**
- **The contract** (`S15_decision_list.md:115`, `:126–127`): "**The region's momentum balance.** From the supplied
  balance: the rate of the region's momentum, the integral of `G^i`, beside item 1's parts at both surfaces and the
  integrated force density, each printed apart"; "The script prints objects and never states conclusions."
- **The engines** (M3):
  - SymPy prints `balance_construction: … Applying divergence and moving-domain transport to the supplied local
    balance produces the printed surface parts and undivided volume force. The balance residual is an equation
    operand, not imposed on arbitrary independent live fields.` (`scripts/S15_two_drain_region.py:310`; transcript
    `:2241`).
  - Wolfram prints `"… Integrating the supplied local balance and applying divergence and moving-domain transport
    gives the separately printed surface parts plus undivided force. The residual is printed as an equation operand,
    not imposed on arbitrary independent live fields."` (`mathematica/S15_two_drain_region.wl:176`; transcript
    `:1365`).
  - Each engine prints `region/component_i/balance/residual` as the rate minus the parts and the integrated force
    density (`.py:297`, `.wl:171`). For general `f` it is not zero. Neither engine reduces it with the supplied local
    balance.
- **The evidence** (M3):
  - Claude: on an instance where the supplied balance holds by construction (`f := ∂_tG + ∂_jΠ`), both engines'
    printed objects close: `rate_minus_parts_minus_f=-1.288e-14` (SymPy) and `-9.769962616701378e-15` (Wolfram),
    with the other components at the same scale.
  - With the supplied balance corrupted (W7, the stress's sign in `Π`), the balance opens by `-6.55`, `-10.49` and
    `5.22`.
  - The `balance_construction` line is byte-identical to the working-tree transcript's in W1, W2 and W7 (Wolfram)
    and in P1 and P2 (SymPy). W1 and P1 flip the inner normal, and W2 and P2 drop `X_2`'s motion from `v_S`. Each of
    these opens the balance.
- **The legs differ.** Grok read the residual as an honest equation operand, which it is ("It is not set to zero";
  M3). Grok did not address the sentence's claim that the printed parts are what the supplied balance produces. That
  claim does not move under ablations that make it false, so it is a typed conclusion (E1). The orchestrator
  adjudicates with Claude.
- **Effect.** Every printed object is unaffected: Claude's evaluator matches all 309 evaluated objects of each engine
  to its own derivation (M5). The record could cite the sentence as the engine's evidence that the balance closes,
  and no engine computation supports that.
- **What must be true.** A transcript line says that the printed surface parts and the integrated force density make
  up the rate under the supplied balance only if a computation in that engine produces that relation, with its
  operands and residual printed. Otherwise the transcript states only what it printed, and it labels the balance
  residual as one that no route reduces. No object of items 1 or 2 changes.

## Notes (not findings)

- **`c_comp`.** Each engine prints a free symbol, `c_comp(x,t)` (SymPy, labelled
  `c_comp_not_present_in_supplied_equations`) and `ccomp[x,y,z,t]` (Wolfram, labelled
  `LIVE_PARAMETER_NOT_IN_SUPPLIED_BALANCE`). It enters no object (M4). S9b defines `c_comp² ≡ dp_br/dρ_br`
  (`directives/S9b_SHARED_PHYSICS.md:367–370`). The record states `c_comp` only as `p_br`'s derivative, and it does
  not cite the free symbol as an object.
- **Claude's process.** The Claude leg stopped two of its own guarded jobs, a stuck derivation and an endless loop in
  its own script. Both ended through the guard's normal exit path. Neither touched the artifact.

## What the legs checked and found sound

Both legs derived the contract's objects first, from the supplied balance, with every field live and both centres and
both radii moving (M6).
- Item 1: the four rest pieces and the two transport pieces at each surface. The inner normal points toward `X_2`. The
  Cartesian components carry no turned frame. Each term is in the group of exactly its labels. `Π_o` and the whole
  pressure are printed whole. The expansion, its groups, its remainder and its condition are printed.
- Items 2 to 4: the region's objects apart, the force density undivided, dimensions from stated signatures, the
  domain printed.
- Live arguments throughout. No force law, normalization, mass or throat identification.
- No value-dependent emission or verdict tokens. The Wolfram engine reads nothing from the SymPy side and is a
  different construction.
- Guarded reruns reproduce both transcripts byte for byte (M5).
- Every FORM ablation moved the output: the inner normal flipped, `X_2`'s motion dropped, `X_1` held fixed, two-label
  terms misfiled, the bulk viscosity dropped, a linear pressure law. So did the one-sided corruption of the supplied
  balance.
- Claude also notes that the printed partition residuals are zero by construction and show completeness only. A
  misfiled term leaves them at `0`; the group-by-group comparison catches it.

## Next

Preserve this reviewed baseline, then repair H1 in the builder's session. Round 1 reviews the repair (fresh Claude +
Grok).
