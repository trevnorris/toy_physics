# S14a part 1 build, review round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_1_build_review_r0_lookups.md`, beside this file
(generator `_scratch/s14a/gen/build_r0_lookups.sh`). "K-n" cites a section of that file.

## The round

- **Artifact.** The Codex build of S14a part 1, to the decision list at `e9eed614`. It has five files, preserved
  unchanged in the commit that adds this file:
  - `scripts/S14a_1_projection.py`;
  - `mathematica/S14a_1_projection.wl`;
  - the two transcripts under `scripts/out/` and `mathematica/out/`;
  - the run record `directives/_measurements/S14a_1_build_runs.md`.

  All five matched the baseline after the round. The working tree was unchanged apart from the lookups file itself
  (K1). The builder exited 0, at about 99k tokens (K1).
- **Prompt and legs.** `_scratch/s14a/build_review/s14a1_build_review_prompt_r0.md`, sha256 `0fbefa30…` (K1). The
  artifact is Codex-written, so the legs were a fresh Claude agent (Opus) and Grok (`grok-4.7`, effort high), each
  given the identical prompt. Both ablate the Wolfram engine, so they ran one after the other: Claude first, then
  Grok. Claude wrote `/var/projects/_review_scratch/s14a1-build-r0-claude/REVIEW.md`, and Grok wrote
  `/var/projects/_review_scratch/s14a1-build-r0-grok/REVIEW.md`.
- **A path in the prompt.** The guard accepts log directories only under `_scratch/s11c`. The prompt's
  `_scratch/s14a/build-review-r0-guard/` was a plain directory, so Claude linked each run into `_scratch/s11c`.
  Before Grok ran, that directory was moved aside to `_scratch/s14a/build-review-r0-guard-claude/` (a rename), and the
  prompt's path was made a link to a fresh `_scratch/s11c/s14a1-build-review-r0-grok/`. So Grok got the same prompt
  and saw none of Claude's run names.
- **Verdicts.** Claude: NEEDS REVISION, 3 findings. Grok: NEEDS REVISION, 2 findings (K2). Each leg derived the
  objects with its own script before opening either engine. Each tested the engines' printed objects on its own exact
  solution of the supplied balances.
- **What both legs found sound.**
  - The field class: live regional fields with jumps, live faces, three in-plane coordinates, no steady state, no
    `J_χ = 0`.
  - The spherical specialization: native fields on radius, `w` and time; faces on radius and time; labelled
    wherever printed; built separately from item 1.
  - The outer flux, printed at a free outer radius with no limit.
  - The projected velocity and its volume flux, which carry the nonzero-density condition.
  - Dimensions, computed from a stated assignment.
  - The supplied inputs. The face factor is printed, never used to build face terms.
  - Emission: no verdict tokens, operands before residuals, nothing conditional on a value.
  - No reads across engines.
  - Byte-identical guarded reruns of both transcripts.
  - Claude also found the core channel, the outer flux, storage and the SymPy order lost matching its derivation on
    its witness.

## Findings

**B1. Both engines leave the inside region's bulk flux out of the face terms.** (Claude 1, Grok 1.) **CONFIRMED.**
- **The list.** In item 1, each projected balance is "a balance for the time rate of a projected density and the
  in-plane divergence of its projected flux, with the faces' motion and tilt kept". It prints "every term that
  leaves at each face". In item 2, the face channel is the net material crossing the faces "measured relative to the
  faces' own motion".
- **The engines.**
  - SymPy builds the endpoint pieces from the time derivative and the in-plane derivatives only. It then takes
    `leaving = atom - endpoint` and adds it to the projected derivative (`:155–160`, `:172–173`, K3).
  - Wolfram builds `pieces` over `Most[flux[[2]]], Most[vars]`, then takes `faceTerms = Expand[atoms - Rest[pieces]]`
    (`:98–99`, `:105`, `:111`, K3).
  - In neither is the bulk-direction flux's own `w`-integral evaluated, so the inside region's bulk flux appears in no
    face term.
- **The evidence.**
  - Each leg's derivation leaves at each face the material crossing relative to the face's motion. Claude's residuals
    against `n a` are `0`, and so is its spherical balance residual. Grok's shell balance minus the integrated balance
    is `0` (K3).
  - On each leg's exact solution, both engines print a nonzero closed-projected-balance residual and a face channel of
    `0`, where the legs' derivations give a nonzero face channel. The engines' shell-balance residual is nonzero
    (K3).
  - SymPy's own closed integral of the material operand prints `0` on Claude's solution, while its "closed projected
    balance" residual does not (K3).
- **Why it changes what is computed.** The face terms, the face channel, `J`, the shell balance and the face channel
  minus the order lost are all built from these terms.
- **What must be true.** Each printed projected balance is the closed-region `w`-integral of its supplied balance,
  with each face's concentrated part at full weight. The face terms are every term that integral leaves at each face,
  with the faces' motion and tilt kept, including the bulk-direction flux. The face channel, `J`, the shell balance
  and the face channel minus the order lost are built from those terms.

**B2. The Wolfram engine computes no face-concentrated part.** (Claude 2, Grok 2.) **CONFIRMED.**
- **The engine.** The region masks use `UnitStep` (`:61–62`). The face parts are taken as
  `Coefficient[Expand[local], DiracDelta[…]]` of the masks' derivatives (`:85–87`, K4).
- **The evidence.**
  - The Wolfram transcript contains `DiracDelta` zero times. Every "face concentrated parts" line prints `{0, 0}`, in
    the general object and in the specialization (K4). The SymPy transcript prints nonzero concentrated parts at the
    same places (K4).
  - Claude's probe shows that `D[UnitStep[w - lo[t]], t]` is a `Piecewise` with no `DiracDelta`, so the coefficient is
    `0`. `HeavisideTheta` gives a `DiracDelta` (K4).
  - Grok's probe prints the engine formula's coefficients as `0`, and the engine's order face parts as `0` against a
    nonzero value from its derivation (K4).
  - Claude's FORM ablation that drops the face parts leaves the Wolfram transcript byte-identical (K4). The object was
    never computed.
  - The run record says "Requested objects not produced: none" for both engines (K4). For Wolfram, that does not hold.
- **What must be true.** The Wolfram engine computes every face's concentrated part at full weight. It enters the
  material and order projections, the projected order source, the order lost and the face channel. A FORM ablation
  that removes it moves the transcript. The run record states any requested object an engine did not produce.

**B3. The two engines share one construction of the face terms.** (Claude 3.) **CONFIRMED.**
- **The list.** "Both engines construct these independently, with the `.wl` blind." The brief adds that the Wolfram
  script "is not a transcription" of the SymPy one.
- **The engines.**
  - The Wolfram script reads nothing from the SymPy side (count 0, K5).
  - Its structure follows the SymPy script function by function: `boundary_differentiation` and `leibnizPieces`, a
    10-argument `supplied`, a `shell` helper. The face terms are the same "atoms minus endpoint pieces" formula
    (py `:160`, wl `:105`), and the tag vocabulary is shared (K5).
- **Why it changes what may be claimed.** B1 is the same defect in both engines. On each leg's solution, both print
  the same wrong face channel and shell residual (K3). A comparator would read that as agreement. Cross-engine
  agreement on the face terms would carry no independent weight.
- **What must be true.** The Wolfram engine reaches its projected balances and face terms by its own construction
  from the supplied balances and the U1 encoding, not by the SymPy engine's route. An error in one engine's projection
  then shows as a disagreement.

## The leg disagreement

None on B1 or B2: both legs found them, each with its own witness. Grok did not raise B3. K5 bears it out.

## Adjudication

All three findings change what is computed (B1, B2) or what may be claimed from agreement (B3). The build goes to a
repair round, with the authors split by engine:
- **SymPy (B1):** the builder's session is resumed.
- **Wolfram (B1, B2, B3):** a fresh Codex author rewrites the engine from the decision list. B3 comes from one author
  constructing both engines. An author who already holds the SymPy construction would bring it to the rewrite.

The two run one after the other, because both update the run record. The repaired build then gets its own two-leg
review, until clear.
