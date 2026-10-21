# S14a part 2 decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S14a_2_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s14a2_dl_p1_lookups.sh`). "K-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** `directives/S14a_2_decision_list.md`, orchestrator-written, sha256 `08da185f…`, with
  `S14a_2_setup_lookups.md`. Both are preserved unchanged in the commit that adds this file. Both matched the
  baseline after the pass, and the working tree was unchanged apart from the lookups file itself (K1).
- **Role.** A pre-builder decision list (G2: one Codex + Grok pass, fold once). It also names the objects the build
  computes, so that content is reviewed until clear (spec row).
- **Prompt and legs.** `_scratch/s14a2/dl_review_prompt_p1.md`, sha256 `2f80f58b…` (K1). Both legs received it
  identically, in parallel. Codex (xhigh, workspace-write in its scratch) wrote
  `/var/projects/_review_scratch/s14a2-dl-p1-codex/REVIEW.md`, and Grok (`grok-4.7`, effort high) wrote
  `/var/projects/_review_scratch/s14a2-dl-p1-grok/REVIEW.md`. Both exited 0 (K1).
- **Verdicts.** Codex: FINDINGS, 3. Grok: FINDINGS, 2 (K2). Each finding below carries its leg's script and filed
  stdout.
- **What Grok checked and found sound:** the displayed balances against stage 006 and the plan, the face factor,
  item 2's material budget as part 1's objects, the readings not chosen kept out of the computed items, the OPEN
  owners, and the labelled specialization.

## Findings

**Q1. "Constant rates" is posed pointwise, on the source densities.** (Codex 1, Grok 1.) **CONFIRMED.**
- **The list.** "The drain and re-ordering source densities, `n Γ_drain` and `n Γ_return`, are independent of
  time. Each stays a live function of position" (l.71–72, K3).
- **The decision it makes operational.** "The throat is at rest, and the drain and re-ordering run at constant
  rates. The pattern stops changing, except that the slab's steady spreading stays a live term" (l.59–61, K3). It
  puts no condition on `n Γ` at each point.
- **The evidence.**
  - Codex: a spreading slab whose return and drain totals are constant (time derivatives `(0, 0)`) while the
    return's first `w` moment grows, `5/16` per unit time in units of `4π`. The support moves, so the source
    densities at a point do change. Conversely, source densities that are independent of time at each point
    (`(1, 0)`) give totals that grow at `7/12` (K4).
  - Grok: homologous spreading with a rate constant along the spreading coordinate satisfies both balances, yet
    `∂t(n Γ_return) = -2*n0*s/h0` at a sample point. The `ω → 0` limit of its modulated family has the same value,
    so the operational premise excludes it. Grok also shows that `∂t Γ = 0` and `∂t(n Γ) = 0` are different
    conditions (K4).
- **Effect.** The build would take a limit other than the one the decision names, adding a pointwise condition
  and still not making the totals constant.
- **What must be true.** "Constant rates" is stated for the rate objects the decision names, with their integration
  domain and measure, and with the frame in which they are constant. No condition is added that the decision does
  not contain.

**Q2. The limit is not defined, so "what the limit drops" is not one object.** (Codex 2, Grok 1.) **CONFIRMED.**
- **The list.** "What steady throughput requires of `n`, `u`, `χ_B`, `J_χ` and the faces is derived by the build
  from the time-dependent balances" (l.76–77). Item 1: "take the limit as defined above", and print "the
  time-dependent object, the limit object and their difference" (l.82–88, K3).
- **The source.** The plan: a static limit "must be **derived** as the `ω→0` limit of the dynamical law, ⛔ never
  posed as an independent frozen solve" (l.599–601, K7).
- **The evidence.**
  - Codex: a solution with source densities independent of time and constant shell totals, whose storage
    `7*cos(t)/12` and core channel `cos(t)/12 + 1` oscillate. Every balance residual is `0` (K5). The balances
    admit it, and the list's definition does not exclude it.
  - Codex: one slow family gives storage `0` at fixed `t` as `ε → 0`, but `7*cos(τ)/12` (divided by `ε`) at
    `t = τ/ε` (K5). The limit depends on the family and its time scaling, and the list gives neither.
  - Grok: within the operational premise, family A has storage `0` and family S has storage `n0*s`. Family E keeps
    a secular change in `χ`, with both source-density derivatives `0` (K5). There is no single persisting part.
- **Effect.** Items 2 and 3 are taken in item 1's limit, so they inherit the gap. The plan's controlled `ω → 0`
  obligation could not be claimed.
- **What must be true.** The list defines the limit itself: the operation that takes any time-dependent solution to
  its limit, and its order relative to the projection and to any far-radius limit. What the limit drops is then one
  object for each printed object. The limit is reached from the time-dependent law, never posed, and leaves S12's
  forms open.

**Q3. The shell totals and the order budget do not name one assembly.** (Codex 3, Grok 2.) **CONFIRMED.**
- **The list.** "the re-ordering total over the shell, from `Γ_return`, and the drain total, from `Γ_drain`"; the
  order budget's "order crossing the faces" and "in-plane order flux through each sphere" (l.93–101, K3).
- **The sources.**
  - Part 1 prints two order-leaving objects at each face, before and after the full concentrated part, and its
    projected order balance uses the second. Its projected order source includes the concentrated parts `C₋`, `C₊`
    (record l.88–91, K6).
  - The S12 decision list: the parts of `n Γ_B` concentrated on a face are part of the source (l.97–98, K6).
  - Part 1's orientations: the outer flux is outward and the core channel inward (record l.126 and l.123, K6).
    Its order lost carries a minus sign and includes the face parts (part 1 decision list l.85–86, K6).
- **The evidence.**
  - Codex: for the same fields, totals over the regular part are `(0, 0)` and totals over the full source are
    `(2, 1)`, in units of `4π`. Adding the separately printed face source to full totals leaves `-1` (K6).
  - Grok: the before-face terms with the full source leave `-2` on a sample. Taking both sphere fluxes outward leaves
    `2*a**2*H_r(a)`, `32` on a sample. The two sign choices for the totals differ by
    `2*int_n_Gamma_drain - 2*int_n_Gamma_return` (K6).
- **What must be true.** Each total names its integrand and measure, and whether face-concentrated parts are in it.
  Each concentrated part is counted exactly once. The face terms in the order budget are named as one of part 1's
  two objects. Each sphere's flux carries its orientation. The signs follow from `Γ_B = Γ_return − Γ_drain`.

## Agreement between the legs

Both legs found Q1 and Q2; Codex reported them as two findings, Grok as one. Both found Q3, each with its own
witness. Neither leg reported a finding the other contradicts.

## Adjudication

Q1 to Q3 each change what is computed. All three are folded once into the list (G2). The fold changes the named
objects, so the folded content is reviewed until clear: a Codex + Grok pass on the folded list runs before any
builder.

**For the orchestrator's own record.** Q1 is a form the orchestrator dictated. "Made operational here" turned the
decision's "constant rates" into a pointwise condition the decision does not contain.
