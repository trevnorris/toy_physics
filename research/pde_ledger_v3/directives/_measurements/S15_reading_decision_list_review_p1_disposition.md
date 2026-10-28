# S15 reading decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_reading_decision_list_review_p1_lookups.md`, beside this
file (generator `_scratch/s15/gen/s15_reading_p1_lookups.sh`). "M-n" cites a section of that file. Every computed
result quoted here is a leg's saved stdout; the legs' scripts were not re-run.

## The pass

- **Artifact.** `directives/S15_reading_decision_list.md`, revision 1, uncommitted on `5a019394`, with its setup
  lookups `directives/_measurements/S15_reading_setup_lookups.md`. Both matched the baseline after the pass (M0).
- **Gate.** G2: a pre-builder decision list written by the orchestrator, so its legs are Codex + Grok, one pass,
  folded once. The list also names the objects the build computes, so that content is reviewed until clear.
- **Prompt and legs.** `_scratch/s15/reading_dl/reading_dl_review_prompt_p1.md`, sha256 `a005df3e…` (M0). Both legs
  received it identically, in parallel:
  - Codex (`gpt-6.1-sol`, effort xhigh): `/var/projects/_review_scratch/s15-readdl-p1-codex/REVIEW.md`;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-readdl-p1-grok/REVIEW.md`.
- **Verdicts.** Both FINDINGS: Codex two, Grok three, all HIGH (M1). Two of them describe the same defect, so there
  are four distinct findings.

## Findings

**F1 (Codex 1): the family conditions are for the geometric part, but the list reads them as conditions on the whole
departure. Confirmed.**
- *The list:* item 4 asks for "the conditions under which that object vanishes", where that object is item 3's
  `E_Q^geo` on the family. The record section reads "the printed conditions under which its departure vanishes" and
  files "gravity's flow requirement for S21, stated as those printed conditions". The intended claim says "the
  conditions under which each closed-form departure vanishes" (M2).
- *The definition:* `E_Q = E_Q^geo + E_Q^resp` (list `:42`). S15's list asks for each dependence "with its
  conditions on every live quantity the piece contains" (`S15_decision_list.md:148`, M2).
- *The evidence:* Codex took the region-momentum group `drain1+drain2` with `b1 = b2 = 0` (M2).
  - With `a1 = −2` and the amplitudes changing in proportion, the parts are `(Egeo/Q, Eresp/Q, E/Q) = (0, 2, 2)`.
    The geometric part vanishes and the departure does not.
  - With `a1 = 1` and other changes, the parts are `(3, −3, 0)`. The departure vanishes and the geometric part does
    not.
  - The all-direction checks repeat both cases.
- *Effect:* conditions for `E_Q^geo = 0` are neither necessary nor sufficient for `E_Q = 0` while the changes are
  live. Read as the list's record section reads them, they would file a flow requirement the objects do not support.
- *What must be true:* nowhere does the list present item 4's conditions as conditions for `E_Q = 0`. Wherever the
  record or the claim uses them, they are conditions for `E_Q^geo = 0`. That is the departure when no quantity
  changes along the family, and `E_Q^resp` is printed beside it as what the changes add.

**F2 (Grok 3): the split into the two parts depends on how a term is written. Confirmed.**
- *The list:* "`E_Q^resp` holds every term that contains a change `F′`; `E_Q^geo` holds the rest." (M3).
- *The evidence:* for `Q = W/d²` written as an integral, `𝒟Q` is one integral whose integrand holds both kinds of
  piece. Read by top-level term, `E_geo = 2W/d²`. Read after setting the changes to zero, `E_geo = 0`. Both readings
  sum to `E_Q` (M3). Codex's stdout also has a re-description of one quantity that moves `−A·K·d` from one part to
  the other while `E` stays unchanged (M3).
- *Effect:* item 2 would let two engines print different parts for the same object, and let the record read a
  vanishing or a nonvanishing geometric departure for an inverse-square integral.
- *What must be true:* the two parts are defined as objects, not by the form in which a term is written. `E_Q^geo` is
  `𝒟Q` with every change `F′` set to zero, plus `2Q`. `E_Q^resp` is `E_Q − E_Q^geo`.

**F3 (Grok 1): the whole pressure, and the whole minus the expansion, need not have a separation derivative on the
condition the list carries. Confirmed.**
- *The list:* it carries S15's pressure conditions unchanged (`:23–24`). "Every object is exact on S15's domain, and
  each pressure object carries its condition." (`:73`) (M4).
- *The source:* the whole pressure is printed "where `p_br` is defined at `ρ_br`". The expansion exists "where
  `p_br` is defined and differentiable at `ρ_2 + ϱ`" (`S15_decision_list.md:110–114`). P6 supplies
  `c_comp(ρ_br)² ≡ dp_br/dρ_br` (`S9b_SHARED_PHYSICS.md:367`) and makes `p_br` "a general function of `ρ_br`"
  (`:370`) (M4).
- *The evidence:* Grok used `p = |ρ − 2|`, with the total density at the kink (M4).
  - The whole pressure is defined there, and its one-sided quotients are `area` and `−area`, so
    `separation_derivative_exists: False`.
  - For the whole minus the expansion, the quotients are `2` and `0`, also `False`.
  - With `p = ρ²`, the derivative exists and matches the chain rule.
- *What must be true:* see F4.

**F4 (Codex 2 = Grok 2): the response part of the first-order expansion needs a second derivative at
`ρ_2 + ϱ`. Confirmed.**
- *The evidence, two independent examples:*
  - Codex used `p = ρ + |ρ − 1|^{3/2}`. It has `p′(1) = 1` from both sides, so the expansion exists on the carried
    condition. With the remainder changing, the outer `first_order_choice` component gives quotients `∓4π(9√h + 5h)/15`
    and a family derivative of `−∞` from both sides (M5).
  - Grok used `p′ = |ρ − 1|`. The response piece has quotients `ρ_1` and `−ρ_1`, so
    `response_derivative_exists: False`, while the geometric piece, with the base held, is defined (M5).
- *Effect (F3 and F4):* item 2 asks for departures that are not defined everywhere on the conditions the list carries,
  and line 73 calls them exact there.
- *What must be true (F3 and F4):* each departure object is stated where the separation derivative of its object
  exists, and it carries the condition under which it does, labelled where it is printed. The measured cases are
  instances of this:
  - the whole pressure, and the whole minus the expansion, need `p_br` differentiable at the total density;
  - the response part of the expansion needs `p_br′` differentiable at `ρ_2 + ϱ`.

## Examined and not filed

- **The operator.** Grok's script confirms that `E_Q = 𝒟Q + 2Q` is zero for degree −2 and for a radial vector. It is
  nonzero for a vector of magnitude `d⁻²` whose direction turns with `d`, and nonzero when the departure vanishes at
  one separation only. An object without `X_1` has `E = 2Q` (M6).
- **The family.** Holding `Ẋ_2` and `Ṙ` keeps them in `Q`; it does not set them to zero (M6).
  - The list's definition differentiates `X_1` and `Ẋ_1` as separate quantities at one instant. A time path gives a
    different derivative, and the list chooses the instant (M6).
  - Scaling the outer radius with the separation adds a term that holding it does not (M6). The list holds the radii.
- **The labelled split of the totals.** The same value placed in drain 1's profile or in a remainder has different
  derivatives along the family (Grok `split/derivatives_differ: True`; Codex's decomposition check, M3). The
  geometric part is therefore relative to S15's labelled split. Neither leg filed this. The list conditions its claim
  on the labelled choices, and `E_Q^resp` carries what a remainder's change adds.
- **The profile family.** The legs found it labelled, limited to items 3 and 4, and elementary on the sphere, including
  the logarithmic special exponent. Its truncated series has zero departure while the exact object does not, so the
  list's refusal of expansions in place of exact objects matters (M6). The legs found no prior-art result on the
  builder's side of the list.
- **Coverage and the decisions.** Grok found that the objects cover S15's reading bullets (`:145–151`). Item 5 keeps
  each drain's velocity, linearly (Grok's report).
- **Stray files.** A Grok run started in the wrong directory left `02_pressure_domain.stdout` (empty) and
  `02_pressure_domain.stderr` in the repository root. The stderr reads "can't open file" (M7). The leg reran the
  script in its scratch (M6). I listed the two files, then moved them to `_scratch/s15/reading_dl/stray_p1/`. After
  that the tree matched the pre-pass status (M0).

## The fold

Revision 2 makes F1–F4 true, once (G2). It changes what is computed (F2, F3, F4) and what may be claimed (F1). Its
spec content therefore gets a second Codex + Grok pass, reviewed until clear.
