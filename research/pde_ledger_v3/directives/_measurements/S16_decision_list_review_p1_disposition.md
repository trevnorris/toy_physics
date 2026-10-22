# S16 decision list, review pass 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_decision_list_review_p1_lookups.md`, beside this file
(generator `_scratch/s14a/gen/s16_dl_p1_lookups.sh`). "N-n" cites a section of that file. Every computed result
quoted here is a leg's saved stdout.

## The pass

- **Artifact.** `directives/S16_decision_list.md` (sha256 `e1fa73f9…`) and its setup lookups (`a0a58f32…`),
  uncommitted on `79511272`. Both matched the baseline after the pass. The working tree was unchanged apart from the
  lookups file (N1).
- **Role.** A pre-builder decision list: one Codex + Grok pass, folded once (G2). The list also names the objects
  the build computes, so that content is reviewed until clear (spec row).
- **Prompt and legs.** `_scratch/s16/dl/dl_review_prompt_p1.md`, sha256 `5b556b41…` (N1). Both legs received it
  identically:
  - Codex (`gpt-6.1-sol`, xhigh, workspace-write, working directory its scratch):
    `/var/projects/_review_scratch/s16-dl-p1-codex/REVIEW.md`, exit 0;
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s16-dl-p1-grok/REVIEW.md`, exit 0.
- **The first Codex launch** ran from the repository root. It stopped before any witness ran, because the guarded
  runner could not reach the systemd bus, and it gave no verdict (N2). It is not counted as a report. Its directory
  was moved aside to `s16-dl-p1-codex-blocked/` before the second launch. The second launch used the working
  directory that earlier decision-list legs used. It received the same prompt and nothing from Grok.
- **Verdicts.** Codex: FINDINGS, 5. Grok: FINDINGS, 4 (N3).

## Findings

Findings from both legs that concern the same object are merged. Each is attributed.

**D1. The split into a regular part and the rest is not determined by the fields on the region.** (Codex 1, Grok 1.)
**CONFIRMED.**
- **The list** (`:70–75`): "the build may write a field on the region as the sum of two parts: a part that is
  regular across the region, the core included; the rest" (N4).
- **The sources.**
  - The list's own region: the core lies inside the inner surface, and its content stays with S22 (`:65–66`).
  - The accepted part 1 record: "The core is excluded, with S22" (`steps/S14a_1_projection.md:52`) (N5).
  - None of the three selected options contains a decomposition.
- **The evidence.**
  - Grok: one field on the annulus `1 ≤ r ≤ 2` has two regular parts, whose gradients at the centre are `(0, 1)`
    and `(1, 0)`. The form holds for one and fails for the other (N6, D1).
  - Codex: two extensions inside the hole agree on the entire shell. Their central values are `(1, 2)`, and their
    contributions through second order are `(7/3, −302/15)`; the remainders compensate. Constant splits `1 + 0` and
    `2 − 1` give the same total (N6, D1).
- **Effect.** The expansions of item 3 and the passive factor of item 4 depend on data the region's fields do not
  contain. They are evaluated at `X_c`, inside the core S22 owns. Two engines can disagree with neither wrong.
- **Cause.** The orchestrator added the split before review, to make item 4 posable.
- **What must be true.** No object of the build depends on data that the supplied fields on the region do not
  determine. Any field of the surroundings is either supplied by a decision or not computed.

**D2. The passive factor is not determined, and its condition holds for every force.** (Codex 2, Grok 2.)
**CONFIRMED.**
- **The list** (`:100–104`): "Suppose its leading term is a scalar of the region times the gradient, at `X_c`, of a
  field of the regular part. The build prints that scalar as the passive factor" (N4).
- **The source.** The selected option: "It also prints the condition under which that force takes the form of a
  gradient times a mass, without assuming it." The guardrail asks that each mass be "Define[d] … independently"
  (`notes/muonium_gravity_corrections.md:21–25`) (N5).
- **The evidence.**
  - Grok: rescaling the field and the scalar leaves the product fixed, `scalars_differ = True`. A second force also
    fits. Flow and stress have shape `(3, 3)`, not that of a gradient of a scalar (N6, D2).
  - Codex: the same integrated force gives passive factors `9317/3000`, `9317/1500` and `−9317/1000` for `z`,
    `z/2` and `−z/3`. A force density with nonzero curl passes (N6, D2).
- **Effect.** The condition constrains no field, and the passive factor and its difference from the inertial
  factor are not unique. The intended claim's "conditions under which the point-particle form would follow"
  (`:162`) would rest on it.
- **What must be true.** Any passive factor the build prints is determined by the supplied data, with its field and
  its normalization named. Where the data do not determine one, none is printed as a factor.

**D3. The "force" includes terms from the motion of the control surfaces.** (Codex 3.) **CONFIRMED.**
- **The list** (`:100–101`): "The force on the region is every term of item 2's relation other than the inertial
  term" (N4).
- **The evidence.**
  - Codex: stationary material with `F = G = Π = f = 0`, and moving surfaces. Then `dP/dt = 0`, and the integrated
    force density and the momentum crossing both surfaces are zero. Yet item 2's remainder is `(9317/3000, 0, 0)`,
    "produced entirely by the moving selection of stationary material" (N6, D3).
  - Codex: with `G = λF` the collected coefficient of the acceleration is `λM`, while the displayed prefactor is `M`.
    Grok: the coefficient of `d²X/dt²` is `M(t)` by identity (N6, D3).
- **Effect.** The "force" would count the change in which material the surfaces enclose as a force. The inertial
  factor is `M` by the definitions, and the balances leave open how momentum relates to motion.
- **What must be true.**
  - The force object is the momentum supplied to the region, as the momentum balance gives it. The terms that come
    from which material the moving surfaces select are kept apart from it, as kinematic terms.
  - The inertial factor is not presented as determined where the relation between momentum and motion is open.

**D4. A constant mass per particle fixes a ratio the accepted O2 spec keeps live.** (Grok 3.) **CONFIRMED.**
- **The list** (`:61–62`): "`m`, a supplied constant mass per unit of `n`" (N4).
- **The sources** (N5).
  - The accepted O2 spec (`72866fcf`) says of `ℐ_br^live`: "No identification with `ρ_br V`, a quadratic live
    kinetic energy, or a constant inertia is supplied" (`directives/O2_SHARED_PHYSICS.md:133`).
  - It says of `ℛ_br`: "General brane-density response … It is input to the density in mass/momentum/energy
    accounting; it is not inferred from the bulk EOS or slab factorization" (`:138`).
  - The part 1 record: "No mass … is delivered" (`steps/S14a_1_projection.md:15`).
- **The evidence.** Grok: on the same number density, a constant and a live ratio give different masses and
  centres. The momentum residual does not contain `m` (N6, D4).
- **Effect.** The list freezes the ratio, a law the selected options do not contain. The inertial factor would then
  be the number content times a constant, by definition.
- **What must be true.** No mass-to-number ratio is fixed. The region's mass, and how its momentum relates to its
  motion, stay as open as O2 keeps them, unless the user decides otherwise.

**D5. A general second moment does not define one width `a_WT`.** (Codex 4, Grok 4.) **CONFIRMED.**
- **The list** (`:95–97`): "`a_WT` is the width that the second moment defines" (N4).
- **The sources** (N5).
  - The source defines the width for "the isotropic Gaussian": `Q_{ij}=\frac{Ma^2}{2}\,\delta_{ij}`
    (`4d_1pn_full.tex:761–765`).
  - The plan states the same relation (`V3_STEP_PLAN.md:732–733`).
- **The evidence** (N6, D5).
  - Grok: three standard scalar reductions differ, and the Gaussian relation fails on an anisotropic moment. Keeping
    the dipole changes the RMS value.
  - Codex: on the density `1 + x²`, no single scalar satisfies all three Gaussian equations.
- **What must be true.** The build prints the second-moment tensor. It defines no scalar width from a general tensor
  without a stated functional, and it does not reduce `a_WT` from the region's moments.

**D6. The normalized objects are requested on a domain that includes zero mass.** (Codex 5.) **CONFIRMED, minor.**
- **The source.** The part 1 record: "The velocity restriction does not restrict the densities, material fluxes or
  balances to a nonzero-density subdomain" (`steps/S14a_1_projection.md:139–142`) (N5).
- **The evidence.** Codex: with no material and live moving surfaces, the mass and the first moment both vanish,
  and the centre is `0/0` (N6, D6).
- **What must be true.** Each normalized object is stated on the domain where its normalizer is nonzero. The
  unnormalized balances and moments keep the full domain.

## Notes, no change

- Both legs read the projected material balance as faithful to the part 1 record (N3 reports; Grok
  `balance_quote.out`, Codex `quote_baseline.stdout`).
- Both legs called the "printed as defined, not as derived" label right for a coefficient fixed by definitions.
- Grok found the owners kept apart, the difference from the plan's `:739–741` stated as the decisions give it, and
  the withdrawn lepton results absent.

## Agreement between the legs

Both legs found D1, D2 and D5. Codex alone found D3 and D6; Grok alone found D4. They do not conflict: Codex's
second D3 witness and Grok's D4 witness are two views of the same open relation between momentum, mass and motion.

## Next

D1, D2, D4 and the inertial part of D3 change what the user's "No potential" and "General balance" decisions can
deliver, so they go to the user first. The rest of D3, D5 and D6 are technical repairs. The list is then folded
once. The fold changes what is computed, so its spec content is reviewed again until clear (pass 2) before any
builder.
