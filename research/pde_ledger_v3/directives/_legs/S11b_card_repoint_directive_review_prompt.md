# Independent review — the S11b card re-point DIRECTIVE

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11b_card_repoint_directive.md`

This is an ORCHESTRATOR-WRITTEN directive that will be handed to Codex to re-point the ledger card
`paper/steps/S11b_interface_coupling_law.tex`. Your review is the gate BEFORE Codex runs (CLAUDE.md rule 7).
The card is physics-bearing (the paper's record of the step), so an error in this directive makes the card
wrong.

## What to check
The card's source of truth is the committed step record
`research/pde_ledger_v3/steps/S11b_interface_coupling_law.md`. Read that step record, the CURRENT card
`paper/steps/S11b_interface_coupling_law.tex`, the decision-list G12 items
(`directives/S11b_unified_decisions.md`), and `paper/macros.tex` FIRST; form your own view of what the
re-point should be; only THEN judge the directive. Check each of the directive's six changes:

1. **Source re-point** — is it right to make `steps/S11b_interface_coupling_law.md` the source of record (A/B
   archival)? Does the step record actually cover the card's content?
2. **G12(c) — `Λ_A⁰` definition + `Λ_p⁰=0` qualifier.** Is the directive's definition of `Λ_A⁰` (DC
   affinity-channel flux coefficient; `Λ_A(ω)=Λ_A⁰/(1−iωτ_A)`) correct and consistent with the card's
   existing `Λ_p`, `Λ_V` kernels and the affinity `𝒜 = μ_s − δp/ρ_m`? Is the slice map `Λ_p⁰ = −Λ_A⁰/ρ_m`
   correct (verify from `𝒜 = μ_s − δp/ρ_m` on the `μ_s=0` slice; check `ZPERM_SLICE_MAP` in the step
   record)? Is `Λ_p⁰=0` the right qualifier for the finite-memory velocity branch?
3. **G12(d) — background-flow limit.** Is the failure condition `|q v₀/ω| ≳ 1` (kc_s0≫ω necessary not
   sufficient) faithful to the step record?
4. **Transverse conditional stability.** Is the correction right — `Im ω=0` only where `μ⊥=μ_R+μ_S/2 ≥ 0`,
   with only the dissipation unconditionally zero, given §5's no-positivity (§0)? Does it match the step
   record?
5. **Verification scope re-point.** Does replacing the `VERDICT: PASS` framing with the T7 comparator result
   (no physics contradiction; format/coefficient-basis/naming/convention + coverage gaps) match the step
   record's adjudication, while keeping the load-bearing "verification comes from reviewer derivations, not
   engine agreement" limit?
6. **Energy-basis note.** Is "ten under §5's total-divergence quotient; representative non-unique; one engine
   over-counted to eleven (X-1), corrected" faithful, and is it right to keep the card's existing `(∇·u)²`
   representative list?

Also: does the directive correctly **preserve** the load-bearing card elements (the "region not prohibition"
framing, the two refuted predictions, the plain-paragraph `Verification scope.` handling since the
`Verification` `\stagefield` is suppressed, the "Known limits" list)? Does it introduce any over-claim or a
result the step record does not carry? Does it wrongly instruct a change that would break the physics?

## What you are handed
- The directive (path above) and its rule-2 twin `directives/_measurements/S11b_card_repoint_directive.md`.
- The step record `steps/S11b_interface_coupling_law.md` (SOURCE OF TRUTH).
- The current card `paper/steps/S11b_interface_coupling_law.tex`.
- The decision list `directives/S11b_unified_decisions.md` (G12) and `paper/macros.tex`.
- `CLAUDE.md`. Full repo; ⛔ no do-not-read list.

## Required method
DOCUMENT review: read the step record + card + G12 + macros first, form your view, then judge the directive.
Quote both the directive and the source for every finding. Where a physics relation is checkable (the slice
map `Λ_p⁰=−Λ_A⁰/ρ_m`, the transverse `μ⊥`), verify it (a short CAS check is welcome; save script + stdout to
named paths if you run one). ⛔ Do not modify the working tree. Physics filter: report a finding only if it
catches the directive specifying something that would make the card wrong, unfaithful to the step record, or
that drops a load-bearing element — not wording preference.

## Output
Per change (1–6) + the preservation check: verdict + exact quotes from the directive and the source. End with
a one-line bottom line: is the directive SAFE TO BUILD THE CARD FROM as written, or the specific edits it
needs first.
