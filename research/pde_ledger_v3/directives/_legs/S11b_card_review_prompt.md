# Independent review — the re-pointed S11b ledger card (physics-bearing `.tex`)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/paper/steps/S11b_interface_coupling_law.tex`

This is a ledger card (paper record of the S11b step), just re-pointed by Codex. It is the INTERPRETATION /
record layer: every physics statement must be faithful to the committed step record, which is the source of
truth. Your job is to catch any claim the card makes that the step record does not support, any internal
inconsistency, and any load-bearing element that was dropped.

## What to check
Read the SOURCE OF TRUTH FIRST — the step record `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md`
and the spec `directives/S11b_SHARED_PHYSICS.md` — form your own view of what the card should say; only THEN
read the card. Check:

1. **Faithfulness.** Every physics statement in the card traces to the step record; ⛔ the card states no
   result the step record does not carry, and no overclaim.
2. **Transverse consistency (the subtle one).** The card's energy list uses the `(∇·u)²` representative (no
   `st²`, no `μ_S`). Therefore its transverse working stiffness must be `μ_⊥ = μ_R` (NOT `μ_R + μ_S/2`), and
   `Im ω = 0` must be **conditional** on `μ_⊥ = μ_R ≥ 0` (§5 assumes no positivity; `μ_⊥ < 0` gives a growing
   root). ⛔ Flag it if the card displays `μ_R + μ_S/2` as the working stiffness (that would use `μ_S`
   undefined and not follow from the displayed energy) — `μ_R + μ_S/2` may appear only as a parenthetical
   cross-engine (SymPy-representative) note. ⛔ Flag it if `Im ω = 0` is stated unconditionally.
3. **`Λ_A⁰` defined + slice map.** Is `Λ_A⁰` now defined (DC affinity-channel coefficient), and the static
   slice map `Λ_p⁰ = −Λ_A⁰/ρ_m` correct (verify from `𝒜 = μ_s − δp/ρ_m` on the `μ_s=0` slice)? Is the
   `Λ_p⁰ = 0` pure-velocity qualifier present on the finite-memory departure? Is `τ_A` narrowed to
   "unrestricted by the passivity conditions (beyond `τ_A ≥ 0`)", not implying negative relaxation times?
4. **Background-flow standing limit.** Is the uncarried `O(v₀|q_n|/ω)` correction recorded, with the failure
   condition `|q v₀/ω| ≳ 1` (`k c_s0 ≫ ω` necessary, not sufficient)?
5. **Verification scope.** Does it describe the frozen-T7 cross-engine comparison (no physics contradiction;
   format / coefficient-basis / naming / convention / denominator-clearing / coverage gaps) rather than the
   old single-engine `VERDICT: PASS` framing, while KEEPING "verification comes from reviewer derivations,
   not engine agreement"?
6. **Energy basis.** Ten invariants under §5's total-divergence quotient; representative non-unique; the X-1
   over-count note; the existing `(∇·u)²` five-term list kept. Faithful?
7. **Preservation.** The "region not prohibition" framing, the two refuted predictions
   (two scalar dof; `ρ_m/(2α)` per face not `ρ_m/α`), the `\paragraph{Verification scope.}` handled as a
   PLAIN paragraph (⛔ NOT `\stagefield{Verification}{...}` — that field is suppressed by default,
   `paper/macros.tex:19`), and all pre-existing "Known limits" items — are they intact (none dropped)?
8. **LaTeX/macros sanity.** Does it use only macros defined in `paper/macros.tex` (`\stagefield`,
   `\claimstatus`, `\StatusText`, `\StageFile`, …)? No obviously broken markup.

## What you are handed
- The card (path above).
- The step record `steps/S11b_interface_coupling_law.md` (SOURCE OF TRUTH) and the spec.
- `paper/macros.tex`. Full repo; ⛔ no do-not-read list.
⚠ You are NOT handed the re-point directive: a card can satisfy its directive and still misrepresent the step
record, and that is what this review catches.

## Required method — DOCUMENT review
Read the step record + spec + macros first; form your view; only then read the card. Quote both the card and
the step record for every finding. Where a relation is checkable (the slice map `Λ_p⁰=−Λ_A⁰/ρ_m`, the
transverse `μ_⊥` in the `(∇·u)²` representative), verify it (a short CAS check is welcome; save script +
stdout to named paths if run). ⛔ Do not modify the working tree. Physics filter: report a finding only if it
catches the card misrepresenting the step record, an internal inconsistency, a dropped load-bearing element,
or broken markup — not wording preference.

## Output
Per item (1–8): verdict + exact quotes from the card and the step record (+ script path/stdout for any CAS).
End with a one-line bottom line: is the card faithful and internally consistent as written, or the specific
corrections it needs.
