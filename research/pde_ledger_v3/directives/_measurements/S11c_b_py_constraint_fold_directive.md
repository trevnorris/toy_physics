# S11c-b PY constraint-fold — decision + build review record (spec-pin B follow-on)

Directive: `directives/S11c_b_py_constraint_fold_directive.md`. Target: `scripts/S11c_b_brane_operator_sympy_audit.py`
(SymPy engine) + a §3b amendment to `directives/S11c_b_SHARED_PHYSICS.md`. This folds the pinned (B) constraint
reduction into PY: the slab U/thickness rows become the constraint-reduced equations (the imported non-uniform
`virtual_constraint` eliminates virtual θ, the same `μ_θ` feeds both reactions), the θ-row becomes the imported
sourced mass-evolution, `μ_θ` stays a separate held-fixed operand, and the jet-depth cascade is raised so the
reaction is retained. Raw `operator_from_density` / `committed_strong_rows` are preserved for #88.

## Decision-leg review (rule 7 — orchestrator-written directive → Codex + Grok)

Prompt `directives/_legs/S11c_b_py_constraint_fold_decision_review.md`. Raw logs
`~/.s11_build/S11c_b_constraint_fold/decision_{codex,grok}.log` (`.log` gitignored). **v1 was REJECTED by BOTH legs,
independently and computation-backed.** Convergent defects (orchestrator-verified against the code, rule 13):

1. **Uniform-constraint freeze (rule 17).** v1 folded `S11c_b_SHARED_PHYSICS.md:143`'s `(uniform linearisation)`
   `δ_vθ + δ_ve_W + ∇_x·δ_vu = 0` as the constraint — freezing `W_bg`. The binding constraint is the MATERIAL
   `δ_vΣ_mat = 0` with live `W_bg`, already imported as S11c-a `virtual_constraint`. Verified via my own last-session
   consult audit (`_measurements/s11c_b_jet_depth_consult_codex/s11cb_advisor_term_audit.py:60`): the constrained U
   row is `E_u,a + ∂_a(μ_θ) − (∂_a W/W)μ_θ`, carrying a non-uniform term the uniform formula drops. Grok's live-
   profile spot-check independently showed the extra `W_0/W_bg` and `∇W_bg/W_bg` coefficients.
2. **θ-row double-count.** `evolution_mass_balance` (S11c-a:1313-1318) is ALREADY the full sourced sum
   (`DENSITY_TIME + VELOCITY_DILATATION + BACKGROUND_ADVECTION + TRUE_AREA_FACE_FLUX`); v1's "add
   `ADVECTIVE_MASS_OPERAND`" would duplicate advection.
3. **Wrong interface / #88 breakage.** `operator_from_density` has no branch/route arg and is #88's raw-EL reference;
   changing it in place silently redefines #88's standing measurement. Fold belongs in a SEPARATE layer in
   `build_operator`.
4. **Depth cascade.** `build_kernel` differentiates the §3b U-row once more (`operator_divergence`), so strong→3
   forces coupling→4 and control→5 with extended jet tables — not the single bump v1 wrote.
5. **Provenance over-reach + rule-5 leaks.** A new top-level `CONSTRAINT_REACTION` origin conflicts with WL (reaction
   lives in `BULK_ENERGY`); v1 also leaked a prior measurement verbatim and stated a typed `−μ_θ` thickness term.

**v2** folds both legs' convergent replacement text (rule 7: one pass, no re-legging). Leak-gated clean.

## Build (Codex, `--sandbox danger-full-access`) — deliverable verification

Log `~/.s11_build/S11c_b_constraint_fold/build_codex.log` (266k tokens). Verified against the DIFF, not the report:
- `git diff --stat`: `S11c_b_SHARED_PHYSICS.md` +16, engine +859/−87. Both target files only.
- Engine `ast.parse` OK; depth constants now `STRONG=3, COUPLING=4, DEPTH_CONTROL=5`.
- **`operator_from_density` and `committed_strong_rows` BYTE-IDENTICAL to HEAD** (checked by extracting both function
  bodies from working tree vs `git show HEAD:` — the #88 raw-reference preservation invariant holds).

## Build-leg review (rule 7 — Codex-written → fresh Claude agent + Grok)

Prompt `directives/_legs/S11c_b_py_constraint_fold_build_review.md`. Evidence
`~/.s11_build/S11c_b_constraint_fold/{claude,grok}_buildleg_evidence/` (the legs' own ablation scripts + literal
stdout). Both legs independently ablated a COPY; working tree untouched.

**The emitted physics rows are CLEAR on both legs** (each with FORM ablations that bite):
- Reaction is COMPUTED from the imported non-uniform `virtual_constraint`, not typed/frozen: swapping the source to
  the uniform three-term relation MOVES the rows (Grok: U `[340,340,340]`, e_W `170`); a wrong constraint moves the
  rows. `live_reaction − ε∇μ_θ ≠ 0` while `uniform_reaction − ε∇μ_θ = 0` — i.e. the v1 freeze was NOT committed.
- θ-row = imported `evolution_mass_balance (branch, DELTA_W, rep)`, no double-counted advection (`(θ+ADVECTIVE)−mass`
  has 3 extra terms, confirming the build did NOT add it). `μ_θ` unchanged/held-fixed/separate.
- Depth cascade load-bearing: cap→2 loses all order-3 atoms; coupling cap→3 loses 174 order-4 terms.
- #88 raw refs intact; provenance stays `BULK_ENERGY` (no new top-level origin); per-term mass row not repeated;
  corruption propagates at the constraint/evolution SOURCE; §5c like-with-like (U/E_W vs S11b `inplane_eom`/
  `thickness_eom`); transverse dispersion free of `μ_θ`; no asserts, no value leak. Full 4-case smoke PASS (564s).

**FINDING (Claude agent, Probe 2) — folded.** The two-route (elimination vs Lagrange-multiplier) residual
`CONSTRAINT_FOLD_ROUTE_RESIDUAL` is **tautological operand theatre**: for a constraint linear in `δ_vθ`,
`λ = −ε·μ_θ/(∂C/∂δ_vθ)` and back-substitution cancels the `δ_vθ` term, so the two routes are algebraically
identical and the residual is `0` by construction — even for a WRONG constraint. My directive's claim that the two
routes are "genuinely independent … not operand theatre" was FALSE. This does NOT corrupt the emitted rows (both
legs confirm the rows are correct and input-driven; the genuine independent check is the blind cross-engine
comparator). The two legs disagreed here — Grok's one-sided CODE corruption moved the residual and read it CLEAR,
but that only tests implementation-separateness, not route-independence; the agent's wrong-input test is decisive.

**Orchestrator CAS verification (rule 13, resolving the disagreement):**
Command: `python3 -c '<generic linear constraint, both routes symbolically>'` (engine-independent).
Literal stdout:
```
elimination - multiplier (should be 0 for ANY a,b,c): 0
WRONG-constraint residual (still 0 => tautology): 0
WRONG rows DIFFER from correct? elim vs elim_w: True
```
⇒ residual is `0` for any linear constraint (tautology confirmed); the reduced ROWS differ under a wrong constraint
(input-driven, correct). Grep confirms nothing asserts the residual `==0` (not in `row_residual`/comparator; engine
has 0 asserts).

**Fix applied (this commit):** the tautological residual is RELABELED
`CONSTRAINT_FOLD_ROUTE_RESIDUAL → CONSTRAINT_FOLD_TRANSCRIPTION_RESIDUAL` with an in-code note that it is a
transcription-consistency check, not an independent physics route; the directive's claim is corrected. Rename is
label/comment-only — the computed rows are unchanged (parse OK, no stale refs, module imports).

**FLAG (Claude agent, Probe 3) — deferred to #88 re-adjudication, no fix now.**
`HESSIAN_FREEZE_STRONG_ROW_RESIDUAL` now compares the FOLDED `live_strong_rows` against the RAW
`committed_strong_rows`, so it is the reaction (nonzero), no longer a zero regression guard. This is a deliberate
consequence of routing controls through the folded path while keeping `committed_strong_rows` raw. Nothing asserts
it `==0` and it is not consumed by the comparator, so it is a printed value; its interpretation belongs to the
DEFERRED #88 KINETIC+θ re-adjudication (`directives/S11c_b_88_blast_radius_build_directive.md`).

## Status

Physics CLEAR (both legs); one verification-hygiene finding folded; one flag deferred to #88. Follow-on: single-case
in-band cross-engine `row_residual` attempt on this box (~16 GB WL peak), then #90, then the ≥64 GB full run.
