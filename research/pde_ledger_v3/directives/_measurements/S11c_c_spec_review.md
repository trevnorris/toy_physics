# S11c-c SHARED_PHYSICS spec — decision-leg review record (2 legs, rule 7 TRIGGER)

Artifact `directives/S11c_c_SHARED_PHYSICS.md` (v1, authored 2026-09-03). Orchestrator-written physics-bearing
spec ⇒ 2 decision legs (Codex + Grok, document review). Prompt `directives/_legs/S11c_c_spec_review.md` (11794 b).
Logs `~/.s11_build/S11c_c_spec_review/{codex,grok}.log`. Both exit 0; each saved its own SymPy derivations +
literal stdout (rule: a prose re-derivation is worth nothing):
- Codex: `directives/_legs/S11c_c_physics_review_derivations.{py,out}`.
- Grok: `directives/_legs/S11c_c_grok_leg_{geometry,closure}.{py,out}`.

## Commands
```
codex exec -c model_reasoning_effort=xhigh --sandbox danger-full-access "$(</abs/prompt)" > codex.log 2>&1   # bg
grok --prompt-file /abs/prompt --cwd /var/projects/toy_physics --model grok-4.6 --effort high \
     --permission-mode bypassPermissions --output-format plain > grok.log 2>&1                                # bg
```

## Verdict split
- **Codex:** NOT safe as-is; needs a physics **re-author + a finer 3-way split** (c1 impermeable curved
  DtN kernel / c2 operator-valued permeable B0c response + resolvent / c3 self-energy fold + re-extraction).
  Driven by (a) review surface + (b) tractability (S11c-b's simpler local operator already OOM'd 30 GB).
- **Grok:** NOT safe as-is; **safe after the 7 listed folds** (F1–F6 local edits, F7 an emit/isolation change);
  does **not** need a full re-author, and the DtN+B0c surface (B0-like, already reviewed in S11b) does **not**
  need splitting **unless** F7's increment-emit is refused.

## CONVERGENT findings (both legs independently → highest confidence)
1. **§2a face-parity leak (Codex #9 MODERATE / Grok F1 top-severity).** My §2a: "faces tilt oppositely ⇒ the
   curvature correction carries a definite parity under `s→−s`." WRONG + leaked (rule 5 + N3). The **lab-w graph
   slope** is odd (`∂h_s/∂x=(s/2)∇W`), but `Z` is defined on the **outward normal**, whose in-plane tilt is
   `−½∇W` for BOTH faces → **even** in `s` (Grok geometry.out OBJ2/4/5; Codex finding 2). ⇒ a "definite parity"
   licenses an ODD correction → both engines manufacture a spurious `δW↔ζ_c` bulk mixing across the DISCONNECTED
   half-spaces and agree. ORCHESTRATOR-VERIFIED (rule 13): `n̂_+ ∝ (−½∇W,+1)`, `n̂_- ∝ (−½∇W,−1)` — in-plane
   part even. **Fold:** delete the implication; parity is whatever the completed per-face operator returns.
2. **DtN is a two-momentum in-plane OPERATOR, not a single-`k` symbol (Codex #1 BLOCKER / Grok F2; rule 17).**
   The first-shape-order curved DtN is a nonlocal kernel `Z_s(ω;k,k′)` with TWO branch legs `q(k),q(k′)`; my
   singular-`q_out` / "three regimes" language invites both engines to diagonalize in `k` (freeze one leg) →
   omit propagating↔evanescent mixing (Codex derivation 1: `δN∝f̂(k−k′)[1−(K²−k·k′)/(q(k)q(k′))]`, depends on
   BOTH legs; Grok OBJ11: `x·e^{ik₀x}` is not an eigenmode of `d/dx`). **Fold:** name it `Z_s(ω;k,k′)` / a
   position-space Schwartz kernel; require the 3×3 regime pairs (grazing on input leg / output leg / both); the
   first-shape-order piece keeps `W_bg(x)` and its jets as `x`-functions (physical-space products or full `(k,k′)`
   kernel); "do not freeze the regime" is insufficient — forbid the single-`k`/locally-constant-`W_bg` freeze.
3. **`Λ_X` misplaced in the θ-row closure block (Codex #6 MAJOR / Grok F4).** My §3c said the θ-row's
   `closure_shape_deriv` carried `Λ_A𝒜+Λ_V V` AND `Λ_X𝒜`. It carries only `Λ_A𝒜+Λ_V V`; `Λ_X𝒜` enters through
   **T-d traction/virtual work** (mechanical rows) — Grok closure.out `T_i_closure_contains_Lambda_X=False`,
   `T_d_traction_contains_Lambda_X=True`; Codex derivation `(M+J)−[J−Λ_A𝒜−Λ_V V]=M+Λ_A𝒜+Λ_V V` (no Λ_X). ⇒ an
   engine that treats the illustration as the fold leaves the mechanical bulk load symbolic → self-energy misses
   the B0c→B2a mechanical channel. **Fold:** substitute closed `(δp_s,J_s,t_s)` into EVERY term-origin-marked
   face/flux/**traction**/closure slot; state `Λ_X` is in T-d traction (thickness+in-plane rows), not the θ-row;
   replace "non-face rows" → "non-face **terms**" (mechanical rows legitimately carry face terms).
4. **Dissipation: operator Hermitian part + independent flux-vs-work, not same-equation Re/Im (Codex #5 MAJOR /
   Grok F6).** For a nonlocal/operator `Z`, dissipation is the Hermitian part under the true-area boundary
   pairing `⟨f,g⟩_∂Ω=Σ_s∫a_s f_s* g_s d³x`, `Z_diss=(Z+Z†ₐ)/2` — NOT the entrywise `Re` (Codex matrix
   counterexamples both directions). AND (Grok F6): inherit S11b's INDEPENDENT pressure-work check
   (`S11b:463-474`) — slab `d(T+U)/dt` off-shell vs independently-computed outgoing bulk flux, real ω /
   propagating / impermeable / Λ_X⁰=0 — because S11c-c is the first curved close and S11c-b's #90 face sign is
   still cross-engine-UNVALIDATED. **Fold:** both — the operator-Hermitian definition AND the independent
   flux-work route; do not state the sign.

## Codex-unique (orchestrator-verified real)
5. **#2 BLOCKER — extract/eliminate DON'T commute.** §3c substitutes closed loads into the ALREADY-EXTRACTED
   `S11CB_COUPLING_KERNEL`; closing `δp` can create a cross-sector term from a term diagonal while `δp` was
   independent (counterexample `R_x=x+p, p=αy`: extract-then-sub → 0; sub-then-extract → α). **Fold:** close
   the FULL `S11CB_SLAB_OPERATOR` first, then RE-EXTRACT the weak off-diagonal blocks from the closed operator
   (repeat S11c-b §3c on the closed operator); substitution into the imported open kernel is at most a residual
   after proving per-channel commutation.
6. **#3 BLOCKER — the §5a face-flattening route is secular at infinity.** The global `w'=(w−ζ_c)/(W_bg+δW)`
   over the infinite half-space gives a transformed Laplacian with coefficients `∝w'∇W_bg` and an outgoing wave
   a secular first variation `i q W_0 w' f(x)` → imposing the flat outgoing condition in `w'` can select the
   wrong radiation branch → route 2 is not a well-defined construction. (Grok's "clean axis" ASSERTED the two
   routes independent but did NOT check infinity — Codex's derivation wins, rule 13.) **Fold:** use a
   cutoff/Hanzawa flattening (= face map near the boundary, identity at infinity) with its transformed radiation
   condition stated, OR a boundary-integral/layer-potential second route.
7. **#7 MAJOR — the algebraic locus protocol can't express operator non-invertibility.** For curved (operator)
   `Z`, "loss of solvability" is a Fredholm/resolvent condition on the function space + profile spectrum — the
   profile-conditioned S11c-d work — not a scalar coefficient locus. **Fold:** emit the FORMAL
   noninvertibility/Fredholm condition in S11c-c; apply the CAS §8 locus protocol only to finite-dimensional
   algebraic degeneracies / the flat Fourier-diagonal regression; defer the profile-conditioned singular set to
   S11c-d.
8. **#8 MAJOR — §5d is a schema mismatch, not a residual.** Comparing a 3-regime keyed map to a 1-regime object
   is a coverage artifact (zero on the propagating case, mismatch elsewhere), and it doesn't test the S11b
   opposite-sheet error. **Fold:** replace with same-domain one-sided corruptions — flip `q_out(k′)→−q_out(k′)`
   on ONE kernel leg and recompute; separately `q_out(k′)→q_out(k)` to detect the momentum freeze (F2); check
   real-axis flux/decay independently, then continue the selected branch without re-selection.
9. **#10 MODERATE — §7 omits the N1/N8 chain output.** No requirement for `scripts/S11c_c_exports.py`,
   `BUILD_INPUT_DIGESTS`, the D3 round-trip, `_RELATIONALS` revival. **Fold:** point verbatim at F1–F9/N8 and
   require the exports file, its exact digest set, round-trip, relational revival, and the own frozen comparator.
10. **#4 BLOCKER-as-stated — grazing vs the drain domain is non-uniform.** §2b conditions every object on
    `|q·v_bulk_normal_0/ω|≪1`, but that → 0 as `q→0` (grazing, §3a) while the impedance/root correction
    DIVERGES (Codex finding 7: `limit_{q→0} stated_smallness = 0`, root correction → ∞). ⚠ This does NOT
    contradict Grok's "convective operator not needed at O(σ_W)" (below) — different questions. **Fold:** state
    grazing as the strict `v_bulk_normal_0=0` result; for nonzero drain exclude the boundary layer
    `|ω v_bulk_normal_0|/(c_s0²|q|)≁≪1`; retain the independent subsonic-speed condition (S11b B9); make the
    order of limits explicit. (Engines still compute the rest-frame object; this refines the STATED domain.)

## Grok-unique (orchestrator-verified real)
11. **F3 — §5a "or" lets the only tilt test of `Z` be skipped (N6).** `Σ_E` doesn't enter `Z` (`ρ_m` constant),
    so an engine choosing the `Σ_E` arm never runs the slope-flip on `Z` → F1's leaked parity goes untested;
    §5c can't catch it (vanishes at η,σ_W=0). **Fold:** require BOTH source mutations, slope-flip MANDATORY on
    `Z`; the `Σ_E` mutation applies to the permeable response / self-energy, not as an alternative to the tilt test.
12. **F5 — no zero-jet (`σ_W→0`, `η` retained) regression on `Z`.** A "bulk cavity between two planes" error
    (treating the disconnected two-face problem as a finite-gap cavity) is `O(η)`, jet-independent → §5b (jet
    ablation) residual 0 AND §5c (η,σ_W→0) vanishes → both engines can share it (N6's vacuous limit on the `η`
    bookkeeper). **Fold:** add an independent zero-jet regression — form the `(σ_W→0, η` retained`)` curved-DtN
    operand and the S11b flat `Z`, emit both + residual; ⛔ do not state the residual is zero.
13. **F7 — the self-energy residual can't isolate S11c-c from S11c-b's DEFERRED residual.** Folding closed loads
    into the imported `S11CB_SLAB_OPERATOR` mixes new bulk physics with S11c-b's cross-engine-unvalidated sign
    conventions (#90 flags). **Fold (prefer over a rewrite):** emit `S11CC_SELF_ENERGY_OPERATOR` as the
    SUBSTITUTION INCREMENT (closed operator − imported-operator-with-symbolic-slots) so S11c-b content cancels
    in the S11c-c residual; if isolation isn't written, defer the fold to a thin sub-step. This also UNBLOCKS
    the fold from the ≥64 GB S11c-b residual.

## Clean axes both legs CONFIRMED (folded corrections aside)
- The two inherited S11b caveats (propagating `Re Z` = bulk radiation; per-face inertial loading against OUTWARD
  acceleration) are legitimate setup caveats, NOT leaked curved answers.
- The two exterior half-spaces are DISCONNECTED; `Z` is per-face; face parity structures the `(δW,ζ_c)`
  ASSEMBLY of two independent impedances, it does not couple the half-spaces.
- The B0b→B0c composition (`δp_s=Z·v_bulk,s` outward-normal → `v_bulk,s=V_s+J_s/ρ_m` → §1d closure) is correct;
  T-a/T-a′/T-a″/T-b/T-c/T-c′/T-d/T-e/T-i consumed at the correct seam (T-i already excludes the bulk solve).
- N11a: dropping the convective operator IS legitimate at first shape order — Grok geometry.out OBJ10:
  `n̂_s·ŵ − s = O(σ_W²)`, so the drain-projection correction is 0 at `O(σ_W¹)`; the tilt×convection cross term
  is higher order inside N11a's domain. (Codex #4/#10 is about the STATED grazing domain, not the operator.)
- N12 truncation (first order ε; first shape order η, σ_W; `W_bg` Hessians kept as first amplitude order;
  `|∇h|²=O(σ_W²)` correctly not first order — Grok OBJ6) is right. N5/N8/N13/N7/N14 correctly stated.
- §5c is correctly a regression only; §5b form-ablation is physical (μ_R,bg ablation structurally absent from
  `Z` — honest zero); the two §5a routes ARE independent in principle once route 2's flattening is repaired (#3).

## The container decision (bring to the user — N2 + rule 11)
Everything in "CONVERGENT" + "Codex-unique" + "Grok-unique" folds the SAME way regardless of container. The ONE
open decision is the CONTAINER: (A) one revised S11c-c spec (Grok — increment-emit isolates the fold; matches
S11c-a/b; tractability handled per-engine as S11c-b did) vs (B) Codex's 3-way split c1/c2/c3 vs (C) a 2-way
split (c1 = DtN + permeable response; c2 = self-energy fold + re-extraction) as a middle ground. Scale/split
count is the user's call (rule 11 / N2). Fold AFTER the container is chosen, so the objects land in the right unit.
