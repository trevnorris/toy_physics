# Decision RE-REVIEW (round 3) — F/G diagnostic BUILD DIRECTIVE. STATIC; ⛔ run nothing, load no `.out`.

## ⚠ Paths
Working dir = `/var/projects/toy_physics` (repo root); all `directives/…`, `scripts/…`, `_measurements/…` paths are
under `research/pde_ledger_v3/`. Artifact: `research/pde_ledger_v3/directives/S11c_c2_FG_diagnostic_build_directive.md`.
⛔ Review THAT F/G diagnostic directive — ⛔ not any WL-engine / T7 / export-chain directive. Physics-bearing →
**review-until-clear**. You are ONE of two independent legs.

**Round 3.** Round 2 (Codex-sol + Grok) confirmed most folds but returned NOT-READY on **7 items**; the orchestrator
folded them. Verify each is now correctly + completely fixed and introduced no new defect. ⛔ No fictional-script
ablation. A real problem is a finding; ⛔ do not rubber-stamp.

## The round-2 items that were folded — verify each in the current directive
1. **F2 entry point.** The uniform map must now reach **both** the slab rows **AND** `build_face`/`kernel_bridge`
   (via `overrides=`/specialized `Inputs`), ⛔ `input_map` alone forbidden (it leaves `W_bg_d*`/profile hats live in
   `χ` — `audit.py:541–546`), and the c1 Fourier hats zeroed. Present and correct?
2. **F3 weak-kernel checker.** Now `Integral.function` **AND** any non-`Integral` remainder (a true weak-zero can
   collapse to no-`Integral` via `audit.py:472–473`; an absent `Integral` is NOT weak-zero). Does F3's checker now
   also handle F6's polynomial (no-`Integral`) sentinel?
3. **F6.** Now emits the checker **result** (baseline vs injected), ⛔ no longer asserts "returns nonzero"?
4. **F8 + F2 citations.** The `SHARED_PHYSICS:385` "increment must vanish" leak removed (F8 cites only decisions N6
   as smoke-test limitation; F2 no longer cites the `:380–388` range)? Any remaining click-through to the F answer?
5. **Per-face pin + `model['faces']`.** Is the per-face set now PINNED (`P_s`, `S_{P,s}=SLAB−SLAB|_{P_s=0}`,
   `I_closed,s=extract(S_{P,s}[χ_s])`, `I_bare,s=−extract(S_{P,s})`, `I_s=I_closed,s+I_bare,s`, `R_add=I_direct−Σ_s I_s`),
   is `model['faces'][s]` now correctly described as `S_{P,s}[χ_s]` (= `I_closed,s`, ⛔ not the increment), and does
   **G2 extract the six blocks from `I_s`** (⛔ not from `model['faces']`)? Cross-check against `audit.py:549`.
6. **G3 adjoint convention PINNED.** Is the convention now literally chosen (bilinear; no conjugation ⇒ `i`/`ω`/`Λ_I(ω)`
   unchanged, no `ω→−ω`; adjoint = swap `extract` trial/test slots + `(X,k_out)↔(Y,k_in)`; face label `s` preserved;
   in-plane IBP boundary 0; `Z`/Fourier not evaluated), ⛔ not left to the builder? Is the pairing correctly c2 §3b /
   S11c-b §3c (⛔ not §3d.4)? Is the one-sided corruption now **named/concrete** (face `s=+`, `u_L→T`/`DIV_U` reverse
   source, relocate `δp_s` slot across sectors + sign flip), with baseline+corrupted for both routes, machine-readable
   provenance, ⛔ no `A−A`? ⚠ Sanity-check the pinned convention is self-consistent (a bilinear no-conjugation adjoint
   for a frequency-dependent `Λ_I(ω)` — is "carried unchanged" the right call, or should the leg flag it?).
7. **Hygiene.** Deleted `reduction/` helpers no longer referenced (generic "use only named machinery")? Header
   (astra-written deliverable → fresh-Claude+Grok legs) coherent?

**No leak / no new defect:** any residual leaked target, manufactured path, or NEW inconsistency the round-3 fold
introduced (e.g. the per-face block, the pinned adjoint convention, or the concrete corruption contradicting §3c or
the audit)?

## Output
For each of 1–7 (+ leak/new-defect): CONFIRMED fixed / still-open, with quoted lines. End: is the directive **now
ready to build**, or the exact remaining decisions/wording that must change.
