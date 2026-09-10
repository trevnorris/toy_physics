# S11c-d SymPy build — STAGE B (part 2): FINISH §2. Focused resume; extend the existing script.

⭐ Staged build (plan: `_measurements/S11c_d_sympy_build_staging_plan.md`). Authority: cleared directive
`directives/S11c_d_sympy_build_directive.md` + spec `399a8516`. Model `gpt-6-astra`. All directive obligations bind
(three script clauses; input-driven; PRINT never assert; builder lane build → run → report → **stop**, ⛔ no
legs/comparator/WL/downstream/commits).

## Read and EXTEND `scripts/S11c_d_mixing_scattering_sympy_audit.py`
Stage A + Stage-B-part-1 are DONE and must be preserved: import/reduction/metadata; `ReducedPencil` sector blocks +
off-diagonal-extraction residual; `ConstantEndPencil` background/abel-limit/phase/integral-symbol/strong+weak
matrices/translate-kernel/curl-gauge operands; computed `K₀,K₋,K₊`; reference/end symbols. The current
`OUTSTANDING_CONSTRUCTIONS` tag lists what is left.

## FINISH §2 this pass (⛔ compute, never type; ⛔ no empty-set/zero substituted for an unexecuted solve)
1. **Full-pencil modes** of each full asymptotic block pencil `𝓛₋^full`, `𝓛₊^full`: incoming/outgoing/evanescent/
   threshold modes (right modes `r_a`, adjoint/left modes `l_a`), with the outgoing/physical-sheet prescription.
2. **Transverse-like / thickness-like classifiers**: computed spectral projectors (or continuous continuation from
   the reference sector basis); emit the classification map; a degeneracy that makes it ill-defined → reported as a
   domain limitation, ⛔ not silently replaced by bare-sector labels.
3. **The transverse gauge quotient** the mode solve needs (the `curl_gauge_operands` groundwork exists — complete it).
4. **Both-row 3-D↔1-D reconstruction round-trips** (close the remaining reduction-metadata gaps the report lists).

Every emitted object carries computed `(ε,η,σ_W)`/`λ` order **and** restored `[L,T,M]` dimension.

## Heavy-CAS + EFFICIENCY (the 7-min full re-run per iteration is the bottleneck — fix it)
- **Iterate cheaply:** the Stage-A reduction is stable. During development, run a **single `(α,ρ)` case** (e.g.
  `LAB_HELD × RHO4_CONSTANT`) to iterate the mode/classifier solve fast; do the **full 4-case run once at the end** of
  the pass to regenerate the complete `.out`. Add a small case/section selector (argv or env) if that helps — ⛔ but
  the final `.out` must cover all four cases and every implemented section.
- **Heavy modes → carrier-first numeric-PIT fingerprint + SHA digest**, ⛔ not a full-symbolic dump (keep the `.out`
  bounded — it is currently ~22 MB, good; keep it so).
- ⛔ Never two memory-heavy CAS jobs concurrently.
- **Do not over-iterate:** implement the four items, run, emit each object + its residual, and once §2 runs
  end-to-end, STOP. Prefer completing all of §2 over polishing one piece — the review legs + orchestrator catch
  issues. ⛔ Do not drift into §3a/§3b/§3c/§3d/§5/export this pass.

## Run + report + stop
Full 4-case run at the end; update the `.out` and `_measurements/S11c_d_sympy_builder_report.md` (mark §2 complete or
name exactly what of §2 remains; per-object which-line-computed-it; §2 residuals literal). ⚠ If §2 still cannot finish
this pass, STOP at a runnable checkpoint and report precisely what remains — ⛔ no fake/empty-substitute. Then STOP.
