# Decision RE-REVIEW (round 2) — the F/G diagnostic BUILD DIRECTIVE. STATIC; ⛔ run nothing, load no `.out`.

## ⚠ Paths
Working dir = `/var/projects/toy_physics` (repo root). Every `directives/…`, `scripts/…`, `_measurements/…` path
is under `research/pde_ledger_v3/`. Artifact under review:
`research/pde_ledger_v3/directives/S11c_c2_FG_diagnostic_build_directive.md`. ⛔ Review THAT F/G diagnostic
directive — ⛔ not any WL-engine / T7 / export-chain directive. This is a **physics-bearing** directive, so it is
**review-until-clear**, ⛔ not a one-pass decision list.

You are ONE of two independent legs. **Round 2**: round 1 (Codex-sol + Grok) returned NOT-READY with ~two dozen
defects; the orchestrator folded them. Verify each is now correctly + completely fixed and introduced no new
defect. ⛔ No fictional-script ablation (defer executable tests to the build review). A real problem is a finding.

## The round-1 defects that were folded — verify each in the current directive
**F:**
1. The uniform map must be **enumerated at the imported operands** (W_bg/μ_R,bg/e_W,bg→consts; ALL their spatial
   jets + w1/m1 profile jets→0; both density representatives via `background_density_map`; density gradient→0;
   η=σ_W=0; η-tied face tilt→0), applied **before** close/extract — and it must ⛔ NOT reuse the audit
   `UNIFORM_LIMIT` override (`audit.py:1081` = the retired 3-symbol proxy leaving `W_bg_d*` live, an M3 freeze), and
   ⛔ NOT be `Z→0`/`Λ→0`/`μ_θ→0`. Is F2 correct and complete?
2. The reconstruction must use an **independent `I_direct = difference(closed_kernel,open_kernel)`**, so `R_split`
   is not an `A−A` tautology (F1). Confirmed?
3. Weak-kernel predicate: inspect `Integral.function` on **generic** trial/test; IBP **only** S11c-b §3c in-plane
   (boundary=0); ⛔ NO IBP/`.doit()` of Z/Fourier integrals (F3). Does this actually exclude `.doit()==0`?
4. F emit **per-face AND assembled** (F4); a **slot-inventory** control (I_closed free of P; I_bare shows the slots)
   (F5); the **sentinel injected into the I_closed path** with a polynomial (no-`Integral`) kernel, same checker
   (F6); FORM ablation = **sector-mixing slot relocation + sign**, ⛔ not a face swap (F7); F flagged as a secondary
   smoke test (F8). All present and correct?
5. `ORDERING_COMMUTATOR` reused via the audit path **including `regression_coordinates`**, kept distinct from
   `I_closed`/`CLOSED_COUPLING_KERNEL`. Confirmed?

**G:**
6. Adjoint pairing is **c2 §3b / S11c-b §3c** stored-energy/kinetic weak variational pairing — ⛔ NOT §3d.4 (G3).
   The **formal-adjoint convention is pinned** (bilinear/sesquilinear; i, ω, memory/outgoing branches, X↔Y, momenta,
   face labels, IBP signs)? The **one-sided corruption is concrete** (dimension/grade-matched tagged FORM mutation
   in the reverse mechanical row before extraction; baseline+corrupted for both routes; machine-readable
   provenance)? ⛔ no `A−A`, ⛔ "not implemented ≠ no route".
7. **Witness table for ALL SIX blocks** (label, coefficient, cleared numerator, denom/domain); empty tables are
   data; the script does NOT self-classify live/blocked (G1). Blocks named to the `extract` keys (DIV_U↔u_L)?
   Per-face + assembled (G2)? Eulerian-only, cross-rep deferred (G4)? Dissipativity/energy-sign forbidden?
8. The G FORM ablation is the **source-level one-sided corruption** (the trial/test **sector-collapse was deleted**
   as a special-field proxy). Confirmed removed?

**Buildability + hygiene:**
9. `close` pinned to `build_case`'s substitution map (no `close()` fn); the `model['faces']` caveat (it is the
   closed carrier, not the per-face increment — define per-face carriers + additivity residual). Buildable from the
   named machinery without a nonexistent object?
10. The nonexistent `reduction/engine_output_checks.py` / `derived_or_declared.py` references are **removed** (they
    were deleted `fb29bba2`). Header no longer says "(Codex-written)"; the malformed path line fixed.

**No leak:** ablations/corruptions emit baseline+mutated+diff and ⛔ never prescribe "must change"/a target;
withheld: whether F weak-vanishes, which G blocks live, whether an adjoint route exists. Any residual leaked target
or manufactured path? Any NEW defect or inconsistency the fold introduced?

## Output
For each of 1–10 (+ leak/new-defect): CONFIRMED fixed / still-open, with quoted lines. End: is the directive **now
ready to build**, or the exact remaining decisions/wording that must change.
