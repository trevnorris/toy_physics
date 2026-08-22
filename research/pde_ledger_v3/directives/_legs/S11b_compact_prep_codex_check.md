# Codex consistency check — S11b compact-prep (can we correctly resume?)

We are about to compact the session. Verify the persisted state is CONSISTENT and that a resume from it would
correctly continue the work — flag any staleness, contradiction, or missing pointer that would mislead a
fresh session. Be adversarial; a rubber-stamp is useless.

## Read (all under /var/projects/toy_physics unless noted)
- `STATUS.md` — the top "CURRENT FRONT — S11b" block (Artifacts #1–#7 + the NEXT block).
- Out-of-repo memory: `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11b_interface_law_result.md`
  (the frontmatter `description` + the top unified-rewrite section).
- `research/pde_ledger_v3/steps/S11b_RUN_CHECKLIST.md` (the live-status sections at the bottom).
- `research/pde_ledger_v3/steps/S11b_wl_engine_review_disposition.md` (the findings + the Codex-consult fold).
- `git log --oneline -20` and the artifacts they name:
  - WL engine `mathematica/S11b_interface_coupling_law_mathematica_audit.wl` + its `.out` under `mathematica/out/`
  - T7 comparator `scripts/S11b_cross_engine_comparator.py` + `scripts/out/S11b_cross_engine_comparison.out`
  - SymPy engine `scripts/S11b_interface_coupling_law_sympy_audit.py` (the X-1 target ~L470-510: `st_squared`
    basis[1]).

## Verify (report each: CONSISTENT / STALE / CONTRADICTION, with the file:line and the fix)
1. **State agreement.** Do STATUS, the memory, and the run-checklist agree on: WL engine DONE, T7 comparator
   DONE+FROZEN+RUN, WL repair DONE (`bd598ae7`, F-WL-1 static `−Λ_A⁰/ρ_m`, all checks genuine), and NEXT =
   the SymPy X-1 repair? Any doc still saying an earlier "NEXT" (e.g. "NEXT = blind Wolfram engine", "NEXT =
   comparator", or "the two SEPARATE repairs")?
2. **X-1 target is actionable.** Does the SymPy engine actually carry `st_squared` as basis[1] at ~L470-510,
   and is the claim ("11 vs correct 10; `st_squared = (2/3)div²+(1/2)curl²` mod a total divergence; §5
   mandates the quotient; the coeff is degenerate with `B_div`/`mu_R`; EOM unaffected") supported by the code
   + the committed comparison output (`ENERGY_BASIS_COUNT` PY 11 vs WL 10)? Is the repair well-enough
   specified for a fresh session to execute it?
3. **WL repair really committed + clean.** Is the WL `.wl` at `bd598ae7` the repaired one (F-WL-1 static, no
   `classifyLeadingBlock[stratumBlock]`, no dead `thresholdUnrelatedClassification`/`unrelatedCausalityAggregation`),
   is its `.out` regenerated (no WolframScript config-warning lines), and is the working tree clean?
4. **Anything a resume would get WRONG** — a claimed-done artifact that is missing/uncommitted, a rule-5 leak
   left in a directive, a pointer to a deleted path, or a contradiction between the docs and the code.

## Output — under ~25 lines
Numbered per check: verdict + file:line + concrete fix. End with a one-line verdict: SAFE TO COMPACT, or fix
these first.
