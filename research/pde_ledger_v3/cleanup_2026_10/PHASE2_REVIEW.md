# Phase 2 review (Claude + Grok) and fix instructions

**Reviewed:** commits `a442df4e`, `b81d87a5` and `8e0a8ee1`, against `c8885c4f`.
**Legs:** Claude (orchestrator) and Grok 4.7. Grok's literal report is [phase2_review_grok.txt](phase2_review_grok.txt),
written in response to [phase2_review_prompt.md](phase2_review_prompt.md).
**Verdict: NEEDS REVISION, with targeted fixes only.** The records are faithful on the main points: the benchmark
numbers, the c2 `[0,2]` wording, the clean-condition scope, the pre-repair status of the original b/c2 verification,
and the paper's visible limits. The banked plan content (S_leak identity, dark-energy postulate, MacCullagh
departure) survived the rewrite. Everything below was verified by Claude with the lookups listed.

## Fixes

**F1. The Wolfram N6 pressure-trace repair is unreviewed.**
- Where: the d record, line 50, and the c2 "Changed after close" section, line 57. Both classify it as
  CONDITIONAL alongside the scoped repairs.
- Source: `_measurements/S11c_wolfram_pressure_trace_repair_report.md:78–80` says "No review, comparator,
  recurring scheduler or push has run."
- Fix: say it is **unreviewed**, as the sheet repair is.
- Also state what it changed: it corrected a reference-pressure double shift by extracting the reference-pressure
  map from the face law and solving its ordered inverse (report lines 3–5;
  `_measurements/S11c_wolfram_repair_audit_report.md:16`).

**F2. The c2 thickness-coordinate row cites the wrong commit.**
- `3b52afcb` edits only `scripts/S11c_b_brane_operator_sympy_audit.py`. The c2 kinetic rewrite is `f618178a`,
  which edits `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`. The c2 export of that family is `537d78fd`
  (`_measurements/S11c_thickness_coordinate_repair_report.md` lines 63, 161 and 270).
- Verified by: `git show --name-only 3b52afcb` and `git show --name-only f618178a`.
- Fix: keep `3b52afcb` in the b record. In the c2 record, cite `f618178a` for the c2 change, plus the c2 export
  commit.

**F3. The inertia clearance domain is missing "gradient-free".**
- The Claude verdict covers a **gradient-free test energy** (`_measurements/S11c_upstream_repair_review_disposition.md:24`,
  and lines 100–101: "valid only for that gradient-free test"). The full material-constraint fold was
  source-inspected only.
- Fix: put both limits into the inertia row of the b and c2 tables.

**F4. The d record upgrades the receiving result.**
- Where: the d record, line 58, says "scoped three-field receiving regularity".
- Source: the handoff (`docs/light_em_investigation_handoff.md:42`) and its result record say something narrower:
  scoped receiving regularity **excludes the recorded C3 real-axis poles**; the three-row response is zero in one
  weighted class only; full five-field transverse poles remain.
- Fix: use the source wording, as FINDINGS line 83 already does.

**F5. The front-door summaries drop the uniform check's main content.**
- The result is: finite modes, nonzero current, and **zero face drives** (light does not push on the bulk) at,
  above and below the speed matches. Source: `_measurements/S11c_d_near_unity_uniform_continue_result.md:3–6, 38`.
- STATUS (lines 13–14), the closeout's "Where we stand" (line 7) and the plan's S11c summary (lines 509–510) keep
  only "finite".
- Fix: add the zero face drives in all three, in plain words, and keep the method-only review limit.

**F6. STATUS dropped open debts that the old STATUS carried.**

Grok found one: the **μ_S shear-normalization debt**. It is listed as OPEN in the old STATUS line 17, and later
reports still say it is open (`_measurements/S11c_inertia_repair_report.md:27`;
`_measurements/S11c_d_end_resolvent_report.md:177–178`). Claude found more. Old versus new keyword counts:
- `KEYING`: 3 → 0;
- `FORMALIZATION_POLICY`: 2 → 0;
- `L control`: 4 → 0;
- `DEFERRED_HEAVY`: 6 → 0;
- `remote compute`: 12 → 0;
- `SUBSTRATE`: 7 → 0.

Command: `git show c8885c4f:STATUS.md | grep -c -i <key>` against the new file.

The items behind those counts:
- the review owed for the CLAUDE.md L control and `lean/FORMALIZATION_POLICY.md`. Both landed unreviewed, in the
  mixed-scope commit `c2bdb663`. The old STATUS also asked to confirm the S10 Lean fidelity reviews were by
  non-authors.
- the ≥64 GB cross-engine runs in `DEFERRED_HEAVY_RUNS.md` and the remote-compute plan;
- the S11c-a full-axis keying debt;
- the substrate requirements (`SUBSTRATE_REQUIREMENTS.md`, never run);
- the state of the cleanup itself: phases 3–5 pending, the branch to be rewritten, nothing pushed since 2026-09-30.

**Fix: do a carry-over sweep of the old STATUS (`git show c8885c4f:STATUS.md`).**
- Every open, owed, deferred or debt item gets exactly one disposition:
  - **carried:** one line in the new STATUS, with its owner, or a pointer to the register that carries it
    (`DEFECT_REGISTER.md`, `DEFERRED_HEAVY_RUNS.md`, `SUBSTRATE_REQUIREMENTS.md`);
  - **resolved:** cite the evidence;
  - **superseded:** cite what superseded it.
- Record the sweep as **Part F of FINDINGS.md**, a table of item, old STATUS line, disposition and evidence.
- Keep STATUS short. Group carried items by owner, and prefer one pointer to a register over many lines.
- Add one line on the cleanup state.

## After the fixes

Stop with a short report: each fix and its commit, and the Part F counts (carried, resolved, superseded).
Claude will verify the fixes against this list. Grok will re-check the STATUS carry-over and Part F only, since that
is new content. Phase 3 starts after that.

## Round 2: check of the F1–F6 fixes (`996eee4b`, `19042324`)

**F1–F5: verified fixed** by Claude:
- The d record (line 50) and c2 record (line 57) now mark the N6 repair as unreviewed.
- The c2 row cites `f618178a` and the export `537d78fd`.
- Both b and c2 tables now say "gradient-free".
- The d record (line 58) uses the source's C3 wording.
- STATUS (lines 13–14), the closeout (line 7) and the plan (lines 509–510) now carry the zero face drives.

**F6: Grok re-checked the carry-over** ([phase2b_review_grok.txt](phase2b_review_grok.txt), prompt
[phase2b_review_prompt.md](phase2b_review_prompt.md)). It found two misses, and Claude verified both:

- **G1.** The FINDINGS Part F row at line 201 marks old STATUS lines 84–103 as **resolved**. But old line 89 carries
  the **three N6 premise caveats** forward:
  1. whether Φ is physically correct;
  2. the `V` transform;
  3. extracted-block leakage.

  These are still open in `steps/S11c_c2_self_energy_fold.md:180–184, 203` and in
  `_measurements/S11c_c2_N6_reconcile_disposition.md` §7.
  **Fix:** split the row. The corrected test and the scoped covariance are **resolved**. The three premise caveats
  are **carried**, owned by S11c-c2/d if reused. Add them to the STATUS open table.
- **G2.** Part F (line 195) and STATUS (line 32) send the **full c2 self-energy cross-engine residual** to
  `DEFERRED_HEAVY_RUNS.md`, but that file has no c2 entry: `grep -c -i "self-energy\|S11c-c2" DEFERRED_HEAVY_RUNS.md`
  returns 0.
  **Fix:** add a short entry to `DEFERRED_HEAVY_RUNS.md` stating what must run (≥64 GB, the assembled self-energy
  operator's full cross-engine residual) and pointing at `steps/S11c_c2_self_energy_fold.md:32–34`. Or carry the
  item in STATUS directly. Use one or the other.

**Next:** apply G1 and G2, then go straight into **Phase 3** without stopping in between. Claude will verify G1 and
G2 at the Phase 3 STOP.

**Phase 3b, AGENTS.md:** the replacement text is in [AGENTS_replacement.md](AGENTS_replacement.md). Copy its content
to the root `AGENTS.md`, unchanged and in its own commit, with the title line kept. This is a process policy, not
physics-bearing: it changes no computation, premise, check or claim. So it gets a Claude review at the Phase 3 STOP,
not a two-leg physics review.
