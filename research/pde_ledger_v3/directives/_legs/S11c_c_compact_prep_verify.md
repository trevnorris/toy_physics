# Verify the S11c-c compact-prep doc/memory updates against the committed reality

I am about to compact this session. Before I do, verify that the status/plan/memory updates I just made are
**accurate and consistent with what is actually committed** — no overclaim, no stale reference, no internal
contradiction. Report every discrepancy; if a file is correct, say so explicitly. This is a correctness audit,
not a rewrite.

## Ground truth — check the docs AGAINST these, do not trust the docs' own claims
- Git: `git -C /var/projects/toy_physics log --oneline -8`. The c1 spec commit should be **`db5cbf88`**
  ("S11c-c START: 2-way split …"). Confirm that SHA exists and that the message matches.
- The committed c1 spec: `/var/projects/toy_projects/../toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md`
  (use `git show db5cbf88:research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md` if the working tree differs).
- The review records: `research/pde_ledger_v3/directives/_measurements/S11c_c_spec_review.md` (the v1 monolith
  review + the 2-way-split decision) and `research/pde_ledger_v3/directives/_measurements/S11c_c1_SHARED_PHYSICS.md`
  (the c1 2-leg review + the union of 10 folds).
- The saved leg derivations under `research/pde_ledger_v3/directives/_legs/S11c_c*` (Codex + Grok scripts + `.out`).

## The updated docs to audit (read each, then check every claim against ground truth)
1. `/var/projects/toy_physics/STATUS.md` — the NEW top clause "S11c-c STARTED — 2-way SPLIT …" (and confirm the
   next clause down, the S11c-b close, is not contradicted). Check: is every SHA, every fold, every "reserved
   for c2" item, and the split description accurate against the committed c1 spec + the review records? Is the
   NEXT (c1 dual-engine build) right? Does it wrongly claim anything is BUILT/committed that is only authored, or
   vice versa?
2. `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md` — the S11c-c
   auto-memory. Check every claim against the committed spec + records: the 2-way split; the c1 folds (the
   two-momentum DtN operator, the 3-object dissipation with the independent traction-vs-far-field-flux route, the
   Σ_E/μ_R,bg-reserved-for-c2, the Hanzawa route, the zero-jet-target, Fredholm-vs-algebraic loci, N11a
   rest-frame); the held-c2 items (extract-then-close non-commutation, θ-row Λ_X/J_s routing, substitution-
   increment). Flag any claim not supported by the committed files. Confirm the YAML frontmatter is valid and the
   description is not truncated.
3. `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/MEMORY.md` — the two S11c one-line hooks
   (S11c-c and S11c-b). Check they are accurate one-liners and their `[[links]]`/paths resolve to real files.

## Specific things to check hard (these are where compact-prep has burned rounds before)
- **Authored vs built vs committed.** The c1 SPEC is committed (`db5cbf88`); NO engine is built yet. Any doc that
  implies the c1 build exists, or that the cross-engine residual ran, is WRONG. Flag it.
- **The fold accuracy.** Does STATUS/memory correctly attribute each fold to c1 vs the held-c2 items? In
  particular: is it correct that (a) the extract/eliminate non-commutation, (b) the θ-row Λ_X routing, and (c)
  the substitution-increment are C2 items (not yet folded into any spec — c2 is unauthored)? And that the 10
  applied folds are all in the committed c1 spec?
- **Stale S11c-b clauses.** STATUS still carries older in-flight S11c-b clauses below the close (STEP-0-FITS-box,
  the integration-pass plan). Do any of them now directly contradict the S11c-b close or the S11c-c start in a
  way a reader would be misled by? (The STEP-0-FITS clause was OVERTURNED — is that clear from ordering, or does
  it read as live?)
- **Any SHA that does not resolve**, any file path that does not exist, any `[[memory-link]]` whose target file
  is absent.

## Output
A short list: for each of the 3 audited docs, either "ACCURATE" or the specific discrepancies (quote the doc line
+ the ground-truth it contradicts). End with a one-line verdict: are we in the clear to compact, or are there
must-fix discrepancies first? Do NOT rewrite the files; just report.
