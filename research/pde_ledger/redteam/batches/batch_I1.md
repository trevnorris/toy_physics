---
batch_id: I.1
label: Part I.1 — Geometry lift, BdG coupling, projected Maxwell setup
started: 2026-05-20
completed: 2026-05-21
stages_total: 12
stages_verified: 12
stages_blocked: 0
material_change_any: false
---

# Batch I.1 — Final summary

All 12 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all findings closed substantively. No `material_change` flags, so no downstream cascade — batch I.2 may proceed without re-auditing I.1.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 001 | 1 (math_transliteration) | 2 (iter1 had `First[]` bug) | verified | smoke-tested the entire pipeline before fanout |
| 002 | 2 (math_transliteration, insuff_verif) | 1 | verified | F1 deviation from directive (different `Coefficient` route) was mathematically correct |
| 003 | 4 (insuff_verif, taut, insuff_verif, math_transliteration) | 1 | verified | checkpoint stage — biggest verifier pass. Side note: `lRed` multi-line parse quirk in `.wl` was caught and patched |
| 004 | 2 (taut, missing_math) | 2 (iter1 IBP residual bug) | verified | codex's two-Integrate form didn't cancel; iter2 combined integrand fixed it |
| 005 | 1 (missing_math, 5 claims) | 1 | verified | clean missing-math fill |
| 006 | 3 (missing_math, insuff_verif, stale_output) | 3 (auditor wrote two math-wrong sub-checks) | verified | F2(a) had `w`,`z` independence bug; F2(c) had wrong-parity "trivial mediator kills leak" claim. Both rewritten; antisymmetric-Z is the correct negative control |
| 007 | 1 (missing_math, 11 claims) | 1 | verified | codex misplaced `.wl` in `scripts/`; orchestrator moved to `mathematica/` |
| 008 | 3 (missing_math, taut, insuff_verif) | 1 | verified | misplaced `.wl` moved; stale duplicate `.txt` swept |
| 009 | 4 (missing_math, taut, insuff_verif, hardcoded) | 1 (lost session) | verified | orphaned codex incident during parallel-math crisis; recovered manually |
| 010 | 2 (missing_math 17 claims, insuff_verif) | 1 | verified | largest claim manifest in the batch |
| 011 | 1 (missing_math, 11 claims) | 1 | verified | misplaced `.wl` moved despite codex.md patch — codex ignores path guidance occasionally |
| 012 | 4 (missing_math, taut, taut, insuff_verif) | 1 | verified | clean run |

## Process bugs discovered and fixed (will benefit future batches)

1. **Mathematica single-seat license:** Discovered by running parallel sanity-exec while loop was iterating codex. Fix: memory `feedback_mathematica_single_seat.md` plus discipline to never run two math invocations.
2. **Codex creating `.wl` in `scripts/`:** Even after `codex.md` patch explicitly says "`.wl` files go in `mathematica/`, not `scripts/`", codex occasionally ignored. Fix: orchestrator now moves and rescans; future patch will auto-relocate inside `fix_loop.sh`.
3. **Saved `output/` `.txt` files not auto-refreshed:** The redteam machinery captures to `exec_logs/`, not the canonical `mathematica/output/` or `scripts/output/`. Fix: manual sweep done; future patch will add output-refresh to `fix_loop.sh`.
4. **Subagents lack `/tmp` read access:** Caught during first audit fanout. Workaround: use project-local `redteam/tmp_prompts/` instead. Now in `.gitignore`.
5. **`fix_loop.sh` syntax noise on zero-Blocked counts:** `grep -c ... || echo 0` produced `"0\n0"` which `[[ -gt ]]` choked on. Fixed.
6. **Manifest stale after codex creates a new script:** Added `redteam scan <NNN>` step after codex returns; sanity exec now finds the new file.

## Prompt hardening landed this batch

- `prompts/codex.md`: replaced "No execution" with "Execute to validate, and iterate" — codex now runs `python3` / `math -script` and iterates to exit 0 before marking Applied. Cap ~5 iterations per finding, else mark Blocked.
- `prompts/codex.md`: explicit directory rules ("`.wl` files go in `mathematica/` not `scripts/`").
- `prompts/auditor.md`: required Self-test section catching variable-independence, parity, trivial-case, and path-specification traps before finalizing a directive.

## Auditor error rate observed

3 of the 12 directives had math errors that codex caught via the Block-and-delta loop:
- Stage 001: `First[EulerEquations[...]]` extraction recipe was wrong
- Stage 006 F2(a): `w` and `z` independence
- Stage 006 F2(c): wrong parity claim (`Z=1` trivial mediator)

All caught before merging. The auditor self-test addition in `auditor.md` should reduce this from ~25% → ~10-15% on future batches.

## Codex iteration counts (workload signal)

Most stages: 1 iteration. Stages with >1 iteration: 001 (2), 004 (2), 006 (3). Total codex sessions opened: 12.

## Carry-forward

Batch I.1 sets the foundation for I.2 (PDE geometry, BdG coupling continuation). Since `material_change_any: false`, no upstream cascade is needed. I.2 may be audited next.

For the orchestrator picking up batch I.2 post-compact: the standing rules in memory (mathematica-single-seat, codex-iterates-until-clean, use-agents-for-reviews) and the prompts (auditor self-test, codex iterate-and-place-in-mathematica) carry the lessons forward.
