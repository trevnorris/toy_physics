# F/G re-grounding — DEFERRED at a tractability wall (2026-09-06)

**Context.** STEP A's physics adjudication resolved F (uniform-limit decoupling) and G (directionality/adjointness)
of the c2 self-energy increment using **orchestrator-authored, now-retired biased** scripts (`verify_F.py`,
`verify_EG.py`). Under the corrected process (orchestrator never authors the CAS instrument; CLAUDE.md `6f8dbd34`)
these were re-grounded the right way. E became the §5c N6 **spec correction** (committed `30d4b72d`). F and G were
carried through the full corrected pipeline; it hit a computational wall.

## What was done (the corrected pipeline, all clean)
1. **Question-vet (Codex-sol):** caught that BOTH my F and G questions were proxies — F needs a **weak-kernel
   identity** (not `.doit()==0`), G directionality is **representation-specific** and the "no adjointness residual"
   disposition needs an **independence audit** (`scratchpad/codex_FG_question_vet.txt`).
2. **Build directive** (`directives/S11c_c2_FG_diagnostic_build_directive.md`), **review-until-clear, 4 rounds**
   (Codex-sol + Grok each round): round 1 ≈24 defects → round 2 (7) → round 3 (few) → round 4 (folded, built). The
   review earned its keep: it caught that my F test pointed at the very 3-symbol `UNIFORM_LIMIT` proxy it was meant
   to replace; the wrong adjoint pairing (§3d.4 vs §3b); the double-map; the un-pinned adjoint convention; etc.
   (leg prompts `directives/_legs/S11c_c2_FG_*`).
3. **astra build** (`scripts/S11c_c2_FG_diagnostic_sympy.py`, 638 lines) — committed here as a REFERENCE only.

## The wall (the finding)
The full-**symbolic** weak-kernel + six-block + adjoint test, over **all retained grades** on the c2 increment
operands (which are enormous — the c2 `.out` is 499 MB), is **computationally impractical**: astra had to iterate
the script ~**78** times; each single-case `TRIAGE` run takes ~**11 min** and emits ~**150 MB** of output (giant
CAS expressions). Even completed, a 150 MB/case deliverable is un-reviewable and un-adjudicable. This is an
**approach** problem, not a directive defect — the build ran clean (empty stderr); it is just far too heavy.

## Disposition — DEFERRED (user decision 2026-09-06)
- F and G are **not blocking**: the self-energy increment **VALUES** (per anchoring) are unaffected — only the F/G
  **interpretations** rest on the retired biased instruments.
- The retired `verify_F.py`/`verify_EG.py` conclusions **do not stand** as re-grounded; they are withdrawn.
- **Correct re-grounding OWED:** a **numeric-probe** diagnostic (Schwartz-Zippel: evaluate the weak-kernels /
  six blocks / pairing residual at random numeric points to test zero/nonzero) — cheap, tiny output. That is a
  separate directive-revision → astra build → build-legs cycle, to be done when F/G are needed (e.g. for the c2
  step record). The reviewed directive above is the symbolic spec; the numeric redesign reuses its object
  definitions (S_P split, per-face pin, pinned adjoint convention, one-sided-corruption requirements) but swaps the
  symbolic weak-zero test for numeric evaluation.
- ⛔ Do NOT run `scripts/S11c_c2_FG_diagnostic_sympy.py` full-symbolic (the wall above); it is committed as a
  reference for the numeric redesign, ⛔ not an accepted/​leg-reviewed deliverable.

## STEP A re-grounding — status
- **E/N6:** RE-GROUNDED → the real finding was a **spec mis-specification** (§5c), corrected + committed `30d4b72d`.
- **F, G:** re-grounding attempted through the full corrected pipeline; **DEFERRED** at the tractability wall;
  numeric re-grounding owed. Increment values unaffected.
