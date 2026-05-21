# Codex — Apply red-team directive for unit {UNIT_ID}

You are applying a red-team directive that was written by an auditor agent. The directive lives at:

  {DIRECTIVE_PATH}

Read it in full. Apply each finding's required change exactly as written. Edit only the files and line ranges the directive names.

## Hard rules

- **Scripts only.** Touch only `.py` and `.wl` script files (and their saved `.txt` outputs if the directive explicitly asks for regeneration). Do NOT touch `paper.tex`, `notes/*.md`, or any other prose document. The red-team is a script-verification loop; doc alignment is handled out-of-band.
- **No new features, no refactors, no stylistic changes.** If the directive says "fix the sign on line 42," fix the sign on line 42. Nothing else.
- **No execution.** Do not run python, mathematica, tests, or any commands beyond what's needed to navigate the repo (cat, grep, ls, etc.). Verification is done in a separate phase, not by you.
- **No guessing.** If a finding's required change is ambiguous, mark it blocked (see below) and skip it. Better to leave it unfixed than to apply the wrong change.
- **Edit in place.** Files are already under git; the orchestrator handles staging/commits.

## Session continuity

This may be a `--resume` of a prior session on the same audit unit. If you recognize the unit ID, prefer applying the *delta* described in the new directive over re-doing earlier work. Earlier `## Applied:` blocks in the directive file tell you what's already been done.

## Per-finding workflow

For each finding `F1, F2, …` in the directive:

1. Read the finding's target file and line range.
2. Apply the required change.
3. Append a block to the directive file under that finding:

```markdown
## Applied: F<n>

- files_changed:
  - `<path>`
- summary: <one sentence describing the edit>
- deviation: <"none" or a one-line explanation of any deviation>
```

If you cannot apply the finding cleanly, append this instead:

```markdown
## Blocked: F<n>

- reason: <what's ambiguous or unsafe>
- question: <a specific question whose answer would unblock you>
```

Continue with the next finding even if one is blocked.

## When the directive includes a missing-script finding

The directive will contain a **claim manifest** listing the physical results the new script must verify. Create the script at the path specified in the directive (or derive it from the existing `paths.sympy` / `paths.mathematica` pattern in `.redteam-config.yaml`).

For SymPy scripts:
- Use `sympy` symbolic primitives. State assumptions explicitly (real, positive, etc.) and only when the physics warrants them.
- For each claim M1, M2, … in the manifest, write a check that re-derives the result from first principles and asserts `simplify(lhs - rhs) == 0` (or equivalent).
- No tautological constructions (`x = expr; assert x == expr`). The assertion must be capable of failing if the math is wrong.
- Print the residual under each check so the saved output shows what's being verified.
- Exit non-zero on failure.

For Mathematica scripts:
- Use `.wl` (script form, not notebook).
- Derive **independently** from physical premises — do NOT transliterate the SymPy script's variable choreography. Different intermediate steps, different function names, different reasoning path. The point is two independent witnesses, not one witness in two languages.
- `Print` each check and its reduced residual.
- Use `Exit[1]` on failure.

Both scripts must be runnable via the runners in `.redteam-config.yaml`.

## When you finish

After all findings have an `## Applied:` or `## Blocked:` block, update the directive's front-matter:

```yaml
applied: true
applied_at: <ISO 8601>
findings_applied: <int>
findings_blocked: <int>
```

That's it. The orchestrator picks up from here — it will run the scripts, capture outputs, and dispatch a verifier agent.
