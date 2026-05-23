# Codex — Apply red-team directive for unit {UNIT_ID}

You are applying a red-team directive that was written by an auditor agent. The directive lives at:

  {DIRECTIVE_PATH}

Read it in full. Apply each finding's required change exactly as written. Edit only the files and line ranges the directive names.

## Hard rules

- **Scripts only.** Touch only `.py` and `.wl` script files (and their saved `.txt` outputs if the directive explicitly asks for regeneration). Do NOT touch `paper.tex`, `notes/*.md`, or any other prose document. The red-team is a script-verification loop; doc alignment is handled out-of-band.
- **No new features, no refactors, no stylistic changes.** If the directive says "fix the sign on line 42," fix the sign on line 42. Nothing else.
- **Execute to validate, and iterate.** After you make edits, RUN the affected scripts and confirm they exit 0. Use `python3 <script-path>` for SymPy and `math -script <script-path>` for Mathematica. If the script fails (non-zero exit, error message, FAIL on an `expectZero`/`assert` check), read the output, diagnose the cause, fix the script, and run again. Do not stop iterating until every script you edited (or created) exits 0 with all its in-file checks passing. The orchestrator dispatches a separate clean-context verifier agent AFTER you confirm the scripts run — that verifier reviews substance, not whether the script runs. Getting the script to run is your job.
- **Iteration cap.** If the same failure persists after ~5 attempts, mark the finding `Blocked` rather than thrash. Better to escalate a hard case than burn the session.
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

The directive will contain a **claim manifest** listing the physical results the new script must verify. Create the script at the path specified in the directive. If the directive doesn't name the path explicitly, derive it from the existing `paths.sympy` / `paths.mathematica` patterns in `.redteam-config.yaml` — DO NOT just put a `.wl` file next to the `.py` file. The conventional layout is:

- **SymPy** `.py` files live in `scripts/`. The matching saved output goes in `scripts/output/`.
- **Mathematica** `.wl` files live in `mathematica/` (NOT `scripts/`). The matching saved output goes in `mathematica/output/`.

So for a missing_mathematica finding on unit `NNN_topic`, the new file path is `mathematica/moving_throat_pde_stage<NNN>_<topic>_mathematica_audit.wl` — placed in the `mathematica/` directory at the project root, NOT in `scripts/`. Read `.redteam-config.yaml`'s `paths.mathematica` glob if uncertain.

For SymPy scripts:
- Use `sympy` symbolic primitives. State assumptions explicitly (real, positive, etc.) and only when the physics warrants them.
- For each claim M1, M2, … in the manifest, write a check that re-derives the result from first principles and asserts `simplify(lhs - rhs) == 0` (or equivalent).
- No tautological constructions (`x = expr; assert x == expr`). The assertion must be capable of failing if the math is wrong.
- Print the residual under each check so the saved output shows what's being verified.
- Exit non-zero on failure.

For Mathematica scripts:
- Use `.wl` (script form, not notebook). PLACED IN `mathematica/`, NOT `scripts/`.
- Derive **independently** from physical premises — do NOT transliterate the SymPy script's variable choreography. Different intermediate steps, different function names, different reasoning path. The point is two independent witnesses, not one witness in two languages.
- `Print` each check and its reduced residual.
- Use `Exit[1]` on failure.

Both scripts must be runnable via the runners in `.redteam-config.yaml`. Run them yourself (`python3 <path>` / `math -script <path>`) and iterate until they exit 0 with all checks passing.

## Common Mathematica pitfalls (defects seen in prior batches)

These have each cost a fix-loop iteration at least once. Check for them proactively when authoring or editing a `.wl`:

1. **Multi-line continuation truncation on `=`.** A statement like

   ```
   lRed =
     1/2 m D[q, t]^2 - 1/2 k q^2
     + 1/2 D[a, t]^2 - 1/2 oA^2 a^2
     + r a ww + gA q a + gW q ww;
   ```

   parses as `lRed = (1/2 m D[q,t]^2 - 1/2 k q^2)` — Mathematica's parser ends the expression at the first line whose trailing token completes a parseable form, silently dropping subsequent `+ ...` lines. **Always parenthesize the RHS** when an assignment spans multiple lines:

   ```
   lRed = (
     1/2 m D[q, t]^2 - 1/2 k q^2
     + 1/2 D[a, t]^2 - 1/2 oA^2 a^2
     + r a ww + gA q a + gW q ww
   );
   ```

   This defect has been caught at stage 003 (batch I.1) and stage 021 (batch I.2). If you are touching any multi-line assignment whose RHS uses leading-`+`/`-` continuation, parenthesize it.

2. **`D[expr, f[t]]` returns 0 when `f[t]` is a head-applied form, not a symbol.** When local variables are bound to function applications:

   ```
   q = qFun[t]; a = aFun[t]; ww = wFun[t];
   lRed = 1/2 m D[q,t]^2 - 1/2 k q^2 + ... + gA q a + gW q ww;
   D[lRed, q]   (* returns -k qFun[t], MISSING the gA aFun[t] + gW wFun[t] terms *)
   ```

   Mathematica's `D` does not differentiate w.r.t. `qFun[t]` for product terms — it treats the head-applied form as not-a-variable for those terms. So a manual Euler-Lagrange operator pattern `D[D[lRed, D[q,t]], t] - D[lRed, q]` silently yields a wrong EOM. The script exits 0, the `expectZero` Prints "FAIL: ...", and you don't notice unless you read the output line.

   **Use `EulerEquations` from the `VariationalMethods` package** instead:

   ```
   Needs["VariationalMethods`"];
   elList = EulerEquations[lRed, {qFun[t], aFun[t], wFun[t]}, t];
   ```

   `EulerEquations` returns a list of `lhs == rhs` equations; coerce to a residual via `/. Equal[lhs_, rhs_] :> lhs - rhs` before passing to `expectZero`. Note that `EulerEquations` may return the opposite sign convention; if so, flip the canonical RHS rather than reintroducing the manual operator. Fallback: substitute true symbols for the head-form aliases before calling `D`.

   This defect has been caught at stage 021 (batch I.2).

3. **`expectZero` Prints "FAIL" but does not `Exit[1]`.** Many existing `.wl` files use a helper that Prints failure but doesn't propagate. The orchestrator's post-fix sanity check now greps the saved output for `^FAIL\b` so this can't pass silently — but if you are writing a new helper, **also call `Exit[1]` on failure**.

4. **`Solve`/`Reduce` returning `ConditionalExpression[0, cond]`, which is not `=== 0`.** When you replace a hand-typed expression with `var /. First[Solve[eqn, var]]`, `Solve[..., Reals]`, or `Reduce[..., var]`, the result is often wrapped in `ConditionalExpression[..., <inequality on declared parameters>]`. Substituting that back into a residual produces `ConditionalExpression[0, cond]`, which is identically zero under `$Assumptions` but does not match the standard `If[TrueQ[res === 0], ...]` check inside `expectZero` — so the residual prints `ConditionalExpression[0, cond]`, the helper Prints "FAIL", and the script exits 1.

   **Strip the wrapper before passing the value forward**, or extend the helper:

   ```
   xEqRaw = xi /. First[Solve[{eqn, xi > 0}, xi, Reals]];
   xEqRaw = xEqRaw /. ConditionalExpression[e_, _] :> e;
   xEq    = FullSimplify[xEqRaw, Assumptions -> $Assumptions];
   ```

   Or, if you control the helper, build the strip in once:

   ```
   expectZero[name_String, expr_] := Module[{res},
     res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
     res = res /. ConditionalExpression[e_, _] :> e;
     res = FullSimplify[res, Assumptions -> $Assumptions];
     Print[name, " = ", fmt[res]];
     If[TrueQ[res === 0], pass[name], fail[name, res]];
   ];
   ```

   Substance-preserving — `ConditionalExpression[0, cond]` is genuinely zero on the declared domain. This defect has been caught at stages 050, 051, and 052 (batch III.2).

5. **`Limit` non-determinism for poles: `Infinity` vs `Infinity/<positive>`.** `Limit[fraction_with_simple_pole, var -> 1, Direction -> "FromBelow"]` returns either `Infinity` or `Infinity/(positive polynomial in declared parameters)` non-deterministically across `math -script` invocations. Strict checks like `If[pi1 =!= Infinity, fail[...]]` are flaky.

   **Test the inverse vanishes instead:**

   ```
   pi1 = Limit[piTr, xi -> 1, Direction -> "FromBelow"];
   If[!TrueQ[Simplify[1/pi1 == 0, Assumptions -> $Assumptions]],
      fail["Pi_tr(xi->1-) is not +infinity", pi1], None];
   ```

   `1/Infinity == 0` and `1/(Infinity/positive) == 0` both reduce to `True`, so this accepts both surface forms while still failing for finite limits or `Indeterminate`. This defect has been caught at stage 051 (batch III.2).

## When you finish

Before marking the directive applied, confirm: every `.py` and `.wl` script you edited or created exits 0 when run, with all in-file `assert` / `expectZero` / `If[..., Exit[1]]` checks passing. If you cannot reach that state for a finding within ~5 iterations, append `## Blocked: F<n>` instead of `## Applied: F<n>`.

After all findings have an `## Applied:` or `## Blocked:` block, update the directive's front-matter:

```yaml
applied: true
applied_at: <ISO 8601>
findings_applied: <int>
findings_blocked: <int>
```

That's it. The orchestrator picks up from here — it dispatches a clean-context verifier agent to assess whether your edits actually closed the *substance* of each finding (not whether they run; that part you already confirmed).
