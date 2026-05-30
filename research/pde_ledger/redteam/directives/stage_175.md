---
unit_id: 175
batch: V.1
created_at: 2026-05-29T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T23:44:35-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 175

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond the edits named. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

---

## F1 — mathematica_transliteration

**Target (Mathematica ONLY):**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl`
  - the `dlog` helper, lines 26-29:
    ```mathematica
    dlog[expr_] := FullSimplify[
      D[Log[FullSimplify[expr, Assumptions -> $Assumptions]], eps] /. eps -> 0,
      Assumptions -> $Assumptions
    ];
    ```
  - the Sigma_N differential block, lines 95-98:
    ```mathematica
    exprPoverDeltaPhys = ((p/delta) /. subsHat) /. subsEps;
    sigmaNDirect = FullSimplify[2*dlog[exprPoverDeltaPhys] - kappa, Assumptions -> $Assumptions];
    sigmaNShape = FullSimplify[dlog[(lambda^2/k) /. subsEps], Assumptions -> $Assumptions];
    expectZero["Sigma_N - dln(Lambda^2/K)", sigmaNDirect - sigmaNShape];
    ```

The SymPy file
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py`
is the REFERENCE engine. Do NOT edit the `.py` at all — its `dlog` route (`sp.diff(sp.log(...))`, lines 100-101; Sigma_N block lines 132-135) is the existing baseline that the new Mathematica route must be independent OF.

**Issue:** The `.wl` Sigma_N differential block is a line-by-line transliteration of the `.py` block. Both engines extract the first-order log-slope via the SAME `dlog` primitive — `D[Log[.], eps] /. eps -> 0` in Mathematica, `sp.diff(sp.log(.), eps).subs(eps, 0)` in SymPy. The differential slope identity `2 dln(P/Delta) - dK = dln(Lambda^2/K)` is therefore SINGLY-ROUTED across the two engines: a wrong target copied from one script into the other would still be likely to pass in both, so the second engine gives no independent slope coverage.

The other Sigma_N coverage does NOT close this gap (per the consult):
- `N0 - Lambda^2` (wl:61 / py:83) is genuinely independent, but it is the STATIC factorization — it does not extract the first-order slope.
- the downstream `Common-shape branch Sigma_N + dK` (wl:113 / py:156) and `Xi_load (all shapes frozen) + dK` (wl:118 / py:164) checks are specializations of `sigmaNDirect` — they inherit the `dlog` slope route rather than adding independence to it.

This is the confirmed structural-independence gap of red-team finding R1 (`codex_reviews/stage_175.md`), and the batch-8 Claude+Codex consult (`codex_reviews/_consult_batch8.md`) CONCURred 4/4 that the slope needs an independent second route.

**Required change (the PRIMARY edit — option B SUPPLEMENT):**

1. Define a NEW, structurally-independent slope extractor near the existing `dlog` helper (immediately after line 29). The load-bearing operation MUST be `Series` + `Coefficient`, NOT `D[Log[...]]`:
   ```mathematica
   dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps];
   ```
   It MAY wrap the series argument with `FullSimplify[..., Assumptions -> $Assumptions]` analogous to `dlog` if needed to land `=== 0`, but the differentiation MUST be the series-coefficient route, not `D[Log[.]]`.

2. ADD exactly ONE new check immediately after the existing `expectZero["Sigma_N - dln(Lambda^2/K)", sigmaNDirect - sigmaNShape];` line (currently line 98), leaving that existing `dlog`-based line UNTOUCHED:
   ```mathematica
   (* Independent second-engine slope route (red-team R1): extract the Sigma_N
      first-order log-slope via Series-coefficient (Series+Coefficient) instead of
      D[Log[.]], so the differential identity no longer relies on a transliteration
      of the SymPy dlog route. Series-route DIRECT (2 dln(P/Delta) - dK) vs the SHAPE
      target dln(Lambda^2/K); -kappa (= -delta_K) kept symbolic. *)
   sigmaNDirectSeries = FullSimplify[2*dlogSeries[exprPoverDeltaPhys] - kappa, Assumptions -> $Assumptions];
   sigmaNShapeSeries  = FullSimplify[dlogSeries[(lambda^2/k) /. subsEps], Assumptions -> $Assumptions];
   expectZero["Sigma_N - dln(Lambda^2/K) [series route]", sigmaNDirectSeries - sigmaNShapeSeries];
   ```
   The variable names (`exprPoverDeltaPhys`, `lambda`, `k`, `subsEps`, `kappa`, `$Assumptions`) are the exact live names in the `.wl` — reuse them as-is; do not rename or redefine them.

**ANTI-GUARDS (do not drift):**
- The new check MUST compare series-route DIRECT (`2*dlogSeries[exprPoverDeltaPhys] - kappa`) vs the SHAPE target (`dlogSeries[(lambda^2/k) /. subsEps]`). It MUST NOT be `dlogSeries[e] - dlog[e]` on the SAME argument — that only checks two extraction methods against each other (a differentiation-method tautology), NOT the Sigma_N physics/factorization identity.
- `-kappa` (= `-delta_K`) MUST stay symbolic; never substitute a numeric value for it.
- Do NOT replace or redefine `sigmaNDirect` or `sigmaNShape`, and do NOT touch the downstream `sigmaNCommon` / `xiLoadFrozen` / `Xi_load` lines (wl:110-118). The supplement leaves the existing `dlog` route in place as corroboration. This is the consult-chosen smaller-blast-radius option B; option A ("replace `sigmaNDirect`/`sigmaNShape` with `dlogSeries`") was explicitly REJECTED for unnecessary blast radius (it would expose the downstream `sigmaNCommon` / `Xi_load` checks to the feasibility risk too).
- Do NOT touch the SymPy `.py` at all. It is the reference engine; its transcript must come out UNCHANGED.
- Do NOT re-touch or re-prescribe the existing `Sigma_N - dln(Lambda^2/K)` `dlog` line or the F1/F2 lines that already pass (see "What this directive does NOT change" below).

**MANDATORY ESCAPE CLAUSE:** If, after running Mathematica, the `dlogSeries` route does NOT robustly reduce the `[series route]` residual to `=== 0` under `FullSimplify` — and the failure is a CAS-normalization / simplification limitation rather than a demonstrated algebraic discrepancy — do NOT force it and do NOT let a previously-passing audit silently change. Instead append `## Blocked: F1` recording the EXACT residual surface form and noting that the series route needs the orchestrator to record a sanctioned MIRROR_POLICY exception (the existing `dlog` block is then accepted as a policy mirror, with R1's structural-independence requirement WAIVED-with-justification, not satisfied — recorded honestly as "waived," NOT "independent series-route coverage achieved"). Apply nothing else in that case.

**Verification command:**
The orchestrator runs `redteam exec-mathematica 175` and `redteam exec-sympy 175`. A new line `Sigma_N - dln(Lambda^2/K) [series route] = 0` must appear in the Mathematica transcript, and both scripts must exit 0. The SymPy transcript must be UNCHANGED (the `.py` is not edited).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl`
- summary: Added an independent Mathematica Series+Coefficient slope route for Sigma_N and checked it against dln(Lambda^2/K).
- deviation: none
