---
unit_id: 050
batch: III.2
created_at: 2026-05-26T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-26T02:57:03-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false  # Q1=(a) applied — see Applied: F2 block below
---

# Codex directive — unit 050

Apply each non-`paper_misalignment` finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit `paper.tex`, `notes/`, or scripts to "fix" the paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch `paper.tex`, `notes/`, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:59-95`

**Issue:**
The Mathematica script reuses the SymPy script's hand-written target expressions byte-for-byte in three load-bearing `expectZero` calls: the derivative form-match (wl:60-64 vs sympy:63-68), the `xEq` closed-form match (wl:56,68-70 vs sympy:74), and the ceiling factored-form match (wl:89-91,92-95 vs sympy:100-107). This violates the second-engine independence policy: a typo or sign error in any of the three SymPy targets would be silently confirmed by the Mathematica script subtracting an identical typo.

**Required change:**
Modify three blocks in `moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl` to use independently-derived reference expressions instead of byte-equivalent hand-typed targets.

(a) Derivative check — replace wl:59-64 (currently `dZetaNdx ... dZetaNdxTarget ... expectZero[..., dZetaNdx - dZetaNdxTarget]`).

Before (wl:59-64):
```
dZetaNdx = FullSimplify[D[zetaN, x], Assumptions -> $Assumptions];
dZetaNdxTarget = -n (n + 1) / ((2 n + 1)^2 (1 + n (n + 1) x)^2);
expectZero[
  "d zeta_n^(twin) / dx + n(n+1)/[(2n+1)^2 (1 + x n(n+1))^2]",
  dZetaNdx - dZetaNdxTarget
];
```

After:
```
dZetaNdx = FullSimplify[D[zetaN, x], Assumptions -> $Assumptions];
expectZero[
  "d zeta_n / dx (denominator structure) : dZetaNdx (2n+1)^2 (1 + n(n+1) x)^2 + n(n+1) = 0",
  dZetaNdx (2 n + 1)^2 (1 + n (n + 1) x)^2 + n (n + 1)
];
```

This anchors the derivative to the denominator structure of `zetaN` itself (multiplying through by the denominator squared and adding `n(n+1)`) rather than to a transliterated target expression. The residual is identically zero by direct differentiation, but the algebraic pathway differs from the SymPy form-subtraction.

(b) `xEq` closed-form match — replace wl:68-70 (currently `expectZero["x_max (from Solve) - ...", xEq - xEqClosedForm]`).

Before (wl:68-70):
```
expectZero[
  "x_max (from Solve) - [1/((2n+1)^2 zeta_req)-1]/[n(n+1)]",
  xEq - xEqClosedForm
];
```

After:
```
expectZero[
  "x_max from Solve satisfies (2n+1)^2 zeta_req (1 + n(n+1) x_max) - 1 = 0",
  (2 n + 1)^2 zetaReq (1 + n (n + 1) xEq) - 1
];
```

This verifies that `xEq` satisfies the defining equation of `x_max` (the same equation `Solve` was given) rather than re-subtracting the hand-typed closed form. The line wl:56 declaration of `xEqClosedForm` can remain so the `Print` statement at wl:65 still references it, or both `xEqClosedForm` and any unused references can be removed if not needed elsewhere; either way, `xEqClosedForm` must NOT appear inside any `expectZero` after this edit.

(c) Ceiling factored-form match — replace wl:87-95 (currently `ceilingDiff ... ceilingDiffTarget ... expectZero[..., ceilingDiff - ceilingDiffTarget]`).

Before (wl:87-95):
```
ceilingDiff = FullSimplify[sNMax - sN, Assumptions -> $Assumptions];
Print["S_n^(max) - S_n^(twin) = ", fmt[ceilingDiff]];
ceilingDiffTarget =
  ((1 - eps) (2 n + 1)^2 n (n + 1) x) /
  (((2 n + 1)^2 - eps) ((2 n + 1)^2 (1 + n (n + 1) x) - eps));
expectZero[
  "S_n^(max) - S_n^(twin) factored form",
  ceilingDiff - ceilingDiffTarget
];
```

After:
```
ceilingDiff = FullSimplify[sNMax - sN, Assumptions -> $Assumptions];
Print["S_n^(max) - S_n^(twin) = ", fmt[ceilingDiff]];
ceilingDiffNumerator = FullSimplify[
  Numerator[Together[ceilingDiff]],
  Assumptions -> $Assumptions
];
expectZero[
  "Numerator of (S_n^(max) - S_n^(twin)) - (1-eps)(2n+1)^2 n(n+1) x",
  ceilingDiffNumerator - (1 - eps) (2 n + 1)^2 n (n + 1) x
];
ceilingDiffDenominator = FullSimplify[
  Denominator[Together[ceilingDiff]],
  Assumptions -> $Assumptions
];
expectZero[
  "Denominator of (S_n^(max) - S_n^(twin)) - ((2n+1)^2 - eps) ((2n+1)^2 (1 + n(n+1) x) - eps)",
  ceilingDiffDenominator - ((2 n + 1)^2 - eps) ((2 n + 1)^2 (1 + n (n + 1) x) - eps)
];
```

This decomposes `ceilingDiff` into numerator and denominator via `Together` + `Numerator`/`Denominator` and checks each piece independently — a derivation pathway not used by the SymPy script.

**Important constraint:** preserve every other line of the Mathematica script as-is (the `Clear`, `$Assumptions`, the `banner`, the doubling and equivalence-criterion blocks at wl:32-50, and the prints / final `Exit[0]`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 050` and confirm the script exits 0 with the three new `expectZero` labels appearing in the output transcript. Diff the new `.wl` against the SymPy `.py` and confirm `dZetaNdxTarget`, `xEqClosedForm`-subtraction, and `ceilingDiffTarget` no longer appear inside `expectZero` arguments.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: replaced the derivative, x_max, and ceiling factored-form Mathematica checks with independently-derived residual checks.
- deviation: none.

## F2 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_050.tex:44` quote: `\stagefield{Output}{Lowest-twin sufficiency \eqref{eq:app-stage050-lowest-success} and higher-harmonic exclusion/softness thresholds.}` — no mention of an enhancement ceiling, and no boxed equation for `S_n^(max)` in the body (stage_050.tex:11-42).
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage050_zeta_threshold_comparison.md:172-174` quote: `S_n^(twin) < S_n^(max) := 1 + (1-eps) / [ (2n+1)^2 - eps ].`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:90` quote: `S_n_max = sp.simplify(1 + (1 - eps) / ((2 * n + 1) ** 2 - eps))`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:94` quote: `expect_zero("S_n^(twin)(x=0) - S_n^(max)", sp.simplify(S_n.subs(x, 0) - S_n_max))`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:104-107` quote: `expect_zero("S_n^(max) - S_n^(twin) factored form", sp.simplify(ceiling_diff - ceiling_diff_target))`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:83` quote: `sNMax = FullSimplify[1 + (1 - eps)/((2 n + 1)^2 - eps), Assumptions -> $Assumptions]`

## Resolve before fix_loop

The stage 050 scripts verify an enhancement-ceiling theorem `S_n^(twin)(x) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` that the notes describe in Section 5 but the paper card's body and `\stagefield{Output}` line do not mention. Should this ceiling theorem be advertised in the paper card for stage 050, or should it be removed from this stage's scripts (and live elsewhere)?

Possible directions (the user picks one):
- (a) Paper card should advertise the ceiling — extend `stage_050.tex` with a fifth boxed equation `S_n^(twin) < S_n^(max) := 1 + (1-eps)/((2n+1)^2 - eps)` and update the Output line to mention "and higher-harmonic enhancement ceiling." No script change.
- (b) Ceiling belongs to a different stage — strip lines sympy:88-112 and wl:82-95 from the stage 050 scripts; relocate the ceiling check to the stage that the paper card identifies as its owner. (If no other stage currently advertises it, option (a) is the more conservative direction.)
- (c) Both stage 050's paper card and the notes Section 5 are derived from a third source that contradicts both the script's chosen scope — flag for deeper review (e.g., the ceiling may belong to a tower-construction stage rather than stage 050).

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## Applied: F2

- direction: (a) — paper card now advertises the enhancement ceiling.
- files_changed:
  - `paper/stages/stage_050.tex:43-51` — added boxed `S_n^{\rm twin}(x;\varepsilon) < S_n^{\rm max}(\varepsilon)` equation with label `eq:app-stage050-Sn-max`
  - `paper/stages/stage_050.tex:50` — Output line updated to reference the new ceiling
- summary: extended paper card with fifth boxed equation; no script change.
- deviation: none.

Codex must NOT re-apply F2 — leave the scripts' existing `S_n_max` checks intact (sympy:88-107, wl:82-95).

## F3 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:46-47`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl:43-44`

**Issue:**
The doubling theorem `S_0 = S(1; eps) = 2` (paper eq:app-stage050-S0) is asserted by direct substitution `zeta = 1` in both scripts, without ever evaluating `twin_support_ratio(0, x)` to confirm that the lowest twin lane indeed has `zeta_0 = 1`. The SymPy script imports `twin_support_ratio` from stage 049 (sympy:17) but does not use it at n=0; the Mathematica script does not import an analogue at all. The upstream chain "n=0 twin gives zeta=1, therefore S_0 = 2" is therefore not anchored to the stage 049 closed form the paper relies on.

**Required change:**

SymPy — insert immediately before the existing `expect_zero("S(1;eps) - 2", ...)` line at sympy:47, so the new anchor sits between sympy:46 (`print("S(zeta;eps) =", S)`) and the existing assertion:

```
expect_zero(
    "zeta_0^(twin) - 1 (anchors doubling to stage 049 import)",
    twin_support_ratio(sp.Integer(0), x) - 1,
)
```

Use `sp.Integer(0)` rather than a bare Python `0` because `twin_support_ratio`'s argument was declared with `positive=True` upstream; substituting an explicit symbolic zero avoids any positivity-assumption interaction during `simplify`.

Mathematica — insert immediately before the existing `expectZero["S(1;eps) - 2", ...]` line at wl:44, so the new anchor sits between wl:43 (`Print["S(zeta;eps) = ", fmt[sEnhance]];`) and the existing assertion:

```
expectZero[
  "zeta_0^(twin) - 1 (anchors doubling)",
  (1/((2 n + 1)^2 (1 + x n (n + 1))) /. n -> 0) - 1
];
```

The `n -> 0` substitution is a literal pattern replacement that fires before `FullSimplify` evaluates anything, so the global `$Assumptions` constraint `n >= 1` (wl:36) does not conflict — the substituted expression `1/((1)^2 (1 + 0)) - 1 = 0` reduces to literal 0 under any assumption.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 050` and `redteam exec-mathematica 050` and confirm:

1. The SymPy output contains a new line `zeta_0^(twin) - 1 (anchors doubling to stage 049 import) = 0`.
2. The Mathematica output contains a new `PASS: zeta_0^(twin) - 1 (anchors doubling)` line.
3. Both scripts still exit 0 with all other PASS lines intact.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage050_zeta_threshold_comparison_mathematica_audit.wl`
- summary: added the zeta_0 twin-lane anchor checks immediately before the existing S(1;eps) doubling assertions in both scripts.
- deviation: none.
