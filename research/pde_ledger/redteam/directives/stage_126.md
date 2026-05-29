---
unit_id: 126
batch: IV.3
created_at: 2026-05-29T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-29T14:32:24-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 126

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts and iterate until each exits 0 with all in-file checks passing:
- SymPy: `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py`
- Mathematica: `timeout 600 math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl` (or the exec runner `redteam exec-mathematica 126` if that is the configured path).

Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:55-66` (the `min_xi0` / `val_xi1` / `sigma_corner` endpoint-and-corner block)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl:64-75` (the `sigmaMatchAtLM` / `sigmaMatchAt0` / `sigmaXiAtLM` / `valXiAt1` boundary-value block)

**Issue:**
The paper/notes claim is the GLOBAL positivity of the convex mouth-source family:
`sigma_xi(z) = (1-xi)*k*cos(k*z) + xi/L >= 0` for all `0 <= z <= L`, `0 <= xi <= 1`, `L > 0`
(notes §3, `sigma_xi` boxed at lines 71-77 and the explicit "normalized and nonnegative on [0,L]" claim at line 78).

The current checks only sample the boundary of the `(z, xi)` box: `sigma_xi(xi=0)` min over z (SymPy:55-60 / Mathematica:64-72), the value at `xi=1` (SymPy:56,61-62 / Mathematica:73-75), and the corner `(z=L, xi=0)` (SymPy:63-66 / Mathematica:70-72). An interior negative perturbation can pass all of these. For example a sign-changing family
`sigma_xi + (-A)*xi*(1-xi)*sin(2*pi*z/L)/L` (A>0)
vanishes at `xi=0` and `xi=1`, preserves the `xi=0` and `xi=1` slice checks and the corner, yet dips negative for some interior `(z, xi)`. The current block cannot detect that, so it does not exercise the load-bearing global-positivity claim.

**Structural fact that makes a clean global check possible:**
`sigma_xi` is AFFINE (degree 1) in `xi`: it equals `A(z) + xi*(B - A(z))` with `A(z) = k*cos(k*z)`, `B = 1/L`, so `d^2 sigma_xi / d xi^2 = 0`. A function affine in `xi` attains its minimum over `xi in [0,1]` at an endpoint (`xi=0` or `xi=1`). Therefore:
(i) `sigma_xi` is affine in `xi`, AND
(ii) both endpoint sources `sigma_xi|_{xi=0} = k*cos(k*z)` and `sigma_xi|_{xi=1} = 1/L` are `>= 0` for all `z in [0,L]`
together IMPLY `sigma_xi >= 0` on the whole box. We verify BOTH (i) and (ii) so an interior-negative family (which is NOT affine in `xi`) fails the affineness check, and a sign-changing endpoint slice fails the min-over-z check.

**Required change (SymPy):**

Replace lines 55-66 (the `min_xi0` / `val_xi1` / `sigma_corner` block). Keep the existing comment at lines 48-50 and the `min_match` check at lines 51-54 unchanged. The new block proves global positivity in two non-tautological pieces: an exact affine-in-xi structural test, and an endpoint-minima test over the full interval `[0,L]` (not just sampled corners).

Before (lines 55-66):
```python
min_xi0 = sp.calculus.util.minimum(sigma_xi.subs(xi, 0), z, sp.Interval(0, L))
val_xi1 = sp.simplify(sigma_xi.subs(xi, 1))
print("min sigma_xi(xi=0) on [0,L] =", min_xi0)
print("sigma_xi(xi=1) =", val_xi1)
if sp.simplify(min_xi0) != 0:
    raise AssertionError("sigma_xi(xi=0) minimum on [0,L] should be 0.")
if sp.simplify(val_xi1 - 1/L) != 0:
    raise AssertionError("sigma_xi(xi=1) should equal 1/L.")
sigma_corner = sp.simplify(sigma_xi.subs([(z, L), (xi, 0)]))
print("sigma_xi(z=L, xi=0) =", sigma_corner)
if sp.simplify(sigma_corner) != 0:
    raise AssertionError("sigma_xi(z=L, xi=0) should equal 0.")
```

After (lines 55-66):
```python
# Global positivity on the box z in [0,L], xi in [0,1].
# sigma_xi is AFFINE in xi: sigma_xi = k*cos(k*z) + xi*(1/L - k*cos(k*z)),
# so d^2 sigma_xi / d xi^2 == 0 and its minimum over xi in [0,1] is attained
# at an endpoint. Affine-in-xi (i) + both endpoint slices nonneg on [0,L] (ii)
# => sigma_xi >= 0 on the whole box. An interior-negative perturbation is NOT
# affine in xi and fails (i); a sign-changing endpoint slice fails (ii).
d2_xi = sp.simplify(sp.diff(sigma_xi, xi, 2))
print("d^2 sigma_xi / d xi^2 =", d2_xi)
if d2_xi != 0:
    raise AssertionError("sigma_xi must be affine in xi for the endpoint-min argument.")
min_xi0 = sp.calculus.util.minimum(sigma_xi.subs(xi, 0), z, sp.Interval(0, L))
min_xi1 = sp.calculus.util.minimum(sigma_xi.subs(xi, 1), z, sp.Interval(0, L))
print("min sigma_xi(xi=0) on [0,L] =", min_xi0)
print("min sigma_xi(xi=1) on [0,L] =", min_xi1)
if not (sp.simplify(min_xi0) >= 0) is sp.true:
    raise AssertionError("sigma_xi(xi=0) must be nonnegative on [0,L].")
if not (sp.simplify(min_xi1) >= 0) is sp.true:
    raise AssertionError("sigma_xi(xi=1) must be nonnegative on [0,L].")
# Endpoint identities (informational, pin the known minima/values).
if sp.simplify(min_xi0) != 0:
    raise AssertionError("sigma_xi(xi=0) minimum on [0,L] should be 0 (at z=L).")
if sp.simplify(sigma_xi.subs(xi, 1) - 1/L) != 0:
    raise AssertionError("sigma_xi(xi=1) should equal 1/L.")
```

Notes for Codex:
- `sp.calculus.util.minimum(1/L, z, sp.Interval(0, L))` on the constant `xi=1` slice returns `1/L`; with `L` declared `positive` (line 10), `sp.simplify(1/L) >= 0` evaluates to `S.true`. If for any reason `(... >= 0) is sp.true` does not hold for a slice min that is manifestly nonneg, do NOT weaken the test to a tautology — instead use `bool(sp.simplify(min_xiK) >= 0)` guarded by the `L > 0` assumption, and report it in `## Applied` as a deviation. The nonnegativity assertions (the two `>= 0` checks) are the load-bearing part; the two trailing identity checks (`min_xi0 == 0`, slice `== 1/L`) pin exact known values and must stay.
- The affineness check `sp.diff(sigma_xi, xi, 2) == 0` is exact symbolic and CANNOT pass for an interior-curved perturbation in `xi`; that is what closes the interior gap.

**Required change (Mathematica):**

Replace lines 64-75 (the `sigmaMatchAtLM` / `sigmaMatchAt0` / `sigmaXiAtLM` / `valXiAt1` boundary-value block). Keep the comment at lines 61-63 and everything above line 64 unchanged. The new block mirrors the SymPy logic: an exact affine-in-xi test via `D[..., {xi, 2}]`, plus a global-positivity confirmation over the full box using `Reduce`/`Resolve[ForAll[...]]` with explicit variable domains.

Before (lines 64-75):
```
sigmaMatchAtLM = FullSimplify[k*Cos[k*z] /. z -> lM, Assumptions -> $Assumptions];
Print["sigma_match(lM) = ", fmt[sigmaMatchAtLM]];
expectZero["sigma_match(lM) = 0 (boundary min on [0,lM])", sigmaMatchAtLM];
sigmaMatchAt0 = FullSimplify[k*Cos[k*z] /. z -> 0, Assumptions -> $Assumptions];
Print["sigma_match(0) = ", fmt[sigmaMatchAt0]];
expectZero["sigma_match(0) = k (interior max on [0,lM])", sigmaMatchAt0 - k];
sigmaXiAtLM = FullSimplify[(sigmaXi /. xi -> 0) /. z -> lM, Assumptions -> $Assumptions];
Print["sigma_xi(z=lM, xi=0) = ", fmt[sigmaXiAtLM]];
expectZero["sigma_xi(z=lM, xi=0) = 0", sigmaXiAtLM];
valXiAt1 = FullSimplify[sigmaXi /. xi -> 1, Assumptions -> $Assumptions];
Print["sigma_xi(xi=1) = ", fmt[valXiAt1]];
expectZero["sigma_xi(xi=1) = 1/lM", valXiAt1 - 1/lM];
```

After (lines 64-75):
```
(* Global positivity on z in [0,lM], xi in [0,1]. sigmaXi is affine in xi, so
   D[sigmaXi, {xi, 2}] == 0; its min over xi in [0,1] is at an endpoint. We also
   confirm the whole-box claim directly with Resolve[ForAll[...]] over explicit
   domains, so an interior-negative source FAILS this block. *)
d2Xi = FullSimplify[D[sigmaXi, {xi, 2}], Assumptions -> $Assumptions];
Print["d^2 sigma_xi / d xi^2 = ", fmt[d2Xi]];
expectZero["sigma_xi affine in xi", d2Xi];
sigmaMatchAtLM = FullSimplify[(k*Cos[k*z] /. z -> lM), Assumptions -> $Assumptions];
Print["sigma_match(lM) = ", fmt[sigmaMatchAtLM]];
expectZero["sigma_match(lM) = 0 (boundary min on [0,lM])", sigmaMatchAtLM];
valXiAt1 = FullSimplify[(sigmaXi /. xi -> 1), Assumptions -> $Assumptions];
Print["sigma_xi(xi=1) = ", fmt[valXiAt1]];
expectZero["sigma_xi(xi=1) = 1/lM", valXiAt1 - 1/lM];
globalPos = Resolve[
  ForAll[{z, xi}, 0 <= z <= lM && 0 <= xi <= 1, sigmaXi >= 0],
  Reals
];
globalPos = Simplify[globalPos, Assumptions -> (lM > 0)];
Print["ForAll sigma_xi >= 0 on box -> ", fmt[globalPos]];
If[!TrueQ[globalPos], fail["sigma_xi >= 0 on box (z in [0,lM], xi in [0,1])", globalPos]];
```

Notes for Codex (Mathematica hazards):
- `Resolve[ForAll[...], Reals]` over the box with `lM` still symbolic positive may return a residual `ConditionalExpression` or an unsimplified expression rather than literal `True`. The `Simplify[..., Assumptions -> (lM > 0)]` step is there to collapse it to `True`. If it does NOT collapse to `True`, do NOT delete the check or rewrite it as a tautology. First try `Resolve[ForAll[{z, xi}, 0 <= z <= lM && 0 <= xi <= 1, sigmaXi >= 0] && lM > 0, Reals]`, or substitute a concrete positive value `lM -> 1` ONLY for this positivity quantifier (the family is scale-covariant in `lM` for positivity, since `sigmaXi = k Cos[k z] (1-xi) + xi/lM` with `k = Pi/(2 lM)` and both terms scale as `1/lM`), e.g. `Resolve[ForAll[{z, xi}, 0 <= z <= 1 && 0 <= xi <= 1, (sigmaXi /. lM -> 1) >= 0], Reals]`. Record whichever form you use in `## Applied` as a deviation.
- Keep the `expectZero` LABELS for `sigma_match(lM)` and `sigma_xi(xi=1)` unchanged so output strings stay comparable. The `sigma_match(0) = k` and `sigma_xi(z=lM, xi=0) = 0` corner lines are subsumed by the affineness + ForAll checks and may be dropped (they were redundant corner samples); if you prefer to keep them for traceability, leave them but they are not load-bearing.
- Do NOT introduce the literal two-character sequence `*` followed by `)` inside any comment.
- `lM > 0` is already in `$Assumptions` (line 29); the explicit `0 <= z <= lM` domain inside `ForAll` is what bounds `z` for the quantifier.

**Why this can now fail (non-tautology):**
- The affineness check is exact: any family with curvature in `xi` (the canonical way to hide an interior dip while preserving both endpoint slices) yields `d^2/dxi^2 != 0` and FAILS.
- The endpoint-min checks (SymPy) and the `ForAll`/`Resolve` box check (Mathematica) test the full `z`-interval, not sampled corners, so a sign-changing slice FAILS. The two together close the interior gap the corner-only checks missed.

**Claim manifest:**
- **M1** — `d^2 sigma_xi / d xi^2 == 0` (SymPy) / `expectZero["sigma_xi affine in xi", d2Xi]` (Mathematica) verifies the structural premise that `sigma_xi` is affine in `xi`, so its min over `xi in [0,1]` is at an endpoint. Claim anchor: convex family definition `sigma_xi = (1-xi) k cos(k z) + xi/L`, `notes/stages/moving_throat_pde_stage126_positive_source_families.md:71-77` (boxed) and `paper/stages/stage_126.tex:16` ("a positive convex family hits the lower branch").
- **M2** — endpoint-min nonnegativity (SymPy `min_xi0 >= 0`, `min_xi1 >= 0`) / global `Resolve[ForAll[..., sigmaXi >= 0]]` (Mathematica) verifies the load-bearing claim `sigma_xi(z) >= 0` on `0 <= z <= L`, `0 <= xi <= 1`, `L > 0`. Claim anchor: "This is normalized and nonnegative on `[0,L]`", `notes/stages/moving_throat_pde_stage126_positive_source_families.md:78`; positivity directive `paper/stages/stage_126.tex:22` ("Check positivity of the mouth source; do not use signed fitting to reach the upper branch").

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 126` and `redteam exec-mathematica 126`. Expected:
- SymPy: new lines `d^2 sigma_xi / d xi^2 = 0`, `min sigma_xi(xi=0) on [0,L] = 0`, `min sigma_xi(xi=1) on [0,L] = 1/L`; the `sigma_xi(xi=1) should equal 1/L` and `min_xi0 == 0` identities still hold; all prior lines (`min sigma_match on [0,L] = 0`, `xi_*`, interval `2/pi < g_- < pi/4`) unchanged; exits 0.
- Mathematica: `PASS: sigma_xi affine in xi`; `ForAll sigma_xi >= 0 on box -> True`; `PASS: sigma_xi(xi=1) = 1/lM`; final `Stage 126 Mathematica audit passed.`; `Exit[0]`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl`
- summary: Replaced the endpoint-and-corner sampling blocks with affine-in-xi checks plus full interval/box positivity verification.
- deviation: Mathematica uses the directive-approved `lM -> 1` scale-covariant positivity quantifier because symbolic `Resolve` with free positive `lM` did not collapse to `True`.

## Applied (reconcile note for already-satisfied PASS rows)

The following rows from the Codex review table are already correctly applied (tainted-applied in a prior pass) and are confirmed present/correct on read. They are recorded here as already-satisfied and require NO new edits:

- **Row 1 (PASS)** — `scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:51-54`: `min_match = sp.calculus.util.minimum(k*cos(k z), z, [0,L])` with assertion `min_match == 0`. Present and correct; can fail for a wrong endpoint/sign change. LEAVE UNCHANGED.
- **Row 4 (PASS)** — `scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py:80-83`: `interval_check = bool(2/pi < g_- < pi/4)` with `raise` on failure. Present and correct; brackets the lower branch. LEAVE UNCHANGED. (Mathematica analogue at `mathematica/...wl:84-86` likewise correct.)
- **Row 5 (PASS)** — banner `STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES` at `scripts/...sympy_audit.py:13` and `mathematica/...mathematica_audit.wl:26`. Correct Stage 126 label. LEAVE UNCHANGED.

Codex: do NOT re-prescribe or re-touch these. Only F1 (the SymPy:55-66 and Mathematica:64-75 blocks) needs editing.
