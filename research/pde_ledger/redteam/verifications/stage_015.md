---
unit_id: 015
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T21:14:30Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 015

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
Codex created the new file `mathematica/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.wl` (200 lines). It defines `expectZero`, `expectNonzero`, and `expectEqual` helpers that each call `Exit[1]` on failure, then walks through the M1-M9 claim manifest with distinct local names (`lagrangian`, `effectiveMass`, `gateMatrix`, `gauntBase`, `realY20Ratio[m_]`) rather than the SymPy names (`L`, `K_eta`, `wall_matrix`, `real_y20_square_ratio`). The saved output `mathematica/output/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.txt` (mtime 2026-05-21 15:01) records `STATUS: PASS` with every M1-M9 PASS line printed, post-dating the script's mtime (13:12), so the file is freshly regenerated.

**Assessment:**
The script is non-tautological and independent of the `.py` choreography in the ways the directive required:

- M3 reads `L2raw = Coefficient[Series[lagrangian, {eps, 0, 2}] // Normal, eps, 2]` rather than the `.py`'s `sp.diff(sp.expand(L), eps, 2).subs(eps, 0) / 2`. Two distinct primitives compute the eps^2 coefficient; the residual check `L2afterIBP - canonicalL2` equates them and the sign-mutation guard (`URR0 + dTwRR0p ...` instead of `URR0 - dTwRR0p ...`) returns `dTwRR0p eta^2`, correctly nonzero.
- M4 builds `D01full = dK - b01 - z01`, `D21full = -(dM + b21 + z21)`, `D41full = -(b41 + z41)` explicitly and then specializes with `wallSpec = {b01 -> 0, b21 -> 0, ...}` — the directive's "do NOT skip the construction" requirement is honored.
- M5 evaluates the Gaussian overlap integrals to closed forms and additionally asserts the *values* `dMoverlap == Sqrt[Pi/3]` and `dKoverlap == 23 Sqrt[Pi] / (3 Sqrt[3])`, which is a stronger substantive guard than the `.py`. I cross-checked: `Integrate[exp(-3w^2), -inf, inf] = Sqrt[Pi/3]`, and the dK breakdown `4*w^2 exp(-3w^2) + 7 exp(-3w^2)` integrates to `(2/3 + 7) Sqrt[Pi/3] = 23 Sqrt[Pi/3] / 3`, matching the closed form. The 6→5 mutation residual `-1/9 * Sqrt[Pi/3]` is exactly what one expects from the integral difference -Sqrt[Pi/3] divided by 9.
- M6 reproduces the `Det = 1/27` and the perturbation shift `2 eps/3` exactly (I checked: `(2/3)(1/9+eps) - 1/27 = 1/27 + 2 eps/3`). Both an equality-to-`2 eps/3` check and a nonzero-guard are present.
- M7 yields `dK = -18 eps`, `dM = -eps` for the perturbed solve, matching the closed-form derivation by hand (Hevenwall = 0 forces `dK = 18 dM`, substituted gives the printed values).
- M8 uses `ThreeJSymbol[{2,0},{2,m},{2,-m}]/ThreeJSymbol[{2,0},{2,0},{2,0}]` with no m=0 short-circuit; the ratio is well-defined because all three j's are 2 so the Gaunt prefactors cancel, and Wigner 3-j values are nonzero. Same-sign cross terms are checked via the `m + m + 0 = 0` selection-rule violation, yielding 0 from `ThreeJSymbol` directly.
- M9 reuses the canonical real-Y20 weights `(1, 1/2, -1)` directly to assemble the grouped trace and `b = 3a` identity.

Engine cross-check is now non-trivially satisfied: independent code paths in two systems agree on every numerical residual.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:129-167` (new block), Codex deleted the two old `wall-only K1 from overlap-generated slots` / `wall-only H_even from overlap-generated slots` assertions and replaced them with concrete Gaussian profiles for `beta`, `delta_mu`, `delta_Tw`, `delta_TO`, `delta_Keta`, computed `dM_overlap_concrete` and `dK_overlap_concrete` via `sp.integrate(..., (wall_w, -sp.oo, sp.oo))`, and added three new assertions:
- `wall-only K1 from concrete Gaussian overlap integrals` (concrete-numerical check)
- `wall-only H_even from concrete Gaussian overlap integrals` (concrete-numerical check)
- `wall-only K1 detects mutated 6*delta_TO coefficient` (substantive coefficient guard via `assert_nonzero` against a 5-instead-of-6 mutated overlap)

**Assessment:**
The new K1 check is no longer the substitution-rename tautology the auditor flagged. With `dM_overlap_concrete = Sqrt[pi/3]` and `dK_overlap_concrete = 23 Sqrt[pi]/(3 Sqrt[3])` (both closed-form sympy integrals), the residual `K1_wall.subs(...) - (-dM_concrete + dK_concrete/9)` exercises concrete `sqrt(pi)` arithmetic. Crucially, the mutation guard exercises the `6 delta_TO` coefficient: changing 6 to 5 shifts `dK_overlap_concrete` by `-Sqrt[pi/3]` and the K1 residual by `-Sqrt[pi/3]/9`, which `assert_nonzero` correctly detects. The Mathematica mirror reports this exact residual `-1/9 Sqrt[Pi/3]` in M5, agreeing with the sympy form.

Side note: the symbolic `dM_overlap`/`dK_overlap` `sp.Integral` definitions at lines 109-113 still exist (unused now) but are inert — they don't paper over the new check.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage015_parent_throat_action_master_sympy_audit.py:25-31`, Codex removed the `if m == 0: return sp.Integer(1)` shortcut. The function now uniformly returns `sp.simplify((sp.Integer(-1) ** m) * gaunt(2, 2, 2, 0, m, -m) / base)` for all m, with the same-sign sanity check gated on `m != 0` (correctly, since m=0 is structurally self-negating).

**Assessment:**
For m=0, the returned value is now `gaunt(2,2,2,0,0,0) / gaunt(2,2,2,0,0,0)`, which sympy reduces to 1 only after evaluating the Wigner 3-j. The assertion `assert_zero("Y20 overlap lane 20", lam20 - 1)` is therefore no longer a `1 - 1 == 0` tautology; the Gaunt machinery actually runs and divides. A mutation to the `base` definition (e.g., multiplying by 2) would now flip `lam20 = 1/2`, surfacing the regression, exactly as the directive specified. The Mathematica mirror in M8 uses the same `gauntBase` denominator structure and produces residual 0 for m=0, confirming the lanes agree across engines.

## Exec log assessment

**SymPy:** exit=0. Notable lines from `redteam/exec_logs/stage_015_sympy.log`:
- `STEP 13 PARENT THROAT ACTION MASTER AUDIT`
- `Boundary operator nonzero sanity check = PASS`
- `Wall-only perturbed-gate solve guard = PASS`
- `Wall-only coefficient determinant guard = PASS`
- `STATUS: PASS`
- `# exit_code: 0`

The sympy `assert_zero`/`assert_nonzero` helpers only raise on failure — they do not print per-assertion PASS lines — so the new `wall-only K1 from concrete Gaussian overlap integrals`, `wall-only H_even from concrete Gaussian overlap integrals`, and `wall-only K1 detects mutated 6*delta_TO coefficient` checks are silently passing (no exception). The clean exit confirms all 23 assertions, including the three new ones, hold.

**Mathematica:** exit=0 (inferred from saved output's final `STATUS: PASS` and the script's structure where any individual claim failure calls `Exit[1]` immediately). No `redteam/exec_logs/stage_015_mathematica.log` file exists in the exec_logs directory — the orchestrator's mathematica log is missing — but the saved output transcript at `mathematica/output/moving_throat_pde_stage015_parent_throat_action_master_mathematica_audit.txt` shows every PASS line through M9 and terminates with `STATUS: PASS`. Notable lines:
- `M3 K_eta canonical quadratic form residual = 0` / `PASS`
- `M3 K_eta dTwRR0p sign mutation residual = dTwRR0p*eta^2` / `PASS` (correctly nonzero)
- `M5 Gaussian dM overlap closed form residual = 0` (i.e., dMoverlap = Sqrt[Pi/3])
- `M5 Gaussian dK overlap closed form residual = 0` (i.e., dKoverlap = 23 Sqrt[Pi]/(3 Sqrt[3]))
- `M5 wall-only K1 detects mutated 6 deltaTO coefficient residual = -1/9*Sqrt[Pi/3]` / `PASS`
- `M6 wall determinant perturbation value residual = 0` (i.e., shift = 2 eps/3 exactly)
- `M7 perturbed solve nonzero dK residual = -18*eps`, `dM residual = -eps`
- `M8 real-Y20 ratio m=0 residual = 0`, `m=1 residual = 0`, `m=2 residual = 0`
- `STATUS: PASS`

**Output freshness:**
- `scripts/.../stage015_sympy_audit.py` mtime 2026-05-21 13:11:45; output mtime 15:00:34 — output is newer.
- `mathematica/.../stage015_mathematica_audit.wl` mtime 2026-05-21 13:12:26; output mtime 15:01:19 — output is newer.
Both saved outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

No derived numerical result quoted forward by downstream units has changed. The substantive content (`K_eta = URR0 - d_TwR_R0p + TwRR0 R0p^2/2`, det = 1/27, Gaunt ratios (1, 1/2, -1)) is unchanged — only the verification surface around it has tightened (added Mathematica mirror, replaced one tautological pair with concrete-integral checks, removed one hardcoded m=0 lane). Downstream units depending on stage 015 results need not re-verify their own assertions.

## Side observations (non-blocking)

- The `redteam/exec_logs/stage_015_mathematica.log` file is absent. The orchestrator may want to capture a fresh `redteam exec-mathematica 015` log alongside the existing sympy log for symmetry; the saved `.txt` output provides equivalent ground truth for now.
- The symbolic `dM_overlap` / `dK_overlap` `sp.Integral` objects at lines 109-113 of the sympy script are now dead code (no assertion references them after F2's removal of lines 129-136). Not blocking — they document the symbolic forms — but they could be deleted in a future janitorial pass.
- The directive's M5 step said `dKoverlap` and `dMoverlap` "should evaluate to closed-form `Sqrt[Pi]`-multiples" without specifying values; Codex went further and added explicit `expectZero` checks against `Sqrt[Pi/3]` and `23 Sqrt[Pi] / (3 Sqrt[3])`. This is a positive deviation beyond the directive — stronger substantive content.

## Verdict justification

All three findings are resolved with non-tautological, substance-bearing edits. The new Mathematica companion independently re-derives every claim using the required primitives (`Series`/`Coefficient`, explicit `D01full/D21full/D41full` construction, `ThreeJSymbol` without short-circuit, Gaussian `Integrate` to closed forms), exits 0, and prints PASS for all M1-M9. The sympy script's old substitution-rename pair is gone, replaced by concrete Gaussian integrals plus a coefficient mutation guard whose `-Sqrt[pi/3]/9` residual matches the Mathematica mirror's `-1/9 Sqrt[Pi/3]` value byte-for-byte. The m=0 Gaunt short-circuit is gone. Exec logs (sympy explicit, mathematica via saved output) confirm exit 0. No regressions in the diff. Verified.
