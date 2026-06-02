---
unit_id: 207
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-02T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 207

## Per-finding outcomes

### F1 — missing_verification_script (subtype `missing_mathematica`)

**Classification:** resolved

**What changed:**
Codex authored the missing second engine at the registered path
`mathematica/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_mathematica_audit.wl`
(7371 bytes, untracked `??` so it does not appear in the diff patch — expected for a brand-new file plus its output `.txt`). No other file was touched: the SymPy script's last commit is 2026-05-11 and `git status` shows no working-tree change to the 207 `.py`. The `.wl` carries an `expectZero`/`fail`+`Exit[1]` harness with the required idioms: `stripConditional` removes `ConditionalExpression[e_, _] :> e` before every zero test (script lines 27, 134, 181), and the run exits 0.

**Assessment:**
Each manifest item M1–M6 appears as an explicit, non-tautological check and PASSes:

- **M1** (lines 95–102): symbolic symmetric 5×5 `H`, contracts `e.H.e - H[[i,i]]` and `(-e).(-e)` form for all five axes → 10 PASS.
- **M2** (lines 104–120): expands `ray.H.ray` for all **10** axis pairs and additionally extracts `Coefficient[form, a b] - 2 H[[i,j]]` (factor-of-2 load-bearing) → 20 PASS. (SymPy checks only the single `(lambda,c)` pair.)
- **M3** (lines 122–125): `-Sign[Gam]*Gam + Abs[Gam]` with conditional stripping → 1 PASS.
- **M4** (lines 127–174): slope-convention coefficient `Coefficient[poly, tau] + k == 0`, plus `Solve`+`SelectFirst` smaller-positive-root selection and an independent `Reduce[... && tau>0, Reals]` set print for both envelopes, then equality with closed form `2H0/(k+Sqrt[k^2-2cH0])` and quadratic residual → 7 PASS.
- **M5** (lines 176–198): `Solve` of the turning quadratic, positive-root `SelectFirst` under `curvature<0`, equality with `Sqrt[-2H0/kappa]`, residual, both labels → 4 PASS.
- **M6** (lines 200–231): native `primitiveDrift[dir]` evaluated on `eps_i*e_i` against the tabulated rows for all five axes → 5 PASS.

None is tautological: M1/M2 contract a fully symbolic matrix against an independently written target; M4/M5 obtain the root by Solve/Reduce and only afterward compare to the closed form (no construct-then-substitute echo). Total 46 checks PASS, 0 FAIL.

## Exec log assessment

**SymPy:** exit=n/a. No `stage_207_sympy.log` present (SymPy side unchanged this iteration; per directive the `.py` may not be edited). The committed SymPy output remains exit-0 PASS per the original audit.

**Mathematica:** exit=0 (`# exit_code: 0`, log line 135). Notable lines:
```
PASS: M2 cross coefficient residual (lambda,c)
M4 tau_lo selected root - closed form = 0
PASS: M4 tau_lo selected root - closed form
PASS: M4 tau_lo quadratic residual
PASS: M5 tau_lo^(tp) - Sqrt[-2H0/kappa_lo]
STAGE 207 MATHEMATICA AUDIT PASSED
```
45 `PASS:` lines + the M4 slope-convention check that PASSes = 46 manifest checks; 0 `FAIL:`. The M4 `Reduce` positive-root sets are printed (lines 90–91) showing the branch logic, not asserted as residuals (informational), consistent with the harness.

**Independent-derivation check (load-bearing):** Confirmed independent, not a transliteration.
- M4: SymPy *constructs* `TauL = 2*H0/(k+sqrt(k**2-2*cL*H0))` (`.py` lines 177–178) and substitutes back; the `.wl` instead `Solve`s the comparison quadratic and `SelectFirst`s the smaller positive root (`.wl` lines 131–160), then checks equality with the closed form. This is the exact decomposition reversal the directive required.
- M5: SymPy writes `sqrt(2*H0/a_turn)` directly (`.py` lines 201–202); the `.wl` `Solve`s the turning quadratic and selects the positive root under `curvature<0` (`.wl` lines 179–190).
- M2: `.wl` adds `Coefficient[form, a b]` extraction and sweeps all 10 pairs vs SymPy's single pair.
M1/M3/M6 are structurally parallel (the table values and Hessian identities are the physics, so the targets must match) but use native idioms (`UnitVector`+dot, `Sign`/`Abs`, `primitiveDrift` index access) rather than re-typed SymPy `Matrix`/`subs` choreography.

**Output freshness:** Confirmed. `.wl` mtime 2026-06-02 10:06:44; output `.txt` mtime 2026-06-02 10:07:49 (newer than the script). Filtered body of the committed `.txt` is byte-identical to the exec-log body (`OUTPUT TXT == LOG BODY`).

## Material-change assessment

`material_change`: false. F1 adds a second verification engine only; no SymPy edit, no change to any derived constant, bracket, or forward result. Downstream units depend on Stage 207's derived quantities, which are unchanged. No `upstream_stale` concern arises from this edit beyond the orchestrator's standard bookkeeping.

## Side observations (non-blocking)

- The original auditor flagged stale display banners in the SymPy script ("STAGE 190" at `.py` lines 35/244 and "from Stage 204" comment at line 104). The new `.wl` correctly uses "STAGE 207" throughout, so no banner staleness in the second engine. The SymPy cosmetic labels remain but were explicitly out of scope (directive forbids editing the `.py`) and carry no math impact. Not blocking.

## Verdict justification

The single finding F1 is fully resolved: the registered `.wl` exists, exits 0, and verifies every claim-manifest item M1–M6 with 46 non-tautological PASS checks (0 FAIL), confirmed against the orchestrator's independent re-run. The derivation is genuinely independent of the SymPy script — most importantly the load-bearing M4/M5 brackets are obtained via Solve/Reduce/SelectFirst rather than the SymPy construct-then-substitute, exactly as the anti-transliteration guard demanded. SymPy was untouched and no derived result changed, so `material_change` is false. Verdict: `verified`.
