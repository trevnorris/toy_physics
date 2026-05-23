---
unit_id: 071
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 071

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

- `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py:65-71` — the tautological `expect_zero("eta - L/ell", eta - L / ell)` line was removed. Two new assertions were added: `Km_expected = sp.pi * a**2 * hbar**2 / (3 * m * rho_w)`, `expect_zero("K_m - pi a^2 hbar^2 / (3 m rho_w)", Km - Km_expected)`, and `expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected * L / Tx) - L / ell)`. The `Km`/`eta` definitions themselves are preserved (still printed for transparency).
- `mathematica/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.wl:77-83` — the tautological `expectZero["eta - L/ell", eta - L/ell]` line was removed. Two new assertions added: `kmExpected = Pi*a^2*hbar^2/(3*m*rhoW)`, `expectZero["K_m - pi a^2 hbar^2 / (3 m rhoW)", km - kmExpected]`, `expectZero["eta - L/ell (from closed-form K_m)", (kmExpected*L/tx) - L/ell]`.

The diff (`stage_071_diff.patch`) shows exactly these two edits and nothing else — no collateral changes to surrounding lines, no refactors elsewhere.

**Assessment:**

The edits match the directive's required-change blocks verbatim (SymPy variable naming `Km_expected` vs Mathematica `kmExpected`, parenthesization, and string labels all align). The new checks are substantive:

1. `Km - Km_expected = 0` requires `Tx/ell = pi*a^2*hbar^2/(3*m*rho_w)`. Since `Tx = pi*a^2*ell*If*hbar^2/(m*rho_w)`, this equality only holds when `If = 1/3`, i.e. when the shape integral has been computed correctly. A wrong factor (sign error, missing factor of `ell`, dropped factor of `hbar`, or mis-evaluated `If`) would cause this to fail — non-tautological.
2. `(Km_expected * L / Tx) - L/ell = 0` requires `Km_expected*L/Tx = L/ell`, i.e. `Km_expected*ell = Tx`. With `Km_expected = pi*a^2*hbar^2/(3*m*rho_w)` (a literal closed form independent of `Tx`'s definition), this exercises the full `Tx` factorization: any error in `Tx`'s `ell`-dependence or `If` value would propagate to a residual. Non-tautological.

Output transcripts confirm both new checks pass (`K_m - pi a^2 hbar^2 / (3 m rho_w) = 0` and `eta - L/ell (from closed-form K_m) = 0`) and that the previously offending line `eta - L/ell = 0` (without the qualifier) no longer appears as an assertion (the standalone `eta = L/ell` print line still appears as a transparency print, but it is not an `expect_zero` call).

## Exec log assessment

**SymPy:** exit=0 (inferred — no `stage_071_sympy.log` file present in `exec_logs/`, but the saved transcript `scripts/output/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.txt` (mtime newer than script) shows the full theorem-ledger banner printed at the end, which only happens if all `expect_zero` calls pass without raising `AssertionError`). Notable lines from the saved transcript:

- `K_m - pi a^2 hbar^2 / (3 m rho_w) = 0`
- `eta - L/ell (from closed-form K_m) = 0`
- `kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2] = 0`
- `W_wall_reduced - Upsilon_w Lambda_ell^2 = 0`

**Mathematica:** exit=0 (inferred — no `stage_071_mathematica.log` file present in `exec_logs/`, but the saved transcript `mathematica/output/moving_throat_pde_stage071_tanh_wall_branch_mathematica_audit.txt` shows all `PASS:` lines and the terminal theorem-ledger banner). Notable lines:

- `PASS: K_m - pi a^2 hbar^2 / (3 m rhoW)`
- `PASS: eta - L/ell (from closed-form K_m)`
- `PASS: T_X exact formula`, `PASS: K_X exact formula`, `PASS: J_1 exact formula`
- `PASS: kappa reduced law`, `PASS: W_wall reduced law`

**Output freshness:** confirmed. Script mtimes are 20:10, transcript mtimes are 20:11 — outputs are newer than scripts. The orchestrator did not write per-stage `stage_071_sympy.log` / `stage_071_mathematica.log` files (only `stage_071_diff.patch` is present in `exec_logs/`), but the `.txt` transcripts under `scripts/output/` and `mathematica/output/` serve the same evidentiary role and are post-fix.

## Material-change assessment

`material_change`: false.

The fix replaced a tautological identity (`eta = L/ell` under the definition `K_m := T_X/ell`) with a non-tautological pin on the closed form of `K_m`. No symbolic quantity (`T_X`, `K_X`, `J_1`, `W_wall`, `kappa`, `eta`) changed value or definition — only the assertion was strengthened. Downstream units that consume `K_m`, `eta`, or any wall-energetics result will see exactly the same closed forms as before; no result that any later stage could depend on has been altered.

## Side observations (non-blocking)

- The theorem-ledger print at the end of both scripts still reads `eta = L/ell under K_m = T_X/ell`. This phrasing is accurate to the script's semantics (the definition `K_m := T_X/ell` is what makes `eta = L/ell`), but a reader scanning only the ledger might still assume the equality is independently checked. The new `K_m - pi a^2 hbar^2 / (3 m rho_w)` assertion does provide the independent check, so the ledger is no longer misleading in substance. Not a blocking concern.
- The exec_logs directory holds only the diff patch for unit 071, no per-engine `.log` files. The saved `.txt` transcripts are sufficient evidence here, but if the orchestrator expects to find `stage_071_{sympy,mathematica}.log` for downstream tracking, that may be worth noting separately.
- The expected `eta` print line (transparency-only) still shows `eta = L/ell` in both transcripts. That is the value of the simplified symbolic `eta` and is a print statement, not an `expect_zero` assertion, so it does not reintroduce the tautology. Confirmed by reading the script source (line 66/68 in SymPy, line 78/80 in Mathematica) — these are pure `print`/`Print` calls.

## Verdict justification

The single finding F1 is fully addressed: the tautological `expect_zero("eta - L/ell", ...)` assertion is removed from both engines and replaced by two substantive checks that pin `K_m` to its literal closed form and reconstruct `eta` from the pinned `K_m`. Both checks pass in the saved transcripts (post-fix mtimes). The diff is minimal (exactly the lines named in the directive, no collateral edits). No regressions are visible in the diff or the transcripts. All other assertions from the original report (A1–A5, A7, A8, B1–B8, B10, B11) continue to pass as before. Verdict: `verified`.
