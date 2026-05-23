---
unit_id: 070
batch: III.3
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 070

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex rewrote `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:26-79` exactly as the directive specified:

- Introduced primitive-only closed forms at the top: `kappaClosed = 4*(m*cSw*L/hbar)^2 + (Ig/IfMoment)*(L/ell)^2` (lines 33-36) and `WwallClosed = 4*rhoW^2*V0^2*L^2/(hbar^2*cSw^2*ell^2)` (lines 37-40).
- Kept the physical-primitive intermediates `Hw, Nphiphi, Gphiphi, Tx, Kx` (lines 42-48) and renamed `kappa -> kappaAssembled` (line 49).
- Replaced `Wwall`'s former reliance on `J1 = IfMoment/Hw` with an inlined form `WwallAssembled = 4*Pi*a^2*L^2*(IfMoment*rhoW/(m*cSw^2))*V0^2/(Tx*ell)` (lines 59-62), i.e. `1/Hw -> rhoW/(m*cSw^2)` substituted in place.
- Replaced `Xi = gphi^2 * I1 * L^2 / Tx` with `Xi = (V0/ell)^2 * (4*Pi*a^2*ell*IfMoment*rhoW/(m*cSw^2)) * L^2 / Tx` (lines 67-70), again inlining `1/Hw` and `V0/ell` rather than naming them.
- Dropped the SymPy-mirror intermediates `J1`, `gphi`, `I1` entirely (confirmed via diff and via search of the post-fix `.wl`).
- Retained the three `expectZero` calls at lines 57, 65, 73 (one each for `kappa`, `W_wall`, `Xi - W_wall`), matching the original assertion count and the directive's claim manifest M1/M2/M3.
- Helpers (lines 1-24), banner, theorem-ledger text, and `Exit[0]` are unchanged.

**Assessment:**
The edit correctly addresses the transliteration finding. The Mathematica script no longer mirrors the SymPy script's variable choreography:

- SymPy uses three named intermediates (`J_1 = If/Hw`, `gphi = V0/ell`, `I_1 = Nphiphi/Hw`) to assemble `Wwall` and `Xi`. The new Mathematica version uses none of these; the `1/Hw` factor is inlined as `rhoW/(m*cSw^2)` at the points of use.
- The verification direction is now "assembled built from physical primitives vs. closed form written directly in primitives," consistent with the directive's "reverse the verification direction" instruction.
- The three new assertions are non-tautological:
  - `kappaAssembled - kappaClosed`: `kappaAssembled` is `Kx*L^2/Tx` with `Kx`, `Tx`, `Hw`, `Nphiphi`, `Gphiphi` built from primitives via four nested algebraic combinations; `kappaClosed` is the direct two-term closed form. Their equality is the non-trivial claim of the unit.
  - `WwallAssembled - WwallClosed`: `WwallAssembled` carries the factor `(IfMoment*rhoW/(m*cSw^2))` divided by `Tx` (which contains `IfMoment*Pi*a^2*ell*hbar^2/(m*rhoW)`); `WwallClosed` is the closed scalar. Cancellation of `IfMoment`, `a^2`, `m`, `Pi` is the substantive content.
  - `Xi - WwallAssembled`: `Xi` is `(V0/ell)^2 * Nphiphi/Hw * L^2 / Tx` (with `1/Hw` inlined) while `WwallAssembled` is `4 pi a^2 L^2 (IfMoment/Hw) V0^2 / (Tx ell)` — structurally different orderings of the same primitives whose equality requires `(V0/ell)^2 * 4 pi a^2 ell IfMoment` to collapse to `4 pi a^2 IfMoment V0^2 / ell`. That residual reduction is the genuine algebraic content; not tautological.
- Output text file `mathematica/output/...stage070...txt` shows all three `expectZero` checks PASSing post-fix:
  - `kappa - expected = 0` → PASS
  - `W_wall - expected = 0` → PASS
  - `Xi - W_wall = 0` → PASS
- Output mtime (May 22 20:09) is newer than script mtime (May 22 20:08), confirming the output was regenerated after Codex's edit.
- Symbol assumptions (positivity / reals) are intact at lines 28-31.

The Mathematica `expectZero` label is "kappa - expected" rather than the directive's example "kappa: assembled - closed", but this is purely descriptive and immaterial to the verification content.

## Exec log assessment

**SymPy:** exit=n/a. No `stage_070_sympy.log` was captured by the orchestrator (the SymPy script was not edited, so re-running it is not required). The saved `scripts/output/...stage070...sympy_audit.txt` was regenerated at May 22 20:09 and shows `expect_zero` lines printing `= 0` for all three checks; the SymPy script raises on non-zero residuals, so the absence of an `AssertionError` confirms exit 0:
- `kappa - expected = 0`
- `W_wall - expected = 0`
- `Xi - W_wall = 0`

**Mathematica:** exit=n/a. No `stage_070_mathematica.log` was captured by the orchestrator. The saved `mathematica/output/...stage070...mathematica_audit.txt` was regenerated at May 22 20:09 (newer than the script's 20:08 mtime). Notable lines from the captured output:
- `kappa - expected = 0` → `PASS: kappa - expected`
- `W_wall - expected = 0` → `PASS: W_wall - expected`
- `Xi - W_wall = 0` → `PASS: Xi - W_wall`
- `Exit[0]` follows the theorem-ledger banner. The `expectZero` helper calls `Exit[1]` on non-zero residuals, so all-PASS implies clean exit.

**Output freshness:** confirmed. Mathematica script mtime is 2026-05-22 20:08; Mathematica output mtime is 2026-05-22 20:09. SymPy script mtime is unchanged (Apr 1) but its output mtime is 2026-05-22 20:09, indicating it was rerun for cross-check. Both outputs post-date the patch.

## Material-change assessment

`material_change`: false.

The edit only changes the algebraic path the Mathematica script takes; the values produced are unchanged (and exactly match the SymPy output: `kappa = L^2*(Ig/(ell^2*IfMoment) + 4*cSw^2*m^2/hbar^2)`, `W_wall = 4*L^2*rhoW^2*V0^2/(cSw^2*ell^2*hbar^2)`, `Xi = W_wall`). The original auditor explicitly noted: "fixing the Mathematica derivation will not change the values of `kappa`, `W_wall`, or `Xi`; it will only re-establish that the agreement is genuine." No downstream unit can depend on a value that did not change.

## Side observations (non-blocking)

- The label string passed to `expectZero` for the `kappa` check is "kappa - expected" while the directive suggested "kappa: assembled - closed". The new closed-form variable is named `kappaClosed` (not `kappaExpected`), so the label is slightly stale relative to the variable name. Purely cosmetic; does not affect the assertion.
- The unused symbol `Hw` is also listed in the `Clear[...]` and `$Assumptions` block (line 28-31). After the rewrite `Hw` is still used (line 42), so this is fine; just noting that the symbol set is unchanged from the pre-fix script.
- The banner still says "STAGE 053" while file names/SymPy say "STAGE 53"/"stage070". This is descriptive text only and was already that way pre-fix; not a regression.

## Verdict justification

The single finding F1 is fully resolved: the Mathematica script now derives `W_wall` and `Xi` by inlining `1/Hw` directly into the primitive product rather than reusing the SymPy-mirror named intermediates `J_1`, `gphi`, `I_1`; it verifies all three claims against closed forms written directly from primitives; and the regenerated output shows all three `expectZero` checks PASS. The diff matches the directive line-for-line, no collateral edits, no helper/banner changes, no `Exit[0]` change. Output freshness confirmed. No regressions. No new findings.
