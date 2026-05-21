---
unit_id: 009
batch: I.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-21T12:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 009

## Per-finding outcomes

### F1 — missing_verification_script

**Classification:** resolved

**What changed:**
A new independent Mathematica audit script was created at
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl`
(1-103). It uses Mathematica idioms (`Integrate`, `Series[..., {x, Infinity, n}]`,
`Normal`, `FullSimplify`, `Assuming`) and a custom `assertZero` guard that calls
`Exit[1]` on failure. It covers all five claim-manifest items (M1a–M1c, M2a–M2b,
M3a–M3b, M4a–M4b, M5a–M5b) and terminates with `Print["STATUS: PASS"]; Exit[0]`.

**Assessment:**
The file lives in the canonical `mathematica/` engine directory (the orchestrator
note confirms it was relocated there from an initial misplacement). The script is
not a transliteration of the SymPy script: it ratios the same polynomial weights
directly via `Series`, rather than calling `sp.integrate` on `W_conc*Z_conc`. For
M5b (the mouth Gaussian) it goes through `Integrate` → `FullSimplify` →
`/. ell -> 1/r` → `Series[..., {r, Infinity, 7}]` — Mathematica produces the
`Erfc` closed form via `Integrate`, not by retyping the formula, satisfying F1's
"derive the erfc form" requirement. Exec log
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_009_mathematica.log`
shows all 11 OK lines and `STATUS: PASS`, exit_code 0.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Lines 190–221 of
`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.py`
replace the four tautological `subs(H_sym, Z_conc)` / `subs(S_sym, C*Z_conc)`
checks with a perturbative profile-alignment block:
- `H_pert = Z_conc + eps*(h0 + h1*w + h2*w**2/2)`,
  `S_pert = C*Z_conc + eps*(s0 + s1*w + s2*w**2/2)`,
- `IH_pert`, `IS_pert` are SymPy-*evaluated* integrals (no
  `sp.Integral`-placeholder substitution),
- two `eps=0` cancellation asserts plus two `sp.series(...,eps,0,2)`
  first-order corrections that compare to `xi - eps*xi*IDh/IZ_conc` and
  `C*mu0 + eps*mu0*IDs/IZ_conc`.

**Assessment:**
The new assertions are genuinely discriminating. The eps^0 checks use
`IH_pert.subs(eps,0)` where `IH_pert` is already an evaluated integral, so the
denominator is the actually-computed `IZ_conc` (not a symbolic `sp.Integral`
placeholder that hides under `.subs(H_sym, Z_conc)`). The eps^1 checks compare
SymPy's first-order Taylor coefficient against a hand-derived ratio of distinct
computed integrals (`IDh/IZ_conc`, `IDs/IZ_conc`) — a typo in any of the
expansions or in the integrate path would break this. The required label changes
("H=Z effective gauge (eps=0 cancellation)" etc.) appear; the old "H=Z effective
gauge"/"S=CZ effective coupling" labels are gone (confirmed in the diff). The
SymPy output txt (mtime newer than the script) shows the perturbative-block
print section with finite expanded expressions, so the script ran end-to-end and
all asserts in this block passed.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
Six new `assert_zero(...)` calls added at the exact directive-specified spots:
- line 111: "half-line Q expansion" — `avg_Q_half` literal match,
- line 151–152: "symmetric mu_eff series",
- line 154–155: "symmetric xi_eff series",
- line 161–162: "mouth mu_eff series",
- line 164–165: "mouth xi_eff series",
- line 231–232: "symmetric Gaussian asymptotic literal".

**Assessment:**
Each new assertion compares a `sp.series(...).removeO()` (or `sp.integrate`
result) against the hand-derived canonical form given in the directive — non-
tautological by construction (the LHS is the SymPy-computed series; the RHS is
the predicted algebraic closed form). The targets exactly match the directive's
six literal expansions. No surplus changes were made. The output txt's expansion
strings ("series = 1 - sigma**2/lambda**2 + 3*sigma**4/(2*lambda**4)" on line
108; the μ_eff/ξ_eff symmetric/mouth expansions on lines 79–84) are consistent
with the asserted literals, and the script reaches `STATUS: PASS` at the bottom
of the output, confirming all six asserts passed.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
Lines 236–246 of the SymPy script replace the hand-typed
`IWZ_mouth_gauss_r = sp.sqrt(sp.pi) * lam * r * sp.erfc(lam*r/2) * sp.exp(lam**2*r**2/4) / 2`
with a derived form
`IWZ_mouth_gauss_r = sp.simplify(IWZ_mouth_gauss.rewrite(sp.erfc).subs(ell, 1/r))`,
and adds an `assert_zero("mouth Gaussian integral equals erfc closed form",
sp.simplify((IWZ_mouth_gauss - IWZ_mouth_gauss_erfc).rewrite(sp.erfc)))` guard.

**Assessment:**
This is the structural fix the finding asked for. The asymptotic series is now
the `series(...)` expansion of the SymPy-`integrate`-evaluated `IWZ_mouth_gauss`,
not of a literal. The new "mouth Gaussian integral equals erfc closed form"
guard verifies that the SymPy integral and the typed erfc closed form differ
by zero, so the typed closed form is now a tested *consequence* of `sp.integrate`
rather than a pre-baked answer. Codex's documented deviation (using
`.rewrite(sp.erfc)` on `IWZ_mouth_gauss` because SymPy's raw `erf` form does
not series-expand at infinity) is reasonable and was disclosed in the
`## Applied: F4` block — the rewrite is a purely cosmetic reformulation of an
already-evaluated integral, not an injection of new information. The existing
A9 ("mouth Gaussian asymptotic from erfc closed form") and A10 ("mouth Gaussian
asymptotic from Taylor integration") remain in place. The output txt shows
`<Z> = sqrt(pi)*lambda*(1 - erf(lambda/(2*ell)))*exp(lambda**2/(4*ell**2))/(2*ell)`
and `series = -120*ell**6/lambda**6 + 12*ell**4/lambda**4 - 2*ell**2/lambda**2 + 1`,
both matching the expected forms.

## Exec log assessment

**SymPy:** exit=n/a. Notable lines:
The file
`/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_009_sympy.log`
does NOT exist — the orchestrator did not capture a SymPy log for this stage.
However, the SymPy output text
`/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage009_projected_maxwell_near_throat_sympy_audit.txt`
has mtime newer than the script (1779384382 > 1779383480) and ends with
`STATUS: PASS` on line 142, which only prints if all `assert_zero(...)` calls
succeeded. Per the user's note that "both engines exit 0; orchestrator
confirmed," I treat this as a passing SymPy run with the log file simply
unsaved.

**Mathematica:** exit=0. Notable lines:
```
OK: M1a half-line IBP recombination residual = 0
OK: M1b half-line dQ expansion residual = 0
OK: M2a Gaussian even-kernel Q residual = 0
OK: M5b mouth Gaussian asymptotic residual = 0
STATUS: PASS
```
All 11 expected `OK: ` labels (M1a, M1b, M1c, M2a, M2b, M3a, M3b, M4a, M4b, M5a,
M5b) appear; exit_code = 0.

**Output freshness:**
- sympy script mtime: 1779383480 (2026-05-21T11:11)
- sympy output txt mtime: 1779384382 (2026-05-21T11:26) — newer, OK.
- mathematica wl mtime: 1779383513 (2026-05-21T11:11)
- mathematica exec log mtime: 1779383803 (2026-05-21T11:16) — newer, OK.
A standalone `..._mathematica_audit.txt` in `scripts/output/` was not generated
(only the exec log captures the Mathematica output), but this is consistent
with how other stages in this batch are handled and the exec log itself
serves as the durable artifact.

## Material-change assessment

`material_change`: false.

The edits change *which assertions* the script enforces and how the mouth
Gaussian asymptotic series is *derived*, but the surviving quantitative results
in the SymPy output txt (the symmetric and mouth Gaussian series, the
μ_eff/ξ_eff expansions, the half-line Q and dQ expansions, and the section-5
closed-form `sqrt(pi)*lambda*(1 - erf(...))*exp(...)/(2*ell)`) are bit-for-bit
identical to the pre-fix outputs. No downstream-visible numeric or symbolic
constant changed; only the verification rigor improved. The newly added
Mathematica script is independent of the SymPy script's outputs and does not
feed any downstream unit.

## Side observations (non-blocking)

- The orchestrator captured a Mathematica exec log but no `mathematica_audit.txt`
  in `scripts/output/`. Stage 008 has both; stage 009 has only the log. This is
  not in scope for any F-finding here, but a future hygiene pass may want to
  align output-capture policy across units.
- The orphaned-codex incident (session id not recorded for stage 009) had no
  observable effect on the edits themselves: the directive's `## Applied: F1–F4`
  blocks are present and the file changes match the directive byte-level. No
  rework signal needed.
- F2's substantive check passes specifically because `IH_pert`/`IS_pert` are
  computed via `sp.integrate` (yielding evaluated polynomial-in-ell results)
  before any `eps=0` substitution. If a future refactor turned them back into
  `sp.Integral(...)` placeholders, the check would silently regress to a
  symbol-substitution tautology again. Worth flagging in the script as a
  comment in a future pass (out of scope here).

## Verdict justification

All four findings (F1, F2, F3, F4) are resolved. The Mathematica script exists,
covers all five manifest claims, and exits 0 with 11 `OK` lines. The SymPy
script's tautological profile-alignment block was replaced with a perturbative
test whose first-order coefficients compare distinct evaluated integrals; the
six previously-printed series now have concrete `assert_zero` anchors; the
mouth-Gaussian asymptotic is derived from the SymPy-evaluated integral and
guarded against the typed erfc closed form. No regressions are visible in the
diff, no downstream-visible constants changed, and the orchestrator confirms
both engines exit 0. Verdict: `verified`.
