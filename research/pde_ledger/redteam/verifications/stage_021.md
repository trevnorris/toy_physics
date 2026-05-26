---
unit_id: 021
batch: I.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T22:30:00Z
verdict: needs_rework
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 0
findings_total: 1
material_change: false
---

# Verification — unit 021

## Per-finding outcomes

### F1 — paper_misalignment (script_missing_paper_claim)

**Classification:** partial

**What changed:**

Codex applied the directive's option (a) (Q6=a / ADD) in both engines but deviated from the directive's suggested assertion shape in a way that makes the new check tautological.

- sympy `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.py:237-245` — declared local `a, c_s` and added:
  ```python
  Dcorr = sp.simplify(-I * Gamma_port * omega**5 * N0)
  expect_zero(
      "delta D_wall^(odd) composed l=2 coefficient",
      Dcorr.subs(Gamma_port, a**5 / (27 * c_s**5))
      - (-I * N0 * a**5 * omega**5 / (27 * c_s**5)),
  )
  ```
- mathematica `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:123-139` — declared `radius, cS` and added the mirror:
  ```mathematica
  dCorr = FullSimplify[-I gammaPort omega^5 n0, ...];
  expectZero[
    "delta D_wall^(odd) composed l=2 coefficient",
    (dCorr /. gammaPort -> radius^5/(27 cS^5)) + I n0 radius^5 omega^5/(27 cS^5)
  ];
  ```

The saved outputs show `delta D_wall^(odd) composed l=2 coefficient = 0` in both engines (sympy line 129; mathematica lines 50–51 with explicit `PASS`).

**Assessment:**

The structural intent is right (assertion added in both engines, mirrors aligned, output prints "= 0"), but the assertion is **tautological / vacuous**:

1. **Sympy.** LHS: `Dcorr.subs(Gamma_port, a^5/(27 c_s^5))` is, by inspection of line 240, `sp.simplify(-I * (a^5/(27 c_s^5)) * omega^5 * N0)` — that is, `-I * N0 * a^5 * omega^5 / (27 c_s^5)`. RHS subtracted: `-I * N0 * a^5 * omega^5 / (27 c_s^5)`. Literally the same string of symbols on both sides; the difference is identically zero before any cancellation in `N0`. The assertion exercises nothing — it is just `X − X = 0`.

2. **Mathematica.** Same shape: `(dCorr /. gammaPort -> radius^5/(27 cS^5))` is `-I * (radius^5/(27 cS^5)) * omega^5 * n0`, and the term added is `+I n0 radius^5 omega^5/(27 cS^5)` = the negation. `X + (−X) = 0` trivially.

The directive's recommended fix (directive lines 47–53) was substantively different and non-tautological: it asked Codex to assert against the expanded `(Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² − R²)²` closed form for N(0). That assertion would have exercised the symbolic equality between (i) the Schur-quotient-derived N0 and (ii) the closed-form positive-square representation — useful cross-bracing between Section III's derivation and the paper's stated form. Codex substituted the symbol `N0` itself on the RHS, which collapses the assertion to a no-op.

In addition, neither script ties `Gamma_port`/`gammaPort` (Section III symbol) to `Gamma5_port`/`gamma5Port` (Section IV symbol, independently shown equal to `a^5/(27 c_s^5)` in assertion A10/B10). The constant `a^5/(27 c_s^5)` is hand-substituted in Section III, not pulled from the Section IV result. So even if the assertion were non-tautological in N0, it would still fail to bridge the two pieces the original finding said were never composed.

A reader/future-edit safety check: if a later edit broke either (a) the Section IV derivation of `Gamma5_port = a^5/(27 c_s^5)`, or (b) the Section III closed form for N0 in terms of `(Ω_A² g_W + R g_A)² / (...)²`, the current "composed" assertion would still pass `0 = 0` and the audit would not catch it. The original finding was precisely about this risk; the fix as written does not remove it.

Verdict for F1: **partial** — assertion is in place, exec passes, mirror parity OK, but the assertion is tautological and does not actually compose the two independently-proven pieces, so the paper deliverable still lacks a meaningful script-side check.

## Exec log assessment

**SymPy:** exit=0 (per orchestrator; canonical output `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage021_reduced_one_port_normal_form_sympy_audit.txt` regenerated 2026-05-25 21:49 after script mtime 21:42). Notable lines:
- `delta D_wall^(odd) composed l=2 coefficient = 0` (line 129)
- `N(omega) compact formula = 0` (line 120)
- `N(0) positive-square form = 0` (line 121)

**Mathematica:** exit=0 (canonical output `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt` regenerated 2026-05-25 21:50 after script mtime 21:43). Notable lines:
- `delta D_wall^(odd) composed l=2 coefficient = 0` and `PASS: delta D_wall^(odd) composed l=2 coefficient` (lines 50–51)

**Output freshness:** confirmed — sympy output mtime 21:49 > script mtime 21:42; mathematica output mtime 21:50 > script mtime 21:43. Both fresh.

Note: the orchestrator-provided exec logs at `redteam/exec_logs/stage_021_*.log` are stale (May 21, pre-fix). I used the canonical `scripts/output/` and `mathematica/output/` files as directed by the prompt.

## Material-change assessment

`material_change`: false.

The only code change is an added assertion (which happens to be tautological). No derived expression downstream of stage 021 is altered — N0, Sigma_cons, transfer factor N(omega), Gamma5_port, and the Y2_hat fingerprint are all unchanged in both engines. Downstream units that depend on the value `a^5/(27 c_s^5)` or N_2(0) see the same numerical/symbolic results as before.

## Side observations (non-blocking)

- The auditor's original report also flagged cosmetic stale labels ("STAGE 004 — …", "FINAL STAGE-4 LEDGER") in the script banners and module docstrings. These remain unchanged in the diff. Not part of any finding; flagging for awareness only.
- Codex correctly declared `radius > 0 && cS > 0` in the `$Assumptions` update on .wl lines 124–126 and re-Cleared the new symbols (`Clear[piOut, gammaPort, radius, cS]`). Mirror-side hygiene is clean.
- Sympy declares `a, c_s` locally inside `outgoing_mixed_dressing_audit()` (line 239), not connected to the `a, c_s` declared inside `outgoing_l2_fingerprint_audit()` (line 264). These are independent `sp.symbols` calls but produce identically-named symbols, so substitutions still work as expected. No bug, but the lack of shared symbol/import means there is literally no Python-level cross-reference between the two sections.

## Delta directive (iteration 2)

### F1 — paper_misalignment (still partial)

The current assertion in both scripts is `X − X = 0` and exercises nothing. Replace it with a non-tautological composition that ties Section III's N0 to the closed-form positive-square representation, exactly as the original directive (lines 47–53) recommended.

**Sympy** — replace lines 240–245 with:

```python
expect_zero(
    "delta D_2^(odd) composed (paper eq:app-stage021-wall-odd)",
    Dcorr.subs(Gamma_port, a**5 / (27 * c_s**5))
    - (-I * (OA**2 * gW + R * gA)**2 / (OA**2 * OW**2 - R**2)**2
        * a**5 / (27 * c_s**5) * omega**5),
)
```

This exercises the equality `N0 == (OA²gW + R gA)² / (OA²OW² − R²)²` as part of the assertion (i.e., subtracting two formulas that are equal only because the Schur-derived N0 simplifies to the positive-square form). If a future edit broke either side, the assertion would fail to simplify to zero.

**Mathematica** — same shape; replace lines 136–139 with:

```mathematica
expectZero[
  "delta D_2^(odd) composed (paper eq:app-stage021-wall-odd)",
  (dCorr /. gammaPort -> radius^5/(27 cS^5))
   - (-I (oA^2 gW + r gA)^2/(oA^2 oW^2 - r^2)^2
        * radius^5/(27 cS^5) * omega^5)
];
```

After the rewrite, both engines should still print `delta D_2^(odd) composed (paper eq:app-stage021-wall-odd) = 0`, but the zero must result from a non-trivial simplification between `N0` (Schur-quotient rational) and `(OA²gW + R gA)² / (OA²OW² − R²)²` (closed-form positive-square). Verify by inspection that the LHS, before simplification, is not literally the negation of the RHS as printed text.

(Optional, not blocking: also bridge `Gamma_port` (Section III) to `Gamma5_port` (Section IV) by hoisting the Section IV result into Section III's scope and asserting `Gamma_port` was set to it before substitution. The auditor did not require this, but it would close the second half of the original finding's risk.)

## Verdict justification

The finding's intent — give the paper's third Output deliverable a script-side check — was acknowledged and the mechanical motion was made: assertions added in both engines, mirror parity preserved, outputs fresh, scripts still pass. But the chosen assertion form is `(−I N0 a⁵ ω⁵/(27 c_s⁵)) − (−I N0 a⁵ ω⁵/(27 c_s⁵)) = 0`, which is identically zero before any sympy/Mathematica simplification touches it. It does not exercise N0's closed form, does not connect to Section IV's independently-derived `Γ_5^port = a^5/(27 c_s^5)`, and would not catch the exact future-edit regression class the original finding warned about. This is the canonical "rubber-stamp" failure mode the verifier prompt explicitly calls out. Classification: `partial`. Overall: `needs_rework`.
