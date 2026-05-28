---
unit_id: 166
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 166

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond the lines named. Do NOT run python or mathematica. Do NOT touch paper.tex, notes/, or any prose document.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py` (insert after line 50)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl` (insert after line 45)

**Issue:** The paper §2 (notes §2 boxed equations) headline results `δlnρ_w = ½ δlnΘ_w` and `δln a = ½ δlnK_s − ¼ δlnΘ_w` are only `print`ed (sympy 47-48; math 42-43), never asserted against their target forms. The only assertions touching `drho`/`da` are the near-tautological forward re-substitution (sympy 53-56; math 48-51 — re-substituting `sp.solve`/`Solve` output into the same equations, guaranteed zero) and the frozen-wall block (sympy 82-83; math 76-77), which sets `dTheta → 0` and thereby deletes exactly the `dTheta` slopes being claimed. A wrong slope on `drho`'s `½ dTheta` or `da`'s `−¼ dTheta` would pass silently.

**Required change:**

SymPy — after line 50 (the four `print(...)` of the solved forms), before `banner("Forward verification")`, insert:

```python
banner("General inversion forms (paper Sec. 2)")
expect_zero("drho general", sol[drho] - sp.Rational(1, 2) * dTheta)
expect_zero("da general", sol[da] - (sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta))
expect_zero("dcs general", sol[dcs] - (sp.Rational(1, 2) * dKs - sp.Rational(1, 4) * dTheta + sp.Rational(1, 5) * dP))
expect_zero("dZ general", sol[dZ] - (dKq - sp.Rational(2, 5) * dP))
```

Mathematica — after line 45 (the four `Print[...]` of `drhoSol`..`dZSol`), before `banner["Forward verification"]`, insert:

```wolfram
banner["General inversion forms (paper Sec. 2)"];
expectZero["drho general", drhoSol - dTheta/2];
expectZero["da general", daSol - (dKs/2 - dTheta/4)];
expectZero["dcs general", dcsSol - (dKs/2 - dTheta/4 + dP/5)];
expectZero["dZ general", dZSol - (dKq - 2*dP/5)];
```

(The `dcs`/`dZ` general checks are included for completeness; the load-bearing additions are `drho general` and `da general`, which the existing assertions never cover.)

**Claim manifest:**
- M1: `δlnρ_w = ½ δlnΘ_w` (general, with `dTheta` present) — notes §2 boxed eq.
- M2: `δln a = ½ δlnK_s − ¼ δlnΘ_w` (general) — notes §2 boxed eq.
- M3: `δln c_s = ½ δlnK_s − ¼ δlnΘ_w + ⅕ δlnP_0` (general) — notes §2 boxed eq.
- M4: `δlnZ_q = δlnK_q − ⅖ δlnP_0` (general) — notes §2 boxed eq.

**Verification command:**
Verifier runs `redteam exec-sympy 166` and `redteam exec-mathematica 166`; confirms the new lines `drho general = 0`, `da general = 0`, `dcs general = 0`, `dZ general = 0` (sympy) and the matching `PASS:` lines (mathematica) appear, and both scripts exit 0.

## Applied: F1
files_changed:
- scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py
- mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl
summary: Inserted a "General inversion forms (paper Sec. 2)" block in both engines that asserts the four solved drifts against their target slopes (drho=½dTheta, da=½dKs−¼dTheta, dcs=½dKs−¼dTheta+⅕dP, dZ=dKq−⅖dP) before the forward-verification banner.
deviation: none

## F2 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py:33` (banner)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:26` (banner) and an independent cross-check block

**Issue:** The `.wl` is a line-by-line port of the `.py` (same ten symbols, same `eq1..eq4`, same `Solve`/`solve` over the same unknowns and order, same bundle/frozen specializations, same target assertions, same numeric, and the identical wrong `"STAGE 149"` banner). The second engine does not derive the inversion by an independent route. Mitigating factor: this is a canonical 4×4 linear inversion, so the change is minimal.

**Required change:**

1. Banner label fix (both engines): change the literal `"STAGE 149 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS"` to `"STAGE 166 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS"`.
   - `.py:33`: `banner("STAGE 149 — ...")` → `banner("STAGE 166 — ...")`.
   - `.wl:26`: `banner["STAGE 149 — ..."]` → `banner["STAGE 166 — ..."]`.

2. Independent second-engine derivation (Mathematica only). After the existing `Solve`-based block, add a distinct matrix-inverse derivation of the same inversion, so the two engines do not share an algebra path. Insert immediately before `banner["Equivalent full-bundle form with P_0 = N_0 / D_0"]` (currently line 53):

```wolfram
banner["Independent matrix-inverse cross-check"];
(* Forward map M: (drho, da, dcs, dZ) -> (dTheta, dKs, dKq, dP) from eq1..eq4. *)
Mmat = {
  {2, 0, 0, 0},   (* dTheta = 2 drho *)
  {1, 2, 0, 0},   (* dKs    = drho + 2 da *)
  {0, -2, 2, 1},  (* dKq    = -2 da + 2 dcs + dZ *)
  {0, -5, 5, 0}   (* dP     = -5 da + 5 dcs *)
};
inv = Inverse[Mmat];
solVec = inv . {dTheta, dKs, dKq, dP};
expectZero["matrix drho", solVec[[1]] - dTheta/2];
expectZero["matrix da", solVec[[2]] - (dKs/2 - dTheta/4)];
expectZero["matrix dcs", solVec[[3]] - (dKs/2 - dTheta/4 + dP/5)];
expectZero["matrix dZ", solVec[[4]] - (dKq - 2*dP/5)];
(* Round-trip: forward map of the matrix solution recovers the observables. *)
expectZero["matrix round-trip", Mmat . solVec - {dTheta, dKs, dKq, dP}];
```

This derives the inversion via `Inverse[Mmat]` rather than echoing the `Solve` call, giving a genuinely distinct path to M1-M4, and the round-trip line confirms `Mmat` was transcribed correctly from `eq1..eq4`.

**Verification command:**
Verifier runs `redteam exec-mathematica 166` and `redteam exec-sympy 166`; confirms both banners read "STAGE 166", the Mathematica transcript shows the `matrix drho/da/dcs/dZ` and `matrix round-trip` PASS lines, and both scripts exit 0.

## Applied: F2
files_changed:
- scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py
- mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl
summary: Fixed the wrong "STAGE 149" banner to "STAGE 166" in both engines, and added an independent matrix-inverse cross-check block (Inverse[Mmat] derivation of M1–M4 plus a round-trip check) before the full-bundle-form banner in the Mathematica script.
deviation: none

## Orchestrator note (F2 round-trip scalarization)
The directive's prescribed `expectZero["matrix round-trip", Mmat . solVec - {dTheta, dKs, dKq, dP}]` produces a length-4 vector residual; `expectZero` tests `res === 0`, which is False for a list, so the first run FAILed with residual `{0,0,0,0}` (every component already zero). Orchestrator scalarized to `Total[(Mmat . solVec - {dTheta, dKs, dKq, dP})^2]` (zero iff every component is zero). Re-run: 19 PASS, 0 FAIL, exit 0.
