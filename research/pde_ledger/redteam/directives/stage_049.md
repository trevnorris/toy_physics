---
unit_id: 049
batch: III.2
created_at: 2026-05-26T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T02:54:44-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 049

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:52-53`

**Issue:**
The v1 F2 fix replaced both the integrator output and the closed-form target with `FullSimplify[Integrate[chiN, {s, 0, l}], ...]` calls, so the subsequent `expectZero["uniform overlap integral", overlapFromIntegral - overlapFormula]` (line 56) now compares two integrator outputs that differ only by whether explicit local `Assumptions` are passed. Under `$Assumptions` (declared at lines 39-43 with `n` integer, `n >= 0`, `l > 0`), both calls reduce to the same form; the assertion is tautological. The notes' boxed closed form `I_n = sqrt(2L)/((n+1/2) pi)` is no longer independently anchored on the Mathematica side, although SymPy still verifies it at lines 70-74. Restore an independent closed-form target without re-introducing the deleted `uniformDnOverlap` helper (which would re-create the structural mirror v1 F2 removed).

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`, replace lines 52-53

Before:
```wolfram
overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}]];
overlapFormula = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
```

After:
```wolfram
overlapFromIntegral = FullSimplify[Integrate[chiN, {s, 0, l}], Assumptions -> Element[n, Integers] && n >= 0 && l > 0];
overlapFormula = Sqrt[2 l]/((n + 1/2) Pi);
```

Leave all other lines unchanged, including the surrounding `Print` statements at lines 54-55, the `expectZero[...]` call at line 56, the `chiN` definition at line 46, and every line from line 58 onward.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 049` and confirm: (a) `overlapFromIntegral` is the integrator call with explicit `Assumptions`, (b) `overlapFormula` is the bare closed-form expression `Sqrt[2 l]/((n + 1/2) Pi)` and no longer a second `Integrate[...]` call, (c) the saved output still shows `uniform overlap integral = 0` with `PASS`, and (d) the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
- summary: Restored the Mathematica overlap target to the independent closed form while keeping the direct integral under explicit assumptions.
- deviation: none
