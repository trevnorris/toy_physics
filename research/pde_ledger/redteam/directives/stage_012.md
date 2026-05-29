---
unit_id: 012
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 012 (v2 paper-grounded re-audit)

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl:50-66`

**Issue:** The "M1" block defines `Z0form`, `Z2form`, `Z4form`, `N0form`, `N2form`, `N4form` on lines 41-48 and then on lines 50-66 calls `expectZero[..., Z0form - Q/Delta]`, `expectZero[..., Z2form - (Q S2 - H Delta)/Delta^2]`, etc. — each residual is identically zero by construction because the RHS of each subtraction is the same expression bound to the LHS form symbol on lines 41-48. The transcript at `mathematica/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.txt` lines 2-13 confirms each residual prints as `0`. Six PASS lines verify nothing beyond `FullSimplify[0] === 0`.

**Required change:**

Delete the six tautological `expectZero` calls and replace with a single labeled Print group that documents the assumed primitive forms as upstream carry-forwards (from Stage 4 / Stage 5, per the SymPy docstring line 55).

Before (lines 50-66):

```mathematica
expectZero["M1 Z0 primitive one-port", Z0form - Q/Delta];
expectZero["M1 Z2 primitive one-port", Z2form - (Q S2 - H Delta)/Delta^2];
expectZero[
  "M1 Z4 primitive one-port",
  Z4form - (Q (S2^2 - Delta) - S2 H Delta)/Delta^3
];
expectZero["M1 N0 primitive one-port", N0form - P^2/Delta^2];
expectZero[
  "M1 N2 primitive one-port",
  N2form - 2 P (P S2 - Delta Gw)/Delta^3
];
expectZero[
  "M1 N4 primitive one-port",
  N4form -
    (Delta^2 Gw^2 - 2 Delta P^2 - 4 Delta P S2 Gw +
        3 P^2 S2^2)/Delta^4
];
```

After (replacing those exact 17 lines):

```mathematica
Print["M1 primitive one-port forms (carried from Stage 4 / Stage 5):"];
Print["  Z0 = ", fmt[Z0form]];
Print["  Z2 = ", fmt[Z2form]];
Print["  Z4 = ", fmt[Z4form]];
Print["  N0 = ", fmt[N0form]];
Print["  N2 = ", fmt[N2form]];
Print["  N4 = ", fmt[N4form]];
```

Do not touch lines 41-48 (the form definitions are still needed by M2..M9). Do not touch any block from line 68 onward (those are the substantive M2..M9 checks).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 012` and confirm:

1. The script exits 0.
2. The output transcript no longer contains the six `M1 Z0 primitive one-port`, `M1 Z2 primitive one-port`, ..., `M1 N4 primitive one-port` `PASS:` lines.
3. The output transcript contains the new `M1 primitive one-port forms (carried from Stage 4 / Stage 5):` header followed by six `Z0 = ...`, `Z2 = ...`, ..., `N4 = ...` print lines.
4. All M2..M9 PASS lines remain unchanged and the final `STATUS: PASS` line still prints.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl`
- summary: Replaced the tautological M1 residual checks with print-only documentation of the upstream primitive one-port forms.
- deviation: none
