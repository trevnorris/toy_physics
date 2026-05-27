---
unit_id: 108
batch: IV.2
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 108

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage108_robustness_classes.md:84-92` quote:
  > "The deformation preserves the canonical outgoing normalization iff `chi_Q = 1`. The exact condition is `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27`. This includes as special cases: pure scale deformation, pure scale+argument deformation with beta = 1, additive core deformations whose odd slot is locked to the static shift."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:50-70` — Class C is built with `Lambda_add = S*Lambda_out + Sigma_0 + Sigma_2 z^2 + Sigma_4 z^4 + I Sigma_5 z^5` (no `beta` substitution into `Lambda_out`); the preservation-locus check only verifies the `beta = 1` reduction `Sigma_5 = -Sigma_0/27`.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:57-81` — same pattern; `lambdaAdd = sNorm*lambdaOut + ...` without `beta z` substitution, locus check only at `beta = 1`.

## Resolve before fix_loop

The notes box a general `beta`-dependent exact preservation submanifold `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` that the scripts do not verify (they only test the `beta = 1` reduction `Sigma_5 = -Sigma_0/27`). Should the scripts be extended to test the general locus, or should the notes be trimmed to match what the scripts actually verify?

Possible directions (the user picks one):
- (a) Notes/paper are correct, scripts under-verify -> extend both scripts: combine scale+argument with additive core (`lambdaAdd = sNorm * lambdaOut(beta z) + sigma_0 + sigma_2 z^2 + sigma_4 z^4 + I sigma_5 z^5`), re-solve `(sigma_2, sigma_4)` against even-moment match (the solutions will now depend on `beta`), then assert `chi_Q = 1` iff `sigma_5 = S(1 - beta^5)/9 - sigma_0/27`. The current `beta = 1` Class B and Class C blocks remain as sanity reductions.
- (b) Scripts are correct, notes overstate -> downgrade `notes/stages/moving_throat_pde_stage108_robustness_classes.md` "Exact preservation submanifold" section to state only the `beta = 1` special case `Sigma_5 = -Sigma_0/27`, and remove the boxed `Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27` formula. No script change.
- (c) Both are correct but verify different things deliberately (notes are aspirational, scripts are the current verified subset) -> add a `\stagefield{Scope}` note in the card and the notes flagging that the general locus is unverified at this audit unit; no script change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:55`

**Issue:** The line `expectZero["chi_arg(beta=1) - 1", chiArg /. beta -> 1 - 1];` is parsed by Mathematica as `chiArg /. (beta -> (1 - 1))` = `chiArg /. (beta -> 0)`, because `Plus`/`Subtract` bind tighter than `Rule` and `ReplaceAll`. Since `chiArg = beta^5`, the substituted residual is `0^5 = 0` regardless of whether `chi_arg(beta=1)` actually equals 1. The check passes by accident, not by computing the intended quantity. The SymPy counterpart on L48 uses `.subs(beta, 1) - 1` and is correct.

**Required change:**
Edit `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl` line 55.

Before:
```
expectZero["chi_arg(beta=1) - 1", chiArg /. beta -> 1 - 1];
```

After:
```
expectZero["chi_arg(beta=1) - 1", (chiArg /. beta -> 1) - 1];
```

This is the only change. Do not touch surrounding lines.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 108` and confirm:
1. The transcript still prints `chi_arg(beta=1) - 1 = 0` and `PASS: chi_arg(beta=1) - 1`.
2. The script exits 0.
3. Line 55 of the `.wl` file now reads `expectZero["chi_arg(beta=1) - 1", (chiArg /. beta -> 1) - 1];` exactly.

## F3 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py:25`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl:26`

**Issue:** Both scripts banner this unit as "STAGE 91" / "STAGE 091", but this is unit 108. The mtime check confirms the saved outputs are fresh; the issue is the *banner string* in the script source, which propagates to the transcripts and misleads downstream grepping.

**Required change:**
Edit `scripts/moving_throat_pde_stage108_robustness_classes_sympy_audit.py` line 25.

Before:
```
banner('STAGE 91 — ROBUSTNESS CLASSES FOR chi_Q')
```

After:
```
banner('STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q')
```

Edit `mathematica/moving_throat_pde_stage108_robustness_classes_mathematica_audit.wl` line 26.

Before:
```
banner["STAGE 091 — ROBUSTNESS CLASSES FOR chi_Q"];
```

After:
```
banner["STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q"];
```

Do not touch any other line. Preserve the em-dash and surrounding whitespace exactly.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 108` and `redteam exec-mathematica 108` and confirm both transcripts now print `STAGE 108 — ROBUSTNESS CLASSES FOR chi_Q` in the banner section. Scripts must continue to exit 0.
