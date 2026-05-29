---
unit_id: 111
batch: IV.2
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 111

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes outside the scope below. Do NOT touch `paper/`, `notes/`, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage111_mixed_sidechannel_pole_mathematica_audit.wl:26-66`

**Issue:**
The Mathematica script reproduces the SymPy script's choreography line-by-line: it builds `lambdaMix` by the same subtraction, extracts the same four `Coefficient[lambdaMix, z, n]` quantities, solves `-l2/l0 == 1/9` then `(l2^2/l0^2 - l4/l0) /. kappa -> kappaMatch == 4/81` in the same order, then forms `chiMix = (-l5/l0)/(1/27)` and its `Series` linearization. The second-engine policy requires that the `.wl` derive at least one of the four paper-side identities by an independent algebraic path, so a choreography error in `.py` cannot silently pass in `.wl`.

**Required change:**

Add an independent cross-check of `chiMix` to the Mathematica script. Insert a new block immediately after the existing `chiMix` line (currently L52) that recomputes `chi_Q^mix` by a path that does not go through `l5`. Then add a fifth `expectZero` comparing the two routes.

Specifically, after line 52 (the existing `chiMix = ...` line), add the following block (keep all existing lines):

```mathematica
(* Independent re-derivation of chi_Q^mix: bypass the L0/L5 extraction *)
(* and compute directly from the geometric-series form of the pole. *)
poleSeries = Series[sigma/(1 - kappa*z^2 - I*gamma*z^5), {z, 0, 5}];
imagPart5 = Coefficient[Normal[poleSeries], z, 5]/I;
chiMixAlt = FullSimplify[
  27 (1/9 - imagPart5)/(3 + sigma),
  Assumptions -> $Assumptions
];
Print["chi_Q^mix (independent route) = ", fmt[chiMixAlt]];
```

Then, after the existing four `expectZero` calls (lines 60-63), add:

```mathematica
expectZero["chi_Q^mix routes agree", chiMix - chiMixAlt];
```

This second route computes `chi_Q^mix = -27 L5 / L0` via the substitutions `L0 = -3 - sigma`, `L5 = 1/9 - imagPart5_of_pole_alone`, without first forming `lambdaMix` or running `Coefficient` on it for `z^5`. Algebraically it must reduce to the same closed form. The new `expectZero` will pass iff both routes give the same result.

Do NOT delete or weaken the existing four `expectZero` calls or the existing `chiMix` definition.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 111` and confirm:
1. The new lines appear at the locations specified.
2. The new transcript line `chi_Q^mix routes agree` prints with residual `0`.
3. `expectZero["chi_Q^mix routes agree", ...]` reports `PASS`.
4. Script still exits 0 with all five `expectZero` calls passing.
