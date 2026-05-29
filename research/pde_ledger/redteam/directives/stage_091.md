---
unit_id: 091
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 091

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1), do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_091.tex:23` quote: "Check \(l=0\) and \(l=2\) orthogonality before applying the geometry firewall."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation.md` (entire file) — the term "orthogonality" and the labels `l=0` / `l=2` do not appear; the notes work purely with the scalar `Kcons(omega)` ansatz.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py` — no angular-mode decomposition; only the scalar `Kcons(omega)` ansatz.
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl` — same; scalar-only.

## Resolve before fix_loop

The stage-091 card's `\stagefield{Checks}` block lists "Check `l=0` and `l=2` orthogonality before applying the geometry firewall" as a verification item, but neither the stage notes nor the audit scripts touch any `l=0`/`l=2` angular decomposition. Is this checklist item supposed to be (i) a hypothesis presumed established by an upstream stage (so the card should reword or relocate the bullet), or (ii) a stage-091 obligation that the scripts must verify here?

Possible directions (the user picks one):
- (a) Hypothesis-from-upstream → reword the card line (e.g., move to `\stagefield{Inputs}` as a prerequisite carried in, or drop it from `\stagefield{Checks}`). Paper-side edit. No script change. A follow-up directive can patch `stage_091.tex:23`.
- (b) Stage-091 obligation → extend `sympy_audit.py` and `mathematica_audit.wl` with an explicit `l=0`/`l=2` orthogonality check. A follow-up directive will specify the integrand (presumably `<Y_0^0, Y_2^0>` on the sphere, or the equivalent radial-projection inner product the program uses) and the assertion form. Script-side edit only.
- (c) Both card and scripts are correct as-is, and the checklist line refers to an audit performed in a separate upstream unit that downstream consumers can locate via a cross-reference → add a `\stagefield{Inputs}` cross-reference to that upstream unit and reword the bullet. Paper-side edit; no script change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl:68` (insert new block after this line, before line 70's `Print[""];`)

**Issue:** The `.wl` script is structurally a line-by-line port of the `.py` script: same Kcons, same series, same coefficient extraction, same `Solve`-on-branch-identity, same eight assertions in the same order. The second-engine policy requires the Mathematica script to derive the main identity by an independent route, not merely echo the SymPy choreography. The algebra here is single-path, but one cheap diversification is available: verify the bottom-line `Yhat_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2)` identity by direct partial-fraction recombination from `K_geom = 3*K_pole` (without going through `Series` + `Coefficient` + `Solve`).

**Required change:**

Insert the following block into the Mathematica script, between the existing `expectZero["zeta_req - 1/3", ...]` (line 68) and the closing `Print[""];` (line 70):

```mathematica
(* Independent derivation: bypass Series + Solve and recombine directly. *)
(* Start from the branch-identity result K_geom = 3*K_pole and verify Yhat via Together, *)
(* not via series expansion. This gives an independent algebraic path to the same identity. *)
kConsBranchDirect = 3*kPole + kPole/(1 - omega^2/omegaQ^2);
k0BranchDirect = 4*kPole;
yHatRecomb = Together[kConsBranchDirect/k0BranchDirect];
yHatTargetRecomb = Together[3/4 + (1/4)/(1 - omega^2/omegaQ^2)];
expectZero["Yhat partial-fraction recombination", yHatRecomb - yHatTargetRecomb];
```

After applying, the script should still end with `Print[""]; Print["Stage 091 Mathematica audit passed."]; Exit[0];` (lines 70-73 unchanged).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 091` and confirm the new check appears in the saved transcript (look for the line `PASS: Yhat partial-fraction recombination`) AND the script exits 0.

---
## Applied: F1 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster A direction (a))

- files_changed: scripts/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_sympy_audit.py
- summary: SymPy docstring now annotates the "l=0 and l=2 orthogonality" Check as an upstream carry-forward from stage 094 (Isotropic Geometry-Decoupling Theorem). No script-side assertion added; no paper edit.
- deviation: none

## Applied: F2
- files_changed: mathematica/moving_throat_pde_stage091_grouped_p2_static_geometry_derivation_mathematica_audit.wl
- summary: Inserted independent partial-fraction recombination block before final banner (PASS: "Yhat partial-fraction recombination"). Plus banner sweep STAGE 074 → STAGE 091.
- deviation: none
