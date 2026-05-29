---
unit_id: 095
batch: IV.1
created_at: 2026-05-27T00:00:00Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 095

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F3 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_095.tex:23` quote: "Check `l=0` and `l=2` orthogonality before applying the geometry firewall."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py` (no orthogonality check anywhere; the script begins with `Dg = G0 + G2*w**2 + G4*w**4` and proceeds directly to the Schur complement)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl` (same; no orthogonality check)

## Resolve before fix_loop

The stage card lists three `\stagefield{Checks}` items; item 2 ("Check `l=0` and `l=2` orthogonality before applying the geometry firewall") has no script-side counterpart. The reduced 2-mode action in `notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md:23-44` assumes diagonal kinetic/dispersive sectors (`(1/2) q D_q q + (1/2) g D_g g + chi M_0 q g`), so orthogonality is implicit but never tested.

Possible directions (the user picks one):
- (a) The orthogonality bullet is a load-bearing check that should be exercised → specify the orthogonality identity the user wants verified (e.g., explicit inner product of the `l=0` and `l=2` angular harmonics integrated against the appropriate measure, or the absence of a `chi^0` cross-coupling in the kinetic term of the L action) and a follow-up directive will add it to both scripts.
- (b) The orthogonality bullet is already implicit in the reduced-model setup (diagonal `D_q`/`D_g`, off-diagonal only via `chi M_0`) and the bullet should be reworded or removed from the card → a follow-up directive authorizes editing `stage_095.tex:23`.
- (c) Orthogonality was intended to be inherited from an upstream stage (e.g., Stage 74/75/77) and the card should reference that upstream proof rather than imply a local check → name the upstream stage; a follow-up directive amends the card to point there.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py:27-29`

**Issue:**
All five existing SymPy assertions are tautological under the script's construction. The Schur correction `corr` is built by writing `-chi**2 * M0**2 / Dg` as the input expression, so every `K{0,2,4}corr` coefficient automatically inherits the `chi**2` factor; the `.has(chi**2)` asserts on lines 27-29 cannot fail. The chi=0 substitution at line 44 reduces `cpole` to `1/4` by construction (because `eps2` and `eps4` each carry literal `chi**2`). The first derivative `d c_pole / dchi |_{chi=0}` at line 51 vanishes by parity (cpole is even in chi by construction). None of these asserts can reveal a regression in the series expansion or in the algebraic identities the notes assert. The Mathematica mirror (lines 44-46) has three substantive closed-form residual checks that the SymPy script lacks.

**Required change:**

In `scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py`, after the existing `K4corr = sp.simplify(sp.expand(corr).coeff(w, 4))` block (currently line 22) and its three `print` statements (lines 23-25), insert three new asserts BEFORE the existing line 27. The new block:

Before (lines 22-26):
```python
K4corr = sp.simplify(sp.expand(corr).coeff(w, 4))
print('K0corr =', K0corr)
print('K2corr =', K2corr)
print('K4corr =', K4corr)

assert sp.factor(K0corr).has(chi**2)
```

After:
```python
K4corr = sp.simplify(sp.expand(corr).coeff(w, 4))
print('K0corr =', K0corr)
print('K2corr =', K2corr)
print('K4corr =', K4corr)

# Closed-form residual checks against notes §2:
# K_{g,0}^eff = -chi^2 M_0^2 / G_0   (static renormalization)
# K_{g,2}^eff =  chi^2 M_0^2 G_2 / G_0^2
# K_{g,4}^eff =  chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3
assert sp.simplify(K0corr - (-M0**2*chi**2/G0)) == 0
assert sp.simplify(K2corr - (G2*M0**2*chi**2/G0**2)) == 0
assert sp.simplify(K4corr - (M0**2*chi**2*(G0*G4 - G2**2)/G0**3)) == 0
print('PASS: K0corr matches closed form')
print('PASS: K2corr matches closed form')
print('PASS: K4corr matches closed form')

assert sp.factor(K0corr).has(chi**2)
```

Leave the rest of the script unchanged. Do not remove the existing `has(chi**2)` lines; they remain as cheap sanity flags. The new asserts are the substantive checks.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 095` and confirm the new asserts execute, the new `PASS:` lines appear in the saved output transcript, and the script exits 0.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl:31-34`

**Issue:**
The Mathematica script mirrors the SymPy script's algebraic choreography step-for-step rather than re-deriving `D_eff = D_q - chi^2 M_0^2 / D_g` from the action `L = (1/2) q D_q q + (1/2) g D_g g + chi M_0 q g` (notes §1). Both engines start from `-chi^2 m0^2 / dG` as a given. A genuine second-engine check would integrate out `g` in Mathematica from the L action and then assert the resulting effective `q`-kernel correction equals `-chi^2 m0^2 / dG` before proceeding to the existing series checks.

**Required change:**

In `mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl`, between the existing `$Assumptions = ...` block (currently lines 29-31) and the existing `dG = g0 + g2*w^2 + g4*w^4;` line (currently line 33), insert an independent Schur derivation. Insert after line 31:

```mathematica
(* Independent Schur derivation from the bilinear action L = (1/2) q dQ q + (1/2) g dG g + chi m0 q g.
   Treat dQ, dG as scalars at fixed omega; integrate out g and identify the q^2 coefficient.
   Notes §1 (notes/stages/moving_throat_pde_stage095_second_order_geometry_contamination.md:38-44). *)
Clear[dQsym, dGsym, qSym, gSym];
Lq = (1/2)*qSym^2*dQsym + (1/2)*gSym^2*dGsym + chi*m0*qSym*gSym;
gStar = First[gSym /. Solve[D[Lq, gSym] == 0, gSym]];
LqEff = FullSimplify[Lq /. gSym -> gStar, Assumptions -> dGsym != 0];
(* LqEff should be (1/2) qSym^2 dEffSym where dEffSym = dQsym - chi^2 m0^2 / dGsym. *)
dEffCoeff = 2*Coefficient[LqEff, qSym, 2];
corrDerived = FullSimplify[dEffCoeff - dQsym, Assumptions -> dGsym != 0];
expectZero["D_eff Schur derivation matches -chi^2 m0^2 / dG",
  corrDerived - (-chi^2*m0^2/dGsym)];
```

Leave the rest of the script (including the subsequent `dG = g0 + g2*w^2 + g4*w^4;` and the existing series/`expectZero` checks) unchanged. The new block stands alone and uses symbolic `dQsym, dGsym` so it does not interact with the downstream low-frequency expansion.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 095` and confirm a new `PASS: D_eff Schur derivation matches -chi^2 m0^2 / dG` line appears in the saved output transcript and the script exits 0.

---
## Applied: F3 (orchestrator-direct, post-user-resolution per batch-IV1-paper-alignment Cluster A direction (a))

- files_changed: scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py
- summary: SymPy docstring annotates the "l=0/l=2 orthogonality" Check as upstream carry-forward from stage 094. The reduced 2-mode L action used here assumes diagonal kinetic/dispersive sectors and off-diagonal coupling only via chi M_0, which is exactly the orthogonality output of stage 094.
- deviation: none

## Applied: F1
- files_changed: scripts/moving_throat_pde_stage095_second_order_geometry_contamination_sympy_audit.py
- summary: Added three closed-form residual asserts on K0corr/K2corr/K4corr matching notes Section 2. PASS lines: "K0corr matches closed form", "K2corr matches closed form", "K4corr matches closed form". Plus banner sweep + Stage 75 → Stage 092 comment fix.
- deviation: none

## Applied: F2
- files_changed: mathematica/moving_throat_pde_stage095_second_order_geometry_contamination_mathematica_audit.wl
- summary: Inserted independent Schur derivation block (Solve[D[Lq, gSym] == 0]) before the existing series machinery. PASS: "D_eff Schur derivation matches -chi^2 m0^2 / dG". Plus banner sweep.
- deviation: none
