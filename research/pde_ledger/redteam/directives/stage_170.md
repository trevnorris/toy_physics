---
unit_id: 170
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 170

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings (F1), do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" F1 unless a follow-up directive explicitly authorizes a direction.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_170.tex:21` quote: "Check the weak-axisymmetric signature \((1,1/2,-1)\) before reducing grouped defects to a scalar."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage170_linear_grouped_outlet_map.md:385-424` quote (boxed): "\(\delta\kappa_W^{(20)}=\epsilon\,\kappa_1,\ \delta\kappa_W^{(21)}=\frac\epsilon2\,\kappa_1,\ \delta\kappa_W^{(22)}=-\epsilon\,\kappa_1\)" with "\(\kappa_1=\frac{3(1-\sigma_*)}{\sigma_* D_0}\left(D_2^{(1)}+\frac19 D_0^{(1)}\right)\)" and the matching \(\gamma_1\) set; "\((\lambda_{20},\lambda_{21},\lambda_{22})=\left(1,\frac12,-1\right)\)."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py` — no block referencing `lambda`, the `(1,1/2,-1)` signature, `kappa_1`, or `gamma_1`; the only `eps` use is the Sec. 1 linearization (L46-54).
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl` — same absence (no `λ`/signature/κ1/γ1 block).

## Resolve before fix_loop

The stage card's Checks list and notes Section 5 enumerate the weak-axisymmetric signature `(λ_20,λ_21,λ_22)=(1,1/2,-1)` and the two scalar amplitudes `κ1`, `γ1` as a deliverable, but neither script verifies it. Should this be verified in-script, or is it an intentional documentation-only corollary that needs no independent check?

Possible directions (the user picks one):
- (a) A script check is required → add a `expect_zero`/`expectZero` block to BOTH engines verifying, under `δD_{A,n}=ε λ_A D_n^(1)` and `δN_{A,0}=ε λ_A N_0^(1)` with `λ=(1,1/2,-1)`: `δκ_W^(21)=(1/2)δκ_W^(20)`, `δκ_W^(22)=-δκ_W^(20)`, the same for `δγ_W`, and that `κ1 = 3(1-σ)/(σ D0)(D2^(1)+D0^(1)/9)`, `γ1 = -(1-σ)/(9σ N0)(N0^(1)-P0 D0^(1))`. Then re-run sympy+mathematica.
- (b) No script check is needed (Sec. 5 is a trivial linearity corollary; the card's Checks line is documentation of a downstream convention) → no script change; leave a note that Sec. 5 is intentionally uncovered.
- (c) Defer until the second full pass.

The orchestrator will not invoke Codex on this unit (for any finding) until the user has chosen a direction for F1.

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:43-127`

**Issue:** The `.wl` is a line-by-line port of the `.py` (`scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py`): identical canonical-branch setup (L43-46), identical `Series`+`Coefficient` linearization mirroring SymPy's `series().coeff()` (L52-57), and the SymPy-only dummy-symbol `Solve` indirection (`du2sym`, `dP0sym` at L67-74) transliterated into Mathematica where it is unnecessary. The second engine therefore reproduces the first engine's algebra path instead of re-deriving independently. (All implemented identities are symbolically correct, so this is a coverage/assurance defect, not a live math bug.)

**Required change:**
Re-derive the load-bearing quantities in the `.wl` by a path distinct from the `.py`, keeping the final `expectZero` targets unchanged (they are the paper's boxed formulas):
1. Replace the `Series`/`Coefficient` extraction (L52-57) with a derivative-based linearization, e.g.:
   - `du2 = FullSimplify[(D[u2Full, eps] /. eps -> 0), Assumptions -> $Assumptions];`
   - `du4 = FullSimplify[(D[u4Full, eps] /. eps -> 0), Assumptions -> $Assumptions];`
   - `dP0 = FullSimplify[(D[P0Full /. N0 -> P0*D0, eps] /. eps -> 0), Assumptions -> $Assumptions];`
   (`D[..., eps] /. eps -> 0` gives the first-order coefficient since each `*Full` is rational and analytic at `eps = 0`; this is a different mechanism than the SymPy `series().coeff(eps,1)` it currently mirrors.)
2. Remove the `du2sym`/`dP0sym` placeholder idiom (L67-74) and invert directly:
   - `dkappaFromdu2 = FullSimplify[dkappa /. First[Solve[du2 == du2Hyb, dkappa]], Assumptions -> $Assumptions];`
   - `dgammaFromdP0 = FullSimplify[dgamma /. First[Solve[dP0/P0 == dP0OverP0Hyb, dgamma]], Assumptions -> $Assumptions];`
Leave every `expectZero[...]` line and every final formula target exactly as-is. Do not change the `.py`.

**Verification command:**
The verifier runs `redteam exec-mathematica 170`. Confirm: the `.wl` no longer contains `du2sym` or `dP0sym`; the linearization uses `D[..., eps]` rather than `Series`/`Coefficient`; every `expectZero` still prints `= 0` / `PASS`; the script exits 0; and the Mathematica output residuals still match the SymPy ones.

## F3 — cosmetic stage-label correction (low)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:30`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:27`

**Issue:** Both `banner(...)` titles read `"STAGE 153 — LINEAR GROUPED-P2 DIRECT OUTLET MAP"`, but this is unit 170. The saved transcripts inherit the wrong label and could mislead a reviewer into thinking the output belongs to a different stage. No effect on any assertion.

**Required change:**
- py L30: change `banner("STAGE 153 — LINEAR GROUPED-P2 DIRECT OUTLET MAP")` → `banner("STAGE 170 — LINEAR GROUPED-P2 DIRECT OUTLET MAP")`.
- wl L27: change `banner["STAGE 153 — LINEAR GROUPED-P2 DIRECT OUTLET MAP"];` → `banner["STAGE 170 — LINEAR GROUPED-P2 DIRECT OUTLET MAP"];`.

**Verification command:**
After `redteam exec-sympy 170` and `redteam exec-mathematica 170`, both transcripts' banner line reads `STAGE 170 — ...` and both scripts still exit 0.

---

## User resolution (F1 — paper_misalignment)

Direction **(a)** chosen by the user on 2026-05-28: add the missing script check (no paper-side edit). Stage 170 is the owner of the weak-axisymmetric Sec. 5 deliverable (the outlet maps δκ_W/δγ_W inheriting the (1,1/2,-1) signature with amplitudes κ1/γ1); stages 171/173 verify the same lane *pattern* only on other quantities (obstructions / load), so a downstream carry-forward (Cluster C) was not a fit. Resolution: orchestrator-direct addition of a Sec. 5 signature block to BOTH engines.

## Applied: F1
files_changed:
- scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py
- mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl
summary: Added a "Weak-axisymmetric signature (1,1/2,-1)" section to both engines — feeds lane-scaled bundle defects δD_(A,n)=eps·λ_A·D_n^(1) (λ=(1,1/2,-1)) through the SAME verified linear outlet maps and asserts δκ_W^(2A)=eps·λ_A·κ1, δγ_W^(2A)=eps·λ_A·γ1 with κ1=3(1-σ)(D2_1+D0_1/9)/(σ D0), γ1=-(1-σ)(N0_1-P0 D0_1)/(9σ N0), plus the (1,1/2,-1) ratio checks (21=½·20, 22=-20). +10 checks per engine.
deviation: none (orchestrator-applied per user direction a).

## Applied: F2
files_changed:
- mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl
summary: Replaced the Series/Coefficient linearization with D[...,eps]/.eps->0, and removed the SymPy-style du2sym/dP0sym placeholder idiom (direct Solve[du2==du2Hyb,dkappa] / Solve[dP0/P0==dP0OverP0Hyb,dgamma]); final expectZero targets unchanged.
deviation: none (orchestrator-applied).

## Applied: F3
files_changed:
- scripts/moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py
- mathematica/moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl
summary: Banner "STAGE 153" -> "STAGE 170" in both engines.
deviation: none (orchestrator-applied).
