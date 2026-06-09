---
unit_id: 225
batch: VII.1
created_at: 2026-06-02T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T16:50:16-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 225

Apply each non-paper_misalignment finding below in order (F1, F2). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F3) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.wl`

(Exact filename — do NOT drop `_mathematica_audit`. `.wl` files live in `mathematica/`.)

**Issue:** Stage 225 (`is_checkpoint: False`, `is_status_only_candidate: False`) has only a SymPy audit. Mathematica can independently verify every claim here (first-order series expansion of dressed rational expressions, exact `2 Sqrt[2]/Pi` arithmetic at the compatibility sample, matrix rank/nullspace, 2x2 determinant), so the dual-engine rule requires a second engine. Write a NEW, INDEPENDENT Mathematica audit — do not transliterate the SymPy choreography.

**Required change:**
Create the `.wl` at the target path. To satisfy the anti-transliteration guard, derive with native Mathematica primitives and a different decomposition than the SymPy script:
- For first-order slopes, use `Series[expr, {eps, 0, 1}]` then `Coefficient[Normal[...], eps]` (NOT `D[...] /. eps->0` mirroring SymPy's `diff(...).subs(eps,0)`).
- For the mechanism sieve, use `NullSpace`, `MatrixRank`, and `Det` on the assembled compensation matrix; obtain the null basis via `NullSpace` directly rather than the SymPy `-A11^{-1} A12` block-elimination route.
- Use exact arithmetic: `kappa = 2 Sqrt[2]/Pi`, rationals for the sample, and `N[..., 16]` only for the final numeric comparisons.
- Guard floating comparisons with a helper like `expectClose[a_, b_, tol_:10^-12] := If[Abs[a-b] > tol, Print["FAIL: ", a, " != ", b]; Exit[1]]` and guard exact identities with `expectZero[e_] := If[FullSimplify[e] =!= 0, Print["FAIL nonzero: ", e]; Exit[1]]`.

**Claim manifest** (the new script must independently verify each):
- **M1** — Arbitrary-base first-order slopes. With `D_A = D0 + eps*lam*D01` (similarly D2A,D4A,N0A), `u2A = -D2A/D0A`, `u4A = (D2A^2 - D0A*D4A)/D0A^2`, `P0A = N0A/D0A`:
  `u2^(1) = (-D0*D21 + D2*D01)/D0^2`,
  `u4^(1) = (-D0*(D0*D41 + D01*D4 - 2*D2*D21) + 2*D01*(D0*D4 - D2^2))/D0^3`,
  `Xi1 = P1/P0 = N01/N0 - D01/D0`.
- **M2** — Conservative compensation surface. Solving `u2^(1)=0` gives `D21 = -u2*D01` with `u2=-D2/D0`; solving `u4^(1)=0` (with that D21) gives `D41 = (D4/D0)*D01 = (u2^2-u4)*D01`. AND the one-pole reduction (this is the substantive content F2 fixes in SymPy too): under the one-pole constraint `D4 = -3*D0*u2^2` (equivalently `u4 = 4*u2^2`), `D41 = -3*u2^2*D01`. Verify the `-3` is load-bearing (a perturbation to `-2` must FAIL).
- **M3** — Primitive log-slope compiler. With base bundle `C=kappa*lamB, GU=lamU, GW=kappa*lamW, R=kappa*lamR`, `Delta=OmU^2*OmW^2-R^2`, `S2=OmU^2+OmW^2`, `H=GU^2+GW^2`, `Q=GU^2*OmW^2+2*GU*GW*R+GW^2*OmU^2`, `P=OmU^2*GW+R*GU`, and moments `B0=C^2/varpi^2`, `B2=C^2/varpi^4`, `B4=C^2/varpi^6`, `Z0=Q/Delta`, `Z2=(Q*S2-H*Delta)/Delta^2`, `Z4=(Q*(S2^2-Delta)-S2*H*Delta)/Delta^3`, `N0=P^2/Delta^2`, dressing each positive primitive `p -> p*Exp[eps*x_p]`, the first-order drifts equal the boxed forms of notes §4.1-4.3:
  `B0_1 = B0*(2 xLB - 2 xV)`, `B2_1 = B2*(2 xLB - 4 xV)`, `B4_1 = B4*(2 xLB - 6 xV)`,
  `Delta_1 = 2 OmU^2 OmW^2 (xOU+xOW) - 2 R^2 xLR`, `S2_1 = 2 OmU^2 xOU + 2 OmW^2 xOW`, `H_1 = 2 GU^2 xLU + 2 GW^2 xLW`,
  `Q_1 = 2 GU^2 OmW^2 (xLU+xOW) + 2 GU GW R (xLU+xLW+xLR) + 2 GW^2 OmU^2 (xLW+xOU)`,
  `P1raw = OmU^2 GW (2 xOU + xLW) + R GU (xLR + xLU)`,
  `Z0_1 = (Q1*Delta - Q*Delta1)/Delta^2`, `Z2_1`, `Z4_1`, `N0_1 = 2 P P1raw/Delta^2 - 2 P^2 Delta1/Delta^3` (use the §4.2 closed forms), and the bundle compiler `D01 = K xK - B0_1 - Z0_1`, `D21 = -(M xM + B2_1 + Z2_1)`, `D41 = -(B4_1 + Z4_1)`, `N01 = N0_1`.
- **M4** — Compatibility sample. With `(lamB,lamU,lamW,lamR,OmU,OmW,varpi,M) = (1/2, 3/10, 2/5, 1/4, 1, 7/5, 2, 1)` and `D0_compat = 3*(M+B2+Z2)^2/(B4+Z4)`, `K_compat = B0+Z0+D0_compat`: `D0 ~ 24.2373099886223`, `D2 ~ -1.18562046858190`, `D4 ~ -0.173991572849491`, `u2 ~ 0.0489171640391802`, `u4 ~ 0.00957155575054425`, `D4/D0 ~ -0.00717866681290820`, `P0_target ~ 0.00206979231806289`, and the one-pole identity `u4 - 4*u2^2 = 0` (exact). These are carried from upstream Stage 223 — verify against these exact targets, do not re-pick new values.
- **M5** — Concrete `\Xi_1` coefficients at the sample (notes §4.4): `xK: -1.00975540977030`, `xM: 0` (vanishes identically on this branch), `xLB: 0.00418038073077834`, `xV: -0.00418038073077834`, `xLU: 0.324464020216766`, `xLW: 1.69086641859305`, `xLR: 0.423379354382463`, `xOU: -0.747843374599229`, `xOW: -4.11424577297551`.
- **M6** — Wall-only no-go: with only `(xK,xM)` active, the compensation system `D21 + u2*D01 = 0`, `D41 - (D4/D0)*D01 = 0` has only the trivial solution `xK=xM=0` on the sample branch (where `D4/D0 != 0`). Pure-BdG no-go: the 2x2 matrix `[[-(B2+u2 B0), 2 B2 + u2 B0],[-(B4 - (D4/D0) B0), 3 B4 - (D4/D0) B0]]` has determinant `Delta_BdG = -B0 B2 (D4/D0) - 2 B0 B4 u2 - B2 B4`, numerically `~ -5.11886996120011e-5 != 0` at the sample.
- **M7** — Mixed/U survivor: with only `(xLU,xLW,xLR,xOU,xOW)` active, the 2x5 compensation matrix (rows from `D21+u2 D01` and `D41-(D4/D0) D01`) has rank 2 / nullity 3; the three null `\Xi_1` values are `~ 1.36026097049402`, `-14.4310278139755`, `-5.01037421295998`. (You may obtain the basis via `NullSpace` directly; the `\Xi_1` values are basis-normalization dependent for v2/v3 but v1 with last-three components `(1,0,0)` gives `1.36026097049402` — match the first-component-normalized basis of notes §5.3, or document the normalization if `NullSpace` returns a different but equivalent basis and verify the *span* and the sigma1 value.)
- **M8** — Transported amplitude windows (carried from upstream Stage 224): with `sigma1 = 1.36026097049402` and budgets `(0.367930328492646, 0.737619063660757, 2.94889585703134, 4.63505472371892)`, the `\Xi_1`-to-amplitude windows `budget/sigma1` are `0.270485102839510`, `0.542262903708006`, `2.16788978070904`, `3.40747461278373`.

**Verification command:**
The verifier will run `redteam exec-mathematica 225` and confirm the new `.wl` appears, every `expectZero`/`expectClose` passes, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M8 with Series-based slopes, exact sample arithmetic, determinant/rank checks, direct NullSpace normalization, and transported window checks.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py:68-69`

**Issue:** The current block defines `one_pole_D41 = sp.simplify((-3 * u2**2) * D01)` and then asserts `one_pole_D41 + 3*u2**2*D01 == 0`, i.e. `0==0` by construction. It does not verify the claimed one-pole reduction `D_{41}=(u_2^2-u_4)D_{01} -> -3 u_2^2 D_{01}` under `u_4 = 4 u_2^2`. The `-3` coefficient is not load-bearing in the current check.

**Required change:**
Replace lines 68-69 with a check that exercises the one-pole constraint. The symbols `u2` and `u4` at lines 32-33 are defined as expressions in `D0,D2,D4` (`u2 = -D2/D0`, `u4 = (D2**2 - D0*D4)/D0**2`), so substituting `u4 -> 4*u2**2` into `u4` is a no-op. Instead impose the one-pole constraint on `D4` directly (one-pole means `D4/D0 = u2**2 - u4 = u2**2 - 4*u2**2 = -3*u2**2`, i.e. `D4 = -3*D0*u2**2`). Use the already-verified general form `D41_comp = D4*D01/D0` (line 62/65) and substitute the one-pole `D4`:

Before (lines 68-69):
```python
    one_pole_D41 = sp.simplify((-3 * u2**2) * D01)
    assert sp.simplify(one_pole_D41 + 3 * u2**2 * D01) == 0
```
After:
```python
    # One-pole reduction: on a one-pole branch D4/D0 = u2^2 - u4 with u4 = 4 u2^2,
    # so D4 = -3 D0 u2^2. Substitute into the general surface D41 = (D4/D0) D01.
    D4_one_pole = -3 * D0 * u2**2
    one_pole_D41 = sp.simplify((D41_comp).subs(D4, D4_one_pole))
    assert sp.simplify(one_pole_D41 - (-3 * u2**2) * D01) == 0
```
Note: `u2 = -D2/D0` does not contain `D4`, so `D4_one_pole = -3*D0*u2**2 = -3*D2**2/D0` is independent of `D4`; substituting `D4 -> D4_one_pole` into `D41_comp = D01*D4/D0` yields `D01*(-3*D2**2/D0)/D0 = -3*D2**2*D01/D0**2`, and `(-3*u2**2)*D01 = -3*(D2**2/D0**2)*D01` — identical. The `-3` is now load-bearing: changing it to `-2` on either side makes the assert fail. Keep the subsequent `print` at line 74 (`one_pole_D41` still holds the reduced form).

**Self-test (performed):** `D41_comp` (line 62) simplifies to `D01*D4/D0`. Substituting `D4 -> -3*D0*u2**2` gives `-3*D0*u2**2*D01/D0 = -3*u2**2*D01` — matches the RHS, so the assert passes for `-3` and fails for any other coefficient. Variable-independence trap cleared: `u2` is independent of `D4`, so the substitution is well-defined and non-trivial.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 225`, confirms the new line-68 block references the one-pole `D4`/`u4=4 u2^2` constraint (not a both-sides restatement), and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py`
- summary: Replaced the tautological one-pole check with a substitution of `D4 = -3 D0 u2^2` into the solved general compensation surface.
- deviation: none

## F3 — paper_misalignment (subtype: notes_contradicts_script)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md:5` quote: "built on the Stage-240 compatibility branch and the transported Stage-241 same-charge ceiling test."
- same notes `:115` quote: "For the concrete Stage-240 compatibility sample"; `:152` "Stage-241 transported same-charge ceilings"; `:50` "after Stage 242"; `:590` cites supporting file `moving_throat_pde_stage242_..._sympy_audit.py` (which does not exist).
- `/var/projects/toy_physics/.../paper/stages/stage_225.tex` — the card's Inputs/Verification-note text draws from these notes (numbering inherited).

**Script side (correct):**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_..._sympy_audit.py:195` quote: "# Concrete Stage 223 compatibility point"
- same script `:240` "Verified concrete Stage 223 compatibility point:"; `:392` "# Transported Stage 224 headroom on the first surviving family"

**Confirmed value provenance (red-team grep):** the compatibility-point literals (`K_compat=24.4737548792910`, `D0=24.2373099886223`, `P0=0.00206979231806289`) are owned by Stage 223 (`scripts/moving_throat_pde_stage223_..._sympy_audit.py:252,254,372`); the four `\Xi_1` ceiling budgets are owned by Stage 224 (`scripts/moving_throat_pde_stage224_..._sympy_audit.py:152-158`). The script's "Stage 223 / Stage 224" comments are correct; the notes' "Stage-240/241/242" and self-reference "Stage 242" are stale (pre-renumbering).

## Resolve before fix_loop

The notes file (and, by inheritance, the paper card's Inputs/Verification-note prose) attributes this unit's carried inputs to "Stage 240 / Stage 241" and self-refers to the unit as "Stage 242", while the script comments and the verified value provenance say the owners are Stage 223 (compatibility point) and Stage 224 (ceiling budgets) and the unit is Stage 225. The values themselves match the script's correct upstream owners; only the stage *numbers* in the prose are stale.

Possible directions (the user picks one):
- (a) Prose is stale → in a follow-up directive, authorize updating the notes (lines 5, 50, 115, 152, 154, and the §8 supporting-file citation at line 590, plus the other "Stage-240/241/242" mentions) and any paper-card text from `240/241/242` to `223/224/225`. No script change.
- (b) Some intentional alternate numbering is in play (e.g. a parallel manuscript using the 24x scheme) → flag for deeper review of the whole 223-226 prose block.
- (c) The card itself is fine and only the notes need the renumber → narrow the prose edit to the notes file only.

F3 is RESOLVED by the user (see the ## RESOLVED block at the end): apply the authorized notes renumber. F1 and F2 are independent script-side fixes and proceed normally.

## RESOLVED — F3 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: renumber notes prose to canonical (values already match; only stage numbers + a filename token are stale). Notes-only; Codex applies, Claude reviews. Do NOT touch scripts, paper.tex, or appendix.
- In `notes/stages/moving_throat_pde_stage225_..._sympy_audit.md`, align all stale stage-number references to the canonical numbering used in the .py script comments + card: cited inputs "Stage 240"/"Stage 241" → "Stage 223"/"Stage 224", and the self-reference "Stage 242" → "Stage 225". Also fix the nonexistent supporting-file citation containing `stage242` to the canonical `stage225` filename token.
- Acceptance: no stale `Stage 240`/`Stage 241`/`Stage 242` or `stage242` reference remains. Append `## Applied: F3`.

## Applied: F3

- files_changed:
  - `notes/stages/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md`
- summary: Renumbered the notes prose from stale Stage 240/241/242 references to canonical Stage 223/224/225 references and fixed the supporting filename token.
- deviation: none
