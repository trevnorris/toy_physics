---
unit_id: 225
batch: VII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-02T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 225 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_225.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at lines 62, 684-745 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py`
- mathematica: (missing)
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.txt`
- mathematica output: (missing)

## What the paper claims

Stage 225 derives the first microscopic formula for the weak-axisymmetric prefactor slope `\Xi_1`. The card's `Output` (stage_225.tex:15) is: "Microscopic compiler $\Xi_1=N_{01}/N_0-D_{01}/D_0$ together with first-order even-preserving compensation conditions and the first surviving mixed-sector family." The notes (§2-§6) enumerate distinct deliverables: (1) exact arbitrary-base first-order formulas for `u_2^{(1)}`, `u_4^{(1)}`, `\Xi_1`; (2) the exact conservative compensation surface `D_{21}=-u_2 D_{01}`, `D_{41}=(D_4/D_0)D_{01}` and its one-pole reduction `D_{41}=-3u_2^2 D_{01}` (using `u_4=4u_2^2`); (3) the full primitive logarithmic-slope compiler from nine `x_p` slopes to `D_{01},D_{21},D_{41},N_{01}`; (4) the concrete compatibility-point values `D_0,D_2,D_4,u_2,u_4,P_{0,target}`; (5) the concrete `\Xi_1` linear-coefficient form; (6) the wall-only and pure-BdG generic/sample no-go; (7) the mixed/U rank-2/nullity-3 nullspace, its null basis, and the corresponding `\Xi_1` values; (8) the transported amplitude windows on the first surviving mixed family. The appendix (lines 684-745) confirms `\Xi_{load}=N_{01}/N_0-D_{01}/D_0=\Xi_1` and that wall-only/pure-BdG do not supply nontrivial same-charge load while a mixed-sector family survives.

## What the script claims to verify

The SymPy script verifies, in order: (a) the arbitrary-base `u_2^{(1)},u_4^{(1)},\Xi_1` formulas by differentiating dressed expressions and asserting against independently-typed closed forms (lines 36-51); (b) the compensation surface by solving `u_2^{(1)}=0`/`u_4^{(1)}=0` and checking the solutions equal `-u_2 D_{01}` and `(D_4/D_0)D_{01}=(u_2^2-u_4)D_{01}` (lines 61-66); (c) the primitive moment-drift compiler `B_{*,1},Z_{*,1},N_{0,1},D_{*1}` against closed forms (lines 140-190); (d) the concrete compatibility point `D_0,D_2,D_4,u_2,u_4,P_0` against the carried numeric literals (lines 230-238); (e) the nine `\Xi_1` linear coefficients (lines 270-278); (f) the wall-only trivial-solution and pure-BdG nonzero-determinant no-go (lines 287-321); (g) the mixed/U rank-2 matrix, null basis, and `\Xi_1` values (lines 326-389); (h) the transported `\epsilon t` amplitude budgets (lines 401-406). This is a substantive, mostly non-tautological battery.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `\Xi_1=N_{01}/N_0-D_{01}/D_0` (Output) | line 51 assert | match |
| `u_2^{(1)},u_4^{(1)}` arbitrary-base | lines 41-50 | match |
| compensation surface `D_{21}=-u_2 D_{01}`, `D_{41}=(D_4/D_0)D_{01}` | lines 64-65 | match |
| one-pole reduction `D_{41}=-3u_2^2 D_{01}` (via `u_4=4u_2^2`) | lines 68-69 | mismatch (tautological — see F2) |
| primitive log-slope compiler | lines 158-190 | match |
| compatibility-point `D_0,D_2,D_4,u_2,u_4,P_0` | lines 230-238 | match |
| concrete `\Xi_1` coefficients | lines 270-278 | match |
| wall-only / pure-BdG no-go | lines 287-321 | match |
| mixed/U nullspace + `\Xi_1` values | lines 326-389 | match |
| transported `\epsilon t` windows | lines 401-406 | match |
| (whole stage) second-engine Mathematica check | none | missing (F1) |
| upstream-stage attribution (notes/card say 240/241/242) | script comments correctly say 223/224 | mismatch (F3, prose-side) |

Overall `paper_alignment: partial` — the mathematical content the script verifies faithfully matches the paper's `Output` and notes deliverables; the partial flag is driven by (i) the stale cross-document stage numbering (F3, prose-side, needs user resolution) and (ii) the tautological one-pole check (F2).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 41-42 | `simplify(u2_1 - …)==0` | u_2^{(1)} | yes |
| A2 | sympy | 43-50 | `simplify(u4_1 - …)==0` | u_4^{(1)} | yes |
| A3 | sympy | 51 | `simplify(Xi1 - (N01/N0 - D01/D0))==0` | Output `\Xi_1` | yes |
| A4 | sympy | 64-65 | `simplify(D21_comp + u2*D01)==0`, `D41_comp - D4*D01/D0` | compensation surface | yes |
| A5 | sympy | 66 | `simplify(D41_comp - (u2**2-u4)*D01)==0` | compensation surface (alt form) | yes |
| A6 | sympy | 68-69 | `one_pole_D41 = (-3*u2**2)*D01`; assert `one_pole_D41 + 3*u2**2*D01 == 0` | one-pole reduction | **no (tautological)** |
| A7 | sympy | 158-182 | `simplify(B*_1/Z*_1/N0_1 - …)==0` | primitive moment compiler | yes |
| A8 | sympy | 188-190 | `simplify(D*_compiler ∓ …)==0` | first-order bundle compiler | yes |
| A9 | sympy | 230-238 | `assert_close(D0/D2/D4/u2/u4/P0, …)` + `u4_s-4*u2_s**2==0` | compatibility point + one-pole | yes |
| A10 | sympy | 270-278 | `assert_close(Xi1_coeffs[…], …)` | concrete `\Xi_1` form | yes |
| A11 | sympy | 295-297 | wall solutions == {xK:0,xM:0} | wall-only no-go | yes |
| A12 | sympy | 312-317 | det form + sample nonzero | pure-BdG no-go | yes |
| A13 | sympy | 341-378 | rank==2, basis, `\Xi_1` values | mixed/U survivor | yes |
| A14 | sympy | 403-406 | `assert_close(t_budgets[…], …)` | transported windows | yes |
| — | mathematica | — | (no script) | all of the above | missing |

## Findings

### F1 — missing_verification_script (subtype: missing_mathematica)

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/` (no `*stage225*.wl` exists)

**What's wrong:**
Stage 225 is `is_checkpoint: False` but `is_status_only_candidate: False` — it does original derivation work (compensation surface, primitive compiler, mechanism sieve). Per the dual-engine rule, a `.wl` is REQUIRED wherever Mathematica CAN independently verify, and the test is "is it possible," not "is it necessary." Every claim here is squarely inside Mathematica's native capability: first-order series/derivative expansion of dressed rational expressions, exact rational + `\sqrt2/\pi` arithmetic at the compatibility sample, matrix rank/nullspace and a 2x2 determinant. There is no genuine obstruction, so single-engine is not justified.

**Why this matters:**
The SymPy script is the sole check on the stage's load-bearing identities (`\Xi_1`, the compensation surface, the mechanism sieve). A second engine deriving the same results from the physical premises via different primitives is the only guard against a systematic SymPy transcription/simplification error.

**Required change:**
Add an independent Mathematica audit script at the exact target path (see directive F1). It must re-derive — not transliterate — the claim manifest M1..M8.

**Verification:**
`redteam exec-mathematica 225` produces the new `.wl`'s output, all in-file checks pass, and the script exits 0.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_microscopic_xi1_compiler_first_order_conservative_compensation_surface_and_mixed_sector_survival_sieve_sympy_audit.py:68-69`

**What's wrong:**
```python
one_pole_D41 = sp.simplify((-3 * u2**2) * D01)
assert sp.simplify(one_pole_D41 + 3 * u2**2 * D01) == 0
```
`one_pole_D41` is *defined* as `(-3*u2**2)*D01`; the assertion then checks `(-3u2²D01) + 3u2²D01 == 0`, which is `0==0` by construction. It cannot fail no matter what the physics is. The paper/notes claim it exercises (notes §3, lines 254-261) is the nontrivial *reduction* of the compensation surface `D_{41}=(u_2^2-u_4)D_{01}` to `D_{41}=-3u_2^2 D_{01}` using the one-pole identity `u_4=4u_2^2`. That reduction step — substituting `u_4=4u_2^2` into `(u_2^2-u_4)` to get `-3u_2^2` — is never symbolically performed; the assertion just restates its own definition. (The one-pole identity `u_4=4u_2^2` itself is checked only numerically at the sample point, line 238, and is never connected to the `D41` reduction.)

**Why this matters:**
A reader scanning the output line `On a one-pole base branch: D41 = -3*D01*D2**2/D0**2` (output line 19) would believe the script proved the one-pole reduction. It did not; it only printed back its own ansatz. The substantive content of paper notes §3 (the `(u_2^2-u_4) → -3u_2^2` step) goes unverified.

**Required change:**
Replace the tautological self-check with a real reduction. On a one-pole branch the symbolic relation is `u_4 = 4 u_2^2`. Verify that the general compensation-surface form `(u2**2 - u4)*D01` reduces to `-3*u2**2*D01` under that substitution, e.g.:
```python
one_pole_D41 = sp.simplify(((u2**2 - u4)*D01).subs(u4, 4*u2**2))
assert sp.simplify(one_pole_D41 - (-3 * u2**2) * D01) == 0
```
where `u2, u4` are the symbolic arbitrary-base expressions already defined at lines 32-33 (`u2 = -D2/D0`, `u4 = (D2**2 - D0*D4)/D0**2`). Note: substituting `u4 -> 4*u2**2` into the *expression* `u4` is a no-op since `u4` is defined in terms of `D0,D2,D4`; the substitution must be applied to the symbolic `u4` *symbol*, OR the reduction must be framed as: impose `D4 = D0*u2**2 - D0*(4*u2**2) = -3*D0*u2**2` (the one-pole constraint `D_4/D_0 = u_2^2 - u_4 = -3u_2^2`) and confirm `D41_comp = (D4/D0)*D01` then equals `-3*u2**2*D01`. See directive F2 for the self-tested concrete form.

**Verification:**
The new assertion at line ~68 should fail if the coefficient `-3` is perturbed (e.g., to `-2`); the verifier confirms the assertion references the one-pole constraint `u_4=4u_2^2` (or equivalently `D_4=-3D_0u_2^2`) rather than restating `-3u2**2*D01` on both sides.

### F3 — paper_misalignment (subtype: notes_contradicts_script)

**Severity:** medium
**Files:**
- paper side: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage225_…_sympy_audit.md:5,50,542,590` and `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_225.tex` (Inputs/Verification note references)
- script side: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage225_…_sympy_audit.py:195,240,392`

**What's wrong:**
The notes file and parts of the paper card attribute the carried inputs to the wrong upstream stages, while the script's comments name the *correct* ones. The carried numeric literals (verified below) are owned by Stage 223 (compatibility point `K_compat=24.4737548792910`, `D_0=24.2373099886223`, `P_0=0.00206979231806289` — confirmed present in `scripts/moving_throat_pde_stage223_…_sympy_audit.py:252,254,372`) and Stage 224 (the four `\Xi_1` ceiling budgets `0.367930328492646`, `0.737619063660757`, `2.94889585703134`, `4.63505472371892` — confirmed in `scripts/moving_throat_pde_stage224_…_sympy_audit.py:152-158`).

- Notes §"Status"/§1.2 (lines 5, 115, 152, 154): "built on the Stage-240 compatibility branch and the transported Stage-241 same-charge ceiling test"; "concrete Stage-240 compatibility sample"; "Stage-241 transported same-charge ceilings."
- Notes lines 50, 542, 590 self-refer to this unit as "Stage 242" (e.g. line 50 "after Stage 242", line 590 cites a nonexistent supporting file `moving_throat_pde_stage242_…_sympy_audit.py`). The file is named `stage225`.
- Script comments (lines 195, 240) correctly say "Concrete Stage 223 compatibility point"; line 392 correctly says "Transported Stage 224 headroom."
- Paper card stage_225.tex `Inputs`/`Downstream` are internally consistent with stage 225 but the notes it draws from carry the 240/241/242 numbering.

This is the classic off-by-(stale-numbering) drift: the notes and card text were authored under an earlier stage-numbering scheme (240/241/242) that was later renumbered to 223/224/225; the script was updated, the prose was not. It is not the `R_Q = script+68` sibling pattern (no `R_Q` quantity appears in this unit), but it is the same family of stale-cross-reference defect the warning anticipated.

**Why this matters:**
A reader following the notes to verify the carried `K_compat`/`P_0`/ceiling literals would look in Stage 240/241 (which own different physics) and either find nothing or the wrong values, breaking the provenance chain. The notes also cite a supporting script `stage242…` that does not exist.

**Required change (USER RESOLUTION — not Codex):**
See directive `## Resolve before fix_loop`. The script is correct and must not be edited. Direction (update notes+card numbering 240/241/242 → 223/224/225, vs. some intended renumbering) is the user's call; Codex does not edit notes/ or paper/.

**Verification:**
After user picks a direction and a follow-up directive authorizes the prose edit, the notes' "Stage-240/241/242" references resolve to "Stage 223/224/225" and the supporting-file citation at notes:590 names the existing `stage225…` script.

## Independent-derivation check (Mathematica)

No `.wl` exists, so transliteration is moot. The directive's claim manifest specifies an independent re-derivation using Mathematica-native primitives (`Series`/`Normal`/`Coefficient`, `NullSpace`, `MatrixRank`, `Det`, exact `2 Sqrt[2]/Pi`) with a different decomposition than the SymPy `diff(...).subs(eps,0)` choreography, to satisfy the anti-transliteration guard.

## Engine cross-check

Only one engine present; `engines_agree: n/a`.

## Verdict justification

The SymPy script is substantively correct and faithfully exercises the paper's `Output` and the notes' enumerated deliverables: I attacked the arbitrary-base formulas (derived by genuine `eps`-differentiation, not hardcoded), the compensation-surface solve (genuine `solve` then independent re-check), the primitive compiler (closed forms typed independently from the dressed-expression derivatives), the sample-point literals (back-derived from the one-pole relation, anchored to upstream Stage 223), and the rank/nullspace/determinant sieve (all computed, none hardcoded) — they hold up. Three findings remain: (F1) no second engine despite Mathematica being fully capable, (F2) a single genuinely tautological self-check at lines 68-69 that fails to verify the claimed one-pole reduction, and (F3) a prose-side stale-numbering misalignment (notes/card say 240/241/242; the actual owners and the script comments say 223/224/225) that requires user resolution. Verdict `findings`; no stop-cold — the math is sound and reconcilable, and F3 is a prose provenance fix, not a math break.

## Self-test notes

I checked the variable-independence trap on F2's prescribed fix: the one-pole reduction must substitute the constraint into a symbol/expression that actually carries `u_4` (or `D_4`); substituting `u4 -> 4*u2**2` into the already-expanded `u4 = (D2**2-D0*D4)/D0**2` would be a no-op, so the directive frames the fix via the `D_4 = -3 D_0 u_2^2` one-pole constraint to guarantee the `-3` coefficient is load-bearing and a perturbation to `-2` would fail. I checked symmetry/parity — no unbounded integrals here. I confirmed the carried literals against their true upstream owners (223 for `K_compat`/`D_0`/`P_0`, 224 for the four ceilings) by direct grep of those scripts, so F3 is a numbering-attribution defect, not a value defect; the values themselves match. No `R_Q` quantity exists in this unit, so the sibling `script+68` typo pattern does not apply.
