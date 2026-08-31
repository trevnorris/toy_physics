# S11c-b #88 — BLAST RADIUS of the §3a frozen-spurion correction (RESULT)

## Headline (settled, quadruple-confirmed)
Correcting the frozen-spurion §3a defect (engine basis 26 → correct 40: complete the spurion
basis to 15/source **and** retain the background curvature — the Hessian ∂²W_bg, ∂²μ_R,bg — in the
Euler–Lagrange total-divergence step) **DISTURBS EVERY strong stored-energy EL row at retained
grade (η,σ_W ≤ 1)**, with genuinely **non-absorbable** new operator structure:

| strong row (SLAB_OPERATOR) | RESIDUAL termcount | RANK_GAIN (const-coeff span) |
|---|---:|---:|
| U_MOMENTUM_1 | 98 | 8 |
| U_MOMENTUM_2 | 98 | 8 |
| U_MOMENTUM_3 | 98 | 8 |
| MU_THETA (θ) | 108 | 6 |
| THICKNESS_EW (e_W) | 186 | 4 |

`RANK_GAIN > 0` on every row ⇒ the correction adds operator directions that **no
constant-coefficient reparametrisation of the frozen row can absorb**. Zeroing the Hessian map
drives `RANK_GAIN → 0` on every row (both build legs) ⇒ the non-absorbable content **is** the
background-curvature (Hessian) coupling; the omitted-invariant first-jet parts are absorbable
re-parametrisations. `NEW_SECOND_JET_ATOMS` (10–12 `w1_profile_dij`/`m1_profile_dij` per row,
absent from every frozen row) is the cheap sufficient witness; `GRADE_SUPPORT = {(0,1),(1,1)}`
(all at σ_W¹, retained).

## Consequence for #89 (the point of #88)
The cross-engine family verdicts were adjudicated on the FROZEN (incomplete) operator. Since the
correct operator differs non-absorbably from the emitted PY operator on these rows:
- **KINETIC family — verdict INVALIDATED, re-adjudicate at #89.** The u-momentum and thickness
  (e_W) inertia rows are disturbed; the settled "WL-correct / representational" kinetic verdict
  rested on the frozen operator.
- **θ-strong row (MU_THETA) — disturbed; reopen any θ-row adjudication at #89.**
- **COUPLING — already GENUINE (#84); also disturbed (the ~118-term bulk residual).**
- ⇒ #89 must re-adjudicate across families, **not only coupling.**

## Scope (this instrument is a PY-side disturbance WITNESS — LAB_HELD stored-energy EL rows only)
⛔ A nonzero row ⇒ the emitted PY row is not §1d-correct there ⇒ its adjudication is invalidated.
A **zero/absent** measurement here does **NOT** clear a family. NOT measured here, each to be
checked at #89:
- **ADVECTIVE / mass row** — built as the constraint `advective_constraint = ε·u_t·∇ρ_br`
  (`sympy_audit.py:1492`), NOT an energy-EL, so the omitted energy invariants do not feed it;
  **likely spared, confirm at #89** (do not assert clearance from this instrument).
- **ADMISSIBILITY ε⁰ operator** (`admissibility_operator`, separate route, already Hessian-aware
  via `background_dx` for the 8 selected) — its omitted-invariant contributions are un-adjudicated.
- **WL frozen rows** — WL's 8 are a different hand-coded subspace (#86); a PY-zero row does not
  imply the WL row is undisturbed.
- **MATERIAL_ADVECTED** rows (frozen `material_pullback` dx, `:1348-1377`).

## How it was established (rule 1/13 — independent constructions that agree)
1. **Orchestrator anchor** (`~/.s11_build/S11c_b_88_probes/`): `probe_blast_radius.py`
   (per-row residual + Control A + driver decomp + new-atom witness), `probe_decomp.py`
   (the Hessian-on-omitted cross-term + monomial non-absorbability). Termcounts 98/98/98/108/186.
2. **Codex-built instrument** `scripts/S11c_b_88_blast_radius.py` (from directive
   `S11c_b_88_blast_radius_build_directive.md`, 2 decision legs folded): exit 0, four structural
   controls assert-pass; termcounts + RANK_GAIN as tabled. Run `~/.s11_build/S11c_b_88_run.out`.
3. **Grok build leg**: independent from-scratch MU_THETA derivation = artifact (108 terms, all 12
   Hessian atoms); Hessian-zeroing → RANK_GAIN 0.
4. **Fresh Claude build leg**: independent from-scratch derivation = artifact (108/6/36/66); every
   control ablation bites; absorbability span correct (`g·q`→gain 0, `H·q`→gain 1).
All four agree. Transcripts `~/.s11_build/S11c_b_88_{decision,build}_{codex,grok}.txt` +
`~/.s11_build/S11c_b_88_probes/ANCHOR.md`.

## OWED — instrument control-hardening (fold at #89's operator rebuild; NOT result-changing)
Grok build leg found two control-robustness gaps (both theoretical here — the result is
quadruple-confirmed and `density_frozen` comes straight from `construct_energy`; the fresh-Claude
leg confirmed the controls that DO exist all bite):
1. `CONTROL_ENGINE` cross-checks two EL recipes on the **same** `density_frozen` — it validates the
   EL recipe but would not catch a **wrong** `density_frozen` (the #86 relabeling class). Harden:
   source the engine row from `build_operator("LAB_HELD")`'s independently reconstructed energy.
2. `CONTROL_JACOBIAN` asserts termcount>0 only; swapped source jets flip the printed identity
   flags but still exit 0. Harden: assert the jet-identity flags (all `CARRIES_SOURCE_JET`).

## Next
#88 SETTLED. NEXT = **#87** (WL undercomplete — code-resolved; quick CAS confirm WL's 8 ⊂ correct 15)
→ **#89** both-engine §3a repair (fixes DIFFER per engine; must re-adjudicate KINETIC + θ, not only
coupling; harden the two controls above) → **#90** PY §3c content fix on the corrected basis.
