---
unit_id: 166
batch: V.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage166_bundle_inversion_four_drifts.md]
  paper_appendix: present
---

# Audit unit 166 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_166.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage166_bundle_inversion_four_drifts.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows: line 63 status table; lines 339-349 bundle-observable narrative)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.txt`

## What the paper claims

`\stagefield{Output}` (stage_166.tex:15): *"Inverts the four observables \((\Theta_w,K_s,K_q,P_0)\) to determine the four remaining microscopic drifts."* The notes give the complete deliverable set. (1) Four exact logarithmic branch laws: `δln Θ_w = 2 δln ρ_w`, `δln K_s = 2 δln a + δln ρ_w`, `δln K_q = δln Z_q + 2 δln c_s − 2 δln a`, `δln P_0 = 5(δln c_s − δln a)`, with `P_0 = (54G/5c⁵) c_s⁵/a⁵`. (2) Exact inversion: `δln ρ_w = ½ δln Θ_w`, `δln a = ½ δln K_s − ¼ δln Θ_w`, `δln c_s = ½ δln K_s − ¼ δln Θ_w + ⅕ δln P_0`, `δln Z_q = δln K_q − ⅖ δln P_0`. (3) Equivalent full-bundle form substituting `δln P_0 = δln N_0 − δln D_0`. (4) Frozen-wall corollary (`δln Θ_w = 0`) giving `δln ρ_w = 0`, `δln a = ½ δln K_s`, etc. (5) Numerics: `Θ_w = 25 λ²ρ²` so `Θ_w^(χ) ≈ 4.06863235008162 λ²` ⟹ `ρ_w^(χ) = √(Θ^(χ)/25) λ⁻¹ ≈ 0.403417022451042 λ⁻¹`.

## What the script claims to verify

The SymPy script (docstring lines 5-12) sets up the four log branch laws as `sp.Eq`, solves the 4×4 linear system for `(drho, da, dcs, dZ)`, then asserts each solved drift equals the notes' boxed closed form via `expect_zero`. It then re-substitutes the solution into the original four laws (forward verification residuals = 0), substitutes `dP → dN0 − dD0` for the bundle form, specializes to the frozen-wall case (`dTheta → 0`), and computes the numeric `rho_w^(chi)`. The Mathematica script does the same and adds an extra "Independent matrix-inverse cross-check" section (lines 59-79): it builds the forward map matrix `Mmat`, inverts it, compares the inverse-derived drifts to the notes' boxed forms, and runs a forward-transcription check confirming `Mmat` encodes the hand-typed boxed laws (`fwdLaws`).

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Four branch laws (notes §1) | eq1–eq4 (py 40-43 / wl 31-34); forward verification (py 59-62 / wl 54-57) | match |
| Inversion forms (notes §2) | `drho/da/dcs/dZ general` (py 53-56 / wl 48-51) | match |
| Bundle form `dP=dN0−dD0` (notes §3) | bundle identities (py 73-80 / wl 86-93) | match |
| Frozen-wall corollary (notes §4) | frozen checks (py 88-91 / wl 104-107) | match |
| `ρ_w^(χ) ≈ 0.403417…` (notes §4) | `rho_chi` (py 94-96 / wl 110-112) | match |
| Matrix-inverse independent route (extra) | wl 59-79 | extra (legitimate cross-check, not in paper) |

The "extra" matrix-inverse block is an additional independent verification route, not a contradiction; no `paper_misalignment`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53-56 | `expect_zero(sol[x] − boxed)` | inversion forms (notes §2) | yes |
| A2 | sympy | 59-62 | `expect_zero(law.subs(sol))` | branch laws round-trip (notes §1) | yes |
| A3 | sympy | 73-80 | `expect_zero(bundle − boxed)` | bundle form (notes §3) | yes |
| A4 | sympy | 88-91 | `expect_zero(frozen − boxed)` | frozen corollary (notes §4) | yes |
| A5 | math | 48-51 | `expectZero[sol − boxed]` | inversion forms | yes |
| A6 | math | 54-57 | `expectZero[law /. sol]` | branch laws round-trip | yes |
| A7 | math | 69-72 | `expectZero[inv.v − boxed]` | inversion via matrix inverse | yes |
| A8 | math | 79 | `expectZero[Total[(Mmat.v − fwdLaws)^2]]` | Mmat ≡ hand-typed boxed laws | yes (non-tautological; fwdLaws typed independently) |
| A9 | math | 86-93 | `expectZero[bundle − boxed]` | bundle form | yes |
| A10 | math | 104-107 | `expectZero[frozen − boxed]` | frozen corollary | yes |

Note on A1/A5: these solve the system independently (`sp.solve` / `Solve`) and compare the *machine-derived* solution against the *hand-typed boxed* RHS from the notes. They are non-tautological because the boxed RHS is not used to build the system — a wrong boxed coefficient would produce a nonzero residual. The output confirms residual = 0 for all.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent re-derivation**, not a transliteration. Three pieces of evidence:

1. **Extra verification route absent from the `.py`.** The Mathematica script adds an entire "Independent matrix-inverse cross-check" section (wl 59-79) — builds `Mmat`, computes `Inverse[Mmat]`, and checks the inverse-derived drift vector against the boxed forms. The SymPy script has no matrix-inverse path at all; it uses only `sp.solve`. A line-by-line port would not invent a second algebraic route the source lacks.

2. **The forward-transcription guard is purpose-built (wl 78-79).** `fwdLaws = {2*drho, drho + 2*da, dZ + 2*dcs - 2*da, 5*(dcs - da)}` is hand-typed from the notes' boxed laws, then compared to `Mmat . {drho,da,dcs,dZ}` via `Total[(...)^2]`. This is a genuine independent check that `Mmat` faithfully encodes the four laws — it has no counterpart in the `.py`. (This is the first-pass batch-7 remediation; it is intact and substantive: the residual would be nonzero for any wrong `Mmat` coefficient, and the squared-sum scalarization is correctly motivated in the comment because `expectZero` tests `res === 0`, which is False for a list.)

3. **Different solve choreography.** SymPy uses `sp.solve((eq1..eq4), (drho,da,dcs,dZ), dict=True)[0]` (py 45); Mathematica uses `Solve[{eq1..eq4}, {drho,da,dcs,dZ}, Reals][[1]]` (wl 35) AND additionally the matrix inverse. The shared *structure* (same four named symbols, same banners) is inevitable because both must verify the same notes — but the algebra is carried out by two different methods, and the `.wl` carries a third method the `.py` does not. Confidence: high that this is independent, not a port.

The boxed branch laws in both files match the notes verbatim: notes §1 `Θ_w = 25 λ²ρ²` ⟹ `δln Θ_w = 2 δln ρ_w` (py 40 `eq1 = sp.Eq(dTheta, 2*drho)` / wl 31 `eq1 = dTheta == 2*drho`); notes §1 `P_0 = (54G/5c⁵) c_s⁵/a⁵` ⟹ `δln P_0 = 5(δln c_s − δln a)` (py 43 `eq4 = sp.Eq(dP, 5*(dcs - da))` / wl 34 `eq4 = dP == 5*(dcs - da)`).

## Engine cross-check

Both engines agree. SymPy output line 5-8 and Mathematica output line 5-8 give identical solved drifts (modulo cosmetic factoring: SymPy `da = dKs/2 - dTheta/4` vs Mathematica `da = (2*dKs - dTheta)/4` — algebraically identical). All `general`, `forward`, `bundle`, and `frozen` residuals are 0 in both. The numeric `rho_w^(chi)`: SymPy `0.403417022451042341`, Mathematica `0.40341702245104232684…` — agree to 18 sig figs and to the notes value `0.403417022451042`. The Mathematica-only matrix-inverse block also yields all-zero residuals (output lines 37-46).

## Verdict justification

`clean`. I attacked: (a) the solve-vs-boxed assertions for circularity — they are not tautological, the boxed RHS is independent of the solved LHS; (b) the matrix-inverse and forward-transcription guards — they are non-trivial and would fail on a wrong coefficient; (c) every coefficient against the notes' boxed laws — all match (25, 54G/5c⁵, ½, ¼, ⅕, ⅖); (d) the `ρ_w^(χ)` numeric — matches notes to full precision; (e) the warned stale constants (`168π²`/`100π²`, Family-1 radius `√(4107−100π²)/(10π)`) — none appear; this stage uses clean rationals only. I confirmed the paper card Output, the notes §1–4 deliverables, and the appendix bundle-observable row all match the script's checks. The first-pass batch-7 vacuous-round-trip remediation is intact and substantive. Both engines independent, both fresh, both pass.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 11 values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δln ρ_w = ½ δln Θ_w` (drho) | py 53 / wl 48; out py5/wl5 | notes.md:102-107 (boxed) | MATCH |
| `δln a = ½ δln K_s − ¼ δln Θ_w` (da) | py 54 / wl 49; out py6/wl6 | notes.md:108-116 (boxed) | MATCH |
| `δln c_s = ½ δln K_s − ¼ δln Θ_w + ⅕ δln P_0` (dcs) | py 55 / wl 50; out py7/wl7 | notes.md:117-127 (boxed) | MATCH |
| `δln Z_q = δln K_q − ⅖ δln P_0` (dZ) | py 56 / wl 51; out py8/wl8 | notes.md:128-136 (boxed) | MATCH |
| branch law `δln Θ_w = 2 δln ρ_w` (coef 2; Θ_w=25λ²ρ²) | py 40 / wl 31 | notes.md:49-57 | MATCH |
| branch law `δln K_s = 2 δln a + δln ρ_w` | py 41 / wl 32 | notes.md:59-67 | MATCH |
| branch law `δln K_q = δln Z_q + 2 δln c_s − 2 δln a` | py 42 / wl 33 | notes.md:70-81 | MATCH |
| branch law `δln P_0 = 5(δln c_s − δln a)` (P_0=54G/5c⁵·c_s⁵/a⁵) | py 43 / wl 34 | notes.md:83-92 | MATCH |
| bundle `δln c_s = ½ δln K_s − ¼ δln Θ_w + ⅕(δln N_0 − δln D_0)` | py 75 / wl 88; out py29/wl51 | notes.md:153-162 (boxed) | MATCH |
| bundle `δln Z_q = δln K_q − ⅖(δln N_0 − δln D_0)` | py 79 / wl 92; out py30/wl52 | notes.md:163-171 (boxed) | MATCH |
| `ρ_w^(χ) ≈ 0.403417022451042 λ⁻¹` (from Θ_w^(χ)=4.06863235008162) | py 94-96 / wl 110-112; out py49/wl77 | notes.md:203-210 (`0.403417022451042`); Θ_w^(χ) at notes.md:180-181 | MATCH |

Frozen-wall corollary forms (`drho|frozen=0`, `da|frozen=½ δln K_s`, `dcs|frozen=½ δln K_s + ⅕ δln P_0`, `dZ|frozen` unchanged) are emitted (out py37-40 / wl61-64) and match notes.md §4 (lines 189-201) — folded into the deliverables above as the `dTheta→0` specialization; all MATCH.

INTERNAL (scaffolding, no finding expected): `Mmat` rows, `inv = Inverse[Mmat]`, `solVec`, `fwdLaws`, the `Total[(...)^2]` transcription residual, all `expect_zero`/`expectZero` residual prints (all 0), PASS flags, banner strings.

## Self-test notes

Checked: (1) Variable independence — no `sp.diff`/`D[]` in this stage (pure linear-algebra inversion), so the zero-derivative trap is N/A. (2) Symmetry/parity — no integrals; N/A. (3) Trivial-case pre-check — substituting the boxed solution back into eq1–eq4 gives identical LHS=RHS (forward residuals 0, confirmed in both outputs); the matrix `Inverse[Mmat].{...}` route reproduces the same boxed forms. (4) Path specs — N/A (no missing-script finding). (5) Paper round-trip — N/A (no fix prescribed). The solve-vs-boxed and forward-transcription checks are non-tautological because the boxed RHS / `fwdLaws` are typed independently of the system being solved.
