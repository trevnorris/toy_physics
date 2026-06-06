---
unit_id: 114
batch: IV.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [notes/stages/moving_throat_pde_stage114_concrete_core_schur.md]
  paper_appendix: present
---

# Audit unit 114 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_114.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage114_concrete_core_schur.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 28: "Stages 114--124: two-channel core realization, parent extraction, and geometric core selection."; stage card `\input` at line 1262)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.txt`

## What the paper claims

The card is terse; the notes are authoritative. Stage 114 introduces a concrete two-channel core model: a 2×2 linear core system `[[K_s, λ],[λ, −K_q·D_W^bare(z)]] (s,q)ᵀ = u·(g_s,g_q)ᵀ`, with mouth feedback `δΛ_core·u = g_s·s + g_q·q`. The card's bottom line (quote line 16): "Concrete shell/mixed core Schur complement reproduces the reduced Robin–mixed outlet." The notes box four distinct deliverables: (1) the exact Schur-complement outlet `δΛ_core = g_s²/K_s − (K_s g_q − λ g_s)²/[K_s(K_s K_q D_W^bare + λ²)]`; (2) its rewrite into the reduced Stage-112 form `ρ_c − σ_c/(1 − κ_c z² − iγ_c z⁵) + O(z⁶)`; (3) the identifications `ρ_c = g_s²/K_s` and `σ_c = (K_s g_q − λ g_s)²/[K_s² K_q (1+r_c)]`; (4) `κ_c = κ_0/(1+r_c)`, `γ_c = γ_0/(1+r_c)`, with `r_c := λ²/(K_s K_q)`. Card `Checks` (lines 21-25) name the Schur-sign check, an L/a normalization check, and an overlap-ratio check; only the Schur reduction is in genuine algebraic scope here (the L/a and overlap items are upstream-geometry sanity reminders, not equations this stage emits).

## What the script claims to verify

Both scripts build the symbolic 2×2 matrix `M = [[K_s, λ],[λ, −K_q·D]]` and vector `c = (g_s, g_q)`, compute `δΛ(D) = Apart(cᵀ·M⁻¹·c, D)` via the engine's built-in matrix inverse, and assert two identities. F-check 1 ("Schur form identity"): `δΛ(D) − (ρ_c − σ̃/(D + r_c)) == 0`, with `ρ_c = g_s²/K_s`, `σ̃ = (K_s g_q − λ g_s)²/(K_s² K_q)`, `r_c = λ²/(K_s K_q)`. F-check 2 ("low-frequency normalized outlet identity"): substitute `D → D_bare = 1 − κ_0 z² − iγ_0 z⁵` into `δΛ(D)` and assert it equals `ρ_c − σ_c/(1 − κ_c z² − iγ_c z⁵)`, with `σ_c = σ̃/(1+r_c)`, `κ_c = κ_0/(1+r_c)`, `γ_c = γ_0/(1+r_c)`. Both scripts then print the four identifications. The assertions are non-tautological: the left side comes from an actual matrix inverse, the right side is an independently-written closed form.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) Exact Schur outlet `δΛ_core(D)` closed form | F-check 1: `delta_D − (ρ_c − σ̃/(D+r_c)) == 0` | match (algebraically equals the notes' boxed form: `σ̃/(D+r_c) = (K_s g_q − λ g_s)²/[K_s(K_s K_q D + λ²)]`) |
| (2) Rewrite into reduced Stage-112 form in z | F-check 2: `delta_z − (ρ_c − σ_c/(1−κ_c z²−iγ_c z⁵)) == 0` | match |
| (3) `ρ_c`, `σ_c` identifications | printed (py L52-53 / wl L56-57); `ρ_c` load-bearing in both asserts; `σ_c` load-bearing in F-check 2 | match |
| (4) `κ_c = κ_0/(1+r_c)`, `γ_c = γ_0/(1+r_c)` | load-bearing in F-check 2 (substituted into target_z); printed py L54-55 / wl L58-59 | match |

`paper_alignment: aligned`. Every notes-boxed deliverable maps to a load-bearing, non-tautological script check.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42 | `expect_zero("Schur form identity", delta_D - target_D)` | claim 1 (exact Schur outlet) | yes |
| A2 | sympy | 49 | `expect_zero("low-frequency normalized outlet identity", delta_z - target_z)` | claims 2/3/4 (reduced z-form + κ_c,γ_c,σ_c) | yes |
| A3 | mathematica | 47 | `expectZero["Schur form identity", deltaD - targetD]` | claim 1 | yes |
| A4 | mathematica | 52 | `expectZero["low-frequency normalized outlet identity", deltaZ - targetZ]` | claims 2/3/4 | yes |

A1/A3 are genuine: `delta_D` is `Apart(cᵀ M⁻¹ c, D)` from the actual inverse, `target_D` is an independent hand-written rational form; their difference vanishing exercises the Schur elimination. A2/A4 substitute the concrete `D_bare(z)` and exercise the `1/(1+r_c)` renormalization of `σ,κ,γ`. None is tautological (no `x := expr; assert x==expr`). No hardcoded numeric answer; all symbolic. Both `expect_zero`/`expectZero` are can-fail (a wrong `target` would leave a nonzero rational residual).

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage114_concrete_core_schur_mathematica_audit.wl:33-52`
- (parallel) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage114_concrete_core_schur_sympy_audit.py:27-49`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`, not an independent re-derivation of the Schur complement. Corresponding sections, in identical order:

- Matrix/vector construction — sympy L27-28: `M = sp.Matrix([[K_s, lam],[lam, -K_q*D]])`, `c = sp.Matrix([g_s, g_q])`; wl L33-34: `m = {{kS, lam}, {lam, -kQ*dSym}}`, `c = {gS, gQ}` — same matrix, same vector.
- Schur reduction — sympy L30: `delta_D = sp.apart((c.T * M.inv() * c)[0], D)`; wl L36: `deltaD = FullSimplify[Apart[c.Inverse[m].c, dSym], ...]` — **same choreography: built-in matrix inverse of the SAME matrix, then `apart`/`Apart` partial-fraction in the SAME variable.** The Mathematica route does not eliminate `q` independently (no `Solve` on the 2×2 system, no explicit `adj/det` Schur formula, no Gaussian row reduction); it re-types the SymPy strategy.
- Definitions and targets — sympy L34-41 define `rho_c, r_c, sigma_tilde, sigma_c, kappa_c, gamma_c, target_D`; wl L39-46 define `rhoC, rC, sigmaTilde, sigmaC, kappaC, gammaC, targetD` in the same order with identical formulas.
- Assertions — sympy L42 then L49 (`Schur form identity`, then `low-frequency normalized outlet identity`); wl L47 then L52 — same two checks, same order, same names.

The only "independence" is that two CAS each implement matrix inversion. Per the IV.2 transliteration-watch mandate (105-175 first-pass orchestrator-direct band), a mirrored matrix-inverse + identical extraction + mirrored assertion order is transliteration even though the result is correct.

**Why this matters:**
The second-engine policy requires the `.wl` to corroborate the Schur reduction by an *independent* path so a shared algebra/transcription error cannot pass both engines. As written, a mistake in the `target_D`/`target_z` construction (the hand-written closed form being checked) would be replicated identically in both files and pass both. The independent value of the Mathematica run is reduced to "Mathematica's `Inverse` agrees with SymPy's `M.inv()`," which is not what the cross-engine check is for.

**Required change:**
Re-derive the core correction in the `.wl` by an independent route, then assert it equals the same `targetD`/`targetZ`. Concretely, replace the built-in `Inverse[m]` Schur path (wl L36) with an explicit elimination of the side-channel variable `q`: `Solve` the second core equation `lam*s - kQ*dSym*q == gQ` for `q`, back-substitute into the first equation `kS*s + lam*q == gS` to get `s`, then form `deltaD = gS*s + gQ*q` (the mouth feedback `δΛ·u = g_s s + g_q q` with `u=1`). Keep the two `expectZero` assertions and the printed identifications unchanged. This makes the Mathematica route derive `δΛ` from the linear-system premise rather than re-typing the matrix-inverse choreography. (If the user prefers, an equivalent independent route is the explicit 2×2 Schur formula `δΛ = g_sᵀ·(K_s)⁻¹·g_s − (Schur term)` written out by hand rather than via `Inverse`.)

**Verification:**
After the rewrite, the `.wl` must contain a `Solve[...]` (or explicit `adj/det`) elimination of `q` instead of `Inverse[m]` at line ~36; both `expectZero` checks must still print `= 0` / `PASS`; the printed `delta_Lambda(D)`, `rho_c`, `sigma_c`, `kappa_c`, `gamma_c` must be byte-identical to the current Mathematica output transcript (lines 5, 13-15). Re-run with `math -script` to confirm Exit[0].

## Independent-derivation check (Mathematica)

Not independent — see F1. The `.wl` mirrors the `.py` step for step: identical matrix `m`, identical `Apart[c.Inverse[m].c, dSym]` versus `sp.apart((c.T*M.inv()*c)[0], D)`, identical six coefficient definitions in identical order, identical two assertions in identical order. The Schur complement is obtained from the engine's built-in inverse in both, so the Mathematica run does not provide an algebraically independent corroboration of the elimination.

## Engine cross-check

Both engines agree at the level they claim. SymPy `delta_Lambda(D) = g_s²/K_s − (K_s g_q − g_s λ)²/[K_s(D K_q K_s + λ²)]` (output lines 5-10) equals the Mathematica `(dSym gS² kQ − gQ² kS + 2 gQ gS λ)/(dSym kQ kS + λ²)` (output line 5) after combining over a common denominator. Both report `Schur form identity = 0` and `low-frequency normalized outlet identity = 0` and pass. The four printed identifications match across engines: `rho_c = g_s²/K_s`; `sigma_c = (K_s g_q − g_s λ)²/[K_s(K_q K_s + λ²)]`; `kappa_c = κ_0 K_q K_s/(K_q K_s + λ²)`; `gamma_c = γ_0 K_q K_s/(K_q K_s + λ²)`. No `engine_disagreement`. (Note these printed `sigma_c/kappa_c/gamma_c` are the `1/(1+r_c)` forms with `1+r_c = (K_q K_s + λ²)/(K_q K_s)` cleared — consistent with the notes' boxed identifications.)

## Verdict justification

`findings`, single finding F1 (`mathematica_transliteration`, medium). The math is correct and fully aligned with the paper/notes: I verified the script's `target_D = ρ_c − σ̃/(D+r_c)` algebraically equals the notes' boxed exact-Schur form, the z-form check exercises the `1/(1+r_c)` renormalization of `σ,κ,γ`, and all four boxed deliverables map to load-bearing non-tautological checks. Attacks that failed: (a) tautology — the LHS comes from a real matrix inverse, not a re-substituted definition, so a wrong `target` would leave a nonzero residual; (b) hardcoded — no numeric literals, fully symbolic; (c) symbol-assumption trap — positivity (`K_s,K_q,κ_0,γ_0>0`, `lam>0`) does not enable a false simplification because the residuals are identically-zero rational functions independent of sign, and `g_s,g_q` are unrestricted real; (d) single-point — checks are full symbolic identities in `D`/`z`, not point evaluations; (e) paper misalignment — none, card+notes+output all reconcile. The lone defect is that the `.wl` is a transliteration of the `.py` rather than an independent re-derivation, which is the specific risk this IV.2 band is being scrutinized for; it does not affect correctness but weakens the second-engine guarantee.

## Self-test notes

Checked: (1) variable independence — no `diff`/`D[]` derivative-against-absent-variable traps; the only substitution is `D → D_bare(z)` which legitimately introduces z-dependence. (2) symmetry/parity — no integrals; n/a. (3) trivial-case — the two `expect_zero` residuals are exact rational-function zeros (confirmed by the printed `= 0` in both transcripts and by hand-checking `σ̃/(D+r_c)` matches the notes' boxed denominator `K_s(K_s K_q D + λ²)`), so they are can-fail, not vacuous. (4) the F1 required change names the full `.wl` path under `mathematica/` and keeps the same `targetD/targetZ` so it cannot introduce a new `paper_misalignment` (paper round-trip clean — the rewrite changes only the derivation route, not any emitted value).

## Value Reconciliation (pass-2 augmentation)

Outputs are fresh (sympy .txt mtime 2026-05-27 15:18 > .py 15:08; mathematica .txt 15:24 > .wl 15:08), so reconciliation uses script source + committed transcripts.

Deliverable values the scripts emit (all symbolic closed forms; no numeric constants in this stage):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δΛ_core(D) = g_s²/K_s − (K_s g_q − λ g_s)²/[K_s(K_s K_q D + λ²)]` | py out L5-10 / wl out L5 (asserted py L42, wl L47) | notes L44-53 (boxed exact Schur) | MATCH |
| reduced z-form `ρ_c − σ_c/(1−κ_c z²−iγ_c z⁵)+O(z⁶)` | asserted py L49 / wl L52 | notes L60-69 (boxed) | MATCH |
| `ρ_c = g_s²/K_s` | py out L15 / wl out L13 | notes L72-75 (boxed) | MATCH |
| `σ_c = (K_s g_q − λ g_s)²/[K_s(K_q K_s + λ²)]` (= `…/[K_s² K_q (1+r_c)]`) | py out L16 / wl out L14 | notes L77-82 (boxed) | MATCH |
| `κ_c = κ_0 K_q K_s/(K_q K_s + λ²)` (= `κ_0/(1+r_c)`) | py out L17 / wl out L15 | notes L84-88 (boxed `κ_c=κ_0/(1+r_c)`) | MATCH |
| `γ_c = γ_0 K_q K_s/(K_q K_s + λ²)` (= `γ_0/(1+r_c)`) | py out L18 / wl out L16 | notes L84-88 (boxed `γ_c=γ_0/(1+r_c)`) | MATCH |
| `r_c = λ²/(K_s K_q)` | py L35 / wl L40 (intermediate, used in 1+r_c) | notes L55-58 (boxed `r_c := λ²/(K_s K_q)`) | MATCH |

INTERNAL scaffolding (no prose expected): `M`/`c` matrix-vector, `sigma_tilde` (= numerator-before-`(1+r_c)`), `target_D`/`target_z`/`targetD`/`targetZ` (comparison RHS), `D_bare`/`dBare`, pass/fail flags ("PASS", "= 0"), banner strings.

reconciliation: complete; 7 values checked, 0 misaligned
