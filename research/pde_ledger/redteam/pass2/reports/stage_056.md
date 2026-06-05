---
unit_id: 056
batch: III.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md]
  paper_appendix: present
---

# Audit unit 056 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_056.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage056_transport_source_asymmetry.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 90; `\input` at 230)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.txt`

## What the paper claims

Stage 056 gives the previously-abstract source-asymmetry parameter `alpha` a physical transport origin. From the conservative 1D drift-diffusion law `∂_t σ + ∂_s J = 0`, `J = -D_σ ∂_s σ + v_σ σ` on the finite throat `s ∈ [0,L]`, the stationary zero-flux branch `J=0` collapses to `D_σ σ' = v_σ σ`, giving `σ ∝ e^{v_σ s/D_σ}`. The two boxed `\stagefield{Output}` deliverables are: (1) the transport interpretation `\boxed{Pe = v_σ L / D_σ}` (eq. app-stage056-Pe) — i.e. the Stage-053 exponential family is exactly this stationary branch with `α = Pe`; and (2) the physical overlap factor `\boxed{Ω_Pe = π Pe(2 Pe e^{Pe} + π) / [(4 Pe² + π²)(e^{Pe} − 1)]}` (eq. app-stage056-Omega-Pe). The notes add further (non-boxed but stated) deliverables: normalized profile `σ_Pe(s) = Pe e^{Pe s/L}/(e^{Pe}−1)`, uniform overlap `I_W = 2√(2L)/π`, transport overlap `I_Pe = 2√(2L)·Pe(2Pe e^{Pe}+π)/[(4Pe²+π²)(e^{Pe}−1)]`, the covariance monotonicity identity `dΩ_Pe/dPe = Cov_Pe(χ₀,s)/I_W`, endpoints `Ω_Pe(0)=1` and `lim_{Pe→∞} Ω_Pe = π/2`, the small-Pe slope `Ω_Pe = 1 + ((4−π)/(2π))Pe + O(Pe²)`, and the large-Pe asymptotic `Ω_Pe = π/2 − π³/(8Pe²) + O(Pe⁻³)`. The appendix row 90 summarizes the stage as "Drift-diffusion zero-flux branch and Peclet-number identification."

## What the script claims to verify

The SymPy docstring lists four checks: (1) the zero-flux drift-diffusion branch gives the normalized exponential family; (2) the exact D/N lowest-mode overlap reproduces the boost formula with `α = Pe`; (3) `dΩ/dPe` equals the exact covariance; (4) small/large-Pe asymptotics. Concretely the assertions verify: zero-flux residual `−D_σ σ' + v_σ σ = 0` for `σ = e^{v_σ s/D_σ}`; normalization `∫₀ᴸ σ_Pe ds = L`; `Ω_Pe` (computed by integrating `σ_Pe·χ₀` and dividing by `I_W`) equals the closed-form expected formula; endpoint limits 1 and π/2; the covariance identity `dΩ/dPe − Cov/I_W = 0`; small-Pe linear coefficient `= (4−π)/(2π)`; and large-Pe asymptotic through `O(Pe⁻²)`. The Mathematica script verifies the identical set independently, additionally pinning the large-Pe coefficient via `Limit[pe²(omegaPe − π/2), pe→∞] = −π³/8`.

## Paper ↔ script cross-check

| Paper / notes deliverable | Script-side check | Status |
|---|---|---|
| Zero-flux branch ⇒ exponential `σ ∝ e^{v_σ s/D_σ}` (eq. transport) | py:38-40 / wl:32-34 `expect_zero(J)` | match |
| `Pe = v_σ L/D_σ` identification (boxed Output 1) | implicit: `σ_Pe = Pe e^{Pe s/L}/(e^{Pe}−1)` py:43 / wl:36 with `α=Pe` | match |
| Normalization `∫σ_Pe = L` | py:44-45 / wl:37,50 | match |
| `I_W = 2√(2L)/π` | py:48-50 / wl:38-39,47 (printed) | match (printed, then used as denominator — exercised) |
| `Ω_Pe` boxed formula (Output 2) | py:57-61 / wl:42-51 `expect_zero(Ω_Pe − expected)` | match |
| Endpoints `Ω_Pe(0)=1`, `lim=π/2` | py:64-71 / wl:53-58 | match |
| Covariance monotonicity `dΩ/dPe = Cov/I_W` | py:73-79 / wl:60-65 | match |
| Small-Pe slope `(4−π)/(2π)` | py:82-88 / wl:67-78 | match |
| Large-Pe `−π³/(8Pe²)` | py:90-95 / wl:69-82 | match |

`paper_alignment: aligned` — every boxed and stated deliverable has a faithful, non-tautological script-side check, in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero(J)` where `J=-D σ' + v σ`, `σ=e^{v s/D}` | zero-flux branch ⇒ exp profile | yes |
| A2 | sympy | 45 | `expect_zero(∫σ_Pe ds − L)` | normalization | yes |
| A3 | sympy | 61 | `expect_zero(Ω_Pe − expected)` (Ω_Pe from integral) | Output 2 (Ω_Pe formula) | yes |
| A4 | sympy | 68-71 | `if Ω0≠1 raise`; `if simplify(ΩInf−π/2)≠0 raise` | endpoints | yes |
| A5 | sympy | 79 | `expect_zero(dΩ/dPe − Cov/I_W)` | covariance monotonicity | yes |
| A6 | sympy | 85-88 | `expect_zero(small-Pe coeff − (4−π)/(2π))` | small-Pe slope | yes |
| A7 | sympy | 92-95 | `expect_zero(Ω_large − (π/2 − π³/(8Pe²)))` | large-Pe asymptotic | yes |
| A8 | math | 34 | `expectZero[jFlux]` | zero-flux branch | yes |
| A9 | math | 50 | `expectZero[norm − ell]` | normalization | yes |
| A10 | math | 51 | `expectZero[omegaPe − omegaExpected]` | Output 2 | yes |
| A11 | math | 57-58 | `expectZero[omega0−1]`, `expectZero[omegaInf−π/2]` | endpoints | yes |
| A12 | math | 65 | `expectZero[D[omegaPe,pe] − cov/iW]` | covariance | yes |
| A13 | math | 75-78 | `expectZero[small coeff diff]` | small-Pe slope | yes |
| A14 | math | 79-82 | `expectZero[largeCoeff + π³/8]` | large-Pe coeff | yes |

No tautological rows: every `expect_zero`/`expectZero` is of the form `derived_quantity − independently_stated_target`, where the derived side is built by integration/differentiation (not by reusing the target). A3/A10 in particular integrate `σ_Pe·χ₀` and divide by an independently integrated `I_W`, then compare against a separately typed-in closed form — a genuine check, not `x==x`.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.txt:11,27` (mtime 2026-05-11 12:43 < script 2026-06-03 15:59)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage056_transport_source_asymmetry_mathematica_audit.txt:11` (mtime 2026-05-11 12:53 < script 2026-06-03 15:59)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py:2,3,31,97` (residual stale `Stage 39`/`STAGE 56` self-labels in current source)

**What's wrong:**
Both saved output `.txt` files predate the current scripts (script mtime 2026-06-03, outputs 2026-05-11), so they are stale by mtime. The captured banners reflect the pre-renumber state: sympy output line 11 prints `STAGE 39 — …` and line 27 `Stage 39 audit passed.`; mathematica output line 11 prints `STAGE 039 — …`. The current `.wl` source banner (line 26) correctly says `STAGE 056`, so the `.wl` output is a pure refresh-lag label mismatch. The current `.py` source, however, still carries stale self-labels: docstring line 2 `Stage 39 SymPy audit`, line 3 begins the stage-39 description, banner line 31 `banner("STAGE 56 — …")` (also non-3-digit, and the output shows the older `STAGE 39`), and line 97 `print("\nStage 39 audit passed.")`. All numeric *results* in both transcripts match the current scripts exactly (I_W, I_Pe, Ω_Pe, endpoints, series, covariance residual all identical), so the staleness is label-only, not a content disagreement.

**Why this matters:**
Refreshing the outputs realigns the committed transcripts with the corrected `056` banners. The residual `Stage 39`/`STAGE 56` self-labels in the `.py` source are the known low-severity numbering-drift class ([[numbering-drift-root-cause]]); they are non-blocking but should be normalized to `056` so the refreshed output no longer reprints `STAGE 39`.

**Required change:**
Self-label normalization in the `.py` source (Reading-2 in-loop scope): line 2 `Stage 39 SymPy audit` → `Stage 056 SymPy audit`; line 31 `banner("STAGE 56 — …")` → `banner("STAGE 056 — …")`; line 97 `print("\nStage 39 audit passed.")` → `print("\nStage 056 audit passed.")`. (Line 3's prose "Stage 39" descriptor likewise → 056.) Do not alter any math. The orchestrator re-runs both scripts to refresh the `.txt` transcripts.

**Verification:**
After refresh, sympy output line 11 reads `STAGE 056 — …` and the closing line reads `Stage 056 audit passed.`; mathematica output line 11 reads `STAGE 056 — …`; all result lines (I_W, I_Pe, Omega_Pe, series, residuals) unchanged from the current transcripts; both scripts exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. It uses Mathematica-native machinery throughout and differs structurally from the `.py` in several load-bearing ways: it pins the large-Pe coefficient via `largeCoeff = Limit[pe^2 (omegaPe - Pi/2), pe -> Infinity]` and asserts `largeCoeff + Pi^3/8 == 0` (wl:70,79-82), whereas SymPy uses `sp.series(..., Pe, oo, 3)` and subtracts the full asymptotic (py:90-95). The `.wl` builds `omegaPe` from `Integrate[sigmaPe chi0, {s,0,ell}]/iW` with explicit `$Assumptions` (Reals, ell>0, dSigma>0, vSigma>0, pe>0) rather than SymPy's per-symbol `positive=True`. Endpoints use `Limit[..., pe->0]`/`pe->Infinity` (wl:53-54) vs `sp.limit` (py:64-65). The covariance block (wl:60-65) and small-Pe `Series[..., {pe,0,2}]` (wl:67) mirror the same mathematics but are expressed with native Mathematica functions and an extra printed coefficient — this is parallel derivation of the same identities, which is exactly what the two-engine policy requires. No `mathematica_transliteration` finding.

## Engine cross-check

The two transcripts agree on every result. SymPy: `I_W = 2√2·√L/π`, `Omega_Pe = π·Pe(2Pe·e^{Pe}+π)/((4Pe²+π²)(e^{Pe}−1))`, `Omega_Pe(0)=1`, `lim=π/2`, small series `1 − Pe/2 + 2Pe/π + …`, large `π/2 − π³/(8Pe²)`, covariance residual 0. Mathematica: `I_W = 2√2·√ell/Pi`, `Omega_Pe = pe·Pi(2 e^{pe} pe + Pi)/((−1+e^{pe})(4pe²+Pi²))` (algebraically identical), `Omega_Pe(0)=1`, `lim=Pi/2`, same small/large series, `largeCoeff = −Pi³/8`, covariance residual 0. The Mathematica `Limit::alimv` warnings are benign (it ignores the `pe>0` assumption while taking the `pe→0`/`pe→∞` limit; the limits are still computed correctly and match SymPy). `engines_agree: true`.

## Verdict justification

The stage holds up. I attacked: (a) the zero-flux assertion — `σ = e^{v s/D}` is defined and `J` is recomputed by differentiation, so A1/A8 are non-tautological (the check would fail for a wrong flux sign or coefficient). (b) The `Ω_Pe` check — the derived side is built by genuine integration of `σ_Pe·χ₀` over `[0,L]` and division by the independently integrated `I_W`, then compared to a separately written closed form; this would catch a wrong constant or sign. (c) Symbol domains — `positive=True`/`pe>0` is justified by the paper's constructive branch `Pe ≥ 0` and `D_σ, v_σ, L > 0`; the script audits only the constructive branch, which matches the paper's `Pe ≥ 0` monotonicity claim, so this is not `missing_branch` (the `Pe<0` destructive branch is descriptive context in the notes, not a script deliverable, and the boxed Output formula is the same analytic expression). (d) The covariance identity is a real differential check, not reusing `Ω`'s closed form. (e) The asymptotic checks pin both the linear small-Pe slope and the `π³/8` large-Pe coefficient. All boxed and stated deliverables reconcile to the docs. The only finding is the label-only `stale_output` (both transcripts predate the `056` banner; the `.py` source still carries `Stage 39`/`STAGE 56` self-labels). No math defect; `material_change: false`.

## Self-test notes

Checked: (1) variable independence — `D[omegaPe, pe]`/`sp.diff(Omega_Pe, Pe)` (A5/A12): `Omega_Pe` genuinely depends on `Pe`, so the derivative is non-trivial; the covariance side is also `Pe`-dependent, so neither side is identically zero — no spurious-zero trap. (2) Symmetry/parity — integrals are over the bounded `[0,L]`, not unbounded, so no odd/even cancellation trap; `I_W = 2√(2L)/π > 0` confirms the denominator is nonzero. (3) Trivial-case — at `Pe→0`, `σ_Pe → 1/L·…` recovers the uniform source and `Ω_Pe(0)=1` (verified in both transcripts), the correct trivial limit. (4)/(5) — the only fix is label normalization plus an output refresh; it touches no math and introduces no new paper_misalignment (no constant changes).

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| zero-flux residual `= 0` | py:40 / wl:34; sympy.txt:13, math.txt:13 | tex:20-24 (J=0 branch); md:65-71 | MATCH |
| normalization `∫σ_Pe ds − L = 0` | py:45 / wl:50; sympy.txt:14, math.txt:18 | md:73-84 (`∫₀ᴸ σ = L`) | MATCH |
| `Pe = v_σ L/D_σ` | py:43 (σ_Pe defn) / wl:36; (def, in σ_Pe) | tex:26-30 boxed; md:84,199 | MATCH |
| `σ_Pe(s) = Pe e^{Pe s/L}/(e^{Pe}−1)` | py:43 / wl:36 | tex:nil; md:79-80,193-195 | MATCH (md) |
| `I_W = 2√(2L)/π` | py:49-50 / wl:39,47; sympy.txt:15, math.txt:15 | md:106 (`I_W = 2√(2L)/π`) | MATCH (md) |
| `I_Pe = 2√(2L)·Pe(2Pe e^{Pe}+π)/[(4Pe²+π²)(e^{Pe}−1)]` | py:52-54 / wl:40,48; sympy.txt:16, math.txt:16 | md:114-115 | MATCH (md) |
| `Ω_Pe = π Pe(2Pe e^{Pe}+π)/[(4Pe²+π²)(e^{Pe}−1)]` | py:53-61 / wl:41-51; sympy.txt:17, math.txt:17 | tex:33-37 boxed; md:119-121,203-204 | MATCH |
| `Ω_Pe(0) = 1` | py:64,68 / wl:53,57; sympy.txt:19, math.txt:26 | md:163,208-209 | MATCH (md) |
| `lim_{Pe→∞} Ω_Pe = π/2` | py:65,70 / wl:54,58; sympy.txt:20, math.txt:27 | md:165,211-212 | MATCH (md) |
| `dΩ/dPe = Cov/I_W` (residual 0) | py:79 / wl:65; sympy.txt:21, math.txt:32 | md:139-143 | MATCH (md) |
| small-Pe slope `(4−π)/(2π)` | py:84-88 / wl:68,75-78; sympy.txt:22-23, math.txt:38,41 | md:175-176 | MATCH (md) |
| large-Pe coeff `−π³/8` (i.e. `−π³/(8Pe²)`) | py:90-95 / wl:70,79-82; sympy.txt:24-25, math.txt:39-40,43 | md:182-183 | MATCH (md) |

INTERNAL (scaffolding, no finding): `sigma = e^{v_σ s/D_σ}` (intermediate), `J` flux residual expression, `p_Pe = σ_Pe/L`, `E_Pe[χ₀]`, `E_Pe[s]`, `E_Pe[χ₀ s]`, `Cov` intermediate, `Omega_expected`/`omegaExpected` (the typed-in comparison target), per-check `PASS`/residual flags, banners.

reconciliation: complete; 12 values checked, 0 misaligned
