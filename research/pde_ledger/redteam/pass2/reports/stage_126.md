---
unit_id: 126
batch: IV.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage126_positive_source_families.md]
  paper_appendix: present
---

# Audit unit 126 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_126.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage126_positive_source_families.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (present; only an `\input{stages/stage_126}` line at :1286 — no separate summary row for this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage126_positive_source_families_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage126_positive_source_families_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage126_positive_source_families_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage126_positive_source_families_mathematica_audit.txt`

## What the paper claims

The stage establishes that the physical lower-compensated Family-1 branch is reached by a simple, exact family of **positive** localized mouth sources — no sign-changing (exotic) profile is needed. Specifically: (1) the self-matched D/N derivative profile `σ_match(z)=k cos(kz)`, `k=π/(2L)`, is normalized to 1 and has exact mouth-bias `𝔤_match = π/4 ≈ 0.785398163397448`; the card body states verbatim "Self-matched derivative profile gives \(\mathfrak g=\pi/4\); a positive convex family hits the lower branch." (2) The exact lower branch from Stages 121–122 is `𝔤_-^{F1} ≈ 0.758035078944663`, giving overshoot `Δ𝔤_match ≈ 0.0273630844527852` and traction ratio `𝔗_-/𝔗_match = (π/4)/𝔤_-^{F1} ≈ 1.036097385480999` (3.61% enhancement). (3) The convex positive family `σ_ξ(z)=(1-ξ)k cos(kz)+ξ/L`, `0≤ξ≤1`, is normalized, nonnegative on `[0,L]`, has bias `𝔤_ξ=(1-ξ)π/4 + ξ·2/π`, and because `2/π < 𝔤_-^{F1} < π/4` there is a unique `ξ_* ∈ (0,1)` with `𝔤_{ξ_*}=𝔤_-^{F1}`; the closed form is `ξ_* = (-37√3 - 5π² + 2√(4107-100π²))/(5(8-π²)) ≈ 0.183918405511538` (an 18.4% admixture). The card `Checks` block demands positivity of the mouth source (no signed fitting), normalization, and that the compensation point be on the lower branch.

## What the script claims to verify

Both engines verify the same chain: (a) `σ_match` normalization `=1` and `𝔤_match=π/4`; (b) the imported lower-branch closed form `𝔤_-^{F1}=(2√(4107-100π²)-37√3)/(20π)`, with derived `Δ𝔤` and traction ratio printed; (c) `σ_ξ` normalization `=1` and `𝔤_ξ=(1-ξ)π/4 + 2ξ/π`; (d) **positivity of `σ_ξ` on the full box** `z∈[0,L], ξ∈[0,1]` — SymPy via affine-in-ξ (`∂²_ξ σ_ξ=0`) plus interval-minimum of both ξ-endpoint slices over `[0,L]`; Mathematica via the same affine-in-ξ check PLUS a direct `Resolve[ForAll[{z,xi}, 0≤z≤1 ∧ 0≤xi≤1, σ_ξ|_{lM→1} ≥ 0], Reals]`; (e) the unique compensation `ξ_*` solved from `𝔤_ξ=𝔤_-^{F1}` with residual `𝔤_ξ(ξ_*)-𝔤_- = 0`; (f) the interval `2/π < 𝔤_-^{F1} < π/4`.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `σ_match=k cos(kz)`, `k=π/(2L)`, normalized to 1 | py:16-17 `norm_match`; wl:32-33; assert `norm_match-1==0` py:22, `normMatch-1` wl:38 | match |
| `𝔤_match = π/4` | py:18,22; wl:34,39 (`gMatch-Pi/4`) | match |
| `𝔤_-^{F1} ≈ 0.758035...` (from 121–122) | py:25-26; wl:41-42 (imported closed form, numeric printed) | match |
| `Δ𝔤_match ≈ 0.027363...` | py:27,31; wl:43,47 (symbolic; numeric implied) | match |
| traction ratio `≈ 1.036097...` | py:28,35; wl:44,51 | match |
| `σ_ξ=(1-ξ)k cos(kz)+ξ/L`, normalized | py:39-40,45; wl:53-54,59 | match |
| `𝔤_ξ=(1-ξ)π/4 + 2ξ/π` | py:41 (`(8xi+π²(1-xi))/(4π)`); wl:55 | match |
| `σ_ξ ≥ 0` on `[0,L]×[0,1]` (no signed fitting) | py:61-77 affine+endpoint-min; wl:68-83 affine + `Resolve[ForAll]` | match |
| unique `ξ_* ≈ 0.183918...` on lower branch | py:79,83-86; wl:85,88 (`g_xi(xi_*)-g_-=0`) | match |
| `2/π < 𝔤_-^{F1} < π/4` (existence/uniqueness) | py:89-94; wl:90-94 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 22 | `norm_match-1==0` and `g_match-π/4==0` (raise on fail) | σ_match normalization + 𝔤_match=π/4 | yes |
| A2 | sympy | 45-46 | `norm_xi-1==0` | σ_ξ normalization | yes |
| A3 | sympy | 53-54 | `min σ_match on [0,L] == 0` | boundary min of σ_match | yes |
| A4 | sympy | 61-64 | `∂²_ξ σ_ξ == 0` (affine) | half of box-positivity argument | yes |
| A5 | sympy | 65-72 | `min σ_ξ(ξ=0) ≥ 0`, `min σ_ξ(ξ=1) ≥ 0` over [0,L] | endpoint-slice positivity (global-in-z) | yes |
| A6 | sympy | 74-77 | `min σ_ξ(ξ=0)==0`, `σ_ξ(ξ=1)==1/L` | pin endpoint minima/values | yes |
| A7 | sympy | 83-86 | `g_xi(xi_*) - g_- == 0` | compensation point on lower branch | yes |
| A8 | sympy | 91-94 | `2/π < g_- < π/4` | existence/uniqueness of ξ_* | yes |
| A9 | mathematica | 38-39 | `expectZero[normMatch-1]`, `expectZero[gMatch-π/4]` | σ_match normalization + 𝔤_match=π/4 | yes |
| A10 | mathematica | 59 | `expectZero[normXi-1]` | σ_ξ normalization | yes |
| A11 | mathematica | 68-70 | `expectZero[D[σ_ξ,{ξ,2}]]` (affine) | half of box-positivity argument | yes |
| A12 | mathematica | 71-73 | `expectZero[σ_match(lM)]` | boundary min of σ_match | yes |
| A13 | mathematica | 74-76 | `expectZero[σ_ξ(ξ=1)-1/lM]` | endpoint value | yes |
| A14 | mathematica | 77-83 | `Resolve[ForAll[{z,xi}, box, σ_ξ|_{lM→1} ≥ 0]] === True` | GLOBAL box-positivity | yes |
| A15 | mathematica | 88 | `expectZero[g_xi(xi_*)-g_-]` | compensation point on lower branch | yes |
| A16 | mathematica | 92-94 | `2/π < g_- < π/4` | existence/uniqueness of ξ_* | yes |

No assertion is tautological by construction: each is computed from the integral / quantifier-elimination engine output, not pre-substituted. Each maps to a stated paper deliverable.

## Findings

None. Adversarial attacks attempted are itemized in the Verdict justification.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration on the load-bearing claim. The two scripts share variable naming and the (unavoidable) same symbolic integrals for normalization / `𝔤_match` / `𝔤_ξ` — those are the physics, not a portable choreography. The positivity claim, which is the stage's adversarial heart ("do not use signed fitting to reach the upper branch"), is verified by **genuinely different methods**:

- SymPy (py:51-77): interval extremum search `sp.calculus.util.minimum(·, z, Interval(0,L))` on each ξ-endpoint slice, combined with the second-derivative affine argument `∂²_ξ σ_ξ = 0` to reduce the box minimum to a ξ-endpoint.
- Mathematica (wl:77-83): in addition to the affine check, a direct quantifier-elimination `Resolve[ForAll[{z,xi}, 0≤z≤1 ∧ 0≤xi≤1, (σ_ξ /. lM→1) ≥ 0], Reals]`, which decides the whole box at once rather than via endpoint slices.

These are different proof strategies for the same fact (one reduces-then-bounds, the other decides the joint quantifier), so the second engine is an independent corroboration. Other steps (`ξ_*` linear solve, interval check) are short enough that any port looks similar regardless of independence; they are not load-bearing for the transliteration question. Call: **not a transliteration**.

## Engine cross-check

Final outputs agree to the printed precision:

| quantity | sympy output | mathematica output |
|---|---|---|
| `𝔤_match` | `pi/4`, `0.78539816339744830962` | `Pi/4`, `0.78539816339744830961566...` |
| `𝔤_-^{F1}` | `(-37√3+2√(4107-100π²))/(20π)`, `0.75803507894466282692` | same closed form, `0.75803507894466282691968...` |
| traction ratio | `1.0360973854809996332` | `1.03609738548099963321630...` |
| `𝔤_ξ` | `(8xi+π²(1-xi))/(4π)` | `-(π(-1+xi))/4 + 2xi/π` (same expansion) |
| `ξ_*` | `(-37√3-5π²+2√(4107-100π²))/(5(8-π²))`, `0.18391840551153962831` | `1+(40+37√3-2√(4107-100π²))/(5(-8+π²))`, `0.18391840551153962830834...` (algebraically identical; verified by rearrangement and matching 20-digit numerics) |
| box-positivity | endpoint minima `0` and `1/L`, both ≥0; affine d²=0 | `Resolve[ForAll]` → `True` |
| compensation residual | `0` | `PASS: g_xi(xi_*) - g_-` (`=0`) |
| interval | `True` | `True` |

No engine disagreement.

## Verdict justification

**clean.** I confirmed I read the paper card, the notes, and the appendix (`\input` line only) and that the script's verified claim matches the paper. Attacks tried and survived: (1) **Box positivity is genuinely global, not corner-sampled.** SymPy proves it via affine-in-ξ (`∂²_ξσ_ξ=0`, py:61) + continuous interval minima of both ξ-endpoint slices (py:65-72) — a sound reduce-then-bound argument that is global in z; Mathematica proves it via `Resolve[ForAll]` over the full `(z,ξ)` box (wl:77-83), output `True` (wl-out:30). Both are global, not sampled. (2) **The `lM→1` substitution in the Mathematica quantifier is legitimately scale-covariant.** Writing `z=lM·u`, `u∈[0,1]`, gives `σ_ξ(z;lM)=(1/lM)[(1-ξ)(π/2)cos(πu/2)+ξ]`, and since `lM>0` the sign equals the `lM→1` bracket — so box-positivity for all `lM>0` ⇔ box-positivity at `lM=1`. The script even re-attaches `lM>0` via `Simplify[..., lM>0]` (wl:81). This matches the directive-authorized historical strengthening. (3) **Not tautological:** every assertion compares an integral/quantifier engine result against an independently-stated target (`π/4`, `1`, `0`, the imported `𝔤_-^{F1}`), none is `x==x`. (4) **The compensation residual `g_xi(xi_*)-g_- = 0` is non-trivial** because `xi_*` is solved from `g_xi=g_-` and then re-substituted into the integral-derived `g_xi`, exercising the full closed form. (5) **Symbol domains are correct:** `z,L` positive real (py:10), `lM>0` real and `z,xi` real (wl:29); positivity assumptions match the physical `[0,L]×[0,1]` setup, no over-strong assumption hides a branch. (6) **Outputs are fresh** (15:27 > 14:31/14:32 script mtimes). (7) **No hardcoded answer**: the only literal closed form (`𝔤_-^{F1}`) is an explicitly-labeled upstream import from Stages 121–122 (notes §2), not an in-stage fabricated number, and it is exercised (it must lie in `(2/π, π/4)` and produce `ξ_*∈(0,1)`). All deliverable values reconcile with the notes (see below).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 9 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `σ_match = k cos(kz), k=π/(2L)` | py:16, wl:32; out implicit | notes:12-16 (boxed) | MATCH |
| `𝔤_match = π/4 ≈ 0.785398163397448` | py:18,33 / wl:34,49; sympy-out:6,11 | notes:36-37 (boxed); tex:16 ("\(\mathfrak g=\pi/4\)") | MATCH |
| `𝔤_-^{F1} ≈ 0.758035078944663` | py:26,34 / wl:42,50; sympy-out:8,12 | notes:47 | MATCH |
| `Δ𝔤_match ≈ 0.0273630844527852` | py:27,31 / wl:43,47; sympy-out:9 (symbolic) | notes:53 | MATCH |
| traction ratio `≈ 1.036097385480999` | py:28,35 / wl:44,51; sympy-out:10,13 | notes:62 (boxed) | MATCH |
| `σ_ξ = (1-ξ)k cos(kz)+ξ/L` | py:39 / wl:53 | notes:71-76 (boxed) | MATCH |
| `𝔤_ξ = (1-ξ)π/4 + 2ξ/π` | py:41 / wl:55; sympy-out:16 | notes:82 | MATCH |
| `ξ_* = (-37√3-5π²+2√(4107-100π²))/(5(8-π²)) ≈ 0.183918405511538` | py:79,81 / wl:85,87; sympy-out:21,22 | notes:100-101 (boxed) | MATCH |
| interval `2/π < 𝔤_-^{F1} < π/4` | py:91-92 / wl:92-93; sympy-out:25-27 | notes:87 | MATCH |

INTERNAL (scaffolding, no finding expected in prose): `norm_match=1`, `norm_xi=1` (normalization residuals → assertion drivers), `min σ_match on [0,L]=0`, `min σ_ξ(ξ=0)=0`, `σ_ξ(ξ=1)=1/L`, `∂²_ξσ_ξ=0`, `g_xi(xi_*)-g_-=0` (residual), `Resolve[ForAll]→True` (positivity verdict flag), all PASS flags. These exist only to drive assertions and are not stage deliverables; though several are also reflected in the notes' boxed `σ_ξ` / normalization statements as MATCH-quality corroboration.

All emitted deliverable values reconcile with the notes; the terse `.tex` card legitimately carries only `𝔤=π/4` plus the prose summary, with the remaining deliverables anchored in the notes (per the augmentation guard, notes coverage = MATCH).

## Self-test notes

Checked the augmentation traps: (1) **Variable independence** — the only `D[·,{ξ,2}]` / `sp.diff(·,xi,2)` is on `σ_ξ`, which genuinely depends on `ξ` (the `(1-ξ)` and `ξ/L` terms), so `=0` is the substantive affine fact, not a degenerate zero; the `Resolve[ForAll]` quantifies over both `z` and `ξ`, both of which `σ_ξ` depends on. (2) **Parity/symmetry** — integrals are over the bounded `[0,L]`, `cos²` (even) and `cos` are non-vanishing positive-area integrands there; the printed nonzero `𝔤_match=π/4`, `g_ξ` are consistent (no spurious-vanishing claim). (3) **Trivial-case** — `σ_ξ` at `ξ=0` gives `k cos(kz)` (min 0 at z=L), at `ξ=1` gives constant `1/L>0`; both nonneg, consistent with the asserted endpoint minima; the box-positivity verdict `True` is the correct global answer since both summands are nonneg on the domain. No directive is written (zero findings).
