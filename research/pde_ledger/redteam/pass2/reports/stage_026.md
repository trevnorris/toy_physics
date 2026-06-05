---
unit_id: 026
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md]
  paper_appendix: present
---

# Audit unit 026 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_026.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage026_concrete_axial_overlaps.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row at line 42)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage026_concrete_axial_overlaps_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage026_concrete_axial_overlaps_mathematica_audit.txt`

## What the paper claims

Stage 026 picks the first concrete finite-throat axial-mode family — a constant N/N zero mode `ν₀=1/√L` and the D/N half-wave ladder `f_n=√(2/L) sin((n+1/2)πs/L)` on `s∈[0,L]` — and computes the exact overlaps. `\stagefield{Output}` reads: "Stage~026 outputs the finite-throat D/N overlap \eqref{eq:app-stage026-kappa0}, the concrete couplings \eqref{eq:app-stage026-couplings}, the finite-throat normalization law \eqref{eq:app-stage026-normalization}, and the required stiffness formula \eqref{eq:app-stage026-Kreq}." The distinct deliverables are: (1) the general overlap `κ_n=√2/[(n+1/2)π]` and the boxed minimal value `κ:=κ₀=2√2/π`; (2) the concrete couplings `C=κλ_B, G_U=λ_U, G_W=κλ_W, R=κλ_R`; (3) the Stage-025-substituted quantities `Δ=Ω_U²Ω_W²−κ²λ_R²`, `Q=λ_U²Ω_W²+2κ²λ_Uλ_Wλ_R+κ²λ_W²Ω_U²`, `P=κ(Ω_U²λ_W+λ_Rλ_U)`, and `X=κ²λ_B²/ϖ²`; (4) the boxed branch-level normalization law `mhat_rad² κ²(Ω_U²λ_W+λ_Rλ_U)² / [Δ(KΔ−Δκ²λ_B²/ϖ²−Q)] = N_Q^target`; (5) the boxed three-term `K_req` formula; (6) `K=K_η+6T_Ω` and the branch condition `K_η+6T_Ω=K_req`. The notes additionally pin the target constant `N_Q^target = 54 G c_s⁵/(5 a⁵ c⁵)` (notes lines 7, 208, 214) and the numeric `κ≈0.900316316` (notes line 130); the `.tex` card carries the target only symbolically as `N_Q^{\rm target}`.

## What the script claims to verify

Both engines verify: the N/N and D/N modes are unit-normalized (`∫u0²=1`, `∫f0²=1`); the general overlap integral `∫₀^L u0 f_n ds` equals the closed form `√2/[(n+1/2)π]`, and at `n=0` equals `2√2/π` (with numeric `0.900316316157106`); the three reused overlap integrals collapse to `κ` and the eta self-overlap is `1`; the reduced coefficients `C,G_U,G_W,R,Δ,Q,P,B0,Z0,N0,D0` follow by substitution; the residual `mhat²·N0/D0 − target` is solved for `K` and the solver's root (a) matches the paper's three-term `K_req` decomposition and (b) back-substitutes to a vanishing residual; and `K=K_η+6T_Ω` is printed for the constant wall branch. The SymPy docstring states the scope verbatim; the comments above each assertion name the claim being exercised.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `κ_n=√2/[(n+1/2)π]` (eq:kappa-n) | py L105 `kappa_n − kappa_n_expected==0`; wl L82 FTC-path vs analytic | match |
| `κ=2√2/π` (boxed eq:kappa0) | py L109 `kappa − 2√2/π==0`; wl L88 | match |
| `C=κλ_B, G_U=λ_U, G_W=κλ_W, R=κλ_R` (boxed couplings) | py L148-151 / wl L116-119 (substitution; printed III.2) | match |
| `Δ,Q,P` (boxed DeltaQP) | py L153-155 / wl L121-123 (printed III.2) | match |
| `X=κ²λ_B²/ϖ²` (eq:X) → B0 | py L157 `B0=C²/varpi²`; wl L125 | match |
| Normalization law (boxed eq:normalization) | py L189 `residual=mhat²·N0/D0 − target`; wl L152 | match |
| `K_req` three-term (boxed eq:Kreq) | py L212 `K_req − K_req_paper==0` & L213 back-sub; wl L173-174 | match |
| `K=K_η+6T_Ω` (eq:Kgeo) | py L216 printed; wl L176 printed | match (printed, not asserted — same as paper, a definition) |
| `N_Q^target=54 G c_s⁵/(5 a⁵ c⁵)` (notes only) | py L188 `target=54 G c_s⁵/(5 a⁵ c⁵)`; wl L151 | match (carried in notes; `.tex` is symbolic) |
| `κ≈0.900316316` (notes) | py L110 numeric `0.900316316157106`; wl L89 | match |

`paper_alignment: aligned`. Every boxed deliverable maps to a non-tautological script assertion; the two prose-only constants (`N_Q^target`, numeric κ) reconcile against the notes.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 83 | `∫u0² − 1 == 0` | mode normalization (eq:u0) | yes |
| A2 | sympy | 84 | `∫f0² − 1 == 0` | mode normalization (eq:DN-ladder) | yes |
| A3 | sympy | 105 | `integrate(u0·f_n) − √2/[(n+½)π] == 0` | claim 1 (κ_n) | yes |
| A4 | sympy | 109 | `κ − 2√2/π == 0` | claim 1 (κ₀, boxed) | yes |
| A5 | sympy | 145 | `integrate(u0·f0) − κ == 0` | claim 2 (overlaps→κ) | yes |
| A6 | sympy | 146 | `integrate(u0·u0) − 1 == 0` | claim 2 (I_ηu=1) | yes |
| A7 | sympy | 212 | `K_req(solve) − K_req_paper == 0` | claim 5 (Kreq, boxed) | yes |
| A8 | sympy | 213 | `residual(K=K_req) == 0` | claim 4 (normalization law) | yes |
| B1 | math | 48 | `∫u0² − 1 == 0` | mode normalization | yes |
| B2 | math | 49 | `∫f0² − 1 == 0` | mode normalization | yes |
| B3 | math | 82 | `κ_n(analytic) − κ_n(FTC) == 0` | claim 1 (κ_n, independent path) | yes |
| B4 | math | 88 | `κ − 2√2/π == 0` | claim 1 (κ₀) | yes |
| B5 | math | 113 | `overlap_u0_f0 − κ == 0` | claim 2 | yes |
| B6 | math | 114 | `overlap_u0_u0 − 1 == 0` | claim 2 | yes |
| B7 | math | 173 | `K_req − K_req_paper == 0` | claim 5 | yes |
| B8 | math | 174 | `residual(k=K_req) == 0` | claim 4 | yes |

No tautological or orphaned assertions: each compares an *independently computed* quantity (integral, or solver root) against an *independently written* target form (closed-form literal or the paper's hand-written decomposition).

## Findings

None. (Stale-output observation noted below in Verdict justification and Self-test; it is non-material and does not rise to a filed finding.)

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration. The load-bearing constant `κ_n` is derived two genuinely different ways from the SymPy route:
- SymPy uses a single definite integral: `kappa_n = sp.integrate(u0 * f_n, (s, 0, L))` (py L97).
- Mathematica computes the *indefinite* integral and applies the fundamental theorem of calculus at the endpoints: `indef = Integrate[u0*fN, s]; kappaNViaFundamentalThm = (indef /. s->l) - (indef /. s->0)` (wl L62-66), and independently writes the analytic short form `Sqrt[2]/((n+1/2)*Pi)` (wl L74-77), then asserts the two agree (wl L82). The inline comment (wl L57-61) explicitly states this exercises a different code path. The analytic short form is anchored because it is checked against the FTC integration path, not merely printed.

The downstream algebra (`Δ,Q,P,B0…D0`, `K_req`) is definitional substitution of the same Stage-025 forms; identical structure across engines is expected there and is not a transliteration concern. Variable naming is fully independent (lowercase `lambdaB, omegaU, k, l, gConst, cSpeed, cs`).

## Engine cross-check

The two outputs agree at every reported quantity:
- `κ = 2√2/π`, numeric `0.900316316157106…` — sympy out L33-35, math out L36-39. Match.
- `Δ = Ω_U²Ω_W² − 8λ_R²/π²`, `Q`, `P`, `B0 = 8λ_B²/(π²ϖ²)`, `Z0`, `N0`, `D0` — sympy out L56-62, math out L62-68. Algebraically identical (Mathematica factors `Q` as `λ_U²Ω_W² + 8λ_W(2λ_Rλ_U+λ_WΩ_U²)/π²`, which expands to SymPy's form).
- `K_req` — sympy out L92-121 (pretty-printed) with assertion `K_req − K_req_paper = 0` (L122); math out L79 closed form `(20 a⁵ c⁵ mhat²(…)²π²)/(27 cs⁵ g(…)²) + Q/Δ + 8λ_B²/(π²ϖ²)`, assertion PASS (L81). The third term reduces correctly: `mhat²·(8/π²)·(…)² / [(54 g cs⁵/(5a⁵c⁵))·(π²Ω_U²Ω_W²−8λ_R²)²/π⁴] = 20 a⁵c⁵ mhat²(…)²π² / [27 cs⁵ g(π²Ω_U²Ω_W²−8λ_R²)²]`. Match.
- Both assert `residual @ K_req = 0` (sympy out L123, math out L83). Match.

`engines_agree: true`.

## Verdict justification

`clean`. I attacked: (1) the boxed `Q` and `P` against the substituted couplings (both expand correctly); (2) the normalization-law denominator `Δ(KΔ−Δκ²λ_B²/ϖ²−Q)` against the script's `Δ²·D0` (identical after factoring Δ²); (3) the three-term `K_req` third-term constant (`54/5 → 20/27` after multiplying through κ²=8/π² and Δ²; verified by hand, matches both engines); (4) the assumption set (`n` integer is required for `cos((n+½)π)=0` and is declared in both engines); (5) tautology (every assertion compares two independently-produced expressions). The target constant `54 G c_s⁵/(5 a⁵ c⁵)` and numeric κ live in the notes and reconcile (the `.tex` card legitimately carries the target symbolically). I confirmed I read the paper card, the notes, and the appendix row, and the scripts' verified claims match all six paper deliverables.

One non-material observation: the SymPy source mtime (2026-06-03) is newer than its saved output (2026-05-26), the literal `stale_output` signal. The cause is commit `e2a4780` (numbering reconciliation Phase 1), a doc-only one-token change ("Stage 9" → "Stage 26" in the docstring); `git diff` against HEAD is empty and no math/code logic changed, so the saved output still reflects the current script's computations. I do not file this as a blocking finding; the orchestrator's independent re-run will refresh the timestamp. (Residual stale "Stage-8"/"Stage 9" labels in the script comments/banners are the known pre-renumber numbering-drift artifacts owned by the separate numbering pass, not math or paper-value findings.)

## Self-test notes

Checked: (1) variable-independence — no `diff`/`D[]` in either script, and `solve(residual,K)` is well-posed since residual depends on K via `D0=K−B0−Z0`; trap N/A. (2) parity — all integrals are over the bounded `[0,L]`, so no symmetric-domain vanishing applies; the assertions check equality to *nonzero* closed forms, not vanishing, and `√2/[(n+½)π]≠0`. (3) trivial-case — at `n=0`, `κ_n` reduces to `2√2/π≈0.9003`, a nonzero literal matching the asserted value; the `K_req` back-substitution `residual(K=K_req)=0` is a genuine zero-residual confirmation, not a trivial pass. No traps tripped.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 10 deliverable values checked, 0 misaligned.

Every result value the scripts emit is reflected in the `.tex` card and/or the `.md` notes. The `.tex` card carries the closed-form symbolic results (κ, couplings, Δ/Q/P, normalization law, K_req); the notes additionally carry the two prose-only constants (numeric κ and the target `54 G c_s⁵/(5 a⁵ c⁵)`), which the terse `.tex` legitimately leaves symbolic. Output files are present; the SymPy `.txt` is timestamp-stale but content-current (see Verdict justification) — reconciliation is based on script source plus the committed outputs.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `κ_n = √2/[(n+½)π]` (=`2√2/(π(2n+1))`) | py L97-98, out L26-27; wl L74-77, out L28 | `.tex:40` eq:kappa-n; `.md:122` | MATCH |
| `κ = κ₀ = 2√2/π` | py L99, out L33; wl L85, out L36 | `.tex:46` boxed eq:kappa0; `.md:126` | MATCH |
| `κ` numeric `0.900316316157106` | py L110, out L35; wl L89, out L39 | `.md:130` (`≈0.900316316`) | MATCH |
| `I_ηφ=κ, I_ηu=1, I_ηw=κ, I_uw=κ` | py L132-135, out L42-45; wl L100-103, out L46-49 | `.tex:64-70` eq:overlaps; `.md:134-140` | MATCH |
| `C=κλ_B, G_U=λ_U, G_W=κλ_W, R=κλ_R` | py L148-151, out L52-55; wl L116-119, out L58-61 | `.tex:75-82` boxed eq:couplings; `.md:168-174` | MATCH |
| `Δ=Ω_U²Ω_W²−κ²λ_R²` (=`Ω_U²Ω_W²−8λ_R²/π²`) | py L153, out L56; wl L121, out L62 | `.tex:89` boxed eq:DeltaQP; `.md:178` | MATCH |
| `Q=λ_U²Ω_W²+2κ²λ_Uλ_Wλ_R+κ²λ_W²Ω_U²` | py L154, out L57; wl L122, out L63 | `.tex:90` boxed; `.md:180` | MATCH |
| `P=κ(Ω_U²λ_W+λ_Rλ_U)` | py L155, out L58; wl L123, out L64 | `.tex:91` boxed; `.md:182` | MATCH |
| `X/B0 = κ²λ_B²/ϖ²` (=`8λ_B²/(π²ϖ²)`) | py L157, out L59; wl L125, out L65 | `.tex:98` eq:X; `.md:186,190` | MATCH |
| `N_Q^target = 54 G c_s⁵/(5 a⁵ c⁵)` | py L188; wl L151 (target) | `.md:7,208,214` (`.tex` symbolic `N_Q^{\rm target}`) | MATCH (notes carrier) |
| `K_req` three-term decomposition | py L206-212, out L92-122; wl L166-173, out L79 | `.tex:115-119` boxed eq:Kreq; `.md:214-216` | MATCH |
| `K = K_η + 6T_Ω` | py L216, out L124; wl L176, out L84 | `.tex:124` eq:Kgeo; `.md:230` | MATCH |

INTERNAL (scaffolding, no prose expected; no finding): `Z0=Q/Δ`, `N0=P²/Δ²`, `D0=K−B0−Z0`, `P0=N0/D0`, the printed "Target residual", and `K_req_paper` (the in-script reconstruction used as the assertion target). `Z0/N0/D0` are intermediate combinations that drive the normalization assertion; `P0` is an unused convenience print; the residual is verification scaffolding.
