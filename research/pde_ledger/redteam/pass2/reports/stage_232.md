---
unit_id: 232
batch: VII.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 232 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_232.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part07.tex` (rows at L76, L867–881; \input at L1453)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage232_known_5pn_data_injection_and_current_branch_verdict_mathematica_audit.txt`

## What the paper claims

Stage 232 injects the numerically-located Family-1 5PN support/source data into the selected-branch same-charge inequalities. `\stagefield{Output}`: "Current branch verdict: the support/source side is safe by a large margin in the injected data; the active unresolved gate is static orbit-lock/coherent placement rather than support starvation." The notes (the authoritative source, much richer than the terse card) enumerate the deliverables: the refreshed geometry (Λ_ℓ=20√2π/x₀₁≈36.9497, κ=(9/5)Λ_ℓ²≈2457.51), the Robin support ceiling (y·tan y=η ⇒ A_K≈1.0000521, ζ_max=A_K·π²/4≈2.46753), the support-drop endpoints (Δ_0≈1.7377e-4, Δ_∞≈2.0172e-2), the fixed-point law Pe=Ξ·Δ(Pe), the two figures of merit Ξ_χ=100·Θ_w^(χ)·Λ_ℓ²≈5.5548e5 and Ξ_J=100·Θ_w^(J)·Λ_ℓ²≈1.2664e5, the two roots Pe_*^(χ)≈11155.73, Pe_*^(J)≈2504.97, the physical ratios ζ_phys/ρ_α,max for both branches (≈2.46753/3.46753), the margins against ζ_req=1/3, ρ_req=4/3, and the near-saturation ceiling gaps (≈9.78e-8, ≈1.94e-6). The appendix (L867–874) carries ζ_req=1/3, ρ_req=4/3, ζ_phys≈2.4675 and the static-first verdict.

## What the script claims to verify

Both engines reconstruct the full 5PN injection chain: (M1) refreshed geometry from the J_0 zero; (M2) the lowest Robin root and the A_K/ζ_max ceiling; (M3) the support-drop kernel and its Δ_0/Δ_∞ endpoints; (M4) the 100-prefactor figures of merit; (M5) the two fixed-point roots from Pe=Ξ·Δ(Pe) with bracket + residual checks; (M6) the physical support ratios ζ_phys, ρ_α,max via the overlap-boost law Ω_Pe; (M7) the injected margins, ratios, and ceiling gaps versus ζ_req=1/3, ρ_req=4/3. Each computed quantity is checked against its pinned decimal AND, where applicable, against an independent route (closed-form vs quadrature in .py; native Integrate + Limit vs uniform-source integral in .wl). Positivity asserts back the bottom-line "safe by a large margin / non-bottlenecked" verdict.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Λ_ℓ, κ refreshed geometry | M1 (py L49–64 / wl L48–60) | match |
| Robin ceiling A_K, ζ_max | M2 (py L69–83 / wl L62–80) | match |
| Δ_0, Δ_∞ endpoints | M3 (py L106–141 / wl L82–116) | match |
| Ξ_χ, Ξ_J figures of merit (100·…) | M4 (py L146–157 / wl L118–128) | match |
| Pe_*^(χ), Pe_*^(J) roots | M5 (py L162–187 / wl L130–153) | match |
| ζ_phys, ρ_α,max ratios | M6 (py L192–212 / wl L155–167) | match |
| margins vs ζ_req=1/3, ρ_req=4/3; ratios; ceiling gaps | M7 (py L217–258 / wl L169–204) | match |
| Output verdict "safe by large margin; gate=static placement" | positivity asserts + printed verdict (py L249–267) | match |

`paper_alignment: aligned` — every notes/appendix deliverable has a faithful, non-tautological script-side check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 62–64 | assert_close Λ_ℓ, χ_s, κ | geometry | yes |
| A2 | sympy | 81–83 | assert_close y, A_K, ζ_max | Robin ceiling | yes |
| A3 | sympy | 113–114,134–139 | assert_close Δ_0/Δ_∞ + quad-vs-closed + endpoint bounds | endpoints | yes |
| A4 | sympy | 156–157 | assert_close Ξ_χ, Ξ_J | figures of merit | yes |
| A5 | sympy | 181–187 | bracket + residual + decimal on Pe_* | roots | yes |
| A6 | sympy | 209–212 | assert_close ζ_phys, ρ_α,max | ratios | yes |
| A7 | sympy | 249–258 | positivity + decimals on margins/ratios/gaps | output verdict | yes |
| M1 | math | 59–60 | expectClose Λ_ℓ, κ | geometry | yes |
| M2 | math | 76–80 | expectTrue range, expectSmall y·tan y−η, decimals | Robin ceiling | yes |
| M3 | math | 97–116 | unevaluated-Integrate guard + formula-vs-Limit + decimals | endpoints | yes |
| M4 | math | 127–128 | expectClose Ξ_χ, Ξ_J | figures of merit | yes |
| M5 | math | 144–153 | bracket + residual + decimals on Pe_* | roots | yes |
| M6 | math | 166–167 | expectClose ζ_phys | ratios | yes |
| M7 | math | 189–204 | expectPositive + decimals on margins/ratios/gaps | output verdict | yes |

No tautological or orphaned rows: every assertion traces to a paper deliverable, and the decimal targets are matched against independently *computed* quantities (not literals re-asserted against themselves).

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine second engine, not a transliteration:
- **M3 (the load-bearing step)**: `.py` hand-derives a closed form `Delta_closed` (py L117–128) and validates it with `mp.quad` numeric quadrature. `.wl` instead performs **native symbolic** `Integrate[kernel·source, {x,0,1}]` (wl L91–96), guards that it actually evaluated (`If[!FreeQ[...,Integrate], fail]`, wl L97–99), and then cross-validates the Δ_0/Δ_∞ formulas via `Limit[...,pVar->0]`, `Limit[...,Infinity]`, AND a uniform-source integral (wl L105–107, residuals 0 at ~1e-56 to ~1e-80). Two structurally different derivation routes (analytic closed-form + numeric quadrature vs symbolic integration + symbolic limits).
- **M5 root-finding**: `.py` uses a custom bisection on `[Ξ·Δ_0, Ξ·Δ_∞]` (py L162–169); `.wl` uses `FindRoot` from a numeric seed (wl L134–141). Different solvers; both verify the fixed-point residual `Pe−Ξ·Δ(Pe)≈0`.
- Shared pinned decimal targets are reference constants, not algebra echoes; cross-engine agreement at ~1e-13 confirms independent reproduction.

Not a `mathematica_transliteration` finding.

## Engine cross-check

Final values agree to high precision (py output vs wl output): Λ_ℓ, κ, A_K, ζ_max, Δ_0, Δ_∞, Ξ_χ=555483.32…, Ξ_J=126637.07…, Pe_*^(χ)=11155.7265863…, Pe_*^(J)=2504.9703142…, ζ_phys^(χ)=2.46752964788…, ζ_phys^(J)=2.46752780517…, margins 2.13419631…/2.13419447…, ceiling gaps 9.78e-8/1.94e-6. All `expectClose`/`assert_close` pass; both outputs end with "All Stage 232 checks passed." No `engine_disagreement`.

## Verdict justification

Clean. Attacks tried that failed: (1) checked for the variable-independence self-test trap — no `sp.diff`/`D[]` derivatives are used as load-bearing checks here, so no vacuous-derivative pass. (2) Checked the hardcoded `Theta_w_chi`/`Theta_w_J` literals (py L146–147, wl L119–120): these are explicitly **imported** carried 5PN/Family-1 wall-depth inputs per notes §1.3 ("two explicit wall-depth extractions … carried Family-1 branch"), not stage-232 derivations, so a literal is appropriate — and they are then fed through Ξ→Pe_*→ζ_phys, not confirmed against themselves. (3) The two seeded `Pe_*` roots are independently re-solved (bisection / FindRoot) and the fixed-point residual is asserted, so they are not tautological. (4) The 100-prefactor (py L149/L154, wl L121) matches the notes' L153/L157 `100·…` and reproduces the notes' quoted Ξ_χ/Ξ_J decimals — the recurring stale "168" family is fully cleared (see below). (5) Positivity asserts genuinely back the "safe by large margin" output. Both engines present, independent, fresh outputs (txt mtime Jun 2 22:01, after both script mtimes). Paper/notes/appendix read and aligned.

## Value Reconciliation (pass-2 augmentation)

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Λ_ℓ ≈ 36.94973154240256 | py L62 / out L5 | notes L44–48 | MATCH |
| χ_s ≈ 18.47486577120128 | py L63 / out L6 | notes L56 | MATCH |
| κ ≈ 2457.508789900114 | py L64 / out L8 | notes L60–64 | MATCH |
| A_K ≈ 1.0000521380385143 | py L82 / out L12 | notes L74–78 | MATCH |
| ζ_max ≈ 2.4675297457259358 | py L83 / out L13 | notes L82–87 (boxed) | MATCH |
| Δ_0 ≈ 1.7377393923469950e-4 | py L113 / out L16 | notes L118–124 | MATCH |
| Δ_∞ ≈ 2.0172162594593645e-2 | py L114 / out L17 | notes L126–131 | MATCH |
| Ξ_χ ≈ 5.5548332017764099e5 | py L156 / out L21 | notes L152–155 (100·) | MATCH |
| Ξ_J ≈ 1.2663707072528143e5 | py L157 / out L22 | notes L156–159 (100·) | MATCH |
| Θ_w^(χ) = 4.06863235008162 | py L146 / wl L119 | notes L144–149 | MATCH (imported) |
| Θ_w^(J) = 0.927552032539308 | py L147 / wl L120 | notes L144–149 | MATCH (imported) |
| Pe_*^(χ) ≈ 11155.7265863205869 | py L186 / out L25 | notes L166 | MATCH |
| Pe_*^(J) ≈ 2504.9703142859238 | py L187 / out L26 | notes L182 | MATCH |
| ζ_phys^(χ) ≈ 2.4675296478814376 | py L209 / out L31 | notes L168–172 | MATCH |
| ρ_α,max^(χ) ≈ 3.4675296478814376 | py L210 / out L32 | notes L174–178 | MATCH |
| ζ_phys^(J) ≈ 2.4675278051675084 | py L211 / out L33 | notes L184–189 | MATCH |
| ρ_α,max^(J) ≈ 3.4675278051675084 | py L212 / out L34 | notes L190–194 | MATCH |
| ζ_req = 1/3, ρ_req = 4/3 | py L217–218 / out | notes L220–224; appendix L870–872 | MATCH |
| margin_ζ^(χ) ≈ 2.1341963145481043 | py L253 / out L37 | notes L230 | MATCH |
| margin_ζ^(J) ≈ 2.1341944718341751 | py L254 / out L39 | notes L240 | MATCH |
| ratio_ζ^(χ) ≈ 7.402588943644313 | py L255 / out L43 | notes L252 | MATCH |
| ratio_ζ^(J) ≈ 7.402583415502525 | py L256 / out L45 | notes L261 | MATCH |
| ratio_ρ^(χ) ≈ 2.600647235911078 | py L257 / out L44 | notes L255 | MATCH |
| ratio_ρ^(J) ≈ 2.600645853875631 | py L258 / out L46 | notes L264 | MATCH |
| ceiling gap χ ≈ 9.784449817674381e-8 | wl L199 / out L49 | notes L274 | MATCH |
| ceiling gap J ≈ 1.9405584273504838e-6 | wl L200 / out L50 | notes L280 | MATCH |

Internal scaffolding (no finding): x01 (J_0 zero, a known transcendental input), y Robin root, brackets lo/hi, residuals, the `Delta_closed`/`deltaIntegralExpr` intermediates, tolerances, Ω_Pe intermediate.

reconciliation: complete; 27 deliverable values checked, 0 misaligned.

**Surviving `168` check:** NONE found as a stale prefactor. Grep of `.py`, `.wl`, card, notes, and the .py output returned zero hits. The only `168` anywhere is an incidental digit-substring inside a residual mantissa in the `.wl` output (`...3821686…`, line 97) — not a quantity, not a prefactor. The recurring stale-168 family is fully cleared at stage 232; notes L153/L157 correctly carry the `100` prefactor and the scripts emit `100` (py L149/L154, wl L121). Notes title correctly reads "Stage 232" (no surviving "Stage 249" drift). Card states "Mathematica audit: none yet" (stage_232.tex L11) while a `.wl` now exists — a stale verification-note lag, not a value misalignment; consistent with the deferred numbering/card-text-lag policy, noted not flagged.

## Self-test notes

Checked the variable-independence trap: no load-bearing `sp.diff`/`D[]` derivative checks exist here, so no vacuous-zero pass. Checked symmetry/parity of the M3 integral: kernel·source on [0,1] is a one-sided integral over a finite domain (not symmetric/unbounded), and both engines confirm Δ_0<Δ(Pe)<Δ_∞, so the nonzero claim is real. Trivial-case pre-check: positivity asserts (margins, gaps) are confirmed by the printed positive decimals; the fixed-point residuals reduce to ~0 (~1e-30/1e-57) at the solved roots. No directive written (zero findings).
