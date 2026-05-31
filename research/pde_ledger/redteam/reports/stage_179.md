---
unit_id: 179
batch: V.2
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-05-30T00:00:00-06:00
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage179_transfer_shape_theorem.md]
  paper_appendix: present
---

# Audit unit 179 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_179.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage179_transfer_shape_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (row at line 89; transfer-shape block at lines 671-715)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage179_transfer_shape_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage179_transfer_shape_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage179_transfer_shape_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage179_transfer_shape_mathematica_audit.txt`

## What the paper claims

Stage 179 is the "Transfer-shape theorem" (Part V, anchor MTDC-T9.2). `\stagefield{Output}` states verbatim: *"Factors each outgoing port as \(N_{A,0}^{(r)}=K_A\mathcal T_{A,r}^2\) and gives \(\Xi_1=2\sum_r\rho_r^{(N)}\tau_r\)."* The notes expand this into five distinct deliverables: (1) the exact wall-normalized factorization `N_0^(r)/K = T_r²` with `T_r = (Ĝ_W + R̂ Ĝ_U)/(1-R̂²)` from the primitive port data `P_r = Ω_U² G_W + R G_U`, `Δ_r = Ω_U²Ω_W² - R²` (§1); (2) the central slope identity `ν_r = κ₁ + 2τ_r` obtained by taking the weak-axisymmetric log-slope of `N_{A,0}^{(r)} = K_A T_{A,r}²` (§3); (3) the explicit transfer-shape slope `τ_r = α̂_r 𝔴_r + β̂_r(𝔲_r+𝔠_r) + (2R̂²/(1-R̂²))𝔠_r` with `α̂_r+β̂_r=1` and the wall-normalized port slopes `𝔴_r,𝔲_r,𝔠_r` defined in §2; (4) exact equivalence to the earlier slippage language, `τ_r = 𝔪_r + I_r/(1+I_r)𝔦_r + H_r/(1-H_r)𝔥_r` (§5); and (5) the grouped-defect collapse `Ξ₁ = 2Σ_r ρ_r^(N) τ_r` given `Σ_r ρ_r^(N)=1`, starting from the Stage-246 relation `Ξ₁ = Σ_r ρ_r^(N)(ν_r-κ₁)` (§4). The appendix row (line 89) restates the same Output; the appendix block (lines 671-686) restates deliverables 1 and 5 identically.

## What the script claims to verify

The SymPy docstring lists four checks that map onto the paper deliverables: the factorization `N0/K = T²`, the slope identity `ν_r = κ₁ + 2τ_r`, the slippage-form equivalence, and the weighted defect identity. The assertions (via `expect_zero`, which simplifies+expands and raises on nonzero) verify exactly these: `N0/K - T² == 0` (line 53), `nu_direct - (kappa1 + 2 tau) == 0` (line 91), `tau - slippage form == 0` (line 102), `(nu-kappa1) - 2*tau_slippage == 0` (line 103), and `Xi_1 - 2 weighted tau == 0` (line 112). Crucially, `nu_direct` (line 73) is computed *independently* by perturbing the primitive port data with first-order slopes and extracting the eps-coefficient of `log(N0A)`, then compared against the transfer-shape-built `nu_expected`; the comparison is non-tautological. The Mathematica `.wl` mirrors the same five checks with `expectZero` and `Exit[1]` on failure.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| (1) `N_0^(r)/K = T_r²`, `T_r=(Ĝ_W+R̂Ĝ_U)/(1-R̂²)` from primitive `P,Δ` | sympy L53 / math L44 `N0/K - T^2 == 0` | match |
| (2) `ν_r = κ₁ + 2τ_r` (log-slope of `N_0=K T²`) | sympy L91 / math L75 `nu_direct - (kappa1+2 tau) == 0` | match |
| (3) `τ_r = α̂𝔴 + β̂(𝔲+𝔠) + (2R̂²/(1-R̂²))𝔠`, `α̂+β̂=1` | folded into `tau` build (sympy L78-84 / math L64-70), exercised by check (2) | match |
| (4) `τ_r = 𝔪 + I/(1+I)𝔦 + H/(1-H)𝔥` | sympy L102 / math L84 `tau - slippage form == 0` (+L103/L85) | match |
| (5) `Ξ₁ = 2Σρ_r τ_r`, `Σρ=1`, from `Ξ₁=Σρ(ν-κ₁)` | sympy L112 / math L93 `Xi_1 - 2 weighted tau == 0` | match |
| (appendix) `T_eff² := Σ_r T_{A,r}² = N_{A,0}/K_A` (eq:app-part05-Teff-def) | (none) | not a stage-179 card/notes deliverable; jointly attributed to "Stages 179--180"; not required here |

`paper_alignment: aligned`. Every variable definition (`Ĝ_W, Ĝ_U, R̂, 𝔴_r, 𝔲_r, 𝔠_r, α̂_r, β̂_r, I_r, H_r, 𝔪_r, 𝔦_r, 𝔥_r`) in both scripts reproduces the corresponding boxed notes definition exactly (verified term-by-term).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 53 | `simplify(N0/K - T^2) == 0` | claim 1 | yes |
| A2 | sympy | 91 | `nu_direct - (kappa1 + 2 tau) == 0` | claims 2 + 3 | yes |
| A3 | sympy | 102 | `tau - tau_slippage == 0` | claim 4 | yes |
| A4 | sympy | 103 | `(nu_expected-kappa1) - 2*tau_slippage == 0` | claims 2+4 | yes |
| A5 | sympy | 112 | `Xi - 2*sum(rho*tau) == 0` | claim 5 | yes (low-substance; see note) |
| B1 | math | 44 | `expectZero[N0/K - T^2]` | claim 1 | yes |
| B2 | math | 75 | `expectZero[nu_direct - (kappa1+2 tau)]` | claims 2 + 3 | yes |
| B3 | math | 84 | `expectZero[tau - tauSlippage]` | claim 4 | yes |
| B4 | math | 85 | `expectZero[(nu-kappa1) - 2*tauSlippage]` | claims 2+4 | yes |
| B5 | math | 93 | `expectZero[xi - 2*sum(rho*tau)]` | claim 5 | yes (low-substance) |

Note on A5/B5: the weighted-defect check uses generic symbolic `tau1,tau2,tau3` and `rho3 = 1-rho1-rho2`, so the identity collapses to `κ₁·1 + 2Στ - κ₁ = 2Στ` — algebraically light, following directly from `Σρ=1`. It is *not* tautological in the prohibited "x=expr then x==expr" sense: it independently constructs `Ξ₁` from the per-port relation `ν_r=κ₁+2τ_r` and confirms the Stage-246 collapse `Σρ(ν-κ₁) → 2Στ`. Deliverable 5 is itself just this algebraic collapse, so a light check is a faithful (not deficient) test. Not a finding.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a genuine independent re-derivation, not a transliteration, on the load-bearing step. Both engines define the same physical variables (necessarily — they are the paper's defined quantities), but the central slope `nu_direct` is obtained by *different mechanisms*: SymPy uses `sp.series(sp.log(N0A), eps, 0, 2).removeO().coeff(eps, 1)/lam` (line 73 — Taylor-series coefficient extraction), whereas Mathematica uses `(D[Log[n0A], eps] /. eps -> 0)/lam` (line 61 — analytic derivative evaluated at eps=0). These are mathematically equivalent extractions implemented via distinct engine machinery; a transliteration would have ported the same `series`/`coeff` call. The simplification back-ends also differ (`sp.simplify(sp.expand(...))` vs. `FullSimplify[Together[Expand[...]]]` under explicit `$Assumptions`). The remaining variable definitions are dictated by the paper and are expected to coincide; the verification mechanism is independent. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines produce identical final residuals on all five checks:

| Check | SymPy output | Mathematica output |
|---|---|---|
| `N0/K - T^2` | `0` (L16) | `0` / `PASS` (L16-17) |
| `nu_direct - (kappa1 + 2 tau)` | `0` (L73) | `0` / `PASS` (L25-26) |
| `tau - slippage form` | `0` (L78) | `0` / `PASS` (L31-32) |
| `(nu-kappa1) - 2*tau_slippage` | `0` (L79) | `0` / `PASS` (L33-34) |
| `Xi_1 - 2 weighted tau` | `0` (L84) | `0` / `PASS` (L39-40) |

The printed `nu_direct` and `tau` expressions also agree symbolically (SymPy pretty-print L21-72 vs. Mathematica InputForm L22-24): both give `tau = gw·ou2(gW-κ₁/2-oW)/(gw·ou2+gu·r) - gu·r(-2gU+κ₁+2(oU+oW-rr))/(2(gw·ou2+gu·r)) + r²(oU+oW-2rr)/(r²-ou2·ow2)`. Engines agree. Both exit 0.

Freshness: sympy `.py` mtime 2026-05-11 11:58:53; sympy `.txt` 12:47:58 (newer). math `.wl` 11:58:53; math `.txt` 13:22:49 (newer). Outputs are fresh (`outputs_fresh: true`). No `stale_output` finding.

## Verdict justification

`clean`. I read the paper card, the §1-§8 notes, and the Part V appendix block, built the five-deliverable model, then attacked the scripts. Attacks attempted and failed to break anything: (a) tried to make the slope identity tautological — it is not, because `nu_direct` is built from perturbed primitive data and compared to the transfer-shape form, with the two sides genuinely distinct; (b) checked the symbol domains — `R` is `real` (more general than `R>0`, correct since `R` can be negative), and no sqrt of `Delta` is ever taken so there is no hidden branch ambiguity, while `1-R̂² = Delta/(Ω_U²Ω_W²)` keeps the factorization exact for any sign of `R`; (c) verified the `N0/K=T²` identity is a pure rational identity unaffected by sign (squared throughout); (d) checked every variable definition term-by-term against the boxed notes — all match; (e) checked deliverable coverage — all five card/notes deliverables have non-tautological checks in both engines, and the only paper expression without a script check (appendix `T_eff²`) is a Stages-179--180 aggregation not listed as a stage-179 deliverable. The script's verified claim matches the paper's stated claim exactly.

One non-blocking observation (not a finding under any of the ten categories, since it gates no assertion and the verified math is correct): the script banners read `"STAGE 162"` (sympy L32 / math L26), the section headers and docstring reference `"Stage 176/160/161"` (sympy L10,L93 / math L77), and these strings propagate into the saved transcripts (sympy `.txt` L11; math `.txt` L11). The notes themselves carry yet another old numbering (`Stage 246/247/244/228/229/245`). All of these are renumbering artifacts left in display strings; the actual physics, variable definitions, and verified identities are unambiguously the Stage-179 transfer-shape theorem and match the current paper. No math/verification defect — purely cosmetic label drift, recorded here for orchestrator awareness, not raised as a numbered finding because it fits none of the ten finding categories and weakens no check.

## Self-test notes

(1) Variable-independence: the only derivative is the eps-slope of `log(N0A)`; `N0A` genuinely depends on `eps` through `KA,OU2A,OW2A,GWA,GUA,RA` (lines 61-66), so the first-order coefficient is non-trivially nonzero (confirmed by the printed multi-term `nu_direct`) — no identically-zero-derivative trap. (2) Symmetry/parity: no integrals in this unit; N/A. (3) Trivial-case: the factorization `N0/K=T²` and the slope identity are exact rational identities verified symbolically in both engines (residual 0), and the weighted identity reduces correctly under `Σρ=1`. (4)/(5): N/A and no fix prescribed — verdict is clean, so no directive is written and no new `paper_misalignment` can be introduced.
