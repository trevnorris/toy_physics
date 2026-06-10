---
unit_id: 248
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-10T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md]
  paper_appendix: present
---

# Audit unit 248 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_248.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md`
- part appendix: cited a wrong project root (mis-rooted) → actual `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.txt`

## What the paper claims

Stage 248 (`\stagefield{Output}`-equivalent is the §7 "What Stage 248 proves" list plus the boxed card equations) promotes the Stage-247 stationary lowered barrier `V := V_eff^(247)` into a dynamic reduced 1D scattering / event-chain compiler. The card deliverables are: (1) exact reduced energy conservation `m_s ddot r = -V'(r)`, `E = (m_s/2) dot r^2 + V`; (2) the finite-radius classical threshold law `v_crit,new = sqrt(2(V_peak - V(r0))/m_s)` and the pure-Coulomb contact threshold `v_contact,Coul = sqrt(2(1/r_contact - 1/r0)/m_s)`, with the window `v_crit,new < v0 < v_contact,Coul`; (3) subbarrier turning points with transport law `dr_±/dE = 1/V'(r_±)`; (4) the WKB action `I_new`, transmission ratio `T_new/T_Coul = exp[-2(I_new - I_Coul)]`, the exact Coulomb outer turning point `r_turn,Coul = 1/E` and the closed-form `I_Coul` (with `arcsin sqrt(E r_contact)` for `E r_contact < 1`); (5) the near-top parabolic normal form `I_top = π ΔE sqrt(m_s/K_peak)/ħ_eff`; (6) the carried diagnostics `Ξ_turn = Ξ1(r_+)`, `λ_th = |E/V'(r_+)|`; and (7) the Session-II benchmark numbers, including the boxed transmission increase `(T_new/T_Coul - 1)×100% ≈ 23.3128%` (notes §6.3, card via transmission-ratio eq).

## What the script claims to verify

Both scripts verify, in five sections: (1) on-shell `dE/dt = 0`; (2) the two threshold-speed closed forms satisfy their defining energy equalities and have positive branches; (3) the Coulomb antiderivative identity `F_coul' = integrand`, the endpoint evaluation `= I_Coul` closed form, both turning-point transport laws, the `λ_th` identity (with an explicit non-vacuity guard), and the `Ξ_turn` carried tag; (4) the near-top parabolic action; (5) the Session-II benchmark numerics — `v_crit,new`, `v_contact,Coul`, `v_0,sub`, `r_turn,Coul`, `T_new/T_Coul`, the `23.3128%` improvement, the cross-speed window ordering, and the `I_Coul` exact-vs-report agreement.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| Energy conservation `dE/dt=0` | py L58 / wl L91 | match |
| `v_crit,new` closed form | py L104 (solve==form), L103 (E@vcrit=Vpeak) / wl L121,L122 (satisfaction) | match |
| `v_contact,Coul` closed form | py L106,L107 / wl L123,L124 | match |
| window `v_crit < v0 < v_contact` | py L308 / wl L272 (Sign product) | match |
| transport `dr_±/dE = 1/V'` | py L178,L179 / wl L193,L194 | match |
| `I_new` action def | py L119 (printed) / wl L? printed | match (symbolic def) |
| `r_turn,Coul = 1/E` | py L125 / wl L141 | match |
| `I_Coul` closed form + antideriv | py L176,L177 / wl L191,L192 | match |
| `T_new/T_Coul = exp[-2(I_new-I_Coul)]` | py L162 / wl L189 (symbolic) | match |
| `I_top` near-top form | py L227 / wl L225 | match |
| `Ξ_turn = Ξ1(r_+)`, `λ_th = E/V'` | py L198-200 / wl L204-206 | match |
| benchmark numerics incl. 23.3128% | py L300-312 / wl L262-274 | match |

`paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `dE_on_shell == 0` | energy conservation | yes |
| A2 | sympy | 103-107 | `E_at_vcrit==0`, solve==form (×2), delta vanishes (×2) | threshold laws | yes |
| A3 | sympy | 176 | `dF_coul == 0` | Coulomb antiderivative | yes |
| A4 | sympy | 177 | `I_coul_endpoints - I_coul_formula == 0` | I_Coul closed form | yes |
| A5 | sympy | 178-179 | `transport_± - 1/V'(r_±) == 0` | transport law | yes |
| A6 | sympy | 198-199 | `λ_th gap raw != 0`, `gap(V=E)==0` | λ_th identity | yes (non-vacuous) |
| A7 | sympy | 227 | `I_top - I_top_expected == 0` | near-top form | yes |
| A8 | sympy | 300-312 | benchmark `abs(...) < tol` ×N | Session-II numerics | yes |
| B1 | math | 91 | `expectZero[dE/dt]` | energy conservation | yes |
| B2 | math | 109-124 | non-vacuity guards + satisfaction `EAtVcrit-Vpeak==0` etc. | threshold laws | yes (satisfaction) |
| B3 | math | 191 | `expectZero[dFcoul]` | Coulomb antiderivative | yes |
| B4 | math | 192 | `IcoulEndpoints - IcoulFormula == 0` | I_Coul closed form | yes |
| B5 | math | 193-194 | transport laws | transport | yes |
| B6 | math | 204-205 | λ_th non-vacuous guard + identity | λ_th identity | yes |
| B7 | math | 225 | near-top normal form | near-top | yes |
| B8 | math | 262-274 | `expectNear` benchmark ×N + Sign flags | Session-II numerics | yes |

No orphaned assertions; every script check traces to a card/notes deliverable.

## Findings

None. Zero findings.

## Independent-derivation check (Mathematica)

The `.wl` is **not** a transliteration; it implements the reframed **satisfaction route** described in the checkpoint heads-up, giving genuinely different information flows from the `.py`:

1. **Threshold speeds (load-bearing).** `.py` L70-76 DERIVES the speed: `sp.solve(sp.Eq(E_launch_new, Vpeak), v0)` then asserts the closed form equals the solve result. `.wl` L97-122 does the opposite — it POSITS `vcritNew = Sqrt[2(Vpeak-V0)/ms]` and verifies it SATISFIES the defining energy equality `EAtVcrit - Vpeak == 0` (back-substitution), gated by a non-vacuity guard (L109-112: fails if `deltaNew` is identically zero OR free of `v0`) and a positive-branch guard (L126-130). One DERIVES what the other POSITS. Non-vacuous: if the posited closed form were wrong, `deltaNew /. v0 -> vcritNew` (L122) would not vanish and `expectZero` would `Exit[1]`; the output L15-21 confirms the gap printed as `(ms*v0^2)/2 + V0 - Vpeak` (genuinely v0-dependent) before collapsing to 0 only after the correct substitution.

2. **Coulomb antiderivative.** `.py` L128-136 builds `F_coul` and checks `diff(F_coul, rv) - coul_integrand == 0`. `.wl` L142-158 does the same identity but through a different Wolfram path: `PowerExpand` → `Refine` under `Esub rv < 1` → an explicit `Sqrt[rv^-1] Sqrt[rv] -> 1` rewrite → `FullSimplify`. The branch-handling choreography differs materially (SymPy needs no PowerExpand/Refine ladder).

3. **λ_th non-vacuity guard.** Both engines guard the identity the same way conceptually (`raw != 0`, then `(V→E) == 0`), but implemented in each CAS's native idiom; outputs L33-34 (py) and L45-48 (wl) both show the raw gap is the non-trivial `(-E + V(r_+))/V'(r_+)` collapsing to 0 only under the turning condition — confirming the identity is not free.

Conclusion: independent re-derivation confirmed; the satisfaction route is present, non-vacuous, and load-bearing. No `mathematica_transliteration` finding.

## Engine cross-check

Both outputs agree at the level claimed. Symbolic residuals are 0 in both (`dE/dt`, antiderivative, endpoints, transport, λ_th, near-top). Benchmark: `.py` prints `improvement (%) = 23.31275...`, `T_new/T_Coul = 1.23312754`, `I_Coul exact = 0.30230579644795263`; `.wl` `expectNear` confirms `improvement` diff `4.57e-5` (tol 1e-3), `ratio` diff `1.7e-8` (tol 5e-8), `I_coul exact vs report` diff `8.28e-5` (tol 1e-3). The two engines' `I_Coul` symbolic forms differ cosmetically (`.py`: `... + π/2` with `asin`; `.wl`: `ArcCos[Sqrt[Esub rContact]]` — note `arccos x = π/2 - arcsin x`, so they are equal), which is the expected sign of independent simplification rather than a port. No `engine_disagreement`.

## Checkpoint-specific verifications (per heads-up)

- **Surviving 168:** `grep` across `.py`, `.wl`, card, notes, both `.txt` outputs returns **NO hits** for `168`. The recurring stale-168 family is fully cleared here. Notes line 506-507 read `\times 100\%` / `23.3128\%` (NOT `168%`); scripts emit `improvement (%) = 23.31...` from `improve_pct = 100.0 * (ratio_num - 1.0)` (py L259) / `improvePct = 100.0 (ratioNum - 1.0)` (wl L253). Confirmed `×100%`, notes match.
- **Banner canonical:** `.wl` L65 banner reads `STAGE 248 — DYNAMIC EVENT-CHAIN COMPILER ...`; output L3 echoes `STAGE 248`. The pass-1 `231→248` banner fix is in place; canonical.
- **Output freshness:** both `.txt` (mtime 2026-06-03 12:33) are newer than `.py` (10:12) and `.wl` (10:21). Not stale.

## Verdict justification

Clean. Attacks tried that failed: (a) tautology hunt on the threshold checks — the `.wl` route is genuinely back-substitution with a real non-vacuity guard that would catch a wrong closed form, and the `.py` route independently solves, so neither is `x = expr; assert x == expr`; (b) variable-independence trap on the transport-law `Solve` and the `λ_th` gap — the raw gap `(-E+V(r_+))/V'(r_+)` provably depends on the turning condition and is nonzero pre-substitution (output confirms), so the `!= 0` guard is real and the `== 0` is not vacuous; (c) Coulomb antiderivative — checked the `arcsin`/`arccos` complementarity and the `E r_contact < 1` branch domain, both consistent across engines; (d) symbol-assumption check — all positivity assumptions (`m_s>0`, `Esub>0`, `r_contact>0`, etc.) are justified by the physical reduced-scattering setup and match between engines; (e) exhaustive 168 sweep — none survive. Paper card, notes, and appendix part08 rows (L94, L250-277) all match the script claims. Both engines present, substantive, independent. Checkpoint bar MET.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 18 deliverable values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `v_crit,new = 2.54139063...` | py out L51 / wl expectNear L262 | notes L450; card §6.1 / appendix eq main-vcrit-new | MATCH |
| `v_contact,Coul = 3.27278339...` | py out L52 / wl L263 | notes L459 | MATCH |
| `v_0,sub = 2.14476202` | py out L53 / wl L264 | notes L475 | MATCH |
| `r_turn,Coul = 0.4` | py out L56 / wl L265 | notes L479 | MATCH |
| `r_peak = 0.23944389` | py out L48 | notes L425 | MATCH |
| `V_peak = 3.42933112` | py out L49 | notes L425 | MATCH |
| `V(r0) = 0.19999794` | py out L50 | notes L426 | MATCH |
| `r_turn,new = 0.39096144` | py out L54 | notes L428 | MATCH |
| `r_inner = 0.19039548` | py out L55 | notes L430 | MATCH |
| `I_new = 0.19744614` | py out L57 | notes L432 | MATCH |
| `I_Coul(report) = 0.30222297` | py out L58 | notes L487 | MATCH |
| `I_Coul(exact) = 0.30230580` | py out L59 / wl L271 | notes L553 | MATCH |
| `T_new = 0.67375262` | py out L60 | notes L491 | MATCH |
| `T_Coul = 0.54637707` | py out L61 | notes L494 | MATCH |
| `T_new/T_Coul = 1.23312756` | py out L62 / wl L268 | notes L501 | MATCH |
| `improvement = 23.3128%` | py out L63 / wl L269 | notes L507; card transmission-ratio | MATCH |
| `Ξ_turn = 0.34437471` | py out L64 | notes L436; card diagnostics eq | MATCH |
| `λ_th = 0.42826825` | py out L65 | notes L438; card diagnostics eq | MATCH |
| `v_cross = 2.59221845` | py out L67 | notes L515 | MATCH |
| `r_turn,Coul@cross = 0.28091...` | py out L68 / wl L270 | notes L534/536 (0.28091705) | MATCH |

INTERNAL (scaffolding, no finding): symbolic residuals (`dF_coul`, `lambda_th_gap`, `I_top - I_top_expected`, transport residuals), `E_launch_new/coul` intermediate forms, `|V'(r_turn)| = 5.8374...` (a back-derived diagnostic check), pass/fail flags, tolerances.

All emitted deliverable values reconcile against the notes (and, where the card carries them, the card/appendix). No MISMATCH, no MISSING-DELIVERABLE.

## Self-test notes

Checked: (1) variable-independence — the transport-law `Solve` differentiates `V(r_±(E)) - E` w.r.t. `E` where `r_±` is an explicit `Function` of `E`, so the chain-rule derivative is genuinely nonzero; the `λ_th` raw gap depends on `V(r_+)≠E` and is nonzero pre-substitution (no zero-derivative trap). (2) Parity/symmetry — the only symmetric-domain integral is `I_top` over `[-yturn, yturn]` of an even integrand `sqrt(2m_s(ΔE - K y²/2))`, correctly nonzero. (3) Trivial-case — the `.wl` satisfaction route would fail (`Exit[1]`) for any wrong closed form, confirmed by the printed non-vacuous gaps. No directive written (zero findings).
