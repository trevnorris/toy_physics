---
unit_id: 248
batch: VIII.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-03T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md"]
  paper_appendix: present
---

# Audit unit 248 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_248.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part08.tex` (rows 54, 94, 250-277 reference this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.txt`

## What the paper claims

Stage 248 promotes the Stage-247 lowered stationary barrier `V := V_eff^(247)` into a reduced one-dimensional dynamic event chain `m_s r-ddot = -V'(r)` with conserved energy `E = (m_s/2) r-dot^2 + V`. The card enumerates these deliverables: (1) energy conservation on the fixed front end; (2) the finite-radius classical threshold speed `v_crit,new = sqrt(2(V_peak - V(r0))/m_s)` and the pure-Coulomb contact threshold `v_contact,Coul = sqrt(2(1/r_contact - 1/r0)/m_s)`, with the classical window `v_crit,new < v0 < v_contact,Coul`; (3) the subbarrier turning-point equations `V(r_±(E)) = E` and the exact transport law `dr_±/dE = 1/V'(r_±(E))`; (4) the WKB action `I_new(E) = (1/ℏ_eff)∫ sqrt(2m_s(V-E)) dr`, the exact Coulomb turning point `r_turn,Coul = 1/E`, the closed-form Coulomb action `I_Coul(E)` (eq stage248-coul-action), and the transmission ratio `T_new/T_Coul = exp[-2(I_new - I_Coul)]`; (5) the near-top parabolic normal form action `I_top(E) = π(V_peak-E)/ℏ_eff · sqrt(m_s/K_peak)`; (6) the carried-forward dynamic diagnostics `Xi_turn(E) = Xi_1(r_+(E))` and `lambda_th(E) = |E/V'(r_+(E))|` (eq stage248-diagnostics). The appendix row (line 94) labels the deliverable set "Energy integral, peak and threshold-speed compiler, turning points, WKB action, transmission ratio, dynamic diagnostics, and Session-II crossing/tunneling benchmark." The notes additionally pin a Session-II numeric benchmark (m_s=1, ℏ_eff=1, r0=5, r_contact=0.18, E_sub=2.5) with reported observables and a ~23.31% transmission enhancement.

## What the script claims to verify

The SymPy docstring lists six items: (1) energy conservation, (2) finite-radius threshold-speed compiler, (3) Coulomb contact-threshold formula, (4) Coulomb turning-point and WKB reference formulas, (5) near-top parabolic action normal form, (6) Session-II benchmark relations. The assertions verify: `dE/dt` vanishes on-shell; `v_crit,new` from the boxed formula equals the `solve`-derived root and the energy at it equals `V_peak`; the same for the Coulomb contact threshold; the proposed Coulomb antiderivative `F_coul` differentiates to the integrand and its endpoint evaluation equals the closed-form `I_Coul`; the transport laws `dr_±/dE = 1/V'(r_±)`; the Gaussian near-top integral equals `π·ΔE·sqrt(m_s/K_peak)/ℏ_eff`; and a numeric benchmark block reproducing the reported Session-II values plus the window ordering. The Mathematica script mirrors this section-by-section. Neither engine symbolically exercises the carried diagnostics `Xi_turn` / `lambda_th` (eq stage248-diagnostics) — they appear only as hardcoded benchmark numbers checked for positivity.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Energy conservation `dE/dt=0` | py:58 / wl:91 `dE_on_shell == 0` | match |
| `v_crit,new` threshold law | py:104,103 / wl:114,115 (solve vs formula, energy=V_peak) | match |
| `v_contact,Coul` Coulomb threshold | py:106,107 / wl:117,118 | match |
| Classical window `v_crit<v0<v_contact` | py:287 numeric / wl:242 sign-flag | match (numeric only; symbolic in §2.3 of notes is just an ordering) |
| Turning-point transport `dr_±/dE=1/V'` | py:178,179 / wl:175,176 | match |
| WKB action `I_new` definition | py:119 (printed Integral, not asserted) | partial (form displayed, not independently verified — but it is a definition, nothing to "prove") |
| Coulomb turning point `r_turn,Coul=1/E` | py:125 (printed) / wl:123 (printed) | match (trivial definition) |
| Coulomb closed-form `I_Coul` | py:176,177 (antiderivative + endpoints) / wl:173,174 | match |
| Transmission ratio `exp[-2(I_new-I_Coul)]` | py:162-164 (printed symbol) / wl:171 | match (it is the definition; numeric ratio checked at py:285) |
| Near-top action `I_top` | py:206 / wl:195 | match |
| Diagnostics `lambda_th = |E/V'(r_+(E))|` | py:290-291 `Xi>0, lambda>0` only / wl:243-244 `Sign-1` only | mismatch — the boxed symbolic identity is never exercised; only positivity of a hardcoded number |
| Diagnostics `Xi_turn = Xi_1(r_+(E))` | py:290 / wl:243 (positivity only) | partial (carried tag; notes downgrade to "not a theorem", but card boxes it) |
| Session-II benchmark numbers | py:279-291 / wl:232-244 | match |

Dominant pattern: most deliverables match; the carried diagnostics are under-exercised → `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `dE_on_shell == 0` | energy conservation | yes |
| A2 | sympy | 103 | `E(v_crit,new)-V_peak == 0` | threshold law | yes |
| A3 | sympy | 104 | `vcrit_new_solved - vcrit_new == 0` | threshold law (solve vs formula) | yes |
| A4 | sympy | 105 | `delta_new.subs(v0,vcrit) == 0` | threshold law (redundant with A2) | partial |
| A5 | sympy | 106 | `vcontact_solved - vcontact == 0` | Coulomb threshold | yes |
| A6 | sympy | 107 | `delta_coul.subs == 0` | Coulomb threshold (redundant with A5) | partial |
| A7 | sympy | 176 | `dF_coul == 0` | Coulomb antiderivative | yes |
| A8 | sympy | 177 | `I_coul_endpoints - I_coul_formula == 0` | Coulomb closed-form action | yes |
| A9 | sympy | 178 | `transport_plus - 1/V'(r_+) == 0` | outer transport law | yes |
| A10 | sympy | 179 | `transport_minus - 1/V'(r_-) == 0` | inner transport law | yes |
| A11 | sympy | 206 | `I_top - I_top_expected == 0` | near-top normal form | yes |
| A12 | sympy | 279-289 | numeric `abs(...) < tol` (×11) | Session-II benchmark + window | yes |
| A13 | sympy | 290-291 | `Xi_turn_num > 0`, `lambda_th_num > 0` | dynamic diagnostics | no (positivity of hardcoded number; symbolic identity untested) |
| B1 | mathematica | 91 | `expectZero dE/dt` | energy conservation | yes |
| B2 | mathematica | 114-118 | `expectZero` ×5 | threshold laws | yes (B with redundancy mirroring A2/A4, A5/A6) |
| B3 | mathematica | 173-176 | `expectZero` ×4 | antiderivative, endpoints, transport | yes |
| B4 | mathematica | 195 | `expectZero near-top` | near-top normal form | yes |
| B5 | mathematica | 232-242 | `expectNear`/`expectZero` ×11 | benchmark + ordering | yes |
| B6 | mathematica | 243-244 | `Sign[Xi]-1`, `Sign[lambda]-1` | dynamic diagnostics | no (positivity only) |

SymPy assert count: 17 (`assert` statements). Docstring claims 6 verification areas — all six are touched, but area 6's "dynamic turning-point diagnostics" sub-claim is only positivity-checked (A13). Mathematica `expectZero`/`expectNear` count: 25 — same coverage profile.

## Findings

### F1 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_sympy_audit.py:226-227,255,290-291`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl:211-212,230,243-244`

**What's wrong:**
The paper card boxes the diagnostic identity (eq stage248-diagnostics):
> `lambda_th(E) = |V(r_+(E))/V'(r_+(E))| = |E/V'(r_+(E))|`

and `Xi_turn(E) = Xi_1(r_+(E))`. This identity carries real content: the second equality uses the turning condition `V(r_+(E)) = E`. The script never exercises it. SymPy only does:
```python
Xi_turn_num = 0.34437471      # line 226
lambda_th_num = 0.42826825    # line 227
...
assert Xi_turn_num > 0        # line 290
assert lambda_th_num > 0      # line 291
```
and prints `Vprime_turn_mag = Esub_num / lambda_th_num` (line 255) without asserting it. Mathematica mirrors this with `Sign[XiTurnNum]-1` / `Sign[lambdaThNum]-1` (lines 243-244). Asserting that a hardcoded number is positive does not verify any Stage-248 relation; it cannot fail for any physics. The boxed identity `lambda_th·V'(r_+) = V(r_+) = E` (the load-bearing equality) is left unverified in both engines.

**Why this matters:**
This is a checkpoint stage (higher bar) and "dynamic diagnostics" is explicitly listed as a deliverable in the appendix row (line 94) and boxed in the card. A `> 0` check on a constant is `insufficient_verification`: it does not exercise the identity the paper states. A downstream stage (notes §8 / Stage 249 helicity compiler) consumes `lambda_th` as the "first dynamic confinement-width trigger," so the relation matters.

**Required change:**
Add a symbolic check (both engines) that the threshold identity is a genuine consequence of the turning condition, not a tautology. With `V` an undefined function, `rp` an outer-turning-point symbol, `Vp = V'(rp)`, and the turning condition `V(rp) = E`:
- Verify `simplify( (V(rp)/Vp) - (E/Vp) ).subs(V(rp), E) == 0` — i.e. the two boxed forms of `lambda_th` agree once `V(rp)=E` is imposed (and confirm it is NOT identically zero before the substitution, so the check is non-vacuous).
- Verify `Xi_turn` is sampled at `r_+`: assert the symbolic definition `Xi_turn == Xi1(rp)` where `Xi1` is an undefined function evaluated at the outer turning-point symbol (this fixes the carry-forward target as `r_+(E)`, matching eq stage248-diagnostics). This is a definitional anchor, but it makes the carried tag explicit rather than a bare positive number.

Do NOT remove the existing benchmark numeric readbacks; add the symbolic checks alongside them.

**Verification:**
New `assert` (SymPy) and `expectZero` (Mathematica) lines appear exercising `lambda_th·V'(r_+) - E` collapsing to 0 under `V(r_+)=E`; both scripts still exit 0. The verifier confirms the new check is non-tautological by noting the pre-substitution residual is nonzero.

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl:95-118,199-244`

**What's wrong:**
The `.wl` is, across its symbolic sections, a step-by-step transliteration of the `.py`: identical section order, identical intermediate-variable choreography, and identical assertion order. Section II is a near-verbatim port:
```python
# py:68-89
E_launch_new = ... ; vcrit_new = sqrt(2(Vpeak-V0)/m_s)
vcrit_new_solved = [sol for sol in sp.solve(Eq(E_launch_new, Vpeak), v0) if not minus][0]
E_at_vcrit = E_launch_new.subs(v0, vcrit_new) - Vpeak
```
```mathematica
(* wl:95-101 *)
ElaunchNew = ...; vcritNew = Sqrt[2(Vpeak-V0)/ms]
vcritNewSolved = (v0 /. First[Solve[ElaunchNew == Vpeak, v0]]) /. ConditionalExpression[val_,_]:>val
EAtVcrit = ElaunchNew /. v0 -> vcritNew
```
Section V (the benchmark block, lines 199-244) is a one-to-one variable-name echo of py:214-291 (`mNum`↔`m_num`, `VpeakNum`↔`Vpeak_num`, etc.) with the same hardcoded literals in the same order. For the benchmark readback layer this echo is acceptable (both must read back the same reported numbers), but the symbolic sections sharing identical choreography weakens the second-engine independence value.

Mitigating: the two load-bearing symbolic checks DO use engine-native routes — the Coulomb antiderivative needed `PowerExpand`/`Refine`/a manual `Sqrt[rv^(-1)] Sqrt[rv] -> 1` rewrite (wl:134-140) absent from SymPy's single `simplify`, the transport solve uses `Derivative[1][rpFun][Eturn] -> drPlus` replacement plumbing different from SymPy's `Function`/`diff`, and `I_Coul` surfaces as `ArcCos` in Mathematica (wl-output:85) vs `π/2 - asin` in SymPy. So this is a partial transliteration, not a total mirror.

**Why this matters:**
Checkpoint stages demand genuine dual-engine independence. The current `.wl` adds little independent confidence in Section II (threshold compiler) and Section V (benchmark), where it merely re-runs SymPy's exact algebra in Mathematica syntax.

**Required change:**
Re-author Section II of the `.wl` to a native, differently-decomposed route: rather than `Solve[ElaunchNew == Vpeak, v0]` then comparing to the hand-written `Sqrt[...]` (mirroring SymPy), derive `v_crit,new` by `Reduce[(1/2) ms v^2 + V0 == Vpeak && v > 0, v]` (native solver, different primitive) and confirm the energy at the resulting root equals `V_peak` symbolically. Likewise for `v_contact,Coul`. Leave Sections I, III, IV as-is (their hard steps are already engine-native) and leave Section V benchmark numbers as-is (legitimate shared readback). Also fix the stale banner (F3).

**Verification:**
Section II of the `.wl` no longer uses the `Solve[...] /. ConditionalExpression` pattern transliterated from SymPy's `solve(...) if not could_extract_minus_sign`; it uses `Reduce`/native root extraction. Script still exits 0 with the same PASS lines.

### F3 — stale_output / mislabel (cosmetic)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage248_dynamic_event_chain_compiler_from_relaxed_stationary_barrier_front_end_turning_point_threshold_speed_and_wkb_mathematica_audit.wl:65`

**What's wrong:**
The top banner reads `banner["STAGE 231 — DYNAMIC EVENT-CHAIN COMPILER ..."]` and the captured output (mathematica output line 11) prints "STAGE 231" for a stage-248 file. This is a leftover label from a cloned-from-231 file — a small but real provenance smell on a checkpoint stage, and it confirms the transliteration origin noted in F2.

**Why this matters:**
A reader scanning the Mathematica transcript sees the wrong stage number. On a checkpoint this should be exact.

**Required change:**
Change `wl:65` banner string from `"STAGE 231 — ..."` to `"STAGE 248 — DYNAMIC EVENT-CHAIN COMPILER FROM THE RELAXED STATIONARY BARRIER FRONT END"`.

**Verification:**
`.wl` line 65 and the regenerated output header read "STAGE 248".

## Independent-derivation check (Mathematica)

Partial independence. The Coulomb antiderivative check is genuinely engine-native: SymPy uses one `sp.simplify(diff(F_coul,rv) - coul_integrand)` (py:136), while Mathematica required a four-step native sequence — `PowerExpand[D[Fcoul,rv]-coulIntegrand]`, `Refine[..., Esub rv < 1]`, the manual rewrite `Sqrt[rv^(-1)] Sqrt[rv] -> 1`, then `FullSimplify` (wl:134-140) — and the closed form surfaces as `ArcCos` rather than `π/2 - ArcSin`. The transport-law check uses different symbolic plumbing (`Derivative[1][rpFun][Eturn] -> drPlus` vs SymPy `Function`/`diff`). However Sections II and V are step-by-step ports with identical variable choreography and (for V) identical literals (see F2). Net: the load-bearing WKB/turning-point checks are independent; the threshold compiler and benchmark layers are transliterated → F2 (low severity).

## Engine cross-check

Engines agree. Both report `dE/dt on-shell = 0`, all five threshold checks 0, the Coulomb antiderivative and endpoint checks 0, both transport laws 0, the near-top action 0, and the benchmark `T_new/T_Coul ≈ 1.2331`, `improvement ≈ 23.3128%`, `I_Coul(exact)` within `1e-3` of the reported `0.30222297`. SymPy reports `I_Coul(exact formula) = 0.302305796` and Mathematica `I_coul exact vs report diff = 8.28e-5` — consistent. The `r_turn,Coul` benchmark agrees at `0.4`. No sign/factor disagreement found.

## Verdict justification

Findings, not clean, not stop-cold. Attacks that FAILED to break the script: the threshold-speed checks are non-tautological (boxed formula compared to an independent `solve` root, A3/A5); the Coulomb antiderivative check genuinely differentiates a proposed closed form against the integrand (A7) and the endpoint evaluation reproduces the boxed `I_Coul` (A8); the transport-law derivatives are over an undefined function `V`, so they are not the variable-independence self-test trap (the chain-rule derivative is a real symbolic object, not identically zero); the near-top Gaussian integral is computed, not assumed (A11); symbol domains (`V0,Vpeak` real; `r0,r_contact,v0,m_s` positive; `y` real on a symmetric parabolic interval) are physically correct. The benchmark numerics reproduce the reported Session-II values within stated tolerances and the window ordering `v_crit < v_cross < v_contact` holds. What does NOT fully hold up: the carried diagnostics `Xi_turn`/`lambda_th` are only positivity-checked on hardcoded numbers, while the card boxes the symbolic identity `lambda_th=|E/V'(r_+)|` (F1, medium); the `.wl` transliterates Sections II/V of the `.py` and carries a stale "STAGE 231" banner (F2/F3, low). I read the paper card, the notes, and the appendix rows for this unit before opening the scripts; the verified content matches the paper on every load-bearing identity, with the one diagnostic-coverage gap noted. No `paper_misalignment` (the boxed values all match the script; the notes' `× 168%` in §6.3 is a prose typo in a notes file the red-team does not edit, and the script correctly uses `× 100`). Hence `verdict: findings`, `paper_alignment: partial`, no stop-cold.

## Self-test notes

Variable-independence trap: checked A9/A10 — `diff(V(r_plus(E_turn)), E_turn)` over the undefined function `V` is a real chain-rule object `V'(r_+)·r_+'(E)`, not identically zero; the transport law is genuine, not the vacuous self-test trap. For the F1 fix I confirmed the prescribed `lambda_th` check is non-vacuous: `(V(rp)/Vp) - (E/Vp)` is nonzero in general and collapses to 0 only after imposing `V(rp)=E`, so the directive instructs Codex to assert both the pre-substitution nonzero and post-substitution zero. Symmetry/parity: the near-top integral integrand `sqrt(2m(ΔE - Kpeak·y²/2))` is even in `y` on `[-yturn, yturn]`, so the nonzero result is consistent; verified the endpoints `y=±yturn` give a zero integrand argument (no domain error). Trivial-case: substituting the Session-II numbers reproduces the boxed `v_crit ≈ 2.5414`, confirming the threshold-law assertions are anchored to real values.
