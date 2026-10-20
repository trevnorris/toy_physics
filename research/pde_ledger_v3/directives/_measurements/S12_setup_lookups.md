# Measurements: S12 setup (generated 2026-10-19 21:58)

Generator: `_scratch/s12/gen/s12_setup_lookups.sh`. Literal retrieval only.

## K1. The plan's S12 entry: the committed law, the ruled-out law, the two inventories
```
$ grep -n 'S12 · The throat' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
579:### S12 · The throat's non-variational SOURCE, plus its boundary data
```

```
$ grep -n -A1 '∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
591:∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n Γ_B ,    Γ_B = Γ_return − Γ_drain ,    [Γ_B] = T⁻¹
592-```
```

```
$ grep -n 'Phase-conversion, \*\*not suction\|RULED OUT\*\* (`:19`\|frozen-wall-premised\. and\|must be \*\*derived\*\*' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
595:material. Phase-conversion, **not suction**.
597:⛔ **RULED OUT** (`:19`, user, 2026-07-24): the G0-card `S_drain` **ρ-mass-sink + remote return**, i.e.
600:up all the calculations."* ⇒ If a static limit is ever wanted it must be **derived** as the `ω→0` limit
```

```
$ grep -n -A5 'Two inventories, kept apart' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
609:⛔ **Two inventories, kept apart:**
610-- **non-variational source / controller functions** — `Γ_return`, `Γ_drain`, `J_χ`, and the return
611-  controllers. ⛔ All `[OPEN]` in form.
612-- **boundary / domain data** — the mouth, collar and IR conditions, and the `𝔅` branch selection.
613-
614-⇒ The throat enters this ledger with its **source form open and its strength supplied**. ⛔ Not "only a
```

```
$ grep -n 'A13 (the conversion balance' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
616:**Register:** ⛔ **A13 (the conversion balance must use the branch chosen at S1)** · **C6**; **C1** (⚠ the drain law is *parameterized, not derived* — and `g_phys` was never
```

## K2. The source of the committed law (v2 note; a lead, not an authority)
```
$ grep -n 'THE COMMITTED CHOICE\|RULED OUT\*\* (frozen-wall' /var/projects/toy_physics/research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md
23:- **✅ (a) Part-I `Γ_B` order-conversion drain — THE COMMITTED CHOICE** (DYNAMICAL wall): `∂_t(χ_B n) + ∇₄·(χ_B n u + J_χ) = n
26:- **❌ (b) G0-card `S_drain` ρ-mass-sink + remote return — RULED OUT** (frozen-wall regime; card §6 freezes `r_B` via `S_hold`
```

```
$ grep -n 'finite-rank return controllers\|The force partition' /var/projects/toy_physics/research/pde_ledger_v2/notes/stage045_nonvariational_block_prep.md
45:- **The finite-rank return controllers `I_ret` / `P_ret`** — projectors enforcing EXACT global mass / momentum / energy closure
51:- **The force partition `F_total = F_var + F_flux + F_𝔅 + F_rad`** — the total inter-defect force split into the variational
```

## K3. Γ_B is signed
```
$ sed -n '14p' /var/projects/toy_physics/docs/oriented_throat_return_source_correction.md
In the recorded equations, `chi_B` measures brane order. Positive `Gamma_B` is bulk-to-brane re-ordering; negative `Gamma_B` is de-structuring. The projected ordered-density source has both `S_flux` and `S_convert`, with `S_convert = Integral(W n Gamma_B dw)`. The older `S_leak` is the boundary/window flux identity; it should not be silently equated with the conversion rate alone. Neither term is the transverse-light
```

## K4. S14a: the bridge the corrected S12 law creates
```
$ grep -n 'conserves total material and converts only order\|Gauss\|Derive, from the dynamical system' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
630:The committed dynamical `Γ_B` law **conserves total material and converts only order**
631:(`ledger_stage006_two_phase_chiB_ontology.md:81`). But the imported force derivation assumes a **Gauss
636:⇒ **Derive, from the dynamical system:** the projected order-loss source · the definition of `J` · the
733:**Gaussian/profile support width** controlling `Q_ij = M a² δ_ij / 2`; the model's **mouth radius** is a
1047:and develops **NO emergent first-class U(1) Gauss constraint** (`parameter_register.md:230`, R66 `:332`;
```

## K5. S_leak and the dark-energy postulate
```
$ grep -n 'Nobody set out to derive a leak\|likely the same physics as S12\|this mechanism I am postulating' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
544:⭐ **Nobody set out to derive a leak; it is what is left over because the window has edges.**
558:> *"this mechanism I am postulating is the reason for dark energy. As the bulk reorders onto the brane,
564:(user-confirmed 2026-08-03; likely the same physics as S12's `Γ_B`, ⚠ unverified).
```

## K6. A13 in the plan: the S1 gate
```
$ grep -n 'A13 — GATE\|A13 gate' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
173:⛔ **Two, not three.** `ψ` plus the order field in whichever form S1's **A13 gate** selects — state the
176:**Register:** ⛔ **A13 — GATE. Resolve before this step banks.** `χ_B` is used here as a **real material
272:⚠ **The equation below is a TEMPLATE, not a pre-decided branch** — S1's A13 gate chooses first, and this step rewrites to match.
```

## K7. Register entries S12 owns or shares
```
$ grep -n '^| \*\*S12\|^| \*\*S22\*\*, with Q2 and S12' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
87:| **S12** *(register inference; no live owner named)* | `R-S12-01` | a reservoir and power budget if a non-passive interface response is adopted | S11b-B / unified S11b, O2 carried condition |
88:| **S12** | `R-S12-02` | background drain/return and separate boundary data | S11b-A/B, S11c uniform, O2 |
90:| **S22**, with Q2 and S12 | `R-S22-01` | how the energy a throat releases splits among light, the extra longitudinal wave and the bulk | S11 departure; photon survey E32–E34, E38 |
```

```
$ grep -n '^- \*\*requirement\*\*' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md | grep -n 'reservoir\|drain/return'
14:502:- **requirement** — a named reservoir and a stated power budget for any adopted non-passive interface
15:531:- **requirement** — the native drain/return functions and their separate boundary data.
```

```
$ grep -n 'O(v₀|q_n|/ω)' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
248:  rest-frame linearisation discards a relative correction of order `O(v₀|q_n|/ω)`, where `v₀` again denotes
536:  `O(v₀|q_n|/ω)` (`steps/S11bB_interface_assembly.md:150–154`), where that record's `v₀` is the steady
```

## K8. O2: the drive premise and the S12 handoffs
```
$ grep -n '| 2 — drive' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
1002:| 2 — drive | Dynamical order-conversion drain; no separate external body force | O2 chooses the absent separate-body-force representation. Source/return and partners S12; S14a/S14 source bridge and S16 response matching remain OPEN, with S16 assignment an inference. |
```

```
$ grep -n '^| O6 `𝒥_map`' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
399:| O6 `𝒥_map` | Exchange/projection and bulk-boundary input; common exchanged material and compatible mass/momentum/traction/work geometry and measures remain required but uncomputed. | S12 dynamical conversion and distinct native source/boundary data; S14a's projected order-loss/far-field bridge remains distinct. |
```

```
$ grep -n 'A13 remains\|A13 branch (S1 gate' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
83:The historical ordered/disordered split of one conserved material is **postulated**; A13 remains
402:Other OPEN operands are A13 branch (S1 gate, propagated S5/S12); conservative momentum/stress and
```

## K9. The S11c closeout's S12 handoff
```
$ grep -n 'S12 is next\|live conversion belongs to' /var/projects/toy_physics/research/pde_ledger_v3/steps/S11c_PARTIAL_CLOSEOUT.md
20:| [Clean-condition symmetry rule](S11c_d_profile_conditioned_scattering.md#clean-condition-result-and-its-owner), measured on SymPy only | No-leak interpretation needs drain frozen; live conversion belongs to **S12**, with no automatic transfer. |
33:**S12 is next:** dynamical bulk-to-brane order conversion and separate boundary data. It does not require a numerical S11c loss factor and does not automatically settle light leakage. [Plan, S12](../V3_STEP_PLAN.md).
```

## K10. The oracle comparison's two questions for S12
```
$ grep -n '^1\. \*\*What drives the flow\|^2\. \*\*Where the material goes\|externally imposed driving potential\|carry the same flow in two ways' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
121:- Visser's force balance (6)–(7), which includes an "externally imposed driving potential" `Φ`;
157:1. **What drives the flow.** Visser's background obeys the force balance (6)–(7), which includes an
158:   externally imposed driving potential `Φ`. The check computes no force comparison
162:2. **Where the material goes.** Visser's continuity has no source. In that setting the §7 velocity forces the
164:   lack (`no_off_sheet_exchange`). So this model can carry the same flow in two ways: with that density and
```

```
$ grep -n 'c0 == Sqrt\|C₁/r^(3/2)' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
78:  nonzero-velocity components, both engines solve the equation as `ρ_br = C₁/r^(3/2)`
83:- **Light speed on this background.** WL's attached flow-only domain prints `c0 == Sqrt[muPerp[r]/rhoBr[r]]`
```

## K11. Step records present in v3 (no S1–S8 record)
```
$ ls /var/projects/toy_physics/research/pde_ledger_v3/steps | grep -c '^S[1-8]_'
0
```

## K12. The distance range S9b and O2 cover: the far field outside one body
```
$ grep -n 'far field, Eulerian steady\|far-field `r>0` domain' /var/projects/toy_physics/research/pde_ledger_v3/steps/O2_steady_brane_balance.md
75:far field, Eulerian steady material profiles, lab time in the brane's far-field rest frame, linear
77:coordinate `w`. Scalar profiles are general live functions of `r=|x|` on the far-field `r>0` domain, with radial
```

```
$ grep -n 'one isolated spherical mass at rest' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
84:The domain is one isolated spherical mass at rest, Eulerian steady draining profiles, far field,
```

## K13. Visser's own statement on constant density with the river speed
```
$ python3 -I -c "import json,sys; d=json.load(open(sys.argv[1])); print(d['V7background']['citation']); print(d['V7background'].get('conditions', d['V7background']))" /var/projects/toy_physics/research/pde_ledger_v3/scripts/S9b_flow_oracle_check_data/oracle_quotes.json | fold -w 400 | grep -i -o '.\{0,160\}incompatible with the continuity.\{0,80\}'
c is position independent. For constant rho and c with that velocity, the paper says: this assignment is incompatible with the continuity equation. Its displayed condition is nabla dot (rho v) != 0. Density-law normal
```

## K14. Earlier dark-matter work: a parked idea and one review objection
```
$ git -C /var/projects/toy_physics grep -n -i -I 'dark matter\|rotation curve\|MOND' -- docs | cut -c1-300
docs/conceptual_foundation.md:243:> curves / a MOND-like modification.** Parked for later.
docs/em_gravity_mined_verdicts.md:27:- **CRACK-7: the MOND/leaky-gravity aside is BACKWARDS-signed** — leaky→extra-dim makes force *steeper* (`1/r³`), not shallower (`1/r` needs effectively 2D). DELETE it. `codex1:2020`.
docs/em_sector_reconsideration.md:63:- **"Long-range Coulomb → `1/r` (MOND-like for EM)" (v1 §2.6):** **REJECTED.** Extra-dimensional leakage makes a force *steeper* (`1/r³`), not shallower; long-range EM deviations are experimentally dead (Faraday rotation, plasmas, lab Coulomb to parts in `10^
```

```
$ sed -n '240,243p' /var/projects/toy_physics/docs/conceptual_foundation.md
>
> **Future item (deferred — do NOT calculate yet, but do not lose):** the same `1/r²`-vs-localization mechanism may, at **large
> separations (galaxy scale, widely-separated bodies)**, turn `1/r²` into `1/r` — a candidate native handle on **galactic rotation
> curves / a MOND-like modification.** Parked for later.
```

## K15. The plan's comparison step for light and gravity
```
$ grep -n -A1 'S21 · What light and gravity share' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
1156:### S21 · What light and gravity share
1157-Both ride the brane; both consume `{ρ_br, μ_R}`. ⛔ **Use the method doc's CORRECTED rule**
```

## K16. Citations in the S12 decision list
```
$ sed -n '78,84p' /var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md
Balance laws:

```
∂_t n + div_4(n u) = 0                                    (total, exact)
∂_t(chi_B n) + div_4(chi_B n u + J_chi) = n Gamma_B ,      [Gamma_B] = T^-1
Gamma_B = Gamma_return - Gamma_drain
```
```

```
$ sed -n '597,601p' /var/projects/toy_physics/research/pde_ledger_v3/V3_STEP_PLAN.md
⛔ **RULED OUT** (`:19`, user, 2026-07-24): the G0-card `S_drain` **ρ-mass-sink + remote return**, i.e.
`∂_t ρ + ∇·(ρv) = S_drain + S_leakage`. It is **frozen-wall-premised**, and the standing principle is
*"the drain cannot be a frozen wall — several instances where we tried to freeze the wall and it screwed
up all the calculations."* ⇒ If a static limit is ever wanted it must be **derived** as the `ω→0` limit
of the dynamical law, ⛔ never posed as an independent frozen solve.
```

```
$ grep -n 'Exchange momentum: the brane' /var/projects/toy_physics/research/pde_ledger_v3/directives/O2_premise_decision_list.md
32:3. **Exchange momentum: the brane's local velocity.** Converted material carries the brane material's local
```

```
$ grep -n 'J_chi = 0. default' /var/projects/toy_physics/research/pde_ledger_v2/notes/stages/ledger_stage006_two_phase_chiB_ontology.md
101:| P9 | `J_chi = 0` default (simplest advective closure); `J_chi ≠ 0` deferred, `[J_chi] = L^-3 T^-1`. |
```

