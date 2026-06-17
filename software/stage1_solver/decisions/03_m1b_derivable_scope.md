# Decision 03 — M1b derivable scope (which V2-22B primitives are derived vs posited vs Path-A-blocked)

**Date:** 2026-06-17
**Mechanism:** Claude+Codex read-only consult ([[feedback_claude_codex_resolve_math]]); a math/derivation
question (what the canonical record actually determines), not a conceptual user gate. Prompt + full log:
`software/stage1_solver/_scratch/m1b_derivability_consult_prompt.md` / `…_consult.log`.
**Context:** the Mathematica second-engine route (M0 ✅, M1a ✅) needs M1b to DERIVE (not posit) the V2-22B
branch packet. Before authoring the M1b build, we resolved exactly what is derivable. Triggered by a schema-
fidelity finding that the V2-22A/V2-21 fixtures CONSUME modes/ports as inputs (Gaps M1/M2).

## Finding — per-primitive classification (D=WP1-derivable, F=free_choice, F†=real derivation gap, S=Path-A, C=convention)

| Primitive | Class | Basis (file:line) |
|---|---|---|
| Wall constitutive `μ_η,T_w,T_Ω,K_η` | **F** | "fixed posited inputs" `docs/stage1_preregistration.md:64` (frozen at GATE A) |
| Wall `K, M` | **D** | weak-form integrals `M_ij,K_ij` `stage_v2_20_…derivation.md:118-136` |
| BdG `varpi_α, φ_α` | **D** | canonical `L_BdG` `notes/moving_throat_pde_program_compact.md` §4.7 ~L1406-1428; V2-20:196-204 |
| BdG overlap `I_ηφ`, coupling `c_α` | **D** | `c_α` = projected wall coupling via `δV_conf` (compact L1080-1085/L1424-1428); V2-20:206-230 |
| `λ_B` (in `c_α=λ_B·I_ηφ`) | **C** | adapter normalization convention; `c_α` is the physical object |
| Mixed profiles `u_r, w_r` | **F†** | no canonical 1-D mixed eigenproblem produces them (V2-20:246-255 "must output") |
| `Ω_U, Ω_W, R` | **F†** | V2-09 ASSUMES the reduced Q/U/W Lagrangian (`:51-68`); does NOT derive the spectral data; "full PDE branch realization open" `stage_v2_09_…:483-493` |
| Mixed overlaps `I_ηu,I_ηw,I_uw` | **D** | adapter integrals `stage_v2_22a_…:53-68` |
| `g_U, g_W` | **S** | wall↔gauge couplings; full derivation needs the open `S_η^(A)` ("named, not formulated", `decisions/02:24-28`) |
| `λ_U,λ_W,λ_R` | **C** | adapter factorization variables; products inherit the physical status above |

`F†` = a REAL derivation gap (the parent action may determine it, but the reduction hasn't been written),
distinct from a clean constitutive free choice like `μ_η`.

## Consequences

1. **`R_norm` is NOT yet a clean falsification test.** It is S_η-independent only at the packet/response
   level (it doesn't use the WP3 low-frequency response) — but its mixed-port inputs `g_U/g_W` are not cleanly
   derivable now. This **refines** the step-8c "R_norm S_η-independent" classification (decisions/02:38): true
   in the narrow sense, but the mixed inputs remain underived. The ANSATZ_LEDGER "branch-determinable" label
   means "the real PDE WOULD select it, posited pending the full solve" (`ANSATZ_LEDGER:17-22`), NOT "already
   derived."

2. **The frontier is localized.** The entire remaining gap to a clean `R_norm` = the **mixed-Maxwell U/W
   eigenproblem** (deriving `Ω/R` from the parent localized-Maxwell action §2.5) **+ the wall↔gauge coupling
   `g_U/g_W`** (forward maybe derivable like `δV_conf`; full return = `S_η`/Path A). That sub-program, not the
   M1b mechanics, is where a breakthrough lives.

## Decision

**Proceed in parallel** (one orchestrator, two workstreams — user-approved 2026-06-17):
- **M1b at the largest honest scope (Q4):** DERIVE wall `K,M` + BdG `varpi_α/φ_α/c_α` (genuine eigensolve) +
  overlaps + V2-21 algebra; carry mixed `Ω_U/Ω_W/R/g_U/g_W` as explicitly-labeled `mixed_ports_status:
  posited_not_derived`; emit the packet + run the V2 chain; label any `R_norm` a **target-blind
  partially-posited diagnostic, NOT a clean test**. Directive: `directives/mt15_02_bdg_wall_derivation.md`.
- **Frontier consult (read-only):** scope whether the mixed-Maxwell sector is derivable from the parent
  action WITHOUT Path A. Prompt: `_scratch/mixed_maxwell_frontier_consult_prompt.md`. Result → a follow-on
  decision (04) on whether/how to attack the mixed-Maxwell derivation.

Still pre-freeze, target-blind; GATE A (freeze free_choice) and the M1c frozen run remain user-gated and
unchanged. Conceptual flag: deciding to ATTACK the mixed-Maxwell/Path-A frontier (vs continue effective
closure) is a user-level strategic call — deferred to decision 04 after the frontier consult.
