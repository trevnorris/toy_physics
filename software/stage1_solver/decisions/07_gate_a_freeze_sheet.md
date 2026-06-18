# Decision 07 — GATE-A freeze sheet (Claude+Codex agreed; PENDING user GO)

**Date:** 2026-06-17
**Mechanism:** Claude+Codex consult ([[claude_codex_resolve_math]]; user directed the GATE-A *values* are a
math/documentation determination for us to agree on — the user's role is the GO + the target-blind guarantee).
Prompt + log: `_scratch/gate_a_values_consult_prompt.md` / `…_consult.log`. Claude reviewed + AGREES with the
determination below.
**Purpose:** the frozen free_choice set for the M1c physical run, with provenance + target-blind attestation.

## Determination
The record pins the **operator form + the freeze discipline**, but does NOT pin the wall constitutive
functions' forms/magnitudes, nor a numeric `R_exit`/`R0(w)` family. Those are genuinely **free_choice → set
target-blind here** (the simplest, most neutral choice = unity, blind to R_norm). Current torch config values
(`r_mouth=1.2, r_exit=0.9`) are engineering-smoke placeholders, NOT physical pins — M1c-run must replace them
with this freeze sheet.

## Freeze sheet (class: bench=external benchmark, D=branch-determinable, F=free_choice posited, C=convention)

| input | FROZEN value/form | class | status |
|---|---|---|---|
| `μ_η(w)` | `1` (const on [0,L]) | F | set target-blind (simplest) |
| `T_w(w)` | `1` | F | set target-blind |
| `T_Ω(w)` | `1` | F | set target-blind |
| `K_η(w)` | `1` | F | set target-blind |
| wall positivity | `μ_η>0, T_w>0, K_η+6T_Ω≥0` (1>0, 1>0, 7≥0 ✓) | C | check |
| `a` = `R_mouth` | `1` (natural units) | C | pinned |
| `L/a` → `L` | `37/20` → `L=1.85` | D | pinned (branch-determinable posit) |
| `ε_r`, `ℓ` | `1/20`, `a/20` (if radial support layer used) | D | pinned |
| `R_exit` | `3/2` (=1.5) | F | set target-blind (STRUCTURED) |
| `R0(w)` family | cubic smoothstep `R0(w)=1+(1/2)(3x²−2x³)`, `x=w/L`, on `[0,1.85]` (a→R_exit=1.5) | F | set target-blind (STRUCTURED) |
| boundary class | mouth Dirichlet; exit open-impedance/Robin; NO hard cap | C | pinned |
| `G, c, c_s` | `1, 1, 1` (natural units; `c_s=c` canonical branch); benchmark stays GR quadrupole `54Gc_s⁵/5a⁵c⁵` | bench/C | pinned |
| `mhat0` | `1` (point-particle natural source-map limit) | C | pinned |
| `S_port` | `1` (N0 exported in gravitationally-normalized port convention) | C | conventional |
| `theta_tail` | `1`, INACTIVE (`tail_sector_active=false`) | C | pinned/inactive |
| `chi_Q, N_Q` | freeze the RELATION `N_Q=chi_Q⁻¹`; **do NOT freeze `chi_Q=1` as an input** (it's the branch-determinable TARGET to extract + compare) | D/C | extract-then-compare |

## Target-blindness (G2/G3) — PROCEDURAL guarantee
Values that can move R_norm (wall functions, R0/R_exit/L/a, c_s/c, mhat0, S_port) are guaranteed target-blind
by the FREEZE DISCIPLINE: freeze + hash BEFORE the physical solve, no residuals/pass-flags in the frozen
packet, no mutation after residuals (prereg §:134-136,176-180; decision-06). **The M1c result is CONDITIONAL:
a miss falsifies THIS effective-closure branch WITH these posited functions/geometry/conventions; a pass =
this branch+posits reproduce the GR-normalized target target-blind.**

## Freeze-hash basis (G4)
Canonical JSON over: `parent_action_status=effective_closure`; branch identity; the 4 wall functions
(domain/units); geometry `{a,R_mouth,L,R_exit,R0_family,ε_r,ℓ,boundary_class}`; gauge `H=Z`; source/port
`{G,c,c_s,mhat0,S_port,theta_tail,tail_sector_active}`; extraction formulas; solver tols/mesh/backend/source
revision; validation protocol. **NO target residuals or pass/fail fields.**

## Conceptual flag (G5)
Freezing posited wall constitutive functions under `effective_closure` is **in-scope Path B, NOT a conceptual
change** (parent-status decision). Escalate to user ONLY if we instead promote `S_Σ[R]`/parent wall dynamics,
activate excluded relaxed/open-system lanes, or change the pre-registered branch/tail scope (= Path A / new
pre-registration).

## GATE-A GO: STRUCTURED branch (user-authorized 2026-06-17; Claude+Codex consult `_scratch/structured_throat_consult*`)
The user authorized GATE A on a **structured throat** (vs the degenerate constant conduit). Codex+Claude agreed
the MINIMAL structured branch: change ONLY the geometry — `R0(w)=1+(1/2)(3x²−2x³)` (cubic smoothstep), round
**`R_exit/a = 3/2`** (principled, modest, NOT fitted); wall functions stay flat (`=1`, the only added free
choice vs trivial is `R_exit` + the standard profile → max target-blind). Valid open throat, non-degenerate,
in-scope Path B. The two rows above are updated to this structured freeze; all other decision-07 pins
unchanged. Set on principle, BEFORE the solve, blind to `R_norm`; covered by the G4 freeze-hash discipline.

**GATE A is GO.** NEXT = M1c-run: write the freeze sheet + hash FIRST → set the torch config to this frozen
branch (NOT the smoke placeholders) → re-solve WP1 → derived chain → tracked byte-reproducible frozen packet
(`software/stage1_solver/frozen/m1c/<hash>/`) → V2 chain → PHYSICAL `R_norm` + §J budget. No value mutated
after any residual.
