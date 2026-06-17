# Decision 04 — mixed-Maxwell frontier: derivable WITHOUT Path A (refines decision 03)

**Date:** 2026-06-17
**Mechanism:** Claude+Codex read-only consult ([[feedback_claude_codex_resolve_math]]). Prompt + full log:
`software/stage1_solver/_scratch/mixed_maxwell_frontier_consult_prompt.md` / `…_consult.log`.
**Question:** can the mixed-Maxwell U/W port sector (`Ω_U,Ω_W,R,g_U,g_W` → `N0` → `R_norm`) be DERIVED from
the parent localized-Maxwell action WITHOUT opening Path A (`S_η^(ψ,A)`)?

## Verdict: (i) DERIVABLE NOW WITHOUT PATH A — but not as the assumed V2-09 U/W eigenproblem

**Path A is NOT the blocker.** The real blocker is a missing, bounded **Maxwell-port reduction**. Key findings:

1. **U/W are a channel-BASIS choice, not eigenmodes (F1).** If `U,W` were exact eigenmodes of the linearized
   localized-Maxwell operator the internal block is diagonal → `R=0`. A nonzero `R` means V2-09's two-
   coordinate `{Ω_U,Ω_W,R,g_U,g_W}` Lagrangian is a chosen diabatic basis. The parent action (`H=Z` localized
   Maxwell, compact:681/688) fixes the full bilinear form + the mixed observables (`E_w=F_{w0}`, `C_a=F_{aw}`)
   but does NOT name a unique two-coordinate basis. ⇒ the clean target is the **basis-INVARIANT
   wall→outgoing mixed transfer (`N0`)**, not basis-dependent port parameters.

2. **The forward coupling needed for `N0` is parent-determined (F2).** Direct `η→δZ` Maxwell-localization
   modulation is NOT in the declared action. But the available forward route IS:
   `η → δV_conf → δψ/δJ_ψ → δA` (wall displaces matter → matter current changes → current sources gauge).
   Every link is fixed: `δV_conf` (compact:1080), the exact current + Maxwell source (compact:655/622),
   already realized in torch (`_matter_number_current` coupled_branch.py:210, Maxwell residual :286).

3. **`S_η^(A)` (the return source / Path A) is NOT required for `N0`/`R_norm` (F2).** `N0 = (Ω_U²g_W+R g_U)²
   /Δ₀²` needs only the FORWARD forced response (`g_U,g_W` from the forced Maxwell/BdG response). `S_η^(A)`
   is needed only for the SELF-CONSISTENT wall equation + reciprocity audit — not the forward transfer number.
   **This refines decision 03's `g_U/g_W = S (Path-A)`: for the `N0` purpose they are forward-derivable; only
   the self-consistent wall closure / reciprocity is Path-A.**

## The bounded subprogram (the "mixed-Maxwell port spike")

What must be built (F3): (1) parent `H=Z` localized Maxwell operator; (2) the REAL ℓ=2 vector-harmonic
reduction (not the torch scalarized engineering-smoke wrapper); (3) open Maxwell BCs + gauge constraint +
physical inner product; (4) the linearized matter current response from `δV_conf`; (5) extract a
**basis-invariant wall→outgoing mixed Green-function transfer** (preferred — sidesteps the U/W basis
ambiguity), or a frozen target-blind U/W basis → `Ω_U,Ω_W,R,g_U,g_W`; (6) feed `Δ, P, N0` to M1b/M1c.

**Size/risk (F4):** a focused research subprogram (NOT 1-2 lines). Smaller than Path A, but HIGHER RISK than
the BdG/wall extraction (gauge fixing, vector harmonics, open BCs, basis ambiguity). This is where a clean
`R_norm` falsification test — the breakthrough — actually comes from, and it is reachable without Path A.

## Decision

- **Strategic implication:** the breakthrough (clean `R_norm`) is reachable WITHOUT the gated Path-A program,
  via the bounded mixed-Maxwell port spike. This is a genuine positive update.
- **PENDING USER DECISION:** whether to commit to the mixed-Maxwell port spike now (it's a real subprogram
  with the listed risks), and whether to target the basis-invariant Green-function transfer (recommended) vs a
  frozen U/W basis. This is a scope commitment = user-level call; teed up after M1b is reviewed.
- **Until/unless the spike is taken:** M1b/M1c keep `mixed_ports_status: posited_not_derived` and `R_norm` as
  a partially-posited diagnostic (decision 03 stands for the posited path).
