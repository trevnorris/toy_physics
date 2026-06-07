---
unit_id: 134
batch: IV.4
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T17:03:44-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 134 (USER-AUTHORIZED independent re-author)

Apply the finding below. After applying, append an `## Applied: F1` block with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is genuinely infeasible (e.g. `DSolveValue` cannot close the
BVP symbolically in this Wolfram version), append `## Blocked: F1` with the exact
error and what you tried, then STOP — do not fall back to re-typing the closed form.

Do NOT introduce unrelated features or refactors. Do NOT touch `paper.tex`, `notes/`,
or any prose document — scripts only. After editing, RUN `math -script <path>` and
iterate until it exits 0 with all in-file checks passing. The orchestrator
independently re-runs afterward.

## F1 — mathematica_transliteration (USER-AUTHORIZED full re-author of the `.wl`)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl`

**Issue:** The `.wl` is a transliteration of the SymPy `.py`. The kernel
`sKernel[p_, k_]` (`.wl:31-34`) is character-for-character the SymPy closed form
`S(Pi, kappa)` (`.py:21-25`), and `S_shell`/`S_q` are then extracted from that SAME
postulated closed form by the same operations in both engines. The second engine
therefore re-derives nothing — it echoes the SymPy algebra. Pass-1 batch-5 removed an
X−X gain-line assert but never made the `.wl` independent
(`notes/MATHEMATICA_MIRROR_POLICY.md` records "134 introduced NO new mirror"). The
user has AUTHORIZED a full independent re-author of the `.wl` for this pass.

**Requirement (WHAT must hold — you design the exact route):**
The `.wl` must obtain the mouth-response kernel `S(Π,κ)` from an INDEPENDENT
derivation, NOT by re-typing the SymPy closed form. The natural independent route is
the one stage 133's `.wl` already uses: solve the underlying scalar mixed
Dirichlet/Neumann boundary-value problem symbolically with `DSolveValue` and read the
kernel off the mouth derivative. Specifically, the kernel is defined by

  −u''(x) + κ² u(x) = G·σ(x),   u(0) = 0,   u'(1) = 0,   on x ∈ [0,1],
  with the localized source σ(x) = Π·Exp[−Π·x]/(1 − Exp[−Π]),
  and  S(Π,κ) = u'(0) / G.

From the INDEPENDENTLY-derived kernel:
  1. Obtain the static-shell value `S_shell = S(Π, 0)` (the κ→0 limit) and assert
     `S_shell = 1` exactly — the static shell channel. (Substantive: a wrong BVP
     coefficient breaks it.)
  2. Obtain the canonical `S_q = S(Π, π/2)`.
  3. Assert the DERIVED kernel equals the boxed paper closed form
     `Π[κ·tanh κ + Π(Exp[−Π]·sech κ − 1)] / [(1 − Exp[−Π])(κ² − Π²)]`
     via an `expectZero` residual — this is the cross-engine check that the
     independent derivation reproduces the paper deliverable (can fail if the closed
     form or the derivation is wrong).
  4. Keep the existing fixed-point-law assembly `Π = M_s + M_q·S_q` and the three
     numeric spot-checks of `S_q` at Π = 1/2, 1, 2 against the EXISTING high-precision
     mpmath literal targets (`.wl:57-62`) — those literals are independent of the
     kernel route and must remain.
  5. Keep the `S_q(Π_*)` print and the "no in-stage gain-line assertion" comment
     block (`.wl:74-77`) — the gain line stays printed-only (it is correctly deferred
     to stages 135/137).

After the re-author, `sKernel` may remain ONLY as the boxed-closed-form comparison
target in the `expectZero` of item 3 (clearly labeled as "paper closed form, checked
against the derived kernel"); it must NOT be the source from which `S_shell`/`S_q`
are computed. If you keep a `sKernel`-style symbol, it is the RHS of the check, never
the route.

**Acceptance criteria (how the verifier confirms independence):**
- A `DSolveValue[...]` (or equivalent from-scratch BVP solve) appears and is the
  source of `S_shell` and `S_q`.
- A new `expectZero` PASSes confirming the derived kernel equals the boxed closed form.
- `S_shell = 1` and the three mpmath spot-checks still PASS; the script exits 0.
- The derived-kernel route is textually distinct from the SymPy `.py` (which still
  uses the postulated `S(Pi,kappa)` — leave the `.py` UNCHANGED).

**Self-test reminder:** `κ = π/2` and `κ = 0` are removable points of the kernel
(the only pole is at `κ = Π`, excluded for generic Π > 0), so the limits exist and
equal the closed form — the `expectZero` is a genuine symbolic 0, not vacuous. No new
numeric constant is introduced, so no new `paper_misalignment` risk. The committed
output should change only by ADDED derivation/PASS lines; the deliverable values
(`S_shell=1`, `S_q` numerics, the fixed-point law) are unchanged.

**Verification command:**
The orchestrator runs `redteam exec-mathematica 134` and confirms (a) a `DSolveValue`
BVP solve is present and feeds `S_shell`/`S_q`, (b) the new derived-kernel-equals-boxed
`expectZero` PASSes, (c) `S_shell=1` and the three spot-checks still PASS, (d) the
script exits 0, and (e) the `.py` is unchanged.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl`
  - `mathematica/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.txt`
- summary: Re-authored the Mathematica audit so `S_shell` and `S_q` are computed from a `DSolveValue` mixed Dirichlet/Neumann BVP mouth-derivative kernel, with the boxed paper closed form used only as a residual check target.
- deviation: none
