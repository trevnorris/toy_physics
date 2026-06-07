---
unit_id: 134
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T22:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 134

## Per-finding outcomes

### F1 — mathematica_transliteration (USER-AUTHORIZED full re-author of the `.wl`)

**Classification:** resolved

**What changed:**
Codex re-authored `mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl` so the mouth-response kernel is DERIVED, not transliterated:

- `.wl:33` defines the localized source `sigma = p*Exp[-p*x]/(1 - Exp[-p])`.
- `.wl:36-40` solves the scalar mixed Dirichlet/Neumann BVP from scratch:
  `uSol = DSolveValue[{-uFun''[x] + k^2*uFun[x] == G*sigma, uFun[0] == 0, uFun'[1] == 0}, uFun[x], x]`.
- `.wl:46` reads the kernel off the mouth derivative:
  `derivedKernel = FullSimplify[(D[u, x] /. x -> 0)/G, ...]` — i.e. `S(Π,κ) = u'(0)/G`, exactly the directive's spec.
- `.wl:61` `sShell` is the `kk -> 0` limit of `derivedKernel /. k -> kk`; `.wl:64` `sQ = derivedKernel /. k -> Pi/2`. Both `S_shell` and `S_q` are now sourced from `derivedKernel`, NOT from the old `sKernel`/`paperKernel` closed form.
- `.wl:49-52` retains the boxed paper closed form as `paperKernel` (renamed from `sKernel`), explicitly commented "Paper closed form, checked against the derived kernel; not used as the source."
- `.wl:54-57` adds four new `expectZero` checks: BVP ODE residual, the two boundary conditions u(0) and u'(1), and the cross-engine check `derivedKernel - paperKernel`.

**Assessment:**
The edit fully resolves the finding and is correct.

(1) Independence — confirmed. `DSolveValue` (`.wl:36`) is the genuine source of `derivedKernel`, and `S_shell`/`S_q` are computed exclusively from `derivedKernel` (`.wl:61,64`). The only surviving appearance of the SymPy-form closed expression (`paperKernel`, `.wl:49-52`) is as the RHS of the `expectZero` at `.wl:57` — it is the comparison target, never the route. The `.wl` is now textually distinct from the `.py`, which still builds `S(Pi,kappa)` from the postulated closed form (`.py:21-25`) and substitutes `pi/2` directly (`.py:30`). Two structurally different routes now arrive at the same kernel — exactly what the dual-engine policy requires. F1 is genuinely resolved, not papered over.

(2) Non-tautology of the new check — the `derivedKernel - paperKernel` `expectZero` is substantive: `derivedKernel` is produced independently by `DSolveValue` + mouth-derivative extraction, while `paperKernel` is the hand-typed boxed form; a sign/factor error in either side fails the residual. The three BVP sanity checks (ODE residual, u(0), u'(1)) further confirm the solved `u` actually satisfies the stated problem before the kernel is read off, so the derivation cannot silently produce a wrong-but-self-consistent kernel. Output shows all four PASS with residual `= 0`.

(3) Preserved deliverables — `S_shell = 1` PASSes (`static shell channel = 0`, PASS), and the three mpmath spot-checks at p=1/2,1,2 still PASS against the unchanged 30-digit literals (`.wl:81-86`). The fixed-point law `Ms + Mq*sQ`, the `S_q(Pi_star)` print, and the print-only/non-asserted gain-line block with its X-X-avoidance comment (`.wl:95-101`) are all retained per directive items 4-5.

No collateral edits beyond the authorized re-author: the only other change is the `expectZero` helper now wrapping `Together[Expand[...]]` (`.wl:21`) instead of bare `Expand[...]` — a necessary and benign simplification-robustness tweak so the rational derived-vs-paper residual collapses to an exact `0` (it does; output shows `= 0`/PASS). This is in-scope support for the new symbolic residual check, not an unrelated refactor.

## Exec log assessment

**SymPy:** exit=0. The `.py` was not edited, so this is unchanged-baseline corroboration. Notable lines:
- `S_shell = 1`
- `OK: S_q matches independent numeric targets at Pi=1/2, 1, 2`
- `OK: S_q(Pi_star) matches notes value 0.658075937605428`

**Mathematica:** exit=0. Notable lines confirming the re-author:
- `BVP ODE residual = 0` / `PASS: BVP ODE residual` (plus `BVP u(0)`, `BVP u'(1)` PASS) — the DSolveValue solution genuinely satisfies the stated BVP.
- `derived mouth kernel equals boxed paper closed form = 0` / `PASS: derived mouth kernel equals boxed paper closed form` — the independent route reproduces the paper deliverable.
- `S_shell = 1` / `static shell channel = 0` / `PASS: static shell channel`.
- `PASS: S_q at p=1/2`, `PASS: S_q at p=1`, `PASS: S_q at p=2` (diffs ~1e-30 against the unchanged literals).

The Mathematica `S_q(p)` is now printed in the DSolveValue-derived algebraic form `(8*E^(Pi/2)*p^2 - 2*E^p*p*(2*(1 + E^Pi)*p + Pi - E^Pi*Pi))/((-1 + E^p)*(1 + E^Pi)*(-4*p^2 + Pi^2))`, which is visibly distinct from (yet provably equal to, via the passing expectZero) the SymPy `cosh`/`tanh` boxed form — direct textual evidence that the kernel was re-derived rather than retyped.

**Output freshness:** confirmed. `.wl` mtime 2026-06-06 17:03:28; its `.txt` output mtime 2026-06-06 21:47:39 (output newer than script → regenerated post-fix by the orchestrator's independent re-run). The `.py` mtime is 2026-05-29 16:41:15 (unchanged), and its `.txt` output mtime is 2026-06-06 21:47:39 (refreshed in the same re-run). `git diff HEAD` and `git status` for the `.py` are both empty → the SymPy script is UNCHANGED, satisfying directive criterion (e).

## Material-change assessment

`material_change`: false.

No deliverable VALUE changed. `S_shell = 1`, `S_q` at the three spot points and at `Pi_star ≈ 0.658075937605428`, the fixed-point law `Π = M_s + M_q S_q(Π)`, and the canonical gain line are all bit-for-bit identical to the pre-fix outputs. The committed Mathematica output changed only by (a) added derivation/BVP/PASS lines and (b) the kernel now being printed in its DSolveValue-derived algebraic form (same mathematical object, different but provably-equal surface syntax — confirmed equal by the passing `expectZero`). No downstream unit can observe a changed numeric. The orchestrator's blanket `upstream_stale` marking for units > 134 is procedurally fine, but there is no substantive concern here: re-audit of downstream stages on account of unit 134 is not warranted by any value change.

## Side observations (non-blocking)

- The `$Assumptions` are deliberately switched twice: the BVP/derivation block (`.wl:29-31`) carries `k != p` so the kernel's only pole (κ=Π) is excluded; then `.wl:59` resets to `p > 0 && kk > 0` for the limit/substitution block. This is correct and necessary (κ=π/2 and κ=0 are removable points, so the limit and substitution are genuine, matching the directive's self-test note). Not a problem — noted only for completeness.
- `paperKernel` retains the symbol legacy of the old `sKernel`; it is now correctly demoted to a check target and clearly commented as such. No action needed.

## Verdict justification

The single finding (F1, the only finding in the original report) is fully resolved by the USER-AUTHORIZED re-author. The `.wl` now derives the mouth kernel independently via `DSolveValue` on the scalar mixed D/N BVP and reads `S(Π,κ) = u'(0)/G` off the solution; `S_shell` and `S_q` are sourced from that derived kernel, with the boxed paper closed form used only as the RHS of a passing cross-check `expectZero`. All four new checks plus the static shell limit and the three mpmath spot-checks PASS, both engines exit 0, the outputs were refreshed post-fix, the `.py` is unchanged, and no deliverable value moved. F1 is genuinely resolved, no regressions appear in the diff or logs, and no material change propagates downstream → `verdict: verified`, `material_change: false`.
