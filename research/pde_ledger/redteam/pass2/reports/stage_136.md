---
unit_id: 136
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: missing
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: unknown
docs_read:
  paper_stage_tex: present
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md]
  paper_appendix: present
---

# Audit unit 136 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_136.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage136_coupled_mouth_status.md` (only one stage136 notes file present)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: (missing)
- mathematica: (missing)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 136 is a **coupled mouth fixed point and gain selection ledger step** — a status/consolidation card, explicitly tagged `is_status_only_candidate: true` and `is_checkpoint: false` in the manifest. The card's `\stagefield{Verification}` states verbatim: *"SymPy audit: none yet. Mathematica audit: none yet."* It does not introduce a new closed-form constant of its own; it records the *state of knowledge* after the upstream coupled solve (Stages 133–135). The card's bottom-line statement (the in-card quote) is: *"The mouth-layer problem is reduced to a gain pair, or one gain under outlet consistency."* The notes carry the substantive content beyond the terse `.tex`: (1) the exact coupled fixed-point law `Π = M_+ S(Π,κ_+) + M_- S(Π,κ_-)`; (2) the Family-1 shell + mixed D/N reduction `Π = M_s + M_q S_q(Π)` with `κ_s=0, κ_q=π/2`; (3) the canonical gain line `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q`; and (4) the outlet-consistent one-gain law `Π = Σ_m[4 − S_q(Π)]` with canonical value `Σ_m* ≈ 0.451485277739090`. The card's `Checks` and `Downstream use` fields explicitly require that numerical fixed points be recorded as *numerically located, not closed-form constants*, and that the open-status tag be carried with the result downstream (feeds Stages 146–153).

## What the script claims to verify

No script exists for this unit. Manifest entry 136 has `files.sympy.path: null / exists: false`, `files.mathematica.path: null / exists: false`, and both output paths null. A direct directory scan of `scripts/`, `mathematica/`, `scripts/output/`, `mathematica/output/` confirms no `stage136` file (the only substring hits — 193, 030, 059 — are coincidental). This is the legitimate status-only configuration the card itself declares ("none yet"). The audit verdict therefore applies to whether the carry-forward chain the card relies on is actually verified by the upstream units the notes cite.

## Paper ↔ script cross-check

| Paper-side deliverable | Carry-forward source (status-only) | Status |
|---|---|---|
| Coupled fixed-point law `Π = Σ M_α S(Π,κ_α)` | appendix eq:app-part04-general-mouth-fixedpoint; upstream 133–135 scripts | match (carried) |
| Family-1 reduction `Π = M_s + M_q S_q(Π)`, `κ_q=π/2` | appendix eq:app-part04-F1-mouth-fixedpoint; stage134/135 .py+.wl | match (carried) |
| Gain line `M_s ≈ 1.50882951349316 − 0.658075937605429 M_q` | stage134 output L25–26 (`M_s = Pi_star − S_q(Pi_star) M_q`) .py+.wl | match (carried) |
| One-gain law `Π = Σ_m[4 − S_q(Π)]`, `Σ_m* ≈ 0.451485277739090` | stage135 output L5,L18 .py+.wl; stage139 .py+.wl | match (carried) |
| Π* ≈ 1.50882951349316 | stage134/135 .py+.wl (and appendix eq:app-part04-Pi-star) | match (carried) |

Every paper-side deliverable on the 136 card is a carry-forward whose underlying value is independently computed and asserted by an upstream unit's scripts in **both** engines. `paper_alignment: aligned`.

## Assertion inventory

No assertions exist in this unit (no script). Status-only — table not applicable. The carry-forward assertions live in upstream units (134/135/139), which are audited under their own units.

## Findings

None.

### Status-only / missing-engine handling

Per the auditor spec and `is_status_only_candidate: true`: a `missing_sympy`/`missing_mathematica` finding is valid **only** if the unit references a result that no upstream unit's scripts actually verify. That is not the case here. The three load-bearing constants the 136 card/notes report are each produced and asserted by upstream scripts in both SymPy and Mathematica:

- **Π* = 1.50882951349316** — `scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py` and `mathematica/...stage134...wl`; also 135.
- **Gain line `1.50882951349315996... − 0.658075937605428494...·M_q`** — `scripts/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.txt:25-26`, with the `.wl` mirror present.
- **Σ_m* = 0.451485277739089696...** — `scripts/output/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.txt:18`, plus stage135/139 `.wl`.

So `missing_verification_script` does **not** fire. The card honestly self-declares "none yet" and routes verification upstream — consistent with status-only policy.

### First-pass attribution-drift check (133–135 vs 184–186)

First-pass history flagged that 136 had cited downstream "184–186" where it should cite 133–135 (Cluster B). The current state is **corrected**: a grep of both the card and the notes for `184/185/186` returns zero hits. The notes line 6 reads *"After Stages 133–135, the mouth-source selection problem has narrowed again."* — the correct upstream attribution. No residual mis-attribution. No finding.

## Independent-derivation check (Mathematica)

Not applicable — no `.wl` for this unit. The Mathematica mirrors of the carried values exist upstream (stage134/135/139 `.wl`) and are audited there.

## Engine cross-check

Not applicable — no engines present in this unit. Upstream engine agreement (e.g., 135's `Sigma_m^* = 0.451485277739089696...` matching across `.py` and `.wl` outputs) is the responsibility of those units' audits; spot-confirmed consistent.

## Verdict justification

`clean`. I read the card, the single notes file, and the relevant appendix section (`subsec:app-part04-coupled-mouth-layer-fixed-point`, eqs Pi-star / F1-mouth-fixedpoint / one-gain-mouth-law / Sigmam-star). The card is a legitimate status-only consolidation ledger entry that explicitly declares "SymPy/Mathematica audit: none yet" and carries forward results from Stages 133–135. Attacks attempted and failed: (1) tried to find an unsupported carry-forward — all three load-bearing constants are independently asserted by upstream scripts in both engines, with exact digit agreement; (2) tried to find the previously-flagged 184–186 mis-attribution — it is corrected to 133–135 with zero residual hits; (3) tried to find a paper↔value mismatch — Π*, the gain-line coefficients, and Σ_m* all match the appendix and the upstream script outputs to all reported digits; (4) checked for a hidden script the manifest missed — none exists. The missing engines are not a finding under status-only policy because the referenced results are all upstream-verified.

## Self-test notes

Trap checks run: (1) variable-independence / derivative-of-constant traps — n/a, no script and no in-unit derivatives. (2) Parity/symmetry of integrals — n/a, no integrals in this unit. (3) Trivial-case pre-check — applied as digit-level reconciliation instead: card values reproduce upstream output digits exactly (Π* 1.50882951349316 vs stage134; Σ_m* 0.451485277739090 vs stage135 L18; slope 0.658075937605429 vs stage134 L26). (4) Path specs — confirmed no `stage136` file exists in `scripts/`, `mathematica/`, or their `output/` dirs, so no `missing_verification_script` directive is warranted. (5) Paper round-trip — no fix prescribed, so no new misalignment introduced. Conclusion: status-only card with a fully upstream-verified carry-forward chain; no directive written.

## Value Reconciliation (pass-2 augmentation)

This is a status-only unit with **no script of its own**; therefore there are no script-emitted result values local to unit 136. The "values" on the 136 card/notes are carry-forwards whose authoritative emitters are the upstream stage 134/135/139 scripts (both engines). I reconcile each carried deliverable value against (a) where it is emitted upstream and (b) where it lives in the 136 docs and the Part IV appendix. Per the augmentation's status-only clause, I note the no-local-script status explicitly and base the reconciliation on the script *outputs* of the carry-forward sources (read, not run).

| value | source (upstream py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Π* ≈ 1.50882951349316 | stage134 `.py`/`.wl` (Pi_star); stage135 | notes L28 (within gain line); appendix L663 `\Pi_*\approx1.50882951349316`; card body (carried) | MATCH |
| gain-line intercept 1.50882951349315996… (= Π*) | `scripts/output/...stage134...sympy_audit.txt:26` | notes L28 `1.50882951349316`; appendix gain context | MATCH |
| gain-line slope 0.658075937605428494… (= S_q(Π*)) | `scripts/output/...stage134...sympy_audit.txt:26` | notes L28 `0.658075937605429` | MATCH |
| Σ_m* ≈ 0.451485277739090 | `scripts/output/...stage135...sympy_audit.txt:18` (0.451485277739089696…); stage139 `.wl` | notes L41 `\Sigma_m^* \approx 0.451485277739090`; appendix L724 | MATCH |
| Family-1 angles κ_s=0, κ_q=π/2 | stage134/135 scripts (S_q := S(Π,π/2)) | notes L22-24; appendix eq F1-mouth-fixedpoint L706 | MATCH |
| 4:−1 weighting → `Π = Σ_m[4 − S_q(Π)]` | stage135 output L5 (`M_s + M_q*S_q - Sigma_m*(4 - S_q) = 0`) | notes L35-38; appendix eq one-gain-mouth-law L718 | MATCH |

INTERNAL / scaffolding (no finding, none local to this unit since there is no script): n/a — there are no local pass/fail flags, residuals, or tolerances to exclude, because no script exists for unit 136.

reconciliation: complete; 6 carried deliverable values checked, 0 misaligned. Single-engine/no-script status: this unit has no script in either engine by design (status-only consolidation card, manifest `is_status_only_candidate: true`); all reconciled values are carry-forwards verified upstream by Stages 134/135/139 in BOTH SymPy and Mathematica, and each agrees to all digits the 136 docs report.
