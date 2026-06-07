---
unit_id: 132
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md]
  paper_appendix: present
---

# Audit unit 132 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_132.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage132_mouth_boundary_layer_status.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (relevant rows: line 29 "Stages 125--132"; the MTDC-T8.4 block, esp. `\subsection{Electrochemical boundary layer}` lines 630-672; result-anchor row line 1178; `\input{stages/stage_132}` line 1298)
- sympy: (missing — status-only, no script of its own)
- mathematica: (missing — status-only, no script of its own)
- sympy output: (missing)
- mathematica output: (missing)

## What the paper claims

Stage 132 is an explicit **status-only consolidation card** ("a positive mouth-source and boundary layer ledger step"; `\stagefield{Verification}{SymPy audit: none yet. Mathematica audit: none yet.}`). Its `\stagefield{Purpose}` says the audit target is "the verification output quoted below," and the quoted body claim is: "The source-shape ambiguity is reduced to the selected bias \(\Pi_m\)." The card carries forward, without re-deriving, the lower compensated core branch, the positivity of the local mouth source, and the explicit GNLS/localized-Maxwell boundary-layer reduction. The notes file enumerate the concrete deliverables it consolidates from Stages 129–131: (1) the boundary-layer source law \(\sigma_\Pi(z)=\Pi e^{-\Pi z/L}/(L(1-e^{-\Pi}))\); (2) the exact Family-1 mouth-bias factor \(\mathfrak g_\Pi = 2\Pi(2\Pi e^\Pi+\pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\); (3) its strict monotone range \(2/\pi\to 1\) as \(\Pi:0^+\to+\infty\); (4) the unique canonical compensation point \(\Pi_*\approx 1.50882951349316\); and (5) the parent threshold \(\partial_z\delta V_{\rm conf}|_m - q_*\partial_z A_0|_m = 1.50882951349316\,\Theta_\sigma/L\). The card explicitly tags downstream use ("feeds Stages 133--145"; "the status tag above must be carried with the result").

## What the script claims to verify

There is no SymPy or Mathematica script for unit 132, which is consistent with the manifest entry (`is_status_only_candidate: true`, `is_checkpoint: false`, both `files.sympy` and `files.mathematica` are `path: null, exists: false`). The deliverables it consolidates are independently verified by the upstream scripts, which I read to confirm the carry-forward is genuinely supported:
- Stage 129 (`scripts/moving_throat_pde_stage129_mouth_boundary_layer_sympy_audit.py`): asserts the source profile \(\sigma_\Pi\) is normalized on [0,L], that the zero-flux current vanishes on the stationary branch, and that the boundary-layer ODE residual is zero (lines 26-31). Mathematica mirror present.
- Stage 130 (`...stage130_mouth_bias_map_sympy_audit.py`): derives \(\mathfrak g_\Pi\) by integrating \(\sigma_\Pi\cos(\pi z/2L)\) and asserts it equals the boxed closed form (lines 15-18); proves the limits \(2/\pi\) and \(1\) (lines 33-36); supplies a GLOBAL strict-monotonicity certificate via the symmetrized covariance identity (lines 53-87); and solves \(\mathfrak g_\Pi=g_-\) for \(\Pi_*\) (lines 102-103). Mathematica mirror present with an independent FindRoot route.
- Stage 131 (`...stage131_parent_mouth_threshold_sympy_audit.py`): asserts \(\Pi_*\) matches the literal `1.50882951349316` to `1e-14` (lines 40-45), the slope \(g'(\Pi_*)\) matches `0.0714453558083195` (lines 47-53), and the parent threshold form \((T_m - q_* A0p)-\Pi\Theta_\sigma/L\) (lines 26-35). Mathematica mirror present, independent root isolation (lines 40-62).

## Paper ↔ script cross-check

Because 132 owns no script, the cross-check is between the card/notes/appendix deliverables and the upstream verifying scripts that the carry-forward depends on.

| Paper-side deliverable (132 notes/appendix) | Verifying engine(s) upstream | Status |
|---|---|---|
| \(\sigma_\Pi(z)=\Pi e^{-\Pi z/L}/(L(1-e^{-\Pi}))\) | sympy+wl 129 (normalization, zero-flux, ODE); 130 line 12 | match |
| \(\mathfrak g_\Pi = 2\Pi(2\Pi e^\Pi+\pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\) | sympy 130:15-18, wl 130:42-44; sympy 131:23, wl 131:40 | match |
| monotone \(2/\pi\to 1\) | sympy 130:33-36 + global certificate 130:53-87; wl 130 mirror | match |
| \(\Pi_*\approx 1.50882951349316\) | sympy 130:102, sympy 131:40-45, wl 130:110, wl 131:62 | match |
| parent threshold \(=\Pi_*\Theta_\sigma/L\) | sympy 131:26-35, wl 131 mirror | match |
| body claim "ambiguity reduced to selected bias \(\Pi_m\)" | qualitative consolidation, no numeric assertion required | match (status-only) |

`paper_alignment: aligned`. Every consolidated value is present and agreeing in the card-supporting notes and the part-04 appendix, and every one is exercised by a non-tautological assertion in at least one upstream engine (most by both).

## Assertion inventory

No assertions exist in this unit (no script). The table below records the upstream assertions that back the carry-forward, for traceability only; they are not 132's own checks.

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| (none) | n/a | — | — | unit has no script | n/a |
| U1 | sympy 129 | 26-31 | `assert simplify(...)==0` (norm, J, ODE) | deliverable 1 (\(\sigma_\Pi\)) | yes |
| U2 | sympy 130 | 17-18 | `g_Pi - boxed == 0` | deliverable 2 (\(\mathfrak g_\Pi\)) | yes |
| U3 | sympy 130 | 33-36 | limits == 2/pi, == 1 | deliverable 3 (range) | yes |
| U4 | sympy 130 | 64-86 | symmetrized Cov + dg/dPi sign | deliverable 3 (monotonicity) | yes |
| U5 | sympy 131 | 42-45 | `\|Pi_*-1.50882951349316\|<1e-14` | deliverable 4 (\(\Pi_*\)) | yes |
| U6 | sympy 131 | 34 | threshold residual form | deliverable 5 | yes |

## Findings

None.

## Independent-derivation check (Mathematica)

Not applicable to unit 132 (no `.wl` of its own). I did spot-check the upstream mirrors that back the carry-forward: stage 131's `.wl` explicitly comments "INDEPENDENT Pi_* route (not a transliteration of SymPy's nsolve)" (line 41) and isolates the root by bracketing rather than echoing SymPy's `nsolve`; stage 130's `.wl` re-derives \(\mathfrak g_\Pi\) via `Integrate` then independently `FindRoot`s. These are not in scope for 132's verdict but confirm the consolidated values rest on dual-engine verification upstream.

## Engine cross-check

Not applicable — unit 132 has no engines. (Upstream 130/131 carry both engines; both pin \(\Pi_* = 1.50882951349316\) to `1e-14`, so the carried value is dual-engine-consistent.)

## Verdict justification

`clean`. I read the paper card, the notes file, and the relevant part-04 appendix block before opening any script, and I confirmed by directory search (`find ... -iname "*stage132*"`) and the pass-2 MANIFEST that unit 132 legitimately owns no SymPy or Mathematica script and is correctly tagged `is_status_only_candidate: true`. Per the status-only handling rule, a missing engine is therefore NOT a finding: every result the card consolidates is referenced from Stages 129–131, and each of those stages' scripts (both engines) actually verifies the referenced result with non-tautological assertions — I checked the source of each rather than trusting comments. Attacks I tried that failed: (a) hunting for a hidden/renamed 132 script under `scripts/` and `mathematica/` (none exists); (b) checking whether the consolidated \(\Pi_*\), \(\mathfrak g_\Pi\), \(\sigma_\Pi\), the \(2/\pi\to 1\) range, and the parent threshold are actually proven upstream or merely asserted in prose (all proven, with explicit tolerance/closed-form/limit checks); (c) checking the first-pass notes-attribution drift — the prior pass flagged 132 citing downstream "180–182" where it should cite 129–131; the current notes line 7 reads "After Stages 129–131" and a `grep` for `180/181/182` over the card and notes returns nothing, so the corrected attribution holds and no residual mis-attribution remains; (d) checking the card body claim against the appendix — the card's "ambiguity reduced to the selected bias \(\Pi_m\)" matches the appendix MTDC-T8.4 narrative and the boxed open question in the notes. Nothing broke.

## Value Reconciliation (pass-2 augmentation)

Unit 132 emits NOTHING of its own (status-only, no script, no saved output). The reconciliation below accounts for each deliverable value the card/notes/appendix REPORT for this stage and confirms each (i) is reflected consistently across the `.tex` appendix and the `.md` notes, and (ii) is actually emitted/verified by the upstream script(s) the carry-forward depends on. Per the augmentation's status-only clause, I reconcile what the present (upstream) engines emit and note the no-script status. I did not run anything; I used the script source plus my reading of the closed forms (no committed 132 output exists).

| value | source (upstream py/wl) | .tex/.md location | status |
|---|---|---|---|
| \(\sigma_\Pi(z)=\Pi e^{-\Pi z/L}/(L(1-e^{-\Pi}))\) | sympy 129:13, 130:12; wl 130:42 (defn) | notes lines 9-12; appendix eq:app-part04-sigmaPi `stage_appendix_part04.tex:646-652` | MATCH |
| \(\mathfrak g_\Pi = 2\Pi(2\Pi e^\Pi+\pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\) | sympy 130:16 (boxed, asserted ==); sympy 131:23; wl 130:43, wl 131:40 | notes lines 14-20; appendix eq:app-part04-gPi `:654-659` | MATCH |
| monotone range \(2/\pi\to 1\) | sympy 130:33-36 (limits asserted); wl 130:54-55 | notes lines 22-26; appendix `:660` "strictly increasing from \(2/\pi\) to \(1\)" | MATCH |
| \(\Pi_*\approx 1.50882951349316\) | sympy 130:102-103, sympy 131:44 (literal, asserted to 1e-14); wl 130:110, wl 131:62 | notes lines 28-31; appendix eq:app-part04-Pi-star `:661-664` | MATCH |
| parent threshold \(=1.50882951349316\,\Theta_\sigma/L\) | sympy 131:26-35; wl 131 mirror | notes lines 33-40; appendix eq:app-part04-parent-mouth-threshold `:665-672` | MATCH |
| selected-bias target \(\Pi_m\approx 1.51\) (open-question framing) | (downstream/qualitative; not a script result) | notes lines 54-68 (boxed open question); appendix MTDC-T8.4 narrative | MATCH (carried qualitatively; consistent with \(\Pi_*\)) |

INTERNAL / scaffolding (no prose reflection expected, raise no finding): none originate in unit 132 (it has no script). The upstream slope `g'(Pi_*) = 0.0714453558083195` and singular-branch separation `Delta g_- = 0.241964921055337` are stage-131 internal anchors, not stage-132 deliverables, and are reconciled in stage 131's own report, not here.

reconciliation: complete; 6 deliverable values checked, 0 misaligned (status-only: no 132-native emissions; all values reconcile against upstream dual-engine scripts + notes + appendix).

## Self-test notes

I checked: (1) variable-independence / spurious-zero-derivative traps — n/a, no new derivatives prescribed (no directive, no script); I did confirm the upstream `sp.diff(gPi, Pi)` is non-trivial since `gPi` genuinely depends on `Pi`. (2) Hidden/renamed-script trap — ran a filename `find` over both script dirs to be sure the "missing" status is real, not a path-mismatch; confirmed absent. (3) Stale-attribution trap — grep-confirmed no residual "180–182" in card/notes and that "129–131" is the live attribution. Conclusion: status-only verdict is sound, every consolidated value is dual-engine-backed upstream and matches the notes and the part-04 appendix; no directive is warranted (zero findings).
