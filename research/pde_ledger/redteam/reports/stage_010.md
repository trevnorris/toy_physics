---
unit_id: 010
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 010 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex`
- notes: (none) — `notes/em_projected/step_NN_*.md` was not committed for the EM-projected stages (004-020); the orchestrator's rendered prompt confirms `notes: (none)`.
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row for stage 010 reads: "Projected Maxwell push into bundle slots / \StatusExactClosure{} / Projected Z_n, N_n slot transport for the grouped response bundle.")
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.txt`

## What the paper claims

The card (`paper/stages/stage_010.tex`) titles the stage "Projected Maxwell push into bundle slots" and labels it the "master transport from projected EM corrections into the grouped bundle slots." It declares the shift convention

> `Z_n \mapsto Z_n + \varepsilon z_n,  N_n \mapsto N_n + \varepsilon n_n,  n \in \{0, 2, 4\}`

with `D_0 = K - B_0 - Z_0`, `D_2 = -(M + B_2 + Z_2)`, `D_4 = -(B_4 + Z_4)`, and then states three explicit first-variation identities (the "transport map"):

1. (eq:stage010-du2) `\delta u_2 = (D_0 z_2 - D_2 z_0)/D_0^2`
2. (eq:stage010-du4) `\delta u_4 = (D_0^2 z_4 - D_0 (2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0)/D_0^3`
3. (eq:stage010-dp0) `\delta P_0 = (D_0 n_0 + N_0 z_0)/D_0^2`

The `\stagefield{Output}` paragraph reads verbatim: "Stage 010 exports the transport map (eq:stage006-projected-shifts)--(eq:stage010-dp0)." So the deliverables are the shift convention plus the three first-variation formulas. Note in particular that `δu_2` and `δu_4` depend on the `z_n` shifts ONLY — they have no `n_n` dependence, marking them as D-side (denominator-only) inversion quantities. By contrast `δP_0` carries an `n_0` contribution and is a full numerator/denominator ratio variation. The part-appendix row reinforces "Projected Z_n, N_n slot transport for the grouped response bundle." but adds no further algebraic detail.

No notes file exists in the repo for the EM-projected stages (004-020). The .tex is therefore the sole paper-side authority here.

## What the script claims to verify

Both scripts (sympy `:38-178` and mathematica `:7-223`) construct perturbed denominators `D_n^p = D_n - \varepsilon z_n` and perturbed numerators `N_n^p = N_n + \varepsilon n_n` for `n \in \{0,2,4\}`, then build three quotients `P_0(\varepsilon) = N_0^p/D_0^p`, `P_2(\varepsilon) = (D_0^p N_2^p - 2 D_2^p N_0^p)/(D_0^p)^2`, and `P_4(\varepsilon) = ((D_0^p)^2 N_4^p - 2 D_0^p (D_2^p N_2^p + D_4^p N_0^p) + 3 (D_2^p)^2 N_0^p)/(D_0^p)^3`. The first-order \varepsilon coefficients `dP_n` are then asserted to match explicit closed forms (sympy `:58-83`, mathematica `:39-63`). Beyond these three, the scripts verify a large secondary structure: one-pole and normalization K-surfaces from a relation `(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0)(T + \varepsilon z_4) = 3(S + \varepsilon z_2)^2` plus `(N_0 + \varepsilon n_0)/(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0) = P_{\mathrm{target}}`; a "transported-target" variant `P_{\mathrm{target,transport}} = (N_0 + \varepsilon n_0)/D_{0,\mathrm{target}}`; first-variation `dcompat`, `dK_one_pole`, `dK_norm`, `dcompat_transport` identities with explicit z_0-cancellation guards; two `assert_nonzero` sign-flip mutation guards; real-Y_{20}-square Gaunt overlap ratios for m=0,1,2 producing lane multipliers (1, 1/2, -1); a grouped trace/anomaly decomposition (xbar, ax, bx) along a "weak-axisymmetric" line b = 3a; and a primitive static Xi formula relating two ratio-perturbations under a specific `N_{0,\mathrm{sym}} = P^2/\Delta^2` substitution.

## Paper ↔ script cross-check

| # | Paper deliverable | Script check | Status |
|---|---|---|---|
| C1 | Shift convention `Z_n -> Z_n + \varepsilon z_n`, `N_n -> N_n + \varepsilon n_n` (paper eq:stage006-projected-shifts) | sympy `:43-48`, mathematica `:24-29` build `D_0^p = D_0 - \varepsilon z_0` etc., `N_0^p = N_0 + \varepsilon n_0` etc. | match |
| C2 | `\delta u_2 = (D_0 z_2 - D_2 z_0)/D_0^2` (paper eq:stage010-du2) | no assertion of this identity anywhere. Closest is sympy `:59-68` / mathematica `:48-53` on `dP_2`, which contains additional `n_0, n_2, N_0, N_2` terms and is a different object. | mismatch |
| C3 | `\delta u_4 = (D_0^2 z_4 - D_0 (2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0)/D_0^3` (paper eq:stage010-du4) | no assertion of this identity. Closest is sympy `:69-83` / mathematica `:55-63` on `dP_4`, again a different object including `N_n` and `n_n`. | mismatch |
| C4 | `\delta P_0 = (D_0 n_0 + N_0 z_0)/D_0^2` (paper eq:stage010-dp0) | sympy `:58` `assert_zero("delta P0", dP0 - (n0/D0 + N0*z0/D0**2))`; mathematica `:44-46` `m1Residual = ... - (n0/D0 + N0 z0/D0^2)`. Algebraically identical to paper form (multiply paper RHS top and bottom by 1; `(D_0 n_0 + N_0 z_0)/D_0^2 = n_0/D_0 + N_0 z_0/D_0^2`). | match |
| C5 | (none) | sympy `:85-139` and mathematica `:65-143` build K-surfaces, compatibility, transported-target compatibility, dK shifts, z_0-cancellation guards, and `assert_nonzero` mutation guards — none of this appears in the paper card. | extra |
| C6 | (none) | sympy `:141-155` and mathematica `:145-189` verify Gaunt-based lane multipliers (1, 1/2, -1) and a trace/anomaly grouping along b = 3a — paper card does not introduce `Y_{20}`, Gaunt overlaps, or a grouped trace decomposition. | extra |
| C7 | (none) | sympy `:157-169` and mathematica `:191-199` verify the "primitive static Xi" formula. Paper card does not introduce Xi, Q, Delta, P, q_1, p_1, d_1, or this substitution. | extra |

Pattern: ONE of three paper-side identities is exercised exactly; TWO are not exercised at all (the scripts verify a related-but-different object); and SEVEN orthogonal claim clusters appear in the scripts that the paper card never mentions. Setting `paper_alignment: partial` reflects this split.

Important context: the script's `δP_n` formulas (sympy `:50-52`, mathematica `:31-37`) are not the paper's `δu_n` even under any obvious specialization. Setting `n_0 = n_2 = 0` in `dP_2` (which would suppress the n-shifts) leaves `N_2 z_0/D_0^2 + 2 N_0 z_2/D_0^2 - 4 D_2 N_0 z_0/D_0^3` — still depends on `N_0, N_2`, whereas `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2` depends only on D's and z's. Hand-derivation of the paper's `δu_2` matches `u_2 = -D_2/D_0` under `D_n -> D_n - \varepsilon z_n` (no N involvement); `δu_4` matches `u_4 = (D_2^2 - D_0 D_4)/D_0^2`. The paper's u_n objects are pure D-side (denominator-polynomial-inversion) coefficients, distinct from the script's full-ratio P_n.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 58 | `assert_zero("delta P0", dP0 - (n0/D0 + N0*z0/D0**2))` | C4 (δP_0) | yes |
| A2 | sympy | 59-68 | `assert_zero("delta P2", dP2 - (...))` | none (script's P_2 differs from paper's u_2; C2 not exercised) | partial (script-internal only) |
| A3 | sympy | 69-83 | `assert_zero("delta P4", dP4 - (...))` | none (script's P_4 differs from paper's u_4; C3 not exercised) | partial (script-internal only) |
| A4 | sympy | 114 | `assert_zero("compatibility surface after eliminating K", ...)` | none | extra |
| A5 | sympy | 115 | `assert_zero("one-pole K shift", ...)` | none | extra |
| A6 | sympy | 116 | `assert_zero("normalization K shift", ...)` | none | extra |
| A7 | sympy | 117 | `assert_zero("compatibility first variation from competing K surfaces", ...)` | none | extra |
| A8 | sympy | 118 | `assert_zero("compatibility first variation from eliminated surface", ...)` | none | extra |
| A9 | sympy | 119 | `assert_zero("compatibility first variation", ...)` | none | extra |
| A10 | sympy | 120 | `assert_zero("transported-target normalization K surface", ...)` | none | extra |
| A11 | sympy | 121-124 | `assert_zero("transported-target compatibility surface", ...)` | none | extra |
| A12 | sympy | 125-128 | `assert_zero("transported-target compatibility first variation", ...)` | none | extra |
| A13 | sympy | 129 | `assert_zero("fixed-target compatibility z0 cancellation", ...)` | none | extra |
| A14 | sympy | 130 | `assert_zero("transported-target compatibility z0 cancellation", ...)` | none | extra |
| A15 | sympy | 131 | `assert_nonzero("normalization K surface still carries z0", ...)` | none | extra |
| A16 | sympy | 132-135 | `assert_nonzero("mutated compatibility transport should fail", ...)` | none | extra (mutation guard) |
| A17 | sympy | 136-139 | `assert_nonzero("mutated transported-target compatibility should fail", ...)` | none | extra (mutation guard) |
| A18 | sympy | 146 | `assert_zero("Y20 overlap lane 20", lam20 - 1)` | none | extra |
| A19 | sympy | 147 | `assert_zero("Y20 overlap lane 21", lam21 - 1/2)` | none | extra |
| A20 | sympy | 148 | `assert_zero("Y20 overlap lane 22", lam22 + 1)` | none | extra |
| A21 | sympy | 152 | `assert_zero("weak-axisymmetric trace", xbar - x0)` | none | extra |
| A22 | sympy | 153 | `assert_zero("weak-axisymmetric a anomaly", ...)` | none | extra |
| A23 | sympy | 154 | `assert_zero("weak-axisymmetric b anomaly", ...)` | none | extra |
| A24 | sympy | 155 | `assert_zero("weak-axisymmetric line b=3a", bx - 3*ax)` | none | extra |
| A25 | sympy | 165-169 | `assert_zero("primitive static Xi formula", ...)` | none | extra |
| B1 | mathematica | 44-46 | `If[FullSimplify[slot0Linear - (n0/D0 + N0 z0/D0^2)] =!= 0, ...]` | C4 (δP_0) | yes |
| B2 | mathematica | 48-53 | `... slot2Linear - m2Target ...` | none (mirror of A2) | partial |
| B3 | mathematica | 55-63 | `... slot4Linear - m3Target ...` | none (mirror of A3) | partial |
| B4-B17 | mathematica | 65-220 | K-surfaces, compat, Gaunt, trace, Xi, mutations | none (mirror of A4-A25) | extra |

ZERO script-side assertions exercise paper claims C2 (`δu_2`) or C3 (`δu_4`). One assertion pair (A1 / B1) exercises C4 (`δP_0`). Everything else is either related-but-different or unrelated extra content.

## Findings

### F1 — paper_misalignment

**Severity:** high
**Subtype:** script_missing_paper_claim (two missing) and target_mismatch (script verifies adjacent objects)
**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex:25-33`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:50-83`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl:31-63`

**What's wrong:**

Paper side (`paper/stages/stage_010.tex:25-33`):
```
\delta u_2 = (D_0 z_2 - D_2 z_0)/D_0^2,
\delta u_4 = (D_0^2 z_4 - D_0 (2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0)/D_0^3.
```
and the `\stagefield{Output}` line declares "Stage 010 exports the transport map (eq:stage006-projected-shifts)--(eq:stage010-dp0)", which by the equation labels means equations `eq:stage006-projected-shifts`, `eq:stage010-du2`, `eq:stage010-du4`, and `eq:stage010-dp0` are all part of the export.

Script side (sympy `:50-52`, mathematica `:31-37`):
```python
    P0p = N0p / D0p
    P2p = (D0p * N2p - 2 * D2p * N0p) / D0p**2
    P4p = (D0p**2 * N4p - 2 * D0p * (D2p * N2p + D4p * N0p) + 3 * D2p**2 * N0p) / D0p**3
```
and the only first-variation closed-form assertion that matches a paper RHS is `dP0` (sympy `:58`, mathematica `:44-46`).

The paper's `u_2`, `u_4` are pure D-side inversion coefficients: setting `u_2 = -D_2/D_0` reproduces the paper's `δu_2` under `D_n -> D_n - \varepsilon z_n`, and `u_4 = (D_2^2 - D_0 D_4)/D_0^2` reproduces the paper's `δu_4` identically. The script's `P_2`, `P_4` are full ratios of N and D pieces (the same structures that appear in higher-order ratio expansions); their first variations are not the paper's `δu_2`, `δu_4` even in the `n_0 = n_2 = n_4 = 0` limit (the N's stay in `δP_n`).

The script verifies one of three paper deliverables and substitutes adjacent-but-distinct content for the other two. Neither `δu_2` nor `δu_4` appears as an assertion target anywhere in the sympy or mathematica scripts.

**Why this matters:**

The paper card's `\stagefield{Output}` line and equation labels promise three first-variation identities. The script certifies one. A downstream consumer of stage 010 (stage 011, 012, ...) reading the paper would expect the algebra `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2` to be machine-checked; it is not. If `δu_2`, `δu_4` are physically the correct paper-side objects, the scripts must be augmented. If `δP_2`, `δP_4` are the correct objects, the paper card must be updated (and the equation labels `eq:stage010-du2`, `eq:stage010-du4` renamed and rewritten with the n-shifts visible). Either way, the divergence is the user's call — Codex must not silently rewrite either side.

**Required change:**
See `## Resolve before fix_loop` block in `/var/projects/toy_physics/research/pde_ledger/redteam/directives/stage_010.md`. Routed to user.

**Verification:**
After the user picks a direction:
- Direction (a) "paper is correct": the scripts must add two new `assert_zero` blocks for `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2` and `δu_4 = (D_0^2 z_4 - D_0(2 D_2 z_2 + D_4 z_0) + 2 D_2^2 z_0)/D_0^3` with appropriate definitions of u_2, u_4 (or by directly asserting the closed forms against `sp.diff` of `(-D_2 + eps z_2)/(D_0 - eps z_0)` etc.). The existing `δP_n` checks may stay if useful as cross-checks but do not satisfy the paper's claim.
- Direction (b) "script is correct, paper should be updated": stage_010.tex equations `eq:stage010-du2` and `eq:stage010-du4` get rewritten with the script's `δP_2`, `δP_4` closed forms (the eight-term and twelve-term expressions in sympy `:60-67` and `:71-82`).
- Direction (c): deeper review by user.

### F2 — paper_misalignment

**Severity:** medium
**Subtype:** paper_missing_script_claim
**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_010.tex` (entire card — no mention of any of the items below)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_sympy_audit.py:85-169`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage010_projected_maxwell_push_bundle_master_mathematica_audit.wl:65-220`

**What's wrong:**

The scripts verify seven content clusters that the paper card never mentions:
1. One-pole K-surface from `(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0)(T + \varepsilon z_4) = 3(S + \varepsilon z_2)^2` (sympy `:85-89`)
2. Normalization K-surface from `(N_0 + \varepsilon n_0)/(K - B_0 - Z_{0,\mathrm{slot}} - \varepsilon z_0) = P_{\mathrm{target}}` (sympy `:90-93`)
3. Compatibility surface and its first variation `dcompat = n_0/P_{\mathrm{target}} - 6 S z_2/T + 3 S^2 z_4/T^2` (sympy `:101-119`)
4. Transported-target compatibility analog with `P_{\mathrm{target,transport}} = (N_0 + \varepsilon n_0)/D_{0,\mathrm{target}}` (sympy `:105-130`)
5. Two negative-control `assert_nonzero` "mutation" guards (sympy `:132-139`)
6. Real-Y_{20}-square Gaunt overlap ratios giving lane multipliers (1, 1/2, -1) and a grouped trace/anomaly decomposition along b = 3a (sympy `:141-155`)
7. A "primitive static Xi" formula using new symbols Q, Delta, P, q_1, p_1, d_1 under the substitution `N_{0,\mathrm{sym}} = P^2/\Delta^2` (sympy `:157-169`)

The paper card (`stage_010.tex`) is short — eight short paragraphs — and mentions none of: S, T, P_target, Z_0_slot, D_0_target, real Y_{20} overlap, Gaunt coefficients, trace/anomaly decomposition, b = 3a line, primitive Xi, Delta, Q, q_1, p_1, d_1. The part appendix row adds no detail beyond "Projected Z_n, N_n slot transport for the grouped response bundle." The notes file that would normally explain the script's intent does not exist in the committed repo (this is one of the EM-projected stages whose per-stage notes were not committed).

**Why this matters:**

If these seven clusters are load-bearing for downstream stages (especially stage 011 "P_2 bridge", stage 012 "primitive bridge", stages 013-014 "mouth-Taylor"), the paper card understates what stage 010 actually does, and a reader cannot reconstruct stage 010's role from the printed paper. Conversely, if these clusters are scaffolding left over from an earlier draft (the script's banner says "step_08_projected_maxwell_push_bundle_master_notes" — the EM-projected step number, not the stage number — suggesting the script was written against an unpublished notes file), they are doing more work than the paper requires and are at risk of bit-rotting against the paper's `Output` line.

The direction of resolution (extend paper vs. trim script vs. file new notes) is the user's call.

**Required change:**
See `## Resolve before fix_loop` block in `/var/projects/toy_physics/research/pde_ledger/redteam/directives/stage_010.md`. Routed to user.

**Verification:**
After the user picks a direction:
- Direction (a) "paper is incomplete, extend it": stage_010.tex grows new paragraphs listing the K-surface compatibility, transported target, Gaunt-overlap multipliers, trace anomaly, and primitive static Xi — and equation labels for each.
- Direction (b) "script over-verifies, trim it": the script blocks at sympy `:85-169` (and mathematica `:65-220`) are removed, leaving only `δP_0` (and possibly `δP_2`, `δP_4`) as audit content.
- Direction (c) "file the missing notes": author publishes `notes/em_projected/step_08_projected_maxwell_push_bundle_master_notes.md` (or equivalent) to anchor these script-side claims, and the stage card optionally cross-refs them via `\StageFile`.

## Independent-derivation check (Mathematica)

The `.wl` derives the same content via different idioms: `Solve[..., K]` and `Coefficient[Normal[Series[...]], eps, 1]` rather than `sp.solve` + `sp.diff(...).subs(eps, 0)`; `ThreeJSymbol` composed into a custom `gauntByThreeJ` rather than calling `sympy.physics.wigner.gaunt`; intermediate naming `den0/num0/slot0` rather than `D0p/N0p/P0p`; K-surfaces solved fresh via `Solve` rather than computed by hand and then verified. The series-coefficient idiom (mathematica `:39-41`) is genuinely Mathematica-native; the SymPy script uses `sp.diff(...).subs(eps, 0)` instead. The Gaunt overlap is constructed from `ThreeJSymbol` (mathematica `:145-156`) in a Wigner-3j formula independent of SymPy's `gaunt()` routine, which uses a different internal recursion. I do not find evidence of line-by-line transliteration; this is not a `mathematica_transliteration` finding.

## Engine cross-check

Both scripts report `STATUS: PASS` and the Mathematica residuals (M1-M15 all zero, M16/M17 mutations equal `6 S^2 z_4/T^2` consistent with the sign-flip definition) match the SymPy assertions one-for-one. Outputs are fresh: sympy script mtime 2026-05-21 11:28, output 11:51; mathematica script 11:29, output 11:51. The engines agree on every shared check, including the two negative-control mutations.

## Verdict justification

Verdict is `findings` because of two `paper_misalignment` items routed to the user. The script's math is internally consistent and the two engines agree (every assertion I attacked — sign on `dcompat`'s `3 S^2 z_4/T^2` term, the z_0-cancellation in transported-target, the Gaunt lane multipliers, the trace coefficients 1+2+2=5 — held up). The defect is paper-grounded: the script verifies one of the paper's three `\stagefield{Output}` deliverables (`δP_0`) faithfully, substitutes a related-but-distinct quantity for the other two (`δP_2`, `δP_4` instead of `δu_2`, `δu_4`), and adds a large amount of extra content (K-surfaces, Gaunt, trace, Xi) that the paper card never introduces. Under v1 (which did not read the paper) the unit looked clean modulo a missing Mathematica engine; under v2 the paper↔script gap is the dominant finding. No `stop_cold` — both findings need user direction, not unilateral patching, and a downstream-impact analysis can only happen after the user has chosen a direction.

## Self-test notes

Checked: (1) sign convention for D_n -> D_n - \varepsilon z_n is consistent between paper (where `Z_n -> Z_n + \varepsilon z_n` and `D_0 = K - B_0 - Z_0`) and both scripts. (2) Hand-derived `u_2 = -D_2/D_0` reproduces the paper's `δu_2 = (D_0 z_2 - D_2 z_0)/D_0^2`; hand-derived `u_4 = (D_2^2 - D_0 D_4)/D_0^2` reproduces the paper's `δu_4` exactly (z_4/D_0 - 2 D_2 z_2/D_0^2 - D_4 z_0/D_0^2 + 2 D_2^2 z_0/D_0^3) — so the paper's u_n have a clear D-only inversion interpretation, distinct from the script's P_n. (3) Confirmed `δP_0 = (D_0 n_0 + N_0 z_0)/D_0^2 = n_0/D_0 + N_0 z_0/D_0^2` is identical to the paper's form, so claim C4 is satisfied. (4) Setting `n_0 = n_2 = 0` in the script's `δP_2` closed form leaves an N-dependent residual `N_2 z_0/D_0^2 + 2 N_0 z_2/D_0^2 - 4 D_2 N_0 z_0/D_0^3`, confirming `δP_2 ≠ δu_2` even in the n-zero limit. (5) Output transcripts are fresh against script mtimes; no `stale_output` finding. No directive ambiguity to self-test — both findings are user-resolution paper_misalignment and the directive does not prescribe Codex edits.
