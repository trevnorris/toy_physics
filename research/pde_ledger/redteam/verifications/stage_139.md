---
unit_id: 139
batch: IV.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 139 (post-rework)

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:32-60` replaces the lone `R_q^comp = 1/4` assertion with a block of nine `assert` statements (six closeness checks against notes literals — `r_F1`, `R_q^nat`, `M_s^nat,*`, `M_q^nat,*`, `M_s^comp,*`, `M_q^comp,*`; two outlet-consistency residuals; one algebraic `R_q^comp - 1/4`). Tolerances are `1e-12` for the literal checks and `1e-25` for the algebraic identities. A trailing `print('all assertions passed')` was added at line 60.
- `mathematica/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:66-79` replaces the single `expectApprox` line with nine `expectApprox` calls covering the same six literal deliverables, both outlet-consistency residuals, and the `R_q^comp - 1/4` identity.

**Assessment:**
The edits match the directive's prescribed block verbatim in both engines. The six literal-deliverable checks are non-tautological: each compares a quantity computed from `Pi_star`, `Sq_star`, and `rF` against an independent notes literal, so a regression in any of those inputs (or in the construction formulas) would now be caught. The outlet-consistency residuals are algebraically built in by the `M_s`/`M_q` definitions, but they still pin the imported `Pi_*` and `S_q` literals satisfy the form, as the directive acknowledged. SymPy exec log ends with `all assertions passed`; Mathematica exec log shows nine numerical-deliverable PASS lines (plus the F2-added `g_minus closed form` PASS) for a total of 10 PASS lines.

### F2 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/.../audit.wl:38-49` replaces the bare `gMinus = N[rF - Sqrt[1 + rF^2]/2, 30]` assignment with a `Solve[(gc - rF)^2 == (1 + rF^2)/4 && gc < rF, gc, Reals]` derivation, a branch-count guard (`Length[gMinusSolutions] =!= 1` triggers FAIL/Exit[1]), and a closed-form cross-check `expectApprox["g_minus closed form", gMinus, rF - Sqrt[1 + rF^2]/2, 10^-25]`.

**Assessment:**
The Mathematica side now independently re-derives `g_c` from the compensated-branch defining equation rather than echoing the SymPy closed form. The branch-selection guard correctly rejects ambiguous cases. The exec log shows `PASS: g_minus closed form` with residual `0``29.099515036788617` (well below `10^-25`), confirming the Solve result agrees with the closed form. The SymPy script is unchanged for this finding, per the directive. With the F3 syntax regression repaired, the Solve-based derivation is now exercised in the transcript.

### F3 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/.../audit.py:5,7` adds two provenance comments: `# r_F1 closed form from Stage 121.` and `# Pi_* and S_q(Pi_*) imported from Stage 134.`
- `mathematica/.../audit.wl:28,30` adds the matching Mathematica comments. The provenance comment for the imported pair was rewritten by the orchestrator as `(* piStar and sQStar (= S_q evaluated at piStar) imported from Stage 134. *)` using ASCII-safe names, avoiding the `Pi_*)` substring that previously triggered Mathematica's block-comment terminator (pitfall #13 from IV.4).

**Assessment:**
Provenance comments are in place in both scripts and now cite Stage 121 (for r_F1) and Stage 134 (for Pi_* and S_q(Pi_*)) — the post-renumber stage IDs per the orchestrator note (not the pre-renumber Stage 223 / Stage 236 numbers the original directive referenced). The ASCII-safe rewrite of the Mathematica comment resolves the prior `Syntax::sntx` regression; the rerun exits 0 and the exec log shows no syntax warning.

### F4 — stale_output (banner typo)

**Classification:** resolved

**What changed:**
`mathematica/.../audit.wl:26` reads `banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];`.

**Assessment:**
The banner is corrected. The Mathematica exec log line 3 confirms the transcript reads `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`. Per orchestrator note, F4 belongs to the Cluster A banner sweep, and the post-cluster output reflects the fix.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `r_F1 = 1.77799353547497811851563491229`
- `R_q^comp = 0.250000000000000000000000000000`
- `all assertions passed`

All printed values match the notes' boxed deliverables to ~15 digits, and the closing `all assertions passed` sentinel confirms all nine SymPy `assert` statements held.

**Mathematica:** exit=0. Notable lines:
- `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS` (F4 banner fix confirmed)
- `PASS: g_minus closed form` (F2 independent-derivation cross-check)
- `PASS: r_F1`, `PASS: R_q^nat`, `PASS: M_s^nat,*`, `PASS: M_q^nat,*`, `PASS: M_s^comp,*`, `PASS: M_q^comp,*` (six numerical-deliverable checks)
- `PASS: outlet consistency nat`, `PASS: outlet consistency comp`, `PASS: R_q^comp - 1/4`
- `Stage 139 Mathematica audit passed.`

Ten PASS lines total, matching the orchestrator's expected count. No `Syntax::sntx` warning appears, confirming the comment-terminator pitfall is resolved.

**Output freshness:**
- SymPy script mtime 2026-05-27 19:49, output mtime 2026-05-27 19:51 (output newer than script).
- Mathematica script mtime 2026-05-27 20:05, output mtime 2026-05-27 20:07 (output newer than script).
Both transcripts are fresh post-fix.

## Material-change assessment

`material_change`: false.

The edits added assertions and provenance comments and replaced a bare assignment with a `Solve`-based derivation. The numerical values printed (`r_F1`, `R_q^nat`, `M_s^nat,*`, etc.) are unchanged to all reported digits compared with the pre-rework run; downstream stages that consume Stage 139's numerical outputs see the same values they always did. The orchestrator's blanket `upstream_stale: true` for units > 139 is a standard procedural marker; no substantive numerical drift to propagate.

## Side observations (non-blocking)

- The SymPy assertion for `R_q^nat` against the notes literal `0.145454452260421` only depends on `rF` (not on `Pi_*` or `S_q`), so this specific check does not pin the imported `Pi_*`/`S_q` literals. The actual pinning of those imports comes from the `Ms_nat`/`Mq_nat`/`Ms_comp`/`Mq_comp` literal checks, which depend on both imports. This is fine in practice; flagging only for completeness.
- The Mathematica exec log prints precision tags like `28.807067524193343` on some values — these are normal Wolfram precision-tracking artifacts, not warnings.

## Verdict justification

All four findings are resolved post-rework. F1's assertion blocks are in place in both engines and non-tautological for the six headline deliverables. F2's Mathematica `Solve`-based independent derivation of `g_minus` now executes and cross-checks against the closed form. F3's provenance comments cite the correct post-renumber stages (121 and 134), with the Mathematica comment rewritten in ASCII-safe form to eliminate the `Pi_*)` block-comment terminator collision. F4's banner reads `STAGE 139`. SymPy exits 0 with `all assertions passed`; Mathematica exits 0 with 10 PASS lines and `Stage 139 Mathematica audit passed.` Output transcripts are fresher than their scripts. Verdict: `verified`.
