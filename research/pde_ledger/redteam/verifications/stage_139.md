---
unit_id: 139
batch: IV.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T18:25:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 139 (REMEDIATION / recovery review)

This verdict SUPERSEDES the prior 2026-05-27 verification (claude-opus-4-7), which described the
pre-rework scripts and erroneously praised the `Solve`-based `gMinus` derivation and the
tautological outlet asserts as "resolved" — exactly the constructs the recovery review later
flagged as R1/R2. Authoritative findings here are `redteam/codex_reviews/stage_139.md`
(verdict FINDINGS, R1–R4). Directive `redteam/directives/stage_139.md` was rewritten to encode
them, with F2 further reworked after the Claude+Codex consult
(`redteam/codex_reviews/_consult_batch5.md` Q4/Q5). Verifier scope is scripts-only; exec logs
read, no scripts run.

## Per-finding outcomes

### R1 (F1) — tautological_check (outlet-consistency residual was X−X)

**Classification:** resolved

**What changed:**
- SymPy `moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py`: the two
  `assert abs(Pi_star - (Ms_* + Mq_* * Sq_star)) < tol_algebraic` asserts are GONE. They are
  replaced by an independent `Sq_recon` reconstruction (py:62-71) and the old residuals survive
  ONLY as labeled structural `print` lines (py:73-79).
- Mathematica `...mathematica_audit.wl`: the two `expectApprox["outlet consistency ..."]` asserts
  are GONE; replaced by `sQRecon` from the kernel (wl:77-85) plus two structural `Print` lines
  (wl:86-88).

**Assessment:**
Correct and non-tautological. The new check recomputes `S_q(Pi_*)` from the Stage-134 closed-form
kernel `S(Pi,kappa) = Pi*(kappa*tanh(kappa) + Pi*(exp(-Pi)/cosh(kappa) - 1)) /
((1-exp(-Pi))*(kappa^2 - Pi^2))` at `kappa = pi/2`. The kernel route references only `Pi_star`
and `tanh`/`cosh` of `pi/2` — it never touches `Ms_*`, `Mq_*`, or `R_q` (confirmed by reading
py:66-71 and wl:81-85). Therefore the assert `|Sq_recon - Sq_star| < tol_literal` can genuinely
fail if either the imported literal or the kernel transcription is wrong; it is NOT X−X. The
denominator uses the `(kappa^2 - p^2)` form mirroring the Stage-134 owner. Exec confirms PASS with
residual ~4.47e-16 (sympy: assertion silent-pass under `all assertions passed`; mathematica log
line 35-36: `S_q recon from Stage 134 kernel residual = 4.468...e-16`, `PASS`). The demoted outlet
residuals print `0` (sympy log 17-18; mathematica log 37-38) and are no longer asserts.

### R2 (F2) — tautological_check (R_q^comp = 1/4 is DEFINITIONAL and branch-blind)

**Classification:** resolved

**What changed:**
- (a) `R_q^comp - 1/4` assert RELABELED as definitional-consistency: sympy py:81-84 comment +
  `assert abs(Rq_comp - sp.Rational(1, 4)) < tol_algebraic, (Rq_comp,)`; mathematica wl:89
  `expectApprox["R_q^comp - 1/4 (definitional-consistency)", rQComp, 1/4, tolAlg];`.
- (b) NEW branch-discrimination anchor `g_-^F1 value` vs the canonical cross-stage literal
  `0.758035078944662826919680890414`: sympy py:49-52; mathematica wl:69-72.
- (c) Mathematica `gMinus` now computed DIRECTLY as the closed form
  `gMinus = N[rF - Sqrt[1 + rF^2]/2, 30];` (wl:45). NO `Solve` statement survives.
- (d) `R_q^nat = 0.145454452260421 ≠ 1/4` intact: defined sympy py:13, asserted py:47;
  defined mathematica wl:36, asserted wl:68.

**Assessment:**
Correct on all four points. (a) The `1/4` check is explicitly labeled definitional-consistency in
both engines and is no longer presented as the anchor; the comments state it is branch-blind (both
`g_±` give 1/4 because the `±` is squared away). (b) The `g_-^F1` value anchor is present and is the
branch DISCRIMINATOR — a sign/branch typo (`g_+ ≈ 2.79`) fails `tol_literal=1e-12` against
`0.758...`. Per the directive's own guidance and the consult orchestrator-refinement, this anchor
is correctly framed as branch-discrimination, NOT an independent re-derivation of `g_-` (g_- has no
independent definition; it IS `rF - sqrt(1+rF^2)/2`), so it is not flagged tautological. (c)
Confirmed by grep: the only occurrence of `Solve` / `(1 + rF^2)/4` in the `.wl` is inside the
explanatory COMMENT at wl:44 ("NOT by solving (gc - rF)^2 == (1 + rF^2)/4, which would re-bake
1/4."), not an executed statement. No `gMinusSolutions`, no `Solve[...]` code remains. (d) The
natural-branch `R_q^nat` is a genuinely different computation `(1-rF)^2/(1+rF^2)` not forced to any
value, proving the branch selection is real. Exec: mathematica log shows `PASS: g_-^F1 value`
(residual 0) and `PASS: R_q^comp - 1/4 (definitional-consistency)` (residual 0); `R_q^nat` residual
-8.88e-16, PASS. The old tautological PASS lines `outlet consistency nat/comp` and `g_minus closed
form` are absent from the log.

### R3 — paper_misalignment (RESOLVED, no script fix)

**Classification:** resolved

**What changed:**
A one-line provenance COMMENT only, in both engines: sympy py:10-11 and mathematica wl:33-34, both
stating the self-matched susceptibility closure `Sigma_0 = M_s` is established at Stage 140 (not
139), per P4-42.

**Assessment:**
Correct. No susceptibility-closure residual check was added to Stage 139 (none appears in the
script body or the exec logs), matching the directive's RESOLVED:R3 instruction. The Mathematica
comment (wl:33-34) is ASCII-safe and contains no `*)`, `_*)`, or `Pi_*)` substring (grep of the
file shows the only `Pi_*)`-shaped text is the legitimate `Print["S_q(Pi_*) = ", ...]` string
literal at wl:51, not a comment hazard). The comment closes cleanly and the rerun shows no
`Syntax::sntx` warning.

### R4 — provenance (RESOLVED, scripts already correct)

**Classification:** resolved

**What changed:**
No change required; scripts retain Stage 121 (`r_F1`, py:5 / wl:28) and Stage 134 (`Pi_*`,
`S_q(Pi_*)`, py:7 / wl:30) provenance.

**Assessment:**
Correct. Grep for `223` and `236` across both scripts returns nothing (exit 1) — the stale
pre-renumber stage numbers were not reintroduced. The 121/134 numbers are the canonical
post-renumber values per the directive's RESOLVED:R4 block.

### F-BANNER — stale_output (banner guard)

**Classification:** resolved

**What changed:** None needed.

**Assessment:** `mathematica/...wl:26` reads `banner["STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS"];`
and the mathematica exec log line 8 prints `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS`. The SymPy
script carries no conflicting stage label. Confirmed.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `g_-^F1 = 0.758035078944662826919680890414` (log line 11) — matches the branch anchor exactly.
- `R_q^comp = 0.250000000000000000000000000000` (line 12) — definitional 1/4, as expected.
- `outlet form residual (nat, structural) = 0` / `(comp, structural) = 0` (lines 17-18) — print-only.
- `all assertions passed` (line 19), `# exit_code: 0` (line 20).

**Mathematica:** exit=0, 9 PASS / 0 FAIL. Notable lines:
- Banner `STAGE 139 — ACTUAL FAMILY-1 MOUTH GAINS` (line 8).
- `PASS: r_F1`, `PASS: R_q^nat`, `PASS: g_-^F1 value`, `PASS: M_s^nat,*`, `PASS: M_q^nat,*`,
  `PASS: M_s^comp,*`, `PASS: M_q^comp,*`, `PASS: S_q recon from Stage 134 kernel`,
  `PASS: R_q^comp - 1/4 (definitional-consistency)` (lines 22-40).
- `S_q recon ... residual = 4.468...e-16` (line 35); `g_-^F1 value residual = 0` (line 25).
- The old tautological PASS lines (`outlet consistency nat/comp`, `g_minus closed form`) are
  ABSENT. `# exit_code: 0` (line 43).

**Output freshness:** confirmed re-generated post-fix. The saved `.txt` outputs live at
`scripts/output/moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.txt` and
`mathematica/output/moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.txt`
(not alongside the scripts). Their mtimes (2026-05-29 18:02:15) are newer than the script mtimes
(2026-05-29 16:51:54 sympy / 16:52:19 mathematica) and newer than the exec logs (17:58). Their
tails match the exec-log content (`all assertions passed`; `Stage 139 Mathematica audit passed.`
with the new `S_q recon` / definitional-consistency PASS lines).

## Material-change assessment

`material_change`: false.

No derived numeric value changed. Every boxed deliverable (`r_F1`, `R_q^nat`, `M_s^nat,*`,
`M_q^nat,*`, `R_q^comp`, `M_s^comp,*`, `M_q^comp,*`, `g_-^F1`, `S_q(Pi_*)`) computes identically to
before. The edits only strengthened the verification surface: two structurally-zero asserts were
demoted to prints and replaced by an independent kernel reconstruction; a branch-discriminating
`g_-^F1` value anchor was added; the `1/4` check was relabeled; the Mathematica branch solve was
swapped for the equivalent direct closed form (same value). Downstream units depending on the 139
gains are unaffected by value, so no narrow re-audit is needed on numerical grounds.

## Side observations (non-blocking)

- The directive's inline line-number references (e.g. "line 41", "line 58", "line 71/72") drifted
  by a few lines after the edits landed (e.g. the `g_-^F1` anchor sits at py:52, the relabeled
  `R_q^comp` assert at py:84). Every required construct is present at its actual location; only the
  citations are slightly stale. Non-blocking.
- The committed `.txt` outputs are under `scripts/output/` and `mathematica/output/`, not the paths
  the verifier prompt's "Output freshness" note implied (alongside the scripts). Located and
  confirmed fresh; non-blocking.
- A stale prior-iteration verdict file (claude-opus-4-7, 2026-05-27) was overwritten by this one; it
  described the pre-rework scripts and should be disregarded.

## Verdict justification

All four authoritative findings (R1–R4) plus the banner guard are resolved. The two formerly
tautological outlet asserts are gone and replaced by a genuinely independent Stage-134 kernel
reconstruction of `S_q(Pi_*)` that can fail; the `R_q^comp = 1/4` check is honestly relabeled as
definitional-consistency and the load-bearing falsifiability now rests on the `g_-^F1` value
(branch discrimination) and `R_q^nat ≠ 1/4`; the Mathematica `gMinus` is the direct closed form
with NO surviving `Solve` (the `(1+rF^2)/4` text exists only as a comment); provenance stays at
Stage 121/134 with no 223/236 regression; and the banner is canonical. Both engines exit 0 with all
assertions passing (9 PASS / 0 FAIL in Mathematica), outputs are freshly regenerated, and no derived
value changed. Verdict: verified.
