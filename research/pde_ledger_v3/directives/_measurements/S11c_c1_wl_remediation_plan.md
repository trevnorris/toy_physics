# S11c-c1 WL engine — FULL REMEDIATION PLAN (user-approved 2026-09-04)

## Why this exists
The repair directive `directives/S11c_c1_wl_repair_directive.md` (orchestrator-written, physics-bearing) was
built against **without its rule-7 decision legs** — the orchestrator rationalized the skip via "SymPy precedent"
+ "repairs are lighter" (both false; [[feedback_directive_design_review]] now closes these). The user flagged it.
The 2 decision legs were run **retroactively** (Codex + Grok, both `EXIT=0`) plus an independent arc-verify (Codex).
They found the directive **not sound**, and the two legs caught DIFFERENT defects — one leg alone would have
shipped it. Evidence: `…/scratchpad/{codex_wl_repairdir.log, grok_wl_repairdir.log, codex_c1_arc_verify.log}`.

The engine's **exported** physics is sound and UNaffected: the two propagated defects are **emit-only** — c2
consumes only the DtN operator / two-momentum kernel / face response, all 2-leg-verified. This remediation fixes
the emit-only propagations + the audit-trail + the missing canonical output, then proceeds to T7.

## The confirmed defects (both decision legs + arc-verify)
| # | Defect | Reached code? | Kind |
|---|---|---|---|
| 1 | R1 invariant LEAKS residual=0 + re-enters at `DTN_KERNEL` (rule 5/3/6) | no (code genuinely unfroze) | directive |
| 2 | `NONINVERTIBILITY_CONDITION` over-protected → input-leg freeze PERSISTS in `fredholmFunctionSpaceOperator` (`.wl:580-597`, both `gZero` on `momentumOutput`; used by `NONINVERTIBILITY_CONDITION` `.wl:1537-1549`) | **YES — PROPAGATED**, re-review structurally blind | directive + CODE (emit-only) |
| 3 | ~~R4b parity labels → DELTA_W/ZETA_C SWAPPED in `PERMEABLE_PORT_HERMITIAN`~~ **CORRECTED** (see note ↓): the real defect is a **dead parity axis** in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (`.wl:1713-1736`) | **YES — PROPAGATED** | directive + CODE (emit-only) |
| 4 | R2 "independent data" wording + no φ-driver requirement | no (code bound `V_s`) | directive |
| 5 | R1 `Ŵ`-in-composition over-spec; R1 probe scope; blanket "no residual" disclaimer inconsistency | no | directive NITs |
| 6 | rule-2 records point to ephemeral `/tmp` logs, not literal output; "+443/−103 one file" is engine-scoped not the 4-file commit | n/a | record |
| 7 | neither c1 engine's canonical `.out` committed (WL nor SymPy) | n/a | repo state |

> **⭐ CORRECTION — defect #3 relocated by the repair-2 decision gate (2026-09-04; Codex + Grok, both exit 0;
> record `_measurements/S11c_c1_wl_repair2_directive_review.md`).** The retroactive legs reported a "parity swap"
> in `PERMEABLE_PORT_HERMITIAN`. The repair-2 decision legs (and my own verification, rule 13) **refute** that:
> `PERMEABLE_PORT_HERMITIAN` is **correct** — its `deltaWHermitian=(plus+minus)/4`, `zetaCHermitian=plus+minus`,
> coupling `(plus−minus)/2` are exactly the **congruence** `Aᵀ·diag(P₊,P₋)·A` blocks under `A=faceToParityMatrix
> =[[1/2,1],[1/2,-1]]`, correct given the outward orientation `V_s=s∂_tζ_s` (`S11c_a:58`). The naive
> `DELTA_W=(ζ₊−ζ₋)` reading (which the retroactive finding and my first repair-2 draft both made) is a **category
> error** that vanishes the thickness port at equal faces — building it would have **corrupted a correct object**.
> The genuine propagated parity defect is a **dead parity axis** in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU`
> (curved per-face memory Hermitian form emitted under both `PARITY_*` keys; should be the parity **combination**,
> spec §3b `:308-320`), with `UNIFORM_LIMIT_*` the same pattern on a flat object (NIT). This is the gate working:
> it caught a wrong fix **before** any build. ⇒ Item A below is updated accordingly.

## THE PLAN — four items, in order. ⛔ Do NOT start T7 until all four are done.

### (A) Repair-2 — fix the 2 PROPAGATED emit-only code defects, FULLY GATED
1. **Author `directives/S11c_c1_wl_repair2_directive.md`** (orchestrator-written). Fixes:
   - **Defect 2 (R1):** the `NONINVERTIBILITY_CONDITION` / Fredholm operator must be the **two-leg DtN** `Z`
     inside `[I+(Λ_A/ρ_m²)Z]` — the SAME both-legs construction as `DTN_OPERATOR`/`operatorCompositionFrom
     Derivation` (carrying `q_out(momentumOutput)`, `q_out(momentumInput)`, and the profile source), ⛔ NOT both
     `N_0` factors frozen to `momentumOutput`. ⛔ **Do NOT add a "diagonal-symbol reduction of the full operator at
     `k=k′`"** — the nonuniform operator is not diagonal, `k=k′` selects a diagonal slice, and a fresh dummy there
     only renames the freeze (repair-2 decision-gate catch); the only legitimate diagonal symbol is the
     already-emitted **flat** `FLAT_DIAGONAL_SYMBOL_RELATION`/`DEGENERATE_LOCI_*` (keep unchanged). Tag count stays
     51. Construction invariant: probes `FREDHOLM_Z_HAS_Q_INPUT`/`_OUTPUT` computed from the emitted operator; the
     re-freeze control is **constructor-level** (`rightLegMomentum→momentumOutput`, as `.wl:1329-1332`) and MOVES
     the probe — ⛔ not a post-hoc `momentumInput→momentumOutput` substitution (the `Count` probe is syntactic).
     ⛔ No leaked residual value / no "minus X is zero" (learn from defect 1).
   - **Defect 3 (R2 — CORRECTED, see the note above):** fix the **dead parity axis** in
     `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (`.wl:1713-1736`, per-face memory Hermitian form under both `PARITY_*`
     keys): emit the **parity-combination** memory Hermitian form per parity (via the SAME change of basis
     `PERMEABLE_PORT_HERMITIAN`/`portParityCombination` already uses), its ωτ_I limits computed on the combined
     form. Construction invariant: the `PARITY_DELTA_W` and `PARITY_ZETA_C` payloads are **distinct computed
     objects** in the curved case, and a one-sided corruption of only the `+`-face memory form MOVES them
     differently. ⛔ Do NOT hand-write block arithmetic and ⛔ do NOT pre-register even/odd. `UNIFORM_LIMIT_*` is
     the same pattern on a flat (parity-independent) object — a lighter NIT, fix only if it preserves
     `UNIFORM_LIMIT_RESIDUAL`. ⛔ **`PERMEABLE_PORT_HERMITIAN` is correct (congruence) — leave it byte-identical.**
   - ⛔ Protect byte-identical ONLY the genuinely-sound core (`DTN_KERNEL`, flat symbol, `DTN_OPERATOR`,
     `DTN_RIGID_SHIFT_*`, the kernel-derived `DTN_HERMITIAN/REGIME/PARITY`, the response inverse `FACE_RESPONSE`,
     **`PERMEABLE_PORT_HERMITIAN`** (correct congruence), the repaired `ENERGY_*` audit, the flat diagonal-symbol
     loci, T-a..T-i, μ_θ) — ⛔ NOT `NONINVERTIBILITY_CONDITION`/`fredholmFunctionSpaceOperator` and ⛔ NOT
     `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (those are what's being fixed).
2. ⛔ **2 DECISION legs on the repair-2 directive BEFORE the build** (Codex + Grok — orchestrator-written). Fold once.
3. Detached Mathematica build (setsid+marker+Monitor; danger-full-access; one kernel; RSS watch).
4. 2 re-review legs (fresh Claude Agent + Grok, serialized) — ablate: the Fredholm operator now carries both legs
   (re-freeze control moves it); the parity labels now bind to the definitions (swap breaks the binding residual);
   the sound core byte-identical. ⛔ **EXPLICIT GATE (a re-review gap the arc-verify flagged): each leg forms its
   own target view / writes its own derivation BEFORE opening the artifact** — the repair-1 re-review legs opened
   the artifact first, which weakened them.
5. Commit after both re-review legs report (rule 9). Baseline before repair-2 = `13f0bd2c`.
6. **Documentary disposition for the directive-only defects (1, 4, 5 — NOT inherited):** the original repair
   directive cannot be un-written, so annotate `directives/S11c_c1_wl_repair_directive.md` with a SUPERSEDED/erratum
   note stating defects 1 (R1 leaked residual=0), 4 (R2 "independent data"/no φ-driver), 5 (R1 `Ŵ`-over-spec,
   probe scope, blanket disclaimer) and that the repair-2 directive corrects the class. Record the retroactive
   repair-directive decision-leg findings in `_measurements/S11c_c1_wl_repair_directive_review.md`.

### (B) Rule-2 record corrections
- Fix the repair record's diff-stat claim: state it is the **engine-scoped** diff (`git diff --stat e139bc61 13f0bd2c -- <wl>` = +443/−103, 1 file); the full commit `e139bc61..13f0bd2c` is **4 files, +691/−103**.
- Fold the **launch/run COMMANDS and the LITERAL key leg outputs** (verdicts + the decisive ablation before/after
  lines) INTO **all three** committed records — `_measurements/S11c_c1_wl_build_directive.md` (which lacks literal
  reviewer output), `_measurements/S11c_c1_wl_build_review.md` and `_measurements/S11c_c1_wl_repair_directive.md`
  (which lack both commands and literal stdout) — ⛔ not just `/tmp` pointers (the `/tmp` logs are ephemeral, rule 2).
- Record the **retroactive repair-directive decision-leg** results (prompt `_legs/S11c_c1_wl_repair_directive_review.md`)
  in a committed `_measurements/S11c_c1_wl_repair_directive_review.md`, with the Codex+Grok verdicts and the
  arc-verify, and their commands + literal key output.

### (C) Generate + commit BOTH canonical `.out` transcripts
- Run the reviewed **WL** engine (post-repair-2) → `mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out`;
  run the reviewed **SymPy** engine `scripts/S11c_c1_bulk_closure_sympy_audit.py` →
  `scripts/out/S11c_c1_bulk_closure_sympy_audit.out`. `datalad save` (annex policy auto-annexes `out/*.out`), then
  `datalad push --to gin` + `git push origin`. ⛔ Never `git add -f` a big `.out`. These are the audit artifacts
  AND the T7 comparator's inputs (see the `project-datalad-gin-out-storage` memory).

### (D) Correct the "DONE" overstatement
- In records/memory: "c1 WL per-engine review complete" is real; "c1 DONE" is not until (A)–(C). Fixed in
  [[project_s11c_c_state]].

## THEN → T7 (only after A–D). See [[project_s11c_c_state]] for the T7 contract + the known representational
residuals to adjudicate after the run (μ_θ, omega-assumption, the two-momentum leg names).

## Process caveats surfaced (apply going forward, not necessarily remediated)
- Repair re-review legs opened the artifact before their own independent derivation (a re-review weakness — form
  the target view first even in a re-review).
- One build-review leg (fresh Claude) printed `False` on its operator-inverse self-check while reporting
  "confirmed"; the verdict survives on Grok + the artifact, but that specific self-claim wasn't cleanly backed.
