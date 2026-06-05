# Batch III.4 (stages 073-084) — pde_ledger SECOND PASS

Date: 2026-06-05

Part III.4 — Family-1 geometry, thresholds, quadrupole (Family-1 geometry map,
healing-length lock / support scale, threshold window, exact n=5 wall-depth lock,
shell-weighted Theta extraction, branch verdict, quadrupole-Pe demand map, zeta
thresholds, product (Pi) thresholds, master quadrupole residual, direct operator
window, full reduced-PDE write-up).

## Method

- v2 paper-grounded auditor PLUS the **exhaustive script→doc value-reconciliation
  augmentation** (`redteam/pass2/RECONCILIATION_AUGMENTATION.md`) — every audit agent
  read both its rendered prompt and the augmentation doc.
- 12 clean per-stage audit agents in parallel (audit agents execute nothing → 0
  Mathematica seats).
- Independent exec reliability gate — orchestrator re-ran BOTH engines for every
  stage (sequential, ≤1 seat; 084 has no SymPy → 23 execs total, all exit 0) and
  refreshed the committed transcripts via `sed '1,/^---$/d;/^# exit_code:/d' <log> >
  <out>`. The exec re-run + an **arbiter grep** on the refreshed committed outputs
  are ground truth. The audit agents again applied an inconsistent threshold to
  numbering labels (080 was flagged `findings` for the SAME cross-ref docstring class
  that sibling 078 was correctly marked CLEAN for), so the orchestrator authored /
  trimmed every directive and made the self-vs-cross call.
- Codex as the sole fix-applier, with the `.py`-only-vs-`.wl` seat distinction: the 3
  `.py`-only fix stages (074, 077, 081) consume 0 seats; the 2 `.wl`-editing math
  stages (075, 083) consume ≤1 seat each — all 5 ran concurrently within the 2-seat
  Mathematica cap. 080 (deferred, cross-refs only) and 076 (refresh-only, canonical
  source) needed no Codex.
- Clean per-stage verify agent for each source-edit stage (074, 075, 077, 081, 083);
  the 7 no-source-edit stages (073, 076, 078, 079, 080, 082, 084) verified directly
  from the audit verdict + exec/arbiter-grep evidence.

## Result

All 12 stages reached `verified` at v2 depth + value-reconciliation augmentation.
**12/12 verified, `material_change=false` on all 12. No stop-cold, no blocked, 0 Codex
deviations, all iter-1 exit 0.** **No genuine `paper_misalignment` anywhere → ZERO
paper/notes edits.** **No checkpoints in range. No EM-projected stages in range.**

| Outcome | Stages |
| --- | --- |
| Real script-side math finding (de-taut + external-literal anchors) → verified (2) | 075, 083 |
| Real script-side verification finding (closed-form-identity assert) → verified (1) | 081 |
| Real script-side symbol fix (latent `positive=True` trap, no result move) → verified (1) | 077 |
| Self-label fix (`.py` docstring) + output refresh → verified (1) | 074 |
| Refresh-only (`stale_output`; source already canonical) → verified (1) | 076 |
| Deferred-clean (audit flags were all CROSS-refs; orchestrator deferred) → verified (1) | 080 |
| Clean → verified (5) | 073, 078, 079, 082, 084 |

084 is **single-engine by design** (Mathematica-only write-up/consolidation skeleton;
card says "SymPy audit: none yet", `is_status_only_candidate: true`) — the absent
SymPy is not a finding.

## Resolutions

### Real script-side math findings (075, 083) — de-tautologization

Both are the SAME defect family: a check that re-asserts a closed form against its own
definition (`X − X ≡ 0`), which cannot catch a paste error. Both were independently
confirmed by the orchestrator before sending to Codex (false-positive guard), and both
fixes anchor the deliverables against fixed EXTERNAL paper/notes literals, so a wrong
factor / Cosh↔Sinh swap would now FAIL.

- **075 — F1 `tautological_check` (both engines):** `Theta_fail := Upsilon_fail/alpha_r²`
  (py:99), so the round-trip `Upsilon_fail − alpha_r²·Theta_fail` reduces to
  `Upsilon_fail − Upsilon_fail ≡ 0` — yet the SymPy comment explicitly (and falsely)
  claimed it is "not the trivial identity 100·Theta == 100·Theta". **Fix:** removed
  both round-trips from both engines; added SymPy `expect_close` anchors of all 8
  deliverables (`Delta_0`, `Delta_inf`, `Upsilon_fail/suff`, `Xi_fail/suff`,
  `Theta_fail/suff`) against fixed external literals. The reduction's real content
  stays locked by `alpha_r² == 100`; the Mathematica side's pre-existing independent
  `expectApprox` numeric battery (which the SymPy side previously lacked) is retained.
- **083 — F1 `tautological_check` + F2 `insufficient_verification` (both engines):**
  `Delta_0_F1 := numer/denom` (py:57-60), so `denom·Delta_0_F1 − numer ≡ 0` by
  construction — the SymPy comment falsely called it "non-tautological … the unique
  solution of a linear equation" and the `.wl` advertised an "Independent BVP
  derivation" it never performs. Two further `.wl` tautologies: `A_F1 independent vs
  closed-form` compared `(Pi/2)² ≡ Pi²/4` (byte-identical operands), and the
  `Omega(Pe)` residual restated the typed definition. F2: the four `Pe` and five
  `zeta` window deliverables were SymPy-PRINTED but never asserted. **Fix:** replaced
  the SymPy `Delta` residuals with external-literal `expect_close` anchors; added the
  4 Pe + 5 zeta window numeric asserts (F2); deleted the tautological `.wl` `A_F1`
  check; and corrected the misleading `.wl` comments to honestly describe the residuals
  as redundant structural checks. The `.wl` numeric battery (the genuine coverage)
  is untouched.

Both `.wl` were EDITED but the verify agents re-confirmed each **remains an independent
engine, NOT a transliteration**. `material_change=false` on both (de-tautologization +
new anchors over already-correct values; no result value moved — confirmed by the
committed-output diff).

**De-taut cross-corroboration:** the Family-1 healing-window endpoints
`Delta_0 = 1.73302079021525e-4` and `Delta_inf = 2.01447565540522e-2` are now
SymPy-anchored in BOTH 075 (threshold-window route) and 083 (direct-operator route) —
two independent derivations of the same window, mutually corroborating the literals.

### Real script-side verification finding (081)

- **081 — F1 `insufficient_verification` (SymPy):** the exact inversion `Q` was pinned
  only at `Q(0)` and `Q(1)` — two points do not pin a rational function; the full form
  rested solely on `sp.solve`. **Fix:** inserted `expect_zero("Q-closedform", Q − (1 +
  zeta − 2·eps_blk·zeta)/(1 − eps_blk·zeta))`, matching the full closed-form identity
  the Mathematica engine (`wl:54`) already asserts. The flagged F2 comment labels
  (`Stage-35`, `Stage 63`) are CROSS-refs → DEFERRED (orchestrator stripped F2 from the
  directive). `material_change=false`.

### Real script-side symbol fix (077)

- **077 — F1 `symbol_assumption_error` (SymPy):** the integration variable `xi` of the
  symmetric full-line integral (cut point `xi_* ≈ −0.3856 < 0`) was declared
  `positive=True`, contradicting the setup and the Mathematica mirror (`xi ∈ Reals`).
  Dormant (the explicit `(-oo,oo)` bounds, not the symbol assumption, set the domain),
  but a latent trap. **Fix:** split the declaration so `xi` is `real=True` only
  (`alpha_r`/`lambda_mu` stay positive). Verify agent confirmed no surviving assertion
  relies on `xi > 0`; the committed output is byte-identical (`I_f = 1/3`, `xi_*`, and
  the four numeric values unchanged). `material_change=false`.

### Self-label fix (074) — numbering-drift interim policy

- **074 — F2 self-label + F1 `stale_output`:** the SymPy module docstring's filename
  line read `moving_throat_pde_stage57_…` (the file's actual name and banner are 074).
  Per the Reading-2 interim policy, this UNAMBIGUOUS self-label was corrected
  NUMBER-only, FORMAT-preserved: `stage57` → `stage074` (3-digit filename-style,
  matching III.1/III.2/III.3 precedent for filename-docstrings). F1 stale_output: the
  committed banners `STAGE 57`/`STAGE 057` → `STAGE 074` via re-run (the orchestrator's
  independent `exec-mathematica` refreshed the `.wl` transcript, so 074's Codex session
  stayed `.py`-only / 0-seat). The `.py` source diff is strip-the-number identical to
  HEAD.

### Refresh-only (076)

- **076 — F1 `stale_output`:** committed transcripts predated the scripts and carried
  the pre-renumber self-banner (`STAGE 59` / `STAGE 059`) while the current source
  already emits `STAGE 076`. No source edit; orchestrator re-run + sed-refresh updated
  both banners. All result/PASS/closed-form lines unchanged.

### Deferred-clean (080) — orchestrator self-vs-cross call

- **080 — audit `findings` overturned to deferred-clean:** the audit flagged five
  comment/docstring labels (`:5,6,27,35,65`) citing `Stage-61`/`Stage 62`. The
  orchestrator made the self-vs-cross call: these are **CROSS-references** to OTHER
  stages (078 and 079), not stage-080 self-labels (its docstring line 3 reads
  `SymPy audit for Stage 080.`; banner is `STAGE 080`). Per Reading-2, cross-refs are
  owned by the dedicated `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (content-keyed,
  never offset-swept) and are LEFT UNTOUCHED. The auditor over-flagged — sibling 078
  was correctly marked CLEAN for the same class this batch. Codex was NOT invoked on
  080; the directive records the deferral for the audit trail.

### Clean stages (073, 078, 079, 082, 084)

`clean` on first read. 073 = Family-1 geometry four-value freeze (`L/a=37/20`,
`ell/a=1/20`, `Lambda_ell=eta=37`); both engines, all reconcile (073's SymPy is the
`..._sympy_audit_refresh.py` variant — committed output `..._refresh.txt`, present and
fresh). 078 = branch-verdict (four boxed Pe windows + ordering/overlap, dual-engine).
079 = quadrupole demand-to-transport map (Mathematica adds non-mirrored Robin-residual
/ series / slope checks). 082 = master quadrupole residual (the delicate `tan y`-near
root handled by the existing `mpmath` bisection route; no finding touched it). 084 =
full reduced-PDE write-up skeleton (Mathematica-only by design; the load-bearing check
independently re-derives the `zeta_max^F1` ceiling via the `Pe→∞` limit + a
`y tan y = 37` root-solve).

## Value reconciliation (pass-2 augmentation)

Applied on all 12 stages; **119 deliverable values checked batch-wide, 0 misaligned,
0 MISSING-DELIVERABLE.** Per stage: 073=5, 074=6, 075=11, 076=8, 077=9, 078=8, 079=6,
080=5, 081=11, 082=15, 083=21, 084=14. No MISMATCH anywhere.

## Numbering deferrals (SCRIPT/OUTPUT-band plan)

All numbering CROSS-refs / variable-names / ambiguous refs in `.py`/`.wl` source
DEFERRED to `redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` (PENDING — content-keyed,
never offset-sweep): 078 docstring `Stage-60`/`Stage-58`; 080 docstring + four comments
`Stage-61/62` (×5); 081 comments `Stage-35`→052 / `Stage 63`→080. Codex held this scope
perfectly (0 cross-ref touched). The arbiter grep confirmed no stale self-epoch
(NNN−17) banner/closing remains on any committed output. Reference memory
`numbering-drift-root-cause`.

## Infra / seat notes

- `$RT exec-*` writes the fresh transcript to `exec_logs/` and does NOT overwrite the
  committed `output/*.txt` — the orchestrator sed-refreshed every committed transcript
  after the re-run (III.2/III.3 lesson re-confirmed). 23/23 execs exit 0 (11 SymPy + 12
  Mathematica; 084 SymPy correctly skipped).
- Seat policy held: the 3 `.py`-only Codex sessions (074, 077, 081) ran at 0 seats; the
  2 `.wl`-editing sessions (075, 083) ran within the 2-seat cap; orchestrator exec
  sequential, after all Codex done.
- Pass-1 `MANIFEST.yaml` untouched (isolation held under `RT_REDTEAM_DIR=redteam/pass2`).

## Trackers

6 prose trackers synced (PAPER_CLEANUP **P5-08** = ZERO paper/notes edits;
MATHEMATICA_MIRROR_POLICY / CHECKPOINT_TRUST_AUDIT / CHECKPOINT_CONSTANT_PROVENANCE /
EM_PROJECTED_INTEGRATION_TRACKER / STAGE_VERIFICATION_COVERAGE each carry a
`## Pass 2 — Batch III.4` entry).

NEXT = III.5 (085-090, "Part III.5 — quadrupole cancellation, loading ratio, verdict";
checkpoints 089 & 090 in range).
