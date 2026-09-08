# Independent build-review leg — S11c-c2 N6 ablation harnesses (executable)

**Verdict: SOUND.** The four astra-built harnesses are faithful, non-fakeable ablation certs that
implement the cleared knife list. Every canonical run and every knife FORM run I recomputed
independently reproduces the committed transcript's SHA-256 digests **exactly** (0 mismatches across
all four harnesses), which simultaneously proves (a) the harness runs the LIVE engine (no hardcoding),
(b) the per-object digests are faithful full-payload SHAs, and (c) each knife is genuinely implemented
and bites as claimed. One **non-blocking finding**: the K_split_route **×2 coefficient companion is
inert** on the WL pin (its FORM knife bites normally).

## Method (how I ablated the harness, not the self-report)

I did NOT trust the transcripts. For each SymPy harness I re-ran its own `worker`/`patch_source` code
path on /tmp trees (fresh subprocess per variant) and compared my independently-computed per-object
`{sha256, nonzero-fingerprint}` to the committed CANONICAL and FORM records. For WL I re-ran the
production Mathematica engine under `timeout --kill-after=5 600` (one kernel at a time) on /tmp
budget+knife-patched copies and compared to committed. Matching CANONICAL digests defeats a
hardcoded/faked cert; matching FORM digests + observing the certified object move defeats an
unimplemented/inert knife. All ablation was on /tmp copies; the working tree was never modified
(confirmed: no harness file shows ` M` in git; all are pre-existing untracked).

Scripts (absolute paths):
- SymPy driver: `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/ablate_sympy.py`
  logs `…/scratchpad/{cov,diag,rec}_ablate.log`
- WL single-variant runner: `…/scratchpad/wl_ablate/wl_one.py` (patched driver `…/wl_ablate/driver.py`,
  transcript redirected to /tmp; per-variant object dumps `…/wl_ablate/obj_*.json`, run logs `…/wl_ablate/run_*.log`)
- transcript parsers: `…/scratchpad/parse_transcript.py`, `…/scratchpad/summarize.py`

## Compact-contract verdict (all four transcripts) — CLEAN

Grep of all four `_measurements/*ablation*.md`: **0** occurrences of `numerator_denominator` (PIT sample
matrix), `ARITHMETIC`(DAG/circuit), `root_nodes`/`_NODES`, `srepr`/`Symbol(`, `SparseArray`,
`PROBE_NUMERATORS`. The only `sp.`/diff text lives inside MANIFEST/RUN_MANIFEST knife-definition fields
and the required unified source diffs — i.e. the harness's own knife text, never leaked engine payloads.
- WL transcript (~1 MB) is **appropriately-compact numeric fingerprints**, not over-verbose: PIT_FINGERPRINT
  entries are `{valid,rejected,nonzero_count,circuit_leaves,numerator_count}` integers keyed by
  `{case,"PROBE",cell,prime}`; the size is many (object×cell×prime×variant) numeric entries plus
  deduplicated COMPONENT_KEYS schemas (structural key names/indices), with zero symbolic content.
- **Error path verified via a real forced error**: my first WL CANONICAL kernel was SIGKILLed (exit 137).
  Its INTEGRITY_GUARDS emitted only bounded diagnostics — `engine.stdout`/`worker.stderr` recorded as
  byte-count + SHA-256 (hashed, never echoed), `production_immutable=1`, no raw object payload anywhere
  (whole failed transcript 12.8 KB). Clause-4 error-path contract holds.
- WL budget patches: exactly the two specified NON-KNIFE edits (driver `Do` iterator → single pinned
  case; `probeCase` drawCount → 4), each `frag_count=1`, classified `BUDGET_RESTRICTION / NON-KNIFE`,
  chained on the production source and hashed separately from the three knives.

## Live-engine / digest-faithfulness (anti-fake)

| Harness | CANONICAL sha match | per-knife FORM sha match (committed) |
|---|---|---|
| covariance | 7/7 | K_junk 7/7, K_circular 7/7, K_rank 7/7 |
| diagnostic | 20/20 | K_EW_rowdrop 16/16, K_slotdrop 20/20 |
| reconcile | 13/13 | K_normal 13/13, K_source_route 13/13, K_operand_swap 13/13 |
| WL (Mathematica) | 15/15 | K_junk 11/11, K_split_route 10/10, K_carrier 12/12 |

All 0 mismatches. A hardcoded/faked cert could not reproduce these from independent runs.
**Digest is a genuine full-payload SHA, not a re-hash of the compact fingerprint** — proven by covariance
**K_circular/FORM**, a same-support bite: my independent run moves the SHA of all 5 certified covariance
objects while the nonzero fingerprint is UNCHANGED (nz 0→0, 146→146), matching committed exactly. If the
digest were a re-hash of the fingerprint it could not move while the fingerprint held.

## Per-knife verdict (bites / inert / faked)

Legend: FORM moves cert fingerprint and/or full-payload digest vs CANONICAL; NO-OP/IDENTITY reproduces
CANONICAL exactly; ×2 companion also moves. All "match" = my independent run == committed record.

### Harness 1 — WL (pin `{MATERIAL_ADVECTED, RHO4_CONSTANT}` = `actualJunkCase`)
- **K_carrier — BITES.** FORM moves N6RC_R_N6 (nz 0→15480), SPLIT_CHECK, SLOT/CLOSURE guards (runtime,
  12/12 match committed). N6COV_R_COV correctly dead (in K_carrier DEAD set) → unmoved. IDENTITY=canonical.
  Code: site `materialNormalKnife=0;` unique (line 17), feeds `materialGeometry` off-block contamination.
- **K_junk — BITES ON THE PIN** (the prior "silently disabled by case pin" defect is RESOLVED). Runtime
  FORM moves N6COV_R_COV (nz 0→1224), N6RC_R_N6 (0→5184), SPLIT_CHECK, CLOSURE_GUARD; SLOT_GUARD unmoved
  (matches committed; 11/11). DEAD (CARRIER_*, SOURCE_PREDICTED, FROZEN_PHI, PHI_DOMAIN_CENSUS) unmoved.
  **Harness FORM-ablation**: neutering K_junk's FORM to a NO-OP reproduces CANONICAL exactly (15/15 identical
  sha; raw stdout byte-identical) → the FORM application is load-bearing. Code: `actualJunkCoefficient=0;`
  unique (line 20); gate `activeJunk=If[case===actualJunkCase, actualJunkCoefficient, 0]` (line 845) with
  `actualJunkCase={"MATERIAL_ADVECTED","RHO4_CONSTANT"}` (line 21) = the pin ⇒ junk live on the pin.
- **K_split_route — FORM BITES; ×2 companion INERT (FINDING).** Runtime FORM moves SPLIT_CHECK
  (nz 0→14688), everything else (incl. R_N6, R_COV, guards, all DEAD) unmoved (10/10 match committed) —
  the intended one-sided admission of the affine family. IDENTITY=canonical. **But K_split_route/×2 moves
  NOTHING** (committed: SPLIT_CHECK digest_equal=1, no fingerprint delta; the ×2 patch DID apply —
  executed_sha ≠ budget_sha, `es[#]-ms[#]`→`2 (es[#]-ms[#])`). Cause: with the canonical `affine=False`
  route the sourceChannel contributes zero to the printed split, so doubling `2×0` is inert; only the FORM
  route-flag (`False→True`) admits a nonzero family. Impact: K_split_route is genuinely implemented (FORM
  bites), but its arithmetic-sensitivity self-test is vacuous on the pin. NON-BLOCKING.

### Harness 2 — covariance (SymPy)
- **K_junk — BITES.** FORM moves R_COV (nz 0→2), R_COV_CONTROL_DELTA (0→2), SOURCE_ACTUAL (same-support
  digest); DEAD (SOURCE_PREDICTED, FROZEN_PHI, PHI_DOMAIN_CENSUS) unmoved; IDENTITY=canonical; ×2 bites (committed).
- **K_circular — BITES (same-support, digest-witnessed).** FORM moves the SHA of all 5 certified with
  fingerprint unchanged (nz 0→0 / 146→146); DEAD SOURCE_ACTUAL only sampler-digest-move (nz unchanged);
  IDENTITY=canonical; ×2 bites R_COV fingerprint (committed).
- **K_rank — BITES.** FORM moves R_COV (nz 0→120), R_COV_INCREMENT (0→18), SOURCE_PREDICTED (146→128),
  FROZEN_PHI; IDENTITY=canonical; ×2 bites (committed). Note: DEAD SOURCE_ACTUAL shows a **column-set
  (MISSING) change** because K_rank shrinks the shared wave-basis (zero columns dropped); shared-column
  nonzero content is unchanged (nz 146→146). This is DEAD-side over-sensitivity, surfaced honestly as
  MISSING (never coerced to zero) — not a physics leak.

### Harness 3 — diagnostic (SymPy)
- **K_EW_rowdrop — BITES.** FORM moves REP_INVARIANCE_MATERIAL_OPERAND (nz 30→15), RESIDUAL (18→24), and
  the native TILT/N4_ADVECTION control groups; DEAD (MU_RECONSTRUCTION_*) only sampler-digest-moves
  (nz unchanged 51→51/0→0); IDENTITY=canonical; ×2 bites; the extra `K_ewsign` coefficient companion moves
  MATERIAL_OPERAND/RESIDUAL (coefficient-class, as designed). 16/16 match.
- **K_slotdrop — BITES.** FORM moves REP_INVARIANCE all three (nz 30→24/30→24/18→14) and the slot/closure
  guards (MISSING columns from the dropped slot); DEAD (MU_AMPLITUDE, FACE_VELOCITY, MU_RECONSTRUCTION_*)
  only sampler-digest-moves (nz unchanged); IDENTITY=canonical; ×2 bites (guards nz+2). 20/20 match.

### Harness 4 — reconcile (SymPy)
- **K_normal — BITES.** FORM moves R_N6 (nz 18→69), CARRIER_MATERIAL (8→12), MATERIAL_OPERAND (30→81),
  CARRIER_CHANNEL (0→51), CARRIER_BRIDGE_RESIDUAL (0→4); DEAD (SOURCE_*, EULERIAN_OPERAND) unmoved;
  IDENTITY=canonical; ×2 bites (same 5). 13/13 match.
- **K_source_route — BITES.** FORM moves SOURCE_CHANNEL (nz 18→30) + SPLIT_CHECK (0→12); R_N6 and DEAD
  (CARRIER_*) unmoved; IDENTITY=canonical; ×2 bites SPLIT_CHECK. 13/13 match.
- **K_operand_swap — BITES.** FORM moves SOURCE_CHANNEL (nz 18→0) + SPLIT_CHECK (0→18); DEAD (CARRIER_*)
  unmoved; IDENTITY=canonical; ×2 bites SPLIT_CHECK. 13/13 match.

## DEAD-control / sampler-awareness (clause 6) — no false negatives, no real leaks

Using `digest_equal` + nonzero-fingerprint as the robust per-object channel: **no DEAD object shows a real
(nonzero-fingerprint) leak in any harness.** DEAD "movement" is confined to (i) SymPy sampler-only digest
shifts with unchanged nonzero fingerprint (documented false-positive-possible on DEAD, never false
negative), (ii) WL `circuit_leaves` — a joint-sample-circuit shared-metadata field that shifts for every
object when the joint circuit grows, with `digest_equal=1` on the DEAD object's own payload, and
(iii) covariance K_rank's shared-wave-basis column drop (zero columns), all reproduced at runtime.

## Operational notes
- Every WL kernel wrapped in `timeout --kill-after=5 600`, run ONE at a time. One K_junk/FORM attempt was
  SIGKILLed at ~411 s **not by my timeout** but by the concurrent second leg's `pkill -x WolframKernel`
  cleanup (kills kernels by exact name); I re-ran it in an idle window and it completed (rc 0). My own
  cleanup is path-specific (`pkill -f …/S11c_c2_N6_ablation_harness.wl`), so it never touches the other
  leg's kernel. No orphan kernel of mine survived any run; working tree untouched.
- Each WL kernel took ~8–10 min (heavier than the committed run due to sharing CPU with the other leg's seat).

## OVERALL: SOUND
The four harnesses faithfully implement the cleared 10-knife list as non-fakeable ablation certs:
live-engine-wrapped (independent digest reproduction), faithful full-payload digests (same-support
K_circular exhibit), all FORM knives bite / NO-OP=IDENTITY collapses to CANONICAL / DEAD sets carry no
nonzero leak, compact-contract clean including error paths, WL K_junk live and biting on the
`actualJunkCase` pin. Single non-blocking finding: **K_split_route ×2 coefficient companion is inert on
the pin** (its FORM knife bites); the arithmetic-sensitivity leg for that one knife is vacuous, though its
physics (FORM) sensitivity is genuine.
