# _measurements — S11c_c1_wl_repair2_directive.md (the FULLY-GATED repair-2 + its 2-leg re-review)

Repair-2 fixes the two propagated emit-only defects the retroactive gate confirmed, **this time with its rule-7
decision legs run BEFORE the build** (record `_measurements/S11c_c1_wl_repair2_directive_review.md`: both legs "do
not build", 3 MUST-FIX including that the first-draft R2 targeted a non-defect — the gate caught it before any
compute). Directive `directives/S11c_c1_wl_repair2_directive.md`; baseline `13f0bd2c`. The `.wl` is Codex-written,
so its re-review legs are **a fresh Claude Agent + Grok** (authorship table), **serialized** (both ablate
Mathematica; 2-seat licence). Re-review prompt: `directives/_legs/S11c_c1_wl_repair2_review.md`.

## Commands (literal)
```bash
# Build (Codex, detached; setsid+marker+Monitor; danger-full-access)
codex exec -c model_reasoning_effort=xhigh --sandbox danger-full-access "$(<…/S11c_c1_wl_repair2_directive.md)" \
  > …/scratchpad/codex_repair2_build.log 2>&1        # EXIT=0; 167,647 tokens
# Re-review leg 1 — fresh Claude Agent (in-process, general-purpose), serialized FIRST
#   (ran wolframscript foreground with timeout 600; artifacts under /tmp/s11cc1_leg/)
# Re-review leg 2 — Grok, serialized SECOND
grok --prompt-file …/_legs/S11c_c1_wl_repair2_review.md --cwd /var/projects/toy_physics \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain \
  > …/scratchpad/grok_repair2_review.log 2>&1        # exit 0
```

## Deliverable verification (orchestrator, rule 13; vs baseline `13f0bd2c`)
- `git diff --stat 13f0bd2c -- <the .wl>` = **+144/−45**, one file; four hunks only:
  `operatorCompositionRecordFromDerivation`/`fredholmFunctionSpaceOperator` (R1), the `NONINVERTIBILITY_CONDITION`
  payload (R1), the memory-parity helpers (R2), the `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` emission (R2).
- Blindness: structural scan clean (0 `Get`/`Import`/`<<`/`ReadString`/`OpenRead`/abs-path); `git diff --check`
  clean; tracked `mathematica/out/` unchanged. Builder in-repo vs isolated-copy runs **byte-identical** (both 51
  tags / 81,881,078 bytes; SHA-256 `d63e616786eebb…`), isolated dir retained exactly its single copied `.wl`.
- Runs: in-repo exit 0 (4m19s), isolated exit 0 (4m21s); peak RSS ~0.8 GiB (light). Tag count **51** unchanged.

## The two repaired controls — each now BITES (BOTH legs, literal one-sided corruption)
- **R1 — Fredholm `Z` is the two-leg DtN operator.** `fredholmFunctionSpaceOperator` delegates to
  `operatorCompositionRecordFromDerivation` (both legs). Base `FREDHOLM_Z_HAS_Q_INPUT = PROVED_TRUE (count 2)`,
  `_Q_OUTPUT = PROVED_TRUE (count 3)`; **constructor-level re-freeze** (`rightLegMomentum → momentumOutput`, DTN
  live left on `momentumInput`) → `FREDHOLM_Z_HAS_Q_INPUT = PROVED_FALSE (count 0)` and `FREDHOLM_SAME_AS_DTN_LIVE
  → False`, while `DTN_OPERATOR COMPOSITION_HAS_Q_INPUT` stays `2` — a genuine constructor change (a payload-wide
  `qOut[…momentumInput]` sprinkle would not split Fredholm from DTN nor move the expanded `q_out(k)q_out(k′)`
  product). `momentumOutput`/`momentumInput` are distinct vectors (`{kOne…}`/`{kPrimeOne…}`), not a shared dummy;
  the flat loci stay the single-`k` symbol (not a `k=k′` slice). Leg-1 `/tmp/s11cc1_leg/r1_ablation.{wl,out}`;
  Grok `PROBE_R1_*`.
- **R2 — dissipation per parity COMBINATION.** `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` now emits the parity-combination
  memory Hermitian via the same `portParityCombination` change of basis. Base `DELTA_W ≠ ZETA_C` (distinct — not
  the pre-repair dead axis); a **one-sided `+`-face memory-Hermitian sign-flip at construction** moves the two
  parity payloads **differently** (`MOVES_IDENTICAL = False`, both individually nonzero), `−`-face unchanged,
  protected tags unchanged. Independent recon of the congruence vs the emitted blocks = `{{0,0},{0,0}}`. Leg-1
  `/tmp/s11cc1_leg/r2_ablation.{wl,stdout}`; Grok `PROBE_R2_*`.

## `PERMEABLE_PORT_HERMITIAN` — INDEPENDENTLY re-adjudicated by BOTH legs (left byte-identical; correct)
Both legs derived, from `V_s = s∂_tζ_s` (`S11c_a:54-59`), that its diagonal blocks are the **congruence**
`Aᵀ·diag(P₊,P₋)·A` (`A = faceToParityMatrix = [[1/2,1],[1/2,-1]]`): both diagonal blocks EVEN
(`ζ_c → P₊+P₋`, `δW → (P₊+P₋)/4`), the odd combination `(P₊−P₋)/2` the coupling. Confirmed on the real artifact
(`baseZC == P₊+P₋`, `baseDW == (P₊+P₋)/4`). The naive `δW↔(ζ₊−ζ₋)` reading is the wrong congruence direction (it
would vanish the thickness port at equal faces). ⇒ the retroactive "parity swap" finding is refuted; keeping
`PERMEABLE_PORT_HERMITIAN` byte-identical is correct.

## Protected core byte-identical (both legs)
Only 3 tags differ vs `13f0bd2c`: `NONINVERTIBILITY_CONDITION` (R1), `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (R2), and
`LOCAL_TAG_NAMES` (bookkeeping — dropped the old Fredholm Module locals, added the R1/R2 helper symbols). All other
48 tags IDENTICAL: `DTN_*` (12), `DTN_FLAT_SYMBOL`, `FACE_RESPONSE`/`_COEFFS`, `PERMEABLE_PORT_HERMITIAN`,
`DEGENERATE_LOCI_*`, `ENERGY_*`, `REP_INVARIANCE_*`, all `CONTROL_*`, `UNIFORM_LIMIT_*`, `ZERO_JET_*`,
`HOMOGENEITY_*`, `DIMENSIONS`, T-a..T-i. Source-identical constructors: `directBoundaryDerivation`,
`operatorCompositionRecordFromDerivation`, `portParityCombination`, `curvedPortResponse`, `faceToParityMatrix`.

## Rule 5 / independence (both legs)
No `Assert`; no stated expected residual; no "residual is zero for the correct labelling." The new identity tests
are `relationalObject[DELTA_W, ZETA_C]` / `relationalObject[move_DW, move_ZC]` — operands, no value gate. R1 freeze
is not `A−A` (Fredholm moves, DTN live does not). Blindness intact (both runs byte-identical isolated-vs-in-repo).

## NITs (both benign; ⛔ NOT MUST-FIX — carry to the c1 step record, rules 10/15)
- **NIT-1 (bookkeeping):** `LOCAL_TAG_NAMES` differs (the expected consequence of the repair's symbol churn; not a
  physics operand).
- **NIT-2 (correct physics):** the two diagonal parity Hermitian blocks are **proportional** (`ζ_c = 4·δW`)
  because both are the EVEN combination `P₊+P₋`; the genuinely independent parity content is the off-diagonal
  coupling `(P₊−P₋)/2` (which IS emitted, `DELTA_W_TO_ZETA_C_BLOCK`). The R2 corruption-move test still
  distinguishes repaired from the dead-axis bug (which gave byte-identical moves), so the control is adequate — but
  the "distinctness" of the two diagonal payloads is a fixed `1/4` vs `1` kinematic ratio, not structural
  independence. Both legs found this independently.

## Verdict
Repair-2 CLEARED — **both re-review legs "no MUST-FIX"**; both repaired controls bite by independent one-sided
corruption; protected core byte-identical; `PERMEABLE_PORT_HERMITIAN` independently re-derived correct. This closes
Item A of the remediation (`_measurements/S11c_c1_wl_remediation_plan.md`). NEXT in the remediation: both `.out`
(Item C), state-memory correction (Item D), then the T7 cross-engine comparator.
