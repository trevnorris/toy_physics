# S11c-c2 N6 ablation harnesses — executable build review record

## Artifacts (four astra-built ablation harnesses + transcripts)
- `mathematica/S11c_c2_N6_ablation_harness.wl` + `scripts/S11c_c2_N6_ablation_harness_wl.py` →
  `_measurements/S11c_c2_N6_ablation_harness_wl.md`
- `scripts/S11c_c2_N6_covariance_ablation_harness.py` → `_measurements/S11c_c2_N6_covariance_ablation_harness.md`
- `scripts/S11c_c2_N6_diagnostic_ablation_harness.py` → `_measurements/S11c_c2_N6_diagnostic_ablation_harness.md`
- `scripts/S11c_c2_N6_reconcile_ablation_harness.py` → `_measurements/S11c_c2_N6_reconcile_ablation_harness.md`

Built by `gpt-6-astra` (high) against the cleared directive revision-2 (`b57420a8`). The reconcile harness's
first recorded run failed on a temp-tree provenance-copy setup fault; astra corrected the setup
(`copy_reconcile_tree`) and a follow-up astra run produced the complete reconcile transcript (11 labels, 143
computed digests, all subprocess exit 0). All four transcripts are compact (no heavy-content leak).

## Legs (Codex/astra-written scripts → fresh Opus agent + Grok)
- **Fresh Opus agent** — report `_legs/S11c_c2_N6_harness_build_review_agent_report.md`.
- **Grok** (`grok-4.6` high) — report `_legs/S11c_c2_N6_harness_build_review_grok_report.md`.
- Identical prompt `_legs/S11c_c2_N6_harness_build_review_prompt.md` (executable ablation, mandatory FORM
  ablation, both build-review nits carried in). Legs ran concurrently (2 Mathematica seats; memory-watched).

## Method (both legs ablated the harness, did NOT trust the transcript)
Both legs copied each harness to /tmp, re-ran its own worker/patch path per variant (fresh subprocess; WL
under `timeout --kill-after=5 600`, one kernel at a time), and compared their INDEPENDENTLY-computed
per-object `{sha256, nonzero fingerprint}` to the committed CANONICAL and FORM records; plus FORM-noop
ablation (neuter the knife → must revert to CANONICAL), one-sided engine corruption, and forced error paths.

## Verdict — BOTH LEGS SOUND (orchestrator concurs, G4)
- **Live-engine-wrapped, not faked.** Independent digest reproduction matched committed EXACTLY — 0
  mismatches across all four harnesses (WL 15/15 CANONICAL + per-knife FORM; covariance 7/7; diagnostic
  20/20; reconcile 13/13). A hardcoded/faked cert could not reproduce these.
- **All 10 knives BITE (FORM); NO-OP/IDENTITY reproduces CANONICAL exactly; FORM-noop reverts to CANONICAL**
  (proves the FORM application is load-bearing, not byte-identical-under-ablation).
- **WL K_junk BITES on the pin** `{MATERIAL_ADVECTED, RHO4_CONSTANT}` = `actualJunkCase` — the round-1 F1
  defect (knife silently disabled by the case pin) is RESOLVED in the built harness (runtime FORM moves
  N6COV_R_COV nz 0→1224, N6RC_R_N6 0→5184).
- **Compact contract clean, incl. error paths.** 0 occurrences of `numerator_denominator` / `ARITHMETIC` /
  `SparseArray` / `PROBE_NUMERATORS` / `srepr` in any transcript. A REAL forced/observed WL kernel SIGKILL
  emitted only bounded diagnostics (byte-count + SHA of stdout/stderr, `production_immutable=1`), no object
  payload. WL ~1 MB is appropriately-compact numeric fingerprints (object×cell×prime×variant integers +
  deduplicated component-key schemas), not over-verbose and not a leak.
- **Digest is a genuine full-payload SHA (same-support witness), not a re-hash of the fingerprint** —
  covariance K_circular/FORM moves all 5 certified SHAs while the nonzero fingerprint is unchanged.
- **DEAD / sampler-awareness (both engines): no false negatives, no real leaks.** Every DEAD "movement" is
  digest-only or `circuit_leaves`-only (shared joint-circuit metadata) or an honest MISSING column-set change
  (K_rank's shared-wave-basis drop) — the nonzero fingerprint (robust channel) never leaks into a DEAD
  object. Confirms the carried-forward nit: sampler-awareness applies on both the WL and SymPy samplers.

## Finding (single, NON-BLOCKING) — WL K_split_route ×2 coefficient companion is INERT
Both legs independently: the WL K_split_route **FORM** knife BITES (runtime SPLIT_CHECK nz 0→14688 via the
`affine False→True` route-flag; 10/10 match committed), IDENTITY=CANONICAL — but the **×2 companion moves
nothing** (SPLIT_CHECK digest_equal=1, no fingerprint delta).
- **Not a harness bug** — orchestrator mechanical check: `scripts/S11c_c2_N6_ablation_harness_wl.py:53`
  applies `SOURCE.replace('es[#] - ms[#]', '2 (es[#] - ms[#])')`, exactly the directive's ×2 spec; both legs
  confirm executed_sha ≠ budget_sha (the patch applied).
- **Mechanism (both legs agree):** under the canonical `affine=False` route the `sourceChannel` contributes
  ZERO to the printed split, so doubling `2×0` is inert; only the FORM route-flag admits a nonzero family.
- **Adjudication:** the FORM knife is genuinely implemented and bites structurally at the same site (proving
  the harness detects changes there); the ×2 inertness is a real physics property (vacuous arithmetic
  self-test for this one knife because the base contribution is zero), CONSISTENT with K_split_route being a
  FORM (family-leaving) knife rather than arithmetic. Documented; not a defect; re-opening the cleared knife
  design to alter the companion is unwarranted. NON-BLOCKING.

## Operational note (not a physics finding)
Running two Mathematica-ablating legs concurrently, one leg's broad `pkill -x WolframKernel` cleanup killed
the OTHER leg's live kernel (SIGKILL at ~411 s, not a timeout); both legs detected it and re-ran in idle
windows to completion, so it did not affect any verdict. Lesson for future concurrent WL reviews: legs must
kill kernels PATH-SPECIFICALLY (`pkill -f <that harness .wl>`), never `pkill -x WolframKernel`.

## Disposition
**SOUND → the four N6 ablation harnesses are the CLEARED ablation cert** for the S11c-c2 N6 knife list.
Stopping rule (G4): nothing outstanding changes what is computed or may be claimed. Commit the reviewed
harnesses + transcripts + reports. The K_split_route/×2 inert property is documented above (no fix).
