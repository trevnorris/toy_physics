# S11c-c2 N6 comparator — fresh Opus review leg report (SOUND)

## VERDICT: SOUND
The astra-written N6 cross-engine comparator is a correct measurement instrument. It joins the two blind engines'
emitted objects on mechanically-typed keys, prints operand_A/operand_B/A_minus_B three-valued residuals, seals the
PIT/digest channels, and decides nothing. No false-agreement, false-disagreement, differenced-PIT, or leaked-value
path could be constructed. Every FORM ablation moved the output in the direction that proves the load-bearing object
is genuinely computed. Tests 23/23.

## Evidence (ablation scripts + stdout under scratchpad/opusleg_*)
- **FORM ablation #1 (residual genuinely computed):** `opusleg_ablation_spurion.py` — injecting `+OPUS_SPURION_FORM`
  into every parsed WL operand B: baseline 160 matched N6COV_R_COV `A_minus_B=Integer(0)` → ablated 0 still-zero, all
  160 moved off zero (80 carry the injected symbol). Proves A_minus_B depends on BOTH operands (comparator :816–824);
  the all-zero R_COV joins are real algebraic cancellations, not a designed-to-agree echo.
- **FORM ablation #2 (typed axis load-bearing):** `opusleg_ablation_wave_collapse.py` — collapsing the WAVE axis in
  make_key: join 160→0, duplicate_key 0→1536. No positional-guess / pre-collapse; block/kernel R_N6 families decode to
  disjoint axis sets → axis_set_mismatch/unmatched_key, never auto-joined (:164–169,:196–198).
- **PIT seal:** `sealed_field` (:623) seals primes/probe_num/den/numerator_denominator/sample_*; support derives from
  the boolean witness, never passed to residual(). 868 pit_sealed SURFACED lines, 0 cross-engine residual.
- **Carrier-rep residual computed (DO-NOT-FOLD):** N6RC_CARRIER_EULERIAN 320 joins → 280 zero + 40 genuinely nonzero
  symbolic residuals surfaced (not sealed/folded).
- **Support one-sided:** both NONZERO_WITNESSED/NO_NONZERO_FOUND appear; support_residual returns UndecidedResidual on
  a one-sided negative; never subtracts native booleans.
- **grade fold not a frozen axis (M3):** `{1,η,σ}→(η,σ)` drops a leading element constant=1 across all 288 WL columns
  (verified); grade[0]≠1 raises rather than silently merges (:110).
- **_NODES/join-set/DoD:** formal_jet srepr round-trip + missing-ARITHMETIC_DAG→parse_failed (:489; tests);
  SHARED set has R_COV_CONTROL_DELTA/SOURCE_CONTROL_DELTA/PHI_DOMAIN_CENSUS/ACTUAL_CONTROL_PARAMETERS (:53–58);
  MEASUREMENT_SCOPE residual_target=none; zero_extract_guard (:768) + exit-2 on parse_failed/zero_extract/all_deferred
  (:1042); join>0 never the exit.

## Non-blocking observations (correct behavior, no change)
1. dimension_record (:844) always KEY_DISAGREE (different schemas SymPy computed/consistent/target vs WL *_SUPPORT) —
   honestly surfaced; the discriminating signal is the sibling dimension_vectors channel. Redundant, not wrong.
2. addition_consistency (:845) always BooleanNotResidualable (PY emits a native bool) — the REQUIRED three-valued
   rejection of a native boolean, not a gap.
No banned per-case verdict token: `*_DISAGREE` are inherited Mismatch.kind descriptors inside the residual object;
`pairing ∈ {matched, sympy_only_column, wl_only_slot}` is descriptive accounting.
