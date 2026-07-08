# II-B3 (ledger_stage010) source map — pathA_29 Check B: slab localization p=2 + NOGO

> Running-start prep captured 2026-07-08 from the stage009 build session (all line refs verified in-session) so a
> fresh session can author the reshape directive without re-discovery. Verify against the cited sources before
> finalizing. Companion: `part2_gravity_atomic_split.md` (row 010) + `stages/ledger_stage009_flat_slab_return_residual.md`
> (the Check-A stage this completes). Build-order id **010**, Part II.
> **Reshape of the Check-B slice** of `software/stage1_solver/tools/pathA_29_brane_bulk_return_sympy.py` (1257 L)
> + a re-authored independent `.wl` (the source `.wl` is digest-bridged — same severing as stage009).
> **Source top-line (verbatim): `RETURN_RESIDUAL_PREDICTION`** — the JOINT Check-A+Check-B verdict. Stage 009
> computed the Check-A component; **stage 010 computes Check B and is the COMPLETING stage** — it consumes 009's
> Check-A classification as a cited input and emits the joint headline (label the composition honestly).

## File inventory
- **Report:** `software/stage1_solver/reports/pathA_29_brane_bulk_return.md` (29 lines; Check-B content :11–17,
  :20–25). Machine ledger `reports/pathA_29_results.yaml`.
- **The Check-B slice of the `.py`:** `solve_dynamic_zero_mode_radial` :164–212, `solve_static_zero_mode_radial`
  :213–250, `solve_static_massive_radial` :251–286, `build_counterfactual_guard` :287–304,
  `classify_dc_sink_gate` :305–316, `branch_verdict_from_p` :317–320, the Check-B part of `compute_symbolics`
  (≈:397–700: compact-cell/Bloch eigenfunctions `f0_abs = 1/√d`, `f1_abs = √(2/d)·cos(πw/d)`, m=0 eigenvalues,
  the NOGO warp construction), `radial_decay_exponent` :127, `radial_operator_residual` :135.
- **Stage-009 territory (do NOT recompute):** `solve_round_trip_phase`, the T_ℓ/ν/p_res/Z/strict-limit slice
  (:321–397), the A-classification :485–486 — consume the Check-A component as a CITED input.

## §1 The Check-B content (report :11–17, :24–25, verbatim facts)
- "Check B was run only on the **admissible DC-sink completions**":
  - `destructuring_absorbing`: solved compact-cell spectrum + a branch-specific 3D radial equation → **derived p=2**;
  - `bloch_stack`: solved q=0 Bloch spectrum + a separate 3D radial equation → **derived p=2**.
- **Counterfactual guard:** "multiplied the solved static zero-mode Green function by r⁻⁴, changing 1/r to 1/r⁵" —
  REJECTED with residual `5/(π·d·r⁷)` (the guard proves the radial solve has teeth; keep it).
- **The mandatory NOGO control:** the anti-localizing half-line warp `μ(w) = exp(2·k_warp·w)`; its zero mode is
  **non-normalizable**; the continuum Green integral gives **p=3**; the SAME classifier returns **`RETURN_NOGO`**
  (proves the classifier able-to-fail — a delocalizing return genuinely kills 1/r² gravity).
- The radiation/Sommerfeld boundary is recorded **`ac_check_a_only`** — NOT a Check-B branch (print as provenance).
- **pde_ledger feed:** "the gravity-range item passes inside the localizing flat-slab family because both DC-sink
  completions give p=2" — 1/r² gravity SURVIVES the slab. Open-item #9 stays not-closed (009 carries the residual).

## §2 The solve mechanics (verified from the `.py`)
- **Dynamic route** (:164–212): seeded by the COMPUTED m=0 transverse eigenvalue (raises if `m_value ≠ 0`); solves
  `g'' + (2/r)g' + ((ω/c_s)² − m²)g = 0` by genuine `dsolve` FIRST, selects the outgoing spherical branch via
  explicit C1/C2 substitution (a **boundary SELECTION** — "normalization fixed by the compact zero-mode overlap
  1/d"; label it as a selection, print the basis), asserts the operator residual ≡ 0, THEN takes `ω → 0` → the
  limit Green function; flow = −d(Green)/dr; exponent from `radial_decay_exponent`.
- **Static route** (:213–250): solves the static radial equation directly. **The static–dynamic CONSISTENCY check
  (the two routes agree) is load-bearing content — KEEP it as a computed check; DROP only the SHA-256 trace-id
  plumbing** (the report's :20–22 "separate traces" hashes = build-reproducibility bookkeeping, not physics).
- **Massive/gapped control** (:251–286): the massive radial solve (Yukawa-class falloff) — the contrast that makes
  p=2 meaningful; keep.
- **Classifier** (:305–320): `classify_dc_sink_gate(branch_ps, quadrupole_survives)` + `branch_verdict_from_p` —
  the verdict is COMPUTED from the derived branch p's; `p=2` both → `RETURN_RESIDUAL_PREDICTION` (with Check A);
  `p=3` → `RETURN_NOGO`. The reshape must keep this a function of the computed exponents (no stamps) and make the
  `quadrupole_survives` input's provenance explicit (cited from stage008's p_raw(ℓ2)=5 export, NOT recomputed ℓ=2).

## §3 Reshape trip-ups (pin in the directive)
1. **Same bridge-severing class as stage009:** strip the `.py`'s JSON emit + expression digest (the `.wl` Imports
   + matches), the runtime YAML reads (`assert_patha28_reuse` reads pathA_28 results; also its own results-yaml
   WRITES + report writer + `validate_results` + `compare_mathematica_if_available`), ALL `sha256_text`/`trace_id`
   bookkeeping. ZERO file I/O in the stage scripts.
2. **The dsolves must stay genuine** — the counterfactual guard (r⁻⁴ rider → residual `5/(πdr⁷)` fires) is the
   existing tooth proving it; keep + verify it fires. No hardcoded 1/r Green anywhere.
3. **Joint-headline honesty:** stage010 emits `RETURN_RESIDUAL_PREDICTION` as the COMPLETED verdict = (Check-A
   component, CITED from stage009 with a citation-integrity tie) + (Check-B p=2, computed here). Print the
   composition; the NOGO branch (classifier → `RETURN_NOGO`) must remain reachable (the warp control proves it).
4. **Boundary-selection labeling:** the C1/C2 outgoing-branch choice + the 1/d overlap normalization are
   SELECTIONS on the postulated slab family — label, don't present as derived.
5. **Do NOT recompute Check A** (009's τ/T_ℓ/ν/p_res/Z) or ℓ=2 content (`T2_applied=false` carried; the
   `quadrupole_survives` input cited from stage008).
6. **`.wl` independence:** the digest bridge dies; the `.wl` solves the radial ODEs natively (`DSolve`), builds its
   own spectra/normalizability checks (`Integrate` for the zero-mode norm; the NOGO non-normalizability must be a
   computed divergence, not a flag), own classifier plumbing. ⚠ Arity discipline + arity self-check block
   (stage008/009 pattern — a def/call arity mismatch silently skips at exit 0).
7. **Dual-corruption anchoring:** every check anchored in-engine to derivations or frozen expected values (the
   stage006/009 lesson — set-then-compare-to-self integrity blocks are the recurring rig class; stage009's
   consumed-kernel block was CAUGHT as exactly this and fixed dual-site. Any consumed value from 008/009 gets the
   dual-site pattern from the start).
8. **W_slab caveat linkage:** the slab width `d` is a register row (stage009, `ACTION` geometry); localization
   holding for the FAMILY ≠ the family being selected by dynamics (R19/`W_slab`, the old L/a item) — print the
   caveat, don't upgrade.

## §4 Teeth candidates (the directive finalizes)
Keep/port: the r⁻⁴ counterfactual guard; NOGO-warp-as-control (classifier flips to `RETURN_NOGO` — able-to-fail
proven); corrupt a radial-operator coefficient (2/r → 3/r) → dsolve residual/exponent asserts fire; corrupt the
m=0 seed (m≠0) → the dynamic solve's guard raises; corrupt the zero-mode normalization → static–dynamic
consistency fires; non-normalizability tooth (make the warp normalizable → the NOGO leg's divergence check fires
the other way); consumed-input corruption (009's Check-A component; 008's p_raw(ℓ2)) → citation-integrity fires;
classification corruption (feed p=3 into the baseline path) → headline assert fires.

## §5 Downstream consumers
- **Stage 023 (pathA_34):** consumes the residual forms + `Z_is_premise` (via 009) — 010's p=2 underpins the
  gravity-range claim it non-regresses.
- **Stage 024 (pathA_43):** consumes the `Φ_ℓ(w,r)` bulk Helmholtz mode + the projected-continuity operator
  lineage — 010 is where the bulk mode's localization is formally earned in the v2 ledger.
- Register: stage010 likely adds NO new knobs (`k_warp` is a control-construction parameter, tracked-not-counted;
  the classifier verdict discharges nothing — R23 stays `PENDING`; note the R19/`W_slab` caveat linkage).

## Verdict tokens + honest scope
Headline: **`RETURN_RESIDUAL_PREDICTION`** (joint, completed here) with `RETURN_NOGO` genuinely reachable.
EARNED: p=2 from real dsolves on BOTH admissible completions (+ the massive contrast), the counterfactual guard,
the NOGO control, static–dynamic consistency. POSTULATED/SELECTED: the slab family (009's rows), the outgoing
boundary selection + 1/d normalization. CITED: 009's Check-A component; 008's `p_raw(ℓ2)=5`. Does NOT close
open-item #9; does NOT select the slab width (R19).

## Process (unchanged, calibrated)
Author reshape directive (§3 pins + §1/§2 faithful cover) → Codex xhigh design-review → fold to `DIRECTIVE_CLEAN`
(no GLM on Parts I–VI) → **pre-exec USER GATE** → Codex builds the two scripts (`--sandbox danger-full-access`,
background, `< /dev/null`, xhigh) → dual-engine exit 0 (repo root + foreign CWD) → arbiter re-run via runners →
tri-review on fresh agents (fidelity + adversarial-with-ablation, incl. the `.wl` arity scan + tally-inflation
spot-check + dual-corruption class) → remediate → fresh-agent re-verify → registration 9→10 + parameter register
update + Codex-verify → note/card/`\input{stages/stage_010}` (Part-II appendix) → PDF → commit + docs/memory sync.
Target stem: `ledger_stage010_slab_localization_p2_nogo` (confirm slug at authoring).
**Codex sessions: `codex exec -c model_reasoning_effort=xhigh`, backgrounded with `< /dev/null`, absolute paths
(a relative-path launch failed once this session), never wrapped in shell timeout.**
