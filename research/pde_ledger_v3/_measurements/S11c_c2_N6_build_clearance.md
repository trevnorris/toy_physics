# S11c-c2 N6 build — CLEARED (both build legs) + residual disposition (2026-09-06)

Companion `scripts/S11c_c2_N6_diagnostic_sympy.py` (918 lines, astra-built per the cleared directive + route-2 spec) +
the audit `ANCHORING_L_MINUS_M` edit (`a24acfb6`). Two build legs (Codex-family author → **fresh Claude agent + Grok**),
review-until-clear. Reports: `scratchpad/{grok_N6_buildreview.txt, tasks/adf923f88c524b2bc.output}`.

## Both legs: BUILD CLEAR (independent derivation + FORM ablation + one-sided route corruption + AST scan)
Convergent evidence (both legs, against /tmp copies; working tree untouched):
- **Independent derivation** (anchor scripts importing nothing from the artifact): reproduced the N4 advection mechanism
  (`μ̄(1)−μ̄(0)=B·g·u`; `g=0`⇒absence; the LAB_HELD e_W shift cancels from the t-difference). Artifact matches.
- **FORM ablation** (drop a δp slot / sign-flip E_W / collapse reverse-onto-forward): operand tables MOVE (Claude
  30→26/28; Grok 63–66/288 cols) ⇒ the increment is computed, ⛔ not a typed payload. Amplitude/reconstruction circuits
  unaffected (ablation localized to `build_increment`).
- **One-sided route corruption** (zero material μ only): Eulerian operand **BYTE-IDENTICAL**; Material + Residual MOVE;
  TILT (route-1) identical ⇒ the two routes are structurally INDEPENDENT; `R_N6` is genuinely `I_E − I_{M→E}`, not
  decoration.
- **N6-specifics all pass:** no `T`/no annihilation (`I_{M→E}` a live carrier); material μ at both sites, face-μ
  verifiably cancels from the δp carrier (`SLOT_GUARD`/`CLOSURE_GUARD` residuals = 0, matching spec §4); shared δp
  symbols (`PRESSURE_IDENTITY` all True); opaque c1 (no θ-sub/Jacobian into `DELTA_P`/`Z`/resolvent); **reverse blocks
  computed-ZERO** (grade-suppression holds in the actual output; the curl is NOT credited as an N4 witness); RHO4 emits
  **computed absence** (no A−A); tilt = injectable Eulerian factory PIT-checked vs the imported carrier (reconstruction
  residual 0), ⛔ not `FLIP_FACE_SLOPE`.
- **PIT sound:** primes {1000000009, 998244353, 1000000033} verified prime & ≡1 mod4; ≥8 draws (11–12); joint singular
  rejection; branch cells before modular reduction; on-shell `q` derived from `k`; honest FN bound with **max** per-prime
  (`conditional_good_prime`), δ≈5.7e-70; residual-zero ⛔ never an exit; **0 `assert`s** (AST scan). Dimensions:
  able-to-fail diagnostic yields nan on the deliberate bad operand; distinct weak dims, no energy-stamping.
- `reduction/derived_or_declared.py` failing to classify the emit = that legacy tool's colon-tagged-srepr grammar vs this
  diagnostic's purpose-built JSONL (no export/LEDGER write); the derived-vs-declared property is instead established by
  the FORM ablation. ⛔ NOT a defect (both legs).
- Line computing `R_N6`: diagnostic `residual(E,M)` at **:853**, `E,M` from `build_increment` **:850-851**.

## The residual DISPOSITION (the computed finding — ⛔ NOT pre-adjudicated)
Mechanical count of the emitted PIT `[num,den]` numerators (baseline outputs, LAB_HELD):
- **`S11CC2_REP_INVARIANCE_RESIDUAL`**: ~**18 certified-nonzero forward-block columns** (`THETA`/`E_W`,
  `TRANSVERSE_TO_THICKNESS`), both `RHOBR_CONSTANT` and `RHO4_CONSTANT`; reverse blocks computed-zero.
- Controls bite: `CONTROL_INDEPENDENCE_RESIDUAL[TILT]` ~34 cols; N4 advection moves (RHOBR live; RHO4 computed absence).

⭐ **Raw nonzero ≠ disagreement** ([[feedback_reconcile_representational_bridge]]): the raw `R_N6` residual is the
computed measurement, ⛔ NOT a verdict that representation invariance fails. Its disposition — a real invariance failure
vs a reconcilable remainder under justified identities — is the **reconcile / c2 step-record** stage, ⛔ not adjudicated
here. Per-engine (SymPy) only; the cross-engine (blind Wolfram) N6 is still owed.

## Owed next (unchanged from the c2 arc)
Adjudicate the R_N6 disposition at reconcile (⛔ don't pre-judge) → blind **Wolfram engine** N6 (imports nothing; 2
decision + 2 build legs) → **c2 T7 comparator + reconcile** (surface, ⛔ don't collapse) → **c2 step record**. Also
still owed program-wide: numeric F/G re-grounding; MEMORY.md compaction; the 499 MB reviewed `.out` is ephemeral.
