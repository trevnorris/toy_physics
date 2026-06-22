# Directive pathA_22b GATE 4 — derive the verdict ratio `g_mhat²/g_G` (action-level cancellation test) → value OR honest blocker

**Status:** READY-TO-EXECUTE (v3). Design-review `SOUND-WITH-FIXES` → 12 fixes applied → Codex confirm-pass `SOUND-AS-IS` →
GLM viability consult returned **(C)-leaning "likely BLOCKED, driven by `α_J`," + recommend a CONFIRMATORY run**. User chose
to run it confirmatory (2026-06-22). v3 folds GLM's framing nit: the `α_J`-doom hypothesis is now the FRONT-AND-CENTER thing
to falsify, and the `W_eff` (Route-a) cancellation is explicitly NECESSARY-BUT-NOT-SUFFICIENT (see §0 "GLM consult" + §2).
**Owner split:** Codex derives + codes + iterates; Claude (orchestrator) reviews and owns this directive + decisions/trackers. Codex does NOT touch git or `decisions/*`.
**Resume context:** `decisions/13` §0 + `decisions/14` §1c/§2/§5. Predecessor gates: pathA_22b Gate 0a (`MHAT_DIMENSIONFUL_CONFIRMED`), Gate 0b (`DOES_NOT_CANCEL (NOT_ESTABLISHED)`), Gate 1 (`Z_int`), Gate 2 (`λγ` = `BETA_GENUINE_GAP` → a genuine INPUT, not "still tunable"; `C_B/C_E=1`), **Gate 3 DONE (`χ_Q≈0.712`)**. This directive SUPERSEDES the one-paragraph GATE-4 stub in `directives/pathA_22b_minimal_combination_xi.md` §1.

---

## 0. Why this gate, framed honestly (read before designing)

The GR-quadrupole verdict is `P0·χ_Q·g_mhat²·λγ⁵/g_G ?=? 54/5`. Four of its factors are settled or classified:
`P0` (branch-determined, extracted), `χ_Q≈0.712` (Gate 3), `λγ` (genuine INPUT, `BETA_GENUINE_GAP`). **The only unsettled
piece is the gravity pair, and the verdict sees it ONLY as the single ratio `g_mhat²/g_G`** (GLM Trap 1 — `g_G` and `g_mhat`
are outputs of the SAME transverse-reduction computation; NEVER calibrate them independently, that would make the quadrupole
match spurious).

Two individual pieces inside that ratio are KNOWN to be underived:
- **`W_eff`** (the reduction width; `N_∞,3 = W_eff·ρ_∞,4`) — `W_EFF_REDUCTION_UNDERIVED`: the parent action defines
  projection `W(w)` (normalized `∫W dw=1`) and localization `Z(w)` but does NOT fix an invariant reduction width
  (canonical `pde.tex:277, 289–295, 496–565`). `Z(w)` (Maxwell kinetic weight) and `W(w)` (observation map) are DECLARED
  DISTINCT and must not be silently identified (pde.tex:277). [`research/pde_ledger/paper/parts/part01_parent_geometry.tex`
  is NON-canonical (historical) — it may be read for context but a derived value's provenance must rest on canonical
  `pde.tex`; if a needed reduction term is ONLY in part01 and not in `pde.tex`, that is itself a `BLOCKED_NEEDS_*` finding.]
- **`α_J`** (mass bridge, inflow `J`→`m_defect`) — `MASS_BRIDGE_FORM_NOT_DERIVED`: pathA_21 SEARCHED for an action-level /
  boundary-source / Noether / Hamiltonian map and found NONE "without restating the target" (pathA_21:28–33, 99–101).

So the direct route (compute `g_G` and `g_mhat` separately) is genuinely blocked on these two underived pieces. **The only
honest, possibly-winnable route — the one Gate 0b explicitly DEFERRED to Gate 4 — is to test whether these underived pieces
CANCEL in the ratio `g_mhat²/g_G`.** Gate 0b found the current sources do not already ESTABLISH a cancellation
(`NOT_ESTABLISHED`), but it is a source-availability verdict, not a derivation: "Before declaring `W_eff`/full transverse
profile irreducible, Gate 4 must test both cancellation routes at the action level" (pathA_22b report:34–35).

**This gate is therefore a DERIVATION with a genuine BLOCKED outcome on the table** — it is NOT a χ_Q-style numeric extraction,
and it is NOT permitted to manufacture `α_J` or `W_eff` by restating the target. A clean `BLOCKED` is a VALID, valuable result.

### GLM viability consult (2026-06-22) — the two doom hypotheses to FALSIFY (run this CONFIRMATORY)
A user-mediated adversarial GLM consult (source-checked against `pde.tex`) judged this gate **likely BLOCKED, driven by
`α_J`**, and recommended running it as a CONFIRMATORY falsification — not an open exploration. Treat these as the PRIMARY
hypotheses; the gate's job is to FALSIFY them with a genuine target-blind action-level derivation, or CONFIRM them as a
*derived* blocked (which hardens the verdict far beyond Gate-0b's source-availability punt):
- **H-α_J (the fatal one): `α_J` does NOT cancel.** `α_J` enters `g_G` via `G_cond ∝ 1/(α_J1·α_J2)` (the stress/gravity lane).
  `g_mhat` is the source-map normalization (`m̂0` = dimensionful scale carrier, `pde.tex:2081`; `m̂=1+O(a²/r²)` dimensionless
  profile, `pde.tex:2097–2099`). GLM reads the action as placing NOTHING carrying `α_J` in the `g_mhat` lane ⇒ `α_J` survives
  in the ratio ⇒ `BLOCKED_NEEDS_alpha_J`. **The ONE genuinely open sub-question that could falsify H-α_J:** does the `m̂0`
  source-map lane *inherit* the inflow/mass-bridge normalization (because the thing it sources is the mass, and
  `m_defect = α_J·ℏ·J/c²`)? Step 4-i must derive `K_source` target-blind and SEE whether it carries that factor — do NOT
  assume either way.
- **H-K_source (strong secondary): `K_source` may not EXIST as an independent target-blind kernel.** The "source-map terms"
  may be only the `Γ5`/Burke–Thorne-facing normalization (which feeds the target), not a genuine transverse kernel derivable
  independently of `54/5`. If so → `BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE` at Step 4-i, before the cancellation test even runs.
- **Independence of blockers (binding on the overall verdict):** H-α_J and the `W_eff` reduction are INDEPENDENT. The `W_eff`
  Route-(a) cancellation is the tractable-feeling part GLM warns we may over-weight — it is **NECESSARY BUT NOT SUFFICIENT**.
  An overall `RATIO_DERIVED` (a real WIN) requires that `K_source` exists target-blind (falsify H-K_source) AND `α_J` cancels
  (falsify H-α_J) AND the kernel functional cancels. If `W_eff` cancels but `α_J` does not, the outcome is STILL
  `BLOCKED_NEEDS_alpha_J` — a `W_eff` cancellation alone MUST NOT be reported as a win.

---

## 1. The target + the structural question (target-blind)

Definitions (carry symbolically; do NOT set any `g_*`, `Θ_Q`, `α_J`, `W_eff` to 1):
```
G   = (a·c_s²/m_GNLS)·g_G          ⟹  g_G    = G·m_GNLS/(a·c_s²)
m̂0  = (c_s/(a²·√m_GNLS))·g_mhat    ⟹  g_mhat = m̂0·a²·√m_GNLS/c_s
```
The conditional `g_G` chain inherited from pathA_21 (cite when used; do NOT re-derive its target-restating steps):
```
G_cond = c_γ⁴·m_GNLS·Θ_Q1·Θ_Q2·I_F,12 / (4π·N_∞,3·α_J1·α_J2·ℏ²),   N_∞,3 = W_eff·ρ_∞,4
```
Note the gravity coupling for the SINGLE radiating defect in the verdict — reduce the two-body `G_cond` to the single-defect
normalization that the `g_G` definition above actually carries. **Required: a written mini-lemma** that tracks the exact
powers of `Θ_Q1·Θ_Q2`, `α_J1·α_J2`, `N_∞,3`, and `I_F,12` from the two-body `G_cond` down to the scalar `g_G`, justifying
each step from the action/definition. FORBIDDEN: square-rooting the two-body expression or dropping one defect's factors "by
convention" to land a clean single-defect form. If the single-defect normalization is not uniquely fixed by the
sources/definition, that is a `BLOCKED_NEEDS_*` finding, not a free choice.

**THE STRUCTURAL QUESTION (the whole gate):** Is the dimensionless ratio `g_mhat²/g_G` independent of BOTH
(i) the transverse reduction kernel (`W_eff` / the `W(w)`,`Z(w)` profile content), AND
(ii) the mass-bridge factor `α_J`?

If yes → the underived pieces cancel and the ratio is a DERIVED O(1) number (or closed form in already-fixed inputs).
If either survives → the ratio is genuinely blocked on that piece.

## 2. The work (requirement + acceptance — Codex designs the route; do NOT pre-script it here)

**Step 4-i — derive the two transverse kernels from the parent action (the load-bearing analytic step).**
From the parent action's stress / source / reduction structure, derive:
- `K_stress(w)` — the transverse kernel carried by `g_G` (the defect-stress→long-range-coupling lane; the Noether-stress /
  reduction content that produces `N_∞,3`, `Θ_Q`, `α_J`).
- `K_source(w)` — the transverse kernel carried by `g_mhat` (the source-map lane; how `m̂0` sources the field, with Gate-0a's
  `m̂=m̂0·g_mhat`, `m̂0` dimensionful reading).
Each kernel: full provenance line-refs, an explicit statement of what is integrated over `w` and against which weight
(`W(w)` projection vs `Z(w)` localization — keep them distinct per pde.tex:277), and a restored-units dimensional check.
Where a kernel piece is itself underived (e.g. `α_J`), carry it as an explicit SYMBOL — do NOT assign it a value.

**BURKE–THORNE FIREWALL (binding, target-blindness).** The `m̂²·Γ5 = 2G/(5c⁵)` Burke–Thorne normalization equations
(`pde.tex:2075–2099`) are TARGET-FACING — they may be used ONLY for dimensional inheritance and the FINAL assembly (Gate 5),
**never to derive `K_source`, `g_mhat`, or the ratio.** Deriving `K_source` from those equations would manufacture it from
the very quantity (`G`/the BT coefficient) we are trying to predict. The `K_source`/source-map STRUCTURE must come from the
action's source-map terms (`pde.tex:2053–2069` and the parent source/current terms), NOT the BT-normalization block.

**STEP 4-i STOP / failure outcomes (before any cancellation algebra).** If `K_stress` or `K_source` CANNOT be derived
target-blind from the action/source content (i.e. only a target-facing or convention-dependent form exists), STOP and report
`BLOCKED_NEEDS_K_STRESS` / `BLOCKED_NEEDS_K_SOURCE` / `BLOCKED_NEEDS_SOURCE_MAP_PROVENANCE`. Do NOT run Step 4-ii on an
undefined or target-derived kernel — a cancellation test between illegitimate kernels proves nothing.

**Step 4-ii — the decisive cancellation test (must be fail-able + non-tautological).**
Test, at the action level, whether `g_mhat²/g_G` is independent of the reduction kernel AND of `α_J`, via the two routes:
- **Route (a) — shared factored scalar:** both `g_G` and `g_mhat²` are the SAME scalar transverse functional times SEPARATE
  `w`-independent field-content factors ⇒ the shared functional cancels in the ratio. Use the PDE-defined measure
  `Z_int = ∫Z(w)dw` (`pde.tex:289–295`) as the default shared functional; do NOT introduce a `√g_w` measure unless Step 4-i
  actually DERIVES it from the action/export (Gate 1 already flagged `√g_w` as not independently exported).
- **Route (b) — pointwise-proportional kernels:** `K_stress(w)=const·K_source(w)` for all `w` (equivalently
  `K_stress(w)K_source(v)=K_stress(v)K_source(w)` ∀ `w,v`) ⇒ the weighted-average ratio is profile-independent.
Also test explicitly whether `α_J` cancels. **`α_J` cancellation rule (binding):** a cancellation is valid ONLY if the SAME,
independently-derived bridge factor demonstrably appears in BOTH lanes with the powers required to divide out — and that
factor's presence in `g_mhat` must be SOURCED from an action / boundary / Hamiltonian relation that is independent of the
`54/5`/BT target. Asserting "`g_mhat` carries `α_J`" by analogy, convenience, or restatement is FORBIDDEN (that is the exact
pathA_21 trap). If `α_J` appears only in `g_G` (the stress lane) and there is no independent reason it sits in `g_mhat`, then
it does NOT cancel → `DOES_NOT_CANCEL`/`BLOCKED_NEEDS_alpha_J`.
**MANDATORY non-tautology guards** ([[feedback-decisive-test-not-tautological]]):
- The two kernels MUST be built from INDEPENDENT provenance (separate derivations / ASTs) BEFORE the comparator sees them —
  never the same symbolic object reused for both. The comparator is applied to independently-derived objects only.
- POSITIVE-comparator control: the negative control is a deliberately non-factorizing / non-proportional kernel pair that
  MUST return `DOES_NOT_CANCEL`. ADDITIONALLY a MUTATED-KERNEL control: perturb one real kernel and confirm the comparator
  then returns a NONZERO cancellation residual (proves the comparator measures the kernels, not a hardcoded pass). A synthetic
  pair alone only proves the comparator works — it does not prove the real derivation.
- The `CANCELS` result must NOT be a consequence of having set the kernels equal by construction, nor of any step that used
  the `54/5` target. State, in prose, the exact structural condition that makes it cancel and why the action forces it.

**Overall-verdict logic (binding — guards against the `W_eff`-cancellation red herring).** Proceeding to a Step-4-iii
`RATIO_DERIVED` (a real WIN) requires ALL THREE: (1) H-K_source FALSIFIED (a genuine target-blind `K_source` exists), AND
(2) H-α_J FALSIFIED (`α_J` provably cancels via the binding α_J rule above, OR is provably absent from `g_G` in the
single-defect reduction), AND (3) the shared kernel functional cancels (Route a) or the kernels are pointwise-proportional
(Route b). If ANY of the three fails, the outcome is the corresponding `BLOCKED_NEEDS_*` — in particular a `W_eff`/Route-(a)
cancellation with `α_J` still present is `BLOCKED_NEEDS_alpha_J`, NOT a win.

**STOP HERE for orchestrator review** (this is the decisive fork). Report the Step-4-i kernels (incl. the H-α_J / H-K_source
findings), the Route-a/Route-b cancellation results, and the overall verdict per the logic above, before doing Step 4-iii.

**Step 4-iii — resolve to a value or a blocker (after the orchestrator gate).**
- If `CANCELS`: derive the residual kernel-free ratio `g_mhat²/g_G` as a closed form / O(1) number in already-fixed inputs.
  Dimensional-check it is dimensionless. Outcome `RATIO_DERIVED_<value>`.
- If `DOES_NOT_CANCEL`: the `RATIO_DERIVED_VIA_WEFF` branch is reachable ONLY under ALL of these conditions (else it is NOT
  reachable — return `BLOCKED_NEEDS_*`):
  (1) the residual functional depends SOLELY on quantities actually present in the frozen export
  (`psi`, `A`, `rho0`, `R0_w`, `Z_w`, grid in `frozen/m1c/…/wp1_background_10x8.json`) — NOT on an arbitrary `W(w)`, an
  undetermined `α_J`, or an assumed infinite-tail continuation;
  (2) the `localization_floor=0.8` issue is honestly resolved — a finite-box `Z_w` with a NONZERO floor is NOT an ideal
  `(-∞,+∞)` reduction (`pde.tex` uses `Z(w)` over `(-∞,+∞)`); a bounded systematic requires EITHER a source-derived finite
  `w`-extent OR a decaying/zero-floor continuation with provenance — an undocumented config floor (`coupled_branch.py:188–192`)
  does not qualify;
  (3) the result carries an explicit `±systematic` from the finite-domain/floor uncertainty.
  Only then → `RATIO_DERIVED_VIA_WEFF_<value> (±systematic)`.
  If the residual piece (`W_eff` kernel and/or `α_J`) is genuinely underivable from the action AND uncomputable from the
  frozen data without an arbitrary choice → `BLOCKED_NEEDS_<W_eff | alpha_J | reduction_kernel>` and STOP (do NOT fake a
  ratio; do NOT calibrate `g_G`/`g_mhat` independently to backfill — that is the GLM-Trap-1 forbidden move).

## 3. Method + discipline (binding)
- **TARGET-BLIND:** derive with NO upstream reference to any target-facing constant — banned upstream of Gate 5:
  `54/5`, `10.8`, `10.8/P0`, `P0_target`, `N_Q`, the outgoing-BT-target / outgoing-factorized-normalization, the
  `pde.tex:2075–2099` BT block, and any imported target constant. Target equations enter ONLY after the ratio is derived or
  blocked. No factor may be reverse-engineered from the target ([[feedback-calibrate-predict-methodology]]).
- **DERIVED-FORM GATE:** no `x==x` posing as a check; no hand-inserted field/coefficient/sign; no `g_*`/`Θ_Q`/`α_J`/`W_eff`
  set to 1 unless DERIVED; no restatement faking closure (the pathA_21 `α_J` trap — do NOT re-commit it). An honest
  `BLOCKED`/`DOES_NOT_CANCEL` is a VALID outcome.
- **Dimensional-check** with units restored (SymPy `M,L,T`) on each kernel, the ratio, and (if reached) the value
  ([[feedback-dimensional-consistency-check]]). Recall `m̂0` is dimensionful (Gate 0a).
- **Dual-engine:** `.wl` via `math -script` (NOT `wolframscript` → exits 255) wherever Mathematica can independently verify
  the symbolic algebra / cancellation ([[feedback-dual-engine-required]]). ≤2 concurrent `math -script`.
- **Transliteration fidelity:** the derived kernels must match the source action term-by-term; a faithful-but-wrong operator
  is the failure mode MMS cannot catch ([[feedback-transliteration-fidelity-audit]]).
- `timeout 600` on every SymPy/`.wl` script (this is symbolic derivation, not a solver run; if anything threatens 600s it is
  mis-scoped — reformulate, do not raise the cap) ([[feedback-script-timeout-policy]]).
- Additive: NEW module(s)/group(s) (e.g. `patha22b_gate4.py` + `.wl`); do NOT modify faithful operators or the
  pathA_18/19/20/20b/21/22a or Gate-0..3 groups. Sources: `research/pde/paper/pde.tex` (NOT `pde_ledger/paper`),
  pathA_21/22a/22b reports, `src/stage1_solver/{patha22b_gate0,patha_extraction,m1c_background_export}.py`.
- Do NOT commit; do NOT touch git or `decisions/*`. Leave outputs for review.

## 4. Deliverables
- A Gate-4 section appended to `reports/pathA_22b_minimal_combination_xi.md`: the two derived kernels (with provenance
  line-refs + dimensional checks), the target-blind attestation, the cancellation test + its negative control, and the gate
  outcome (§2 classes) + residual ledger (what, if anything, remains blocked and the minimal next step to close it).
- New SymPy module/group + `.wl` cross-check; tests including (1) a target-blindness guard (no banned target constant is read
  upstream of assembly), (2) the negative-control test (non-factorizing kernel pair ⇒ `DOES_NOT_CANCEL`), (3) a
  mutated-kernel control (perturbing one real kernel ⇒ the comparator returns a NONZERO cancellation residual), and (4) an
  independence check that `K_stress` and `K_source` are built from separate provenance (not the same reused symbolic object).
- Run logs; do NOT commit.

## 5. Review plan (after Codex exit 0)
Two CLEAN agents ([[feedback-review-agents]], [[feedback-offload-to-agents]]):
1. **Dimensional-fidelity audit** — each derived kernel vs the parent-action source, term by term; dimensions restored;
   `m̂0` dimensionful; `W(w)`/`Z(w)` not conflated.
2. **Adversarial review** — is the cancellation a REAL proof or assumed? Was `α_J`/`W_eff` restated to manufacture closure?
   Is the negative control genuine (can the test actually return `DOES_NOT_CANCEL`)? Any `g_*`-set-to-1? Any `54/5` leakage?
   Distrust an all-clean `CANCELS` result hardest ([[feedback-negative-verdict-short-circuit]] applied symmetrically — a
   too-convenient WIN needs the same scrutiny as a blocker).
Route the DECISIVE cancellation test (Step 4-ii) through the GLM tertiary consult before accepting `CANCELS`
([[feedback-decisive-test-not-tautological]] — GLM previously caught a tautological "decisive" test that Claude+Codex passed).
Claude synthesizes; bring the gate outcome to the user. If `RATIO_DERIVED`, flag the GATE-A freeze implication (a derived
value un-pins `m̂0²·S_port`) as a user methodology call before Gate 5 assembly.
