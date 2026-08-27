# S11c-a T7 comparator — CROSS-ENGINE RECONCILIATION plan (post-compact roadmap)

## ⭐⭐⭐ DONE — COMPARATOR BUILT + RUN, ENGINES AGREE (2026-08-27, `50f43123`) — this plan is CLOSED; see below
The reviewed comparator (`scripts/S11c_a_cross_engine_comparator.py`, 2 build legs PASS) diffed all 39
families. THE TWO ENGINES AGREE: after rule-13 classification every residual is REPRESENTATIONAL (WL writes
face-eval bulk quantities μ_θ^α/δj as APPLIED FUNCTIONS of x, PY as BARE SYMBOLS; δρ=ρ_br·θ; CONORMAL §3c
verdict-A), NO genuine physics disagreement. ⭐ Built via the rule-15 pivot — the prose directive failed THREE
leg rounds, so rev-4 was a build brief that DELEGATED per-family extraction to the builder with mandatory
accounting. PY transcript committed `afdc8158`. Detail: SCOUT `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§23–24;
STATUS.md front; memory `project_s11c_a_state`. OPEN = (1) from-spec ADJUDICATE the representational identities
are benign (are WL's face-eval args inert per §3a/§1b? is δρ=ρ_br·θ?) — the #3-pattern name-map trap; (2)
characterize the control battery (bite/hold); (3) step record + family card. The blocks below are the
superseded planning history.


## ⭐ PROGRESS (updated 2026-08-26 — Step 5 coverage + CONORMAL/PROJECTION adjudications committed aef0ed36)
- **Step 0 (feasibility/adjudication matrix): DONE, committed `3c7f9137`** —
  `S11c_a_T7_adjudication_matrix.md` (+ twin + census `_measurements/s11ca_t7_census/`), 2 legs + fold.
- **Step 1 (adjudicate vs spec): DONE, committed `3491a376`** — `S11c_a_T7_adjudication_verdicts.md`
  (+ twin + `_measurements/s11ca_t7_adjudication/`), 2 legs + fold. **User chose FULL RECONCILE.**
  Verdicts (fixes on BOTH engines): B density → PY correct; C virtual-work → WL correct; H.1 coverage →
  PY correct; BG → WL correct. See the verdicts doc for the per-engine fix list.
- **Step 2 (ENGINE PATCHES): DONE — BOTH engines committed.**
  - **WL patch committed `a7459cb8`** (`S11c_a_wl_patch_directive.md` + `S11c_a_wl_patch_repair_directive.md`
    + twins + `_legs/S11c_a_wl_patch_review_prompt.md`): B projection density axis + H.1 §5b/§5c coverage +
    kinematic/flux axis drop; needed a repair for a patch-introduced §5a `materialShape` blocker; 2 fresh legs
    (Agent + Grok) PASS; Grok's F1 (PROJECTION_STATIC form control) resolved by computed counts → working
    control (fires 6/24, honest invariant elsewhere).
  - **PY patch committed `2d0f0055`** (`S11c_a_py_patch_directive.md` + twin + `_legs/S11c_a_py_patch_review_prompt.md`):
    C full physical×virtual virtual-work grid (16 cases) + BG boundary loads + density-map branch drop; 2 fresh
    legs (Agent + Grok) PASS; virtual-work off-diagonal crux confirmed genuine 3 ways.
  - ⚠ Shared non-defect note carried forward: BACKGROUND_STATE DIMENSIONS tuple is declarative (skips the
    supplied zeros V⁰/J⁰/𝒜⁰) — both engines may differ in that bookkeeping ⇒ a DIMENSIONS comparator-bridge.
- **Step 3 (SCOUT + MEASURE-FIRST): DONE 2026-08-26 — measurement only, NO commits.** Instrument
  `~/.s11_build/S11c_a_measure_reconcile.py`; findings `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md`.
  Namespace joins 39↔39. ⛔ WL `.out` CONTROL_FORM_* are multi-line (Write-newline) — comparator loader
  reassembles blocks (user's call). ⛔⛔ **No translation dictionary** (S10/S11b = mechanical transliteration +
  join-on-name only; the trap is smuggling an algebraic identity into naming — S10 `omega2→omega**2`).
  Measure-first per the user: a SMALL EXPLICIT map (~11 params + ~18 fields + WL-Derivative→PY-jet decode +
  mixed-partial canon + the ONE spec relation `e_W≡δW/W₀` spec:45). **4 primaries AGREE EXACTLY** (T-0
  FACE_NORMAL, T-b FACE_MEASURE, T-c FACE_VELOCITY, T-g VIRTUAL_CONSTRAINT — 32 cases, residual 0). **ONE
  candidate FINDING (adjudicate vs §3c/§379):** traced-bulk shape-derivatives (RELATIVE_FLUX, TRACTION,
  EVOLUTION, KINEMATIC OPERAND_A) — PY `(∂_n background)·shift` vs WL `∂_n(perturbation)`. ⇒ **NOW RESOLVED
  — see Step 4.**
- **Step 4 (ADJUDICATE + FIX the shifted-trace finding): DONE, committed `c36beac4` (2026-08-26).** The
  dual-engine method caught TWO real PY bugs; WL was correct. Adjudicated by 2 from-spec CAS derivations (fresh
  Agent + Grok, each residual 0) + a source check: **(D1)** PY fabricated background-normal-jet PREMISE symbols
  the spec never supplies (`d_w_v_bulk_0`/`d_w_delta_p_0`/`d_w_j_0`/`d_w_rho_4D_0`) → a spurious §3c term-2;
  **(D2)** PY froze the traced perturbation at the flat face → dropped the mandated ε·η term-1 correction. Fixed
  SPEC-QUESTION-FIRST (§3c clarified in the shared spec, structural premises only; then the PY engine fix):
  directive `S11c_a_py_shifted_trace_fix_directive.md` (+ twin), 2 directive-review legs (Codex+Grok, folded
  once — they caught a leaky, under-scoped rev-1) → build → 2 ablation legs (fresh Agent+Grok, both PASS, form
  ablations confirm term-2 derives from the real background). CONFIRM (fixed PY vs WL via
  `~/.s11_build/S11c_a_reconcile_fixed.py`): **RELATIVE_FLUX 8/8→0 + TRACTION 16/16→0** (TRACTION also needs a
  mechanical `mu_theta_L→mu_theta` rename + `sp.cancel` for the λ_X complex denominators). ⚠ EVOLUTION residual
  = a COMPARATOR single-face-canonicalization gap (its case key has no FACE axis; raw WL complete, both faces),
  NOT a finding. Engines: WL `a7459cb8` UNCHANGED (was correct), PY corrected in `c36beac4`.
- **Step 5 (comparator + coverage): IN PROGRESS (2026-08-26).** Measure-first coverage via
  `~/.s11_build/S11c_a_cov_all.py`: **11 tag-families reconcile EXACTLY to 0** (FACE_NORMAL/MEASURE/VELOCITY,
  VIRTUAL_CONSTRAINT, RELATIVE_FLUX, TRACTION, CLOSURE, EVOLUTION [implemented+verified face-detect],
  KINEMATIC A+B, VIRTUAL_WORK 16). Declared map += 4 params + the EVOLUTION face-detect (`{±1/2*W0}` →
  `XFACEX` tokens). **CONORMAL (T-a′) ADJUDICATED VERDICT A — no finding** (2 from-spec legs + WL-source
  check). ⛔ **PROJECTION "window drops a face" was a FALSE finding — rule-13 catch** (MY truncated-fragment
  prompt; the raw WL window is 2-arg `𝒪(w−W0/2,−w−W0/2)`, both-face derivatives, identical to PY).
  ⭐⭐ **BACKGROUND-CURRENT — was a REAL finding, NOW FIXED (WL `6fae82b8`).** WL carried the rest-frame
  background bulk current as FREE PREMISE symbols (`currentWBackground`/`currentXBackground{i}`), violating §3c
  ("none may be introduced as a free premise"), §1b, and its own `bulkVelocityZero→0` (survived 1660 hits in the
  dynamic projection; PY had `j⁰=ρ⁰·v⁰=0`). 2 from-spec adjudication legs confirmed divergence; on the user's
  rule-6 pivot (symbolic-carry is MORE informative than PY's hardcoded `v⁰=0`), 2 physics consults (Codex+Grok,
  own CAS+stdout) COMPUTED the survivor = 0 under static continuity `∇₄·j⁰=0` (`∫δΩ·∇₄·j⁰` → total `w`-deriv →
  boundary term → 0) — benign, but WL internally inconsistent. USER DECISION: build `j⁰=ρ⁰·v⁰` in WL + record the
  continuity-cancellation in the step record (not a silent hardcode). Directive rev 2 `1610b9e9` (2 legs; Codex
  caught 3 real defects — a leaked downstream result, an acceptance that couldn't distinguish construction from
  `:=0` → the velocity-probe ablation, weak §6 clauses; Grok clean; folded once). Fix `6fae82b8` (2 build legs
  Agent+Grok PASS — velocity-probe ablation shows the bg current TRACKS a nonzero probe ⇒ genuine `ρ·v`, not
  `:=0`; a `Part::pkspec1` symbolic-index bug that ballooned the run to 14 GB was fixed via the `_Integer`
  pattern restriction). Regenerated `.out`: 40 tags, **0** bg-current (was 667), only 10 bg-current-consumer tags
  changed, 30 byte-identical. Detail `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§10–16d.
  **NET: THREE confirmed T7 findings, ALL FIXED — shifted-trace (PY `c36beac4`) + free-premise bg-current
  (WL `6fae82b8`) + current-freezing (PY `49b5c525`).** ⛔ CORRECTION (this was WRONG above):
  the PROJECTION integrand was NOT mechanical. Issue (1) the "perturbation current map" was the exact NAME-MAP
  TRAP the user flagged — a freeze-diagnostic + a CAS consult (Codex+Grok, Q1 = DIFFERENT with concrete
  witnesses) proved `currentWPerturbation(w)` (WL field) and `delta_j_bulk_4` (PY constant) are DIFFERENT
  objects: PY froze the current at the face, so `WINDOW_NORMAL_CURRENT ≡ 0` and the `∂_wδj_w` term §1b requires
  was ABSENT = finding #3 (unanimous from-spec adjudication; fixed `49b5c525`, 2 build legs ablation-PASS).
  Issue (2) the window shape-derivative IBP form (Q2) is genuinely benign/SAME. **NEXT:** the DEFERRED comparator
  (PLAN below — REBUILD rev-1 grounded in reality) → the full cross-engine T7 result (should now show PROJECTION
  AGREEING post-fix) + the FACE_SHIFT/origins/bookkeeping sweep → step record + family card. Detail §§17–22.

- **Step-5 PLAN (Codex comparator build — ⛔ DEFERRED; user chose "lean measurement first, then batch fixes"):**
  ⚠ A rev-1 build directive was written and drew ~14 defects from its own 2-leg gate (uncommitted WIP
  `directives/S11c_a_comparator_build_directive.md` +twin +review leg): WRONG PY input (`exports.py` is a
  LEDGER; the real PY tag stream is the uncommitted FIXED transcript `~/.s11_build/S11c_a_py_fixed_run2.out` [47 tags, post-`49b5c525`]),
  false "S11b reassembles multi-line WL" (that lives in the ephemeral scratch loader), dropped FACE_SHIFT,
  `_RESIDUAL_OPERAND` (engines emit `_RESIDUAL`), **2 SMUGGLING folds** — `mu_theta_L/M→mu_theta` (branch
  collapse) and the CONORMAL §3c Taylor identification — which must move to an explicit REGISTRY, NOT the
  naming map, `canon_key` CASE-LOSS (needs FULL-AXIS keying + extracted/unpaired/duplicate counts),
  DIMENSIONS/BACKGROUND_DENSITY_MAP mislabeled supplied (they are COMPUTED). The controls were instead measured
  by hand and are CLEAN (§§18,20). The REBUILT comparator must: commit a PY `.out` + port the ephemeral loaders,
  then fold the FULL
  measured declared map into the frozen contract (`scripts/S11b_cross_engine_comparator.py` — reuse its
  parsers/transliterate/residual + the multi-line WL reassembly). The map already includes (measured this
  step): params/fields/jet-decode/`e_W`/`d_w_X↔X_dw`/`mu_theta`/`sp.cancel` + the `XFACEX` EVOLUTION
  face-detect [DONE, verified] + the CONORMAL §3c-flat-face identification [DONE, verdict A]; still to add from
  the remaining sweep: the bulk-current map + the IBP window canonicalization + the origin/control/bookkeeping
  bridges. Codex-written ⇒ fresh Agent + Grok legs, synthetic fixtures, rule 5 → freeze + run → THE full
  cross-engine result → step record + family card (pin the reconciliation map forward for S11c-b…e). ⚠
  FACE_SHIFT carries the density rep (16 cases). ⚠ the reconciliation is a spec-grounded per-field map (not a
  byte relabel); the residual after it is the physics. (CLOSURE is ALGEBRAIC traced-bulk, already 16/16→0 —
  NOT a window integral; only PROJECTION_* are Inactive[Integrate].)
- ⚠ The "Revised workflow" below is the original roadmap; steps 0-2 are now complete. The case-structure
  divergences it anticipated were adjudicated in step 1 and fixed on both engines in step 2. Tag granularity
  differs (PY 47 vs WL 40) — the comparator joins per matching key; bookkeeping/control tags bridged or excluded.


## The headline (Codex plan-review, orchestrator-verified): the engines DISAGREE on CASE STRUCTURE
The re-emit was chosen to get a full mechanical cross-check. Codex's plan review then found — and I verified
against the two real `.out`s — that **the two engines diverge deeper than serialization: they compute
different CASE SETS and store different semantic content.** So the "emit-only, physics-preserved" assumption
is FALSE. These divergences are **potential physics findings** (rule 1 — the disagreement IS the measurement),
NOT comparator plumbing. They must be **adjudicated against the spec BEFORE any comparator**; fixing the wrong
engine is a COMPUTATION change with a FULL review, not a relabel.

### Case-structure / semantic disagreements to ADJUDICATE (verified)
1. **VIRTUAL_WORK_SHAPE_DERIV — 8 (PY) vs 16 (WL).** PY ties the virtual DOF to the physical DOF → 8
   `(branch,dof,density)` diagonal cases (sympy 675-687, 919-924). WL loops physical DOF × virtual DOF
   independently → 16 `(branch,density,DOF,VIRTUAL_DOF)` cases (WL 1037-1051; `.out` shows the 8 off-diagonal
   `DOF_x|VIRTUAL_DOF_y`, x≠y, that DON'T EXIST in PY). Spec §4 T-d/§5: is the virtual work the full
   physical×virtual matrix (WL) or the diagonal contraction (PY)? Infects the rep/independence/form/uniform
   controls that consume it.
2. **Projection cluster density-decomposition — PY yes, WL no.** PY loops both density representatives
   (`rho_shape`/`rho0` change; sympy 1002-1039, 1077-1087) → keys carry RHO4/RHOBR; WL projection source has
   NO density parameter (WL 758-786), keys are branch|DOF only. ⚠ Codex ran it: the two PY density cases are
   NOT equal → not a droppable duplicate. Spec: is the projection density-decomposed?
3. **CONTROL_FORM_* coverage** — PY ablates all 5 projection objects + EVOLUTION_TERM_ORIGINS (sympy
   1406-1424); WL ablates only PROJECTION_SHAPE_DERIV + evolution balance (WL 1457-1495). Missing cases.
4. **UNIFORM_LIMIT_* coverage** — PY iterates every RAW_PRIMARY (sympy 1636-1644); WL emits only projection
   shape, omits static/dynamic/residual/origins + evolution origins (WL 1613-1640).
5. **Semantic leaf: graded-term vs coefficient.** PY stores the FULL graded term with explicit `epsilon_shape`
   (e.g. FACE_NORMAL/conormal/measure store a `(background, ε·derivative, total)` tuple, sympy 774-783); WL's
   `shapeDerivative` strips `waveOrder` and records the order only in MULTIGRADE (WL 42-53), so WL `EXPRESSION`
   is the derivative COEFFICIENT. Reconcilable ONLY by a declared algebraic reconstruction identity, ⛔ not
   byte-relabeling. FACE_SHIFT (PY aggregate vs WL FIELD-exploded) and origin inventories are similar.
6. **Envelope + dimension rep** — PY flat `VALUE`/`MULTIGRADE`/col-matrix `DIMENSION` (sympy 188-212) vs WL
   nested `SHAPE_DERIVATIVE.EXPRESSION` + `EXACT_SOURCE` + list `DIMENSION` (WL 943-952); DIMENSIONS +
   HOMOGENEITY have incompatible special envelopes (sympy 1647-1709 vs WL 1750-1808). WL flatten is
   mechanical; PY needs per-family leaf extraction; some are genuine mismatches (rows 1-5).

### Serialization-only divergences (reconcilable by emit-patch OR a trivial-comparator bridge) — from the legs
Case-key encoding (PY positional `Tuple(Str,Integer,Str)`, integer 1 = face AND direction; WL pipe-string with
49 extra tokens DIRECTION/FIELD/ORIGIN/VIRTUAL_DOF + FACE_MAP `PLUS`/`MINUS`; PY unmentioned `BOTH_FACES` axis);
field names (VALUE↔EXPRESSION, MULTIGRADE↔MULTIGRADE_EPSILON_ETA_SIGMAW); CAS reps (window 2-arg BOTH sides;
`Inactive[Integrate]`↔`Integral`, build relationals UNEVALUATED — `parse_mathematica("Equal[0,0]")`→native
True; nested `Derivative[o,{a,b}]`; integer Association keys `<|1->,-1->|>` in BACKGROUND_STATE → S11b parser
rejects). Full evidence: `~/.s11_build/S11c_a_comparator_directive_{codex,grok}.log`.

## State
PY `49b5c525` (projection current-freezing fix; earlier §3c fix `c36beac4`), WL `6fae82b8` (all committed+verified; the WL `.out` under `mathematica/out/` regenerated at `6fae82b8` — 40 tags, 0 background-current). PY tag stream from the FIXED engine = `~/.s11_build/S11c_a_py_fixed_run2.out` (47 tags, NOT committed — the comparator must commit a PY `.out`). WL committed `.out` under `mathematica/out/`. ⛔ COMPARATOR DEFERRED — rev-1 directive drew ~14 defects (see Step-5 PLAN); rebuild grounded in reality. Full session record `~/.s11_build/S11c_a_T7_SCOUT_FINDINGS.md` §§17–22.
Superseded comparator-side directive `S11c_a_comparator_build_directive.md` (+twin+review prompt) kept for the
record. ⚠ Running the SymPy engine rewrites committed `S11c_a_exports.py` (Dummy-index counters) AND the
emitted payload feeds `export_candidates` (sympy 333-343 → 1852-1878), so ANY PY emit reformat mutates
exports — decouple/​preserve/​review+commit the regenerated export deliberately. ⚠ another session works
untracked in repo `memory/` — commit EXPLICIT paths only, never `git add -A`.

## Revised workflow (post-compact)
**0. FEASIBILITY / ADJUDICATION MATRIX (the real first step; reviewed).** One systematic pass, every tag:
`(tag, semantic leaf, PY producer path, WL producer path, exact case axes+count, bookkeeper convention,
verdict ∈ {1:1-serialization | identical-duplicate-axis | MISSING-case | SEMANTIC-mismatch})`. This IS the
cross-engine audit. Author it, then 2 legs (it decides physics).
**1. ADJUDICATE each non-serialization row against the spec** (§4 T-a…T-i, §5 controls, §7 grammar): which
engine is correct, or did the spec under-specify? A MISSING-case or SEMANTIC-mismatch is a FINDING → spec
decision → **computation-layer patch to the wrong engine (a one-engine physics fix ⇒ FULL review, and it is a
spec question first — `[[feedback_one_engine_fix_is_a_spec_question]]`)**. A spec under-specification → a
narrow normative §7/§4 addendum (2 legs+fold; lighter than reopening the whole closed doc).
**2. Reconcile the SERIALIZATION-only rows** — either a per-engine emit-patch to a pinned shared schema
(labeled typed-axis keys — `FACE→FACE_PLUS` vs `DIRECTION→DIRECTION_n` disjoint, handle `BOTH_FACES` as its
own scope axis not a face sign; flat `{EXPRESSION,MULTIGRADE_EPSILON_ETA_SIGMAW,DIMENSION_L_T_M}` with one
ordered-triple dimension rep; per-family canonical leaf), OR keep them as trivial-comparator bridges. Choose
per row in the matrix.
**3. PROVENANCE-MANIFEST gate** (replaces "byte-compare values") for every emit-patch: `(tag, canonical case,
semantic slot) → (old path, new path, canonical CAS hash)`; require COMPLETE coverage; many-to-one only if
hashes prove identical; drops only via an explicit diagnostic allowlist; for the ε-coefficient normalization
verify a declared algebraic reconstruction identity (not byte equality); validate case sets, leaf counts,
dims, multigrades, duplicate-axis rejection. **Light re-review is valid ONLY for a proven pure adapter under
this gate; any computation change (rows 1-5) needs a FULL engine review.**
**4. Re-run** (WL ~23 min, watch orphan kernel; PY ~3 min) → commit aligned engines + WL `.out`; regenerate
PY `.out`; deliberately handle the export.
**5. TRIVIAL comparator** — join, promote both keyed maps to Associations, residual per matching key, frozen
T7 contract, + the CAS-syntactic bridges. Build (Codex) + 2 legs (Agent+Grok), synthetic fixtures (rule 5).
Freeze + run → cross-engine result.
**6. Step record** (carry: T-c′ definitional-identity residual; §5a vacuous for T-d LAB_HELD; WL dead code
`widthGradientAnsatz` 797-800; AND the case-structure findings + how each was adjudicated) → family card.
⭐ Pin the shared §7 schema + the adjudicated case-structure FORWARD for S11c-b…e.

## ⚠ Scope note for the user
This is bigger than the "emit-only re-emit" I quoted: several families need SPEC ADJUDICATION and a possible
COMPUTATION-layer patch to one engine (with full review), because the engines genuinely disagree on case
structure. That is the cross-engine method paying off — it caught real divergences — but the first work is
adjudication, not serialization. Cost is therefore higher and front-loaded on physics decisions.

## Governing reminders
Codex `-c model_reasoning_effort=xhigh`; `--sandbox danger-full-access` for Mathematica; ONE grok/user;
background via `run_in_background`, ⛔ never `&`; logs OUTSIDE repo, prompts UNDER project; rule-2 twin +
leak-gate every directive; commit reviewed baseline before overwriting; no commit before both legs report.
