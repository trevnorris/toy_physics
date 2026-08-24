# S11c-a T7 comparator — CROSS-ENGINE RECONCILIATION plan (post-compact roadmap)

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
SymPy `9b6438fa`, WL `ddecdbc2` (both committed+verified). PY tag stream = fresh run to
`~/.s11_build/S11c_a_sympy_engine.out` (NOT committed). WL committed `.out` under `mathematica/out/`.
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
