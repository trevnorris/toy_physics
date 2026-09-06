Reading additional input from stdin...
OpenAI Codex v0.153.4
--------
workdir: /var/projects/toy_physics
model: gpt-5.6-sol
provider: openai
approval: never
sandbox: danger-full-access
reasoning effort: xhigh
reasoning summaries: none
session id: 01a0756b-4874-7ab0-adfc-9de6cbe018d0
--------
user
# Verify — the S11c-c2 compact-prep state (docs accurate? next-prompt correct? clear to compact?)

The orchestrator is about to compact context. Verify, against the REAL repo + git history, that (1) the state
docs are accurate and not overstated, (2) the carry-forward set is right, and (3) the hand-back prompt is correct
and complete. ⛔ Do not rubber-stamp; if something is wrong or overstated, say so plainly with file:line / commit
evidence. This is a document + light-git verification; ⛔ do not modify the working tree.

## What just happened (claimed)
Two things closed this session, both committed on `ledger-v3-rebuild`:
- **c2 fold physics** was independently 2-leg reviewed (fresh Claude agent + Grok) → **SOUND, 0 confirmed defects**,
  committed `8f3a017f`. Grok flagged three "defects" (F uniform-limit, G adjointness, E N6) that the orchestrator
  adjudicated as FALSE POSITIVES by its own computation (rule 13).
- **c2 export repair** (publication-only, astra build) shrank `scripts/S11c_c2_exports.py` **60,516,900 →
  22,441,522 bytes**, directive gated `a5f7a06c` (2 decision legs), committed `aa76105a`, re-reviewed clear.

## 1. STATUS + adjudication accuracy (verify against git + the files)
Read the new top clause of `STATUS.md` (2026-09-06, "c2 STEP A + STEP B/C DONE") and the two adjudication records
`_measurements/S11c_c2_physics_review_adjudication.md` and `_measurements/S11c_c2_export_repair_rereview_adjudication.md`.
Check against `git log --oneline -8`, `git show --stat 8f3a017f aa76105a`, and the artifacts:
- Are commits `8f3a017f` / `a5f7a06c` / `aa76105a` / `d2befb7c` real and as described?
- Is `scripts/S11c_c2_exports.py` actually ~22.4 MB, does it import, is `s11cc2SelfEnergyIncrement` absent, and are
  both `s11cc2ClosedSlabOperator` + `s11cc2ClosedCouplingKernel` present (4-case trees)?
- Is the diff `8f3a017f`→working of `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` genuinely publication-only
  (only `EXPORT_ROOTS`, the `export_key` map in `run()`, `publish`, and the new `publication_compact` helper — NO
  construction function changed)?
- Is "physics SOUND, 0 defects" appropriately supported, or overstated? In particular: are the F/G/E dispositions
  (F genuine coupling decouples / G directional self-energy / E leading-order rep-invariance holds, σ_W deferred)
  defensible from the adjudication's own evidence, or does any read as the orchestrator rationalizing away a real
  leg finding?

## 2. The carry-forward set (right + complete?)
The step record is claimed to owe: (F) "genuine coupling decouples", not "increment vanishes" (+ a light §5e/§3c
spec-wording fix); (E) σ_W-sector rep-invariance remnant deferred; (G) directional self-energy; the two S11c-b
sign conventions surfaced by the WL comparator. Against the physics adjudication + `directives/S11c_c2_SHARED_PHYSICS.md`
(§5e, §3c, §7): is this set correct and complete, or is a real carry-forward missing / miscategorized (e.g. is any
of these actually a build defect rather than a record note)?

## 3. The hand-back prompt (correct + complete + clear to compact?)
Read `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md`.
Verify: the STATE line matches reality; the NEXT sequence is right (STEP 1 light spec-wording fix = review-until-clear;
STEP 2 blind WL engine with 2 decision legs BEFORE + 2 build legs SERIALIZED; STEP 3 T7 comparator + reconcile →
step record); the model policy + ops notes are correct; nothing outstanding is being skipped. Is anything in the
prompt wrong, missing, or an overstatement? Is there any reason we are NOT clear to compact (an uncommitted
load-bearing artifact, an open review leg, a claimed-committed thing that isn't)?

## Output
For each of 1–3: your finding + evidence (commit / file:line). End with: **CLEAR TO COMPACT** or the exact list to
fix first, and (if any) the precise wording corrections for the STATUS clause / next-prompt.
codex
I’ll verify this as a read-only audit: first the repo/history and exact artifacts, then the review logic and directive wording, and finally the external hand-back prompt plus working-tree state.
exec
/bin/bash -lc "git status --short --branch && git log --oneline --decorate -8 && git show --stat --oneline --summary 8f3a017f && git show --stat --oneline --summary aa76105a && git show --stat --oneline --summary a5f7a06c && git show --stat --oneline --summary d2befb7c && stat -c '%n %s bytes' scripts/S11c_c2_exports.py && rg -n '"'^#|c2 STEP A|STEP B/C|8f3a017f|a5f7a06c|aa76105a|d2befb7c|SOUND|carry|sigma|σ|directional|genuine coupling'"' STATUS.md _measurements/S11c_c2_physics_review_adjudication.md _measurements/S11c_c2_export_repair_rereview_adjudication.md | head -n 500" in /var/projects/toy_physics
 exited 1 in 0ms:
## ledger-v3-rebuild
?? research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md
d2befb7c (HEAD -> ledger-v3-rebuild) STATUS: c2 STEP A (physics SOUND) + STEP B/C (export repair 60->21 MB) DONE + committed
aa76105a S11c-c2 export repair COMMITTED — 60 MB -> 21.4 MB (publication-only, re-reviewed clear)
a5f7a06c S11c-c2 export-repair directive — GATED (2 decision legs, one pass, folded)
8f3a017f S11c-c2 self-energy fold — physics review (2 legs) ADJUDICATED: fold physics SOUND (0 defects)
42428835 S11c-c2 compact-prep: fold Codex verify — export = BOTH closed operators (not coupling-only); transparent-factored + equiv-check; STATUS hedges
eef7c269 STATUS: STEP 1 c1 record corrections DONE; c2 build directive gated + astra build ran (overstepped + over-exported) → repair pending
1ae6c336 S11c-c2 build directive — fix 2 defects the first astra build exposed (pre-repair)
c1de32b0 S11c-c2 SymPy build directive — GATED (2 decision legs) + folded one pass
8f3a017f S11c-c2 self-energy fold — physics review (2 legs) ADJUDICATED: fold physics SOUND (0 defects)
 .../S11c_c2_adjudication_verify_EG.out             |   33 +
 .../S11c_c2_adjudication_verify_EG.py              |   51 +
 .../S11c_c2_adjudication_verify_F.out              |   36 +
 .../_measurements/S11c_c2_adjudication_verify_F.py |   60 +
 .../_measurements/S11c_c2_physics_claude_leg.md    |   74 +
 .../_measurements/S11c_c2_physics_grok_leg.txt     |  169 ++
 .../S11c_c2_physics_review_adjudication.md         |  114 +
 .../_measurements/S11c_c2_sympy_builder_report.md  |  370 ++++
 .../S11c_c2_sympy_object_provenance.json           | 2339 ++++++++++++++++++++
 .../_measurements/S11c_c2_sympy_stdout_index.json  | 1395 ++++++++++++
 .../_legs/S11c_c2_physics_review_prompt.md         |  149 ++
 .../scripts/S11c_c2_selfenergy_fold_sympy_audit.py | 1053 +++++++++
 12 files changed, 5843 insertions(+)
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_physics_claude_leg.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_physics_grok_leg.txt
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json
 create mode 100644 research/pde_ledger_v3/directives/_legs/S11c_c2_physics_review_prompt.md
 create mode 100644 research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
aa76105a S11c-c2 export repair COMMITTED — 60 MB -> 21.4 MB (publication-only, re-reviewed clear)
 .../S11c_c2_export_repair_adjudication_verify.out  |    7 +
 .../S11c_c2_export_repair_adjudication_verify.py   |   83 +
 .../S11c_c2_export_repair_builder_report.md        |  523 +++++
 .../S11c_c2_export_repair_rereview_adjudication.md |   54 +
 ...t_repair_rereview_claude_recompute_artifact.out |   15 +
 .../S11c_c2_export_repair_rereview_grok.txt        |  173 ++
 .../S11c_c2_sympy_guard_evidence.json              | 2028 ++++++++++++++++++++
 .../_measurements/S11c_c2_sympy_progress.json      |   10 +
 .../_legs/S11c_c2_export_repair_rereview_prompt.md |   72 +
 research/pde_ledger_v3/scripts/S11c_c2_exports.py  |   83 +
 .../scripts/S11c_c2_selfenergy_fold_sympy_audit.py |  168 +-
 11 files changed, 3210 insertions(+), 6 deletions(-)
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.out
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.py
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_grok.txt
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_sympy_progress.json
 create mode 100644 research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_rereview_prompt.md
 create mode 100644 research/pde_ledger_v3/scripts/S11c_c2_exports.py
a5f7a06c S11c-c2 export-repair directive — GATED (2 decision legs, one pass, folded)
 .../S11c_c2_export_repair_decision_codex.md        | 120 ++++++++++++++++
 .../S11c_c2_export_repair_decision_gate.md         |  36 +++++
 .../S11c_c2_export_repair_decision_grok.txt        | 153 +++++++++++++++++++++
 .../directives/S11c_c2_export_repair_directive.md  | 119 ++++++++++++++++
 .../_legs/S11c_c2_export_repair_decision_prompt.md |  58 ++++++++
 5 files changed, 486 insertions(+)
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_codex.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_gate.md
 create mode 100644 research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_grok.txt
 create mode 100644 research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md
 create mode 100644 research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_decision_prompt.md
d2befb7c STATUS: c2 STEP A (physics SOUND) + STEP B/C (export repair 60->21 MB) DONE + committed
 STATUS.md | 8 +++++++-
 1 file changed, 7 insertions(+), 1 deletion(-)
stat: cannot statx 'scripts/S11c_c2_exports.py': No such file or directory

codex
The named commits exist in the expected order. One immediate caveat: the tree is not literally clean—there is one untracked verification-prompt file. I’ll determine whether it is merely this audit’s input or an outstanding load-bearing artifact before treating that as a compaction blocker.
exec
/bin/bash -lc "stat -c '%n %s bytes' research/pde_ledger_v3/scripts/S11c_c2_exports.py && git ls-files --stage research/pde_ledger_v3/scripts/S11c_c2_exports.py research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && nl -ba STATUS.md | sed -n '1,90p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md && nl -ba /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md" in /var/projects/toy_physics
 succeeded in 0ms:
research/pde_ledger_v3/scripts/S11c_c2_exports.py 22441522 bytes
100644 210812ec4b63b1a44d83cf60bd55f9a4c74ec6b7 0	research/pde_ledger_v3/scripts/S11c_c2_exports.py
100644 43218c15643605271de5ce4cefad7081af4a11ce 0	research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
     1	# STATUS — where the Path-A program is (single front door)
     2	
     3	> ⚠ **SHA NOTE (2026-09-01 DataLad rewrite):** the git history of `ledger-v3-rebuild` was rewritten — **every commit SHA cited in a clause dated before 2026-09-01 is PRE-rewrite and no longer resolves.** Find those commits by their message, not the SHA. Post-rewrite anchors: #89 directive `660a3055`, #88 `05cb1ea5`, #87 `bab2b828`, #89 engine checkpoint `f655ea65`, #89 clearance `9f40c18e`, **#89a WL basis `d4adbd99`**. `.out` files are now git-annex pointers (`datalad get` to read).
     4	
     5	## ⭐⭐⭐ c2 STEP A (physics review) + STEP B/C (export repair) DONE — fold physics SOUND; export 60 MB → 21.4 MB; both committed; NEXT = light spec fix → blind Wolfram engine (2026-09-06)
     6	**(A) c2 PHYSICS 2-leg review DONE — fold physics SOUND, 0 confirmed defects, committed `8f3a017f`.** Real independent legs (Codex-written script → **fresh Claude agent + Grok**, identical prompt, blindness-by-absence). Both agree the fold WIRING + A/C/D1–D6 are SOUND (δp not J_s, operator-inverse, kernel bridge, V_s→face_velocity, computed jets, ε-strip; FORM ablations respond). **Grok flagged F/G/E as CONFIRMED defects — ALL FALSE POSITIVES**, refuted by my own computation (rule 13; `_measurements/S11c_c2_physics_review_adjudication.md` + `…_adjudication_verify_{F,EG}.py/.out`): **F** uniform-limit — the genuine closure-induced coupling **DECOUPLES** (Integral integrand ≡0, `.doit()`→0 all 4 cases; Grok saw `Integral(...)` without evaluating it — its recurring pattern, cf. c1 "tangential freeze"); **G** adjointness — the builder never claimed adjoint-by-construction, both blocks emitted, the self-energy is **DIRECTIONAL/one-way** (reverse block ≡0, expected since δp lives only in the θ/mech rows) — dissipative directional, not a defect; **E** N6 — leading O(ε)/O(εη) rep-invariance **HOLDS** (σ_W→0 annihilates the whole residual), σ_W-sector remnant DEFERRED. **Carry-forwards (step-record, ⛔ NOT defects):** F wording "genuine coupling decouples" NOT "increment vanishes" (→ light §5e/§3c spec fix owed); E σ_W-sector deferred; G directional self-energy; the 2 S11c-b sign conventions → the WL comparator.
     7	**(B/C) EXPORT REPAIR DONE + committed `aa76105a` — 60,516,900 → 22,441,522 bytes (63% smaller), publication-only, physics untouched.** Directive gated `a5f7a06c` (2 decision legs Codex sol xhigh + Grok, one pass, the R2 guard-spec blocker folded — `_measurements/S11c_c2_export_repair_decision_{gate,codex.md,grok.txt}`). astra build (`gpt-6-astra` high, setsid-detached) **BEHAVED — no overstep** (the fixed §7 fence held) but SLOW (~50 min / 3 full runs). Repair: drop `SelfEnergyIncrement` from `EXPORT_ROOTS`+the export-key map → EMIT-only; new `publication_compact` (transparent-factors **only** the per-case VALUE leaves, keeps the cased tree, preserves Integrals + the pole set — **reverts** if factoring would change it — and raises if any temp Dummy escapes); a **separate** strict-recursive Integral-aware semantic guard (srepr roundtrip kept alongside); `display`→root-name. **2 re-review legs** (astra-written → fresh Claude + Grok): **Grok CLEAR** (guards bite — `out+1` compaction + membership re-add/drop all hard-stop). The **fresh-Claude leg WEDGED** on reaped background waiters and produced **2 FALSE POSITIVES** — both artifacts of its **cross-process** re-run (fresh Integral bound-dummies vs the stored compact: `k·d1·X − k·d2·X`). **My rule-13 verify** (same-process, `_measurements/S11c_c2_export_repair_adjudication_verify.py/.out`): FAITHFUL **5/5**, GUARD BITES **5/5** (doubling, marker-add, in-integrand) → repair SOUND; adjudication `_measurements/S11c_c2_export_repair_rereview_adjudication.md`.
     8	**Lessons this session:** (1) Integral-bearing operators MUST be compared **same-process** (cross-process re-run mismatches bound dummies → spurious FAIL). (2) **Every `run_in_background` task was reaped** (~30 s, even trivial loops; NOT OOM — 23 GB free); the astra build survived only via **setsid-detached launch + DONE-marker + re-armed watcher** [[feedback_background_tasks_can_die_spuriously]]. (3) **astra assessment** (user asked): KEEP astra for CODE — this build's code was correct + defect-free under adversarial review; the earlier c1 overstep was a directive-fence defect (fixed), not astra; only weakness is speed.
     9	**NEXT:** (a) light **§5e/§3c spec-wording clarification** ("genuine coupling decouples", not "increment vanishes"; review-until-clear); (b) the **blind Wolfram engine** (imports nothing, re-derives) → its **2 legs** → the **c2 T7 comparator + reconcile** → the **c2 step record** (NO per-substep card — one S11c roll-up after S11c-e). Live handoff: [[project_s11c_c_state]]. ⛔ c1 STANDS. ⚠ `MEMORY.md` ~20 KB — compact next fresh session.
    10	
    11	## ⭐⭐⭐ STEP 1 (c1 record corrections) DONE; c2 BUILD DIRECTIVE GATED; astra build RAN — OVERSTEPPED + OVER-EXPORTED → repair; NEXT = real legs + export repair (2026-09-05) [SUPERSEDED by the clause above]
    12	**(1) STEP 1 — c1 RECORD corrections DONE + committed `efae0337`.** Review-until-clear (2 rounds, Codex sol xhigh + Grok, both computation-backed). Round 1 found 2 defects (grazing-inverse **mis-attribution** — the `~1/η` singular object is the DtN inverse `N⁻¹`/`Z`, ⛔ NOT the permeable resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹` which is REGULAR at grazing for Λ_A≠0; energy-parenthetical sign), folded with the legs' wording; round 2 both **CLEAR**. c1 STANDS confirmed **4×** (density opaque-recoverable, residual 0). Gate record `_measurements/S11c_c1_record_corrections_review_adjudication.md`.
    13	**(2) STEP 2 — c2 SymPy BUILD DIRECTIVE (thin; physics=pointers to spec `16849fc6`) GATED + committed `c1de32b0`.** 2 decision legs (both ran `ledger_fold`): wiring SOUND (18 closure-covering roots + `coupling_kernel` 19th; exact-lookup `IMPORT_KEYS` rule; bare-vs-prefixed `face_response` hazard; `v_bulk_normal_0` excluded), **7 symbol-map defects folded one-pass** (source=`DELTA_P` not `PRESSURE`; w-jets COMPUTED not imported; **ε-normalization** avoid O(ε²); `V_s`→`face_velocity` not `kinematic_balance`; `dtn_operator`→`dtn_kernel` bridge; `coupling_kernel` promoted; clauses verbatim). Gate record `_measurements/S11c_c2_build_directive_decision_review.md`.
    14	**(3) astra SymPy build RAN** (`gpt-6-astra` high, detached, 19 min, EXIT 0) — the **builder REPORTS** implementing the intended construction (`_measurements/S11c_c2_sympy_builder_report.md`; ⚠ physics **NOT independently verified** — that report even carries `DERIVED_OR_DECLARED: FAIL` + a parser-limitation, and — now stale — claims no reviewer was launched) — **BUT two problems, both from MY directive defects (now fixed `1ae6c336`, refined post-Codex-verify):** (a) **OVERSTEP** — astra read the `/build` skill and ran its own INVALID self-reviews + ~40 orchestrator files (self-check trap; §7 had embedded the orchestrator's review-launch); **quarantined (37 files)** to `/var/projects/toy_physics_ext_logs/c2_builder_overstep/`; [[feedback_builder_directive_no_orchestrator_process]]. (b) **OVER-EXPORT** — 60 MB = 3 fully-expanded closed objects (`ClosedCouplingKernel` 24.7 / `ClosedSlabOperator` 19.9 / `SelfEnergyIncrement` 15.9); §5 didn't split EMIT (→`.out`) from EXPORT (→delta). **Directive fixed:** §7 bounds the builder to build→verify→report; §5 exports **BOTH closed operators** (`CLOSED_SLAB_OPERATOR` full closed op — d needs the closure-modified diagonal for its spectrum — + `CLOSED_COUPLING_KERNEL` off-diagonal, per d's declared scope `S11c_decisions.md:47,52,83`), in a **transparent factored** encoding with a semantic-equivalence check, and the increment (a comparison/emit artifact) + term_origins + §3d + §5 controls EMIT-only (Codex-verify corrected an earlier "coupling-kernel-only" guess). **Artifacts SALVAGED** (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py`, `scripts/S11c_c2_exports.py`, the builder report — **UNcommitted**, pending review + export repair).
    15	**NEXT:** (a) run the **REAL independent legs** (Codex-written → **fresh Claude agent + Grok**, orchestrator-launched) on the **physics** (fold map, ε, kernel bridge, N6 rep-invariance, the uniform-limit `delta_p_minus` concern astra self-flagged, the 6 re-adjudications, mandatory FORM ablation + one-sided corruption); (b) **export repair** (membership + held/factored — a script re-run); (c) adjudicate → re-review until clear → commit; then the **blind Wolfram engine** → **T7 comparator + reconcile** → **c2 step record**. Live handoff: [[project_s11c_c_state]]. ⚠ `MEMORY.md` ~20 KB — compact next fresh session. **The clause below is SUPERSEDED — kept for the artifact map.**
    16	
    17	## ⭐⭐⭐ c2 SPEC GATED + c1 RETRO-REVIEW DONE (c1 STANDS) + CLAUDE.md RESTRUCTURED (2026-09-05) [SUPERSEDED by the clause above]
    18	Three things landed after the c1 close. **(1) c2 SHARED_PHYSICS spec (the self-energy fold) GATED + committed `16849fc6`** — 2 decision legs (Codex sol xhigh + Grok) both ⚙-computational + CONVERGED, v1 NOT sound (central **isolation claim FALSE** — the substitution increment does NOT cancel S11c-b's two unvalidated SIGN CONVENTIONS; they multiply the substituted δp_s slots, both legs computed `2·carrier·incr`), folded ONE pass → **v2** (9 findings: real θ-row has `J_s` ELIMINATED ⇒ substitute **closed δp_s+w-jets** NOT J_s; `Λ_X` traction-only; **assembled two-face** operator; **3→6 re-adjudications**; positional `load_model` + `s11c_c1_` keys; `Z→0≠Λ_A=0`). Record `_measurements/S11c_c2_spec_review.md`. **(2) c1 SPEC RETRO-REVIEW DONE** (user's *'full correct path'* — the c1 spec was folded once, never re-legged): 2 legs DISAGREED (Grok CLEAR / Codex BLOCK-reopen), adjudicated (rule 13, `_measurements/S11c_c1_spec_retro_review_adjudication.md`) — ⭐ **c1's ENGINES/EXPORTS STAND, ⛔ do NOT reopen** (PY carries density OPAQUELY, 0 derivatives; re-binding to `background_density_map` recovers the O(εη) channel exactly → c2 is **specified** to recover it, c2 not yet built; Codex-sol verify CONFIRMED the override), but record corrections OWED — **3 MUST-level** (seal-5 "harmless O(η²)"→**O(εη)-recoverable-representational**; grazing→`NOT_ESTABLISHED`; independence-scoping) **+ lower-severity carry-forward** (energy-residual orientation, `h_s`/DtN-vs-impedance terminology, density-as-mult-operator, `K_a` Hermitian, η²/ησ_W/σ_W² caveat, drain-projection, flat-`Z₀` leakage; full list in the adjudication record). **(3) CLAUDE.md RESTRUCTURED + committed `3da095c6`** (17 rules→M/E/G/S + an **artifact→review-discipline table**; the fix: **specs=review-until-clear**, not fold-once; ⚠ the preservation legs caught + removed a self-introduced **FABRICATION** in the rewrite). **NEXT = (a) the 3 c1 RECORD corrections (2 legs, review-until-clear; ⛔ NOT engine reopens), THEN (b) c2 SymPy BUILD DIRECTIVE** (own 2 decision legs BEFORE the astra build; determine `IMPORT_KEYS` root set + fold symbol map vs the real files). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = c2) is SUPERSEDED — kept for the artifact map.**
    19	
    20	## ⭐⭐⭐ S11c-c1 DONE — T7 comparator + CROSS-ENGINE RECONCILE + step record committed `072d0b75`; NEXT = c2 (2026-09-05) [SUPERSEDED by the clause above]
    21	Both c1 engines per-engine SOUND (clause below). The **T7 comparator** was gated (`84686a54`, 2 decision legs, 12 folds) → built → reviewed baseline `7141e6ad` → scoped repair (astra) `704308af`, all re-review legs CLEAN = **SOUND** (measurement-only, prints/decides nothing). ✅ **CROSS-ENGINE RECONCILE RUN** (family-scoped, 46/50 families, 30 GB, peak ~317 MB) + **c1 STEP RECORD committed `072d0b75`** (`steps/S11c_c1_curved_bulk_closure.md`; record `_measurements/S11c_c1_comparator_reconcile.md`; reproducible `_measurements/S11c_c1_reconcile_reproduce.py`). Method = the **staged representational bridge** ([[feedback_reconcile_representational_bridge]]): comparator PRINTS raw (nonzero≠disagree), I apply identifications read from the engine sources (⛔ not tuned to zero; adversarial corruptions bite) and test collapse. ⭐ **EARNED cross-engine AGREE:** the two-momentum **DtN KERNEL** — all 4 (anchoring×face) cases EXACT ZERO off-diagonal (Stage-2 residual collapsed to the dispersion relation itself) — and the **face-response COEFFICIENTS** (LAB_HELD off-diagonal + MATERIAL on-diagonal); seals 1–4 AGREE. **2 step-record legs (Grok + Codex sol xhigh) SCOPED an earlier "closure-wide AGREE" over-claim** — folded once (`_measurements/S11c_c1_step_record_review_{grok,codex_sol}.txt`). ⛔⛔ **UNDECIDED — ⛔ NOT cross-engine-closed (c2 must NOT inherit these):** (5) **bg-density = a SURFACED rule-17 freeze** (PY bare constant vs WL live field `(ρbr/W₀)WBg(x)`; functional dependence agrees + 0 derivatives, but a scalar≠a field globally — c2 re-adjudication MANDATORY); the raw **`dtn_operator`** whole-form (kernel-AGREE does NOT extend to it); **ENERGY** (closed-form vs far-field integral); an off-diagonal **flat-resolvent leg-labeling**. ⚠ **4 giants + full symbolic residual UNMEASURED — DEFERRED ≥64 GB** (`DEFERRED_HEAVY_RUNS.md`, ⛔ not pre-adjudicated). No cross-engine DISAGREEMENT found. **The kernel is cross-engine closed; the rest is UNDECIDED/deferred.** ⚠ Codex-sol correction-verify pass was run on the fold (in flight/complete). **NEXT = c2** (self-energy fold, a CODE build → gpt-6-astra, own gated directive with 2 decision legs BEFORE the build; carry the mandatory seal-5 + `dtn_operator` re-adjudication). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = comparator) is SUPERSEDED — kept for the artifact map.**
    22	
    23	## ⭐⭐⭐ S11c-c1 BOTH ENGINES DONE + per-engine reviewed; NEXT = the c1 T7 CROSS-ENGINE COMPARATOR (2026-09-04) [SUPERSEDED by the clause above]
    24	The export migration + delta rebuild (clause below) is **DONE** (`f90e7630`); the c1 SymPy engine is **DONE** (`d6e16471`, 906 KB own-rows delta export); the c1 **blind Wolfram** engine is **DONE** (`c2b0f9dd`→`e139bc61`→`13f0bd2c`→ **repair-2 `dd34d564`**, fully gated, both re-review legs CLEAR). ✅ **A FULL REMEDIATION is COMPLETE** (`_measurements/S11c_c1_wl_remediation_plan.md`): the WL repair DIRECTIVE had skipped its rule-7 decision legs; the remediation re-ran repair-2 **WITH** its decision legs, which **caught that a first-draft fix targeted a non-defect** (would have corrupted the correct `PERMEABLE_PORT_HERMITIAN`) — the gate working. Both c1 `.out` are committed (`4a14100a`, annex+GIN; WL SHA-256 `d63e616786…`). ⛔ **Both c1 engines are BUILT + reviewed — do NOT rebuild.** **NEXT = the c1 T7 cross-engine comparator** — frozen contract `directives/S11c_c1_SHARED_PHYSICS.md:580-587` (N8 inherited verbatim; ⛔ point at it, don't paraphrase): it **computes and prints, deciding nothing** (join by object name with the S11c-a axis-typed keys, injective `mechanical_lower_camel`, paired residual operands + per-family accounting, repoint ablation); the three-valued AGREE/DISAGREE/UNDECIDED characterization is the POST-RUN adjudication, ⛔ not a script verdict; ⛔ no representational fold pre-registered. Reuse `scripts/S11c_a_cross_engine_comparator.py`. Orchestrator-written → Codex+Grok decision legs BEFORE the build → build → 2 legs. Known representational residuals to adjudicate AFTER the run (never pre-registered): μ_θ, the ω-assumption artifact, the two-momentum leg names, and the new proportional-4 parity NIT. ⚠ The full cross-engine residual may need the ≥64 GB box (this box is 30 GB) — build per-engine here, defer the heavy residual if it OOMs. Live handoff: [[project_s11c_c_state]]. **The two clauses below (c1 build PAUSED / NEXT = migration) are SUPERSEDED — kept for the artifact map.**
    25	
    26	## ⭐⭐⭐ EXPORT ARCHITECTURE REDESIGN (bind-closure LEDGER) — machinery BUILT + VERIFIED; c1 build PAUSED for it; NEXT = c1 MIGRATION + delta rebuild (2026-09-04) [SUPERSEDED by the clause above]
    27	The user's *"why is the ledger still that big"* question **mid-c1-build** triggered a full export-chain redesign. **Measured:** the LEDGER over-exports **~56%** by the bind test (`S11c_b_exports.py` = 58.75 MiB, 2441 rows → ~26 MiB bind-closure; **~69% / 18 MiB after dropping `coupling_kernel`** per D4) — **~29% is `*_term_origins` diagnostics** (29.5% = all origin/provenance rows incl. `dissipation_origin`) + dead accumulated S9/S10/S11/S11b history nothing binds. Root cause: membership was "every §4 MAIN primary" (`D1`) **and** each write ACCUMULATES the whole inherited file forward; **`F10` (2026-09-03) re-asserted a topology the 2026-08-11 `chain_accumulate_or_generate_decision.md` had already SUPERSEDED (generate, not accumulate) and that was never built.**
    28	- **Process:** two-engine **deliberation** ("what is the minimum export?", Codex+Grok, convergent) → the **bind-closure design** `directives/export_ledger_bind_closure_design.md` → **two-leg design review** (Codex+Grok, both "not sound as-is", convergent must-fixes) → v2 folded → committed **`c04e071f`**, §D3 refined **`fd8c89d0`**. (Earlier: the first trim `d1fe1bf0`; made lawful via F10 `3ee34ebc` — F10 now SUPERSEDED `4274b403`.)
    29	- **The design (committed):** **D1 Membership = the BIND-CLOSURE** — export a row iff a later step binds its **F9 write-key** as a construction/regression operand, recursively closed over symbols + `dimension_key`; **by binding, not label**; missing manifest = **hard error**. **D2 Topology = GENERATE over a frozen base** — F9 is write-time, the fold is chronological **last-wins (never re-applies F9)**, **no write-time pass-through**; the existing files **through `S11c_b_exports.py` are FROZEN** ⇒ no chain re-run, no `BUILD_INPUT_DIGESTS` cascade. **D3 Guard** = manifest + **edge resolution by full SYMBOL IDENTITY (srepr = name+assumptions)** + bidirectional smoke-test; **`D3`-round-trip provably does NOT catch under-export**. **D4** = representation deferred; **`coupling_kernel` DROPPED** from the minimum (c2 re-extracts, doesn't bind it) ⇒ minimum ~319 rows/18 MiB.
    30	- **F9c-severity scan** (`_measurements/f9c_pair_scan.py`): 16 routed pairs, **11 identity-resolvable, 5 genuine same-srepr — and 0 of the 5 on the critical path** (only `omega`, resolvable). ⇒ the "routed-key contract" is **not** a c1/c2 blocker; the 5 off-path pairs are a documented **deferred residual**.
    31	- **Fold+guard MODULE `scripts/ledger_fold.py` (+ `test_ledger_fold.py`): BUILT + 2-leg-reviewed + 2 repairs + VERIFIED** (committed **`1be64084`**). `load_model(base ⊔ deltas)` merges by exact key (no F9 re-apply); `check_consumer` resolves edges by identity (`omega` resolves, structural `O_window`/coordinates skipped, symbol *users* aren't mistaken for *producers*), smoke-test = the under-export backstop. **Verified: `check_consumer` resolves all 12 tested roots** (the 11 critical-path roots + `coupling_kernel`, resolved cleanly though c2 doesn't bind it) against the 2441-row frozen base (the mandatory `real_ledger_roots` test — its absence let repair-1's false-positives ship). Directives: `8041799f` (build), `57c4e57a`/`4cc4073d` (repairs), baseline `05f86947`.
    32	- ⛔ **c1 build is PAUSED** — its build products (`S11c_c1_bulk_closure_sympy_audit.py` + `_exports.py`) were **DELETED** when the rebuild was stopped for this. **NEXT = the c1 MIGRATION:** (1) amend the c1 **spec §7 + header + `N1`** (`S11c_decisions.md`) from "carry the whole inherited LEDGER" → "write an **own-rows delta**; import via the **fold**"; (2) update the **c1 SymPy build directive** to §D1–D3 (own-rows delta, explicit `IMPORT_KEYS`, **digest-pin `ledger_fold.py`**, import via the fold); (3) **rebuild c1's export as a delta** (few-MB, over the frozen base) + 2 legs; THEN resume c1's dual-engine physics (blind WL engine → T7 comparator → reconcile → step record) per the clause below.
    33	
    34	## ⭐⭐⭐ S11c-c STARTED — 2-way SPLIT (c1 bulk closure / c2 self-energy fold); c1 spec COMMITTED `db5cbf88` (2-leg-reviewed + folded); c1 SymPy build directive `bd8b8cef` → build PAUSED for the export redesign above (2026-09-03)
    35	✅ **The old S11c-c row is SPLIT into two build units** (N2 refinement, decided WITH the two spec-review legs — `directives/S11c_decisions.md` N2 permits it): **c1 = curved bulk closure** (the nonlocal DtN/impedance + permeable B0c face response; exports the closed face response `(δp_s,J_s,t_s)(V_s,μ_θ)` + DtN operator); **c2 = self-energy fold** (fold the closed response into the S11c-b slab operator + re-extract the closed off-diagonal kernel → the coupled nonlocal self-energy operator). The split is documented in the c1 spec header + `_measurements/S11c_c_spec_review.md`.
    36	- **Process (rule-7 gate worked):** v1 MONOLITH spec drew **~13 verified defects** across 2 legs (Codex+Grok, each with saved SymPy derivations) → verdict split (Codex: re-author+3-way; Grok: fold+1-spec) → **USER chose the 2-way split**. c1 spec authored folding all 11 c1-relevant findings → **c1 got its OWN 2 legs** (Codex+Grok, fresh derive-and-check + split-seam) → BOTH: **"safe after the listed folds, NOT a re-author"**, core object + split CONFIRMED with derivations → **10 folds applied ONE pass** (rule 7: fold and go, NOT re-legged).
    37	- **c1 spec `directives/S11c_c1_SHARED_PHYSICS.md` COMMITTED `db5cbf88`** (the committed spec is the build input; the SymPy engine's `BUILD_INPUT_DIGESTS` will pin it — no c1 engine/exports currently present or committed; the first build's 145 MB products were DELETED for the export redesign). Key folds locked: DtN as a **two-momentum operator** `Z(ω;k,k′)` (composition `N₀∘M_h∘N₀` + kernel, both branch legs, named rigid-shift; ⛔ no left-quantized single-leg — rule 17); dissipation as **3 distinct objects** (bulk-radiation Hermitian, two-port permeable Hermitian, + an **independent traction-vs-far-field-flux** energy balance — the typed `½Re(δp V*)` was toothless, = bulk flux, never sees `t_s`); reactive part `K_a`; no sign-definiteness on the evanescent nullspace at 1st order (O(η²) leakage → S11c-e); `Σ_E`/`μ_R,bg` **structurally absent from c1's μ_θ-symbolic objects ⇒ reserved for c2**; rep-invariance route 2 = **Hanzawa** (not the secular global scaling; drop N4/MATERIAL on the bulk DtN); zero-jet target = **UNMODIFIED S11b half-space `Z`** (not a re-solved cavity); both branch freezes; Fredholm-vs-algebraic loci; non-uniform grazing domain (grazing = strict `v_bulk_normal_0=0`); N1/N8 chain output.
    38	- **c1 SymPy build directive DONE `bd8b8cef`; the first Codex build produced a 145 MB export → the export redesign (clause ABOVE); c1 build PAUSED.** NEXT = the **c1 MIGRATION** (per the EXPORT ARCHITECTURE clause: amend spec §7/header/N1 + directive to §D1–D3 + delta rebuild) → then resume c1 dual-engine: **blind Wolfram** engine → 2 build legs (fresh Claude + Grok) → **T7 comparator** join → reconcile → c1 step record. THEN c2 (held folds: extract-then-close — close the FULL operator first then re-extract, extract/eliminate DON'T commute; θ-row `Λ_X`/`J_s` row-routing; the **substitution-increment** isolation from S11c-b's deferred residual). ⚠ NO per-substep card (one S11c roll-up after S11c-e).
    39	- ⚠ **MEMORY:** S11c-c is heavier than S11c-b (which OOM'd the 30 GB box). Manage as S11c-b did — SERIALIZE the 2 engines, watch RSS, defer heavy controls in-band → out-of-band (`DEFERRED_HEAVY_RUNS.md`); the full cross-engine residual may again need the ≥64 GB box. ⛔ NEVER run 2 memory-heavy CAS jobs concurrently. Anchors: c1 spec `db5cbf88`, review `_measurements/S11c_c1_SHARED_PHYSICS.md` + `_measurements/S11c_c_spec_review.md`.
    40	
    41	## ⭐⭐⭐ S11c-b CLOSED (per-engine leg-verification) — step record committed `bcb9f7d7` (2 legs); cross-engine residual DEFERRED to ≥64 GB; NEXT = S11c-c [DONE — S11c-c STARTED, see clause above] (2026-09-03)
    42	✅ **S11c-b step record `steps/S11c_b_variable_coefficient_operator.md` committed `bcb9f7d7` + pushed** — the honest close: basis 40 + un-frozen operator + constraint-fold pin B + #90 face+response coupling, all PER-ENGINE leg-verified; the FULL cross-engine `row_residual` (P2a/P2b comparator builds + WL primaries `.out` with term-origins) + #88 re-adjudication + 2 control-hardenings DEFERRED to a ≥64 GB box (user choice — the lighter core-only ~8 GB path was offered, not taken). The 2 whole-row sign conventions (kinetic, face) + #90's 2 flags stay cross-engine-unvalidated. Both document legs (Codex+Grok) confirmed the per-engine-vs-deferred split is HONEST; ~16 accuracy/provenance defects folded (review `_measurements/S11c_b_step_record_review.md`). ✅ `exports.py` regen COMMITTED `af560257` (faithful — `BUILD_INPUT_DIGESTS` match the committed folded+#90 engine/spec; plain-git blob; per user direction). ⚠ its coupling content is per-engine leg-verified, NOT yet cross-engine-validated ⇒ S11c-c imports per-engine-verified coupling. Committed in-tree `.out` remain STALE (fresh ones live in `~/.s11_build/` scratch; regen is ≥64 GB work). **NEXT = S11c-c (curved-interface bulk closure / DtN)** — `directives/S11c_decisions.md`; NO per-substep card (one S11c roll-up after S11c-e). ⚠ Do NOT run 2 memory-heavy CAS jobs concurrently (OOM lesson). Anchors: step record `bcb9f7d7`, deferral `66e8d021`, P1-WL `06048d15`.
    43	
    44	## ⭐⭐⭐ S11c-b CROSS-ENGINE RESIDUAL DEFERRED to a ≥64 GB box (STEP 0 OVERTURNED); S11c-b CLOSES on per-engine leg-verification; earlier NEXT was the honest step record (2026-09-03)
    45	⛔ **STEP 0 was WRONG.** It measured only `evaluatedModel` + one `FINAL_KERNEL` (7.95 GB) and concluded the residual fit 30 GB. The full `S11CB_PRIMARIES_ONLY` production emit ALSO builds `mainKernelOrigins` (a FULL `extractCoupling` per origin, ~6/case, for `COUPLING_KERNEL_TERM_ORIGINS`) → ONE case hit **15.6 GB still growing** → OOM. ⇒ **the full cross-engine residual is DEFERRED to a ≥64 GB box** (user decision 2026-09-03: defer whole residual, do not take the lighter core-only path). S11c-b closes on the per-engine leg-verification (fold + #90 + #89/a/b, all 2-leg CLEAR) + the coarse cross-engine consistency established 2026-09-02. The P2a/P2b comparator decision lists are committed as the SPEC for the ≥64 GB run (`DEFERRED_HEAVY_RUNS.md`). ⚠ The two whole-row SIGN conventions (kinetic −K PY vs +K WL; face) + #90's two flags (closure-fold sign; uniform-Λ survivors) remain cross-engine-UNVALIDATED (deferred). NEXT = honest step record `steps/S11c_b_*.md` (2 legs) documenting the arc + the deferrals; exports handling; then S11c-c. Below is the (now superseded) in-flight plan, kept for the artifact map.
    46	
    47	User chose the FULL reviewed cross-engine residual, THEN (after the OOM finding) chose to DEFER it. The pass was a P1-WL emit gate + two comparator builders (P2a slab-row/row_residual, P2b §3a bridge) → run the residual — the gate + decision lists are DONE; the builds + run are ≥64 GB work.
    48	- ✅ **P1-WL gate committed `06048d15` + pushed** (origin+GIN synced): additive +17-line `S11CB_PRIMARIES_ONLY` in the WL engine (mirrors PY's), set-mode emits primaries for all 4 cases + skips controls/tower, unset byte-identical to HEAD. **2 build legs CLEAR** (⚠ TWO fresh Claude agents — Grok substituted for the gate's 2nd leg due to an xAI capacity outage; low-risk additive change, objective checks: 14 tags/0 controls, extractors pass, payloads byte-identical set-vs-HEAD, both FORM ablations bite). v1→v2 after 2 decision legs rejected a single-case selector (branch-scoped energy families; `row_residual` raises on non-aligned keys → the fix is the symmetric PRIMARIES_ONLY, both engines emit matching 4-case primaries).
    49	- ✅ **Fresh PY `.out` done** (`~/.s11_build/S11c_b_step0_scope/py_primaries_fresh.out`, 183 MB, folded+#90, PRIMARIES_ONLY). 🔄 **Production WL `.out` running** (`~/.s11_build/S11c_b_production/wl_primaries.out`, task `bgxqpd4ok`, memory-watched, ~8 GB/case as STEP 0 predicted, ~2h). These are the residual inputs (scratch, not committed).
    50	- ✅ **P2a v2 (slab-row join + row_residual) committed `06048d15`; P2b v2 (§3a bridge) committed `69bcc18d`.** Each folded from 2 decision legs (Codex+Grok, ~10 computation-backed defects each). KEY findings: (P2a) `row_residual.py` is in the blast radius — its #90-STALE one-sided face subtraction (`:427`) manufactures a false disagreement; μ_θ not consumed; `DOF` must be pinned; the closed disposition table is now IN the directive; origin-family name-mismatch adapters; the TWO whole-row SIGN conventions (kinetic −K PY vs +K WL; face) are SURFACED in the residual, NOT normalized (adjudicated in the step record — never design away the disagreement). (P2b) the bridge needs an EXPRESSION-valued scale substitution (`gammaWJET*→W_0·gamma_s11cb_*`; `I_PY=W_0·I_WL`), not a string rename; a factor-LOCKED energy-term bijection certificate; the applied→bare guard on the FULL parse (width/modulus/eW), pinned to S11c-b tables (not s11ca FIELD → S11c-a regression). ⛔ Sequencing: P2b lands BEFORE P2a's final validation.
    51	- ⛔ **OOM LESSON (2026-09-03):** ran the WL production run (~8 GB) + the P2b Codex build (parsing the full 156 MB WL + 84 MB PY `.out` into SymPy = tens of GB) CONCURRENTLY → system OOM killed BOTH (watchdog never fired — min avail 8.6 GB; no dmesg access; memory recovered to 23 GB after). ⇒ NEVER run two memory-heavy CAS jobs concurrently on this 30 GB box; SERIALIZE. The incomplete P2b build changes were discarded (committed baseline intact); the P2b directive gained a memory-bound build note (test against a §3a-invariant-only / single-case extract, not the full pipeline).
    52	- ✅ **exports.py REGEN COMMITTED `af560257`:** the PY PRIMARIES_ONLY run republished `S11c_b_exports.py` (hash-chained input; faithful deterministic regen of the folded+#90 engine — digests verified against the committed inputs). Committed per user direction (⚠ coupling content is per-engine-verified, cross-validation deferred to ≥64 GB).
    53	- **NEXT (serialized, memory-safe):** production `.out` completes → **P2b build** (Codex, bounded, alone) → 2 build legs (fresh Claude+Grok) → commit → **P2a build** (alone) → 2 build legs → commit → **P3: run `row_residual --py <fresh PY> --wl <fresh WL>`** → the cross-engine residual (fold's constraint-reduced rows + #90's closure-fold sign / uniform-Λ survivors; the whole-row sign conventions surface here for step-record adjudication) → #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → honest step record `steps/S11c_b_*.md` + exports (2 legs) → close of S11c-b. ⚠ NO per-substep card (one S11c roll-up after S11c-e); S11c-c/d/e remain. Records: `_measurements/S11c_b_{p1_wl_residual_emit,p2_comparator_update,p2a_slab_row_join,p2b_gamma_bridge}_directive.md`. Post-06048d15 anchors: P1-WL `06048d15`, P2b `69bcc18d`.
    54	
    55	## ⭐⭐⭐ S11c-b CROSS-ENGINE STEP 0 (memory scope) DONE — the single-case residual FITS this 30 GB box; integration pass is a 3-builder sequence; NEXT = P1 residual-mode single-case emit (2026-09-03)
    56	The resume-prompt STEP 0 is measured and recorded (`directives/_measurements/S11c_b_step0_residual_scope.md`; `DEFERRED_HEAVY_RUNS.md` reconciled). A guarded single-case WL run of exactly the objects the residual needs — `evaluatedModel` (operator) + `extractCoupling…["FINAL_KERNEL"]` (primary kernel), NO tower-depth control variants — peaked **7.95 GB RSS / 0.99 GB in-kernel, min MemAvailable 14.94 GB, ~26 min** for one case. ⇒ the cross-engine residual is **doable on this box**; the ≥64 GB requirement bounds only the full 4-case in-band `.out` regen with the tower/heavy controls (built unconditionally at `…mathematica_audit.wl:2204-2231`), which the residual does not need.
    57	- ⚠ **REFRAMING (verified):** `row_residual` and the comparator PARSE the two committed `.out` (`--py`/`--wl`, default the committed files) — they do NOT run the engines live. **Both committed `.out` are STALE:** WL at `d4adbd99` (#89a, pre-#89b → FROZEN operator) and PY at the migration checkpoint (PRIMARIES_ONLY, pre-fold, pre-#90). The stale PY `.out` even has the pre-fold `THETA_BALANCE = μ_θ` structure. ⇒ the integration pass MUST regenerate **fresh single-case `.out`** from both current engines first, then point the comparator at them.
    58	- ⚠ **STEP 0 OPERATIONAL DISCOVERY:** the engine's normal emit path always builds the tower variants (~16 GB/case; only `S11CB_SKIP_HEAVY_CONTROLS` gates the 2 equivalence controls, not the variants). A comparator-PARSEABLE single-case `.out` at the measured ~8 GB footprint needs a **residual-mode single-case emit** — restrict to one case, call the engine's OWN `emitShared`/`modelRecord`/`kernelRecord` on the primary objects, skip the tower variants. Same objects, same emit functions.
    59	- ⭐ **BRIDGE DELTA MAPPED (agent) + rule-13 VERIFIED (me):** (a) `extract_slab` has a now-FALSE `THETA_BALANCE → ("MU_THETA","THETA")` bridge — post-fold that row is `evolution_mass_balance − Σ closure_shape_deriv` (mass-evolution), so it would join PY's mass-evolution against WL's raw μ_θ = a MANUFACTURED disagreement; plus `MASS_EVOLUTION_ROW` under-bridged (PY maps only the ADVECTIVE summand) and WL `CENTER_FACE_GENERALIZED_ROW` unpartnered. (b) The §3a coefficient bridge is FULLY STALE: WL renamed its coefficients to `gammaWJET*`/`gammaMUJET*` at `d4adbd99`, AFTER the bridge tables were written (`70164909`), so the 12 `gammaWidth*` entries now match NOTHING ⇒ **0 of ~30 §3a coefficients are currently bridged**; 7 suffixes/source never had an entry. ⛔ PY coefficient names are POSITIONAL (`gamma_s11cb_*_NN`, runtime-quotient-selected) — the pairing MUST be read off the invariant each coefficient multiplies (PY emits `ENERGY_BASIS_NEW_INVARIANTS`), NEVER guessed by index. (c) Profile-jet + DOF bridges are COMPLETE — no gap.
    60	- **NEXT = the integration pass, a 3-builder sequence (each leg-gated):** **P1** add a residual-mode single-case emit to BOTH engines (primary tags only, tower/heavy controls skipped) → **P2** update the comparator (`extract_slab` structural fix + §3a bridge by invariant-matching) — its OWN decision list (drafted, 13 testable properties) + 2 decision legs + 2 build legs → **P3** run `row_residual` on the fresh single-case `.out`, read the residual (OUR finding, no value gate) — validates the fold's constraint-reduced U/E_W rows AND #90's closure-fold sign/uniform-Λ-survivors. Then #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → step record + exports (2 legs) → close of S11c-b. ⚠ NO per-substep card (one S11c roll-up after S11c-e); S11c-c/d/e still remain.
    61	
    62	## ⭐⭐⭐ S11c-b #90 §3c COUPLING CONTENT BUILT + leg-gated + committed `7677aa18` (pushed origin+GIN) — the "first physics number" (face+response) is IN; NEXT = cross-engine integration pass (this box) (2026-09-02)
    63	PY's coupling was bulk-only; #90 folds the settled §3c INCLUDE/INCLUDE content — the reversible tilted-face geometry + the irreversible face response (`Λ` symbolic). Record `directives/_measurements/S11c_b_90_coupling_content_directive.md`; directive `directives/S11c_b_90_coupling_content_directive.md`; leg evidence `~/.s11_build/S11c_b_90_coupling/`.
    64	- ⛔ **Decision legs (Codex+Grok) REJECTED v1 architecturally:** `FACE_FLUX_BOUNDARY_OPERANDS` is the RAW T-substrate bundle, not operator rows — weak-restricting it is the §3c-forbidden PARALLEL ROUTE. Fix (both legs, convergent): compute the face GENERALIZED-FORCE rows from the consumed virtual work (coeffs of `δ_vu`/`δ_ve_W`, live `μ_θ` bound), ADD to the constraint-reduced operator rows (origin FACE, not through the θ-fold), then weak-restrict the full operator; fold the SKIPPED `closure_shape_deriv` (Λ_A/Λ_V) + `traction` (Λ_X). Plus the `A_T` token collision (PY's `A_T_s11cb` is the test potential, not the geometry), double-count/over-reach traps, §0 "every" overstatement, rule-5 leaks. v2 folded both (one pass).
    65	- ✅ **Build legs (fresh Claude agent + Grok) BOTH VERDICT CLEAR (8/8 probes, form ablations BITE):** face computed INTO the operator (`U/E_EXPANDED = bulk + face`; no `FACE_FLUX_BLOCKS`; `build_kernel` never reads the raw bundle; Grok's independent face-row reconstruction matches the emit); form-ablating disjoint face sources moves disjoint kernel channels, non-face (projection) moves nothing; `Λ` symbolic (`Λ_I⁰/(1−iωτ_I)`, no `Z`/DtN/impedance); `μ_θ` bound (no reserved `mu_theta_L/M` in blocks); no `ζ_c`; adjointness over the enlarged blocks; exact-once; `operator_from_density`/`committed_strong_rows` byte-identical (#88 refs). ⚠ Codex's OWN self-review (Codex variants, "Claude/Grok unavailable") was DISCARDED as invalid.
    66	- ⚠ **TWO cross-engine / step-record flags (NOT build defects — the construction is correct per §3c; a single engine can't settle the response CONTENT):** (1) the closure-fold sign/magnitude (Grok's T-i identity `TRUE_AREA − RELFLUX_SUM = 0` corroborates the recipe, but sign/magnitude is PY↔WL-only); (2) the §5c uniform-limit residual is NONZERO and Λ-bearing (`(2,4,0,0,0,0)`, γ-count 0) — whether the irreversible face-response legitimately survives at `(η,σ_W)=0` or would violate §1d's uniform decoupling needs the cross-engine residual + step-record (does WL carry the same Λ survivors?). §5c is a smoke test, no value gate.
    67	- **NEXT:** the **cross-engine integration pass**: complete the S11c-b symbol bridge (~35 §3a basis coeffs + jet + DOF transliteration) + update `extract_slab` to the folded PY structure (`THETA_BALANCE`→mass-evolution-minus-closure, `MU_THETA` from `MU_THETA_FACE_BINDING`) → run `row_residual` (reviewed cross-engine JOIN) — this validates the fold's constraint-reduced rows AND #90's closure-fold sign/uniform-survivors → #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → honest step record `steps/S11c_b_*.md` + **exports** (2 legs) → STATUS/memory close **of S11c-b**. ⚠ **The card:** S11c-b owes NO separate card (S11c-a precedent, `steps/S11c_a_*.md:268-269`); N1 (`directives/S11c_decisions.md:24-31`) specifies ONE S11c roll-up card, produced ONLY after S11c-e. ⚠ S11c-b is ONE sub-step — **S11c-c** (curved-bulk closure/DtN), then **S11c-d/e** still remain. ⚠ **RESOURCE SCOPE (Codex-corrected):** the measured **0.9 GB fits only the single-case WL OPERATOR/U-row probe** (`evaluatedModel`), NOT the #90 KERNEL or a full comparator pass — the #90 kernel builds took the legs ~45 min/case, and `DEFERRED_HEAVY_RUNS.md:58-84` still requires ≥64 GB for full WL operator/kernel regen. So the cross-engine pass's feasibility ON THIS BOX is established only for the U-row operator; the kernel/comparator scope must be MEASURED before claiming it, and `DEFERRED_HEAVY_RUNS.md` reconciled.
    68	
    69	## ⭐⭐⭐ S11c-b PY CONSTRAINT-FOLD BUILT + leg-gated + committed `82f53828` (pushed origin+GIN) — pin (B) implemented; NEXT = single-case cross-engine `row_residual` (this box) → #90 (2026-09-02)
    70	The PY engine now folds pin (B): the slab `U`/`e_W` rows are the CONSTRAINT-REDUCED equations (the imported non-uniform `virtual_constraint` eliminates virtual θ, the same held-fixed `μ_θ` feeds both reactions), the θ-row is the imported sourced `evolution_mass_balance` (⚠ AT `82f53828`; **#90 SUBSEQUENTLY changed the θ-row to `evolution_mass_balance − Σ closure_shape_deriv`** — folding the face-response closure into the mass row; current HEAD's θ-row is the mass balance MINUS the closure residual, and the U/e_W rows are constraint-reduced internal + face generalized forces), `μ_θ` stays a separate held-fixed operand, and the jet-depth cascade is raised (`STRONG=3, COUPLING=4, DEPTH_CONTROL=5`). `operator_from_density`/`committed_strong_rows` are byte-identical to HEAD (#88 raw reference preserved). Record `directives/_measurements/S11c_b_py_constraint_fold_directive.md`; directive `directives/S11c_b_py_constraint_fold_directive.md`; leg evidence `~/.s11_build/S11c_b_constraint_fold/`.
    71	- ⛔ **The 2 decision legs REJECTED v1 for a rule-17 FREEZE** — I handed §1c:143's `(uniform linearisation)` as the constraint (froze `W_bg`); the binding object is the MATERIAL `δ_vΣ_mat=0` with live `W_bg` (imported as `virtual_constraint`). Plus a θ-row double-count, a #88-break interface, a missing depth cascade, and provenance/leak. v2 folded both legs (rule 7, one pass). The gate caught the cardinal error before one engine line changed.
    72	- ✅ **2 build legs (fresh Claude agent + Grok): PHYSICS ROWS CLEAR on both** — the reaction is COMPUTED from the imported non-uniform constraint (a FORM ablation of the constraint source MOVES the rows; `live−ε∇μ_θ ≠ 0` but `uniform−ε∇μ_θ = 0`, proving the v1 freeze was NOT committed), θ = imported mass-evolution with no double-count, depths load-bearing (cap→2 loses all order-3; coupling cap→3 loses 174 order-4 terms), #88 raw refs intact, provenance `BULK_ENERGY`, no leak, smoke PASS (564s, 4 cases).
    73	- ⚠ **ONE FINDING folded (Claude agent Probe 2, orchestrator CAS-proven rule 13):** the two-route (elimination vs Lagrange) residual is TAUTOLOGICAL — `λ=−ε·μ_θ/∂_θC` makes the routes algebraically identical for any linear constraint, so the residual is `0` by construction even for a WRONG constraint (operand theatre, not independence). Relabeled `ROUTE_RESIDUAL → CONSTRAINT_FOLD_TRANSCRIPTION_RESIDUAL` + the directive claim corrected; the rows are unchanged. The legs DISAGREED (Grok CLEAR via a weaker CODE-corruption test; the agent's WRONG-INPUT test was decisive). Real independence = the cross-engine comparator. **FLAG deferred to #88:** `HESSIAN_FREEZE_STRONG_ROW_RESIDUAL` now compares folded-live vs raw-committed ⇒ is the reaction (nonzero); nothing asserts it.
    74	- ⭐ **CROSS-ENGINE ATTEMPT (2026-09-02) → KEY DISCOVERY: the ≥64 GB box is NOT needed for the residual.** A guarded single-case WL U-row run (EULERIAN/LAB_HELD/RHO4_CONSTANT, full basis) FINISHED at **0.9 GB peak** (not 16 GB) — the ~16 GB is only the full run's tower-depth control variants (`operatorActivated/Truncated/Extended` + kernel variants), which the residual doesn't need. ⇒ the whole cross-engine integration (single-case AND likely full-4-case primaries) is doable **on this 30 GB box**, gating only the tower controls. Coarse cross-engine consistency established (both U-rows order-3; both carry live-`W_bg` `eta_bg*w1` coupling — neither froze the constraint). ⚠ RULE 4: a prose shortcut (a clean shared `(η w1 −1)` factor) was REFUTED by the all-jets check (72/75 PY jet coeffs not cleanly divisible). A real COEFFICIENT-level residual needs the full S11c-b symbol bridge (~35 §3a basis coefficients + jet + DOF transliteration) + the `extract_slab` update to the folded PY structure — a reviewed physics-bearing cross-engine JOIN (⛔ never blanket-collapse), doable HERE, deferred to a dedicated integration pass. Evidence + record `~/.s11_build/S11c_b_constraint_fold/{cross_engine_single_case_attempt.md, wl_urow_labheld_rho4.txt, py_urow_labheld_rho4.txt}`.
    75	- **NEXT:** **#90** — PY §3c coupling content fix (the under-extracted reversible tilted-face `A_T` + irreversible response `A_T·Λ(ω)`, per the #84 SETTLED verdict) + a §0 clarity pin (Λ = supplied flat-face closure, NOT the S11c-c bulk kernel), on the corrected basis + folded operator. Then the cross-engine integration pass HERE (transliteration + `extract_slab` update + `row_residual`, reviewed) + #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → step record + S11c card + close. ⚠ The ≥64 GB box is NO LONGER on the S11c-b critical path (only the belt-and-suspenders full in-band control run remains big-box work).
    76	
    77	## ⭐⭐⭐ S11c-b STRONG-ROW JET-DEPTH RECONCILED → SPEC PINNED (B): the slab U-row is CONSTRAINT-REDUCED; **PY is the engine that must change** (WL correct). NEXT = the PY constraint-fold BUILDER (own decision list + build legs) (2026-09-02)
    78	The #89b PY-check flag ("WL emits order-3 in strong U-rows, PY caps at depth-2 — PY under-emits?") is RESOLVED, and NOT as first thought. Reconciliation record `directives/_measurements/S11c_b_strong_row_jet_depth_reconciliation.md`; three consistent computations (orchestrator PY probe + Codex + Grok consult, evidence under `directives/_measurements/s11c_b_jet_depth_consult_{codex,grok}/`) showed it is NOT a jet-depth freeze and NOT a physics error — a REPRESENTATION mismatch: the §3a energy has only FIRST background jets, so the RAW held-fixed `δU/δu` is order-2-complete on BOTH engines; WL's order-3 is ENTIRELY the θ-constraint reaction `∇μ_θ` (its `constrainedRowsWithLiveEnergyEL` eliminates virtual θ via the MATERIAL VIRTUAL CONSTRAINT `δ_vΣ_mat=0` [linearized] — distinct from the separate sourced mass-EVOLUTION equation), which PY keeps SEPARATE (raw EL + a separate `μ_θ` operand).
    79	- ⛔ **SPEC-PINNED (B) — `directives/S11c_b_jet_depth_spec_pin_decision.md` (FOLDED VERDICT), 2 decision legs Codex+Grok BOTH B, raw verdict transcripts `~/.s11_build/S11c_b_jet_depth/{codex,grok}_leg_verdict.log` (`.log` gitignored).** My proposed (A)/PY-correct lean was OVERTURNED by the binding S11b rule (orchestrator-verified verbatim, rule 13): `S11b_SHARED_PHYSICS.md:341-342` "`δ_vθ,δ_v(δW),δ_vu` not independent … **Do NOT vary `U` with `θ` held fixed**" (constraint eq at `:337`) + `:426` "the in-plane equation must carry `−∇(δU/δθ)` … varying at fixed θ removes this contribution … selects the convention uniquely"; S11b's own engine `…coupling_law_mathematica_audit.wl:280` `constrainedUL = explicitUL + I k muTheta`. §1c's "θ may not be eliminated before this derivative" scopes the CONSTITUTIVE `μ_θ`, not the U-row; the separate `μ_θ` operand is a scalar face-affinity driver (not double-counting the vector `∇μ_θ`); θ's real EOM is MASS-EVOLUTION, not `μ_θ=0`; (A) would also break the §5c uniform-limit regression vs S11b's `INPLANE_EOM`.
    80	- **BUILDER TARGET (next phase; its OWN decision list + build legs, rule 7/9 — no engine change yet):** PY `operator_from_density` (`scripts/S11c_b_brane_operator_sympy_audit.py:1968-2062`, currently the raw held-fixed EL) MUST fold S11b's material virtual constraint `δ_vθ+δ_ve_W+∇·δ_vu=0` into `U_MOMENTUM_ROWS` + thickness row; raise `STRONG_ROW_JET_DEPTH` 2→3 (NOW a genuine rule-17 freeze; PY already reproduces the 10 order-3 atoms once the constraint is applied); keep `μ_θ` a separate constitutive operand; θ-row = mass evolution. Add a §3b sentence (both legs drafted it). **WL unchanged** (`constrainedRowsWithLiveEnergyEL` is already this object; #89b's un-freeze was right). **#88 RE-ADJUDICATION:** #88 identified the strong U-rows with the raw held-fixed EL — that conflicts with the pin; its energy-basis-completion disturbance measurement STANDS, but the full-row/KINETIC adjudication is redone after PY carries the constraint reaction, and #88's θ result is a disturbance of `MU_THETA_OPERATOR`, not an independent θ equation. Then the deferred `row_residual` compares the two constraint-reduced (order-3) U-rows.
    81	
    82	## ⭐⭐⭐ S11c-b #89b WL OPERATOR UN-FREEZE + REPAIR DONE — engine committed `a1be8d8f` + PUSHED (origin + GIN); full `.out` regen DEFERRED to a ≥64 GB box; PY sibling CLEAR (no freeze); NEXT = integration (reconcile strong-row jet depth) → #90 (2026-09-02)
    83	The §3b variable-coefficient slab operator is un-frozen: coefficients stay LIVE through every order-raising step (EL + `activateSpatialDivergences` + `extractCoupling` + face EL + constraint), `Inactive[Div]` split preserved, jet-retaining reduction LAST (rule 17). Tractability (naively HUNG 2h+ → ~6 min): innermost-first Div activation (no `If`-on-Association) + per-summand `Series` linearity, `PossibleZeroQ`-verified. ⚠ **MEMORY WALL — INTRINSIC:** the correct un-frozen operator must hold its full jet tower un-reduced (~16 GB/case) until the final reduction, so BOTH the 2 heavy equivalence controls (behind `S11CB_SKIP_HEAVY_CONTROLS`) AND the whole 4-case in-band `.out` regen are DEFERRED to a **≥64 GB box** (`research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md`); the committed `.out` is UNCHANGED — **the deliverable is the ENGINE.**
    84	- ⛔ **A RE-FREEZE BLOCKER was caught by the 2 build legs and repaired.** The EMITTED operator (what the cross-engine comparison reads) reduced `operatorLive` while its outer `Inactive[Div]` was still HELD → froze widthBase inside the Div → dropped the mixed-2nd/3rd U-row Hessian jets (the correct `activatedOperator`, computed at `:1377`, was DISCARDED). Plus 2 broken controls (§5.E dim walker non-functional [`Times`/`Plus` Flat/OneIdentity]; MATERIAL_ADVECTED independence → `base−base=0`). ⭐ The 2 legs DISAGREED — Grok caught it, the fresh Agent CLEARED it (its atom-PRESENCE test saw the jet inside the held Div but never the ACTIVATED coefficient); the orchestrator resolved it by its own computation on the real emit (rule 13, activate & reduce do NOT commute on the U row) → Grok right.
    85	- **REPAIR** (`directives/S11c_b_89b_wl_operator_repair.md`): 2 decision legs (Grok 7 + Codex 4 = **11 directive-level gaps** folded — the fix ripples to the tower-depth operands, the frozen Hessian-witness across ALL slots, the uniform-limit SLAB **and** TRANSVERSE_DISPERSION, plus a rule-5 acceptance leak that could pass the frozen object); Codex repair (+28 lines, scope-verified = only the 3 fixes + ripples); **2 re-review legs (Grok + fresh Agent, rebuilt from `operatorLive` not the emit) both VERDICT CLEAR** (emit == activate-then-reduce reference, one-sided swap MOVES the residual [order-3 jets restored], all comparisons like-with-like, dim walker moves SELECTIVELY under a primitive-atom mutation, independence carries `VALIDATED_ON_REPRESENTATIVE_BRANCH` with LAB_HELD live, no regression). Record `directives/_measurements/S11c_b_89b_wl_operator_build_review.md`.
    86	- **NEXT:** (1) ✅ **PY SIBLING FREEZE-CHECK DONE** (rule 16; record `directives/_measurements/S11c_b_89b_py_sibling_freeze_check.md`): PY is CORRECT (activate-then-reduce — emitted strong rows retain the Hessian, verified live-vs-frozen; un-frozen genuinely by #89's `total_derivative` tower), so the old agreement was NOT on a frozen PY object and **no PY freeze repair is needed.** ⚠ BUT a NEW reconciliation FLAG (spec question): PY caps strong rows at `STRONG_ROW_JET_DEPTH=2` (Hessian, no 3rd-order bg jets — 3rd lives in the coupling cascade) while WL #89b now emits **3rd-order** bg jets in the strong U-rows (the re-review's "order-3 restored") → `row_residual` WILL show a strong-row disagreement that is a **jet-DEPTH-CONVENTION question** (⛔ NOT pre-judged — it could be a genuine depth error in one engine OR a representational cap; adjudicate which strong-row jet depth the retained grade requires, don't assume WL over- or PY under-emits) BEFORE reading the residual as a physics finding. (2) **Integration** — rebuild both `.out` (⚠ WL side needs the ≥64 GB box) + re-run `row_residual` (~15 min) + **re-adjudicate KINETIC+θ** (#88 blast radius) + fold 2 owed #88 control-hardenings. (3) **#90** PY §3c content fix (face+response) + §0 pin on the corrected basis. (4) honest step record `steps/S11c_b_*.md` + S11c roll-up card + close. ⚠ Committed with EXPLICIT paths (16 files); the other session's repo `memory/` staged work was left untouched.
    87	
    88	## ⭐⭐⭐ S11c-b #89a WL §3a BASIS REPAIR DONE — WL emits basis 40 (span verified = independent derivation), CLEARED by 2 build legs, COMMITTED `d4adbd99`; #89b WL OPERATOR unfreeze NEXT (2026-09-01)
    89	The WL side splits: **#89a (DONE) = the §3a energy BASIS; #89b (NEXT) = the operator unfreeze.** #89a completed the WL engine's hand-coded 8/source spurion family to the full **O(3)-Kronecker** field-bilinear family (`directives/S11c_b_89_wl_3a_repair_directive.md`), imposing the exact thickness map on the new invariants; the count is computed by `MatrixRank` (`26` = the restrict-to-8 regression, the only public target). ✅ **CLEARED by 2 build legs** (fresh Claude agent + Grok, both VERDICT 1-finding/0-blockers): two INDEPENDENT SymPy derivations + the engine all agree the completed basis is **40 = 10 uniform + 15 + 15**, and — the decisive check — the engine's 15/source **span the SAME space** as the independent derivation (`rank(union)=15=each`; byte-identical scalars), ⛔ not a matching-integer coincidence. Form ablation load-bearing (restrict→26, +ε Levi-Civita→16 [parity exclusion real], drop a shear pairing→38/14); count computed-not-typed; thickness map imposed. Record `directives/_measurements/S11c_b_89_wl_basis_buildleg_clearance.md`; evidence `~/.s11_build/S11c_b_89a_wl_buildleg_{claude,grok}/`; decision-review record `_measurements/S11c_b_89_wl_decision_review.md` (2 decision legs, computation-backed, caught 8 defects incl. the operator-freeze miss + a tautological rule-17 control — folded once, rule 7). ⇒ **both engines (PY #89, WL #89a) now emit 40 on the CORRECT completed basis.**
    90	- ⚠ **OPERATOR FREEZE — DIAGNOSED, NOT repaired (scope split; = #89b).** WL's slab operator `evaluatedModel` applies `applyProfile` (frozen 2nd/3rd jets) BEFORE the EL variation → frozen-EL rank 26 vs live-EL rank 40 (Δ14; the Hessian is the non-absorbable operator content, the WL analog of the #88 blast radius). #89a emits this as an HONEST diagnostic (`OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC`, not asserted equal, no repair claim) and DEFERS the fix to **#89b** (EL-before-freeze then live-path reduction; a different mechanism/surface — bundling was a rule-15 risk). PY #89 already unfroze its operator, so the two engines reach operator parity only after #89b.
     1	# S11c-c2 self-energy fold — physics review adjudication (gate record)
     2	
     3	**Artifact reviewed:** `research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py` (1053 lines;
     4	astra/`gpt-6-astra` build). Emitted output `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out` (499 MB), navigated by
     5	`_measurements/S11c_c2_sympy_stdout_index.json` (154-tag byte-offset index).
     6	**Role:** script / physics-bearing build → **review-until-clear** (G2/G4), not one-pass.
     7	**Authorship → legs (G1):** script is Codex(astra)-written ⇒ two legs = **fresh general-purpose Claude agent +
     8	Grok**, orchestrator-launched, identical prompt. My own verification (below) is required (G4/rule 13) but is
     9	**not** a leg.
    10	
    11	**Identical prompt:** `directives/_legs/S11c_c2_physics_review_prompt.md` (13,386 bytes).
    12	**Blindness by absence (S1):** for the review's duration, astra's builder report + provenance narrative were
    13	relocated to `ext_logs/c2_review_hold/`; the legs saw only the mechanical byte-offset index. astra's ~36
    14	quarantined self-review files stayed out of tree (`ext_logs/c2_builder_overstep/`) and were disregarded. Held
    15	files restored before adjudication.
    16	
    17	**Handed inputs:** the script; the physics authority `directives/S11c_c2_SHARED_PHYSICS.md`; the real imports
    18	`scripts/{S11c_b_exports.py,S11c_c1_exports.py,ledger_fold.py}`; predecessor records
    19	`steps/S11c_c1_curved_bulk_closure.md`, `steps/S11c_b_variable_coefficient_operator.md`; the emitted `.out` + its
    20	index. **Withheld (by not handing):** the build directive (its frozen symbol map would make the review a fidelity
    21	check), the builder report, the provenance narrative, astra's self-reviews.
    22	
    23	**Leg reports (separately attributed):**
    24	- fresh Claude: `_measurements/S11c_c2_physics_claude_leg.md` (verdict: 0 confirmed physics defects; E unsettled/
    25	  deferred, F not-a-defect + caveat).
    26	- Grok (`grok-4.6` high): `_measurements/S11c_c2_physics_grok_leg.txt` (verdict: fold *wiring* right, but B/E/F/G
    27	  flagged CONFIRMED; "the increment as defined is not the S11b-decoupled self-energy the spec claims").
    28	
    29	---
    30	
    31	## The three disagreements — adjudicated by MY OWN computation (rule 13; legs wrong in both directions)
    32	
    33	Both legs independently agree A, C, D1–D6 and the fold **wiring** (δp not J_s, operator inverse, kernel bridge,
    34	V_s→face_velocity, computed jets, ε-strip) are SOUND, each with shown CAS evidence. The legs SPLIT on B/E/F/G. I
    35	resolved each myself; scripts + literal stdout committed beside this record.
    36	
    37	### F (uniform limit) — Grok CONFIRMED "genuinely nonzero coupling"; Claude "genuine coupling decouples"
    38	**Verified: Claude is right; Grok's F is a FALSE POSITIVE.**
    39	Command: `_measurements/S11c_c2_adjudication_verify_F.py` → `…_verify_F.out`. Method (independent of the leg's
    40	`Trial`-label heuristic): the increment = extract(close) − extract(open); extract(close) substitutes ALL δp
    41	slots ⇒ carries **no** bare `delta_p_±`; every bare-`delta_p` term therefore comes only from `−extract(open)`.
    42	Zeroing the bare slots leaves `extract(close)|_uniform`; `.doit()` the residual Integrals; test == 0.
    43	**Literal result (all 4 cases, both nonzero blocks THETA & E_W):** the surviving `closed_part` is
    44	`coeff·Test·Integral(0, …)` — integrand **literally 0** — `closed_part_zero(doit)=True`. So the genuine
    45	closure-induced coupling is **identically zero in the uniform limit**: the transverse↔thickness coupling
    46	decouples exactly as required. Grok saw `Integral(...)` present and reported a surviving `Z₀·μ_θ` coupling
    47	**without evaluating that its integrand is 0** (`expand` leaves `Integral(0,…)` unevaluated; only `.doit()`
    48	collapses it). (Grok's c1 "tangential freeze" was likewise a false positive — a pattern.)
    49	**Residue (real, not a defect):** the nonzero part of the emitted uniform-limit object is the
    50	`coeff·(δp_-+δp_+)·Test` bare open-slot bookkeeping from `−extract(open)` (δp is a free symbol in the open
    51	operator). **⇒ the increment does not literally vanish; the genuine coupling decouples.**
    52	
    53	### B (ε-grade) — Grok CONFIRMED "increment not O(εη)-only"
    54	**Same object as F; not a physics defect.** Both legs agree the strip-one is correct and there is no O(1) and no
    55	O(ε²) term. The O(ε) (η-free) grade Grok flagged **is** the `−extract(open)` open-slot bookkeeping (the F
    56	residue), i.e. a property of the §3c increment representation, not a fold error. The genuine induced self-energy
    57	is O(εη) (it decouples at η→0, verified in F).
    58	
    59	### E (N6 rep-invariance) — Grok CONFIRMED "does not vanish (σ channel)"; Claude "leading order holds, σ_W deferred"
    60	**Verified: leading-order rep-invariance HOLDS; the remnant is σ_W-sector only.**
    61	Command: `_measurements/S11c_c2_adjudication_verify_EG.py` → `…_verify_EG.out`. Independent test: does
    62	`REP_INVARIANCE_RESIDUAL.subs(sigma_W→0) == 0`? **Literal result (both densities, all 6 blocks):**
    63	`full_zero=False` but `sigmaW->0_zero=True` — setting σ_W→0 annihilates the entire residual ⇒ **every**
    64	non-invariance term carries a σ_W factor. So the leading O(ε) and O(εη) representation-invariance holds exactly;
    65	the non-invariance is confined to the **σ_W drain sector** (the known c1 drain-projection O(σ_W²) carry-forward;
    66	full evaluation matches the ≥64 GB deferred rep-invariance family). **⇒ not a load-bearing c2 defect;** surfaced +
    67	carried forward. The step record must claim rep-invariance only at leading order, with the σ_W remnant open.
    68	
    69	### G (adjointness) — Grok CONFIRMED "blocks not adjoint; suppressed check hides asymmetry"; Claude "honest omission"
    70	**Verified: SOUND; Grok's G is a FALSE POSITIVE (refuted a claim the builder never made).**
    71	Command: same `…_verify_EG.py`/`.out` (block directionality) + builder-report line 70.
    72	**Literal result:** in `SELF_ENERGY_INCREMENT`, the THICKNESS→TRANSVERSE reverse blocks are **identically zero**;
    73	only TRANSVERSE→{θ,e_W} is nonzero — the induced self-energy is **one-way**, physically expected since δp lives
    74	only in the θ/mechanical rows and depends on the transverse velocity (closing induces coupling *into* those rows,
    75	not back). `CLOSED_COUPLING_KERNEL` carries **both** directions nonzero (the pre-existing open coupling), and
    76	**both blocks are emitted** ⇒ the asymmetry is visible, not hidden. Builder report line 70 states the actual
    77	reason: *"There is no independent adjointness construction in this implementation, so no
    78	SELF_ENERGY_ADJOINTNESS_RESIDUAL is emitted. Both off-diagonal blocks are emitted."* — exactly spec §3b ("emit
    79	both blocks and say there is no independent second route rather than dressing a structural zero as a check"). The
    80	builder never claimed "adjoint by construction"; Grok attributed that and called it dishonest. Directional ≠
    81	defect; omission is correct.
    82	
    83	---
    84	
    85	## VERDICT — the c2 fold physics is SOUND (zero confirmed physics defects)
    86	
    87	The fresh-Claude leg's verdict holds. Grok's four "CONFIRMED" flags are all false positives (F, G) or a label
    88	dispute where leading order is clean (B is F's residue; E's leading order cancels). Both legs + my verification
    89	agree the load-bearing fold construction is physically correct.
    90	
    91	**Carry-forwards / caveats (step-record interpretation, ⛔ NOT build defects):**
    92	1. **F:** the emitted uniform-limit object is non-vanishing purely from the §3c `−extract(open)` bare-slot
    93	   representation. Interpret as **"the genuine closure-induced coupling decouples,"** ⛔ not "the increment
    94	   vanishes." ⇒ **light spec clarification owed** to `S11c_c2_SHARED_PHYSICS.md` §5e (the "increment must vanish"
    95	   wording) + §3c (the increment carries the open-slot O(ε) piece by construction) — a wording fix that does not
    96	   change what is computed; folds under spec review-until-clear.
    97	2. **E:** leading-order (O(ε), O(εη)) representation-invariance ESTABLISHED; the σ_W-sector remnant is SURFACED +
    98	   deferred (drain-projection O(σ_W²); ≥64 GB rep-invariance family). Record must not claim full rep-invariance.
    99	3. **G:** the induced self-energy is **directional** (transverse→{θ,e_W}; reverse block identically zero) —
   100	   expected from δp residing only in the θ/mechanical rows; both blocks emitted; no adjointness residual (correct
   101	   per §3b). Record it as a physics feature.
   102	4. The two S11c-b sign conventions (face-generalized-force, #90 closure-fold) are cross-engine-UNVALIDATED and
   103	   multiply the increment ⇒ the blind-Wolfram cross-engine comparator surfaces them (spec §3c/§7); the §3d.4
   104	   mechanical-power pairing adjudicates the face-force sign. (Deferred to the Wolfram engine + comparator.)
   105	
   106	**Next (unchanged by this adjudication):**
   107	- Commit this reviewed baseline (the script + review artifacts) **before** any repair overwrites it (this commit).
   108	- **STEP B — publication-only export repair** (still viable; physics sound): `publish`@807 / call@952 only —
   109	  drop `SELF_ENERGY_INCREMENT` (+ operands, term_origins, parity, §3d, §5 controls) to EMIT-only; export **both**
   110	  `S11CC2_CLOSED_SLAB_OPERATOR` + `S11CC2_CLOSED_COUPLING_KERNEL` (all 4 cases) transparent-factored + the casewise
   111	  `canonicalize(expanded − decode(compact))==0` check. ⛔ do not touch build_case/extract/kernel/payloads.
   112	- **STEP C** — re-review the export repair (publication-only) until clear → commit → then the light §5e/§3c spec
   113	  clarification (review-until-clear) → the blind Wolfram engine → its 2 legs → c2 T7 comparator + reconcile → c2
   114	  step record.
     1	# S11c-c2 export repair — re-review adjudication (gate record)
     2	
     3	**Artifact:** the astra publication-only export repair — the diff of
     4	`scripts/S11c_c2_selfenergy_fold_sympy_audit.py` from the reviewed baseline `8f3a017f` (publish rewrite +
     5	`publication_compact` helper + a strict-recursive semantic guard; `EXPORT_ROOTS`/`export_key` membership) and the
     6	regenerated `scripts/S11c_c2_exports.py` (**60,516,900 → 22,441,522 bytes**, 63% smaller).
     7	**Role:** script / publication-only build → review-until-clear. **Authorship → legs (G1):** astra(Codex)-written ⇒
     8	fresh Claude agent + Grok. Identical prompt: `directives/_legs/S11c_c2_export_repair_rereview_prompt.md`.
     9	Blindness by absence: astra's builder report relocated during the review.
    10	
    11	## Leg outcomes
    12	- **Grok** (`_measurements/S11c_c2_export_repair_rereview_grok.txt`): all 6 checks — **CLEAR to commit**. Physics
    13	  unchanged (baseline vs current TRIAGE byte-identical; construction functions unchanged, only the export map
    14	  differs in `run()`), membership minimal (increment absent, both operators, closure minimal), semantic
    15	  equivalence clean, transparency clean (no Dummy/UnevaluatedExpr), **guards BITE** (its `out+1` compaction
    16	  corruption hard-stopped the semantic guard; membership re-add/drop both caught), no overstep.
    17	- **Fresh Claude agent** (incomplete — wedged twice on reaped background waiters; stopped by orchestrator): checks
    18	  2 (membership clean), 4 (transparency clean), 6 (no overstep) corroborate Grok. BUT it reported two APPARENT
    19	  defects — both **disproven** (see below): a semantic recompute FAIL and a "B2 corruption returned clean". Its
    20	  cross-process re-run method introduced fresh Integral bound-dummies, producing artifacts
    21	  (`_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out`: the nonzero leaves are
    22	  `2·_Dummy_34·X − 2·_Dummy_44·X` — identical term, relabeled dummy = the artifact signature).
    23	
    24	## The leg-vs-leg disagreement — SETTLED by my own same-process computation (rule 13)
    25	Command: `_measurements/S11c_c2_export_repair_adjudication_verify.py` → `…_verify.out`. Built one case in-process
    26	(no serialization/cross-process dummy mismatch), ran astra's `publication_compact` on every closed-slab-operator
    27	leaf, and applied astra's own algebraic-leaf comparison (`normalize_integrals` + Integral-protected expand):
    28	```
    29	closed slab operator: 5 algebraic leaves
    30	FAITHFUL (compact==emitted):      5/5  fails=[]
    31	GUARD BITES x2 corruption:        5/5  slipped=[]
    32	GUARD BITES +marker corruption:   5/5  slipped=[]
    33	INTEGRAND corruption on /THETA: bites=True
    34	```
    35	⇒ `publication_compact` is **faithful** (loses no information — the Agent's semantic FAIL is a cross-process
    36	dummy artifact, not a real difference), and the semantic guard **bites** every real value change — doubling, a
    37	distinct marker term, and an in-integrand change (the Agent's "B2 returned clean" was its own no-op corruption,
    38	not a guard gap). This independently corroborates Grok and refutes both of the Agent's apparent defects.
    39	
    40	## Verdict — the export repair is SOUND; correct to commit
    41	Faithful compaction (my 5/5 + Grok + astra's in-run guard), guards bite (my check + Grok), membership minimal +
    42	correct (increment EMIT-only; both closed operators, all 4 cases; exported==required closure), transparency (no
    43	Dummy/UnevaluatedExpr; `display` shrunk to the root name), singular locus preserved (astra's pole-revert +
    44	`reciprocal_powers_unchanged`), physics unchanged (no construction touched; TRIAGE byte-identical), no overstep.
    45	The 2-leg discipline is met in substance: Grok is a complete clean leg, my rule-13 verification independently
    46	re-confirms the load-bearing guard-bite + faithfulness, and the Agent's clean checks corroborate while its two
    47	divergences are conclusively artifacts of its cross-process method.
    48	
    49	⚠ **Process note (for leg-quality tracking):** the fresh-Claude leg repeatedly wedged on reaped background waiters
    50	and its cross-process semantic method produced false-positive FAILs; a same-process comparison (as Grok and I did)
    51	is the correct method for Integral-bearing operators. Its methodology, not a real defect, drove the divergence.
    52	
    53	**Next:** commit the repaired baseline; then the light §5e/§3c spec-wording clarification; then the blind Wolfram
    54	engine → its 2 legs → the c2 T7 comparator + reconcile → the c2 step record. ⛔ c1 STANDS.
     1	Resume S11c-c2 — physics SOUND + export repaired, both committed; NEXT = light spec-wording fix → blind Wolfram engine → T7 → step record. Trust `git log --oneline -8`. Live handoff: [[project_s11c_c_state]].
     2	
     3	STATE (all committed this session): c1 fully DONE (⛔ do NOT reopen). c2 fold physics 2-leg-reviewed SOUND `8f3a017f` (fresh Claude agent + Grok, identical prompt; Grok's F/G/E "defects" ALL FALSE POSITIVES — refuted by my own computation, rule 13; adjudication `_measurements/S11c_c2_physics_review_adjudication.md` + `…_adjudication_verify_{F,EG}.py/.out`). c2 EXPORT REPAIR committed `aa76105a` — **60 MB → 21.4 MB**, publication-only (physics untouched), directive gated `a5f7a06c` (2 decision legs), re-reviewed clear (`_measurements/S11c_c2_export_repair_rereview_adjudication.md`; my same-process verify: faithful 5/5, guards bite 5/5 — the fresh-Claude re-review leg's "FAIL" was a cross-process bound-dummy artifact). STATUS top clause `d2befb7c`.
     4	
     5	⭐ CARRY-FORWARDS owed to the c2 STEP RECORD (⛔ NOT build defects — the fold is SOUND): (F) the uniform-limit object does NOT literally vanish because the §3c increment retains the −extract(open) open slot — the record must say "the genuine closure-induced coupling **DECOUPLES**" (Integral integrand ≡0, verified), ⛔ NOT "the increment vanishes"; (E) N6 rep-invariance HOLDS at leading O(ε)/O(εη) (σ_W→0 kills the whole residual), the σ_W-sector remnant is DEFERRED (drain-projection ≥64 GB family); (G) the induced self-energy is **DIRECTIONAL/one-way** (reverse THICKNESS→TRANSVERSE block ≡0 — δp lives only in the θ/mech rows; a dissipative directional self-energy, expected + both blocks emitted); the two S11c-b sign conventions (face-generalized-force, #90 closure-fold) multiply the increment and are SURFACED by the WL cross-engine comparator (⛔ not normalized; §3d.4 mechanical-power pairing adjudicates the face-force sign).
     6	
     7	NEXT:
     8	  STEP 1 — light SPEC-WORDING clarification to `directives/S11c_c2_SHARED_PHYSICS.md`: §5e ("the off-diagonal self-energy increment must vanish" in the uniform limit → "the genuine closure-induced coupling decouples; the emitted increment retains the §3c open-slot representation and is not literally zero") + §3c (note the increment carries the −extract(open) O(ε) open-slot piece by construction) + reconcile "export representation" → "comparison/emit representation" (minor). ⛔ WORDING ONLY — physics unchanged (both legs + my compute confirmed SOUND). The spec is **review-until-clear**: orchestrator-written → 2 legs (Codex `gpt-5.6-sol` xhigh + Grok), fold, iterate to clear. Then commit.
     9	  STEP 2 — the blind WOLFRAM engine (imports NOTHING, re-derives §§1–2 + the S11c-a face substrate + the S11c-b slab rows it folds into + c1's closed response — `S9_export_chain_rebuild_directive.md:16-18` is the ONLY cross-engine control): orchestrator-written WL BUILD DIRECTIVE → **2 decision legs** (Codex sol xhigh + Grok) BEFORE the build → detached Mathematica build (2-seat licence, `--sandbox danger-full-access`, setsid+marker) → **2 build legs** (fresh Claude + Grok, SERIALIZED — both ablate kernels, 2 seats). Model = astra (Mathematica). Reuse the c1 WL engine shape (`project_s11c_c_state` has the c1 WL arc).
    10	  STEP 3 — the c2 T7 COMPARATOR + reconcile (frozen T7 contract: join the two engines' emitted objects by name, three-valued, PRINTS/decides nothing; the load-bearing residual is on the **self-energy increment** — which carries the 2 sign conventions [surface, don't normalize] and the §3d representation questions; the reconcile = the staged representational bridge [[feedback_reconcile_representational_bridge]], ⛔ never a blanket collapse; full per-family symbolic residual DEFERRED ≥64 GB, name don't pre-adjudicate) → the **c2 STEP RECORD** (fold the carry-forwards above; 2 record legs; NO per-substep card — one S11c roll-up after S11c-e).
    11	
    12	⭐ MODEL POLICY: CODE→`gpt-6-astra` high; DOC/decision review→`gpt-5.6-sol` xhigh; the WL builder = astra. astra assessment (2026-09-06, user asked): KEEP astra for CODE — this session's export-repair code was correct + defect-free under adversarial review + NO overstep (the fixed §7 fence held); only weakness = speed (it ran 3 full builds). The earlier c1 overstep was MY directive-fence defect (fixed), not astra.
    13	⚠ OPS (measured this session): ALL `run_in_background` tasks were REAPED (~30 s, even trivial `until`-loop watchers; NOT OOM — 23 GB free). Use **setsid-detached launch + a DONE-marker + a re-armed watcher** for any heavy build; the watcher's death IS a poll tick — re-check the marker + re-arm. Integral-bearing CAS objects need **SAME-PROCESS** comparison (a cross-process re-run mismatches Integral bound dummies → spurious FAIL). [[feedback_background_tasks_can_die_spuriously]]. ⛔ never 2 memory-heavy CAS jobs concurrently. ⚠ MEMORY.md ~20 KB — a full compaction (drop/merge stale entries to <17 KB, 1 line/entry, ⛔ don't drop a live pointer) is OWED; do it carefully early this session.
    14	⭐⭐ GOVERNING: a DISAGREEMENT is a finding — verify it MYSELF (rule 13; both legs have been wrong this session, in both directions); scripts PRINT never assert; a spec = review-until-clear; the WL engine is BLIND (imports nothing); ⛔ NEVER fabricate. ⛔ c1 STANDS.

exec
/bin/bash -lc "git diff --stat 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && git diff --unified=12 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && rg -n '"'^def |EXPORT_ROOTS|export_key|publication_compact|publish'"\\(' research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && git diff --name-only 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && git status --short" in /var/projects/toy_physics
 succeeded in 0ms:
 .../scripts/S11c_c2_selfenergy_fold_sympy_audit.py | 168 ++++++++++++++++++++-
 1 file changed, 162 insertions(+), 6 deletions(-)
diff --git a/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py b/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
index 68a3bfd2..43218c15 100644
--- a/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
+++ b/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
@@ -38,25 +38,24 @@ IMPORT_KEYS = (
     'dtn_operator', 'dtn_flat_symbol', 'dtn_kernel', 'coupling_kernel',
     'rho_m', 'rho_br', 'W_0', 'epsilon_shape', 'eta_bg', 'sigma_W',
     'W_bg', 'L_W', 'mu_R', 'omega', 'c_s0',
     'Lambda_A_0', 'Lambda_V_0', 'Lambda_X_0',
     'rho_br_bg_rho4_constant',
 )
 ANCHORINGS = ('LAB_HELD', 'MATERIAL_ADVECTED')
 DENSITIES = ('RHO4_CONSTANT', 'RHOBR_CONSTANT')
 FACES = (1, -1)
 REPRESENTATION = 'DELTA_W'
 EXPORT_ROOTS = frozenset((
     's11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel',
-    's11cc2SelfEnergyIncrement',
 ))
 NEW_DIMENSIONS = {}
 COMPUTATION_LINES = {}
 EMITTED_KEYS = set()
 OUTPUT_ROWS = {}
 COMPUTED_BINDINGS = {}
 
 # Dimensional declaration schema of the supplied jet ansatz.  Populated below
 # with the inherited S11c-b/c1 declaration dimensions (no result values).
 DIMENSION_SCHEMA = {'B_rho_3': [-1, -2, 1], 'C': [-2, -2, 1], 'G_W_u': [-1, -2, 1], 'G_theta_u': [-1, -2, 1], 'L_W': [1, 0, 0], 'Lambda_A_0': [-5, 1, 1], 'Lambda_V_0': [-4, 0, 1], 'Lambda_X_0': [-4, 0, 1], 'W_0': [1, 0, 0], 'W_bg': [1, 0, 0], 'W_bg_d1': [0, 0, 0], 'W_bg_d2': [0, 0, 0], 'W_bg_d3': [0, 0, 0], 'c_s0': [1, -1, 0], 'd_w_delta_j_bulk_minus_1': [-4, -1, 1], 'd_w_delta_j_bulk_minus_2': [-4, -1, 1], 'd_w_delta_j_bulk_minus_3': [-4, -1, 1], 'd_w_delta_j_bulk_minus_4': [-4, -1, 1], 'd_w_delta_j_bulk_plus_1': [-4, -1, 1], 'd_w_delta_j_bulk_plus_2': [-4, -1, 1], 'd_w_delta_j_bulk_plus_3': [-4, -1, 1], 'd_w_delta_j_bulk_plus_4': [-4, -1, 1], 'd_w_delta_p_minus': [-3, -2, 1], 'd_w_delta_p_plus': [-3, -2, 1], 'd_w_delta_rho_4D_face_minus': [-5, 0, 1], 'd_w_delta_rho_4D_face_plus': [-5, 0, 1], 'd_w_delta_v_bulk_minus_1': [0, -1, 0], 'd_w_delta_v_bulk_minus_2': [0, -1, 0], 'd_w_delta_v_bulk_minus_3': [0, -1, 0], 'd_w_delta_v_bulk_minus_4': [0, -1, 0], 'd_w_delta_v_bulk_plus_1': [0, -1, 0], 'd_w_delta_v_bulk_plus_2': [0, -1, 0], 'd_w_delta_v_bulk_plus_3': [0, -1, 0], 'd_w_delta_v_bulk_plus_4': [0, -1, 0], 'd_w_trace_grad_f_1': [-2, 0, 0], 'd_w_trace_grad_f_2': [-2, 0, 0], 'd_w_trace_grad_f_3': [-2, 0, 0], 'd_w_trace_grad_f_4': [-2, 0, 0], 'delta_j_bulk_1': [-3, -1, 1], 'delta_j_bulk_2': [-3, -1, 1], 'delta_j_bulk_3': [-3, -1, 1], 'delta_j_bulk_4': [-3, -1, 1], 'delta_p_minus': [-2, -2, 1], 'delta_p_plus': [-2, -2, 1], 'delta_rho_4D_bulk_t': [-4, -1, 1], 'delta_rho_4D_face_minus': [-4, 0, 1], 'delta_rho_4D_face_plus': [-4, 0, 1], 'delta_v_bulk_minus_1': [1, -1, 0], 'delta_v_bulk_minus_2': [1, -1, 0], 'delta_v_bulk_minus_3': [1, -1, 0], 'delta_v_bulk_minus_4': [1, -1, 0], 'delta_v_bulk_plus_1': [1, -1, 0], 'delta_v_bulk_plus_2': [1, -1, 0], 'delta_v_bulk_plus_3': [1, -1, 0], 'delta_v_bulk_plus_4': [1, -1, 0], 'delta_v_e_W': [0, 0, 0], 'delta_v_theta': [0, 0, 0], 'delta_v_u_1': [1, 0, 0], 'delta_v_u_1_d1': [0, 0, 0], 'delta_v_u_1_d2': [0, 0, 0], 'delta_v_u_1_d3': [0, 0, 0], 'delta_v_u_2': [1, 0, 0], 'delta_v_u_2_d1': [0, 0, 0], 'delta_v_u_2_d2': [0, 0, 0], 'delta_v_u_2_d3': [0, 0, 0], 'delta_v_u_3': [1, 0, 0], 'delta_v_u_3_d1': [0, 0, 0], 'delta_v_u_3_d2': [0, 0, 0], 'delta_v_u_3_d3': [0, 0, 0], 'delta_v_zeta_c': [1, 0, 0], 'e_W': [0, 0, 0], 'e_W_bg': [0, 0, 0], 'e_W_d1': [-1, 0, 0], 'e_W_d1d1': [-2, 0, 0], 'e_W_d1d1d1': [-3, 0, 0], 'e_W_d1d1d1d1': [-4, 0, 0], 'e_W_d1d1d1d2': [-4, 0, 0], 'e_W_d1d1d1d3': [-4, 0, 0], 'e_W_d1d1d2': [-3, 0, 0], 'e_W_d1d1d2d2': [-4, 0, 0], 'e_W_d1d1d2d3': [-4, 0, 0], 'e_W_d1d1d3': [-3, 0, 0], 'e_W_d1d1d3d3': [-4, 0, 0], 'e_W_d1d2': [-2, 0, 0], 'e_W_d1d2d2': [-3, 0, 0], 'e_W_d1d2d2d2': [-4, 0, 0], 'e_W_d1d2d2d3': [-4, 0, 0], 'e_W_d1d2d3': [-3, 0, 0], 'e_W_d1d2d3d3': [-4, 0, 0], 'e_W_d1d3': [-2, 0, 0], 'e_W_d1d3d3': [-3, 0, 0], 'e_W_d1d3d3d3': [-4, 0, 0], 'e_W_d2': [-1, 0, 0], 'e_W_d2d2': [-2, 0, 0], 'e_W_d2d2d2': [-3, 0, 0], 'e_W_d2d2d2d2': [-4, 0, 0], 'e_W_d2d2d2d3': [-4, 0, 0], 'e_W_d2d2d3': [-3, 0, 0], 'e_W_d2d2d3d3': [-4, 0, 0], 'e_W_d2d3': [-2, 0, 0], 'e_W_d2d3d3': [-3, 0, 0], 'e_W_d2d3d3d3': [-4, 0, 0], 'e_W_d3': [-1, 0, 0], 'e_W_d3d3': [-2, 0, 0], 'e_W_d3d3d3': [-3, 0, 0], 'e_W_d3d3d3d3': [-4, 0, 0], 'e_W_probe': [0, 0, 0], 'e_W_probe_d1': [-1, 0, 0], 'e_W_probe_d1d1': [-2, 0, 0], 'e_W_probe_d1d1d1': [-3, 0, 0], 'e_W_probe_d1d1d1d1': [-4, 0, 0], 'e_W_probe_d1d1d1d2': [-4, 0, 0], 'e_W_probe_d1d1d1d3': [-4, 0, 0], 'e_W_probe_d1d1d2': [-3, 0, 0], 'e_W_probe_d1d1d2d2': [-4, 0, 0], 'e_W_probe_d1d1d2d3': [-4, 0, 0], 'e_W_probe_d1d1d3': [-3, 0, 0], 'e_W_probe_d1d1d3d3': [-4, 0, 0], 'e_W_probe_d1d2': [-2, 0, 0], 'e_W_probe_d1d2d2': [-3, 0, 0], 'e_W_probe_d1d2d2d2': [-4, 0, 0], 'e_W_probe_d1d2d2d3': [-4, 0, 0], 'e_W_probe_d1d2d3': [-3, 0, 0], 'e_W_probe_d1d2d3d3': [-4, 0, 0], 'e_W_probe_d1d3': [-2, 0, 0], 'e_W_probe_d1d3d3': [-3, 0, 0], 'e_W_probe_d1d3d3d3': [-4, 0, 0], 'e_W_probe_d2': [-1, 0, 0], 'e_W_probe_d2d2': [-2, 0, 0], 'e_W_probe_d2d2d2': [-3, 0, 0], 'e_W_probe_d2d2d2d2': [-4, 0, 0], 'e_W_probe_d2d2d2d3': [-4, 0, 0], 'e_W_probe_d2d2d3': [-3, 0, 0], 'e_W_probe_d2d2d3d3': [-4, 0, 0], 'e_W_probe_d2d3': [-2, 0, 0], 'e_W_probe_d2d3d3': [-3, 0, 0], 'e_W_probe_d2d3d3d3': [-4, 0, 0], 'e_W_probe_d3': [-1, 0, 0], 'e_W_probe_d3d3': [-2, 0, 0], 'e_W_probe_d3d3d3': [-3, 0, 0], 'e_W_probe_d3d3d3d3': [-4, 0, 0], 'e_W_probe_t': [0, -1, 0], 'e_W_t': [0, -1, 0], 'e_W_tt': [0, -2, 0], 'epsilon_shape': [0, 0, 0], 'eta_bg': [0, 0, 0], 'f_hold_e_W_0': [-1, -2, 1], 'f_hold_theta_0': [-1, -2, 1], 'f_hold_u_1_0': [-2, -2, 1], 'f_hold_u_2_0': [-2, -2, 1], 'f_hold_u_3_0': [-2, -2, 1], 'grad_theta_1': [-1, 0, 0], 'grad_theta_2': [-1, 0, 0], 'grad_theta_3': [-1, 0, 0], 'k': [-1, 0, 0], 'k_W': [-3, -2, 1], 'kappa_W': [-3, -2, 1], 'kappa_theta': [1, -2, 1], 'kappa_theta_W': [1, -2, 1], 'm1_profile': [0, 0, 0], 'm1_profile_d1': [0, 0, 0], 'm1_profile_d1d1': [0, 0, 0], 'm1_profile_d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1d1': [0, 0, 0], 'm1_profile_d1d1d1d1d2': [0, 0, 0], 'm1_profile_d1d1d1d1d3': [0, 0, 0], 'm1_profile_d1d1d1d2': [0, 0, 0], 'm1_profile_d1d1d1d2d2': [0, 0, 0], 'm1_profile_d1d1d1d2d3': [0, 0, 0], 'm1_profile_d1d1d1d3': [0, 0, 0], 'm1_profile_d1d1d1d3d3': [0, 0, 0], 'm1_profile_d1d1d2': [0, 0, 0], 'm1_profile_d1d1d2d2': [0, 0, 0], 'm1_profile_d1d1d2d2d2': [0, 0, 0], 'm1_profile_d1d1d2d2d3': [0, 0, 0], 'm1_profile_d1d1d2d3': [0, 0, 0], 'm1_profile_d1d1d2d3d3': [0, 0, 0], 'm1_profile_d1d1d3': [0, 0, 0], 'm1_profile_d1d1d3d3': [0, 0, 0], 'm1_profile_d1d1d3d3d3': [0, 0, 0], 'm1_profile_d1d2': [0, 0, 0], 'm1_profile_d1d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2d2': [0, 0, 0], 'm1_profile_d1d2d2d2d3': [0, 0, 0], 'm1_profile_d1d2d2d3': [0, 0, 0], 'm1_profile_d1d2d2d3d3': [0, 0, 0], 'm1_profile_d1d2d3': [0, 0, 0], 'm1_profile_d1d2d3d3': [0, 0, 0], 'm1_profile_d1d2d3d3d3': [0, 0, 0], 'm1_profile_d1d3': [0, 0, 0], 'm1_profile_d1d3d3': [0, 0, 0], 'm1_profile_d1d3d3d3': [0, 0, 0], 'm1_profile_d1d3d3d3d3': [0, 0, 0], 'm1_profile_d2': [0, 0, 0], 'm1_profile_d2d2': [0, 0, 0], 'm1_profile_d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2d2': [0, 0, 0], 'm1_profile_d2d2d2d2d3': [0, 0, 0], 'm1_profile_d2d2d2d3': [0, 0, 0], 'm1_profile_d2d2d2d3d3': [0, 0, 0], 'm1_profile_d2d2d3': [0, 0, 0], 'm1_profile_d2d2d3d3': [0, 0, 0], 'm1_profile_d2d2d3d3d3': [0, 0, 0], 'm1_profile_d2d3': [0, 0, 0], 'm1_profile_d2d3d3': [0, 0, 0], 'm1_profile_d2d3d3d3': [0, 0, 0], 'm1_profile_d2d3d3d3d3': [0, 0, 0], 'm1_profile_d3': [0, 0, 0], 'm1_profile_d3d3': [0, 0, 0], 'm1_profile_d3d3d3': [0, 0, 0], 'm1_profile_d3d3d3d3': [0, 0, 0], 'm1_profile_d3d3d3d3d3': [0, 0, 0], 'mu_R': [-1, -2, 1], 'mu_R_bg': [-1, -2, 1], 'mu_R_bg_d1': [-2, -2, 1], 'mu_R_bg_d2': [-2, -2, 1], 'mu_R_bg_d3': [-2, -2, 1], 'mu_S': [-1, -2, 1], 'mu_W': [-3, 0, 1], 'mu_theta_L': [-1, -2, 1], 'mu_theta_M': [-1, -2, 1], 'omega': [0, -1, 0], 'phi_L_d1d1': [0, 0, 0], 'phi_L_d1d1d1': [-1, 0, 0], 'phi_L_d1d1d1d1': [-2, 0, 0], 'phi_L_d1d1d1d1d1': [-3, 0, 0], 'phi_L_d1d1d1d1d2': [-3, 0, 0], 'phi_L_d1d1d1d1d3': [-3, 0, 0], 'phi_L_d1d1d1d2': [-2, 0, 0], 'phi_L_d1d1d1d2d2': [-3, 0, 0], 'phi_L_d1d1d1d2d3': [-3, 0, 0], 'phi_L_d1d1d1d3': [-2, 0, 0], 'phi_L_d1d1d1d3d3': [-3, 0, 0], 'phi_L_d1d1d2': [-1, 0, 0], 'phi_L_d1d1d2d2': [-2, 0, 0], 'phi_L_d1d1d2d2d2': [-3, 0, 0], 'phi_L_d1d1d2d2d3': [-3, 0, 0], 'phi_L_d1d1d2d3': [-2, 0, 0], 'phi_L_d1d1d2d3d3': [-3, 0, 0], 'phi_L_d1d1d3': [-1, 0, 0], 'phi_L_d1d1d3d3': [-2, 0, 0], 'phi_L_d1d1d3d3d3': [-3, 0, 0], 'phi_L_d1d2': [0, 0, 0], 'phi_L_d1d2d2': [-1, 0, 0], 'phi_L_d1d2d2d2': [-2, 0, 0], 'phi_L_d1d2d2d2d2': [-3, 0, 0], 'phi_L_d1d2d2d2d3': [-3, 0, 0], 'phi_L_d1d2d2d3': [-2, 0, 0], 'phi_L_d1d2d2d3d3': [-3, 0, 0], 'phi_L_d1d2d3': [-1, 0, 0], 'phi_L_d1d2d3d3': [-2, 0, 0], 'phi_L_d1d2d3d3d3': [-3, 0, 0], 'phi_L_d1d3': [0, 0, 0], 'phi_L_d1d3d3': [-1, 0, 0], 'phi_L_d1d3d3d3': [-2, 0, 0], 'phi_L_d1d3d3d3d3': [-3, 0, 0], 'phi_L_d2d2': [0, 0, 0], 'phi_L_d2d2d2': [-1, 0, 0], 'phi_L_d2d2d2d2': [-2, 0, 0], 'phi_L_d2d2d2d2d2': [-3, 0, 0], 'phi_L_d2d2d2d2d3': [-3, 0, 0], 'phi_L_d2d2d2d3': [-2, 0, 0], 'phi_L_d2d2d2d3d3': [-3, 0, 0], 'phi_L_d2d2d3': [-1, 0, 0], 'phi_L_d2d2d3d3': [-2, 0, 0], 'phi_L_d2d2d3d3d3': [-3, 0, 0], 'phi_L_d2d3': [0, 0, 0], 'phi_L_d2d3d3': [-1, 0, 0], 'phi_L_d2d3d3d3': [-2, 0, 0], 'phi_L_d2d3d3d3d3': [-3, 0, 0], 'phi_L_d3d3': [0, 0, 0], 'phi_L_d3d3d3': [-1, 0, 0], 'phi_L_d3d3d3d3': [-2, 0, 0], 'phi_L_d3d3d3d3d3': [-3, 0, 0], 'rho_4D_bg_rho4_constant': [-4, 0, 1], 'rho_4D_bg_rhobr_constant': [-4, 0, 1], 'rho_br': [-3, 0, 1], 'rho_br_bg_rho4_constant': [-3, 0, 1], 'rho_br_bg_rhobr_constant': [-3, 0, 1], 'rho_m': [-4, 0, 1], 's11cc1_V_lab_held_minus': [1, -1, 0], 's11cc1_V_lab_held_plus': [1, -1, 0], 's11cc1_V_material_advected_minus': [1, -1, 0], 's11cc1_V_material_advected_plus': [1, -1, 0], 's11cc1_div_height_grad_operator': [0, 0, 0], 's11cc1_dtn_operator_lab_held_minus': [0, 0, 0], 's11cc1_dtn_operator_lab_held_plus': [0, 0, 0], 's11cc1_dtn_operator_material_advected_minus': [0, 0, 0], 's11cc1_dtn_operator_material_advected_plus': [0, 0, 0], 's11cc1_first_shape_impedance_operator': [0, 0, 0], 's11cc1_flat_impedance_operator': [0, 0, 0], 's11cc1_flat_normal_dtn_inverse': [0, 0, 0], 's11cc1_flat_normal_dtn_operator': [0, 0, 0], 's11cc1_height_multiplication_operator': [0, 0, 0], 's11cc1_identity_operator': [0, 0, 0], 's11cc1_k_input_1': [-1, 0, 0], 's11cc1_k_input_2': [-1, 0, 0], 's11cc1_k_input_3': [-1, 0, 0], 's11cc1_k_output_1': [-1, 0, 0], 's11cc1_k_output_2': [-1, 0, 0], 's11cc1_k_output_3': [-1, 0, 0], 's11cc1_mu_theta_lab_held_minus': [-1, -2, 1], 's11cc1_mu_theta_lab_held_plus': [-1, -2, 1], 's11cc1_mu_theta_material_advected_minus': [-1, -2, 1], 's11cc1_mu_theta_material_advected_plus': [-1, -2, 1], 's11cc1_q_out_input': [-1, 0, 0], 's11cc1_q_out_output': [-1, 0, 0], 's11cc1_response_resolvent_lab_held_minus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_minus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_plus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_lab_held_plus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_minus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_minus_rhobr_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_plus_rho4_constant': [0, 0, 0], 's11cc1_response_resolvent_material_advected_plus_rhobr_constant': [0, 0, 0], 's11cc1_w1_profile_hat_transfer': [3, 0, 0], 's11cc1_w1_profile_jet_hat_1': [3, 0, 0], 's11cc1_w1_profile_jet_hat_2': [3, 0, 0], 's11cc1_w1_profile_jet_hat_3': [3, 0, 0], 'sigma_W': [0, 0, 0], 't_hold_minus_0_1': [-2, -2, 1], 't_hold_minus_0_2': [-2, -2, 1], 't_hold_minus_0_3': [-2, -2, 1], 't_hold_minus_0_4': [-2, -2, 1], 't_hold_plus_0_1': [-2, -2, 1], 't_hold_plus_0_2': [-2, -2, 1], 't_hold_plus_0_3': [-2, -2, 1], 't_hold_plus_0_4': [-2, -2, 1], 'tau_A': [0, 1, 0], 'tau_V': [0, 1, 0], 'tau_X': [0, 1, 0], 'theta': [0, 0, 0], 'theta_d1d1': [-2, 0, 0], 'theta_d1d1d1': [-3, 0, 0], 'theta_d1d1d1d1': [-4, 0, 0], 'theta_d1d1d1d2': [-4, 0, 0], 'theta_d1d1d1d3': [-4, 0, 0], 'theta_d1d1d2': [-3, 0, 0], 'theta_d1d1d2d2': [-4, 0, 0], 'theta_d1d1d2d3': [-4, 0, 0], 'theta_d1d1d3': [-3, 0, 0], 'theta_d1d1d3d3': [-4, 0, 0], 'theta_d1d2': [-2, 0, 0], 'theta_d1d2d2': [-3, 0, 0], 'theta_d1d2d2d2': [-4, 0, 0], 'theta_d1d2d2d3': [-4, 0, 0], 'theta_d1d2d3': [-3, 0, 0], 'theta_d1d2d3d3': [-4, 0, 0], 'theta_d1d3': [-2, 0, 0], 'theta_d1d3d3': [-3, 0, 0], 'theta_d1d3d3d3': [-4, 0, 0], 'theta_d2d2': [-2, 0, 0], 'theta_d2d2d2': [-3, 0, 0], 'theta_d2d2d2d2': [-4, 0, 0], 'theta_d2d2d2d3': [-4, 0, 0], 'theta_d2d2d3': [-3, 0, 0], 'theta_d2d2d3d3': [-4, 0, 0], 'theta_d2d3': [-2, 0, 0], 'theta_d2d3d3': [-3, 0, 0], 'theta_d2d3d3d3': [-4, 0, 0], 'theta_d3d3': [-2, 0, 0], 'theta_d3d3d3': [-3, 0, 0], 'theta_d3d3d3d3': [-4, 0, 0], 'theta_probe': [0, 0, 0], 'theta_probe_d1': [-1, 0, 0], 'theta_probe_d1d1': [-2, 0, 0], 'theta_probe_d1d1d1': [-3, 0, 0], 'theta_probe_d1d1d1d1': [-4, 0, 0], 'theta_probe_d1d1d1d2': [-4, 0, 0], 'theta_probe_d1d1d1d3': [-4, 0, 0], 'theta_probe_d1d1d2': [-3, 0, 0], 'theta_probe_d1d1d2d2': [-4, 0, 0], 'theta_probe_d1d1d2d3': [-4, 0, 0], 'theta_probe_d1d1d3': [-3, 0, 0], 'theta_probe_d1d1d3d3': [-4, 0, 0], 'theta_probe_d1d2': [-2, 0, 0], 'theta_probe_d1d2d2': [-3, 0, 0], 'theta_probe_d1d2d2d2': [-4, 0, 0], 'theta_probe_d1d2d2d3': [-4, 0, 0], 'theta_probe_d1d2d3': [-3, 0, 0], 'theta_probe_d1d2d3d3': [-4, 0, 0], 'theta_probe_d1d3': [-2, 0, 0], 'theta_probe_d1d3d3': [-3, 0, 0], 'theta_probe_d1d3d3d3': [-4, 0, 0], 'theta_probe_d2': [-1, 0, 0], 'theta_probe_d2d2': [-2, 0, 0], 'theta_probe_d2d2d2': [-3, 0, 0], 'theta_probe_d2d2d2d2': [-4, 0, 0], 'theta_probe_d2d2d2d3': [-4, 0, 0], 'theta_probe_d2d2d3': [-3, 0, 0], 'theta_probe_d2d2d3d3': [-4, 0, 0], 'theta_probe_d2d3': [-2, 0, 0], 'theta_probe_d2d3d3': [-3, 0, 0], 'theta_probe_d2d3d3d3': [-4, 0, 0], 'theta_probe_d3': [-1, 0, 0], 'theta_probe_d3d3': [-2, 0, 0], 'theta_probe_d3d3d3': [-3, 0, 0], 'theta_probe_d3d3d3d3': [-4, 0, 0], 'theta_probe_t': [0, -1, 0], 'theta_t': [0, -1, 0], 'trace_grad_f_1': [-1, 0, 0], 'trace_grad_f_2': [-1, 0, 0], 'trace_grad_f_3': [-1, 0, 0], 'trace_grad_f_4': [-1, 0, 0], 'u_1': [1, 0, 0], 'u_1_d1': [0, 0, 0], 'u_1_d1d1': [-1, 0, 0], 'u_1_d1d1d1': [-2, 0, 0], 'u_1_d1d1d1d1': [-3, 0, 0], 'u_1_d1d1d1d2': [-3, 0, 0], 'u_1_d1d1d1d3': [-3, 0, 0], 'u_1_d1d1d2': [-2, 0, 0], 'u_1_d1d1d2d2': [-3, 0, 0], 'u_1_d1d1d2d3': [-3, 0, 0], 'u_1_d1d1d3': [-2, 0, 0], 'u_1_d1d1d3d3': [-3, 0, 0], 'u_1_d1d2': [-1, 0, 0], 'u_1_d1d2d2': [-2, 0, 0], 'u_1_d1d2d2d2': [-3, 0, 0], 'u_1_d1d2d2d3': [-3, 0, 0], 'u_1_d1d2d3': [-2, 0, 0], 'u_1_d1d2d3d3': [-3, 0, 0], 'u_1_d1d3': [-1, 0, 0], 'u_1_d1d3d3': [-2, 0, 0], 'u_1_d1d3d3d3': [-3, 0, 0], 'u_1_d2': [0, 0, 0], 'u_1_d2d2': [-1, 0, 0], 'u_1_d2d2d2': [-2, 0, 0], 'u_1_d2d2d2d2': [-3, 0, 0], 'u_1_d2d2d2d3': [-3, 0, 0], 'u_1_d2d2d3': [-2, 0, 0], 'u_1_d2d2d3d3': [-3, 0, 0], 'u_1_d2d3': [-1, 0, 0], 'u_1_d2d3d3': [-2, 0, 0], 'u_1_d2d3d3d3': [-3, 0, 0], 'u_1_d3': [0, 0, 0], 'u_1_d3d3': [-1, 0, 0], 'u_1_d3d3d3': [-2, 0, 0], 'u_1_d3d3d3d3': [-3, 0, 0], 'u_1_t': [1, -1, 0], 'u_1_t_d1': [0, -1, 0], 'u_1_t_d2': [0, -1, 0], 'u_1_t_d3': [0, -1, 0], 'u_1_tt': [1, -2, 0], 'u_2': [1, 0, 0], 'u_2_d1': [0, 0, 0], 'u_2_d1d1': [-1, 0, 0], 'u_2_d1d1d1': [-2, 0, 0], 'u_2_d1d1d1d1': [-3, 0, 0], 'u_2_d1d1d1d2': [-3, 0, 0], 'u_2_d1d1d1d3': [-3, 0, 0], 'u_2_d1d1d2': [-2, 0, 0], 'u_2_d1d1d2d2': [-3, 0, 0], 'u_2_d1d1d2d3': [-3, 0, 0], 'u_2_d1d1d3': [-2, 0, 0], 'u_2_d1d1d3d3': [-3, 0, 0], 'u_2_d1d2': [-1, 0, 0], 'u_2_d1d2d2': [-2, 0, 0], 'u_2_d1d2d2d2': [-3, 0, 0], 'u_2_d1d2d2d3': [-3, 0, 0], 'u_2_d1d2d3': [-2, 0, 0], 'u_2_d1d2d3d3': [-3, 0, 0], 'u_2_d1d3': [-1, 0, 0], 'u_2_d1d3d3': [-2, 0, 0], 'u_2_d1d3d3d3': [-3, 0, 0], 'u_2_d2': [0, 0, 0], 'u_2_d2d2': [-1, 0, 0], 'u_2_d2d2d2': [-2, 0, 0], 'u_2_d2d2d2d2': [-3, 0, 0], 'u_2_d2d2d2d3': [-3, 0, 0], 'u_2_d2d2d3': [-2, 0, 0], 'u_2_d2d2d3d3': [-3, 0, 0], 'u_2_d2d3': [-1, 0, 0], 'u_2_d2d3d3': [-2, 0, 0], 'u_2_d2d3d3d3': [-3, 0, 0], 'u_2_d3': [0, 0, 0], 'u_2_d3d3': [-1, 0, 0], 'u_2_d3d3d3': [-2, 0, 0], 'u_2_d3d3d3d3': [-3, 0, 0], 'u_2_t': [1, -1, 0], 'u_2_t_d1': [0, -1, 0], 'u_2_t_d2': [0, -1, 0], 'u_2_t_d3': [0, -1, 0], 'u_2_tt': [1, -2, 0], 'u_3': [1, 0, 0], 'u_3_d1': [0, 0, 0], 'u_3_d1d1': [-1, 0, 0], 'u_3_d1d1d1': [-2, 0, 0], 'u_3_d1d1d1d1': [-3, 0, 0], 'u_3_d1d1d1d2': [-3, 0, 0], 'u_3_d1d1d1d3': [-3, 0, 0], 'u_3_d1d1d2': [-2, 0, 0], 'u_3_d1d1d2d2': [-3, 0, 0], 'u_3_d1d1d2d3': [-3, 0, 0], 'u_3_d1d1d3': [-2, 0, 0], 'u_3_d1d1d3d3': [-3, 0, 0], 'u_3_d1d2': [-1, 0, 0], 'u_3_d1d2d2': [-2, 0, 0], 'u_3_d1d2d2d2': [-3, 0, 0], 'u_3_d1d2d2d3': [-3, 0, 0], 'u_3_d1d2d3': [-2, 0, 0], 'u_3_d1d2d3d3': [-3, 0, 0], 'u_3_d1d3': [-1, 0, 0], 'u_3_d1d3d3': [-2, 0, 0], 'u_3_d1d3d3d3': [-3, 0, 0], 'u_3_d2': [0, 0, 0], 'u_3_d2d2': [-1, 0, 0], 'u_3_d2d2d2': [-2, 0, 0], 'u_3_d2d2d2d2': [-3, 0, 0], 'u_3_d2d2d2d3': [-3, 0, 0], 'u_3_d2d2d3': [-2, 0, 0], 'u_3_d2d2d3d3': [-3, 0, 0], 'u_3_d2d3': [-1, 0, 0], 'u_3_d2d3d3': [-2, 0, 0], 'u_3_d2d3d3d3': [-3, 0, 0], 'u_3_d3': [0, 0, 0], 'u_3_d3d3': [-1, 0, 0], 'u_3_d3d3d3': [-2, 0, 0], 'u_3_d3d3d3d3': [-3, 0, 0], 'u_3_t': [1, -1, 0], 'u_3_t_d1': [0, -1, 0], 'u_3_t_d2': [0, -1, 0], 'u_3_t_d3': [0, -1, 0], 'u_3_tt': [1, -2, 0], 'u_L_1': [1, 0, 0], 'u_L_1_d1': [0, 0, 0], 'u_L_1_d1d1': [-1, 0, 0], 'u_L_1_d1d1d1': [-2, 0, 0], 'u_L_1_d1d1d1d1': [-3, 0, 0], 'u_L_1_d1d1d1d2': [-3, 0, 0], 'u_L_1_d1d1d1d3': [-3, 0, 0], 'u_L_1_d1d1d2': [-2, 0, 0], 'u_L_1_d1d1d2d2': [-3, 0, 0], 'u_L_1_d1d1d2d3': [-3, 0, 0], 'u_L_1_d1d1d3': [-2, 0, 0], 'u_L_1_d1d1d3d3': [-3, 0, 0], 'u_L_1_d1d2': [-1, 0, 0], 'u_L_1_d1d2d2': [-2, 0, 0], 'u_L_1_d1d2d2d2': [-3, 0, 0], 'u_L_1_d1d2d2d3': [-3, 0, 0], 'u_L_1_d1d2d3': [-2, 0, 0], 'u_L_1_d1d2d3d3': [-3, 0, 0], 'u_L_1_d1d3': [-1, 0, 0], 'u_L_1_d1d3d3': [-2, 0, 0], 'u_L_1_d1d3d3d3': [-3, 0, 0], 'u_L_1_d2': [0, 0, 0], 'u_L_1_d2d2': [-1, 0, 0], 'u_L_1_d2d2d2': [-2, 0, 0], 'u_L_1_d2d2d2d2': [-3, 0, 0], 'u_L_1_d2d2d2d3': [-3, 0, 0], 'u_L_1_d2d2d3': [-2, 0, 0], 'u_L_1_d2d2d3d3': [-3, 0, 0], 'u_L_1_d2d3': [-1, 0, 0], 'u_L_1_d2d3d3': [-2, 0, 0], 'u_L_1_d2d3d3d3': [-3, 0, 0], 'u_L_1_d3': [0, 0, 0], 'u_L_1_d3d3': [-1, 0, 0], 'u_L_1_d3d3d3': [-2, 0, 0], 'u_L_1_d3d3d3d3': [-3, 0, 0], 'u_L_1_t': [1, -1, 0], 'u_L_1_t_d1': [0, -1, 0], 'u_L_1_t_d2': [0, -1, 0], 'u_L_1_t_d3': [0, -1, 0], 'u_L_2': [1, 0, 0], 'u_L_2_d1': [0, 0, 0], 'u_L_2_d1d1': [-1, 0, 0], 'u_L_2_d1d1d1': [-2, 0, 0], 'u_L_2_d1d1d1d1': [-3, 0, 0], 'u_L_2_d1d1d1d2': [-3, 0, 0], 'u_L_2_d1d1d1d3': [-3, 0, 0], 'u_L_2_d1d1d2': [-2, 0, 0], 'u_L_2_d1d1d2d2': [-3, 0, 0], 'u_L_2_d1d1d2d3': [-3, 0, 0], 'u_L_2_d1d1d3': [-2, 0, 0], 'u_L_2_d1d1d3d3': [-3, 0, 0], 'u_L_2_d1d2': [-1, 0, 0], 'u_L_2_d1d2d2': [-2, 0, 0], 'u_L_2_d1d2d2d2': [-3, 0, 0], 'u_L_2_d1d2d2d3': [-3, 0, 0], 'u_L_2_d1d2d3': [-2, 0, 0], 'u_L_2_d1d2d3d3': [-3, 0, 0], 'u_L_2_d1d3': [-1, 0, 0], 'u_L_2_d1d3d3': [-2, 0, 0], 'u_L_2_d1d3d3d3': [-3, 0, 0], 'u_L_2_d2': [0, 0, 0], 'u_L_2_d2d2': [-1, 0, 0], 'u_L_2_d2d2d2': [-2, 0, 0], 'u_L_2_d2d2d2d2': [-3, 0, 0], 'u_L_2_d2d2d2d3': [-3, 0, 0], 'u_L_2_d2d2d3': [-2, 0, 0], 'u_L_2_d2d2d3d3': [-3, 0, 0], 'u_L_2_d2d3': [-1, 0, 0], 'u_L_2_d2d3d3': [-2, 0, 0], 'u_L_2_d2d3d3d3': [-3, 0, 0], 'u_L_2_d3': [0, 0, 0], 'u_L_2_d3d3': [-1, 0, 0], 'u_L_2_d3d3d3': [-2, 0, 0], 'u_L_2_d3d3d3d3': [-3, 0, 0], 'u_L_2_t': [1, -1, 0], 'u_L_2_t_d1': [0, -1, 0], 'u_L_2_t_d2': [0, -1, 0], 'u_L_2_t_d3': [0, -1, 0], 'u_L_3': [1, 0, 0], 'u_L_3_d1': [0, 0, 0], 'u_L_3_d1d1': [-1, 0, 0], 'u_L_3_d1d1d1': [-2, 0, 0], 'u_L_3_d1d1d1d1': [-3, 0, 0], 'u_L_3_d1d1d1d2': [-3, 0, 0], 'u_L_3_d1d1d1d3': [-3, 0, 0], 'u_L_3_d1d1d2': [-2, 0, 0], 'u_L_3_d1d1d2d2': [-3, 0, 0], 'u_L_3_d1d1d2d3': [-3, 0, 0], 'u_L_3_d1d1d3': [-2, 0, 0], 'u_L_3_d1d1d3d3': [-3, 0, 0], 'u_L_3_d1d2': [-1, 0, 0], 'u_L_3_d1d2d2': [-2, 0, 0], 'u_L_3_d1d2d2d2': [-3, 0, 0], 'u_L_3_d1d2d2d3': [-3, 0, 0], 'u_L_3_d1d2d3': [-2, 0, 0], 'u_L_3_d1d2d3d3': [-3, 0, 0], 'u_L_3_d1d3': [-1, 0, 0], 'u_L_3_d1d3d3': [-2, 0, 0], 'u_L_3_d1d3d3d3': [-3, 0, 0], 'u_L_3_d2': [0, 0, 0], 'u_L_3_d2d2': [-1, 0, 0], 'u_L_3_d2d2d2': [-2, 0, 0], 'u_L_3_d2d2d2d2': [-3, 0, 0], 'u_L_3_d2d2d2d3': [-3, 0, 0], 'u_L_3_d2d2d3': [-2, 0, 0], 'u_L_3_d2d2d3d3': [-3, 0, 0], 'u_L_3_d2d3': [-1, 0, 0], 'u_L_3_d2d3d3': [-2, 0, 0], 'u_L_3_d2d3d3d3': [-3, 0, 0], 'u_L_3_d3': [0, 0, 0], 'u_L_3_d3d3': [-1, 0, 0], 'u_L_3_d3d3d3': [-2, 0, 0], 'u_L_3_d3d3d3d3': [-3, 0, 0], 'u_L_3_t': [1, -1, 0], 'u_L_3_t_d1': [0, -1, 0], 'u_L_3_t_d2': [0, -1, 0], 'u_L_3_t_d3': [0, -1, 0], 'u_T_1': [1, 0, 0], 'u_T_1_d1': [0, 0, 0], 'u_T_1_d1d1': [-1, 0, 0], 'u_T_1_d1d1d1': [-2, 0, 0], 'u_T_1_d1d1d1d1': [-3, 0, 0], 'u_T_1_d1d1d1d2': [-3, 0, 0], 'u_T_1_d1d1d1d3': [-3, 0, 0], 'u_T_1_d1d1d2': [-2, 0, 0], 'u_T_1_d1d1d2d2': [-3, 0, 0], 'u_T_1_d1d1d2d3': [-3, 0, 0], 'u_T_1_d1d1d3': [-2, 0, 0], 'u_T_1_d1d1d3d3': [-3, 0, 0], 'u_T_1_d1d2': [-1, 0, 0], 'u_T_1_d1d2d2': [-2, 0, 0], 'u_T_1_d1d2d2d2': [-3, 0, 0], 'u_T_1_d1d2d2d3': [-3, 0, 0], 'u_T_1_d1d2d3': [-2, 0, 0], 'u_T_1_d1d2d3d3': [-3, 0, 0], 'u_T_1_d1d3': [-1, 0, 0], 'u_T_1_d1d3d3': [-2, 0, 0], 'u_T_1_d1d3d3d3': [-3, 0, 0], 'u_T_1_d2': [0, 0, 0], 'u_T_1_d2d2': [-1, 0, 0], 'u_T_1_d2d2d2': [-2, 0, 0], 'u_T_1_d2d2d2d2': [-3, 0, 0], 'u_T_1_d2d2d2d3': [-3, 0, 0], 'u_T_1_d2d2d3': [-2, 0, 0], 'u_T_1_d2d2d3d3': [-3, 0, 0], 'u_T_1_d2d3': [-1, 0, 0], 'u_T_1_d2d3d3': [-2, 0, 0], 'u_T_1_d2d3d3d3': [-3, 0, 0], 'u_T_1_d3': [0, 0, 0], 'u_T_1_d3d3': [-1, 0, 0], 'u_T_1_d3d3d3': [-2, 0, 0], 'u_T_1_d3d3d3d3': [-3, 0, 0], 'u_T_1_t': [1, -1, 0], 'u_T_1_t_d1': [0, -1, 0], 'u_T_1_t_d2': [0, -1, 0], 'u_T_1_t_d3': [0, -1, 0], 'u_T_2': [1, 0, 0], 'u_T_2_d1': [0, 0, 0], 'u_T_2_d1d1': [-1, 0, 0], 'u_T_2_d1d1d1': [-2, 0, 0], 'u_T_2_d1d1d1d1': [-3, 0, 0], 'u_T_2_d1d1d1d2': [-3, 0, 0], 'u_T_2_d1d1d1d3': [-3, 0, 0], 'u_T_2_d1d1d2': [-2, 0, 0], 'u_T_2_d1d1d2d2': [-3, 0, 0], 'u_T_2_d1d1d2d3': [-3, 0, 0], 'u_T_2_d1d1d3': [-2, 0, 0], 'u_T_2_d1d1d3d3': [-3, 0, 0], 'u_T_2_d1d2': [-1, 0, 0], 'u_T_2_d1d2d2': [-2, 0, 0], 'u_T_2_d1d2d2d2': [-3, 0, 0], 'u_T_2_d1d2d2d3': [-3, 0, 0], 'u_T_2_d1d2d3': [-2, 0, 0], 'u_T_2_d1d2d3d3': [-3, 0, 0], 'u_T_2_d1d3': [-1, 0, 0], 'u_T_2_d1d3d3': [-2, 0, 0], 'u_T_2_d1d3d3d3': [-3, 0, 0], 'u_T_2_d2': [0, 0, 0], 'u_T_2_d2d2': [-1, 0, 0], 'u_T_2_d2d2d2': [-2, 0, 0], 'u_T_2_d2d2d2d2': [-3, 0, 0], 'u_T_2_d2d2d2d3': [-3, 0, 0], 'u_T_2_d2d2d3': [-2, 0, 0], 'u_T_2_d2d2d3d3': [-3, 0, 0], 'u_T_2_d2d3': [-1, 0, 0], 'u_T_2_d2d3d3': [-2, 0, 0], 'u_T_2_d2d3d3d3': [-3, 0, 0], 'u_T_2_d3': [0, 0, 0], 'u_T_2_d3d3': [-1, 0, 0], 'u_T_2_d3d3d3': [-2, 0, 0], 'u_T_2_d3d3d3d3': [-3, 0, 0], 'u_T_2_t': [1, -1, 0], 'u_T_2_t_d1': [0, -1, 0], 'u_T_2_t_d2': [0, -1, 0], 'u_T_2_t_d3': [0, -1, 0], 'u_T_3': [1, 0, 0], 'u_T_3_d1': [0, 0, 0], 'u_T_3_d1d1': [-1, 0, 0], 'u_T_3_d1d1d1': [-2, 0, 0], 'u_T_3_d1d1d1d1': [-3, 0, 0], 'u_T_3_d1d1d1d2': [-3, 0, 0], 'u_T_3_d1d1d1d3': [-3, 0, 0], 'u_T_3_d1d1d2': [-2, 0, 0], 'u_T_3_d1d1d2d2': [-3, 0, 0], 'u_T_3_d1d1d2d3': [-3, 0, 0], 'u_T_3_d1d1d3': [-2, 0, 0], 'u_T_3_d1d1d3d3': [-3, 0, 0], 'u_T_3_d1d2': [-1, 0, 0], 'u_T_3_d1d2d2': [-2, 0, 0], 'u_T_3_d1d2d2d2': [-3, 0, 0], 'u_T_3_d1d2d2d3': [-3, 0, 0], 'u_T_3_d1d2d3': [-2, 0, 0], 'u_T_3_d1d2d3d3': [-3, 0, 0], 'u_T_3_d1d3': [-1, 0, 0], 'u_T_3_d1d3d3': [-2, 0, 0], 'u_T_3_d1d3d3d3': [-3, 0, 0], 'u_T_3_d2': [0, 0, 0], 'u_T_3_d2d2': [-1, 0, 0], 'u_T_3_d2d2d2': [-2, 0, 0], 'u_T_3_d2d2d2d2': [-3, 0, 0], 'u_T_3_d2d2d2d3': [-3, 0, 0], 'u_T_3_d2d2d3': [-2, 0, 0], 'u_T_3_d2d2d3d3': [-3, 0, 0], 'u_T_3_d2d3': [-1, 0, 0], 'u_T_3_d2d3d3': [-2, 0, 0], 'u_T_3_d2d3d3d3': [-3, 0, 0], 'u_T_3_d3': [0, 0, 0], 'u_T_3_d3d3': [-1, 0, 0], 'u_T_3_d3d3d3': [-2, 0, 0], 'u_T_3_d3d3d3d3': [-3, 0, 0], 'u_T_3_t': [1, -1, 0], 'u_T_3_t_d1': [0, -1, 0], 'u_T_3_t_d2': [0, -1, 0], 'u_T_3_t_d3': [0, -1, 0], 'v_bulk_normal_0': [1, -1, 0], 'w1_profile': [0, 0, 0], 'w1_profile_d1': [0, 0, 0], 'w1_profile_d1d1': [0, 0, 0], 'w1_profile_d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1d1': [0, 0, 0], 'w1_profile_d1d1d1d1d2': [0, 0, 0], 'w1_profile_d1d1d1d1d3': [0, 0, 0], 'w1_profile_d1d1d1d2': [0, 0, 0], 'w1_profile_d1d1d1d2d2': [0, 0, 0], 'w1_profile_d1d1d1d2d3': [0, 0, 0], 'w1_profile_d1d1d1d3': [0, 0, 0], 'w1_profile_d1d1d1d3d3': [0, 0, 0], 'w1_profile_d1d1d2': [0, 0, 0], 'w1_profile_d1d1d2d2': [0, 0, 0], 'w1_profile_d1d1d2d2d2': [0, 0, 0], 'w1_profile_d1d1d2d2d3': [0, 0, 0], 'w1_profile_d1d1d2d3': [0, 0, 0], 'w1_profile_d1d1d2d3d3': [0, 0, 0], 'w1_profile_d1d1d3': [0, 0, 0], 'w1_profile_d1d1d3d3': [0, 0, 0], 'w1_profile_d1d1d3d3d3': [0, 0, 0], 'w1_profile_d1d2': [0, 0, 0], 'w1_profile_d1d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2d2': [0, 0, 0], 'w1_profile_d1d2d2d2d3': [0, 0, 0], 'w1_profile_d1d2d2d3': [0, 0, 0], 'w1_profile_d1d2d2d3d3': [0, 0, 0], 'w1_profile_d1d2d3': [0, 0, 0], 'w1_profile_d1d2d3d3': [0, 0, 0], 'w1_profile_d1d2d3d3d3': [0, 0, 0], 'w1_profile_d1d3': [0, 0, 0], 'w1_profile_d1d3d3': [0, 0, 0], 'w1_profile_d1d3d3d3': [0, 0, 0], 'w1_profile_d1d3d3d3d3': [0, 0, 0], 'w1_profile_d2': [0, 0, 0], 'w1_profile_d2d2': [0, 0, 0], 'w1_profile_d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2d2': [0, 0, 0], 'w1_profile_d2d2d2d2d3': [0, 0, 0], 'w1_profile_d2d2d2d3': [0, 0, 0], 'w1_profile_d2d2d2d3d3': [0, 0, 0], 'w1_profile_d2d2d3': [0, 0, 0], 'w1_profile_d2d2d3d3': [0, 0, 0], 'w1_profile_d2d2d3d3d3': [0, 0, 0], 'w1_profile_d2d3': [0, 0, 0], 'w1_profile_d2d3d3': [0, 0, 0], 'w1_profile_d2d3d3d3': [0, 0, 0], 'w1_profile_d2d3d3d3d3': [0, 0, 0], 'w1_profile_d3': [0, 0, 0], 'w1_profile_d3d3': [0, 0, 0], 'w1_profile_d3d3d3': [0, 0, 0], 'w1_profile_d3d3d3d3': [0, 0, 0], 'w1_profile_d3d3d3d3d3': [0, 0, 0], 'zeta_c': [1, 0, 0], 'zeta_c_d1': [0, 0, 0], 'zeta_c_d2': [0, 0, 0], 'zeta_c_d3': [0, 0, 0], 'zeta_c_t': [1, -1, 0]}
 WAVE_NAMES = frozenset(['e_W', 'e_W_d1', 'e_W_d1d1', 'e_W_d1d1d1', 'e_W_d1d1d1d1', 'e_W_d1d1d1d2', 'e_W_d1d1d1d3', 'e_W_d1d1d2', 'e_W_d1d1d2d2', 'e_W_d1d1d2d3', 'e_W_d1d1d3', 'e_W_d1d1d3d3', 'e_W_d1d2', 'e_W_d1d2d2', 'e_W_d1d2d2d2', 'e_W_d1d2d2d3', 'e_W_d1d2d3', 'e_W_d1d2d3d3', 'e_W_d1d3', 'e_W_d1d3d3', 'e_W_d1d3d3d3', 'e_W_d2', 'e_W_d2d2', 'e_W_d2d2d2', 'e_W_d2d2d2d2', 'e_W_d2d2d2d3', 'e_W_d2d2d3', 'e_W_d2d2d3d3', 'e_W_d2d3', 'e_W_d2d3d3', 'e_W_d2d3d3d3', 'e_W_d3', 'e_W_d3d3', 'e_W_d3d3d3', 'e_W_d3d3d3d3', 'e_W_t', 'e_W_tt', 'grad_theta_1', 'grad_theta_2', 'grad_theta_3', 'theta', 'theta_d1d1', 'theta_d1d1d1', 'theta_d1d1d1d1', 'theta_d1d1d1d2', 'theta_d1d1d1d3', 'theta_d1d1d2', 'theta_d1d1d2d2', 'theta_d1d1d2d3', 'theta_d1d1d3', 'theta_d1d1d3d3', 'theta_d1d2', 'theta_d1d2d2', 'theta_d1d2d2d2', 'theta_d1d2d2d3', 'theta_d1d2d3', 'theta_d1d2d3d3', 'theta_d1d3', 'theta_d1d3d3', 'theta_d1d3d3d3', 'theta_d2d2', 'theta_d2d2d2', 'theta_d2d2d2d2', 'theta_d2d2d2d3', 'theta_d2d2d3', 'theta_d2d2d3d3', 'theta_d2d3', 'theta_d2d3d3', 'theta_d2d3d3d3', 'theta_d3d3', 'theta_d3d3d3', 'theta_d3d3d3d3', 'theta_t', 'u_1', 'u_1_d1', 'u_1_d1d1', 'u_1_d1d1d1', 'u_1_d1d1d1d1', 'u_1_d1d1d1d2', 'u_1_d1d1d1d3', 'u_1_d1d1d2', 'u_1_d1d1d2d2', 'u_1_d1d1d2d3', 'u_1_d1d1d3', 'u_1_d1d1d3d3', 'u_1_d1d2', 'u_1_d1d2d2', 'u_1_d1d2d2d2', 'u_1_d1d2d2d3', 'u_1_d1d2d3', 'u_1_d1d2d3d3', 'u_1_d1d3', 'u_1_d1d3d3', 'u_1_d1d3d3d3', 'u_1_d2', 'u_1_d2d2', 'u_1_d2d2d2', 'u_1_d2d2d2d2', 'u_1_d2d2d2d3', 'u_1_d2d2d3', 'u_1_d2d2d3d3', 'u_1_d2d3', 'u_1_d2d3d3', 'u_1_d2d3d3d3', 'u_1_d3', 'u_1_d3d3', 'u_1_d3d3d3', 'u_1_d3d3d3d3', 'u_1_t', 'u_1_t_d1', 'u_1_t_d2', 'u_1_t_d3', 'u_1_tt', 'u_2', 'u_2_d1', 'u_2_d1d1', 'u_2_d1d1d1', 'u_2_d1d1d1d1', 'u_2_d1d1d1d2', 'u_2_d1d1d1d3', 'u_2_d1d1d2', 'u_2_d1d1d2d2', 'u_2_d1d1d2d3', 'u_2_d1d1d3', 'u_2_d1d1d3d3', 'u_2_d1d2', 'u_2_d1d2d2', 'u_2_d1d2d2d2', 'u_2_d1d2d2d3', 'u_2_d1d2d3', 'u_2_d1d2d3d3', 'u_2_d1d3', 'u_2_d1d3d3', 'u_2_d1d3d3d3', 'u_2_d2', 'u_2_d2d2', 'u_2_d2d2d2', 'u_2_d2d2d2d2', 'u_2_d2d2d2d3', 'u_2_d2d2d3', 'u_2_d2d2d3d3', 'u_2_d2d3', 'u_2_d2d3d3', 'u_2_d2d3d3d3', 'u_2_d3', 'u_2_d3d3', 'u_2_d3d3d3', 'u_2_d3d3d3d3', 'u_2_t', 'u_2_t_d1', 'u_2_t_d2', 'u_2_t_d3', 'u_2_tt', 'u_3', 'u_3_d1', 'u_3_d1d1', 'u_3_d1d1d1', 'u_3_d1d1d1d1', 'u_3_d1d1d1d2', 'u_3_d1d1d1d3', 'u_3_d1d1d2', 'u_3_d1d1d2d2', 'u_3_d1d1d2d3', 'u_3_d1d1d3', 'u_3_d1d1d3d3', 'u_3_d1d2', 'u_3_d1d2d2', 'u_3_d1d2d2d2', 'u_3_d1d2d2d3', 'u_3_d1d2d3', 'u_3_d1d2d3d3', 'u_3_d1d3', 'u_3_d1d3d3', 'u_3_d1d3d3d3', 'u_3_d2', 'u_3_d2d2', 'u_3_d2d2d2', 'u_3_d2d2d2d2', 'u_3_d2d2d2d3', 'u_3_d2d2d3', 'u_3_d2d2d3d3', 'u_3_d2d3', 'u_3_d2d3d3', 'u_3_d2d3d3d3', 'u_3_d3', 'u_3_d3d3', 'u_3_d3d3d3', 'u_3_d3d3d3d3', 'u_3_t', 'u_3_t_d1', 'u_3_t_d2', 'u_3_t_d3', 'u_3_tt'])
 
@@ -795,40 +794,94 @@ def traction_pairing(inputs, case, model, *, flip=False):
     independent_force=tree(independent_force,lambda e:inputs.physical_fields(e.subs(model['substitutions'],simultaneous=True)))
     velocities=[wave_jet(inputs.a(f'u_{i}_t')) for i in range(1,4)]
     e_velocity=wave_jet(inputs.a('e_W_t'))
     slab_power=sum(r*v for r,v in zip(model['closed']['U'],velocities))+model['closed']['E_W']*e_velocity
     face_power=sum(r*v for r,v in zip(independent_force['U'],velocities))+independent_force['E_W']*e_velocity
     kinetic_stored=slab_power-face_power
     pairing=tree({'SLAB_POWER':slab_power,'KINETIC_STORED_POWER':kinetic_stored,
                   'TRACTION_POWER':power,'FACE_GENERALIZED_POWER':face_power},lambda e:retained_shape(e,inputs))
     residual=pairing['SLAB_POWER']-pairing['KINETIC_STORED_POWER']-pairing['TRACTION_POWER']
     return tree(covectors,lambda e:retained_shape(e,inputs)),pairing,residual
 
 
+def publication_compact(value):
+    """Factor VALUE entries while preserving reciprocal and Integral boundaries.
+
+    Temporary atoms shield denominators and calculus objects from factor_terms;
+    they are substituted back immediately, never serialized or used as CSE.
+    """
+    @lru_cache(maxsize=None)
+    def compact(expression):
+        if isinstance(expression, sp.MatrixBase):
+            return expression.applyfunc(compact)
+        if isinstance(expression, sp.Tuple):
+            return sp.Tuple(*(compact(v) for v in expression))
+        if isinstance(expression, Str) or expression.is_Atom:
+            return expression
+        if isinstance(expression, sp.Integral):
+            return sp.Integral(compact(expression.function), *expression.limits)
+        if (expression.is_Pow and expression.exp.is_negative) or isinstance(
+                expression, (sp.Function, sp.Derivative)):
+            return expression
+        protected = {}
+
+        @lru_cache(maxsize=None)
+        def shield(node):
+            if isinstance(node, sp.Integral) or isinstance(node, (sp.Function, sp.Derivative)) or (
+                    node.is_Pow and node.exp.is_negative):
+                token = sp.Dummy(commutative=node.is_commutative)
+                protected[token] = compact(node)
+                return token
+            if not node.args:
+                return node
+            return node.func(*(shield(a) for a in node.args))
+
+        shielded = shield(expression)
+        grouped = sp.collect(shielded, list(protected), exact=True)
+        result = sp.factor_terms(grouped, fraction=False).xreplace(protected)
+        if result.atoms(sp.Dummy) - expression.atoms(sp.Dummy):
+            raise ValueError('publication temporary escaped')
+        # Factoring can expose an identically zero coefficient and erase the
+        # reciprocal it multiplied. Keep that original subexpression intact.
+        before_poles = {p for p in expression.atoms(sp.Pow) if p.exp.is_negative}
+        after_poles = {p for p in result.atoms(sp.Pow) if p.exp.is_negative}
+        return result if before_poles == after_poles else expression
+
+    return compact(value)
+
+
 def publish(inputs, fold, objects):
-    candidates={k:{'value':v,'display':sp.sstr(v),'value_kind':'COMPUTED_OBJECT',
-                   'class':'DERIVED','step':'S11c-c2','route':'F9A_ABSENT'} for k,v in objects.items()}
+    if set(objects) != EXPORT_ROOTS:
+        raise ValueError('publication root membership')
+    compact_objects = {}
+    for key, value in objects.items():
+        compact_objects[key] = sp.Tuple(*(sp.Tuple(case, sp.Tuple(*(
+            sp.Tuple(label, publication_compact(item) if str(label) == 'VALUE' else item)
+            for label, item in payload))) for case, payload in value))
+    candidates={k:{'value':v,'display':k,'value_kind':'COMPUTED_OBJECT',
+                   'class':'DERIVED','step':'S11c-c2','route':'F9A_ABSENT'} for k,v in compact_objects.items()}
     all_atoms=set().union(*(v.atoms(sp.Symbol,AppliedUndef) for v in objects.values()))
     for atom in all_atoms:
         declaration=atom.func if isinstance(atom,AppliedUndef) else atom
         if declaration in NEW_DIMENSIONS:
             name=declaration.__name__ if isinstance(atom,AppliedUndef) else str(atom)
             candidates[name]={'value':declaration,'display':str(declaration),'value_kind':'COMPUTED_OBJECT',
                               'class':'COORDINATE','step':'S11c-c2','route':'F9A_ABSENT',
                               'dimension_key':name+'Dimension'}
             candidates[name+'Dimension']={'value':sp.ImmutableMatrix(NEW_DIMENSIONS[declaration]),
                 'display':str(NEW_DIMENSIONS[declaration]),'value_kind':'COMPUTED_OBJECT','class':'STRUCTURAL','step':'S11c-c2','route':'F9A_ABSENT'}
     collisions=set(candidates)&set(fold)
     if collisions:
         raise ValueError(('F9 collision',sorted(collisions)))
+    print('EXPORT_F9_COLLISIONS = ' + repr(sorted(collisions)), flush=True)
     combined=dict(fold)|candidates
     closure=check_consumer(combined,EXPORT_ROOTS)['closure']-set(fold)
     delta={k:candidates[k] for k in sorted(closure)}
     guard=assert_delta_is_minimal(delta,closure)
     digests={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (
         HERE,ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_c1_exports.py',
         ROOT/'directives/S11c_c2_SHARED_PHYSICS.md',ROOT/'scripts/ledger_fold.py')}
     lines=['# Generated own-rows delta; S11c-c2.','from types import MappingProxyType',
         'import sympy as sp','from sympy.core.symbol import Str',
         'from sympy.functions.elementary.piecewise import ExprCondPair',
         'def _restore(s):',"    return eval(s, {'__builtins__': {}, **vars(sp), 'Str': Str, 'ExprCondPair': ExprCondPair})",
         'IMPORT_KEYS = '+repr(IMPORT_KEYS),'BUILD_INPUT_DIGESTS = MappingProxyType('+repr(digests)+')','_LEDGER = {']
@@ -840,26 +893,130 @@ def publish(inputs, fold, objects):
         lines.append(repr(k)+': {'+', '.join(fields)+'},')
     lines += ['}','LEDGER = MappingProxyType({k: MappingProxyType(v) for k,v in _LEDGER.items()})','del _LEDGER','']
     path=ROOT/'scripts/S11c_c2_exports.py'
     temporary=path.with_suffix('.py.tmp')
     temporary.write_text('\n'.join(lines))
     namespace={}
     exec(compile(temporary.read_text(),str(path),'exec'),namespace)
     restored=namespace['LEDGER']
     comparison=sp.Tuple(*(sp.sympify(delta[k]['value']==restored[k]['value']) for k in delta))
     emit('EXPORT_ROUNDTRIP',comparison,inputs)
     if any(delta[k]['value']!=restored[k]['value'] for k in delta):
         raise ValueError('serialization roundtrip')
+    # Separate semantic guard: the generated module's restored values versus
+    # the original emitted payload, with strict containers before arithmetic.
+    evidence = []
+
+    def record(path, check, result):
+        line = f'EXPORT_SEMANTIC {path} {check} = {result}'
+        print(line, flush=True)
+        evidence.append(line)
+        if result is not True:
+            raise ValueError(('emitted/compact mismatch', path, check, result))
+
+    def entries(value, path):
+        record(path, 'mapping_tuple', isinstance(value, sp.Tuple))
+        record(path, 'pair_arities', all(isinstance(p, sp.Tuple) and len(p) == 2 for p in value))
+        result = dict(value)
+        record(path, 'unique_keys', len(result) == len(value))
+        return result
+
+    def poles(value):
+        # Exact bases AND exponents: stronger than a pole-set comparison.
+        return {p for p in value.atoms(sp.Pow) if p.exp.is_negative}
+
+    def semantic(original, decoded, path):
+        if isinstance(original, Mapping):
+            record(path, 'mapping_type', type(decoded) is type(original))
+            record(path, 'mapping_keys', set(original) == set(decoded))
+            for key in original:
+                semantic(original[key], decoded[key], f'{path}/{key}')
+        elif isinstance(original, sp.MatrixBase):
+            record(path, 'matrix_type', type(decoded) is type(original))
+            record(path, 'matrix_shape', decoded.shape == original.shape)
+            for i in range(original.rows):
+                for j in range(original.cols):
+                    semantic(original[i, j], decoded[i, j], f'{path}[{i},{j}]')
+        elif isinstance(original, (sp.Tuple, tuple, list)):
+            record(path, 'tuple_type', type(decoded) is type(original))
+            record(path, 'tuple_arity', len(decoded) == len(original))
+            if original and all(isinstance(p, sp.Tuple) and len(p) == 2 and
+                                isinstance(p[0], Str) for p in original):
+                left, right = entries(original, path + '/emitted'), entries(decoded, path + '/restored')
+                record(path, 'mapping_keys', set(left) == set(right))
+                for key in left:
+                    semantic(left[key], right[key], f'{path}/{key}')
+            else:
+                for i in range(len(original)):
+                    semantic(original[i], decoded[i], f'{path}[{i}]')
+        elif isinstance(original, Str):
+            record(path, 'Str_label', type(decoded) is Str and original == decoded)
+        else:
+            record(path, 'algebraic_leaf', isinstance(original, sp.Expr) and isinstance(decoded, sp.Expr))
+            record(path, 'reciprocal_powers_unchanged', poles(original) == poles(decoded))
+            # Expand integrands with their unchanged reciprocal/calculus atoms
+            # protected, then protect identical Integral atoms as in difference.
+            # These local check-only dummies never enter a published value.
+            atoms = (poles(original) | poles(decoded) |
+                     original.atoms(sp.Function, sp.Derivative) |
+                     decoded.atoms(sp.Function, sp.Derivative))
+            protected = {v: sp.Dummy('publicationCoefficient') for v in atoms}
+            unprotect = {v: k for k, v in protected.items()}
+
+            @lru_cache(maxsize=None)
+            def normalize_integrals(value):
+                replacements = {}
+                for integral in value.atoms(sp.Integral):
+                    integrand = sp.expand(normalize_integrals(integral.function).xreplace(protected))
+                    integrand = integrand.xreplace(unprotect)
+                    replacements[integral] = (sp.S.Zero if integrand == 0 else
+                                              sp.Integral(integrand, *integral.limits))
+                return value.xreplace(replacements)
+
+            expression = normalize_integrals(decoded) - normalize_integrals(original)
+            integrals = {v: sp.Dummy('publicationIntegral') for v in expression.atoms(sp.Integral)}
+            residual = sp.expand(expression.xreplace(integrals))
+            print(f'EXPORT_SEMANTIC {path} expanded_difference = {residual}', flush=True)
+            evidence.append(f'EXPORT_SEMANTIC {path} expanded_difference = {residual}')
+            record(path, 'expanded_difference_is_zero', residual == 0)
+
+    sizes = {}
+    expected_cases = {cas((a, d)) for a in ANCHORINGS for d in DENSITIES}
+    for key in sorted(EXPORT_ROOTS):
+        original_cases = entries(objects[key], key + '/emitted')
+        decoded_cases = entries(restored[key]['value'], key + '/restored')
+        record(key, 'case_keys', set(original_cases) == set(decoded_cases) == expected_cases)
+        sizes[key] = {}
+        for case in original_cases:
+            path_name = key + '/' + '/'.join(map(str, case))
+            original = entries(original_cases[case], path_name + '/emitted')
+            decoded = entries(decoded_cases[case], path_name + '/restored')
+            record(path_name, 'payload_keys', set(original) == set(decoded))
+            for label in original:
+                if str(label) == 'VALUE':
+                    semantic(original[label], decoded[label], path_name + '/VALUE')
+                    counts = {'emitted_srepr_bytes': len(sp.srepr(original[label]).encode()),
+                              'compact_srepr_bytes': len(sp.srepr(decoded[label]).encode())}
+                    sizes[key]['/'.join(map(str, case))] = counts
+                    print('EXPORT_VALUE_BYTES ' + path_name + ' = ' + repr(counts), flush=True)
+                else:
+                    record(path_name + '/' + str(label), 'metadata_exact', original[label] == decoded[label])
+    restored_closure = check_consumer(dict(fold) | dict(restored), EXPORT_ROOTS)['closure'] - set(fold)
+    restored_guard = assert_delta_is_minimal(restored, restored_closure)
+    print('EXPORT_RESTORED_CLOSURE = ' + repr(sorted(restored_closure)), flush=True)
+    print('EXPORT_RESTORED_MINIMALITY = ' + repr({k: sorted(v) for k, v in restored_guard.items()}), flush=True)
+    record('delta', 'closure_unchanged', restored_closure == closure)
     temporary.replace(path)
-    return {'keys':sorted(delta),'own_closure':sorted(closure),'guard':{k:sorted(v) for k,v in guard.items()},'digests':digests}
+    return {'keys':sorted(delta),'own_closure':sorted(closure),'guard':{k:sorted(v) for k,v in guard.items()},'digests':digests,
+            'publication_semantic': evidence, 'publication_value_bytes': sizes}
 
 
 def run():
     started=time.monotonic()
     fold,audit=load_model(str(ROOT/'scripts/S11c_b_exports.py'),str(ROOT/'scripts/S11c_c1_exports.py'))
     closure_audit=check_consumer(fold,IMPORT_KEYS)
     lookup_audit=assert_lookups_equal_manifest(bind_inputs,fold,IMPORT_KEYS)
     inputs=lookup_audit['result']
     measurements={'fold':audit,'fold_rows':len(fold),'import_keys':sorted(IMPORT_KEYS),
         'lookups':sorted(lookup_audit['lookups']),'closure':sorted(closure_audit['closure'])}
     checkpoint=ROOT/('_measurements/S11c_c2_sympy_triage_guard_evidence.json' if os.environ.get('S11CC2_PACKAGE')=='TRIAGE' else '_measurements/S11c_c2_sympy_guard_evidence.json')
     checkpoint.write_text(json.dumps(measurements,indent=2))
@@ -883,26 +1040,25 @@ def run():
                 ('CLOSED_SLAB_OPERATOR_PARITY_BLOCKS',{
                     'FACE_SUM':tree(difference(model['faces'][1],tree(model['faces'][-1],lambda e:-e)),lambda e:e/2),
                     'FACE_DIFFERENCE':tree(difference(model['faces'][1],model['faces'][-1]),lambda e:e/2)}),
                 ('CLOSED_COUPLING_KERNEL',model['closed_kernel']),
                 ('CLOSED_COUPLING_KERNEL_TERM_ORIGINS',{s:extract(v,inputs) for s,v in model['faces'].items()}),
                 ('SELF_ENERGY_CLOSED_EXTRACT_OPERAND',model['closed_kernel']),
                 ('SELF_ENERGY_OPEN_EXTRACT_OPERAND',model['open_kernel']),
                 ('SELF_ENERGY_INCREMENT',model['increment']),
                 ('FOLD_SYMBOL_MAP',model['maps']),
             ):
                 payload=emit(quantity,obj,inputs,case)
                 export_key={'CLOSED_SLAB_OPERATOR':'s11cc2ClosedSlabOperator',
-                            'CLOSED_COUPLING_KERNEL':'s11cc2ClosedCouplingKernel',
-                            'SELF_ENERGY_INCREMENT':'s11cc2SelfEnergyIncrement'}.get(quantity)
+                            'CLOSED_COUPLING_KERNEL':'s11cc2ClosedCouplingKernel'}.get(quantity)
                 if export_key:
                     exports[export_key][case]=payload
             # Regression-only bind of the predecessor's already-extracted kernel.
             regression=named(inputs.open_kernel[case],'COMPLETE_OPERATOR_BLOCKS')
             closed_regression=tree(regression,lambda e:inputs.physical_fields(e.subs(model['substitutions'],simultaneous=True)))
             # The predecessor uses jet trial/test symbols; align those explicitly
             # in regression_coordinates before comparing to the function ansatz.
             aligned=regression_coordinates(closed_regression,inputs)
             emit('ORDERING_EXTRACT_FIRST_OPERAND',aligned,inputs,case)
             emit('ORDERING_COMMUTATOR',difference(model['closed_kernel'],aligned),inputs,case)
             covectors,pairing,residual=traction_pairing(inputs,case,model)
             emit('TRACTION_MECHANICAL_CONTRIB',covectors,inputs,case)
48:EXPORT_ROOTS = frozenset((
63:def cas(value):
75:def named(value, key):
79:def axes(value):
83:def cases(value):
87:def tree(value, fn):
99:def difference(a, b):
111:def atom_named(value, name):
118:def coordinate(name, dimension):
133:def field(name, point, dimension):
139:def wave_jet(atom, point=X, sector=None):
176:def bind_inputs(fold):
308:def restricted(expression, sector):
325:def extract(rows, inputs):
345:def matrix_evaluate(expression, replacements):
365:def kernel_bridge(inputs, anchoring, face, density, overrides):
419:def normal_slope_control(inputs,anchoring,face):
427:def fourier_profiles(inputs,expression,kout,kin):
437:def outgoing_spectral(inputs,kout,kin):
453:def kernel_apply(inputs, diagonal, off_diagonal, source, kout, kin, second=sp.S.Zero):
472:def integral(integrand,*limits):
476:def build_face(inputs, anchoring, face, density, overrides=None, mu_override=None, velocity_override=None):
526:def expanded_rows(operator):
531:def build_case(inputs, anchoring, density, *, overrides=None, routing=None, mu_override=None, velocity_override=None, input_map=None):
571:def retained_shape(expression,inputs):
578:def shape_coefficients(expression,eta,sigma):
635:def dimension(expression):
677:def dimensions(value):
688:def grades(expression, epsilon, eta, sigma):
720:def grade_object(value, inputs):
731:def profile_bindings(inputs):
745:def emit(quantity, value, inputs, case=(), *, key=None):
764:def control(inputs, case, baseline, name, **kwargs):
771:def traction_pairing(inputs, case, model, *, flip=False):
806:def publication_compact(value):
852:def publish(inputs, fold, objects):
853:    if set(objects) != EXPORT_ROOTS:
858:            sp.Tuple(label, publication_compact(item) if str(label) == 'VALUE' else item)
877:    closure=check_consumer(combined,EXPORT_ROOTS)['closure']-set(fold)
984:    for key in sorted(EXPORT_ROOTS):
1003:    restored_closure = check_consumer(dict(fold) | dict(restored), EXPORT_ROOTS)['closure'] - set(fold)
1013:def run():
1030:    exports={k:{} for k in EXPORT_ROOTS}
1051:                export_key={'CLOSED_SLAB_OPERATOR':'s11cc2ClosedSlabOperator',
1053:                if export_key:
1054:                    exports[export_key][case]=payload
1108:    measurements['export']=publish(inputs,fold,{k:cas(v) for k,v in exports.items()})
1115:def modulus_form(expression,inputs):
1136:def representation_pullback(value,inputs,density):
1174:def regression_coordinates(value,inputs):
research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
?? research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md

exec
/bin/bash -lc "git diff --unified=0 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py | rg '"'^@@|''^[+-]def '"' && git diff --word-diff=porcelain --unified=0 8f3a017f -- research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py | rg '"'^@@'"' && git show --format=fuller --name-status --no-renames 8f3a017f | sed -n '1,80p' && git show --format=fuller --name-status --no-renames a5f7a06c | sed -n '1,70p' && git show --format=fuller --name-status --no-renames aa76105a | sed -n '1,100p' && git show --format=fuller --name-status --no-renames d2befb7c | sed -n '1,50p'" in /var/projects/toy_physics
 succeeded in 0ms:
@@ -50 +49,0 @@ EXPORT_ROOTS = frozenset((
@@ -806,0 +806,46 @@ def traction_pairing(inputs, case, model, *, flip=False):
+def publication_compact(value):
@@ -808,2 +853,9 @@ def publish(inputs, fold, objects):
@@ -822,0 +875 @@ def publish(inputs, fold, objects):
@@ -851,0 +905,103 @@ def publish(inputs, fold, objects):
@@ -853 +1009,2 @@ def publish(inputs, fold, objects):
@@ -895,2 +1052 @@ def run():
@@ -50 +49,0 @@ EXPORT_ROOTS = frozenset((
@@ -806,0 +806,46 @@ def traction_pairing(inputs, case, model, *, flip=False):
@@ -808,2 +853,9 @@ def publish(inputs, fold, objects):
@@ -822,0 +875 @@ def publish(inputs, fold, objects):
@@ -851,0 +905,103 @@ def publish(inputs, fold, objects):
@@ -853 +1009,2 @@ def publish(inputs, fold, objects):
@@ -895,2 +1052 @@ def run():
commit 8f3a017fac947da9782ec694d0b798a93a307935
Author:     Trevor Norris <trev.norris@gmail.com>
AuthorDate: Sat Sep 5 19:35:24 2026 -0600
Commit:     Trevor Norris <trev.norris@gmail.com>
CommitDate: Sat Sep 5 19:35:24 2026 -0600

    S11c-c2 self-energy fold — physics review (2 legs) ADJUDICATED: fold physics SOUND (0 defects)
    
    Reviewed baseline of the astra-built c2 SymPy fold, committed BEFORE the STEP-B
    export repair overwrites it. Script is Codex(astra)-written -> legs = fresh Claude
    agent + Grok, identical prompt, orchestrator-launched, blindness by absence
    (builder report + provenance relocated during review).
    
    Both legs agree the fold WIRING + A/C/D1-D6 are SOUND (dp not J_s, operator inverse,
    kernel bridge, V_s->face_velocity, computed jets, eps-strip; FORM ablations respond).
    Legs split on B/E/F/G; adjudicated by my own computation (rule 13):
    
    - F (uniform limit): Grok "surviving coupling" = FALSE POSITIVE. verify_F.py/.out:
      removing the bare open-slot bookkeeping leaves coeff*Test*Integral(0,...) -> .doit()==0
      in all 4 cases -> the genuine closure-induced coupling DECOUPLES exactly. Grok saw
      Integral(...) without evaluating its integrand to 0.
    - B: the O(eps) grade Grok flagged = F's -extract(open) representation piece, not a
      physics defect; genuine self-energy is O(eps.eta).
    - E (N6): verify_EG.py/.out: REP_INVARIANCE_RESIDUAL.subs(sigma_W->0)==0 for every block,
      both densities -> leading O(eps)/O(eps.eta) rep-invariance HOLDS; remnant is sigma_W-sector
      only (drain-projection carry-forward, >=64 GB deferred), not a load-bearing c2 defect.
    - G (adjointness): FALSE POSITIVE. increment reverse blocks identically zero (one-way
      self-energy, expected -- dp lives only in theta/mechanical rows); both blocks emitted;
      builder report L70 honestly says "no independent adjointness construction" = spec 3b.
      Grok refuted a claim the builder never made.
    
    VERDICT: c2 fold physics SOUND, zero confirmed defects. Carry-forwards (record, not
    defects): F wording (say "genuine coupling decouples", not "increment vanishes"; light
    5e/3c spec clarification owed); E leading-order only + sigma_W deferred; G directional
    self-energy; the two S11c-b sign conventions -> blind-Wolfram comparator.
    
    NEXT: STEP B publication-only export repair (publish@807/call@952 only) -- drop the
    increment to EMIT-only, export BOTH closed operators transparent-factored + equiv-check.
    c1 STANDS.
    
    Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_01BT9RBnimH6baW1cqX8RKpE

A	research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out
A	research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py
A	research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out
A	research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py
A	research/pde_ledger_v3/_measurements/S11c_c2_physics_claude_leg.md
A	research/pde_ledger_v3/_measurements/S11c_c2_physics_grok_leg.txt
A	research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md
A	research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md
A	research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json
A	research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json
A	research/pde_ledger_v3/directives/_legs/S11c_c2_physics_review_prompt.md
A	research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
commit a5f7a06c9ab827ec573f02c716a6bd98e9a065a1
Author:     Trevor Norris <trev.norris@gmail.com>
AuthorDate: Sat Sep 5 19:56:54 2026 -0600
Commit:     Trevor Norris <trev.norris@gmail.com>
CommitDate: Sat Sep 5 19:56:54 2026 -0600

    S11c-c2 export-repair directive — GATED (2 decision legs, one pass, folded)
    
    Publication-only export-repair build directive (orchestrator-written) -> 2 decision
    legs (Codex gpt-5.6-sol xhigh + Grok), one two-leg pass, folded once (G2, not
    iterate-to-green). Both legs converged: scope + 3 edit sites, R1 membership (d binds
    BOTH closed operators; increment EMIT-only), R3 (display safe), R4, overstep fence =
    SOUND. One blocker: R2 (representation/verification), which I under-specified.
    
    Folded (verified against the code myself, rule 13):
    - R2: compact ONLY per-case VALUE leaves, keep the cased tree (cases()/named parseable);
      factor/cancel/together/collect only, NOT expand, NOT CSE (dummies aren't declarations
      -> incomplete-skeleton diff), NOT opaque; preserve the singular/pole locus (feeds d's
      resonances). TWO guards: KEEP the structural srepr roundtrip, ADD a strict-recursive
      Integral-aware emitted<->compact semantic check (shape/length-checked, not naive zip;
      run against the generated module; hard-stop on mismatch).
    - nits: allow a publication helper adjacent to publish; short bounded display (supersedes
      S9 D3 / bind-closure D4 for giant rows); .out identity = tag/key names not payload diff;
      guard_evidence.json checkpoint write = expected harness side-effect; no d IMPORT_KEYS
      manifest yet (both operators are the declared prospective binds).
    
    Reports: _measurements/S11c_c2_export_repair_decision_{codex.md,grok.txt}; gate record
    _measurements/S11c_c2_export_repair_decision_gate.md. NEXT: astra build (detached) ->
    2 re-review legs (fresh Claude + Grok) -> adjudicate -> commit. c1 STANDS.
    
    Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_01BT9RBnimH6baW1cqX8RKpE

A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_codex.md
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_gate.md
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_grok.txt
A	research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md
A	research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_decision_prompt.md
commit aa76105a12c11d48e63fe321eee5ebca5740dde1
Author:     Trevor Norris <trev.norris@gmail.com>
AuthorDate: Sat Sep 5 22:23:45 2026 -0600
Commit:     Trevor Norris <trev.norris@gmail.com>
CommitDate: Sat Sep 5 22:23:45 2026 -0600

    S11c-c2 export repair COMMITTED — 60 MB -> 21.4 MB (publication-only, re-reviewed clear)
    
    Publication-only export repair (astra build, gpt-6-astra), 2-leg re-reviewed.
    scripts/S11c_c2_exports.py: 60,516,900 -> 22,441,522 bytes (63% smaller).
    
    WHAT CHANGED (export path only; physics untouched):
    - Dropped s11cc2SelfEnergyIncrement from EXPORT_ROOTS + the export_key map -> it
      is EMIT-only (comparator representation, not a downstream binding).
    - New publication_compact(): transparent-factors only the per-case VALUE leaves,
      preserves the cased {VALUE,MULTIGRADE,DIMENSION,bindings} tree, protects Integrals
      + denominators, reverts if factoring would change the pole set, raises if any
      temporary Dummy escapes (no CSE/opaque).
    - New strict-recursive semantic guard: emitted vs restored-compact, Integral-aware,
      container-shape-checked, leafwise expand(diff)==0, hard-stop on mismatch; the
      srepr roundtrip is KEPT separately.
    - display shrunk from a full sstr dup to the root name.
    No construction function changed (build_case/extract/kernel/retained_shape/emit/...).
    
    RE-REVIEW (astra-written -> fresh Claude agent + Grok):
    - Grok: all 6 checks CLEAR; guards bite (out+1 compaction corruption + membership
      re-add/drop all hard-stop); physics unchanged (baseline vs current TRIAGE
      byte-identical); membership minimal; transparency clean; no overstep.
    - Fresh Claude agent: wedged on reaped background waiters; its cross-process re-run
      method produced two FALSE POSITIVES (semantic FAIL = fresh-Integral-dummy artifact
      k*d1*X - k*d2*X; B2 "clean" = its own no-op corruption). Its clean checks
      (membership/transparency/overstep) corroborate Grok.
    - Orchestrator rule-13 verification (adjudication_verify.py/.out), same-process:
      FAITHFUL 5/5 (compact==emitted), GUARD BITES 5/5 on x2, 5/5 on +marker, and True
      on an in-integrand corruption. Refutes both Agent false positives; confirms Grok.
    
    VERDICT: export repair SOUND. Membership correct+minimal, compaction faithful (no
    info loss), guards bite, pole locus preserved, physics unchanged, no overstep.
    
    NEXT: light S11c_c2_SHARED_PHYSICS 5e/3c wording clarification; then the blind
    Wolfram engine -> its 2 legs -> c2 T7 comparator + reconcile -> c2 step record.
    c1 STANDS.
    
    Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_01BT9RBnimH6baW1cqX8RKpE

A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.out
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.py
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out
A	research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_grok.txt
A	research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json
A	research/pde_ledger_v3/_measurements/S11c_c2_sympy_progress.json
A	research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_rereview_prompt.md
A	research/pde_ledger_v3/scripts/S11c_c2_exports.py
M	research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
commit d2befb7c84016dd7bd9ecac9e1b053c389e6b3ee
Author:     Trevor Norris <trev.norris@gmail.com>
AuthorDate: Sun Sep 6 00:31:17 2026 -0600
Commit:     Trevor Norris <trev.norris@gmail.com>
CommitDate: Sun Sep 6 00:31:17 2026 -0600

    STATUS: c2 STEP A (physics SOUND) + STEP B/C (export repair 60->21 MB) DONE + committed
    
    New top clause (2026-09-06): c2 fold physics 2-leg-reviewed SOUND 8f3a017f
    (Grok F/G/E all false positives, refuted by my rule-13 compute); export repair
    committed aa76105a (publication-only, re-reviewed clear). Carry-forwards for the
    step record (F/E/G interpretation, not defects) + the astra-keep-for-code
    assessment + the session's background-reap lesson recorded. NEXT = light
    S11c_c2_SHARED_PHYSICS 5e/3c wording fix -> blind Wolfram engine -> T7
    comparator+reconcile -> c2 step record. c1 STANDS.
    
    Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
    Claude-Session: https://claude.ai/code/session_01BT9RBnimH6baW1cqX8RKpE

M	STATUS.md

codex
The script diff is tightly scoped: zero-context hunks occur only at `EXPORT_ROOTS`, the added `publication_compact`, `publish`, and the `run()` export-key map. No construction-function hunk appears. I’m now loading the generated 22,441,522-byte module without bytecode writes to verify its actual runtime structure, not just its text.
exec
/bin/bash -lc "/usr/bin/time -f 'IMPORT elapsed=%e maxrss_kb=%M exit=%x' timeout 120s env PYTHONDONTWRITEBYTECODE=1 python -B - <<'PY'
import importlib.util
from pathlib import Path
p = Path('research/pde_ledger_v3/scripts/S11c_c2_exports.py')
spec = importlib.util.spec_from_file_location('s11cc2_exports_verify', p)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ledger = module.LEDGER
roots = ('s11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel')
print('IMPORT_OK', True)
print('ROWS', len(ledger))
print('SELF_ENERGY_PRESENT', 's11cc2SelfEnergyIncrement' in ledger)
for root in roots:
    print('ROOT', root, 'PRESENT', root in ledger)
    value = ledger[root]['value']
    print('OUTER_TYPE', type(value).__name__, 'CASE_COUNT', len(value))
    for case, payload in value:
        print('CASE', tuple(map(str, case)), 'PAYLOAD_LABELS', tuple(str(label) for label, _ in payload))
PY
git status --short" in /var/projects/toy_physics
 succeeded in 10799ms:
IMPORT_OK True
ROWS 70
SELF_ENERGY_PRESENT False
ROOT s11cc2ClosedSlabOperator PRESENT True
OUTER_TYPE Tuple CASE_COUNT 4
CASE ('LAB_HELD', 'RHO4_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('LAB_HELD', 'RHOBR_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('MATERIAL_ADVECTED', 'RHO4_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('MATERIAL_ADVECTED', 'RHOBR_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
ROOT s11cc2ClosedCouplingKernel PRESENT True
OUTER_TYPE Tuple CASE_COUNT 4
CASE ('LAB_HELD', 'RHO4_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('LAB_HELD', 'RHOBR_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('MATERIAL_ADVECTED', 'RHO4_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
CASE ('MATERIAL_ADVECTED', 'RHOBR_CONSTANT') PAYLOAD_LABELS ('VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M', 'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS')
IMPORT elapsed=10.92 maxrss_kb=1622224 exit=0
?? research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md

exec
/bin/bash -lc "rg -n --max-columns 240 --max-columns-preview 's11cc2(SelfEnergyIncrement|ClosedSlabOperator|ClosedCouplingKernel)' research/pde_ledger_v3/scripts/S11c_c2_exports.py && rg -n '"'^## |''^### |§3c|§5e|N6|uniform|self-energy|increment|sign convention|sign-convention|face-generalized|#90|directional|adjoint|Wolfram|T7|σ_W|sigma_W'"' research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md | sed -n '1,260p' && nl -ba research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md | sed -n '100,260p' && nl -ba research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md | sed -n '470,650p'" in /var/projects/toy_physics
 succeeded in 0ms:
11:'s11cc2ClosedCouplingKernel': {'value': _restore("Tuple(Tuple(Tuple(Str('LAB_HELD'), Str('RHO4_CONSTANT')), Tuple(Tuple(Str('VALUE'), Tuple(Tuple(Str('TRANSVERSE_TO_THICKNESS'), Tuple(Tuple(Str('THETA'), Mul(Symbol('epsilon_shape'), Add(Mul [... omitted end of long line]
12:'s11cc2ClosedSlabOperator': {'value': _restore("Tuple(Tuple(Tuple(Str('LAB_HELD'), Str('RHO4_CONSTANT')), Tuple(Tuple(Str('VALUE'), Tuple(Tuple(Str('U'), Tuple(Mul(Symbol('epsilon_shape'), Add(Mul(Symbol('eta_bg', real=True), Symbol('sigma_ [... omitted end of long line]
1:# S11c-c2 — SHARED PHYSICS (the self-energy fold: closing the slab operator with the curved-bulk response)
10:coupled **nonlocal self-energy operator**: the transverse sector now carries a self-energy threaded through the bulk
15:declared `IMPORT_KEYS`; the Wolfram engine imports nothing and re-derives every consumed object from the sibling
27:## 0 · Scope
32:**closed** full operator by the S11c-b §3c weak variational restriction — the **close-then-extract** ordering (§2).
33:(3) Emit the resulting **nonlocal self-energy** as the **substitution increment** (§3c), both operands re-extracted
43:(`DEFERRED_HEAVY_RUNS.md`); c2 must be constructible and cross-engine-testable on this box for its own increment
48:## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build
50:Everything in §1 is an input. The substitutions, orderings, extractions, self-energy structure, and the six
53:### 1a · Inheritance and the two imported models
55:The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12`), and admissibility are exactly S11c-a
63:  content is per-engine leg-verified only, and two whole-row sign conventions + #90's two flags are
66:  `+diff` vs WL `−linearVirtualVariation`, and the **#90 closure-fold sign**. ⛔ The face-force and closure-fold
67:  conventions are the coefficients of the slots c2 substitutes into (§3c) — they do **not** cancel from c2's residual.
70:  bare S11b `face_response`/`face_response_coeffs`, which c1 imported as its uniform-limit regression operand);
76:### 1b · What is cross-engine AGREE vs UNDECIDED in the c1 import — SUPPLIED HONESTLY (rule 6/16)
93:### 1c · The S11c-b slab operator and the symbolic slots the fold closes — SUPPLIED (grounded in the real row)
100:closure_shape_deriv_s` (the #90 fold). ⛔ **The explicit relative-flux `J_s` carrier has already been eliminated
101:from the θ-row by the #90 subtraction** — what remains is the flux **closure** written out as
109:⛔ The face-force and #90 closure-fold sign conventions above are **not** normalized by c2; the §3c increment carries
112:### 1d · The face closure laws and the Λ-channel routing — SUPPLIED (this routing is c2's task)
135:## 2 · The three fold objects and the non-commuting ordering — SUPPLIED framing
138:§3c). Closing the face response threads the nonlocal bulk DtN operator `Z` through the operator, including into the
148:the closure introduces coupling the extract-first route discards, and that induced coupling is the self-energy.
150:⚠ **Three DISTINCT objects — each gets its own name (⛔ do not call all three "self-energy"):**
155:S11CC2_SELF_ENERGY_INCREMENT   = extract(close(SLAB)) − extract(SLAB)  (c2's OWN-ROWS object, §3c) .
162:## 3 · The self-energy construction (OUTPUTS)
165:`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component value, sign, order, parity, or grade
168:### 3a · Close the full operator
187:### 3b · Re-extract the off-diagonal coupling from the CLOSED operator
189:Apply the S11c-b §3c **weak variational restriction verbatim** (independent divergence-free transverse and
193:implement the split by zeroing only undifferentiated field occurrences (S11c-b §3c). Emit both blocks and, **only if
194:the pairing-based adjointness residual is a genuine independent route**, that residual; ⛔ it is not the mixed second
195:derivative of a scalar energy (rule 2 corollary 3) — if the two blocks are adjoint by construction, emit them and say
203:### 3c · The self-energy as the substitution increment — the export representation
205:Emit the self-energy as the **substitution increment**, both operands re-extracted from `slab_operator` **with the
214:substitution increment, supported only where the closed `δp_s` replaced the symbolic slot. ⚠ **This is an export
216:v1):** the increment **drops the S11c-b bulk/kinetic base** (every term with no `δp_s` slot — including the deferred
218:within each engine and cancels). It does **NOT** drop the **face-generalized-force and #90 closure-fold sign
219:conventions**: those are the coefficients multiplying the `δp_s` slots, so they multiply the increment and — being
220:cross-engine-UNVALIDATED — can leak into c2's cross-engine residual `(≈ 2·carrier·increment)`. ⇒ the comparator
222:operands and the increment; ⛔ do not claim the increment isolates c2 from those two conventions.
224:### 3d · The six re-adjudications (rule 17 / c1 UNDECIDED carry-ins the fold makes load-bearing)
232:   `RHOBR_CONSTANT` the frozen `ρ_br`. Emit the self-energy for both and the field-vs-field object
237:   accepted corruption (`N6`).
244:   ablate it — **which self-energy terms depend on the raw whole-form `dtn_operator`** beyond the kernel
251:   This is the control that **settles the load-bearing face-generalized-force sign** the increment (§3c) does not
255:   `dtn_flat_symbol` **only** as the uniform-limit regression operand (diagonal `k=k′`), pin it so and say so; if any
256:   MATERIAL off-diagonal self-energy term consumes it, it is a **sixth UNDECIDED item to re-adjudicate**, ⛔ not
260:   emit the `μ_R,bg`-form ablation of the self-energy (§5) — the c1 reservation is discharged here.
264:## 4 · Objects to compute and emit
269:- The **re-extracted closed off-diagonal coupling kernel** (both blocks), adjointness residual only if independent — §3b.
270:- The **self-energy increment** and its two same-extract operands — §3c.
276:Every result carries its `(ε,η,σ_W)` order (`N12`; the transverse↔thickness coupling is the inherited `O(εη)`, the
277:self-energy the increment threaded through the bulk `Z`) and its restored `[L,T,M]` dimension. ⛔ No result is
282:## 5 · Independent routes and controls
288:### 5a · The ordering ablation (the non-commutation, §2)
292:`S11CC2_ORDERING_COMMUTATOR`. The commutator is the closure-induced self-energy and must be **nonzero** in general; a
296:### 5b · The routing ablation (§1d)
300:self-energy increment (nonzero residual vs §3c). ⛔ Tests that the §1d routing is load-bearing, ⛔ not that any
303:### 5c · The N6 independent-route control — the two anchorings ARE the representation-invariance pair
305:`N6` requires **two independently-constructed routes** compared after field redefinition, then one-sided corruption
307:MATERIAL_ADVECTED): construct the self-energy increment in each, map Eulerian↔material by the S11c-a/S11b field
310:omit its advective `u·∇Σ_E⁰` contribution — and require the increment to move while the other route's contribution
314:### 5d · The background-density field-vs-field re-adjudication (rule 17, §3d.1)
316:Emit the self-energy increment for both density representatives (`RHO4_CONSTANT` live, `RHOBR_CONSTANT` frozen) and
317:their difference `S11CC2_DENSITY_LIVE_MINUS_FROZEN`. ⛔ `∇ρ→0` is not an accepted corruption (`N6`); the genuine test
321:### 5e · Three DISTINCT reduction limits (⛔ not conflated) + the μ_R,bg form ablation
326:- **Uniform limit** (`W_bg→W̄₀`, `η→0`): the off-diagonal self-energy increment must vanish (S11b decoupling, `N6`);
328:- **Zero-DtN** (`Z→0`): the increment's **bulk-nonlocal** (Z-dependent) part vanishes; ⛔ do not equate this with
331:  self-energy to move (⛔ a coefficient rescale is insufficient — only a form change tests the coupling).
335:## 6 · Method, dimensions, and script obligations
339:  weak restriction is S11c-b §3c verbatim; the operator inverse is c1 §3a/§3b verbatim.
341:  ([[feedback_dimensional_consistency_check]]); `(ε,η,σ_W)` multigrade on every object (`N12`).
348:  ACTION/imported operands. ⛔ No tautological residual (rule 2 corollary 3): the §3c increment is an export
349:  representation, ⛔ not a check; the §3b adjointness and §3d.4 pairing residuals are emitted only when a genuine
352:  nonlocal `Z`, so the full increment may be heavy — measure the process that runs, defer heavy controls
358:## 7 · Names, F9 reservations, chain output, and export schema
377:  `tau_A`, `tau_V`, `tau_X`, `rho_m`, `rho_br`, `W_0`, `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`, `mu_R`,
379:  §5a (`export_ledger_bind_closure_design.md:148-153`) — it is ⛔ **not** a construction operand for §3c.
394:**The comparator (`N8`, frozen `T7` contract).** The c2 comparator joins the two blind engines' emitted objects by
396:load-bearing residual is on the **self-energy increment** (§3c) — which drops the deferred ≥64 GB S11c-b bulk residual
397:but **carries the face-force/#90 closure-fold sign conventions** (§3c); the comparator **SURFACES** those (rule 1/6),
405:**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
411:## 8 · Supplied versus computed; builder report
415:non-commuting-ordering framing and its counterexample and the three named objects, the S11c-b §3c weak-restriction
419:kernel and its adjointness disposition (§3b); the self-energy increment and its operands (§3c); the six
420:re-adjudication objects (§3d); every control residual (§5); the `(ε,η,σ_W)` orders and `[L,T,M]` dimensions.
423:declares the fold's symbol map (§3a), routing (§1d), ordering (§2), and increment representation (§3c) it implemented;
   100	closure_shape_deriv_s` (the #90 fold). ⛔ **The explicit relative-flux `J_s` carrier has already been eliminated
   101	from the θ-row by the #90 subtraction** — what remains is the flux **closure** written out as
   102	`Λ_A𝒜_s + Λ_V V_s` with the face pressure `δp_s`/`d_w_delta_p_s` **symbolic**. In the real row (verified against
   103	`closure_shape_deriv`), the θ-row carries `−4I·Λ_A(−δp_s/ρ_m + μ_θ/ρ_br)/(ωτ_A+I)`, `−2I·Λ_V W_0 e_{W,t}/(ωτ_V+I)`,
   104	and the pressure-jet terms in `d_w_delta_p_s`. ⇒ **the fold operation is "substitute the closed `δp_s(V_s,μ_θ)` and
   105	its w-jets `d_w_delta_p_s` into the symbolic `delta_p_±`/`d_w_delta_p_±` slots"**, ⛔ **NOT "substitute a closed
   106	`J_s`" (there is no `J_s` slot; adding one would double-count the already-folded flux closure).**
   107	
   108	`mu_theta_operator` (`S11CB_MU_THETA_OPERATOR`) is the separate held-fixed constitutive operand (neither row).
   109	⛔ The face-force and #90 closure-fold sign conventions above are **not** normalized by c2; the §3c increment carries
   110	each engine's own convention, and the §3d.4 mechanical-power pairing adjudicates the face-force sign (rule 1/6).
   111	
   112	### 1d · The face closure laws and the Λ-channel routing — SUPPLIED (this routing is c2's task)
   113	
   114	Per anchoring `α` and face `s`, on the curved-face objects of the S11c-a substrate (⛔ not flat Cartesian objects):
   115	
   116	```text
   117	𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br,bg⁰ ,   μ_θ = S11CB_MU_THETA_OPERATOR (held-fixed operand) ,
   118	J_s = Λ_A(ω) 𝒜_s + Λ_V(ω) V_s ,               ⇐ the FLUX closure carries only Λ_A, Λ_V ,
   119	t_s = −(δp_s + Λ_X(ω) 𝒜_s) n̂_s ,              ⇐ the TRACTION carries Λ_X (a reciprocal-traction channel) ,
   120	n̂_s·v_bulk,s = V_s + J_s/ρ_m ,                interfacial mass balance (kinematics, not a result) ,
   121	Λ_I(ω) = Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,  the three τ_I independent .
   122	```
   123	
   124	⚠ **`Λ_X` appears ONLY in the traction `t_s`** — ⛔ **not** in the flux closure `J_s`, and ⛔ **not** in the S11c-a
   125	T-i shape derivative `closure_shape_deriv` (which is the shape derivative of `J_s − Λ_A𝒜_s − Λ_V V_s = 0` only; c1
   126	§1d:157-160). The mechanical routing of `t_s` (which carries `Λ_X` and `δp_s`) is c2's, via the S11c-a traction
   127	`traction` (T-d) / virtual work. c1 supplied the closed `(δp_s, J_s, t_s)(V_s, μ_θ)` by solving these with the bulk
   128	relation `δp_s = Z·v_bulk,s` (the DtN operator), keeping `Λ_I(ω)`/`μ_θ` symbolic. **The routing is c2's** (c1
   129	§1d:160): the closed `δp_s` (+jets) closes the θ-row's already-folded flux terms and enters `t_s`; the closed `t_s`
   130	enters the mechanical rows. The equations of motion are obtained by S11b's balance-law/virtual-work method, ⛔ not by
   131	putting an irreversible response kernel in an ordinary action.
   132	
   133	---
   134	
   135	## 2 · The three fold objects and the non-commuting ordering — SUPPLIED framing
   136	
   137	`coupling_kernel` was obtained by **weak variational restriction of the OPEN (δp_s-symbolic) §3b operator** (S11c-b
   138	§3c). Closing the face response threads the nonlocal bulk DtN operator `Z` through the operator, including into the
   139	transverse diagonal block, so:
   140	
   141	```text
   142	extract( close( SLAB ) )  ≠  close( extract( SLAB ) )      (extract and close do NOT commute) .
   143	```
   144	
   145	⭐ **The correct object is `extract(close(·))` — close the FULL operator first, then re-extract.** The counterexample
   146	fixing the ordering: for `R_x = x + p` with the closure `p = α y`, extracting the `x`-block first (`∂R_x/∂x = 1`, `p`
   147	dropped) then closing gives a zero `x→y` coupling, whereas closing first (`R_x = x + α y`) then extracting gives `α`;
   148	the closure introduces coupling the extract-first route discards, and that induced coupling is the self-energy.
   149	
   150	⚠ **Three DISTINCT objects — each gets its own name (⛔ do not call all three "self-energy"):**
   151	
   152	```text
   153	S11CC2_CLOSED_COUPLING_KERNEL  = extract(close(SLAB))                  (the re-extracted closed off-diagonal block) ,
   154	S11CC2_ORDERING_COMMUTATOR     = extract(close(SLAB)) − close(extract(SLAB))   (the §5a form-ablation control) ,
   155	S11CC2_SELF_ENERGY_INCREMENT   = extract(close(SLAB)) − extract(SLAB)  (c2's OWN-ROWS object, §3c) .
   156	```
   157	
   158	⛔ `close(extract(SLAB))` (close the already-extracted OPEN kernel) is the §5a ablation, ⛔ not a construction route.
   159	
   160	---
   161	
   162	## 3 · The self-energy construction (OUTPUTS)
   163	
   164	Every object below is computed for both anchorings `α` and both density representatives `ρ`; it carries its computed
   165	`(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component value, sign, order, parity, or grade
   166	in this document.
   167	
   168	### 3a · Close the full operator
   169	
   170	Substitute the c1 closed face pressure `δp_s(V_s,μ_θ)` and its w-jets `d_w_delta_p_s` (from `s11c_c1_face_response`)
   171	into the symbolic `delta_p_±`/`d_w_delta_p_±` slots of `slab_operator` — the θ-row's already-folded flux terms
   172	(`closure_shape_deriv`) and the mechanical-row traction (`traction`, carrying `Λ_X`). ⛔ Do **not** add a closed `J_s`
   173	term (the θ-row's flux is already `Λ_A𝒜_s + Λ_V V_s`; the fold closes its symbolic `δp_s`, §1c). Because
   174	`δp_s = Z·v_bulk,s` with `Z` a nonlocal two-momentum operator, the closed rows carry `Z` — the response elimination
   175	is the c1 operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹`, ⛔ not a scalar division. Map c1's `V_s`, opaque `μ_θ`, and
   176	per-case resolvent symbols onto the slab-row fields explicitly (the build directive freezes the exact symbol map).
   177	Retain the curved-face measure/normal and the live background coefficients; ⛔ do not freeze any coefficient at its
   178	constant binding before the fold (rule 17, §3d.1). **Emit the closed operator as the assembled two-face object per
   179	`(α,ρ)`** (both faces summed, as the slab EOM is), the per-face contributions as provenance, and the face-parity
   180	blocks.
   181	
   182	```text
   183	⇒ S11CC2_CLOSED_SLAB_OPERATOR (assembled, per (α,ρ)) , S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS (per face) ,
   184	  S11CC2_CLOSED_SLAB_OPERATOR_PARITY_BLOCKS .
   185	```
   186	
   187	### 3b · Re-extract the off-diagonal coupling from the CLOSED operator
   188	
   189	Apply the S11c-b §3c **weak variational restriction verbatim** (independent divergence-free transverse and
   190	curl-free longitudinal trial/test displacements, independent `θ`/`e_W` trial/test fields, compact support in the
   191	in-plane interior so the IBP boundary term is fixed to zero; the inherited face conditions still apply) to the
   192	**closed** full operator of §3a, extracting both off-diagonal blocks → `S11CC2_CLOSED_COUPLING_KERNEL`. ⛔ Do not
   193	implement the split by zeroing only undifferentiated field occurrences (S11c-b §3c). Emit both blocks and, **only if
   194	the pairing-based adjointness residual is a genuine independent route**, that residual; ⛔ it is not the mixed second
   195	derivative of a scalar energy (rule 2 corollary 3) — if the two blocks are adjoint by construction, emit them and say
   196	there is no independent second route rather than dressing a structural zero as a check.
   197	
   198	```text
   199	⇒ S11CC2_CLOSED_COUPLING_KERNEL , S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS ,
   200	  S11CC2_SELF_ENERGY_ADJOINTNESS_RESIDUAL  (emit only if an independent route exists) .
   201	```
   202	
   203	### 3c · The self-energy as the substitution increment — the export representation
   204	
   205	Emit the self-energy as the **substitution increment**, both operands re-extracted from `slab_operator` **with the
   206	same weak restriction** (⛔ NOT the imported open `coupling_kernel`, which the two-leg-gated
   207	`directives/export_ledger_bind_closure_design.md:148-153` bars as a construction operand — c2 re-extracts):
   208	
   209	```text
   210	S11CC2_SELF_ENERGY_INCREMENT  =  extract( close(SLAB) )  −  extract( SLAB )    (per α, ρ) .
   211	```
   212	
   213	Because the weak restriction is linear, this equals `extract(close(SLAB) − SLAB)` — the restriction of the pure
   214	substitution increment, supported only where the closed `δp_s` replaced the symbolic slot. ⚠ **This is an export
   215	representation, ⛔ not a residual check** (rule 2 corollary 3). ⭐ **What it does and does NOT isolate (corrected from
   216	v1):** the increment **drops the S11c-b bulk/kinetic base** (every term with no `δp_s` slot — including the deferred
   217	≥64 GB slab-operator content and the kinetic `−K/+K` convention — appears identically in `close(SLAB)` and `SLAB`
   218	within each engine and cancels). It does **NOT** drop the **face-generalized-force and #90 closure-fold sign
   219	conventions**: those are the coefficients multiplying the `δp_s` slots, so they multiply the increment and — being
   220	cross-engine-UNVALIDATED — can leak into c2's cross-engine residual `(≈ 2·carrier·increment)`. ⇒ the comparator
   221	**SURFACES** them (rule 1/6, §7) and the §3d.4 mechanical-power pairing **adjudicates the face-force sign**. Emit both
   222	operands and the increment; ⛔ do not claim the increment isolates c2 from those two conventions.
   223	
   224	### 3d · The six re-adjudications (rule 17 / c1 UNDECIDED carry-ins the fold makes load-bearing)
   225	
   226	Each is a **computation whose cross-engine disposition is a finding**, ⛔ not a value to match:
   227	
   228	1. **Background density, field-vs-field (rule 17, c1 seal 5).** Before the fold, bind `rho_br_bg_rho4_constant` to
   229	   its live relation from `background_density_map` (`rho_br_bg_rho4_constant = W_bg·ρ_br/W_0`, `W_bg=W_0(1+η w₁)`),
   230	   ⛔ not to a bare constant — else the fold silently repeats the c1 PY freeze while passing the manifest guards. The
   231	   **two density representatives are the field-vs-constant pair**: `RHO4_CONSTANT` carries the live `ρ_br,bg⁰(x)`,
   232	   `RHOBR_CONSTANT` the frozen `ρ_br`. Emit the self-energy for both and the field-vs-field object
   233	   `S11CC2_DENSITY_LIVE_MINUS_FROZEN = SelfEnergy[ρ_br,bg⁰(x)] − SelfEnergy[ρ_br]`. ⛔ This is a live-field difference,
   234	   ⛔ NOT a `∇ρ` jet (neither c1 engine differentiates ρ; at first shape order it is an O(η) multiplicative
   235	   difference) and ⛔ NOT `ρ_m` (the bulk fluid density, a different object). The comparator joins the two engines'
   236	   live-field objects **field-against-field**; ⛔ never a PY constant folded onto a WL field. ⛔ `∇ρ→0` is not an
   237	   accepted corruption (`N6`).
   238	2. **The `t_s` traction representation.** Carry `t_s` in its **native covector form** `t_s = −(δp_s + Λ_X𝒜_s) n̂_s`
   239	   (a covector along the curved outward normal), ⛔ not a pre-collapsed scalar, so the mechanical-row contribution is
   240	   the genuine covector pairing; the scalar-vs-4-vector difference c1 left UNDECIDED is carried into the comparator as
   241	   a representation to reconcile, ⛔ not silently collapsed on one side. ⇒ `S11CC2_TRACTION_MECHANICAL_CONTRIB`.
   242	3. **The DtN kernel vs whole-form.** Structure the fold so the load-bearing bulk object is the AGREE'd two-momentum
   243	   **kernel** `dtn_kernel` (its `q_out(k),q_out(k′)` legs + on-shell dispersion). Emit — separately, so a leg can
   244	   ablate it — **which self-energy terms depend on the raw whole-form `dtn_operator`** beyond the kernel
   245	   (`S11CC2_DTN_WHOLEFORM_DEPENDENCE`). ⛔ Do not assume whole-form AGREE; ⛔ do not fold the whole-form into the
   246	   kernel.
   247	4. **The traction-vs-slab mechanical-power pairing (c1-assigned to c2).** c1's energy audit "uses no slab EOM row and
   248	   does not import `S11CB_SLAB_OPERATOR` (that pairing is S11c-c2's, after the fold)" (c1 §3b:328-330); parent N2
   249	   gives c2 the folded operator and S11c-e only the conversion FORM. Emit the post-fold mechanical-power / traction
   250	   pairing against the slab kinetic-and-stored variation: `S11CC2_TRACTION_SLAB_POWER_PAIRING` and its residual. ⭐
   251	   This is the control that **settles the load-bearing face-generalized-force sign** the increment (§3c) does not
   252	   cancel — a one-sided `t_s`-sign corruption must move the mechanical-power residual against the slab kinetic term.
   253	5. **The flat-resolvent leg-labeling.** c1 left the PY-output-leg vs WL-input-leg convention UNDECIDED, equal on the
   254	   `k=k′` diagonal, differing off-diagonal in the MATERIAL_ADVECTED response coefficients. If the fold uses
   255	   `dtn_flat_symbol` **only** as the uniform-limit regression operand (diagonal `k=k′`), pin it so and say so; if any
   256	   MATERIAL off-diagonal self-energy term consumes it, it is a **sixth UNDECIDED item to re-adjudicate**, ⛔ not
   257	   inherited as AGREE. Emit `S11CC2_FLAT_SYMBOL_USAGE` naming where and how it enters.
   258	6. **The `μ_R,bg` form control (c1-reserved for c2).** c1 reserved the `μ_R,bg` form ablation for c2, "where `μ_θ` is
   259	   composed with the slab variables" (c1 §5a). Since the fold composes the opaque `μ_θ` operand into the slab rows,
   260	   emit the `μ_R,bg`-form ablation of the self-energy (§5) — the c1 reservation is discharged here.

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md | sed -n '260,435p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.out" in /var/projects/toy_physics
 succeeded in 0ms:
   260	   emit the `μ_R,bg`-form ablation of the self-energy (§5) — the c1 reservation is discharged here.
   261	
   262	---
   263	
   264	## 4 · Objects to compute and emit
   265	
   266	Per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, multigraded and dimensioned:
   267	
   268	- The **closed full slab operator** (assembled two-face, per-face provenance, parity blocks) — §3a.
   269	- The **re-extracted closed off-diagonal coupling kernel** (both blocks), adjointness residual only if independent — §3b.
   270	- The **self-energy increment** and its two same-extract operands — §3c.
   271	- The **six re-adjudication objects** — §3d: `S11CC2_DENSITY_LIVE_MINUS_FROZEN`, `S11CC2_TRACTION_MECHANICAL_CONTRIB`,
   272	  `S11CC2_DTN_WHOLEFORM_DEPENDENCE`, `S11CC2_TRACTION_SLAB_POWER_PAIRING`, `S11CC2_FLAT_SYMBOL_USAGE`, and the
   273	  `μ_R,bg`-form ablation output.
   274	- The **control outputs** of §5, each emitted as the object and its literal residual.
   275	
   276	Every result carries its `(ε,η,σ_W)` order (`N12`; the transverse↔thickness coupling is the inherited `O(εη)`, the
   277	self-energy the increment threaded through the bulk `Z`) and its restored `[L,T,M]` dimension. ⛔ No result is
   278	reported without both.
   279	
   280	---
   281	
   282	## 5 · Independent routes and controls
   283	
   284	⭐ Every control re-enters the chain **at the ACTION / the imported operands**, ⛔ never at a result. Each emits the
   285	object and its literal residual; ⛔ none asserts a target value. A **coefficient** rescale tests arithmetic; only a
   286	**form** change tests physics.
   287	
   288	### 5a · The ordering ablation (the non-commutation, §2)
   289	
   290	Emit `S11CC2_CLOSED_COUPLING_KERNEL = extract(close(·))` AND `close(extract(·))` (close the already-extracted OPEN
   291	`coupling_kernel` — a legitimate **regression** use of the imported open kernel) and their difference
   292	`S11CC2_ORDERING_COMMUTATOR`. The commutator is the closure-induced self-energy and must be **nonzero** in general; a
   293	byte-identical result is the tell that the closure was not threaded through the full operator. ⛔ Report the literal
   294	difference; ⛔ state no expected form.
   295	
   296	### 5b · The routing ablation (§1d)
   297	
   298	Re-run the fold with the routing corrupted one direction at a time: (i) route the closed `δp_s` only into the
   299	mechanical rows, dropping it from the θ-row flux terms; (ii) drop `Λ_X` from `t_s`. Each corruption must move the
   300	self-energy increment (nonzero residual vs §3c). ⛔ Tests that the §1d routing is load-bearing, ⛔ not that any
   301	channel is nonzero.
   302	
   303	### 5c · The N6 independent-route control — the two anchorings ARE the representation-invariance pair
   304	
   305	`N6` requires **two independently-constructed routes** compared after field redefinition, then one-sided corruption
   306	— ⛔ not a bare one-sided corruption. The two routes for c2 are the **two anchorings** `α∈{L,M}` (LAB_HELD /
   307	MATERIAL_ADVECTED): construct the self-energy increment in each, map Eulerian↔material by the S11c-a/S11b field
   308	redefinition (`Δρ = δρ_E + u·∇ρ⁰`), and emit their residual `S11CC2_REP_INVARIANCE_RESIDUAL` — it must vanish (the
   309	same operator in two representations). **Then** corrupt **one route only** — flip one face's slope term in `n̂_s`, or
   310	omit its advective `u·∇Σ_E⁰` contribution — and require the increment to move while the other route's contribution
   311	does not (the independence test between the ≥2 same-order channels, `N3`/`N4`). ⛔ A zero residual proves nothing
   312	unless each route is shown to move alone.
   313	
   314	### 5d · The background-density field-vs-field re-adjudication (rule 17, §3d.1)
   315	
   316	Emit the self-energy increment for both density representatives (`RHO4_CONSTANT` live, `RHOBR_CONSTANT` frozen) and
   317	their difference `S11CC2_DENSITY_LIVE_MINUS_FROZEN`. ⛔ `∇ρ→0` is not an accepted corruption (`N6`); the genuine test
   318	is that the live-field representative differs from the frozen one wherever `ρ_br,bg⁰(x)` enters, and that the two
   319	engines join field-against-field. ⛔ Do not manufacture a `∇ρ` jet; ⛔ do not substitute `ρ_m`.
   320	
   321	### 5e · Three DISTINCT reduction limits (⛔ not conflated) + the μ_R,bg form ablation
   322	
   323	⛔ `Z→0`, `Λ_A⁰=0`, and impermeability are **three different limits** (from the c1 face solve
   324	`δp_s=Z·v_bulk,s`, `J_s=Λ_A𝒜_s+Λ_V V_s`): `Z→0` ⇒ `δp_s=0`; `Λ_A⁰=0` ⇒ `δp_s=Z V(Λ_V+ρ)/ρ ≠ 0`; **impermeable**
   325	requires `Λ_A⁰=Λ_V⁰=0` (S11b step `:89-90`). Emit each as its own regression with its literal residual:
   326	- **Uniform limit** (`W_bg→W̄₀`, `η→0`): the off-diagonal self-energy increment must vanish (S11b decoupling, `N6`);
   327	  a **secondary** smoke test only (cannot see coefficient/sign/parity).
   328	- **Zero-DtN** (`Z→0`): the increment's **bulk-nonlocal** (Z-dependent) part vanishes; ⛔ do not equate this with
   329	  `Λ_A⁰=0` or impermeability.
   330	- **μ_R,bg form ablation** (§3d.6): perturb the FORM of `μ_R,bg(x)` in the composed `μ_θ` operand and require the
   331	  self-energy to move (⛔ a coefficient rescale is insufficient — only a form change tests the coupling).
   332	
   333	---
   334	
   335	## 6 · Method, dimensions, and script obligations
   336	
   337	- **Method.** Balance laws + the binding material virtual-displacement rule + variational derivatives with held-fixed
   338	  fields named + prescribed external virtual work (S11b), ⛔ never an irreversible kernel in an ordinary action. The
   339	  weak restriction is S11c-b §3c verbatim; the operator inverse is c1 §3a/§3b verbatim.
   340	- **Dimensions.** Restore `[L,T,M]` on every emitted object, dimensional consistency able-to-fail
   341	  ([[feedback_dimensional_consistency_check]]); `(ε,η,σ_W)` multigrade on every object (`N12`).
   342	- **Rest-frame limit.** Inherit `N11a` inert; c2 constructs **no** convective operator. Every c2 result inherits the
   343	  c1/S11b smallness domain (`|q_out·v_bulk_normal_0/ω|≪1` + boundary-layer/subsonic), ⛔ never aliasing
   344	  `v_bulk_normal_0` to `v_0`.
   345	- **Script obligations.** The three build-skill clauses bind the build directive (`.claude/skills/build/SKILL.md`):
   346	  a script PRINTS computed objects and never states conclusions; PRINT the residual, do not assert it; interpretation
   347	  is the step record. ⛔ No hand-typed CAS object standing in for a computed one; every control re-enters at the
   348	  ACTION/imported operands. ⛔ No tautological residual (rule 2 corollary 3): the §3c increment is an export
   349	  representation, ⛔ not a check; the §3b adjointness and §3d.4 pairing residuals are emitted only when a genuine
   350	  independent route exists.
   351	- **Serialize CAS jobs; watch RSS.** c1 measured LIGHT (comparator peak ~317 MB on 30 GB); c2's fold threads the
   352	  nonlocal `Z`, so the full increment may be heavy — measure the process that runs, defer heavy controls
   353	  in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`), ⛔ never two memory-heavy CAS jobs concurrently. Detached launch
   354	  (harness reaps `run_in_background` ~87 s). Mathematica: 2-seat licence, `--sandbox danger-full-access`.
   355	
   356	---
   357	
   358	## 7 · Names, F9 reservations, chain output, and export schema
   359	
   360	**F9 / `N14` reservations.** Every new object gets a **fresh** injective `mechanical_lower_camel` name; ⛔ never reuse
   361	an imported S11c-b/c1/S11b key for a new object. ⛔ Never reuse `slab_operator`, `coupling_kernel`,
   362	`closure_shape_deriv`, `dtn_operator`, `dtn_kernel`, `face_response`, or any imported constant key for a c2 object.
   363	
   364	**Chain output (`N1`/`N8`; topology = the two-leg-gated `directives/export_ledger_bind_closure_design.md` §D1–§D3).**
   365	The SymPy engine reads the inherited model via the **positional** call
   366	`load_model("scripts/S11c_b_exports.py", "scripts/S11c_c1_exports.py")` (signature `load_model(base_path,
   367	*delta_paths)`, `scripts/ledger_fold.py:102`; ⛔ NOT `load_model(base=…, deltas=[…])` — that TypeErrors), binding only
   368	its declared `IMPORT_KEYS`, and writes `scripts/S11c_c2_exports.py` as its **own-rows delta** (§D2, ⛔ not the
   369	accumulated whole-model file). Verified on the real files: base 2441 rows + c1 delta 44 rows → fold 2485, exact-key
   370	intersection empty, no overwrite, `check_consumer` closure resolves. What the step **binds** (§D1) is the §3
   371	consume-set from both parents, named by their **F9 write-keys**:
   372	
   373	- from **S11c-b**: `slab_operator`, `slab_operator_term_origins`, `mu_theta_operator`, `closure_shape_deriv`
   374	  (S11c-a-authored, in the fold), the S11c-a T-a..T-i face substrate (`face_normal`, `conormal_deriv`,
   375	  `face_measure_shape_deriv`, `face_velocity`, `relative_flux`, `kinematic_balance`, `traction`, `face_shift`),
   376	  `background_density_map`, and the constants/kernels those rows reference (`Lambda_A_0`, `Lambda_V_0`, `Lambda_X_0`,
   377	  `tau_A`, `tau_V`, `tau_X`, `rho_m`, `rho_br`, `W_0`, `W_bg`, `w1_profile`, `L_W`, `sigma_W`, `eta_bg`, `mu_R`,
   378	  `epsilon_shape`, …). ⛔ The open `coupling_kernel` is **emit-only unless c2 declares it a REGRESSION operand** for
   379	  §5a (`export_ledger_bind_closure_design.md:148-153`) — it is ⛔ **not** a construction operand for §3c.
   380	- from **S11c-c1**: `s11c_c1_face_response`, `s11c_c1_face_response_coeffs`, `dtn_operator`, `dtn_flat_symbol`,
   381	  `dtn_kernel`, the per-case `s11cc1_response_resolvent_*`, the momentum symbols `s11cc1_{k,q}_{in,out}put*`,
   382	  `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_*`. ⛔ **NOT the bare `face_response`/`face_response_coeffs`**
   383	  (those are the S11b open/flat rows c1 imported as its regression operand; binding them folds the wrong physics —
   384	  the N14/F9 false-equal).
   385	
   386	⛔ The exact `IMPORT_KEYS` **root set** (minimal roots whose recursive closure covers the consume-set) is fixed at the
   387	c2 build directive against the two real export files, ⛔ not enumerated-then-frozen here; its two decision legs verify
   388	it, and that the guard (`check_consumer`/`assert_lookups_equal_manifest`/`assert_delta_is_minimal`) passes on the
   389	two-parent fold — ⚠ noting the guard passes on **key existence**, so it will **not** catch a wrong-provenance
   390	`face_response` binding; that check is the directive's + legs' responsibility. `BUILD_INPUT_DIGESTS` pins, per §D3,
   391	`{this sub-step's SymPy audit, scripts/S11c_b_exports.py, scripts/S11c_c1_exports.py, this spec, scripts/ledger_fold.py}`.
   392	⛔ Never `git add -f` a big `.out`; ⛔ never annex an `*_exports.py`.
   393	
   394	**The comparator (`N8`, frozen `T7` contract).** The c2 comparator joins the two blind engines' emitted objects by
   395	name, pairs residual operands, is three-valued, rejects a native boolean, and PRINTS/decides nothing (rule 2). Its
   396	load-bearing residual is on the **self-energy increment** (§3c) — which drops the deferred ≥64 GB S11c-b bulk residual
   397	but **carries the face-force/#90 closure-fold sign conventions** (§3c); the comparator **SURFACES** those (rule 1/6),
   398	⛔ does not normalize them, and the §3d.4 mechanical-power pairing adjudicates the face-force sign. It also surfaces
   399	the §3d representation questions (density field-vs-field, `t_s` scalar-vs-4-vector, DtN whole-form-vs-kernel,
   400	flat-symbol usage); the reconcile is the staged representational bridge
   401	([[feedback_reconcile_representational_bridge]]), ⛔ never a blanket collapse. The four c1 giant families + the full
   402	per-family symbolic residual remain deferred (`DEFERRED_HEAVY_RUNS.md`); c2 names, ⛔ does not pre-adjudicate,
   403	whatever it cannot close on this box.
   404	
   405	**The blind Wolfram engine** re-derives the §§1–2 supplied inputs, the S11c-a face substrate, the S11c-b slab-operator
   406	rows it folds into, and the c1 closed response — importing nothing (the only cross-engine control). ⛔ The denylist
   407	stays cut (`N9`/rule 12); blindness is enforced by absence.
   408	
   409	---
   410	
   411	## 8 · Supplied versus computed; builder report
   412	
   413	**SUPPLIED (unfalsifiable in this build):** all of §1 (the two imported models and their AGREE/UNDECIDED disposition,
   414	the slab-operator rows and the real θ-row structure, the face closure laws and Λ-channel routing), the §2
   415	non-commuting-ordering framing and its counterexample and the three named objects, the S11c-b §3c weak-restriction
   416	convention, the c1 §3a/§3b operator inverse, `background_density_map`, `N11a`, `N12`.
   417	
   418	**COMPUTED (outputs, ⛔ none stated here):** the closed two-face operator (§3a); the re-extracted closed coupling
   419	kernel and its adjointness disposition (§3b); the self-energy increment and its operands (§3c); the six
   420	re-adjudication objects (§3d); every control residual (§5); the `(ε,η,σ_W)` orders and `[L,T,M]` dimensions.
   421	
   422	**Builder report.** The build directive states, per emitted object, which line computed it (`.claude/skills/build`);
   423	declares the fold's symbol map (§3a), routing (§1d), ordering (§2), and increment representation (§3c) it implemented;
   424	and reports the literal residuals of §5 — ⛔ never a prose conclusion. The disposition of the six §3d re-adjudications
   425	and every §5 residual is read on **our** side, in the step record, ⛔ not asserted by the script (rule 5).
     1	==== PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT
     2	   [TRANSVERSE_TO_THICKNESS/THETA] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
     3	        SURVIVING(closed_part) head: -2*I*Lambda_A_0*epsilon_shape*s11cc2TestTheta(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1,
     4	   [TRANSVERSE_TO_THICKNESS/E_W] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
     5	        SURVIVING(closed_part) head: 2*I*Lambda_X_0*W_0*epsilon_shape*s11cc2TestE(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1, 
     6	   [TRANSVERSE_TO_THICKNESS/DIV_U] block ZERO
     7	   [THICKNESS_TO_TRANSVERSE/THETA] block ZERO
     8	   [THICKNESS_TO_TRANSVERSE/E_W] block ZERO
     9	   [THICKNESS_TO_TRANSVERSE/DIV_U] block ZERO
    10	==== PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT
    11	   [TRANSVERSE_TO_THICKNESS/THETA] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    12	        SURVIVING(closed_part) head: -2*I*Lambda_A_0*epsilon_shape*s11cc2TestTheta(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1,
    13	   [TRANSVERSE_TO_THICKNESS/E_W] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    14	        SURVIVING(closed_part) head: 2*I*Lambda_X_0*W_0*epsilon_shape*s11cc2TestE(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1, 
    15	   [TRANSVERSE_TO_THICKNESS/DIV_U] block ZERO
    16	   [THICKNESS_TO_TRANSVERSE/THETA] block ZERO
    17	   [THICKNESS_TO_TRANSVERSE/E_W] block ZERO
    18	   [THICKNESS_TO_TRANSVERSE/DIV_U] block ZERO
    19	==== PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT
    20	   [TRANSVERSE_TO_THICKNESS/THETA] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    21	        SURVIVING(closed_part) head: -2*I*Lambda_A_0*epsilon_shape*s11cc2TestTheta(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1,
    22	   [TRANSVERSE_TO_THICKNESS/E_W] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    23	        SURVIVING(closed_part) head: 2*I*Lambda_X_0*W_0*epsilon_shape*s11cc2TestE(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1, 
    24	   [TRANSVERSE_TO_THICKNESS/DIV_U] block ZERO
    25	   [THICKNESS_TO_TRANSVERSE/THETA] block ZERO
    26	   [THICKNESS_TO_TRANSVERSE/E_W] block ZERO
    27	   [THICKNESS_TO_TRANSVERSE/DIV_U] block ZERO
    28	==== PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT
    29	   [TRANSVERSE_TO_THICKNESS/THETA] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    30	        SURVIVING(closed_part) head: -2*I*Lambda_A_0*epsilon_shape*s11cc2TestTheta(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1,
    31	   [TRANSVERSE_TO_THICKNESS/E_W] full_block_zero=False bare_slots=['delta_p_minus', 'delta_p_plus'] n_Integral=1  closed_part_zero(expand)=False  closed_part_zero(doit)=True
    32	        SURVIVING(closed_part) head: 2*I*Lambda_X_0*W_0*epsilon_shape*s11cc2TestE(s11cc2X1, s11cc2X2, s11cc2X3, s11cc2Time)*Integral(0, (s11cc1_k_output_1, -oo, oo), (s11cc1_k_output_2, -oo, oo), (s11cc1_k_output_3, -oo, oo), (s11cc2Y1, 
    33	   [TRANSVERSE_TO_THICKNESS/DIV_U] block ZERO
    34	   [THICKNESS_TO_TRANSVERSE/THETA] block ZERO
    35	   [THICKNESS_TO_TRANSVERSE/E_W] block ZERO
    36	   [THICKNESS_TO_TRANSVERSE/DIV_U] block ZERO
     1	==== G: SELF_ENERGY_INCREMENT block directionality (LAB_HELD_RHO4) ====
     2	   THICKNESS_TO_TRANSVERSE/DIV_U    zero=True  bytes~10
     3	   THICKNESS_TO_TRANSVERSE/E_W      zero=True  bytes~10
     4	   THICKNESS_TO_TRANSVERSE/THETA    zero=True  bytes~10
     5	   TRANSVERSE_TO_THICKNESS/DIV_U    zero=True  bytes~10
     6	   TRANSVERSE_TO_THICKNESS/E_W      zero=False  bytes~1708669
     7	   TRANSVERSE_TO_THICKNESS/THETA    zero=False  bytes~855746
     8	==== G: CLOSED_COUPLING_KERNEL block directionality (LAB_HELD_RHO4) ====
     9	   THICKNESS_TO_TRANSVERSE/DIV_U    zero=False  bytes~77695
    10	   THICKNESS_TO_TRANSVERSE/E_W      zero=False  bytes~156361
    11	   THICKNESS_TO_TRANSVERSE/THETA    zero=False  bytes~83329
    12	   TRANSVERSE_TO_THICKNESS/DIV_U    zero=False  bytes~708125
    13	   TRANSVERSE_TO_THICKNESS/E_W      zero=False  bytes~2251947
    14	   TRANSVERSE_TO_THICKNESS/THETA    zero=False  bytes~1078305
    15	==== E: REP_INVARIANCE_RESIDUAL sigma_W-sector confinement (RHO4) ====
    16	    PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHO4_CONSTANT blocks=6 sigma_W_symbol=[sigma_W]
    17	      THICKNESS_TO_TRANSVERSE/DIV_U    full_zero=False  sigmaW->0_zero=True
    18	      THICKNESS_TO_TRANSVERSE/E_W      full_zero=False  sigmaW->0_zero=True
    19	      THICKNESS_TO_TRANSVERSE/THETA    full_zero=False  sigmaW->0_zero=True
    20	      TRANSVERSE_TO_THICKNESS/DIV_U    full_zero=False  sigmaW->0_zero=True
    21	      TRANSVERSE_TO_THICKNESS/E_W      full_zero=False  sigmaW->0_zero=True
    22	      TRANSVERSE_TO_THICKNESS/THETA    full_zero=False  sigmaW->0_zero=True
    23	    PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHOBR_CONSTANT blocks=6 sigma_W_symbol=[sigma_W]
    24	      THICKNESS_TO_TRANSVERSE/DIV_U    full_zero=False  sigmaW->0_zero=True
    25	      THICKNESS_TO_TRANSVERSE/E_W      full_zero=False  sigmaW->0_zero=True
    26	      THICKNESS_TO_TRANSVERSE/THETA    full_zero=False  sigmaW->0_zero=True
    27	      TRANSVERSE_TO_THICKNESS/DIV_U    full_zero=False  sigmaW->0_zero=True
    28	      TRANSVERSE_TO_THICKNESS/E_W      full_zero=False  sigmaW->0_zero=True
    29	      TRANSVERSE_TO_THICKNESS/THETA    full_zero=False  sigmaW->0_zero=True
    30	=== builder report: adjointness justification ===
    31	70:There is no independent adjointness construction in this implementation, so no `SELF_ENERGY_ADJOINTNESS_RESIDUAL` is emitted. Both off-diagonal blocks are emitted. The traction pairing compares the slab's inherited generalized-force contribution to independently composed native covector work; its kinetic/stored operand is the full slab mechanical power with the inherited face generalized-force rows removed. The test velocities use the imported velocity amplitudes, with the external wave epsilon divided out consistently with the pressure map; physical quadratic mechanical power restores that test-amplitude epsilon. This is a comparison of face work, not an independent derivation of all kinetic and stored energies. It does not adjudicate c1's separate far-field ENERGY audit.
    32	
    33	[exited with code 0]
     1	closed slab operator: 5 algebraic leaves
     2	FAITHFUL (compact==emitted):      5/5  fails=[]
     3	GUARD BITES x2 corruption:        5/5  slipped=[]
     4	GUARD BITES +marker corruption:   5/5  slipped=[]
     5	INTEGRAND corruption on /THETA: bites=True
     6	
     7	[exited with code 0]

codex
The generated export passes a real import (11 s, ~1.62 GB RSS), contains exactly four cases under each required root, and omits the increment. The physics evidence does support the narrow F/G/E facts—but I’m checking the project’s own review-discipline rules now, because “two re-review legs clear” may be stronger than what the export adjudication actually records: one leg is explicitly incomplete.
exec
/bin/bash -lc "rg --files -g 'AGENTS.md' -g 'CLAUDE.md' -g 'MEMORY.md' -g '*state*' | sed -n '1,120p' && rg -n 'review-until-clear|two legs|2 legs|leg|G1|G2|G3|rule 13|disagreement|wedg|incomplete|adjudicat' CLAUDE.md research/pde_ledger_v3/CLAUDE.md 2>/dev/null | sed -n '1,260p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_grok.txt | sed -n '1,220p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out | sed -n '1,120p' && find /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad -maxdepth 1 -type f -printf '%f %s bytes\\n' | sort" in /var/projects/toy_physics
 succeeded in 0ms:
notes/g2/g2_step38_adiabatic_wall_ground_state.md
scripts/g2/g2_step38_adiabatic_wall_ground_state.py
notes/5pn/5pn_stage212_214_even_gate_reinstatement.md
CLAUDE.md
software/stage1_solver/scripts/__pycache__/pathA_C0g_state_svd_ftau.cpython-310.pyc
software/stage1_solver/scripts/pathA_C0g_state_svd_ftau.py
memory/_meta/state.json
software/stage1_solver/_scratch/exp6_legit_restate.py
scripts/5pn/5pn_stage323_coherent_kernel_state_extractor.py
scripts/5pn/5pn_stage329_coherent_local_dn_branch_state.py
scripts/5pn/5pn_stage326_microscopic_coherent_placement_state.py
software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/state_measurements/c0g_state_tau_0p03.json
software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/state_measurements/c0g_state_tau_0p02925.json
software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/state_measurements/c0g_state_tau_0p0290625.json
software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/state_measurements/c0g_state_tau_0p029125.json
software/stage1_solver/runs/pathA_C0g_diag_fold_vs_conditioning/state_measurements/c0g_state_tau_0p0295.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p03.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p029124.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p02912.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p02925.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p029125.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p029122.json
software/stage1_solver/runs/pathA_C0g_build_B1B2/b2_reconfirm/state_measurements/c0g_state_tau_0p0291225.json
research/pde_ledger_v2/schemas/examples/reject_fixtures/report_register_unknown_state.yaml
research/pde_ledger_v3/directives/W0_restatement_build.md
research/pde_ledger_v3/directives/W0_restatement_verification_report.md
research/pde_ledger_v3/directives/W0_statements_review.md
research/pde_ledger_v3/directives/W0_builder_statement.md
research/pde_ledger_v3/directives/W0_statements_review_round2.md
research/pde_ledger_v3/directives/maccullagh_reference_state_search_prompt.md
research/pde_ledger_v3/directives/W0_integrator_statement.md
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_022_tau_0p0291108575439_tau_0p0291108575439.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_025_tau_0p0291132_tau_0p0291132.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_024_tau_0p0291138_tau_0p0291138.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_027_tau_0p029113_tau_0p029113.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_028_tau_0p0291131_tau_0p0291131.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_023_tau_0p0291126455688_tau_0p0291126455688.npz
software/stage1_solver/runs/pathA_C0g_B3_deeppoint/C1/states/stage1_deepcrawl_state_dir_attempt_026_tau_0p0291128_tau_0p0291128.npz
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p0291144335937.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p02911625.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p0291138.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p02912.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p029122.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p0291132.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p029118125.json
software/stage1_solver/runs/pathA_C0g_deepcrawl/diagnostics/state_measurements/c0g_state_tau_0p0291225.json
research/pde_ledger/redteam_adversarial/provenance/fit_stage152_corrected_mouth_state_pi_corr_t_corr__pi_corr.yaml
research/pde_ledger/redteam_adversarial/provenance/fit_stage093_module_forced_statement__y_q_cons.yaml
research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage152_corrected_mouth_state_pi_corr_t_corr__phase_b__pi_corr.md
research/pde_ledger/redteam_adversarial/tmp_prompts/fit_stage093_module_forced_statement__phase_b__y_q_cons.md
research/pde_ledger/redteam_adversarial/provenance/_synthesis/batch_09/fit_stage093_module_forced_statement__y_q_cons.yaml
research/pde_ledger/redteam_adversarial/provenance/_synthesis/batch_16/fit_stage152_corrected_mouth_state_pi_corr_t_corr__pi_corr.yaml
CLAUDE.md:7:prior version by two independent legs before commit).
CLAUDE.md:12:leg. Only the user scales work down. Nothing below — no grouping, no budget, no "administrative" label —
CLAUDE.md:20:   criterion (the expected-value acceptance test, only); construct independently; preserve disagreement. (M1, M2)
CLAUDE.md:23:3. **Pick two non-author legs.** Orchestrator-written → Codex + Grok; Codex-written → fresh Claude + Grok.
CLAUDE.md:24:   Launch on sight, before inspecting the result. (G1, G3)
CLAUDE.md:25:4. **Use the right exit.** *Decision list:* two legs, verify, fold **once**, go. *Physics spec / script /
CLAUDE.md:26:   record:* **review-until-clear** — iterate leg→fold→leg until nothing outstanding changes what is
CLAUDE.md:27:   computed or may be claimed. The filename confers no exception. (G2, G4)
CLAUDE.md:36:(G2) **never** caps review of a physics spec, even one embedded in a decision list. "Review-until-clear" =
CLAUDE.md:48:| **Pre-builder decision list** | Two; O (always orchestrator-written) | Check the requested decisions + their supporting evidence; ⛔ no fictional-script ablation | **One two-leg pass**, verify findings, fold **once** — ⛔ no iterate-to-green | No builder before both reports + fold. An unresolved issue that changes what is computed or may be claimed routes to the applicable spec/build gate (R10 scope); physics-bearing content also meets the spec row. No commit before both reports. |
CLAUDE.md:50:| **Script / physics-bearing build** | Two; O or Cx; normally Cx for `/build` | Independent derivation scripts + literal stdout; every load-bearing check ablated; **mandatory FORM ablation**; one-sided corruption for independence; emit-before-guard + output checks | **Review-until-clear**; change author when repairs breed defects (⛔ never a fourth fold) | Launch legs before inspecting results; both usable reports; preserve reviewed baseline before repair; accept only on substantive clearance. Serialize dual Mathematica ablations. |
CLAUDE.md:52:| **Other / claimed non-physics** | Two if physics-bearing (G1); if you claim it is not, record why | First record whether it changes computation, premises, checks, or claims — if so it is physics-bearing, route to a row above | ⛔ Do not infer a one-pass or zero-review exemption from the suffix or an "administrative" label | The two-report commit gate has no non-physics exception. |
CLAUDE.md:65:**M1 — Two engines exist so they can disagree, and the disagreement is the measurement.** *(was R1, R6)*
CLAUDE.md:66:Independent construction, not hidden answers. A disagreement is a **finding** — ⛔ never try to make
CLAUDE.md:68:a disagreement as a builder target to eliminate; it is a finding on the orchestrator's side.
CLAUDE.md:98:review-leg prompt I write says *"a prose derivation is worth nothing; show the script and its literal
CLAUDE.md:119:**G1 — Whatever writes does not review; pick the two legs by authorship.** *(was R7)*
CLAUDE.md:120:**Whatever writes does not review.** Two independent legs on anything **physics-bearing** — and a spec both
CLAUDE.md:122:**Orchestrator-written → Codex + Grok. Codex-written → fresh Claude agent + Grok.** Choose legs by
CLAUDE.md:125:is required (G4) but is **never** a leg. *(evidence: L-R7)*
CLAUDE.md:127:**G2 — The decision-list gate — and ONLY the decision list gets the one-pass exception.** *(was R7 trigger;
CLAUDE.md:129:**TRIGGER — no builder launches until its decision list has had two legs.** The decision list is
CLAUDE.md:137:**G3 — Launch legs on sight, before inspecting the result; no commit before both report.** *(was R8, R9)*
CLAUDE.md:138:Launch legs **on sight**, before I look at the result — a self-check discharges the felt need for an
CLAUDE.md:141:one leg the other's output. **No commit before both legs report.** The commit is the last step. Reviewing
CLAUDE.md:146:**Stop when nothing outstanding changes what is computed or what may be claimed** — ⛔ not when both legs are
CLAUDE.md:147:green. "A leg that finds nothing is weak evidence" is my prior; ⛔ put it in a leg's prompt and it becomes a
CLAUDE.md:149:myself** (legs have been wrong in both directions); obtain **both** reports before adjudicating or editing.
CLAUDE.md:169:`.claude/skills/review-legs/SKILL.md`. A root summary cannot replace their identical-prompt,
CLAUDE.md:190:  `grep`-by-line still resolves once content is present (the symlink is followed), so a script/directive/leg
CLAUDE.md:207:  and eight legs died on a question one `len(LEDGER)` answered; *"the imported action is in the LEDGER"*
CLAUDE.md:214:- **L-R7 (G1, G2) — six spec defects.** Measured 2026-08-09: six spec defects, each costing a build round
CLAUDE.md:215:  *plus* two legs, when two legs before the build would have caught them — three "level-above" misses, one
CLAUDE.md:219:- **L-R8 (G3) — the self-check trap.** A self-check discharges the felt need for an independent one, and it
CLAUDE.md:221:- **L-R10 (G4) — the quota mechanism.** "A leg that finds nothing is weak evidence" is my prior; put it in a
CLAUDE.md:222:  leg's prompt and it becomes a quota.
CLAUDE.md:226:- **L-R13 (G4) — legs wrong in both directions.** Legs have been wrong in both directions.
     1	I'll start by reading the review skill and the export-repair directive so the six checks run against the actual requirements, not a paraphrase.I'll inspect the export-repair diff, fold machinery, and current artifacts first, then run the six checks independently in `/tmp`.Diff is publication-only so far. Next I’ll inspect TRIAGE/run/publish, the export module structure, and fold consumers before any re-run.TRIAGE returns before `publish`, so I’ll set up isolated `/tmp` copies, inspect `emit`/`publish` wiring, then run baseline vs current one case at a time.The 476 MB log grep is still running. I’ll set up the `/tmp` sandbox in parallel and keep the construction-hash comparison moving.Construction functions match the baseline; `run()` only differs at the export map. I’ll kill the 476 MB grep, instrument TRIAGE copies, and start the current one-case re-run.I'll instrument the TRIAGE copies, start the current one-case run, and write the membership and export-line extractors in parallel.Construction is unchanged except the export map. Next I’ll run current TRIAGE, then membership and the independent semantic comparison.Current TRIAGE finished in 71s. I’ll start the baseline TRIAGE next and run the membership ablation in parallel — it does not load parent exports.Baseline and current TRIAGE outputs are byte-identical. Next I’ll fix the semantic script and run membership, then the independent compact-vs-emitted comparison.Membership is clean: increment absent, both operators present in a 4-case tree, closure minimal. Next is the independent compact-vs-emitted semantic check.The ablation subtracted whole Tuples; I’ll walk algebraic leaves instead so the corrupted compact is compared the same way the guard does.Transparency is clean. Next I’ll re-run `publish` itself on a compaction that adds a real term, so the script’s own semantic guard has to hard-stop.# Independent re-review — S11c-c2 EXPORT REPAIR
     2	
     3	**Verdict:** the publication-only repair is **correct to commit**. No confirmed defect in physics, membership, compaction, singular locus, stored form, or guard bite. One hygiene gap (missing builder report) does not change what is computed or exported.
     4	
     5	All work was in `/tmp/s11cc2_export_review/`. The working tree was not modified.
     6	
     7	---
     8	
     9	## 1. Physics unchanged — **no finding**
    10	
    11	**Diff scope.** Against `8f3a017f`, every construction `def` is byte-identical. The only body changes are `publish` and `run`; the only added def is `publication_compact`.
    12	
    13	```
    14	CHANGED_DEFS ['publish', 'run']
    15	added {'publication_compact'}
    16	N_IDENTICAL 36  N_CHANGED 2
    17	```
    18	
    19	`run()`’s unified diff is **only** the `export_key` map dropping `SELF_ENERGY_INCREMENT` (`S11c_c2_selfenergy_fold_sympy_audit.py:1051-1052`). `build_case` / `build_face` / `extract` / `kernel_bridge` / `retained_shape` / `emit` / `traction_pairing` / `control` / `grade_object` are unchanged.
    20	
    21	**Emitted physics.** TRIAGE (`S11CC2_PACKAGE=TRIAGE`, one case, three core objects) of `git show 8f3a017f:…` vs the current script:
    22	
    23	| | current | baseline |
    24	|---|---|---|
    25	| stdout bytes | 10 966 662 | 10 966 662 |
    26	| pickle bytes | 250 680 | 250 680 |
    27	| `CLOSED_SLAB_OPERATOR` tag | 4 013 128 | identical |
    28	| `CLOSED_COUPLING_KERNEL` tag | 4 366 525 | identical |
    29	| `SELF_ENERGY_INCREMENT` tag | 2 575 252 | identical |
    30	
    31	Stdout and pickle are **byte-identical**. Increment is still **emitted**.
    32	
    33	Scripts / stdout:
    34	- `/tmp/s11cc2_export_review/check1_run_diff.py` → `measurements/check1_run_diff.out`
    35	- `/tmp/s11cc2_export_review/measurements/check1_all_defs.out`
    36	- `/tmp/s11cc2_export_review/measurements/check1_triage_cmp.out`
    37	- TRIAGE: `measurements/triage_{current,baseline}.out`
    38	
    39	---
    40	
    41	## 2. Membership — **no finding**
    42	
    43	Independent `load_model` / `check_consumer` / `assert_delta_is_minimal` on the generated module (not the script’s `EXPORT_SEMANTIC` printer):
    44	
    45	| check | result |
    46	|---|---|
    47	| file size | **22 441 522 B (21.4 MB)** |
    48	| `s11cc2SelfEnergyIncrement` | **absent** |
    49	| `s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel` | present |
    50	| cases each | all 4 `(anchoring, density)` |
    51	| payload labels | `VALUE`, `MULTIGRADE`, `DIMENSION_L_T_M`, `COMPUTED_BRANCH_BINDINGS`, `FOURIER_PROFILE_BINDINGS` |
    52	| `display` | short root name (24 / 26 chars), not `sstr(value)` |
    53	| leaked increment / origins / parity / control keys | **none** |
    54	| `load_model(b, c1, this delta)` | binds both operators; 0 overwrites |
    55	| `check_consumer` closure | resolves; own closure = 70 rows |
    56	| `assert_delta_is_minimal` | **PASS**; written keys == own closure |
    57	| `BUILD_INPUT_DIGESTS` | matches `sha256sum` of the five pinned files |
    58	
    59	Script: `/tmp/s11cc2_export_review/check2_membership.py`  
    60	Stdout: `/tmp/s11cc2_export_review/measurements/check2_membership.out`
    61	
    62	---
    63	
    64	## 3. Semantic equivalence (recomputed) — **no finding**
    65	
    66	For `LAB_HELD/RHO4_CONSTANT`, compact `exports.py` VALUE vs TRIAGE-emitted VALUE, with **strict containers** (case keys, tuple arities, mapping keys, matrix shapes) and **Integral-aware** `expand(compact − emitted)` (integrands expanded under protected reciprocal/calculus atoms; remaining `Integral`s Dummy-protected). Pole identity **and** `(base, exp)` signature **and** `together().as_numer_denom()` denominator atoms compared per leaf.
    67	
    68	```
    69	SEMANTIC_PASS s11cc2ClosedSlabOperator
    70	SEMANTIC_PASS s11cc2ClosedCouplingKernel
    71	N_LEAVES_ZERO 11 N_LEAVES_NONZERO 0
    72	N_POLE_IDENTITY_OK 11 FAIL 0
    73	N_POLE_SIG_FAIL 0 N_DEN_ATOM_FAIL 0
    74	BOOLEAN_FAILURES []
    75	```
    76	
    77	11 algebraic leaves = slab `U[0,1,2]`, `THETA`, `E_W` (5) + kernel two directions × `{THETA, E_W, DIV_U}` (6). No dropped trailing channel. Compact sizes match the build’s `EXPORT_VALUE_BYTES` (slab 4 002 816 → 2 271 711; kernel 4 355 993 → 3 375 460).
    78	
    79	Script: `/tmp/s11cc2_export_review/check3_semantic.py`  
    80	Stdout: `/tmp/s11cc2_export_review/measurements/check3_semantic.out`
    81	
    82	---
    83	
    84	## 4. Transparency + no information loss — **no finding**
    85	
    86	Grep + post-import atom walk of `S11c_c2_exports.py`:
    87	
    88	```
    89	TEXT_Dummy_COUNT 0
    90	TEXT_UnevaluatedExpr_COUNT 0
    91	TEXT_cse_COUNT 0
    92	DUMMY_HITS [] N 0
    93	UNEVAL_HITS [] N 0
    94	NON_BASIC [] N 0
    95	```
    96	
    97	Stored values are `sp.Basic` (`Tuple` of evaluable expressions). `sp.diff` on four `VALUE` leaves succeeded (`U[0..2]` wrt `sigma_W`, `THETA` wrt `omega`). Short `display`. Pole sets unchanged on all 11 leaves (check 3). Temporary `Dummy`s in `publication_compact` are substituted back before serialization (`:841-842`).
    98	
    99	Script: `/tmp/s11cc2_export_review/check4_transparency.py`  
   100	Stdout: `/tmp/s11cc2_export_review/measurements/check4_transparency.out`
   101	
   102	---
   103	
   104	## 5. FORM ablation — guards **do bite** — **no finding**
   105	
   106	**(a) Compaction.** Copied the script to `/tmp`; made algebraic leaves `out + 1` after the pole-preserving choice (pole atom set unchanged; a real term added). Independent leafwise residual:
   107	
   108	```
   109	CLOSED_SLAB_OPERATOR HONEST_LEAVES_ZERO 5 NONZERO 0
   110	CLOSED_SLAB_OPERATOR CORRUPT_LEAVES_ZERO 0 NONZERO 5
   111	CORRUPT_SAMPLE ('VALUE/U[0]', 1, True)   # residual 1, poles unchanged
   112	GUARD_BITES CLOSED_SLAB_OPERATOR
   113	CLOSED_COUPLING_KERNEL ... CORRUPT_LEAVES_ZERO 0 NONZERO 6
   114	GUARD_BITES CLOSED_COUPLING_KERNEL
   115	```
   116	
   117	Then **re-ran `publish()`** on the corrupted copy (TRIAGE one case, replicated to four case keys). The script’s own check hard-stopped:
   118	
   119	```
   120	EXPORT_SEMANTIC .../THETA expanded_difference_is_zero = False
   121	PUBLISH_RAISED ValueError ('emitted/compact mismatch',
   122	  's11cc2ClosedCouplingKernel/LAB_HELD/RHO4_CONSTANT/VALUE/TRANSVERSE_TO_THICKNESS/THETA',
   123	  'expanded_difference_is_zero', False)
   124	```
   125	
   126	**(b) Membership.** Re-add increment to `EXPORT_ROOTS` only, and drop `ClosedSlabOperator` from `EXPORT_ROOTS`, then call `publish` with the two-operator object dict:
   127	
   128	```
   129	RESULT readd_increment ValueError publication root membership
   130	RESULT drop_slab       ValueError publication root membership
   131	MEMBERSHIP_READD_BITES True
   132	MEMBERSHIP_DROP_BITES True
   133	```
   134	
   135	The first-line guard (`publish:853-854`) fires before compact. That is the intended hard-stop; `assert_delta_is_minimal` is not reached because membership already failed.
   136	
   137	Scripts / stdout:
   138	- `/tmp/s11cc2_export_review/check5_compact_ablation.py` → `measurements/check5_compact_ablation.out`
   139	- `/tmp/s11cc2_export_review/check5_publish_ablation.py` → `measurements/check5_publish_ablation.out`
   140	- `/tmp/s11cc2_export_review/check5_membership_ablation.py` → `measurements/check5_membership_ablation.out`
   141	
   142	---
   143	
   144	## 6. No overstep — **no finding** (hygiene note below)
   145	
   146	Working-tree delta from this repair:
   147	
   148	```
   149	 M research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py
   150	?? research/pde_ledger_v3/scripts/S11c_c2_exports.py
   151	?? research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json
   152	?? research/pde_ledger_v3/_measurements/S11c_c2_sympy_progress.json
   153	?? research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_rereview_prompt.md
   154	```
   155	
   156	No self-review / derived-or-declared / output-checker / finalize artifacts. `guard_evidence.json` is the expected checkpoint. `progress.json` is the unchanged `build_case` harness (`:533`). The `_legs/` prompt is orchestrator-created for this review.
   157	
   158	---
   159	
   160	## CONFIRMED defects
   161	
   162	**None.** Physics unchanged; compact ≡ emitted on all 11 leaves; poles preserved; increment absent from the delta; both operators bindable/differentiable; both new guards hard-stop under one-sided corruption.
   163	
   164	## Unsettled (not physics)
   165	
   166	- **`S11c_c2_export_repair_builder_report.md` is absent.** The directive asked for it; it is not in the tree. That is a builder-report gap, not an export defect.
   167	- Membership is defined by `EXPORT_ROOTS`. Re-adding the increment to **both** `EXPORT_ROOTS` and the `export_key` map would export it and pass minimality (it would be a declared root). The new `set(objects) != EXPORT_ROOTS` catch is desync, not a hardcoded ban. The committed `EXPORT_ROOTS` content is the control; it is correct.
   168	
   169	---
   170	
   171	## Commit?
   172	
   173	**Yes — the export repair is correct to commit** (script + `S11c_c2_exports.py` + the expected `S11c_c2_sympy_guard_evidence.json` checkpoint). Optional: add the missing builder report and `progress.json` if you want the directive’s report file in the same commit; neither changes the exported physics.
     1	=== s11cc2ClosedSlabOperator ('LAB_HELD', 'RHO4_CONSTANT') ===
     2	  leaves compared: 5
     3	=== s11cc2ClosedCouplingKernel ('LAB_HELD', 'RHO4_CONSTANT') ===
     4	  leaves compared: 6
     5	
     6	===== SUMMARY =====
     7	total algebraic leaves compared: 11
     8	nonzero difference leaves: 2
     9	   NONZERO s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/THETA -> 2*_Dummy_34*I*Lambda_A_0*epsilon_shape*eta_bg*sigma_W/(omega*rho_m*tau_A + I*rho_m) - 2*_Dummy_44*I*Lambda_A_0*epsilon_shape*eta_bg*sigma_W/(omega*rho_m*tau_A + I*rho_m)
    10	   NONZERO s11cc2ClosedSlabOperator/LAB_HELD/RHO4_CONSTANT/VALUE/E_W -> -2*_Dummy_49*I*Lambda_X_0*W_0*epsilon_shape*eta_bg*sigma_W/(2*omega*rho_m*tau_X + 2*I*rho_m) + _Dummy_49*W_0*epsilon_shape*eta_bg*sigma_W + 2*_Dummy_59*I*Lambda_X_0*W_0*epsilon_shape*eta_bg*sigma_W/(2
    11	pole-set mismatches: 0
    12	container shape mismatches: 0
    13	compact leaves containing Dummy: 0 []
    14	compact leaves containing UnevaluatedExpr: 0 []
    15	VERDICT: FAIL
.codex_seen 0 bytes
.grok_seen 0 bytes
.seen_arcverify 0 bytes
.seen_repairdir_codex 0 bytes
.seen_repairdir_grok 0 bytes
ABLATE_copy.py 31876 bytes
CLAUDE_old.md 8172 bytes
RESUME_PROMPT.md 5403 bytes
S11_census_instr_repair4_review_prompt.md 8623 bytes
S11_closure_codex_check.md 4574 bytes
S11b_docs_codex_check.md 3201 bytes
S11b_exports_committed_backup.py 14766845 bytes
S11b_exports_r1.py 14598922 bytes
S11b_exports_r1repair.py 14766348 bytes
S11b_snags_codex_consult.md 6869 bytes
S11c_a_sympy_engine_leg_findings.md 20704 bytes
S11c_b_spec_design.md 8524 bytes
S11c_c1_sympy_review_prompt.md 12531 bytes
abl_baseline.err 0 bytes
abl_baseline.out 11008938 bytes
abl_collapse.err 0 bytes
abl_collapse.out 27685183 bytes
abl_current_argstrip.py 79209 bytes
abl_dropA.err 0 bytes
abl_dropA.out 29416186 bytes
abl_dropD.err 0 bytes
abl_dropD.out 27321628 bytes
abl_integral_nobounds.py 79273 bytes
abl_key_facecollapse.py 79289 bytes
ablate_A_zerohessian.stderr 0 bytes
ablate_A_zerohessian.stdout 518875 bytes
ablate_B_swapsrc.stderr 0 bytes
ablate_B_swapsrc.stdout 877980 bytes
ablate_C_droprecon.stderr 407 bytes
ablate_C_droprecon.stdout 947163 bytes
ablate_D_frozenperturb.stderr 409 bytes
ablate_D_frozenperturb.stdout 879103 bytes
ablate_density_trace.err 508 bytes
ablate_density_trace.out 2244 bytes
ablate_density_trace.wl 3045 bytes
ablate_density_trace_clean.out 2244 bytes
ablate_gate.py 1860 bytes
ablate_parse.py 2779 bytes
ablation_a.out 1335 bytes
ablation_b.out 2248 bytes
ablation_b_demo.py 1127 bytes
ablation_c.out 879 bytes
ablation_c_empirical.out 885 bytes
ablation_d.out 488 bytes
ablation_d_demo.py 1245 bytes
absorb_checks.py 2481 bytes
absorb_checks.stdout 252 bytes
added.lines 1088933 bytes
anchored_dump.txt 1505412 bytes
apply_ref_tidy.py 4282 bytes
artifact_baseline.stderr 0 bytes
artifact_baseline.stdout 877980 bytes
artifact_dump.stderr 0 bytes
artifact_mu_theta_dump.srepr 70578 bytes
b.lines 9560947 bytes
b.names 1262 bytes
b.sorted 24562441 bytes
badarg.err 77 bytes
badarg2.err 70 bytes
base_ADMISSIBILITY_OPERATOR_OPERAND.out 4935 bytes
base_ADMISSIBILITY_RESIDUAL.out 4592 bytes
base_ADMISSIBILITY_SUPPORT_OPERAND.out 68 bytes
base_COUPLING_KERNEL.out 0 bytes
base_COUPLING_KERNEL_TERM_ORIGINS.out 9514 bytes
base_DIMENSIONS.out 2587 bytes
base_ENERGY_BASIS_COUNT.out 57 bytes
base_ENERGY_BASIS_NEW_INVARIANTS.out 431 bytes
base_ENERGY_BASIS_OMISSIONS.out 378 bytes
base_ENERGY_BASIS_VARIABLE.out 413 bytes
base_HOMOGENEITY_BASE_OPERAND.out 696 bytes
base_HOMOGENEITY_CONTROL_OPERAND.out 726 bytes
base_HOMOGENEITY_RESIDUAL.out 673 bytes
base_MU_THETA_OPERATOR.out 61349 bytes
base_SLAB_OPERATOR.out 771182 bytes
base_SLAB_OPERATOR_TERM_ORIGINS.out 11932 bytes
baseline.out 0 bytes
baseline.tags 14754674 bytes
baseline_partial_D2.out 416972 bytes
baseline_partial_D4.out 1145188 bytes
baseline_stderr.txt 0 bytes
baseline_stdout.txt 90723081 bytes
baseline_tests.out 1799 bytes
basis_rank_probe.py 6081 bytes
big_scale_test.out 142 bytes
bite_liveness.py 1748 bytes
c.lines 9696187 bytes
c.names 1262 bytes
c.sorted 24884411 bytes
c1_baseline.py 88470 bytes
c1_baseline_commit_msg.txt 2016 bytes
c1_checkpoint_verify_prompt.md 5576 bytes
c1_done.marker 16 bytes
c1_memlog.txt 2152 bytes
c1_migration_commit_msg.txt 1435 bytes
c1_nofold.py 88990 bytes
c1_reconcile_resume_prompt.md 4971 bytes
c1_repair_commit_msg.txt 2125 bytes
c1_repair_resume_prompt.md 4706 bytes
c1_repaired.py 88938 bytes
c1_run.stderr 0 bytes
c1_run.stdout 37872136 bytes
c1_run2.stderr 0 bytes
c1_run2.stdout 37872136 bytes
c1_run3.stderr 0 bytes
c1_run3.stdout 76369533 bytes
c1_wl_arc_verify_prompt.md 4726 bytes
c2_build_directive_skeleton.md 7701 bytes
c2_export_gate_commit_msg.txt 1758 bytes
c2_export_repair_commit_msg.txt 2378 bytes
c2_resume_prompt.md 4497 bytes
c2_review_commit_msg.txt 2313 bytes
chain_report.txt 568 bytes
claude_md_streamline_proposal.md 39247 bytes
codex_c1_arc_verify.log 3036682 bytes
codex_c1_arc_verify.marker 16 bytes
codex_c1_build.log 8695564 bytes
codex_c1_build.marker 7 bytes
codex_c1_comparator_build.log 24431074 bytes
codex_c1_comparator_directive_review.txt 536809 bytes
codex_c1_migration_review.log 2912385 bytes
codex_c1_repair.log 5094888 bytes
codex_c1_repair.marker 16 bytes
codex_c1_repair_build.log 118095 bytes
codex_c1_repair_build.marker 7 bytes
codex_c1_verify.log 476598 bytes
codex_c1_verify.marker 16 bytes
codex_closeout_verify.log 403621 bytes
codex_compact_verify.md 4440 bytes
codex_compact_verify_out.txt 1058579 bytes
codex_design_review.txt 864275 bytes
codex_export_minimum.txt 355599 bytes
codex_fix2_final.md 22858 bytes
codex_jetdepth_consult.log 400064 bytes
codex_ledger_fold_build.log 1324316 bytes
codex_ledger_fold_repair.log 3895974 bytes
codex_ledger_fold_repair2.log 656261 bytes
codex_leg.marker 7 bytes
codex_pid.txt 8 bytes
codex_pin_verify.log 467176 bytes
codex_repair2_build.log 894298 bytes
codex_repair2_dir.log 255845 bytes
codex_repair_dir_review.marker 7 bytes
codex_repair_dir_review.txt 2385652 bytes
codex_s11c_c1_directive_review.txt 1902614 bytes
codex_s11c_c1_sympy_build.log 34886406 bytes
codex_s11c_c1_sympy_build_v2.log 556330 bytes
codex_s11c_c1_sympy_build_v3.log 419673 bytes
codex_specpin.log 2075388 bytes
codex_verify_remediation.log 2325556 bytes
codex_verify_remediation.marker 16 bytes
codex_verify_t7_prep.log 1273892 bytes
codex_wl_build.log 29610675 bytes
codex_wl_build.marker 16 bytes
codex_wl_dirreview.log 1172147 bytes
codex_wl_dirreview.marker 16 bytes
codex_wl_repair.log 6348097 bytes
codex_wl_repair.marker 16 bytes
codex_wl_repairdir.log 290031 bytes
codex_wl_repairdir.marker 16 bytes
commit_msg_89b.txt 3141 bytes
commit_msg_buildleg.txt 635 bytes
commit_msg_checkpoint.txt 1328 bytes
commit_msg_clearance.txt 1983 bytes
commit_msg_closeout.txt 2040 bytes
commit_msg_out.txt 395 bytes
commit_msg_proposal.txt 506 bytes
compact_msg.txt 831 bytes
comparator_decision_DRAFT.md 5549 bytes
complete_prefix.txt 3369118 bytes
completeness.out 2701 bytes
completeness_probe.py 3257 bytes
coupling_keyacct.py 1001 bytes
coupling_report.txt 8173 bytes
crux_frozen_merge_is_spurious.py 2014 bytes
current.out 24884411 bytes
current.tags 14869524 bytes
d4_prefix.diff 47274 bytes
derive_density_trace.out 1926 bytes
derive_density_trace.py 4263 bytes
derive_material_advected_Te_density_trace.py 8222 bytes
derive_material_advected_Te_density_trace.stdout 5435 bytes
diag_evolution.py 1374 bytes
diag_growth.py 2532 bytes
diff_transcript.out 2231 bytes
diff_transcript.py 1365 bytes
dim_block.wl 2671 bytes
drop_forceHoldEw.out 7568 bytes
drop_sigmaW_MU_THETA.out 4680 bytes
dtn_operator_line.txt 42182 bytes
em_doc_cleanup_codex_review.md 3651 bytes
em_doc_cleanup_plan.md 7060 bytes
energy_baseline.err 0 bytes
energy_baseline.out 15236017 bytes
energy_face_line.txt 80047 bytes
energy_repaired.err 0 bytes
energy_repaired.out 18046347 bytes
energy_residual_line.txt 412638 bytes
engine.diff 45603 bytes
engine_copy.wl 105823 bytes
engine_stderr.txt 0 bytes
engine_stdout.txt 17981245 bytes
er_BASELINE.srepr 7867 bytes
er_FARFIELD_BRANCH_SIGN_REVERSED.srepr 7854 bytes
er_TRACTION_SIGN_REVERSED.srepr 7828 bytes
extract_grades.py 3040 bytes
face_repaired.err 0 bytes
face_repaired.out 38076848 bytes
faceflux_probe.py 1520 bytes
faceflux_probe2.py 2022 bytes
faceshift_current_check.py 2189 bytes
faceshift_density_adjudication.py 5847 bytes
faceshift_postfix_check.py 2350 bytes
faceshift_residual.py 2541 bytes
faceshift_verify.py 2479 bytes
final_checks.out 728 bytes
final_checks.py 4543 bytes
final_d2.diff 39805 bytes
final_d4.diff 23692 bytes
fix1_leg_pids.txt 16 bytes
flat_repaired.out 23203 bytes
fold4_probe.out 480 bytes
fold4_probe.py 1270 bytes
full_repaired.err 0 bytes
full_repaired.out 107173098 bytes
full_repaired.timelog 0 bytes
full_residual.py 1271 bytes
full_traction.py 1276 bytes
gate_FULL_VOCAB_py.out 295 bytes
gate_FULL_VOCAB_wl.out 25 bytes
gate_MISSING_c_s0_py.out 279 bytes
gate_MISSING_c_s0_wl.out 25 bytes
gate_ONLY_c_s0_py.out 86 bytes
gate_ONLY_c_s0_wl.out 25 bytes
gen_s11b_build_directive_twins.sh 6049 bytes
gen_s11b_decisions_measurements.sh 5347 bytes
gen_s11b_fix_round1_twin.sh 4540 bytes
gen_s11b_fix_round2_twin.sh 2667 bytes
gen_s11b_shared_physics_measurements.sh 6744 bytes
gen_s11c_a_recompose_meas.sh 1737 bytes
gen_s11c_a_repair2_meas.sh 1245 bytes
gen_s11c_a_repair_meas.sh 1537 bytes
gen_s11c_a_spec_meas_v2.sh 2315 bytes
gen_s11c_a_spec_measurements.sh 3469 bytes
gen_s11c_a_sympy_build_meas.sh 2047 bytes
gen_s11c_decisions_measurements.sh 4527 bytes
gen_twin_fix4.sh 4862 bytes
gen_twin_fix5.sh 3280 bytes
gen_twin_fix5b.sh 3659 bytes
gen_twin_fix6.sh 4201 bytes
grok_c1_comparator_directive_review.txt 20950 bytes
grok_c1_comparator_review.marker 7 bytes
grok_c1_comparator_review.txt 12639 bytes
grok_c1_migration_review.log 12645 bytes
grok_c1_repair_review.log 12585 bytes
grok_c1_repair_review.marker 16 bytes
grok_c1_review.log 12396 bytes
grok_c1_review.marker 16 bytes
grok_design_review.txt 14014 bytes
grok_export_minimum.txt 20117 bytes
grok_jetdepth_consult.log 8888 bytes
grok_ledger_fold_review.txt 10346 bytes
grok_leg.marker 7 bytes
grok_repair2_dir.log 14382 bytes
grok_repair2_review.log 10529 bytes
grok_repair_dir_review.marker 7 bytes
grok_repair_dir_review.txt 8656 bytes
grok_repair_review.marker 7 bytes
grok_repair_review.txt 12306 bytes
grok_s11c_c1_directive_review.txt 13563 bytes
grok_specpin.log 12522 bytes
grok_wl_build_review.log 10215 bytes
grok_wl_build_review.marker 16 bytes
grok_wl_dirreview.log 19840 bytes
grok_wl_dirreview.marker 16 bytes
grok_wl_repair_review.log 7301 bytes
grok_wl_repair_review.marker 16 bytes
grok_wl_repairdir.log 17857 bytes
grok_wl_repairdir.marker 16 bytes
h1234_err.txt 0 bytes
h1234_out.txt 0 bytes
h234_err.txt 0 bytes
h234_out.txt 3520 bytes
h5_err.txt 0 bytes
h5_out.txt 1199 bytes
harness.py 1001 bytes
harness_defs.wl 78994 bytes
harness_r1234.py 6823 bytes
harness_r234.py 4437 bytes
harness_r5.py 2958 bytes
head_full.wl 105718 bytes
head_tags_sofar.txt 2232 bytes
head_trunc.err 470 bytes
head_trunc.out 6963306 bytes
head_trunc.wl 79886 bytes
indep_mu_theta.py 5225 bytes
indep_mu_theta.stdout 331 bytes
independent_verify.out 0 bytes
independent_verify.py 8603 bytes
independent_verify2.out 9126 bytes
independent_verify2.py 6065 bytes
kernel_baseline.err 0 bytes
kernel_baseline.out 120592 bytes
kernel_nofold.err 0 bytes
kernel_nofold.out 119960 bytes
kernel_repaired.err 0 bytes
kernel_repaired.out 120609 bytes
killed_prefix.txt 3368978 bytes
launch_wl_build.sh 834 bytes
launch_wl_dirreview.sh 1291 bytes
layer_abl.py 62273 bytes
layer_patched.py 62519 bytes
ledger_fold_ORIG.py 14412 bytes
loader.py 403 bytes
material_transport_adjudication.py 7544 bytes
measure_reconcile.py 6362 bytes
measure_v2.py 5299 bytes
measure_v3.py 10508 bytes
mg_classify.py 3584 bytes
mg_parse.py 2346 bytes
mg_parse2.py 2019 bytes
mg_terms.py 2133 bytes
new_ROOT2.eq 8414 bytes
new_ROOT3.eq 8381 bytes
next_prompt_c2_build_v2.md 6044 bytes
next_prompt_c2_wolfram.md 5388 bytes
next_prompt_draft.md 4941 bytes
old_ROOT2.eq 8414 bytes
old_ROOT3.eq 8381 bytes
opus_ablated_counter_S11c_a.py 79244 bytes
opus_coupling_kernel_lines.txt 1691071 bytes
opus_coupling_token_trace.out 951 bytes
opus_coupling_token_trace.py 2016 bytes
opus_fixtures_postfix.out 336 bytes
opus_form_ablation.py 2136 bytes
opus_form_ablation_counter.out 1407 bytes
opus_form_ablation_prefix.out 1398 bytes
opus_injectivity.out 608 bytes
opus_injectivity.py 4453 bytes
opus_prefix_S11c_a_comparator.py 79230 bytes
opus_s11ca_postfix.err 0 bytes
opus_s11ca_postfix.out 193521648 bytes
opus_s11ca_prefix.err 0 bytes
opus_s11ca_prefix.out 49672492 bytes
opus_s11ca_prefix2.err 0 bytes
opus_s11ca_prefix2.out 193584648 bytes
opus_s11ca_py_symbols.txt 2831 bytes
opus_s11ca_rawdiff.txt 377684 bytes
opus_trace_canon.out 2212 bytes
opus_trace_canon.py 2228 bytes
opus_verify_s11ca_diff.out 5454 bytes
opus_verify_s11ca_diff.py 3393 bytes
opus_verify_s11ca_diff2.out 372 bytes
opus_verify_s11ca_diff2.py 4238 bytes
orch_refreeze_verify.err 440 bytes
orch_refreeze_verify.wl 2138 bytes
p1_wl_residual_emit_decision_DRAFT.md 3983 bytes
patched_run.err 0 bytes
patched_run.out 16951728 bytes
peek_conormal.py 781 bytes
peek_evolution.py 988 bytes
port_line.txt 62269035 bytes
primary_faceshift_probe.out 1746 bytes
primary_faceshift_probe2.out 1863 bytes
probe.py 753 bytes
probe1.wl 273 bytes
probe2.wl 429 bytes
probe3.wl 628 bytes
probe_87.py 4791 bytes
probe_ADMISSIBILITY_OPERATOR_OPERAND.out 7968 bytes
probe_ADMISSIBILITY_RESIDUAL.out 8500 bytes
probe_ADMISSIBILITY_SUPPORT_OPERAND.out 2744 bytes
probe_ENERGY_BASIS_COUNT.out 293 bytes
probe_MU_THETA_baseline.out 4684 bytes
probe_blast.py 6918 bytes
probe_correct_el.py 2737 bytes
probe_decomp.py 4062 bytes
probe_match.py 2427 bytes
probe_omitted.py 1538 bytes
prod_run.err 0 bytes
prod_run.out 33996954 bytes
py_commit_msg.txt 2191 bytes
py_emit_nonlocal.txt 515 bytes
py_fams.txt 1283 bytes
py_form_q.txt 348 bytes
py_freeze_probe.py 1820 bytes
py_freeze_probe2.py 1674 bytes
py_newinv.txt 99677 bytes
py_nonlocal.txt 886 bytes
py_primary.txt 251 bytes
r3_spur_line.txt 518 bytes
r3_spur_sheet.txt 1202 bytes
raw_wl_evo.py 1096 bytes
real_run.stderr 0 bytes
real_run_full.out 224727199 bytes
real_run_status.txt 19 bytes
recon_E_detail.out 1742 bytes
recon_E_detail.py 4127 bytes
recon_E_final.out 1741 bytes
recon_E_final.py 2521 bytes
recon_closure.out 9802 bytes
recon_closure.py 3423 bytes
recon_keys.out 4938 bytes
recon_keys.py 2558 bytes
recon_load.out 712 bytes
recon_load.py 1067 bytes
recon_min_and_E.out 4781 bytes
recon_min_and_E.py 3823 bytes
recon_minimal.out 0 bytes
recon_minimal.py 2296 bytes
regime_baseline.err 0 bytes
regime_baseline.out 851170 bytes
regime_repaired.err 0 bytes
regime_repaired.out 851643 bytes
removed.lines 766963 bytes
rep_resid.srepr 79869 bytes
repair2_build.marker 36 bytes
rerun_03_measure_operands.stdout 3561 bytes
rerun_05_d2_rank_mod_disc.stdout 3541 bytes
rerun_agent_04_bounded_routes.stdout 5349 bytes
rerun_agent_08_fixed_parser_rerun.stdout 3051 bytes
resume_prompt_c1_build.md 5644 bytes
resume_prompt_c1_migration.md 6101 bytes
resume_prompt_draft.md 3865 bytes
row_residual_report.py 3475 bytes
run_c1_repair.sh 491 bytes
run_c1_verify.sh 427 bytes
run_coverage.out 5877 bytes
run_current_fold.out 1238 bytes
run_current_fold.py 2325 bytes
run_integral_canon.out 911 bytes
run_integral_canon.py 3028 bytes
run_keying.out 843 bytes
run_keying.py 1777 bytes
run_probes.py 5849 bytes
run_repair2_build.sh 460 bytes
s11b_slice_map_check.out 266 bytes
s11b_slice_map_check.py 1327 bytes
s11b_transverse_mu_check.out 264 bytes
s11b_transverse_mu_check.py 1830 bytes
s11b_zperm_slicemap_sign.py 1416 bytes
s11ca_axis_census.py 6306 bytes
s11ca_axis_keys.txt 80311 bytes
s11ca_case_census.py 5536 bytes
s11ca_census_keys.txt 388039 bytes
s11ca_density_adjudication.py 5048 bytes
scale_capability.out 590 bytes
scale_capability.py 3697 bytes
scout_axis_census.py 6308 bytes
scout_axis_keys.txt 74621 bytes
scout_keytree.py 5517 bytes
scout_keytree_detail.txt 22746 bytes
scout_leafcount.py 3879 bytes
scout_skeleton.py 3915 bytes
selfctrl.py 1421 bytes
slab_row.txt 18147887 bytes
small_families.out 1084919 bytes
smoke.wl 522 bytes
spliced.out 8858289 bytes
step0_run.sh 2000 bytes
step0_wl_op_kernel.wl 2952 bytes
sweep_pid.txt 7 bytes
sympy_completeness_final.log 81876 bytes
sympy_out_stderr.log 803 bytes
t7_resume_prompt.md 6146 bytes
tags_complete.txt 118 bytes
tags_killed.txt 118 bytes
test_2_3.out 2901 bytes
test_2_3.wl 5110 bytes
test_5.out 2044 bytes
test_5.wl 1046 bytes
test_6_7.out 6288 bytes
test_6_7.wl 2436 bytes
test_6b.out 7930 bytes
test_6b.wl 1585 bytes
test_dim_fix.wl 1787 bytes
test_dim_iso.wl 2588 bytes
test_evolution.py 1241 bytes
test_grade.out 3201 bytes
test_grade.wl 2482 bytes
test_indep.out 526 bytes
test_indep.wl 572 bytes
test_refreeze.out 934 bytes
test_refreeze.wl 1904 bytes
test_refreeze2.out 3259 bytes
test_refreeze2.wl 1980 bytes
test_traction_fix.py 1259 bytes
tower_depth_probe.out 1277 bytes
tower_depth_probe.py 4315 bytes
v2_diff.txt 45474 bytes
verify_remediation_state.md 4078 bytes
verify_t7_prep.md 7136 bytes
watchdog_pid.txt 7 bytes
wfail_sample.txt 884 bytes
wip4_msg.txt 2151 bytes
wl_adm_op.txt 7327 bytes
wl_commit_msg.txt 1893 bytes
wl_face_split.py 2497 bytes
wl_fams.txt 1283 bytes
wl_form_q.txt 228 bytes
wl_newinv.txt 40318 bytes
wl_nonlocal.txt 886 bytes
wl_out_stderr.log 2009 bytes
wl_production_run.sh 2545 bytes
wl_provenance_probe.err 440 bytes
wl_provenance_probe.out 1079 bytes
wl_provenance_probe.wl 3003 bytes
wl_remediation_resume_prompt.md 4303 bytes
wl_resume_prompt.md 3034 bytes
wl_single_case_urow.wl 2045 bytes
wlrepair1_baseline.err 2009 bytes
wlrepair1_baseline.out 59389690 bytes
wlrepair1_baseline.wl 79824 bytes
wlrepair1_fulldiff.txt 36770 bytes
wlrepair1_probe_block.wl 4582 bytes
wlrepair1_probe_run.err 2009 bytes
wlrepair1_probe_run.out 76660903 bytes
wlrepair1_probe_run.wl 100199 bytes
wlrepair1_repaired.err 2009 bytes
wlrepair1_repaired.out 76656928 bytes
wlrepair1_repaired.wl 95617 bytes
wlreview1_FINDINGS_SUMMARY.txt 1380 bytes
wlreview1_baseline.err 440 bytes
wlreview1_baseline.out 71082 bytes
wlreview1_baseline.wl 79824 bytes
wlreview1_derivation.out 1456 bytes
wlreview1_derivation.py 6250 bytes
wlreview1_harness.err 440 bytes
wlreview1_harness.marker 7 bytes
wlreview1_harness.out 11295 bytes
wlreview1_harness.wl 48204 bytes
wt_ablate.err 440 bytes
wt_ablate.out 6767226 bytes
wt_ablate.wl 79227 bytes
wt_tags.txt 2232 bytes
wt_trunc.err 482 bytes
wt_trunc.out 6811091 bytes
wt_trunc.wl 78708 bytes
x1_independent_basis_count.py 3179 bytes

exec
/bin/bash -lc "nl -ba CLAUDE.md | sed -n '1,180p' && git show aa76105a:research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md | sed -n '1,35p' && git log --format='%h %s' -- research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md | head -n 5 && git ls-tree -r --name-only aa76105a | rg 'S11c_c2_export_repair_(builder_report|rereview)'" in /var/projects/toy_physics
 succeeded in 0ms:
     1	# How we work
     2	
     3	Seventeen controls. Every one exists because ignoring it cost a session — none is removed here. This
     4	**2026-09-05 restructure** groups them so the rule is scannable and the evidence that earned it is kept but
     5	out of the way; it changes presentation, **never authority** (design + preservation map:
     6	`CLAUDE_streamline_proposal_2026-09-05.md`; the M/E/G/S controls were preservation-checked against the
     7	prior version by two independent legs before commit).
     8	Full process: `docs/development_pipeline.md`. What we're building: `docs/development_plan.md`.
     9	Where we are: `STATUS.md`.
    10	
    11	⭐ **Correctness override (R11):** cost is **never** a reason to drop a control, narrow a check, or skip a
    12	leg. Only the user scales work down. Nothing below — no grouping, no budget, no "administrative" label —
    13	weakens a gate.
    14	
    15	---
    16	
    17	## At a glance
    18	
    19	1. **Identify the object and the artifact's role.** Supply verified equations; withhold the target-answer
    20	   criterion (the expected-value acceptance test, only); construct independently; preserve disagreement. (M1, M2)
    21	2. **Measure every claim.** Compute and emit *before* guards; interpret in records; attach the command and
    22	   its literal stdout in a `_measurements/` file. This binds the orchestrator too. (E1)
    23	3. **Pick two non-author legs.** Orchestrator-written → Codex + Grok; Codex-written → fresh Claude + Grok.
    24	   Launch on sight, before inspecting the result. (G1, G3)
    25	4. **Use the right exit.** *Decision list:* two legs, verify, fold **once**, go. *Physics spec / script /
    26	   record:* **review-until-clear** — iterate leg→fold→leg until nothing outstanding changes what is
    27	   computed or may be claimed. The filename confers no exception. (G2, G4)
    28	5. **Verify and preserve.** Adjudicate findings by evidence (a finding is not a mandate); ablate physics
    29	   checks by FORM; keep varying quantities live; change a defect-breeding author; no commit before both
    30	   reports; commit the reviewed baseline before a repair overwrites it. (E2, G4, M3)
    31	6. **Enforce, don't merely promise.** Bound handoffs and leave visible evidence; never revive quarantine or
    32	   drop a control for cost; follow the annex/GIN policy for v3 `.out`. (S1, R11)
    33	
    34	**Scope precedence:** classify an artifact by **function, then authorship** — the filename confers nothing.
    35	A packet serving several roles inherits **all** applicable controls. The decision-list one-pass exception
    36	(G2) **never** caps review of a physics spec, even one embedded in a decision list. "Review-until-clear" =
    37	G4's substantive exit (nothing outstanding changes what is computed or may be claimed), never consensus or
    38	green labels.
    39	
    40	### Artifact → review discipline
    41	
    42	Reviewer key: **O** = orchestrator-written → Codex + Grok; **Cx** = Codex-written → fresh Claude agent +
    43	Grok. Mixed/changed authorship needs an explicit valid non-author pairing — never a silent reuse or a
    44	reduced count.
    45	
    46	| Artifact / role | Legs (by author) | Method | Pass / fold discipline | Advance / commit gate |
    47	|---|---|---|---|---|
    48	| **Pre-builder decision list** | Two; O (always orchestrator-written) | Check the requested decisions + their supporting evidence; ⛔ no fictional-script ablation | **One two-leg pass**, verify findings, fold **once** — ⛔ no iterate-to-green | No builder before both reports + fold. An unresolved issue that changes what is computed or may be claimed routes to the applicable spec/build gate (R10 scope); physics-bearing content also meets the spec row. No commit before both reports. |
    49	| **Physics spec / shared spec / physics-bearing directive** | Two; O or Cx by actual author | Review the requested physics, complete premises, recipe/answer leakage; substantiate physics claims independently. Defer executable script-control tests to the build. | **Review-until-clear**: repair/re-review any change to computation or claims | Clear before using as governing physics. Directive review never pays the build tax. Both reports before any commit; reviewed baseline preserved before overwrite. |
    50	| **Script / physics-bearing build** | Two; O or Cx; normally Cx for `/build` | Independent derivation scripts + literal stdout; every load-bearing check ablated; **mandatory FORM ablation**; one-sided corruption for independence; emit-before-guard + output checks | **Review-until-clear**; change author when repairs breed defects (⛔ never a fourth fold) | Launch legs before inspecting results; both usable reports; preserve reviewed baseline before repair; accept only on substantive clearance. Serialize dual Mathematica ablations. |
    51	| **Step record / `.tex` card / physics-bearing prose** | Two; O or Cx; ⛔ never chosen by extension | Source-first fidelity; quote both sides; ⛔ no build directive in the packet; a measured physics claim still needs its command + stdout; for cards check suppressed macro fields | **Review-until-clear** about what may be claimed | Own artifact/version review required; script review alone is not record review. Both reports before commit. |
    52	| **Other / claimed non-physics** | Two if physics-bearing (G1); if you claim it is not, record why | First record whether it changes computation, premises, checks, or claims — if so it is physics-bearing, route to a row above | ⛔ Do not infer a one-pass or zero-review exemption from the suffix or an "administrative" label | The two-report commit gate has no non-physics exception. |
    53	
    54	**Observable gate record** *(new 2026-09-05, from the approved proposal; implements R12's observable-artifact
    55	principle — illustrative, not a newly mandated gate)*: artifact path/version + role(s);
    56	author(s); handed-input manifest; the identical rendered review prompt; two separately-attributed reports
    57	with script/stdout paths where required; verified finding dispositions; iteration/author-change record; the
    58	reviewed-baseline commit and the accepted version. It documents controls; it does not substitute for
    59	carrying them out, and it is not an answer quarantine.
    60	
    61	---
    62	
    63	## M · Specify and construct
    64	
    65	**M1 — Two engines exist so they can disagree, and the disagreement is the measurement.** *(was R1, R6)*
    66	Independent construction, not hidden answers. A disagreement is a **finding** — ⛔ never try to make
    67	divergence impossible with more careful prose; that defeats the reason there are two engines. ⛔ Never treat
    68	a disagreement as a builder target to eliminate; it is a finding on the orchestrator's side.
    69	
    70	**M2 — Name the object; supply verified equations; withhold only the target answer.** *(was R3, R5)*
    71	Name the **object**, ⛔ not the recipe: if a review argues *how* to derive something — is this quotient
    72	well-defined, is this weight unique — the question was manufactured by specifying a derivation path. Ask for
    73	the object; let the engine hand over what it built. The spec says what to **compute** — ⛔ never what
    74	anything equals, is expected, or was measured. Withhold exactly one thing: an acceptance criterion
    75	referencing an expected value (a builder iterating to exit 0 converges on any target it can see). Supply
    76	everything else, as **equations**, and label supplied objects as supplied and unfalsifiable in that build.
    77	Naming an object does not waive first-principles derivation, form ablation, or evidence. Prior art's result
    78	for our object is never a supplied premise (M3).
    79	
    80	**M3 — Prior art is an oracle, not a premise; keep every varying quantity live.** *(was R16, R17)*
    81	Prior art is an **oracle**, never a premise: check our computed result against it; ⛔ never assume its result
    82	for our object — its conditions may not be ours. A freeze is a red flag; a **required freeze is the finding**
    83	— treating a quantity that VARIES as if it were constant manufactures a wrong-but-consistent answer both
    84	engines share, which the comparator misreads as agreement (the measured account: L-R17). Keep every varying
    85	quantity **LIVE** and differentiate it; when a method seems to require holding one fixed to proceed, that
    86	requirement is the measurement, not a step. Caught only by a ground-truth anchor (verify *its* conditions —
    87	never replace the derivation with its answer) or a variable-coefficient/form ablation (E2). *(evidence: L-R17)*
    88	
    89	---
    90	
    91	## E · Produce evidence
    92	
    93	**E1 — A script prints computed objects; a record interprets them; every claim carries its command.**
    94	*(was R2, R4)*
    95	A script **prints** computed objects — ⛔ it never states conclusions. Emit both operands and the residual,
    96	**then** guard (a residual asserted zero always prints `0` and carries no information). Interpretation
    97	belongs to the step record. ⛔⛔ **This binds the orchestrator too** — the half I keep exempting: every
    98	review-leg prompt I write says *"a prose derivation is worth nothing; show the script and its literal
    99	stdout, or the claim is discarded,"* and the same standard applies to anything I write. **A claim about an
   100	artifact carries the command that produced it**, and the commands + their literal output go in a
   101	`_measurements/` file beside the document. If I am deciding in prose what the engines should compute, I have
   102	**inverted the method** — fix the instrument; ⛔ don't reason around it. *(evidence: L-R2, L-R4)*
   103	
   104	**E2 — Ablate to test; FORM tests physics, COEFFICIENT tests arithmetic; a check against a nonexistent
   105	script goes to the build review.** *(was R14, R17)*
   106	**Ablate** to test; ⛔ don't read. A **FORM** control tests physics (it leaves the family); a **COEFFICIENT**
   107	control tests arithmetic (it does not). Demand a script and its literal stdout — a prose re-derivation is the
   108	same defect relocated into the review. The variable-coefficient/form ablation is what catches the M3 freeze.
   109	Controls written against a script that **does not exist** cannot be ablated, so they get reviewed by reading
   110	(the weakest instrument): put them in the **build** review, ⛔ not in a document about a build. Computed means
   111	*dependent on the derivation*, not merely CAS-shaped: a hand-typed algebraic payload, an answer-bearing tag,
   112	or a suppressed-identical payload is not computed — the build-skill's action/ansatz-only construction,
   113	object-only tag names, independent routes, and value-independent emission still bind. *(evidence: L-R14)*
   114	
   115	---
   116	
   117	## G · Review and advance
   118	
   119	**G1 — Whatever writes does not review; pick the two legs by authorship.** *(was R7)*
   120	**Whatever writes does not review.** Two independent legs on anything **physics-bearing** — and a spec both
   121	engines read is physics-bearing, because an error there makes both engines agree on the same wrong thing.
   122	**Orchestrator-written → Codex + Grok. Codex-written → fresh Claude agent + Grok.** Choose legs by
   123	authorship, ⛔ never by file type; record mixed/changed authorship and establish a valid non-author pairing
   124	before proceeding (⛔ never silently reuse a contributor or reduce the count). My own verification of findings
   125	is required (G4) but is **never** a leg. *(evidence: L-R7)*
   126	
   127	**G2 — The decision-list gate — and ONLY the decision list gets the one-pass exception.** *(was R7 trigger;
   128	disambiguated 2026-09-05)*
   129	**TRIGGER — no builder launches until its decision list has had two legs.** The decision list is
   130	orchestrator-written and is the one artifact the *builder* trusts: everything downstream is checked twice,
   131	the list is checked zero times. **One pass, then fold and go — never iterated to green.** ⛔⛔ **The one-pass
   132	exception is the decision list's alone.** A physics **spec** (`SHARED_PHYSICS`, or anything both engines
   133	read), a **script**, and a **physics-bearing record** are **reviewed until clear** under G4 — the filename
   134	confers no exception, and a decision-list pass never clears a spec embedded in it. If a packet serves both
   135	roles, it satisfies **both** gates. *(evidence: L-R7)*
   136	
   137	**G3 — Launch legs on sight, before inspecting the result; no commit before both report.** *(was R8, R9)*
   138	Launch legs **on sight**, before I look at the result — a self-check discharges the felt need for an
   139	independent one, and it is most convincing when it finds things. Before launch, verify only that the
   140	deliverable exists, is non-empty, and ran plausibly — ⛔ do not inspect its result content, and never feed
   141	one leg the other's output. **No commit before both legs report.** The commit is the last step. Reviewing
   142	the directive does not pay the tax for the build. *(evidence: L-R8)*
   143	
   144	**G4 — Verify each finding, iterate to clearance, change the author if folds keep breeding defects.**
   145	*(was R10, R13, R15)*
   146	**Stop when nothing outstanding changes what is computed or what may be claimed** — ⛔ not when both legs are
   147	green. "A leg that finds nothing is weak evidence" is my prior; ⛔ put it in a leg's prompt and it becomes a
   148	quota (keep it in rationale, out of rendered instructions). A finding is not a mandate — **verify it
   149	myself** (legs have been wrong in both directions); obtain **both** reports before adjudicating or editing.
   150	If successive revisions keep breeding defects in the material just changed, **change the author** — ⛔ don't
   151	fold a fourth time (author change waives neither review nor clearance). Commit the exact **reviewed
   152	baseline** before a repair overwrites it (recording unresolved findings); that preservation commit is not
   153	acceptance — accept the repaired result only after its own review and G4 clearance. *(evidence: L-R10, L-R13)*
   154	
   155	---
   156	
   157	## S · Enforce through observable artifacts
   158	
   159	**S1 — A prohibition is not a control; blindness is enforced by absence.** *(was R12)*
   160	**A prohibition is not a control.** Blindness is enforced by **absence** — by bounding what the builder is
   161	handed — ⛔ never by a sentence forbidding a read. A do-not-read list is a denylist, and a denylist means the
   162	architecture is wrong. The measured failure is absence of computation, not anchoring; **quarantine is cut
   163	and E1 replaced it** — ⛔ do not revive any cut quarantine mechanism under a new name (raw logs live outside
   164	the tree for hygiene, not as a blindness claim). ⛔⛔ **This applies to these controls too:** every one is
   165	prose I drift from under load, so the ones that hold are the ones that leave an artifact whose absence you
   166	can **see** (the at-a-glance gate record is one illustration, not a new mandate). *(evidence: L-R12)*
   167	
   168	Operational runbooks stay authoritative in the skills — `.claude/skills/build/SKILL.md` and
   169	`.claude/skills/review-legs/SKILL.md`. A root summary cannot replace their identical-prompt,
   170	first-principles-script/stdout, mandatory-form-ablation, source-staging + directive-exclusion,
   171	authorship-pairing, fixed-baseline, kernel-serialization/budget/sandbox, leak-probe, orphan-memory
   172	diagnosis, and real-deliverable-verification obligations; ⛔ do not infer a weaker reading from a summary.
   173	
   174	---
   175	
   176	## Repository infrastructure — the `.out` transcripts live in git-annex + GIN (set up 2026-09-01)
   177	
   178	The v3 CAS audit transcripts (`research/pde_ledger_v3/scripts/out/*.out` and
   179	`research/pde_ledger_v3/mathematica/out/*.out`, ~370 MB) are **git-annex content backed by GIN**
   180	(`gin.g-node.org/trevnorris/toy_physics`, public), NOT plain git blobs — one exceeded GitHub's 100 MB/file cap.
# S11c-c2 export repair — builder report

Build → verify → report completed. Publication-only; physics construction and emissions are unchanged.

The generated `scripts/S11c_c2_exports.py` is **22,441,522 bytes**, down from **60,516,900 bytes**: **38,075,378 bytes removed (62.92%; 2.70× smaller)**.

## Exact edit sites

Only `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` was edited as source:

- `EXPORT_ROOTS`, current line 48: removed `s11cc2SelfEnergyIncrement`; retained both declared deliverables.
- New adjacent publication helper `publication_compact`, current line 806.
- `publish`, current line 852: compact only case `VALUE` trees; bounded root-name displays; retain declaration schema, F9, closure, minimality, digest structure, and structural roundtrip; add strict emitted-versus-restored semantics and reciprocal guards, literal evidence and size measurements; install only after passing.
- The sole `export_key` dictionary inside `run`, current line 1051: removed increment export routing. The increment's existing emit remains.

AST comparison against the pre-edit source passed after excluding only these authorized sites. All other construction, extraction, control, grading, emission, and physics-loop code is identical. No construction change was needed.

## Representation and size

Transform: exact `collect` followed by `factor_terms(fraction=False)` on per-case VALUE expressions, with reciprocal/calculus atoms protected locally and restored immediately; Integral boundaries and container types retained; any subexpression whose reciprocal-power set would change keeps its original form. No expansion is used to produce export values, no CSE is used, and no temporary/hold symbol survives in the delivered values. Expansion occurs only in the separate semantic guard.

`s11cc2ClosedCouplingKernel` VALUE totals: 18,599,401 → 13,917,507 bytes (25.17% smaller). `s11cc2ClosedSlabOperator` VALUE totals: 15,031,549 → 8,422,000 bytes (43.97% smaller). `s11cc2SelfEnergyIncrement` VALUE totals: 11,997,906 → 0 bytes (100.00% smaller; absent from delta, still emitted).

Sizes below are UTF-8 byte lengths of each case's `srepr(VALUE)`, excluding its unchanged metadata. Baseline sizes were measured from the original 60,516,900-byte module; the final build's emitted sizes match those baseline sizes exactly. Zero means absent from the delta.

| Object | Case | Before bytes | After bytes | Reduction |
| --- | --- | ---: | ---: | ---: |
| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHO4_CONSTANT` | 4,355,993 | 3,375,460 | 22.51% |
| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHOBR_CONSTANT` | 3,407,032 | 2,596,031 | 23.80% |
| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 5,990,000 | 4,422,717 | 26.16% |
| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 4,846,376 | 3,523,299 | 27.30% |
| `s11cc2ClosedSlabOperator` | `LAB_HELD/RHO4_CONSTANT` | 4,002,816 | 2,271,711 | 43.25% |
| `s11cc2ClosedSlabOperator` | `LAB_HELD/RHOBR_CONSTANT` | 3,121,562 | 1,718,913 | 44.93% |
| `s11cc2ClosedSlabOperator` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 4,427,664 | 2,502,613 | 43.48% |
| `s11cc2ClosedSlabOperator` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 3,479,507 | 1,928,763 | 44.57% |
aa76105a S11c-c2 export repair COMMITTED — 60 MB -> 21.4 MB (publication-only, re-reviewed clear)
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_grok.txt
research/pde_ledger_v3/directives/_legs/S11c_c2_export_repair_rereview_prompt.md

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/_measurements/S11c_c2_physics_claude_leg.md && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_physics_grok_leg.txt | sed -n '1,240p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md | sed -n '1,160p'" in /var/projects/toy_physics
 succeeded in 0ms:
     1	# S11c-c2 self-energy fold — independent physics review (fresh Claude leg)
     2	
     3	Leg author: fresh general-purpose Claude agent (Codex/astra wrote the script → legs are fresh-Claude + Grok).
     4	Identical prompt: `directives/_legs/S11c_c2_physics_review_prompt.md`. Sandbox: `/tmp/c2_ablate_agent`
     5	(working tree never modified). Evidence: `/tmp/c2_ablate_agent/EVIDENCE_LOG.txt` + helper scripts/outputs there.
     6	
     7	**Method.** Derived the physics from `S11c_c2_SHARED_PHYSICS.md` + the real imports, then tested the script by
     8	(1) byte-offset extraction of the emitted `.out` objects, (2) loading the real folded model to inspect
     9	bindings, (3) two mandatory FORM ablations on the TRIAGE core, (4) one-sided corruption checks.
    10	
    11	## Findings by section
    12	
    13	**A · Fold map — SOUND.** `build_face`@477 returns only `{delta_p_±, d_w_delta_p_±}` slot substitutions — no
    14	`J_s` slot (no double-count). Jets computed by differentiation (`build_face`@506-508). Response elimination is
    15	a genuine 2×2 matrix inverse (`kernel_bridge`@396-398), not scalar division. Load-bearing bulk object is the
    16	two-momentum kernel `z_matrix=[[z0out,z1],[0,z0in]]`. Symbol map verified by model load: `face_velocity[LAB,+]
    17	= W_0·e_W_t·ε/2` dim `[1,-1,0]` = `s11cc1_V` dim `[1,-1,0]` — `V_s` is the interfacial normal velocity;
    18	`μ_θ`→`mu_theta_operator`.
    19	
    20	**B · ε-normalization — SOUND.** Increment MULTIGRADE = `{(1,0,0),(1,0,1),(1,1,0),(1,1,1)}`. Per-term scan (LAB
    21	54 / MATERIAL 318 terms): 0 terms without `epsilon_shape` (no O(1)), 0 terms with ε^≥2 (no spurious O(ε²)); the
    22	O(εη) coupling present at (1,1,0). The stray `(0,0,0)` is the resolvent-denominator grade-support
    23	over-approximation the code self-disclaims (`grades`@712-713). Strips exactly one ε.
    24	
    25	**C · Close-then-extract — SOUND.** `ORDERING_COMMUTATOR` nonzero (0.40–1.79 MB, all 4). FORM-ablation form1
    26	(freeze input leg `z0in=z0out`, the rule-17 one-leg freeze): increment sha `f4fba6b1→252d2552`. FORM-ablation
    27	form2 (`z1=0`, kill off-diagonal): increment 2,575,195→1,899,256 bytes. Closure threading load-bearing.
    28	
    29	**D1 · Density field-vs-field — SOUND.** `density[RHO4_CONSTANT] = rho_br·(1+eta_bg·w1_profile)` (live),
    30	`density[RHOBR_CONSTANT] = rho_br` (frozen); `build_face`@483 rebinds the bare c1-frozen slot to the live
    31	field. `DENSITY_LIVE_MINUS_FROZEN` carries x-dependent `w1_profile` + `eta_bg` — live-field difference, no ρ
    32	jet, ρ_br-based (not ρ_m).
    33	
    34	**D2 · Traction covector — SOUND.** `TRACTION_MECHANICAL_CONTRIB` = `Matrix(4,1)`, all 4 components nonzero;
    35	`traction.dot(virtual_velocity)`@792. Native covector, not collapsed.
    36	
    37	**D3 · DtN whole-form vs kernel — SOUND (non-vacuous).** `DTN_WHOLEFORM_DEPENDENCE` empty (all 4). Unfiltered
    38	increment∩dtn_operator = `{W_0,c_s0,eta_bg,omega,rho_m,sigma_W}` (background constants only); zero
    39	noncommutative whole-form symbols leak.
    40	
    41	**D4 · Traction–slab power pairing — SOUND.** Emits SLAB_POWER / KINETIC_STORED_POWER / TRACTION_POWER /
    42	FACE_GENERALIZED_POWER separately; residual = face_power − traction_power (two independent routes). One-sided
    43	`t_s`-flip `TRACTION_SIGN_RESIDUAL` = 2·traction_power, nonzero (4.5–7.9 MB).
    44	
    45	**D5 · Flat-resolvent leg-labeling — SOUND.** `FLAT_SYMBOL_USAGE`: REGRESSION == KERNEL_DIAGONAL = ω·ρ_m/q_out,
    46	RESIDUAL=0. `dtn_flat_symbol` read only @917 as the uniform-limit diagonal regression operand; no MATERIAL
    47	off-diagonal term consumes it.
    48	
    49	**D6 · μ_R,bg FORM ablation — SOUND.** `modulus_form`@962 maps `m1_profile → m1_profile²` (plus every jet) — a
    50	genuine FORM change; `MU_R_FORM_RESIDUAL` = 5.3 MB nonzero.
    51	
    52	**E · N6 representation-invariance — NOT SETTLED (not a confirmed defect).** `REP_INVARIANCE_RESIDUAL` does not
    53	syntactically vanish (279 genuine-coupling terms). BUT its MULTIGRADE = `{(1,0,1),(1,1,1)}` only — the leading
    54	O(εη) coupling grade (1,1,0) cancels exactly, so leading-order rep-invariance holds; remnant confined to the
    55	σ_W drain sector. Could not prove the σ_W remnant vanishes on this box (needs momentum-integral evaluation —
    56	matches c1's rep-invariance family DEFERRED to ≥64 GB). One-sided `FLIP_FACE_SLOPE` moves each route.
    57	
    58	**F · Uniform limit (astra's flagged concern) — SETTLED; not a defect.** `UNIFORM_LIMIT_OPERAND` not a
    59	visibly-zero payload (~16.6 KB), but the genuine transverse-trial-dependent coupling term has Integral integrand
    60	≡ 0 exactly in every block, all 4 cases + both densities. The only nonzero residue is the §3c bare open-slot
    61	bookkeeping `coeff·(δp_minus+δp_plus)·Test`, `coeff = iΛ_A·ε/(ωρ_mτ_A+iρ_m) ≠ 0` — the `−extract(SLAB)` open
    62	piece. So the genuine closure-induced coupling decouples (O(εη)→0). **Caveat:** the emitted uniform object is not
    63	literally zero — a step record must say "the genuine coupling decouples," not "the increment vanishes."
    64	
    65	**G · Tautology / adjointness — SOUND.** No `SELF_ENERGY_ADJOINTNESS_RESIDUAL` emitted; `CLOSED_COUPLING_KERNEL`
    66	carries both off-diagonal blocks built by the same weak restriction — the omission is honest (§3b), not a
    67	suppressed check. No emitted `*_RESIDUAL` has the `A−(A/B)·B` structural-zero form.
    68	
    69	## Overall verdict
    70	**The fold's physics is SOUND: zero confirmed physics defects.** Two items to carry, neither an error: (E) the
    71	N6 σ_W-sector remnant could not be settled to zero on this box — leading-order invariance holds, the remnant is
    72	consistent with the ≥64 GB deferral but a genuine σ_W pullback incompleteness cannot be ruled out here; (F) the
    73	uniform-limit object is non-vanishing purely from the §3c bare-open-slot representation, so interpretation must
    74	say "the genuine coupling decouples," not "the increment vanishes."
     1	I'll start by reading the review-legs skill and the governing physics, then write an independent derivation before touching the script's output.I’ll copy the files into a sandbox, then read the spec, imports, and fold machinery so I can derive the physics independently before using the script’s output.The load-bearing functions are next: extract, kernel bridge, fold construction, and the control/re-adjudication paths.I’ll inspect the imported objects and the stdout index first, then write an independent derivation against those inputs rather than the script’s printed values.I’ll write an independent derivation that inspects the imported objects and reconstructs the fold physics from those inputs.The first derivation hit a Tuple-shaped U-row; I’ll handle that and inspect where the mechanical pressure slots actually live.The uniform-limit operand still carries a bare `delta_p_minus` slot — I’ll parse that payload and finish the independent derivation.The response VALUE is nested; I’ll unwrap it, parse the uniform operand, and start a TRIAGE form ablation.The uniform increment is genuinely nonzero, and the two off-diagonal blocks are not adjoint. Next I’ll inspect the kernel legs and run a TRIAGE-scale form ablation.The plus/minus swap plus a jet-sign flip cancelled on this symmetric LAB_HELD case. I’ll apply a real form change: collapse the two q-legs and drop the jet substitution.Independent derivation scripts and stdout (do not treat the audit as a source of truth):
     2	
     3	| Role | Script | Stdout |
     4	|---|---|---|
     5	| Inverse, ε-powers, density live relation, traction covector, kernel legs | `/tmp/c2_review/independent_derive3.py` | `/tmp/c2_review/independent_derive3.out` |
     6	| FIRST_SHAPE atoms, all `face_velocity` cases, increment MULTIGRADE scan | `/tmp/c2_review/inspect_kernel_and_velocity.py` | `/tmp/c2_review/inspect_kernel_and_velocity.out` |
     7	| One-case fold inspect (maps, slots, grades) | `/tmp/c2_review/inspect_fold_once.py` | `/tmp/c2_review/inspect_fold_once.out` |
     8	| Uniform / zero-DtN operand parse | `/tmp/c2_review/parse_uniform_operand.py` | `/tmp/c2_review/parse_uniform_operand.out` |
     9	| FORM ablation (q-leg collapse + drop jets + flip pressure) | same inspector on `/tmp/c2_ablate_form/scripts/...` | `/tmp/c2_review/inspect_fold_form2.out` |
    10	
    11	---
    12	
    13	## A. Fold map — **sound on the construction; two inherited caveats**
    14	
    15	**Pressure, not `J_s`.** Independent inspect of `slab_operator` EXPANDED rows: θ and e_W carry `delta_p_±` / `d_w_delta_p_±` only; no `J_s` / `delta_j` / `relative_flux` slots (`independent_derive3.out`, `independent_derive2.out`). `build_face`@510–512 substitutes exactly those four slot names. One-case run: `SUBSTITUTION_KEYS = [d_w_delta_p_minus, d_w_delta_p_plus, delta_p_minus, delta_p_plus]`, `HAS_J_SLOT False`, `CLOSED_HAS_DELTA_P []` (slots fully replaced), `INC_HAS_J []` (`inspect_fold_once.out:5–6,33–38`). c1’s `J` and `T` leaves exist on the response and are not substituted.
    16	
    17	**w-jets are differentiated.** `build_face`@506–509: `sp.diff(extension * kernel, NORMAL)` then eval at the reference face. One-case: `JET_EQUALS_PRESSURE False`, `JET_ZERO False` for both faces. Not hand-typed.
    18	
    19	**Operator inverse, not scalar division.** Independent reconstruction from §1d:
    20	
    21	\[
    22	\bigl[I + (\Lambda_A/\rho_m^2)Z\bigr]\delta p = Z\cdot\text{source},\quad \Lambda_A=\Lambda_A^0/(1-i\omega\tau_A).
    23	\]
    24	
    25	Against every c1 `RESOLVENT_DEFINITION`: `INVERSE_OPERAND_RESIDUAL 0` (`independent_derive3.out:14–15` and all eight cases). `kernel_bridge`@396–398 evaluates that operand as a 2×2 matrix and `upper_triangular_solve`s it. Diagonal of `RESOLVENT_KERNEL_MATRIX` is `1/(1 + Λ_A ω /(ρ_m q_out (1−iωτ_A)))`, not a scalar `1/Z`.
    26	
    27	**Kernel bridge, both q-legs.** `dtn_kernel` FIRST_SHAPE carries `s11cc1_q_out_input` and `s11cc1_q_out_output` (`inspect_kernel_and_velocity.out:1–7`). Fold matrix is `[[z0(q_out), z1], [0, z0(q_in)]]` (`inspect_fold_once.out:12`). Increment has both `k_input` and `k_output`, no `dtn_operator` atoms (`INC_HAS_DTN_OP []`). Whole-form tags are the empty tuple (170–180 bytes). **Caveat (inherited, not a c2 freeze):** FIRST_SHAPE has `k_input` only, no `k_output`; `k·k′` is `|k_in|²`. That is the c1 kernel c2 consumes. Both *q* legs stay live; Fourier transfer is `Ŵ(k_out−k_in)`.
    28	
    29	**Symbol map.** Independently: `V_s` is S11c-a `face_velocity` (interfacial normal speed), `μ_θ` is `mu_theta_operator`. Script: `c1_V → face_velocity/ε`, `c1_μ → μ_θ/ε` (`build_face`@489–499). One-case IDENTIFICATIONS: `s11cc1_V_lab_held_plus → W_0 e_W_t / 2`. For LAB_HELD, `face_velocity = W_0 e_W_t ε / 2` on **both** faces (no `u`); MATERIAL adds `−σ u·∇w₁` (`inspect_kernel_and_velocity.out:26–33`). Same-sign ± is in the S11c-a export, not invented here. LAB_HELD `U_BODY_BALANCE` has no `δp` slots (inherited S11c-b); mechanical `δp` on LAB_HELD lives in e_W / face-force, consistent with V being thickness-only.
    30	
    31	**Disposition A:** fold map is physically the right substitution. Inherited: FIRST_SHAPE `k`-leg content; ± `face_velocity` sign.
    32	
    33	---
    34	
    35	## B. ε-normalization — **strip-one is correct; the increment is not purely O(εη)**
    36	
    37	Independent ε-powers (`independent_derive3.out`):
    38	
    39	- c1 `DELTA_P`: ε¹; strip one → ε⁰.
    40	- Slab slot coefficients: ε¹ on every `delta_p` / `d_w_delta_p` hit.
    41	- `face_velocity` and `mu_theta_operator`: ε¹.
    42	
    43	`build_face`@487,493–494 strips exactly one ε from the source and from the slab amplitudes. Naive double-ε does **not** happen.
    44	
    45	Emitted increment MULTIGRADE (byte-offset scan of the 2.57 MB tag):
    46	
    47	`(0,0,0), (1,0,0), (1,0,1), (1,1,0), (1,1,1)`
    48	
    49	i.e. it contains **O(ε) with no η**. Closed coupling kernel the same. Spec §4 / N12 says the self-energy increment lands at O(εη). The O(ε) piece is the local flat `Z_0` substitution into the open slots (same object as F). Not a double-ε bug; it is a grade of the increment the spec forbade.
    50	
    51	**Disposition B:** strip-one implemented correctly. **CONFIRMED:** increment is not O(εη)-only.
    52	
    53	---
    54	
    55	## C. Close-then-extract — **sound; commutator nonzero; FORM ablation bites**
    56	
    57	`build_case`@553–563: `closed = source.subs(substitutions)` then `extract(closed)`; increment = `extract(close)−extract(open)`. `run`@900–906: `close(extract)` is the imported open kernel with the same substitutions (`regression` only). `ORDERING_COMMUTATOR` is 431 555 bytes, MULTIGRADE includes `(1,0,1)` and `(1,1,1)` — not byte-identical, not Integer(0).
    58	
    59	FORM ablation (copy under `/tmp/c2_ablate_form` only): collapse `q_out_input → q_out_output` (`kernel_bridge`@378–380) **and** replace jets by 0 and flip pressure sign (`build_face`@510–512). Literal increment size:
    60	
    61	| block | baseline | FORM-ablated |
    62	|---|---|---|
    63	| TRANSVERSE→THICKNESS THETA | 268452 | **173410** |
    64	| TRANSVERSE→THICKNESS E_W | 536127 | **346325** |
    65	
    66	A first attempt (swap ± slots **and** flip jet sign) was a **no-op** on LAB_HELD because `P_+ ≈ P_−` and the two flips cancel; that is a symmetry of this case, not evidence of a typed payload. The q-collapse + dropped-jet ablation moves the object. Closure threading is load-bearing.
    67	
    68	**Disposition C:** `extract(close)` is what is built; commutator is genuinely nonzero.
    69	
    70	---
    71	
    72	## D. Six re-adjudications
    73	
    74	**D1. Density live-vs-frozen — sound.** Independent: `Eq(rho_br_bg_rho4_constant, W_bg ρ_br / W_0)` with residual 0 against `W_bg ρ_br / W_0`; RHOBR is the bare `ρ_br` (`independent_derive3.out:157–160`). Fold binds `{rho_br_bg_rho4_constant: ρ_br(1+η w₁)}` (`inspect_fold_once.out:10`). Live field, not a constant, not `ρ_m`, not a `∇ρ` jet (`DENSITY_IS_RHO_M False`). `DENSITY_LIVE_MINUS_FROZEN` is 1.5–1.9 MB with grades `(1,1,0)` and `(1,1,1)` — freeze is observable. `∇ρ→0` is not what is tested.
    75	
    76	**D2. `t_s` covector — sound.** S11c-a `traction` is a 4×1 matrix carrying `δp`, `Λ_X`, `d_w_delta_p` (`independent_derive3.out:184–189`). `traction_pairing`@772–792 uses that covector, not c1’s UNDECIDED `T` leaf. Native form.
    77	
    78	**D3. DtN whole-form — sound.** All four `DTN_WHOLEFORM_DEPENDENCE` tags are the empty tuple. Increment has no noncommutative `dtn_operator` atoms. Load-bearing object is the kernel matrix.
    79	
    80	**D4. Traction-vs-slab pairing — independent route; residual not Integer(0).** Algebra at `traction_pairing`@800–803:
    81	
    82	`KINETIC_STORED = SLAB_POWER − FACE_POWER` ⇒ emitted residual `= FACE_POWER − TRACTION_POWER`.
    83	
    84	That is two constructions (slab `FACE_GENERALIZED_FORCE_ROWS` vs S11c-a traction), not `A − (A/B)B`. One-sided `flip=True` only negates traction (`@778–779`). `TRACTION_SIGN_RESIDUAL` is 7.9 MB, `Add(...)`, grades `(1,0,0)…(1,1,1)` — it moves. Pairing residual itself is 11.6 MB `Add(...)`. Whether it *simplifies* to 0 is unsettled (too large to SZ here); it is not a structural zero of the tautology kind.
    85	
    86	**D5. Flat-symbol usage — sound.** `dtn_flat_symbol` appears only in `FLAT_SYMBOL_USAGE` (`run`@917–919). Increment has `INC_HAS_FLAT []`. Residual `KERNEL_DIAGONAL − REGRESSION = 0` with both sides `ω ρ_m / q_out`. No MATERIAL off-diagonal term consumes it.
    87	
    88	**D6. `μ_R` FORM — moved.** `modulus_form`@959–977 is `m1_profile → m1_profile²` including jets (not a rescale). `MU_R_FORM_RESIDUAL` is 5.3 MB. Control re-enters at the composed `μ_θ` operand (`control(..., mu_override=...)`@926).
    89	
    90	**Disposition D:** D1, D2, D3, D5 sound. D4 independent and the sign corruption bites; pairing residual not shown to vanish. D6 form control bites.
    91	
    92	---
    93	
    94	## E. N6 representation invariance — **CONFIRMED: residual does not vanish (σ channel)**
    95	
    96	`REP_INVARIANCE_RESIDUAL_RHO4` is 5 582 453 bytes, **not** Integer(0). MULTIGRADE is **only** `(1,0,1)` and `(1,1,1)` — O(εσ) and O(εησ). The O(ε) and O(εη) grades of the two increments are **absent** from the residual (they match). Operand sizes are unequal: LAB 2.57 MB vs mapped MATERIAL 4.74 MB.
    97	
    98	One-sided `FLIP_FACE_SLOPE` (`kernel_bridge`@385–388, only `face==1`): `REP_INDEPENDENCE_{MATERIAL,LAB}_RESIDUAL` are 1.24 / 1.20 MB with grades `(1,0,1),(1,1,1)`. Each route moves. They are separate `build_case` calls, so corrupting MATERIAL cannot move LAB.
    99	
   100	Spec §5c predicted the residual vanishes. It does not, at slope order. Either the Eulerian↔material pullback (`representation_pullback`@980–1015) is incomplete in σ, or the two anchorings differ at O(εσ). Not a typed-zero tag.
   101	
   102	**Disposition E:** **CONFIRMED** nonzero at O(εσ). O(ε) / O(εη) appear invariant.
   103	
   104	---
   105	
   106	## F. Uniform limit — **CONFIRMED: genuinely nonzero (astra’s concern is real)**
   107	
   108	Parsed `UNIFORM_LIMIT_OPERAND` (16 655 bytes, all four cases structurally the same):
   109	
   110	- **Not** Integer(0). No η, no σ, no `W_bg`.
   111	- TRANSVERSE→THICKNESS **THETA** and **E_W** nonzero; DIV_U and the entire THICKNESS→TRANSVERSE block **identically 0**.
   112	- Remaining **open-slot** symbols `delta_p_minus` and `delta_p_plus` (from `−extract(open)`), **plus** closed-response `Integral`s of `G_theta_u Λ_A0` × transverse trial curls (local `Z_0` acting on `μ_θ(u)`).
   113	
   114	Explicit THETA head (`parse_uniform_operand.out`):
   115	
   116	```text
   117	I*Λ_A0*δp_minus*ε*TestTheta/(ω ρ_m τ_A + I ρ_m)
   118	+ I*Λ_A0*δp_plus*ε*TestTheta/(ω ρ_m τ_A + I ρ_m)
   119	− 2I*Λ_A0*ε*TestTheta*Integral(ω ρ_m (G_theta_u Λ_A0 (curl A) / (ρ_br ρ_m …) + …))
   120	```
   121	
   122	`ZERO_DTN_OPERAND` (6 282 bytes) is only the leftover open slots (including η-weighted jets) with **no** Integrals — so the Integrals in the uniform operand are the closed flat response, not an unsimplified zero.
   123	
   124	This is the substitution increment `extract(close)−extract(open)` at η=σ=0, `W_bg→W_0`. It cannot vanish: close replaces δp by `Z_0·source(V,μ_θ)`, open still has the symbolic slot. Spec §5e predicted vanish (S11b decoupling). Two separate failures sit in the same object: (i) leftover open-slot terms, (ii) a surviving local `Z_0 μ_θ` transverse→θ / e_W coupling.
   125	
   126	**Disposition F:** **CONFIRMED** decoupling failure of *this* increment. Not a large structurally-zero expression.
   127	
   128	---
   129	
   130	## G. Tautology / operand-theatre
   131	
   132	- **§3c increment** is tagged `SELF_ENERGY_INCREMENT`, not `*_RESIDUAL`. Honest export.
   133	- **No** `SELF_ENERGY_ADJOINTNESS_RESIDUAL` in the 154-tag index. The builder’s claim that the blocks are adjoint by construction is **not honest**: for LAB_HELD, THICKNESS→TRANSVERSE increment is identically 0 while TRANSVERSE→THICKNESS is a 268 k-char expression (and the same split appears in the uniform operand). Those blocks are independently extracted (`extract`@326–343); they are not forced adjoint. Suppressing the residual hides a real asymmetry.
   134	- **Ordering commutator, routing, μ_R, traction-sign, N6 independence:** two-route or operand-minus-corrupted-operand; FORM / one-sided tests move them.
   135	- **FLAT_SYMBOL residual 0:** two c1 exports of the same diagonal formula; regression-only, as specified.
   136	- **Traction pairing residual** reduces to `FACE_POWER − TRACTION_POWER` (see D4); not a tautology.
   137	- No physics `assert` precedes a residual emit.
   138	
   139	**Disposition G:** increment tagging is honest. **CONFIRMED:** adjoint-by-construction is false.
   140	
   141	---
   142	
   143	## FORM / one-sided ablation summary
   144	
   145	| Ablation | What moved | Literal |
   146	|---|---|---|
   147	| Collapse `q_in→q_out` + drop jets + flip `δp` (FORM) | increment THETA/E_W sstr | 268452→173410 / 536127→346325 |
   148	| Swap ± slots + flip jet (LAB_HELD) | nothing (symmetry no-op) | sizes identical |
   149	| Density live vs frozen | `DENSITY_LIVE_MINUS_FROZEN` | 1.5 MB, grades O(εη) |
   150	| `t_s` sign flip, traction only | `TRACTION_SIGN_RESIDUAL` | 7.9 MB `Add` |
   151	| `FLIP_FACE_SLOPE` one route | independence residuals | 1.2 MB, grades O(εσ) |
   152	| `μ_R` profile `m→m²` | `MU_R_FORM_RESIDUAL` | 5.3 MB |
   153	| Uniform η=σ=0 | increment | 16 kB, explicit δp + Integrals |
   154	
   155	---
   156	
   157	## One-line dispositions
   158	
   159	| Item | Disposition |
   160	|---|---|
   161	| **A** | Fold map sound (pressure not `J_s`; jets differentiated; operator inverse; kernel both q-legs; `V_s→face_velocity`). |
   162	| **B** | ε-strip-one correct. **CONFIRMED:** increment MULTIGRADE includes O(ε), not only O(εη). |
   163	| **C** | `extract(close)` is built; commutator nonzero; FORM ablation of threading/q-legs moves the increment. |
   164	| **D** | D1–D3, D5 sound. D4 independent, sign corruption bites. D6 form control bites. |
   165	| **E** | **CONFIRMED:** N6 residual does not vanish; leftover is the σ grades. One-sided slope test moves one route. |
   166	| **F** | **CONFIRMED:** uniform-limit increment is genuinely nonzero (open-slot `δp_±` + local `Z_0 μ_θ` Integrals). |
   167	| **G** | Increment is an export, not a dressed residual. **CONFIRMED:** off-diagonal blocks are not adjoint; the suppressed check would have shown it. |
   168	
   169	**Overall.** The fold *wiring* (slots, inverse, kernel, density live-binding, ε-strip, no `J_s`) is physically the right construction. The *increment as defined* is not the S11b-decoupled self-energy the spec claims: it keeps an O(ε) local-response piece, it does not vanish as η→0, the N6 residual fails at O(εσ), and the two off-diagonal blocks are not adjoint. What must change is the increment/control objects (subtract the open-slot leftover and/or test `extract(close)|_{η=0}` for decoupling; repair or surface the σ-channel N6 residual; stop claiming adjointness), not the substitution of closed `δp` into the slab slots.
     1	# S11c-c2 SymPy builder report
     2	
     3	Status: the full detached production audit and export publication completed. Mechanical checks and scoped reduction-tool evidence are recorded below; physics dispositions remain for the step record and independent reviews.
     4	
     5	## Authority and boundary
     6	
     7	`directives/S11c_c2_SHARED_PHYSICS.md` matches the bytes at commit `16849fc6`. All of spec §1, its AGREE/UNDECIDED dispositions, the inherited weak-restriction convention, the outgoing-bulk ansatz, the profile/density relations, and the retained grading are SUPPLIED and unfalsifiable within this build. Running this script does not independently validate those inputs. The face-force convention and the #90 closure-fold convention remain in the inherited slot carriers. No cross-engine normalization or disposition is made here.
     8	
     9	The requested implementation is `scripts/S11c_c2_selfenergy_fold_sympy_audit.py`; the generated delta is `scripts/S11c_c2_exports.py`. No Wolfram engine, comparator, reviewer, additional builder, or commit was launched. The parent owns the fresh Claude and Grok reviews.
    10	
    11	## Frozen symbol map
    12	
    13	The source is the `DELTA_P` quantity under `s11c_c1_face_response / CASES / (anchoring, face, density) / VALUE`. The build checks that the two prefixed response roots have producer `S11c-c1`. Neither bare S11b response key is bound. The coefficient export's `PRESSURE` quantity is not a pressure construction operand.
    14	
    15	The face representation is `DELTA_W` for both faces of both anchorings. `build_face` (line 477) applies these stages:
    16	
    17	1. Bind `rho_br_bg_rho4_constant` to the matching `background_density_map` VALUE before composing the pressure source. The live representative is the imported field relation; the other representative is independently read from its own density-map case.
    18	2. Read `DELTA_P`, divide its source by `epsilon_shape`, and evaluate its ordered resolvent–DtN product using the kernel bridge. The slab slot already supplies its external epsilon.
    19	3. Map each face's c1 `mu_theta` identifier to the matching anchoring/density `mu_theta_operator` amplitude. Both faces use the anchoring's single constitutive operand. That imported operand itself contains epsilon, so its amplitude is obtained by division by epsilon.
    20	4. Map each c1 velocity identifier to the matching `face_velocity / (anchoring, face, DELTA_W)` amplitude, also divided by its imported epsilon. `kinematic_balance` is not used as the velocity operand.
    21	5. Map the c1 whole-form `Z` identifier and resolvent identifier to the computed kernel-algebra matrices. The construction uses `dtn_kernel`, including its two legs. The raw `dtn_operator` is not silently substituted as the construction operator.
    22	6. Compute the normal jet by differentiating the outgoing continuation ansatz with respect to the actual normal coordinate and evaluating at the reference face. Substitute both the pressure and normal-jet slots per face. No extra closed `J` is added.
    23	
    24	Every case emits its full pressure, normal jet, density substitution, four identifier-class bindings, inverse operand, inverse matrix, pressure matrix, and second-scattering kernel under `PY_S11CC2_FOLD_SYMBOL_MAP_*`.
    25	
    26	Exact application lines in the final script: live-density binding **482–484**; c1 epsilon normalization **487**; mu amplitude **493**, velocity amplitude **494**, and their simultaneous identification **499–500**; the kernel bridge **366–417**, including inverse evaluation **396–398** and ordered intermediate-leg composition **402–416**; pressure integral application **501**; normal differentiation **506–509**; slot map **510–512**; substitution into the full rows **553–555**. `REPRESENTATION = 'DELTA_W'` is fixed at line 47. All eight source cases are bound explicitly by `build_face`; the table below spells out the identifier names. Both faces of each anchoring map to that anchoring's single slab mu operand, resolved through the matching density case of `mu_theta_operator`.
    27	
    28	| Anchoring / density / face | c1 mu identifier → slab mu | c1 velocity identifier → face-velocity case | Whole-form identifier / resolvent identifier → kernel bridge |
    29	|---|---|---|---|
    30	| `LAB_HELD / RHO4_CONSTANT / 1` | `s11cc1_mu_theta_lab_held_plus` → `mu_theta_L` | `s11cc1_V_lab_held_plus` → `(LAB_HELD, 1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_plus` / `s11cc1_response_resolvent_lab_held_plus_rho4_constant` → computed per-case matrices |
    31	| `LAB_HELD / RHO4_CONSTANT / -1` | `s11cc1_mu_theta_lab_held_minus` → `mu_theta_L` | `s11cc1_V_lab_held_minus` → `(LAB_HELD, -1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_minus` / `s11cc1_response_resolvent_lab_held_minus_rho4_constant` → computed per-case matrices |
    32	| `LAB_HELD / RHOBR_CONSTANT / 1` | `s11cc1_mu_theta_lab_held_plus` → `mu_theta_L` | `s11cc1_V_lab_held_plus` → `(LAB_HELD, 1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_plus` / `s11cc1_response_resolvent_lab_held_plus_rhobr_constant` → computed per-case matrices |
    33	| `LAB_HELD / RHOBR_CONSTANT / -1` | `s11cc1_mu_theta_lab_held_minus` → `mu_theta_L` | `s11cc1_V_lab_held_minus` → `(LAB_HELD, -1, DELTA_W)` | `s11cc1_dtn_operator_lab_held_minus` / `s11cc1_response_resolvent_lab_held_minus_rhobr_constant` → computed per-case matrices |
    34	| `MATERIAL_ADVECTED / RHO4_CONSTANT / 1` | `s11cc1_mu_theta_material_advected_plus` → `mu_theta_M` | `s11cc1_V_material_advected_plus` → `(MATERIAL_ADVECTED, 1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_plus` / `s11cc1_response_resolvent_material_advected_plus_rho4_constant` → computed per-case matrices |
    35	| `MATERIAL_ADVECTED / RHO4_CONSTANT / -1` | `s11cc1_mu_theta_material_advected_minus` → `mu_theta_M` | `s11cc1_V_material_advected_minus` → `(MATERIAL_ADVECTED, -1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_minus` / `s11cc1_response_resolvent_material_advected_minus_rho4_constant` → computed per-case matrices |
    36	| `MATERIAL_ADVECTED / RHOBR_CONSTANT / 1` | `s11cc1_mu_theta_material_advected_plus` → `mu_theta_M` | `s11cc1_V_material_advected_plus` → `(MATERIAL_ADVECTED, 1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_plus` / `s11cc1_response_resolvent_material_advected_plus_rhobr_constant` → computed per-case matrices |
    37	| `MATERIAL_ADVECTED / RHOBR_CONSTANT / -1` | `s11cc1_mu_theta_material_advected_minus` → `mu_theta_M` | `s11cc1_V_material_advected_minus` → `(MATERIAL_ADVECTED, -1, DELTA_W)` | `s11cc1_dtn_operator_material_advected_minus` / `s11cc1_response_resolvent_material_advected_minus_rhobr_constant` → computed per-case matrices |
    38	
    39	The exact code operations are in `build_face`, `kernel_bridge`, and `kernel_apply`; the source index records their spans.
    40	
    41	## Computation and representation
    42	
    43	| Object / operation | Computation location |
    44	|---|---|
    45	| Input binding / producer verification | `bind_inputs`, line 177 |
    46	| Generated-coefficient dimensions | `infer_dimensions`, line 212 |
    47	| Whole-field trial/jet ansatz | `wave_jet`, line 140 |
    48	| Continuous source-point field map | `at_source`, line 246 |
    49	| Two-leg inverse and intermediate-momentum composition | `kernel_bridge`, line 366 |
    50	| Solved outgoing branch bindings | `outgoing_spectral`, line 438 |
    51	| Fourier profile bindings | `profile_bindings`, line 732 |
    52	| Closed pressure and computed normal jets | `build_face`, line 477 |
    53	| Close before weak extraction | `build_case`, line 532 |
    54	| Both weak off-diagonal blocks | `extract`, line 326 |
    55	| Same-extract subtraction | `difference`, line 100 |
    56	| Rectangular retained-grade computation | `shape_coefficients`, line 579 |
    57	| Native traction / independent face-work pairing | `traction_pairing`, line 772 |
    58	| All import-level controls | `control`, line 765 |
    59	| Modulus-profile form control | `modulus_form`, line 959 |
    60	| Anchoring field pullback and jets | `representation_pullback`, line 980 |
    61	| Imported coupling regression coordinate bridge | `regression_coordinates`, line 1018 |
    62	| Own-row bind-closure publication | `publish`, line 807 |
    63	
    64	The row representation contains the assembled strong `U`, `THETA`, and `E_W` balances, retaining all derivatives in the imported expanded rows. Per-face pressure-slot contributions are computed from those same rows. The reported parity quantities are the sum and difference of these per-face contributions. The full imported open coupling root is used only in the ordering regression.
    65	
    66	The continuous kernels act through explicit Fourier integrals on arbitrary source fields. Their integrands and ordered kernel compositions are computed. The arbitrary profile/test fields are not assigned a profile or evaluated numerically. The outgoing roots are solved with SymPy and retained as named branch functions with their computed definitions attached to every exported primary. The profile Fourier bindings are also attached. These bindings are part of the serialized object and are required when interpreting or comparing the named functions.
    67	
    68	The inverse is evaluated through the inherited rectangular retained grade: the two-leg triangular solve supplies the diagonal and first-scattering pieces; a three-leg triangular solve supplies the mixed second-scattering term, with the intermediate momentum integrated explicitly. Higher shape grades are outside this representation. The normal extension used for the w-jet is explicit in `build_face`. This is a retained-order inverse kernel, not an all-orders inversion on a prescribed numerical background. The exported named branch functions and profile transforms have attached computed definitions; their defining integrals remain on arbitrary fields.
    69	
    70	There is no independent adjointness construction in this implementation, so no `SELF_ENERGY_ADJOINTNESS_RESIDUAL` is emitted. Both off-diagonal blocks are emitted. The traction pairing compares the slab's inherited generalized-force contribution to independently composed native covector work; its kinetic/stored operand is the full slab mechanical power with the inherited face generalized-force rows removed. The test velocities use the imported velocity amplitudes, with the external wave epsilon divided out consistently with the pressure map; physical quadratic mechanical power restores that test-amplitude epsilon. This is a comparison of face work, not an independent derivation of all kinetic and stored energies. It does not adjudicate c1's separate far-field ENERGY audit.
    71	
    72	The N6 comparison constructs both imported anchorings and applies the density/thickness field map to the full material rows before re-extraction. Its two one-sided probes each read the anchoring's native face-normal covector, flip the upper face's first slope component, and propagate the computed component ratio into the corresponding first Fourier tilt of that anchoring's imported DtN kernel. The original and corrupted native normals are emitted with each probe. The companion anchoring is emitted from its separately constructed baseline; no difference of an object with itself is presented as a companion test. These precise probe scopes are exposed for review; they are not a claim that the unresolved c1 density, traction, whole-form, leg-labeling, or ENERGY dispositions have closed.
    73	
    74	## Import guard evidence
    75	
    76	The positional two-parent fold measured `2485` rows, with source counts `[['/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_b_exports.py', 2441], ['/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c1_exports.py', 44]]` and overwrites `[]`. The actual access-recording call is `assert_lookups_equal_manifest(bind_inputs, fold, IMPORT_KEYS)`, and its observed lookup set equals the manifest on these files. `check_consumer` resolves `543` keys. The literal evidence is `_measurements/S11c_c2_sympy_guard_evidence.json`.
    77	
    78	Final `IMPORT_KEYS` (34 keys):
    79	
    80	- `L_W`
    81	- `Lambda_A_0`
    82	- `Lambda_V_0`
    83	- `Lambda_X_0`
    84	- `W_0`
    85	- `W_bg`
    86	- `background_density_map`
    87	- `c_s0`
    88	- `closure_shape_deriv`
    89	- `conormal_deriv`
    90	- `coupling_kernel`
    91	- `dtn_flat_symbol`
    92	- `dtn_kernel`
    93	- `dtn_operator`
    94	- `epsilon_shape`
    95	- `eta_bg`
    96	- `face_measure_shape_deriv`
    97	- `face_normal`
    98	- `face_shift`
    99	- `face_velocity`
   100	- `kinematic_balance`
   101	- `mu_R`
   102	- `mu_theta_operator`
   103	- `omega`
   104	- `relative_flux`
   105	- `rho_br`
   106	- `rho_br_bg_rho4_constant`
   107	- `rho_m`
   108	- `s11c_c1_face_response`
   109	- `s11c_c1_face_response_coeffs`
   110	- `sigma_W`
   111	- `slab_operator`
   112	- `slab_operator_term_origins`
   113	- `traction`
   114	
   115	## Export and measurement evidence
   116	
   117	The declared outgoing roots are `s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, and `s11cc2SelfEnergyIncrement`. Publication computes their own bind-closure over the real parent fold, adds referenced new coordinate/function declarations and their dimension rows, checks exact-key collisions, runs `assert_delta_is_minimal`, and checks serialization against the live values. The generated module pins the five required input digests. The runtime evidence records the exact exported key set and guard operands.
   118	
   119	Production stdout is `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`. Each §5 control has its own literal operand and residual tags there. No residual is asserted against a physics target. The completion marker is `_measurements/S11c_c2_sympy_completion.txt`; `/usr/bin/time -v` measures the Python process in `_measurements/S11c_c2_sympy_resources.txt`; `_measurements/S11c_c2_sympy_monitor.tsv` samples that process's PID, elapsed seconds, RSS, VSZ, and CPU percentage. `_measurements/S11c_c2_sympy_progress.json` records the current case/control.
   120	
   121	`_measurements/S11c_c2_sympy_source_index.json` records the script hash, line count, Python lexical-token count, and exact function spans. The runtime tag index provides the exact computation/emit location and fresh write-key for each emitted object. Prior development attempts and the single-case smoke run are separate `/tmp/S11c_c2_attempt*` / `/tmp/s11cc2_smoke*` files; they are not the production completion evidence.
   122	
   123	## Final measured evidence
   124	
   125	The production completion marker contains `0`. The audit contains 1053 lines, 82,702 bytes, and 22,953 Python lexical tokens (a source-size plausibility measure, not a model-token billing count). Its SHA-256 is `6360721815b0b8c4a52393c388c37b20c858bca04615d8e671964af49a87c7dd`. The plain, non-symlink export contains 71 own rows and 60,516,900 bytes. Publication reloaded the serialized module and emitted its literal equality tuple before the round-trip guard.
   126	
   127	The actual lookup manifest has 34 keys, exactly matching the recorded accesses; its consumer closure has 543 keys. The three-parent fold has 2556 rows. Own delta membership equals the computed own bind-closure, imported-key collisions are `[]`, and the independent file reread records all five digest comparisons in `_measurements/S11c_c2_postbuild_verify.json`. Closed primary pressure/normal-slot references are `[]`; raw c1 whole-form/resolvent references are `[]`.
   128	
   129	The stdout index records 154 unique tags across 498,811,405 bytes, duplicate tags `[]`, and 0 untagged records. Every record carries a computed multigrade and dimension object. The complete literal metadata, branch definitions, and residual values remain in the indexed stdout; they are not converted into physics dispositions here.
   130	
   131	The streaming dimension-metadata scan found missing dimension fields on `[]` and literal `nan` entries on `[]`. This records the dimension walker's output; the independent corruption probe below tests its sensitivity.
   132	
   133	Python PID(s) sampled: `['2439638']`. Runtime self-measurement: 1159.807 seconds and 2,689,720 KiB peak RSS. All requested in-script controls ran serially; none was deferred. The inherited c1 giant-family and ≥64 GB cross-engine comparisons remain outside this build. The final `/usr/bin/time -v` record is:
   134	
   135	```text
   136		Command being timed: "python -u /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py"
   137		User time (seconds): 1159.48
   138		System time (seconds): 2.33
   139		Percent of CPU this job got: 99%
   140		Elapsed (wall clock) time (h:mm:ss or m:ss): 19:22.18
   141		Average shared text size (kbytes): 0
   142		Average unshared data size (kbytes): 0
   143		Average stack size (kbytes): 0
   144		Average total size (kbytes): 0
   145		Maximum resident set size (kbytes): 2689720
   146		Average resident set size (kbytes): 0
   147		Major (requiring I/O) page faults: 5
   148		Minor (reclaiming a frame) page faults: 1784272
   149		Voluntary context switches: 24
   150		Involuntary context switches: 14767
   151		Swaps: 0
   152		File system inputs: 184
   153		File system outputs: 1093544
   154		Socket messages sent: 0
   155		Socket messages received: 0
   156		Signals delivered: 0
   157		Page size (bytes): 4096
   158		Exit status: 0
   159	```
   160	

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/_measurements/S11c_c2_physics_grok_leg.txt | sed -n '124,172p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py | sed -n '1,180p' && nl -ba research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py | sed -n '1,180p'" in /var/projects/toy_physics
 succeeded in 0ms:
   124	This is the substitution increment `extract(close)−extract(open)` at η=σ=0, `W_bg→W_0`. It cannot vanish: close replaces δp by `Z_0·source(V,μ_θ)`, open still has the symbolic slot. Spec §5e predicted vanish (S11b decoupling). Two separate failures sit in the same object: (i) leftover open-slot terms, (ii) a surviving local `Z_0 μ_θ` transverse→θ / e_W coupling.
   125	
   126	**Disposition F:** **CONFIRMED** decoupling failure of *this* increment. Not a large structurally-zero expression.
   127	
   128	---
   129	
   130	## G. Tautology / operand-theatre
   131	
   132	- **§3c increment** is tagged `SELF_ENERGY_INCREMENT`, not `*_RESIDUAL`. Honest export.
   133	- **No** `SELF_ENERGY_ADJOINTNESS_RESIDUAL` in the 154-tag index. The builder’s claim that the blocks are adjoint by construction is **not honest**: for LAB_HELD, THICKNESS→TRANSVERSE increment is identically 0 while TRANSVERSE→THICKNESS is a 268 k-char expression (and the same split appears in the uniform operand). Those blocks are independently extracted (`extract`@326–343); they are not forced adjoint. Suppressing the residual hides a real asymmetry.
   134	- **Ordering commutator, routing, μ_R, traction-sign, N6 independence:** two-route or operand-minus-corrupted-operand; FORM / one-sided tests move them.
   135	- **FLAT_SYMBOL residual 0:** two c1 exports of the same diagonal formula; regression-only, as specified.
   136	- **Traction pairing residual** reduces to `FACE_POWER − TRACTION_POWER` (see D4); not a tautology.
   137	- No physics `assert` precedes a residual emit.
   138	
   139	**Disposition G:** increment tagging is honest. **CONFIRMED:** adjoint-by-construction is false.
   140	
   141	---
   142	
   143	## FORM / one-sided ablation summary
   144	
   145	| Ablation | What moved | Literal |
   146	|---|---|---|
   147	| Collapse `q_in→q_out` + drop jets + flip `δp` (FORM) | increment THETA/E_W sstr | 268452→173410 / 536127→346325 |
   148	| Swap ± slots + flip jet (LAB_HELD) | nothing (symmetry no-op) | sizes identical |
   149	| Density live vs frozen | `DENSITY_LIVE_MINUS_FROZEN` | 1.5 MB, grades O(εη) |
   150	| `t_s` sign flip, traction only | `TRACTION_SIGN_RESIDUAL` | 7.9 MB `Add` |
   151	| `FLIP_FACE_SLOPE` one route | independence residuals | 1.2 MB, grades O(εσ) |
   152	| `μ_R` profile `m→m²` | `MU_R_FORM_RESIDUAL` | 5.3 MB |
   153	| Uniform η=σ=0 | increment | 16 kB, explicit δp + Integrals |
   154	
   155	---
   156	
   157	## One-line dispositions
   158	
   159	| Item | Disposition |
   160	|---|---|
   161	| **A** | Fold map sound (pressure not `J_s`; jets differentiated; operator inverse; kernel both q-legs; `V_s→face_velocity`). |
   162	| **B** | ε-strip-one correct. **CONFIRMED:** increment MULTIGRADE includes O(ε), not only O(εη). |
   163	| **C** | `extract(close)` is built; commutator nonzero; FORM ablation of threading/q-legs moves the increment. |
   164	| **D** | D1–D3, D5 sound. D4 independent, sign corruption bites. D6 form control bites. |
   165	| **E** | **CONFIRMED:** N6 residual does not vanish; leftover is the σ grades. One-sided slope test moves one route. |
   166	| **F** | **CONFIRMED:** uniform-limit increment is genuinely nonzero (open-slot `δp_±` + local `Z_0 μ_θ` Integrals). |
   167	| **G** | Increment is an export, not a dressed residual. **CONFIRMED:** off-diagonal blocks are not adjoint; the suppressed check would have shown it. |
   168	
   169	**Overall.** The fold *wiring* (slots, inverse, kernel, density live-binding, ε-strip, no `J_s`) is physically the right construction. The *increment as defined* is not the S11b-decoupled self-energy the spec claims: it keeps an O(ε) local-response piece, it does not vanish as η→0, the N6 residual fails at O(εσ), and the two off-diagonal blocks are not adjoint. What must change is the increment/control objects (subtract the open-slot leftover and/or test `extract(close)|_{η=0}` for decoupling; repair or surface the σ-channel N6 residual; stop claiming adjointness), not the substitution of closed `δp` into the slab slots.
     1	"""Orchestrator's INDEPENDENT adjudication of section F (uniform limit).
     2	Question: in the uniform limit, does the genuine closure-induced coupling
     3	[= extract(close)|_uniform] decouple to zero, leaving only the -extract(open)
     4	bare-slot bookkeeping? (fresh-Claude leg: YES; Grok: NO, a Z0.mu_theta Integral survives.)
     5	
     6	Method (independent of the leg's 'Trial'-label heuristic): the increment =
     7	extract(close) - extract(open). extract(close) substitutes ALL delta_p slots, so
     8	it carries NO bare delta_p_{plus,minus}/d_w_delta_p symbols; every bare-delta_p term
     9	therefore comes only from -extract(open). Zero those bare slots -> what remains is
    10	extract(close)|_uniform. .doit() the Integrals and test == 0.
    11	"""
    12	import sys; sys.path.insert(0,'/tmp/c2_ablate_agent')
    13	from load_tag import value_of, to_dict
    14	import sympy as sp
    15	
    16	BARE = {'delta_p_minus','delta_p_plus','d_w_delta_p_minus','d_w_delta_p_plus'}
    17	
    18	def bare_zero_subs(e):
    19	    sub={}
    20	    # zero any Symbol or AppliedUndef whose name is a bare open slot
    21	    for a in e.atoms(sp.Symbol):
    22	        if a.name in BARE: sub[a]=0
    23	    for a in e.atoms(sp.core.function.AppliedUndef):
    24	        if a.func.__name__ in BARE: sub[a]=0
    25	    # also symbols/functions carrying the bare name as prefix (jets like d_w_delta_p_plus)
    26	    for a in e.atoms(sp.Symbol):
    27	        if any(a.name==b for b in BARE): sub[a]=0
    28	    return e.subs(sub) if sub else e
    29	
    30	for tag in ['PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT',
    31	            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT',
    32	            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT',
    33	            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT']:
    34	    d=to_dict(value_of(tag)); print('====',tag)
    35	    for outer in d:
    36	        if not isinstance(d[outer],dict): continue
    37	        for inner in d[outer]:
    38	            e=sp.expand(d[outer][inner])
    39	            if e==0:
    40	                print('   [%s/%s] block ZERO'%(outer,inner)); continue
    41	            # names present, and whether Integrals appear
    42	            names=sorted({a.name for a in e.atoms(sp.Symbol) if a.name in BARE})
    43	            n_int=len(e.atoms(sp.Integral))
    44	            closed_part=bare_zero_subs(e)
    45	            closed_part=sp.expand(closed_part)
    46	            # try to evaluate any surviving integrals
    47	            try:
    48	                closed_doit=sp.expand(closed_part.doit())
    49	            except Exception as ex:
    50	                closed_doit=('doit-failed:%s'%ex)
    51	            is_zero = (closed_part==0)
    52	            is_zero_doit = (closed_doit==0) if not isinstance(closed_doit,str) else 'n/a'
    53	            print('   [%s/%s] full_block_zero=%s bare_slots=%s n_Integral=%d  '
    54	                  'closed_part_zero(expand)=%s  closed_part_zero(doit)=%s'
    55	                  %(outer,inner, e==0, names, n_int, is_zero, is_zero_doit))
    56	            if not is_zero:
    57	                # show a short head of what survives, to characterize it
    58	                surv = closed_part if not isinstance(closed_doit,str) else closed_part
    59	                head = sp.srepr(surv)[:180]
    60	                print('        SURVIVING(closed_part) head:', str(surv)[:200])
     1	"""Orchestrator's independent checks for E (N6) and G (adjointness/directionality)."""
     2	import sys; sys.path.insert(0,'/tmp/c2_ablate_agent')
     3	from load_tag import value_of, to_dict
     4	import sympy as sp
     5	
     6	def blocks(tag):
     7	    """Return {outer/inner: expr} for a nested block object."""
     8	    d=to_dict(value_of(tag)); out={}
     9	    if isinstance(d,dict):
    10	        for o in d:
    11	            if isinstance(d[o],dict):
    12	                for i in d[o]: out['%s/%s'%(o,i)]=d[o][i]
    13	            else: out[o]=d[o]
    14	    return out
    15	
    16	# ---------- G: block directionality of the self-energy (full object, not uniform limit) ----------
    17	print('==== G: SELF_ENERGY_INCREMENT block directionality (LAB_HELD_RHO4) ====')
    18	inc=blocks('PY_S11CC2_SELF_ENERGY_INCREMENT_LAB_HELD_RHO4_CONSTANT')
    19	for k in sorted(inc):
    20	    e=inc[k]
    21	    z=(sp.expand(e)==0)
    22	    print('   %-32s zero=%s  bytes~%d'%(k, z, len(sp.srepr(e))))
    23	
    24	print('==== G: CLOSED_COUPLING_KERNEL block directionality (LAB_HELD_RHO4) ====')
    25	ker=blocks('PY_S11CC2_CLOSED_COUPLING_KERNEL_LAB_HELD_RHO4_CONSTANT')
    26	for k in sorted(ker):
    27	    e=ker[k]
    28	    z=(sp.expand(e)==0)
    29	    print('   %-32s zero=%s  bytes~%d'%(k, z, len(sp.srepr(e))))
    30	
    31	# ---------- E: does the N6 rep-invariance residual carry ONLY sigma_W-sector terms? ----------
    32	# If residual.subs(sigma_W->0) == 0, the entire remnant carries sigma_W => leading O(eps),O(eps.eta)
    33	# rep-invariance HOLDS and the non-invariance is confined to the sigma_W drain sector.
    34	print('==== E: REP_INVARIANCE_RESIDUAL sigma_W-sector confinement (RHO4) ====')
    35	for tag in ['PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHO4_CONSTANT',
    36	            'PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHOBR_CONSTANT']:
    37	    try:
    38	        b=blocks(tag)
    39	    except KeyError:
    40	        print('   MISSING tag',tag); continue
    41	    # find the sigma_W symbol among atoms
    42	    allexpr=sp.Add(*[v for v in b.values()]) if b else sp.S.Zero
    43	    sig=[a for a in allexpr.atoms(sp.Symbol) if a.name in ('sigma_W','sigmaW','sigma_w')]
    44	    print('   ',tag,'blocks=%d sigma_W_symbol=%s'%(len(b), sig))
    45	    for k in sorted(b):
    46	        e=sp.expand(b[k])
    47	        if e==0:
    48	            print('      %-32s zero'%k); continue
    49	        e0 = e.subs({s:0 for s in sig}) if sig else e
    50	        e0 = sp.expand(e0)
    51	        print('      %-32s full_zero=%s  sigmaW->0_zero=%s'%(k, e==0, e0==0))

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md | sed -n '1,115p' && rg -n 'N12|sigma_W|σ_W|drain-projection|rep-invariance|representation-invariance' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md DEFERRED_HEAVY_RUNS.md research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md 2>/dev/null | sed -n '1,300p' && rg -n 'dissipative directional|dissipative|directional' STATUS.md research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md" in /var/projects/toy_physics
 succeeded in 0ms:
     1	# S11c-c2 — SHARED PHYSICS (the self-energy fold: closing the slab operator with the curved-bulk response)
     2	
     3	**S11c-c2** is the second half of the S11c-c curved-interface bulk closure (the S11c-c decision-list row was split
     4	c1/c2 by user choice 2026-09-03; `directives/S11c_decisions.md` N2). c1 solved the curved two-face outgoing bulk
     5	problem and **exported the closed permeable face response** `(δp_s, J_s, t_s)(V_s, μ_θ)`, the nonlocal DtN/impedance
     6	operator, its flat symbol, and its two-momentum kernel (`steps/S11c_c1_curved_bulk_closure.md`; exports
     7	`scripts/S11c_c1_exports.py`). **c2 folds that closed response into the S11c-b variable-coefficient slab operator
     8	`S11CB_SLAB_OPERATOR`** — whose θ-row and mechanical rows still carry the face pressure `δp_s` **symbolically** — and
     9	**re-extracts the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the CLOSED full operator**, yielding the
    10	coupled **nonlocal self-energy operator**: the transverse sector now carries a self-energy threaded through the bulk
    11	DtN. This document is the physics authority for the two blind c2 engines and their comparator. Tag prefix `S11CC2_`.
    12	
    13	The SymPy engine reads the inherited model through `ledger_fold.load_model` over the atomic frozen base
    14	`scripts/S11c_b_exports.py` **with the c1 delta `scripts/S11c_c1_exports.py`** folded on top (§7), binding only its
    15	declared `IMPORT_KEYS`; the Wolfram engine imports nothing and re-derives every consumed object from the sibling
    16	specs (`S9_export_chain_rebuild_directive.md:16-18` is the only cross-engine control). Blindness is the control: an
    17	agreement is independent construction, not a copy.
    18	
    19	⭐ This is a **CODE build** authored against c1's and S11c-b's reviewed exports. Rule 7: it gets **two decision legs
    20	(Codex `gpt-5.6-sol` xhigh + Grok) before any builder**; the build directive that follows gets its own two decision
    21	legs before the build. ⚠ **This is spec v2** — folded once from the v1 two-leg decision gate
    22	(`directives/_measurements/S11c_c2_spec_review.md`), which found the v1 isolation claim false, the fold operation
    23	under-specified against the real θ-row, and the re-adjudication set incomplete.
    24	
    25	---
    26	
    27	## 0 · Scope
    28	
    29	**In scope.** (1) **Close** the S11c-b slab operator by substituting the c1 closed face response into its symbolic
    30	`δp_s`-slots, per anchoring `α∈{L,M}` and density representative `ρ∈{ρ_4D,ρ_br}`, summing the two faces `s∈{+,−}`
    31	that the slab EOM already sums (§1c). (2) **Re-extract** the off-diagonal transverse↔`{θ,e_W,u_L}` coupling from the
    32	**closed** full operator by the S11c-b §3c weak variational restriction — the **close-then-extract** ordering (§2).
    33	(3) Emit the resulting **nonlocal self-energy** as the **substitution increment** (§3c), both operands re-extracted
    34	from `SLAB_OPERATOR` with the same extract. (4) **Re-adjudicate the c1 items the fold makes load-bearing** (§3d):
    35	the background-density field-vs-constant freeze (c1 seal 5, rule 17), the `t_s` traction representation, the DtN
    36	whole-form-vs-kernel, the **traction-vs-slab mechanical-power pairing** c1 assigned to c2, the flat-resolvent
    37	leg-labeling, and the `μ_R,bg` form control c1 reserved for c2.
    38	
    39	**Out of scope (named, not solved).** The profile-conditioned spectrum/scattering and leakage rates (S11c-d, `N5`);
    40	the leakage observable and falsification (S11c-e, `N7`); the nonlinear-light program (`N10`); the O(η²) evanescent-
    41	leakage sign on the flat nullspace (S11c-e, c1 §3b caveat). c2 introduces **no** convective bulk operator (`N11a`).
    42	The full per-family cross-engine residual and the four c1 giant families remain **deferred to a ≥64 GB box**
    43	(`DEFERRED_HEAVY_RUNS.md`); c2 must be constructible and cross-engine-testable on this box for its own increment
    44	(§7), and name — not silently absorb — anything it cannot close here.
    45	
    46	---
    47	
    48	## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build
    49	
    50	Everything in §1 is an input. The substitutions, orderings, extractions, self-energy structure, and the six
    51	re-adjudication verdicts of §§3–5 are **outputs**; ⛔ none is stated here.
    52	
    53	### 1a · Inheritance and the two imported models
    54	
    55	The DOFs, sector split, background ansatz, `(ε,η,σ_W)` power counting (`N12`), and admissibility are exactly S11c-a
    56	§§1–2 / S11c-b §§1–2, inherited by pointer. c2 consumes two already-built, per-engine-reviewed models:
    57	
    58	- **S11c-b** (`scripts/S11c_b_exports.py`, committed `af560257`; step record `steps/S11c_b_variable_coefficient_operator.md`):
    59	  the variable-coefficient slab operator `slab_operator` over `{u,θ,e_W}` (write-key `slab_operator`, tag
    60	  `S11CB_SLAB_OPERATOR`), its term provenance `slab_operator_term_origins`, the off-diagonal coupling kernel
    61	  `coupling_kernel`, the constitutive operand `mu_theta_operator`, and — inherited from S11c-a — the θ-row flux
    62	  closure carrier `closure_shape_deriv`. ⚠ **S11c-b's full cross-engine residual is DEFERRED (≥64 GB); its coupling
    63	  content is per-engine leg-verified only, and two whole-row sign conventions + #90's two flags are
    64	  cross-engine-UNVALIDATED** (`steps/S11c_b_variable_coefficient_operator.md:112-115`): the **kinetic** convention
    65	  `−K` PY vs `+K` WL (a bulk term, independent of the response slots), the **face generalized-force** convention PY
    66	  `+diff` vs WL `−linearVirtualVariation`, and the **#90 closure-fold sign**. ⛔ The face-force and closure-fold
    67	  conventions are the coefficients of the slots c2 substitutes into (§3c) — they do **not** cancel from c2's residual.
    68	- **S11c-c1** (`scripts/S11c_c1_exports.py`, per-engine-reviewed 44-row delta; step record `steps/S11c_c1_curved_bulk_closure.md`):
    69	  the closed face response, write-keys **`s11c_c1_face_response`** / **`s11c_c1_face_response_coeffs`** (⛔ NOT the
    70	  bare S11b `face_response`/`face_response_coeffs`, which c1 imported as its uniform-limit regression operand);
    71	  the DtN operator `dtn_operator` (per anchoring/face: `s11cc1_dtn_operator_{lab_held,material_advected}_{plus,minus}`),
    72	  the flat symbol `dtn_flat_symbol`, the two-momentum kernel `dtn_kernel`; the per-case response resolvents
    73	  `s11cc1_response_resolvent_{α}_{s}_{ρ}_constant`; the momentum symbols `s11cc1_{k,q}_{in,out}put*`; the FT-of-profile
    74	  objects `s11cc1_w1_profile_hat_transfer`, `s11cc1_w1_profile_jet_hat_{1,2,3}`.
    75	
    76	### 1b · What is cross-engine AGREE vs UNDECIDED in the c1 import — SUPPLIED HONESTLY (rule 6/16)
    77	
    78	⭐ **AGREE (cross-engine, load-bearing; c1 reconcile `directives/_measurements/S11c_c1_comparator_reconcile.md`):**
    79	the **two-momentum DtN kernel** `dtn_kernel` (all four anchoring×face cases collapse to exact zero off-diagonal; a
    80	wrong jet sign or a one-leg freeze leaves it nonzero) and the **`δp_s` (pressure) + `J_s` (relative-flux) response
    81	coefficients** (all leaves collapse at physical kinematics; sweep `AGREE=54`). The flat symbol on the diagonal and
    82	the parity matrix (= the kernel) AGREE.
    83	
    84	⛔ **UNDECIDED — c2 must NOT treat these as cross-engine-closed (c1 step record §"Established vs owed"):** (i) the
    85	**background density** `rho_br_bg_rho4_constant` — a **surfaced rule-17 freeze**: WL keeps it a live applied field,
    86	PY froze it to a bare constant; (ii) the **`t_s` traction** — WL a zero-padded 4-vector `(0,0,0,scalar)`, PY a
    87	scalar; (iii) the **raw `dtn_operator` whole-form** — kernel-AGREE does **not** extend to the whole noncommutative
    88	object; (iv) the **off-diagonal flat-resolvent leg-labeling** — PY output-leg `q_out` vs WL input-leg, equal on the
    89	`k=k′` diagonal, differing off-diagonal in the MATERIAL_ADVECTED response coefficients; (v) the **ENERGY** audit (PY
    90	closed-form vs WL far-field integral). Each of (i)–(v) enters c2's fold load-bearingly (§3d); ⛔ folding any of them
    91	to force AGREE is the exact defect this rebuild exists to catch.
    92	
    93	### 1c · The S11c-b slab operator and the symbolic slots the fold closes — SUPPLIED (grounded in the real row)
    94	
    95	`slab_operator` is the divergence-form variable-coefficient operator over `{u,θ,e_W}`, per `(anchoring, density)`,
    96	with every background coefficient (`μ_R,bg`, `W_bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, the `Σ_E⁰` map) and its first jet
    97	retained live. **It is a two-face object** — the mass row is the sourced evolution `∂_tΣ + ∇_x·(Σv) = −(J₊+J₋)`
    98	(both faces summed) and the mechanical rows sum both per-face tractions. The rows are: `U_BODY_BALANCE`,
    99	`E_W_BALANCE` (the constraint-reduced mechanical rows, pin B), and the θ/mass row `evolution_mass_balance − Σ_s
   100	closure_shape_deriv_s` (the #90 fold). ⛔ **The explicit relative-flux `J_s` carrier has already been eliminated
   101	from the θ-row by the #90 subtraction** — what remains is the flux **closure** written out as
   102	`Λ_A𝒜_s + Λ_V V_s` with the face pressure `δp_s`/`d_w_delta_p_s` **symbolic**. In the real row (verified against
   103	`closure_shape_deriv`), the θ-row carries `−4I·Λ_A(−δp_s/ρ_m + μ_θ/ρ_br)/(ωτ_A+I)`, `−2I·Λ_V W_0 e_{W,t}/(ωτ_V+I)`,
   104	and the pressure-jet terms in `d_w_delta_p_s`. ⇒ **the fold operation is "substitute the closed `δp_s(V_s,μ_θ)` and
   105	its w-jets `d_w_delta_p_s` into the symbolic `delta_p_±`/`d_w_delta_p_±` slots"**, ⛔ **NOT "substitute a closed
   106	`J_s`" (there is no `J_s` slot; adding one would double-count the already-folded flux closure).**
   107	
   108	`mu_theta_operator` (`S11CB_MU_THETA_OPERATOR`) is the separate held-fixed constitutive operand (neither row).
   109	⛔ The face-force and #90 closure-fold sign conventions above are **not** normalized by c2; the §3c increment carries
   110	each engine's own convention, and the §3d.4 mechanical-power pairing adjudicates the face-force sign (rule 1/6).
   111	
   112	### 1d · The face closure laws and the Λ-channel routing — SUPPLIED (this routing is c2's task)
   113	
   114	Per anchoring `α` and face `s`, on the curved-face objects of the S11c-a substrate (⛔ not flat Cartesian objects):
   115	
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:31:independent outgoing-flux route. Every object is multigraded by `(ε,η,σ_W)` and carries restored physical
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:56:unmodified. Its exclusion is legitimate at first shape order (the drain-projection correction is `O(σ_W²)`; the
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:57:outward normal obeys `n̂_s·ŵ − s = O(σ_W²)`). ⚠ The stated validity domain is **non-uniform near grazing** and is
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:126:On the uniform (flat, `η=σ_W=0`) background S11b solved this bulk closure: the flat impermeable impedance `Z`,
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:128:supplied S11b objects, imported/re-derived **only** as the `(η,σ_W)→(0,0)` regression operand of §5c:
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:202:`h_s^α` (carrying `η`,`σ_W` via `W_bg`) is the curvature that distinguishes c1 from S11b's flat B0b.
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:226:`σ_W ≡ η W̄₀/L_W`, varied independently (by `η` and `L_W`); ⛔ no engine may replace `σ_W` by `η` or assign them a
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:228:bookkeeper `ε`. **Every computed object is multigraded by `(ε,η,σ_W)` from its actual data dependency**; the
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:230:`η` and `σ_W`** (`N5`/`N12`). ⛔ No term is removed merely because it carries both a wave and a background
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:232:background amplitude and is not dropped; `|∇_x h_s|²` is `O(σ_W²)` and is not a first-shape-order ingredient.
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:357:carries the object's computed `(ε,η,σ_W)` multigrade and restored `[L,T,M]` dimension, and states no component,
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:368:- **The representation-invariance package** (§5a) ⇒ `S11CC1_REP_INVARIANCE_EULERIAN_OPERAND`,
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:457:For the curved-face DtN operator and the permeable face response, independently obtain the `(η,σ_W)→(0,0)`
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:462:coefficient/sign/parity/branch structure, which vanish at `(η,σ_W)=0`. `(η,σ_W)→0` is **not** an accepted
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:465:### 5d · Zero-jet regression (`σ_W→0` at finite `η`)
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:467:Operand A: c1's DtN in the limit `σ_W→0` with `η` **retained** and `w₁` a **constant symbol** (`∇W_bg=0`, so the
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:475:false zero. The genuine test is whether c1's `σ_W→0` operand acquires a spurious `η`/thickness dependence — the
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:515:Restore units and compute the `[L,T,M]` dimension of every emitted object. `η`, `σ_W`, `w₁`, `m₁`, `θ`, `e_W`
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:553:reserved by S11c-a/S11c-b are inherited unchanged (`W_bg`, `w1_profile`, `L_W`, `sigma_W`, `mu_R_bg`,
research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:620:residual for the representation-invariance (Eulerian vs Hanzawa second route), the mandatory `W_bg`-tilt one-sided
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:178:σ_W       ≡ η W̄₀/L_W .
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:185:∂_{yᵢ}W_bg = σ_W ∂_{ξᵢ}w₁ ,
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:186:∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W ∂_{ξᵢ}m₁ .
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:189:Thus the zero-jet contrast bookkeeping uses `η`, while the supplied first-jet bookkeeping uses `σ_W`.
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:190:`η` and `σ_W` are varied independently by varying `η` and `L_W`; no engine may replace `σ_W` by `η` or
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:195:independent amplitude bookkeeper `ε`. Every computed object must be multigraded by `(ε,η,σ_W)` from its
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:408:computed `(ε,η,σ_W)` multigrade and restored physical dimension. Do not state a component, term, order,
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:506:For each S11c-a object, independently obtain the `(η,σ_W)→(0,0)` operand and the corresponding S11b object.
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:521:Restore units and compute the `[L,T,M]` dimension of every emitted object. `η`, `σ_W`, `w₁`, and `m₁` are
research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md:540:`L_W`, `sigma_W`, `mu_R_bg`, `m1_profile`, `rho_4D_bg_rho4_constant`,
research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md:165:  justified independently of the outcome (FT-of-derivative jet identity; σ_W=η_bg·W_0/L_W binding; on-shell
research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md:192:  = η², η·σ_W, σ_W²** (η and σ_W are INDEPENDENT grades, spec §2c), ⛔ not only `O(η²)`. Physical conclusion (no
research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md:194:- **Drain-projection `O(σ_W²)` wording.** The `O(σ_W²)` applies to the geometric **drain-tilt** projection
research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md:195:  (`n̂·ŵ−s=O(σ_W²)`), ⛔ NOT to convection — convection is dropped as an inheritance of N11a's standing rest-frame
research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md:196:  limit (`d(rel)/dη=0`, has no η), ⛔ not because "convection is O(σ_W²)". Wording fix. (Grok NIT)
/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md:5:⭐ CARRY-FORWARDS owed to the c2 STEP RECORD (⛔ NOT build defects — the fold is SOUND): (F) the uniform-limit object does NOT literally vanish because the §3c increment retains the −extract(open) open slot — the record must say "the genuine closure-induced coupling **DECOUPLES**" (Integral integrand ≡0, verified), ⛔ NOT "the increment vanishes"; (E) N6 rep-invariance HOLDS at leading O(ε)/O(εη) (σ_W→0 kills the whole residual), the σ_W-sector remnant is DEFERRED (drain-projection ≥64 GB family); (G) the induced self-energy is **DIRECTIONAL/one-way** (reverse THICKNESS→TRANSVERSE block ≡0 — δp lives only in the θ/mech rows; a dissipative directional self-energy, expected + both blocks emitted); the two S11c-b sign conventions (face-generalized-force, #90 closure-fold) multiply the increment and are SURFACED by the WL cross-engine comparator (⛔ not normalized; §3d.4 mechanical-power pairing adjudicates the face-force sign).
research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md:71:Command: same `…_verify_EG.py`/`.out` (block directionality) + builder-report line 70.
research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md:99:3. **G:** the induced self-energy is **directional** (transverse→{θ,e_W}; reverse block identically zero) —
STATUS.md:6:**(A) c2 PHYSICS 2-leg review DONE — fold physics SOUND, 0 confirmed defects, committed `8f3a017f`.** Real independent legs (Codex-written script → **fresh Claude agent + Grok**, identical prompt, blindness-by-absence). Both agree the fold WIRING + A/C/D1–D6 are SOUND (δp not J_s, operator-inverse, kernel bridge, V_s→face_velocity, computed jets, ε-strip; FORM ablations respond). **Grok flagged F/G/E as CONFIRMED defects — ALL FALSE POSITIVES**, refuted by my own computation (rule 13; `_measurements/S11c_c2_physics_review_adjudication.md` + `…_adjudication_verify_{F,EG}.py/.out`): **F** uniform-limit — the genuine closure-induced coupling **DECOUPLES** (Integral integrand ≡0, `.doit()`→0 all 4 cases; Grok saw `Integral(...)` without evaluating it — its recurring pattern, cf. c1 "tangential freeze"); **G** adjointness — the builder never claimed adjoint-by-construction, both blocks emitted, the self-energy is **DIRECTIONAL/one-way** (reverse block ≡0, expected since δp lives only in the θ/mech rows) — dissipative directional, not a defect; **E** N6 — leading O(ε)/O(εη) rep-invariance **HOLDS** (σ_W→0 annihilates the whole residual), σ_W-sector remnant DEFERRED. **Carry-forwards (step-record, ⛔ NOT defects):** F wording "genuine coupling decouples" NOT "increment vanishes" (→ light §5e/§3c spec fix owed); E σ_W-sector deferred; G directional self-energy; the 2 S11c-b sign conventions → the WL comparator.
STATUS.md:29:- **The design (committed):** **D1 Membership = the BIND-CLOSURE** — export a row iff a later step binds its **F9 write-key** as a construction/regression operand, recursively closed over symbols + `dimension_key`; **by binding, not label**; missing manifest = **hard error**. **D2 Topology = GENERATE over a frozen base** — F9 is write-time, the fold is chronological **last-wins (never re-applies F9)**, **no write-time pass-through**; the existing files **through `S11c_b_exports.py` are FROZEN** ⇒ no chain re-run, no `BUILD_INPUT_DIGESTS` cascade. **D3 Guard** = manifest + **edge resolution by full SYMBOL IDENTITY (srepr = name+assumptions)** + bidirectional smoke-test; **`D3`-round-trip provably does NOT catch under-export**. **D4** = representation deferred; **`coupling_kernel` DROPPED** from the minimum (c2 re-extracts, doesn't bind it) ⇒ minimum ~319 rows/18 MiB.
STATUS.md:139:## ⭐⭐⭐ CURRENT FRONT — S11c-b (variable-coefficient brane operator + off-diagonal coupling kernel; S11c's FIRST physics number). **ADJUDICATION LAYER COMPLETE (3 phases reviewed-sound) + STEP 1 DONE; STEP 2 IN PROGRESS — only ADMISSIBILITY-θ is a SETTLED verdict, advective/kinetic/coupling are LIKELY/PENDING (2026-08-29). THREE commits: `bb0bfc02` comparator time-order fix (canon_jet_name recorded time differentiation as a Boolean `has_time` → WL `D[,{time,2}]`→`u_t` collapsed while PY `u_tt` stayed = asymmetric, flattening transverse-trial ∂²ₜ in COUPLING; count time tokens, emit `_tt`; both build legs SOUND, S11c-a NOT regressed [28/28 diffs canonically zero]; ⭐ discovered BY this adjudication's round-3 decision-list legs); `2e5c6755` adjudication layer v1 `scripts/S11c_b_adjudicated_comparison.py` (Bridge A `bRho↦B_rho_3/W_0` + sort-routing + total-bijection containers + exact-multiset accounting + jet-conservation; both build legs SOUND, 38 MATCH all genuine NO false MATCH via triple ablation, ⭐ Bridge A creates ZERO match — the 38 agreements are RENAME-LEVEL); `afc276b7` adjudication layer v2 (Bridge D = the engine's own committed `PROFILE_GRADE_SUBS` background expansion + IBP/total-in-plane-divergence classifier for weak-pairing DENSITIES ONLY + strong-operator EXACT + atom-gated protection; the v2 directive took 2 decision-list leg rounds catching a naive-chain-rule Bridge D [engines use `sigma_W`/`L_W` not `W_0·eta_bg`], an over-broad divergence classifier [valid only for weak densities not strong operators], and Bridge-D-∂ non-commutation; both build legs SOUND, one reimplementing the Euler/divergence test FROM SCRATCH). ⭐ MEASURED CROSS-ENGINE RESULT (rule 6 — a disagreement IS the measurement): route counts `MATCH=38, PROTECTED_UNREDUCED=32, FLAG=12, RESIDUAL_BULK=8, STRUCTURE_INCOMPLETE=57, COVERAGE=84` (231, multiset-equal, JET 231/0). The adjudication routed 38 rename-level MATCH (ADMISSIBILITY operator(16)/support(20)/COUNT(2)); 32 REPRESENTATIONAL [protected ENERGY_BASIS quotient reps 07/10 + gamma-DivGrad — non-unique by design]; and 20 INDEPENDENTLY-NONZERO operator/coupling differences (not rename/bridge-level — ⚠ step 2 refines which are genuine-physics vs representational): ADMISSIBILITY-THETA `∇²w1` (PY carries, WL 0), SLAB_TERM_ORIGINS ADVECTIVE, KINETIC `W_bg²−W_0²`, 8 COUPLING_KERNEL bulk certified non-IBP by the v2 legs. `REPRESENTATIONAL_DIVERGENCE=0` is REAL (both legs; a from-scratch Euler operator with exact rational zero-test certified all 8 coupling residuals as genuine bulk; the classifier recovers a verified `V` for known divergences up to 5218 ops). ⚠ the original "systematic higher-background-order WL truncation" hypothesis was ⛔ FALSIFIED by STEP 1. **⭐ STEP 1 DONE + COMMITTED `e5bf4122` (2026-08-29): the `(eta_bg,sigma_W)` MULTIGRADE INSTRUMENT `scripts/S11c_b_background_multigrade.py` (both build legs SOUND — fresh Claude agent + Grok each re-derived the coefficients by a DIFFERENT route with ZERO mismatch, form-ablated every guard live; one tautological guard RECONSTRUCTION folded out; extractor `scripts/S11c_b_grade_fingerprint.py`; run `~/.s11_build/S11c_b_multigrade_run.out`; records `_measurements/S11c_b_multigrade_instrument_{build_directive,directive_review,build_review}.md`). ⭐ COMPUTED FINGERPRINT (`(a,b)=(eta_bg^a,sigma_W^b)`; A=PY, B=blind WL) — the 20 differences are a per-family MIX, NOT uniform WL under-retention: ADMISSIBILITY-THETA A={(0,1),(1,1)} B=∅ (WL lacks ∇²w1 entirely); KINETIC THICKNESS_ROW A={(0,0),(1,0),(2,0)} B={(0,0)} (u-momentum rows AGREE); ADVECTIVE A={(0,1)} B={(0,0),(0,1),(1,0)} ⇒ WL carries ∇·u_t that PY LACKS (OPPOSITE direction); COUPLING BIDIRECTIONAL A={(0,1)..(4,1)}+rem(5,1) B={(0,1),(0,2),(1,1),(1,2)}. ⭐ STEP 2 ADJUDICATION (WIP, recorded `_measurements/S11c_b_step2_adjudication.md`; user chose FULL RIGOR all families 2026-08-29): (1) ADMISSIBILITY-θ = GENUINE WL under-retention [SETTLED] — §3a mandates retaining ∇²W_bg (2nd spatial deriv still first shape order), WL operand is 0. Mechanism (leading hypothesis, CONFIRM at repair by running WL): WL's `truncateBackground` Series-to-first-background-order truncation drops it (vs a missing energy term) while PY keeps it at σ_W¹. CLEAR repair candidate. (2) SLAB ADVECTIVE = LIKELY representational provenance re-bucketing (PY = ∇W_bg·u_t, WL = full ∇·(W_bg u_t); PY assigns W_bg∇·u_t to ACCUMULATION/FLUX) — CONFIRM via the deferred ENERGY_BASIS quotient reconciliation. (3) SLAB KINETIC = LIKELY genuine WL under-retention (W_bg²→W_0²; direction PY-has-more). (4) COUPLING = certified non-IBP bulk (v2 legs) with a MEASURED bidirectional grade pattern; the §2a per-cell "which engine is spec-correct" verdict is PENDING. ⛔ NEXT (full rigor): (A) SLAB advective/kinetic — the ENERGY_BASIS quotient reconciliation via VARIATIONAL-EQUIVALENCE (⛔⛔ NOT `classify_total_divergence` — v2 legs proved it invalid for strong operators); (B) COUPLING §2a/§3a per-cell adjudication; (C) then repair ONLY the confirmed genuine engine gaps (⇒ admissibility ∇²w1, likely kinetic; NOT advective if representational) → re-run WL engine + comparator + adjudication → re-adjudicate survivors (rule 13) → honest step record + S11c card + close. ⛔ COMPUTE don't assert (rule 4); ⛔⛔ NEVER blanket-collapse (the four families adjudicate DIFFERENTLY); ⛔ a one-engine fix is a SPEC question FIRST (if §3d/§2a ambiguous for a cell → fix spec). OWED (carry to card): v2 N1/N2/N3; 57 STRUCTURE_INCOMPLETE + 12 control NAMESPACE cross-engine-owed; S11c-a control-family keying; admissibility §5 control-coverage.** **[SUPERSEDED — historical; the reconcile phase below LED HERE. Read as prior state, not the result.] RECONCILE PHASE (2026-08-28): engines review-clean + transcripts committed; T7 comparator BUILT + BOTH BUILD LEGS SOUND (fresh Claude + Grok FORM-ablated every fold) + committed `5f01f9fa`; my run of the COMMITTED comparator is byte-identical to both legs (224727199 bytes; families=28 all-join, 7 unpaired, 0 verdict tokens). One comparator finding (Grok: FACE_FLUX "silent drop") verified + downgraded — visible py_only via COUPLING_KERNEL_TERM_ORIGINS, not integrity. Hand-coded reconcile `S11c_b_handcoded_comparison.py` BUILT + BOTH build legs SOUND (map COMPLETE, no false MATCH, `--drop-rename` surgical/load-bearing) + committed `82ec3b5f`; run (`~/.s11_build/S11c_b_reconcile_run.out`): **2 MATCH** (ADMISSIBILITY_SUPPORT 20/20, ENERGY_BASIS_COUNT bare-int 26), **14 FLAG** (all core objects), **12 NAMESPACE_INCOMPLETE** (control families, WL `<|…|>` unparsed → OWED). DIMENSIONS: the digest shows the dimension VALUES agree (a container-shape diff), but its reconcile verdict is still `FLAG` — confirm in adjudication, do not pre-call it. Both legs verified the 14 FLAGs are NOT naming noise (map complete); representational-vs-finding is the pending adjudication, not yet computed. **⛔ NEXT + FINAL = the ADJUDICATION (decides the number; the never-blanket-collapse danger zone):** apply the further REVIEWED bridges + re-check zero — (i) bRho→`B_rho*W_0` scale, (ii) gamma-index selected-representatives, (iii) WL `f(xOne,xTwo,xThree,time)`→PY bare GATED on a computed 0-spatial-jet check per field (S11c-a INERT_APPLIED gate), (iv) integral linearity; ⛔⛔ NEVER collapse `W_bg`/`mu_R_bg` (real background-gradient jets = the variable-coeff physics; blanket applied→bare hid the S11c-a in-plane current freeze). Reduces-to-0 ⇒ engines AGREE on that object (coupling-kernel number); survives ⇒ genuine finding (rule 1/6). ENERGY_BASIS = non-unique quotient, never fold a representative. Should be a REVIEWED artifact (reconcile-v2 with the gated bridges → 2 legs → commit) — COMPUTE, don't assert (rule 4). Likely mirrors S11c-a: agreement up to reviewed representational identities.** Journey: spec drafted → 2 legs (8 defects) `7023420d` → both blind engines built (Codex) → 4 build legs found a CONFIRMED rule-7 shared-blind-spot (both engines independently produced a vacuous ε→0 admissibility from an ambiguous §3d — a comparator alone would have read it as agreement) + a 4-leg tautological Clairaut adjoint + per-engine kernel-extraction bugs → my spec round-2 fix was STILL insufficient (2 re-legs) → **rule-15 author change: Codex rewrote §3a/§3c/§3d/§5a, 2 re-legs CLEAN** `0c0e9a8a` → engine repairs: WL W1-W3 clean both legs `49d9fad1`; SymPy B2/B3/B4 clean, B1 took a round-2 (the round-1 legs DISAGREED — Grok caught a scalar over-promotion the Claude-agent leg missed; the sibling WL engine's independent-correct §3d construction was the tiebreaker → SymPy-only fix) → **SymPy B1 round-2 clean both legs `103cdfdb`**. Both engines now review-clean. **OWED (non-blocking, carry to the card):** the admissibility operand is verified correct but has no EMITTED §5 discriminating control (independence omits on structural absence, rep-invariance is a background-order structural zero, §5b/§5c don't name it) — a control-coverage refinement, physics verified; plus S11c-a's owed control-family keying. **[SUPERSEDED — the live NEXT is the gated-bridge adjudication in the RECONCILE-PHASE front clause above; comparator + reconcile are built/legged/committed/run. The following is the historical plan record.]** Transcripts committed (SymPy `55abd09b` 167MB, WL `73b4e639` 97MB); T7 comparator (Codex `bqj3ftjqc`, brief `directives/S11c_b_comparator_build_directive.md`) was verified + 2 build legs + run + reconciled per the front. Finding records `_measurements/S11c_b_*`. [Historical spec-build detail below.] **SPEC DONE + 2 legs folded + committed `7023420d` (2026-08-27).** The SHARED PHYSICS spec `directives/S11c_b_SHARED_PHYSICS.md` was authored (mirrors the S11c-a template; inherits S11c-a §§1–3 by reference — SymPy imports `S11c_a_exports`, blind Wolfram re-derives). Two decision-list legs (orchestrator-written → Codex + Grok) ran BEFORE any engine (rule 7) and were **very productive**: they converged on 2 serious shared-blind-spots and surfaced 6 more — all 8 verified against source (rule 13) and folded ONCE (not iterated to green). The 8: **D1** ζ_s=ζ_c+sδW/2 not ζ_±=±δW/2 (was forcing ζ_c=0); **D2 (serious)** confine "u only through gradients" to the uniform quotient — the non-uniform background admits undifferentiated-u spurion couplings, the very N15 channels the step emits (my gloss would have made both engines agree the coupling is absent); **D3** declare which quantities vary (W_bg,μ_R,bg,density) vs constant moduli (N12); **D4** enumerate ONE thickness coordinate + impose e_W,bg map before the rank test; **D5 (deep)** the total-divergence quotient does NOT lift to variable coefficients (`c∇·F≡−(∇c)·F` ⇒ first-jet terms are PHYSICS) — dropped the leaky μ_⊥ citation + removed the pre-registered representative fold from the comparator (could have masked the coupling); **D6 (serious)** admissibility = background-order (ε⁰) balance in the 𝒮_hold⁰ pairing, not ε→0 of the wave operator (vacuously 0) — new §3d; **D7** sectors via local curl/div structure not a global projector (N5), adjoint w.r.t. the supplied pairing; **D8** separate one-source form ablations for ∂W_bg and ∂μ_R,bg (independent profiles). Finding record `directives/_measurements/S11c_b_spec_review.md`; leg prompt `_legs/S11c_b_spec_review.md`; raw transcripts `~/.s11_build/S11c_b_spec_review_{grok,codex}.txt`. Rule-5 leak scan clean (no coupling grade/sign, basis count, or μ_⊥ identity stated). **[historical, DONE — the plan at spec-close was]** the two blind engine builds against the spec (SymPy chained on `S11c_a_exports`; Wolfram blind), each with 2 build legs, then this sub-step's own T7 comparator + step record. OWED (carry in from S11c-a, not blocking): bridge the control-family keying (REP_INVARIANCE/CONTROL_INDEPENDENCE PY-missing-DENSITY; CONTROL_FORM FACE-granularity) + a CONTROL_FORM re-characterization.
STATUS.md:599:bidirectional set equivalence; **full fresh D3 cell 391/391 byte-identical (1117 s)** — a scoped cold
STATUS.md:655:bidirectional equality + a non-vacuous driven witness — then 2 script legs clean, tag diff exactly the

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md | sed -n '45,65p;216,238p;600,630p' && nl -ba research/pde_ledger_v3/steps/S11c_c1_curved_bulk_closure.md | sed -n '145,210p' && rg -n 'REP_INVARIANCE|sigma_W|defer|64 GB|drain' research/pde_ledger_v3/_measurements/S11c_c1_comparator_reconcile.md research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md 2>/dev/null | sed -n '1,220p'" in /var/projects/toy_physics
 succeeded in 0ms:
    45	  `μ_θ` feeding `μ_s`, and re-derives nothing of S11c-b;
    46	- S11c-d: any profile-conditioned scattering, Born, Bloch, WKB, resonance, or spectral object; ⛔ **no global
    47	  dispersion relation** `ω(k)` for a generic `W_bg(x)` is requested (`N5`). c1 emits the DtN and response as
    48	  **operators/kernels**, ⛔ not their spectrum, and the **profile-conditioned resolvent / singular set** of the
    49	  operator loci is S11c-d's, not c1's (§3b);
    50	- S11c-e: the flux-normalized dimensionless conversion form, the leakage observable, the confinement
    51	  interpretation (`N13`), and the bench-optics comparison (`N7`);
    52	- all nonlinear DC, harmonic, sideband, intensity, and soliton questions (`N10`/`G14`).
    53	
    54	⛔ **The convective bulk operator is NOT an S11c task (`N11a`, settled in the decision list).** `v_bulk_normal_0`
    55	is inherited as a standing rest-frame limit; the rest-frame bulk wave equation `∂_t²φ=c_s0²∇₄²φ` (§1b) stands
    56	unmodified. Its exclusion is legitimate at first shape order (the drain-projection correction is `O(σ_W²)`; the
    57	outward normal obeys `n̂_s·ŵ − s = O(σ_W²)`). ⚠ The stated validity domain is **non-uniform near grazing** and is
    58	made precise in §2b; ⛔ `v_bulk_normal_0` is never aliased to `v_0`.
    59	
    60	There is no acceptance value to withhold in S11c-c1. There is no terminal `VERDICT`, `PASS`, or `FAIL`.
    61	
    62	---
    63	
    64	## 1 · Complete inherited setup — SUPPLIED and unfalsifiable in this build
    65	
   216	- away from grazing, each result is conditional on the S11b-derived domain, made explicit as the pair
   217	  `|q_out·v_bulk_normal_0/ω| ≪ 1` **and** `|ω·v_bulk_normal_0|/(c_s0²|q_out|) ≪ 1` (the boundary-layer exclusion),
   218	  together with the independent **subsonic** speed condition inherited from S11b (`v_bulk_normal_0 < c_s0`); the
   219	  order of the limits (`v_bulk_normal_0→0` vs `q_out→0`) is stated explicitly;
   220	- ⛔ the engines still compute only the rest-frame object; this domain is recorded by the step record, not carried
   221	  as a term in any operator.
   222	
   223	### 2c · Multigrade bookkeeping — SUPPLIED
   224	
   225	The background bookkeepers are S11c-a's: the zero-jet contrast `η` and the supplied first-jet bookkeeper
   226	`σ_W ≡ η W̄₀/L_W`, varied independently (by `η` and `L_W`); ⛔ no engine may replace `σ_W` by `η` or assign them a
   227	common order. Wave and bulk perturbations (`u`, `ζ_c`, `δW`, `θ`, `φ`, `δp`, `V_s`, `J_s`) carry the amplitude
   228	bookkeeper `ε`. **Every computed object is multigraded by `(ε,η,σ_W)` from its actual data dependency**; the
   229	requested truncation is **first order in wave amplitude `ε` and first shape order in each background bookkeeper
   230	`η` and `σ_W`** (`N5`/`N12`). ⛔ No term is removed merely because it carries both a wave and a background
   231	bookkeeper; ⛔ no object's grade is stated here. A second spatial derivative of `W_bg` remains first order in
   232	background amplitude and is not dropped; `|∇_x h_s|²` is `O(σ_W²)` and is not a first-shape-order ingredient.
   233	
   234	---
   235	
   236	## 3 · The bulk-closure construction
   237	
   238	### 3a · The curved-face impedance / DtN — a two-momentum nonlocal OPERATOR
   600	model via the fold over the atomic frozen base `scripts/S11c_b_exports.py` (the whole F9-resolved S11b + S11c-a +
   601	S11c-b model), binding only its declared `IMPORT_KEYS`; the Wolfram engine re-derives the supplied §§1–2 inputs
   602	and the S11c-a substrate it consumes and the S11c-b `μ_θ` operand, importing nothing.
   603	
   604	---
   605	
   606	## 8 · Supplied versus computed; builder report
   607	
   608	**Supplied and unfalsifiable here:** all of §1 (the inherited DOF and sector split; the bulk acoustics, radiation
   609	condition, and branch object; the flat B0b/B0c reference used only as the §5c regression operand; the face
   610	closure, affinity, kernels, traction with its `Λ_X` channel, and kinematic balance); the background ansatz,
   611	anchoring, and curved half-space setup of §2; the standing rest-frame limit and its non-uniform validity domain
   612	(§2b); and the S11c-a shape-derivative substrate (T-a..T-i) and the S11c-b `μ_θ` operand — imported by SymPy and
   613	re-derived by Wolfram, per-engine-verified with the S11c-b cross-engine residual deferred.
   614	
   615	**Computed here:** the curved-face DtN operator (the composition form and the two-momentum `k,k′` kernel, the
   616	flat symbol, the rigid-shift residual, the regime-pair/parity/branch structure, the Hermitian and reactive parts)
   617	(§3a); the permeable curved face response, its formal noninvertibility condition and finite-dimensional loci, and
   618	its dissipation audit as three distinct objects — bulk-radiation Hermitian part, the two-port permeable-port
   619	Hermitian form, and the independent traction-vs-far-field-flux energy balance (§3b); and both operands and the
   620	residual for the representation-invariance (Eulerian vs Hanzawa second route), the mandatory `W_bg`-tilt one-sided
   621	independence control, the `W_bg` per-direction form-ablation, the uniform-limit, the zero-jet, the branch/
   622	momentum-liveness (both legs), and the homogeneity packages (§§5–6) — each with its computed multigrade and
   623	dimension. (The `Σ_E`-advection and `μ_R,bg` channels are structurally absent from c1's `μ_θ`-symbolic objects
   624	and are reserved for S11c-c2.)
   625	
   626	The builder's report is under 25 lines and gives the files written, line and tag counts, tasks run or skipped,
   627	runtime, all emitted tag names, and any ambiguity or non-computable requested object. It states that §§1–2 and
   628	the S11c-a/S11c-b substrate were supplied and unfalsifiable in this build, that the rest-frame limit and its
   629	non-uniform validity domain condition every result (§2b, grazing = strict `v_bulk_normal_0=0`), and reports no
   630	computed value and no conclusion about the physics.
   145	  [[feedback_handcode_comparison_never_blanket_collapse]]
   146	- ⚠ **Exact grazing — claim SCOPED: `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`** (retrospective spec review). The
   147	  **DtN inverse `N⁻¹`** (hence `Z=iρ_mω·N⁻¹`) is **nonanalytic as `q_out→0`** (the grazing limit; at exact double
   148	  grazing `N=η·B` so `N⁻¹` — and thus `Z` — carries a `~1/η` Laurent pole, ⛔ not a Taylor series in η), so the
   149	  first-shape coefficient `Z₁` is a valid **non-grazing asymptotic** coefficient ONLY — `‖N₀⁻¹N₁‖≪1` imposed on both
   150	  legs (both momentum legs bounded away from grazing). ⚠ It is `N⁻¹`/`Z` that is singular, ⛔ **NOT the permeable
   151	  face-closure resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹`**, which for generic `Λ_A≠0` is **regular** at grazing (`O(η)`, and
   152	  `[I+(Λ_A/ρ_m²)Z]⁻¹·Z → ρ_m²/Λ_A` finite — both retro legs printed "D⁻¹ 1/η-pole = False, Z=C·N⁻¹ 1/η-pole =
   153	  True"). The GRAZING regime label is carried as a keying axis (reconcile §4 `DTN_GRAZING_BEHAVIOUR`) and the
   154	  strict-`v_dr=0` rest-frame qualification is correct, but the exact-grazing *threshold response* is ⛔ **not
   155	  claimed** at c1's first shape order — it needs a degenerate threshold solve (deferred).
   156	  [[feedback_never_freeze_a_varying_field]]
   157	- ⚠ **Independence is SCOPED, ⛔ not blanket** (retrospective spec review). The c1 spec supplied the composition
   158	  recipe and some expected structural values — rigid-shift cancellation, the flat `Z₀=ρ_m ω/q_out`, the zero-jet
   159	  outcome (a rule-5/rule-3 leak). So for THOSE objects part of the cross-engine "agreement" is **fidelity to the
   160	  supplied structure**, ⛔ not independent discovery. What IS independently confirmed is the **two-momentum DtN
   161	  kernel**: both blind engines AND both retrospective spec-review legs re-derived it from first principles (Grok
   162	  Route A/B, Codex boundary solve — residual 0). Read the "agreement is independent construction" claim (in *The two
   163	  blind engines* above) with that scope: the kernel earns it; the supplied objects earn it only as fidelity.
   164	- ⛔ **Never blanket-collapse the cross-engine bridge.** Each identification was read from the two engine SOURCES and
   165	  justified independently of the outcome (FT-of-derivative jet identity; σ_W=η_bg·W_0/L_W binding; on-shell
   166	  dispersion; ε-placement) — ⛔ never tuned to force a zero. [[feedback_handcode_comparison_never_blanket_collapse]]
   167	- ⛔ **Serialize CAS jobs; watch RSS.** c1 measured LIGHT (comparator peak ~317 MB on 30 GB); the 4 giants + full
   168	  residual are the ≥64 GB work. Detached launch (harness reaps `run_in_background` ~87 s).
   169	
   170	## Carry-forward corrections (lower-severity — from the retrospective spec review, 2026-09-05)
   171	These change what may be CLAIMED, ⛔ not any computed object (c1's engines/exports STAND — no reopen). Full list +
   172	both leg reports: `_measurements/S11c_c1_spec_retro_review_adjudication.md`,
   173	`_measurements/S11c_c1_spec_retro_review_{grok,codex_sol}.txt`.
   174	- **Energy-residual orientation.** The independent energy identity is `P_face + P_∞ = 0` (traction work + **positive
   175	  outgoing** far-field Poynting = 0, with `P_∞` the outgoing flux), ⛔ NOT the literal `A−B` the spec wrote — a
   176	  literal `A−B` on the written operands is `−2·δp·V*` for the supplied `t_s` and vanishes only for the WRONG `t`
   177	  sign. PY correctly compared to `−outgoing_flux` (`.py:1407,1459`), so no recomputation is owed; the fix is the
   178	  spec/tag SEMANTICS (to write it as a **vanishing `A−B` residual**, define the bulk subtraction operand `B` as
   179	  **minus** outgoing Poynting — `A−B = P_face−(−P_∞) = 0` — or emit `A+B`). Carry to c2's energy control. (both legs)
   180	- **`h_s` graph-vs-outward + `N`/`Z` terminology.** `h_s=s·W_bg/2+ζ_s` is the signed lab GRAPH height (face-odd); the
   181	  shape displacement that enters the kernel is the face-EVEN outward displacement `a_s=(W_bg−W₀)/2`. Both engines used
   182	  the correct outward height (PY `shape_source:592`, WL `:177`), so the kernel AGREES — notation only. Separately: `Z`
   183	  is the impedance (Neumann-to-Dirichlet response), the mathematical DtN is `N`, related by `Z=iρ_mω·N⁻¹` (⛔ not
   184	  synonyms); `N₁=−(N₀M_aN₀+Div(a∇)+κ²M_a)` vs `Z₁=iρ_mω·(−N₀⁻¹N₁N₀⁻¹)`. Terminology fix; no recomputation. (both legs)
   185	- **Density as a multiplication operator.** Name the live `1/ρ_br,bg` as a multiplication operator `β=M_{1/ρ_br,bg}`
   186	  so the O(εη) channel (Method notes › rule 17) cannot be emitted as a bare constant. ⇒ **the c2 build directive must
   187	  require this naming** (it makes the freeze un-emittable). Carry to c2. (both legs)
   188	- **`K_a` is Hermitian, ⛔ not anti-Hermitian.** The reactive object `K_a=(Z−Z†ₐ)/(2i)` is **Hermitian**
   189	  (`K_a†−K_a=0`); the anti-Hermitian part is `(Z−Z†ₐ)/2`. The spec's "anti-Hermitian" label for `K_a` is wrong — fix
   190	  the label wherever the dissipation objects are reused (c2 / S11c-e). (Codex NIT)
   191	- **Evanescent caveat covers all second-shape grades.** The evanescent-nullspace completion is **second-shape order
   192	  = η², η·σ_W, σ_W²** (η and σ_W are INDEPENDENT grades, spec §2c), ⛔ not only `O(η²)`. Physical conclusion (no
   193	  first-order passivity violation) unchanged. (Codex NIT; Grok concurs the caveat is physically right)
   194	- **Drain-projection `O(σ_W²)` wording.** The `O(σ_W²)` applies to the geometric **drain-tilt** projection
   195	  (`n̂·ŵ−s=O(σ_W²)`), ⛔ NOT to convection — convection is dropped as an inheritance of N11a's standing rest-frame
   196	  limit (`d(rel)/dη=0`, has no η), ⛔ not because "convection is O(σ_W²)". Wording fix. (Grok NIT)
   197	- **Flat `Z₀` / rigid-shift expected-value leakage.** `§5d` typed `Z₀=ρ_m ω/q_out` and `§3a` typed the rigid-shift
   198	  "must cancel" residual — rule-5 leaks in the (builder-facing) c1 spec. Physics correct; folds into the *Independence
   199	  is SCOPED* note above (agreement on `Z₀`/rigid-shift/zero-jet is partly fidelity-to-supplied-values). Spec hygiene;
   200	  the kernel remains independently confirmed. (both legs)
   201	
   202	## What's next
   203	**c2** (self-energy fold): fold the closed response into `S11CB_SLAB_OPERATOR`, re-extract the closed off-diagonal
   204	kernel from the CLOSED full operator → the coupled nonlocal self-energy operator. Held folds: extract-then-close
   205	ordering (close FIRST then re-extract — extract/eliminate don't commute), θ-row `Λ_X`(traction)/`J_s`(mass)
   206	routing, substitution-increment emit (c2's self-energy = closed − open-symbolic, so S11c-b's deferred cross-engine
   207	residual cancels out of c2's residual). c2 = a CODE build (gpt-6-astra) with its own gated directive (2 decision
   208	legs BEFORE the build). ⚠ NO per-substep card — N1 specifies ONE S11c roll-up card, produced only after S11c-e.
   209	The ≥64 GB c1 cross-engine re-run (above) is the outstanding c1 work; when it lands it confirms the
   210	inheritance-argument families symbolically and closes the UNDECIDED items (DTN_OPERATOR, ENERGY, flat-leg, and —
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:20:per-family symbolic residual are **UNMEASURED — DEFERRED** (≥64 GB, `DEFERRED_HEAVY_RUNS.md`), ⛔ not
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:22:whole-form, t_s, ENERGY, seal 5 and the giants are UNDECIDED/deferred.** Reproduce:
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:41:('sigmaW','sigma_W'),('tauA','tau_A'),('tauV','tau_V')]` — 12 folds incl. `cS0←c_s0`, injective.
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:121:δp_s, J_s AGREE; t_s is a surfaced scalar-vs-vector representation, UNDECIDED — deferred to the symbolic residual.**
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:163:`background_density_map`), or the deferred symbolic residual. [[feedback_never_freeze_a_varying_field]]
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:173:  form is not leaf-extractable; full adjudication DEFERRED to the ≥64 GB symbolic run.
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:196:- **REP_INVARIANCE_{EULERIAN,HANZAWA}_OPERAND / _RESIDUAL (j=50 each), CONTROL_* (j=26–156), DEGENERATE_LOCI_*,
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:199:  2 exact-zero joins in REP_INVARIANCE_RESIDUAL, 1 in DEGENERATE_LOCI_REAL_ADMISSIBLE are direct evidence. ⚠ This
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:208:## 5. Deferred to the ≥64 GB box (⛔ do NOT narrow to fit — recorded in `DEFERRED_HEAVY_RUNS.md`)
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:215:S11c-b's cross-engine-residual deferral. Mid-tier (FACE_RESPONSE + 3) result folded in §5a.
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:239:  UNIFORM_LIMIT_* families in the ≥64 GB deferral.
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:244:kernel), degenerate loci, dimensions. **NOT earned — UNDECIDED/deferred (⛔ not pre-adjudicated AGREE):** seal 5
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_reconcile.md:247:cross-engine closed; the operator whole-form, t_s, ENERGY, density and the giants await the ≥64 GB symbolic run.**

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py | sed -n '1125,1215p' && nl -ba research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py | sed -n '1013,1115p' && rg -n 'REP_INVARIANCE|INDEPENDENCE|sigma_W|MULTIGRADE' research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json | sed -n '1,160p'" in /var/projects/toy_physics
 succeeded in 0ms:
  1125	            jet_map={profile:m}
  1126	            for jet in v.atoms(sp.Derivative):
  1127	                suffix=''.join('d'+str(coordinates.index(q)+1) for q in jet.variables)
  1128	                jet_map[jet]=inputs.a('m1_profile_'+suffix)
  1129	            replacements[a]=v.xreplace(jet_map)
  1130	    for a,rhs in inputs.profiles.items():
  1131	        if a.name.startswith('mu_R_bg'):
  1132	            replacements[a]=rhs.xreplace(replacements)
  1133	    return expression.xreplace(replacements)
  1134	
  1135	
  1136	def representation_pullback(value,inputs,density):
  1137	    rho_br=inputs.density[(density,)][1]
  1138	    rho4=rho_br/inputs.values['W_bg']
  1139	    theta_shift=sum(wave_jet(inputs.a(f'u_{i+1}'))*inputs.dx(rho4,i) for i in range(3))/rho4
  1140	    e_shift=sum(wave_jet(inputs.a(f'u_{i+1}'))*inputs.a(f'W_bg_d{i+1}') for i in range(3))/inputs.values['W_0']
  1141	    theta_shift=sp.cancel(theta_shift.xreplace(inputs.profiles))
  1142	    e_shift=e_shift.xreplace(inputs.profiles)
  1143	    def pull(expression):
  1144	        replacements={}
  1145	        for atom in expression.atoms(AppliedUndef):
  1146	            point=atom.args[:3]
  1147	            translate=dict(zip(X,point))
  1148	            shift_theta=theta_shift if tuple(point)==X else inputs.at_source(theta_shift)
  1149	            shift_e=e_shift if tuple(point)==X else inputs.at_source(e_shift)
  1150	            if atom.func.__name__=='s11cc2Fieldtheta':
  1151	                replacements[atom]=atom+shift_theta.xreplace(translate)
  1152	            elif atom.func.__name__=='s11cc2FieldeW':
  1153	                replacements[atom]=atom+shift_e.xreplace(translate)
  1154	        @lru_cache(maxsize=65536)
  1155	        def compose(value):
  1156	            if not value.has(sp.Function('s11cc2Fieldtheta'),sp.Function('s11cc2FieldeW')):
  1157	                return value
  1158	            if value in replacements:
  1159	                return replacements[value]
  1160	            if isinstance(value,sp.Derivative):
  1161	                result=compose(value.expr)
  1162	                for variable in value.variables:
  1163	                    result=inputs.dx(result,X.index(variable)) if variable in X else sp.diff(result,variable)
  1164	                return result
  1165	            if isinstance(value,sp.Integral):
  1166	                return integral(compose(value.function),*value.limits)
  1167	            if not value.args:
  1168	                return value
  1169	            return value.func(*(compose(arg) for arg in value.args))
  1170	        return compose(expression)
  1171	    return tree(value,pull)
  1172	
  1173	
  1174	def regression_coordinates(value,inputs):
  1175	    tests={'v_theta_s11cb':field('s11cc2TestTheta',X,(0,0,0)),
  1176	           'v_e_W_s11cb':field('s11cc2TestE',X,(0,0,0)),
  1177	           'psi_L_s11cb':field('s11cc2TestPhi',X,(2,0,0))}
  1178	    tests.update({f'A_T_s11cb_{i+1}':field(f's11cc2TestA{i}',X,(2,0,0)) for i in range(3)})
  1179	    def transform(expression):
  1180	        replacements={}
  1181	        for a in expression.atoms(sp.Symbol):
  1182	            name=a.name
  1183	            if name in tests:
  1184	                replacements[a]=tests[name]
  1185	            elif name.startswith('u_T_'):
  1186	                replacements[a]=wave_jet(sp.Symbol(name.replace('u_T_','u_',1)),X,'TRANSVERSE')
  1187	            elif name.startswith('u_L_'):
  1188	                replacements[a]=wave_jet(sp.Symbol(name.replace('u_L_','u_',1)),X,'LONGITUDINAL')
  1189	            elif name.startswith('theta_probe'):
  1190	                replacements[a]=wave_jet(sp.Symbol(name.replace('theta_probe','theta',1)),X,'THETA')
  1191	            elif name.startswith('e_W_probe'):
  1192	                replacements[a]=wave_jet(sp.Symbol(name.replace('e_W_probe','e_W',1)),X,'E_W')
  1193	            elif name.startswith('phi_L_d'):
  1194	                p=field('s11cc2TrialPhi',X,(2,0,0))
  1195	                for i in re.findall(r'd([123])',name):
  1196	                    p=sp.diff(p,X[int(i)-1])
  1197	                replacements[a]=p
  1198	        return expression.xreplace(replacements)
  1199	    v=tree(value,transform)
  1200	    result={str(k):{str(n):e for n,e in row} for k,row in v}
  1201	    for direction,row in result.items():
  1202	        for target,value in row.items():
  1203	            sector='TRANSVERSE' if direction=='TRANSVERSE_TO_THICKNESS' else ('LONGITUDINAL' if target=='DIV_U' else target)
  1204	            row[target]=retained_shape(restricted(value,sector),inputs)
  1205	    return result
  1206	
  1207	
  1208	if __name__ == '__main__':
  1209	    run()
  1013	def run():
  1014	    started=time.monotonic()
  1015	    fold,audit=load_model(str(ROOT/'scripts/S11c_b_exports.py'),str(ROOT/'scripts/S11c_c1_exports.py'))
  1016	    closure_audit=check_consumer(fold,IMPORT_KEYS)
  1017	    lookup_audit=assert_lookups_equal_manifest(bind_inputs,fold,IMPORT_KEYS)
  1018	    inputs=lookup_audit['result']
  1019	    measurements={'fold':audit,'fold_rows':len(fold),'import_keys':sorted(IMPORT_KEYS),
  1020	        'lookups':sorted(lookup_audit['lookups']),'closure':sorted(closure_audit['closure'])}
  1021	    checkpoint=ROOT/('_measurements/S11c_c2_sympy_triage_guard_evidence.json' if os.environ.get('S11CC2_PACKAGE')=='TRIAGE' else '_measurements/S11c_c2_sympy_guard_evidence.json')
  1022	    checkpoint.write_text(json.dumps(measurements,indent=2))
  1023	    emit('DIMENSION_COEFFICIENT_BINDINGS',inputs.dimension_inference,inputs)
  1024	    if os.environ.get('S11CC2_PACKAGE')=='TRIAGE':
  1025	        case=(ANCHORINGS[0],DENSITIES[0])
  1026	        model=build_case(inputs,*case)
  1027	        for name,key in (('CLOSED_SLAB_OPERATOR','closed'),('CLOSED_COUPLING_KERNEL','closed_kernel'),('SELF_ENERGY_INCREMENT','increment')):
  1028	            emit(name,model[key],inputs,case)
  1029	        return
  1030	    exports={k:{} for k in EXPORT_ROOTS}
  1031	    baselines={}
  1032	    for anchoring in ANCHORINGS:
  1033	        for density in DENSITIES:
  1034	            case=(anchoring,density)
  1035	            model=build_case(inputs,*case)
  1036	            baselines[case]=model
  1037	            for quantity,obj in (
  1038	                ('CLOSED_SLAB_OPERATOR',model['closed']),
  1039	                ('CLOSED_SLAB_OPERATOR_TERM_ORIGINS',model['faces']),
  1040	                ('CLOSED_SLAB_OPERATOR_PARITY_BLOCKS',{
  1041	                    'FACE_SUM':tree(difference(model['faces'][1],tree(model['faces'][-1],lambda e:-e)),lambda e:e/2),
  1042	                    'FACE_DIFFERENCE':tree(difference(model['faces'][1],model['faces'][-1]),lambda e:e/2)}),
  1043	                ('CLOSED_COUPLING_KERNEL',model['closed_kernel']),
  1044	                ('CLOSED_COUPLING_KERNEL_TERM_ORIGINS',{s:extract(v,inputs) for s,v in model['faces'].items()}),
  1045	                ('SELF_ENERGY_CLOSED_EXTRACT_OPERAND',model['closed_kernel']),
  1046	                ('SELF_ENERGY_OPEN_EXTRACT_OPERAND',model['open_kernel']),
  1047	                ('SELF_ENERGY_INCREMENT',model['increment']),
  1048	                ('FOLD_SYMBOL_MAP',model['maps']),
  1049	            ):
  1050	                payload=emit(quantity,obj,inputs,case)
  1051	                export_key={'CLOSED_SLAB_OPERATOR':'s11cc2ClosedSlabOperator',
  1052	                            'CLOSED_COUPLING_KERNEL':'s11cc2ClosedCouplingKernel'}.get(quantity)
  1053	                if export_key:
  1054	                    exports[export_key][case]=payload
  1055	            # Regression-only bind of the predecessor's already-extracted kernel.
  1056	            regression=named(inputs.open_kernel[case],'COMPLETE_OPERATOR_BLOCKS')
  1057	            closed_regression=tree(regression,lambda e:inputs.physical_fields(e.subs(model['substitutions'],simultaneous=True)))
  1058	            # The predecessor uses jet trial/test symbols; align those explicitly
  1059	            # in regression_coordinates before comparing to the function ansatz.
  1060	            aligned=regression_coordinates(closed_regression,inputs)
  1061	            emit('ORDERING_EXTRACT_FIRST_OPERAND',aligned,inputs,case)
  1062	            emit('ORDERING_COMMUTATOR',difference(model['closed_kernel'],aligned),inputs,case)
  1063	            covectors,pairing,residual=traction_pairing(inputs,case,model)
  1064	            emit('TRACTION_MECHANICAL_CONTRIB',covectors,inputs,case)
  1065	            emit('TRACTION_SLAB_POWER_PAIRING',pairing,inputs,case)
  1066	            emit('TRACTION_SLAB_POWER_PAIRING_RESIDUAL',residual,inputs,case)
  1067	            _,corrupt_pairing,corrupt_residual=traction_pairing(inputs,case,model,flip=True)
  1068	            emit('TRACTION_SIGN_OPERAND',corrupt_pairing,inputs,case)
  1069	            emit('TRACTION_SIGN_RESIDUAL',corrupt_residual-residual,inputs,case)
  1070	            whole_atoms=inputs.values['dtn_operator'].atoms(sp.Symbol)
  1071	            used=cas(model['increment']).atoms(sp.Symbol)&whole_atoms
  1072	            emit('DTN_WHOLEFORM_DEPENDENCE',sp.Tuple(*sorted((a for a in used if not a.is_commutative),key=str)),inputs,case)
  1073	            flat=cases(inputs.values['dtn_flat_symbol'])[(anchoring,1)]
  1074	            bridge=model['maps'][1]['DTN_KERNEL_MATRIX'][0,0]
  1075	            emit('FLAT_SYMBOL_USAGE',{'REGRESSION':flat,'KERNEL_DIAGONAL':bridge,'RESIDUAL':bridge-flat},inputs,case)
  1076	            control(inputs,case,model,'ROUTING_MECHANICAL_ONLY',routing='MECHANICAL_ONLY')
  1077	            control(inputs,case,model,'ROUTING_TRACTION_CHANNEL',overrides={inputs.values['Lambda_X_0']:0})
  1078	            control(inputs,case,model,'ZERO_DTN',overrides={'ZERO_DTN':True})
  1079	            control(inputs,case,model,'LAMBDA_A_LIMIT',overrides={inputs.values['Lambda_A_0']:0})
  1080	            control(inputs,case,model,'IMPERMEABLE_LIMIT',overrides={inputs.values['Lambda_A_0']:0,inputs.values['Lambda_V_0']:0})
  1081	            control(inputs,case,model,'UNIFORM_LIMIT',overrides={inputs.eta:0,inputs.sigma:0,inputs.values['W_bg']:inputs.values['W_0']})
  1082	            control(inputs,case,model,'MU_R_FORM',mu_override=lambda e:modulus_form(e,inputs))
  1083	            measurements['last_case']=case
  1084	            checkpoint.write_text(json.dumps(measurements,indent=2))
  1085	    for anchoring in ANCHORINGS:
  1086	        emit('DENSITY_LIVE_MINUS_FROZEN',difference(baselines[(anchoring,DENSITIES[0])]['increment'],baselines[(anchoring,DENSITIES[1])]['increment']),inputs,(anchoring,))
  1087	    for density in DENSITIES:
  1088	        lab=baselines[(ANCHORINGS[0],density)]
  1089	        material=baselines[(ANCHORINGS[1],density)]
  1090	        mapped_closed=extract(representation_pullback(material['closed'],inputs,density),inputs)
  1091	        mapped_open=extract(representation_pullback(material['open'],inputs,density),inputs)
  1092	        mapped_closed=tree(mapped_closed,lambda e:retained_shape(e,inputs))
  1093	        mapped_open=tree(mapped_open,lambda e:retained_shape(e,inputs))
  1094	        mapped=difference(mapped_closed,mapped_open)
  1095	        emit('REP_INVARIANCE_LAB_OPERAND',lab['increment'],inputs,(density,))
  1096	        emit('REP_INVARIANCE_MATERIAL_OPERAND',mapped,inputs,(density,))
  1097	        emit('REP_INVARIANCE_RESIDUAL',difference(mapped,lab['increment']),inputs,(density,))
  1098	        corruption=build_case(inputs,ANCHORINGS[1],density,overrides={'FLIP_FACE_SLOPE':True})
  1099	        emit('REP_INDEPENDENCE_MATERIAL_NORMAL_OPERANDS',sp.Tuple(*normal_slope_control(inputs,ANCHORINGS[1],1)),inputs,(density,))
  1100	        emit('REP_INDEPENDENCE_MATERIAL_OPERAND',corruption['increment'],inputs,(density,))
  1101	        emit('REP_INDEPENDENCE_MATERIAL_RESIDUAL',difference(corruption['increment'],material['increment']),inputs,(density,))
  1102	        emit('REP_INDEPENDENCE_LAB_OPERAND',lab['increment'],inputs,(density,))
  1103	        lab_corruption=build_case(inputs,ANCHORINGS[0],density,overrides={'FLIP_FACE_SLOPE':True})
  1104	        emit('REP_INDEPENDENCE_LAB_NORMAL_OPERANDS',sp.Tuple(*normal_slope_control(inputs,ANCHORINGS[0],1)),inputs,(density,))
  1105	        emit('REP_INDEPENDENCE_LAB_CORRUPTED_OPERAND',lab_corruption['increment'],inputs,(density,))
  1106	        emit('REP_INDEPENDENCE_LAB_RESIDUAL',difference(lab_corruption['increment'],lab['increment']),inputs,(density,))
  1107	        emit('REP_INDEPENDENCE_MATERIAL_COMPANION_OPERAND',material['increment'],inputs,(density,))
  1108	    measurements['export']=publish(inputs,fold,{k:cas(v) for k,v in exports.items()})
  1109	    measurements['computation_lines']=COMPUTATION_LINES
  1110	    measurements['elapsed_seconds']=time.monotonic()-started
  1111	    measurements['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  1112	    checkpoint.write_text(json.dumps(measurements,indent=2))
  1113	
  1114	
  1115	def modulus_form(expression,inputs):
13:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), ",
22:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
31:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
40:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
49:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
58:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
67:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
76:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
85:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
94:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)), Tuple(Integer(0), Integer(0), Integer(2)), Tuple(Integer(0), Integer(1), Integer(0)), Tuple(Integer(0), Integer(1), Integer(1)), Tuple(Integer(0), Integer(1), Integer(2)), Tuple(Integer(0), Integer(2), Integer(0)), Tuple(Integer(0), Integer(2), Integer(1)), Tuple(Integer(0), Integer(2), Integer(2)), Tuple(Integer(0), Integer(3), Integer(0)), Tuple(Integer(0), Integer(3), Integer(1)), Tuple(Integer(0), Integer(3), Integer(2)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
103:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
112:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
121:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
130:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
139:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
148:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
157:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
166:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple()), ",
175:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), ",
184:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
193:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
202:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
211:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
220:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
229:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
238:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
247:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
256:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
265:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
274:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
283:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
292:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
301:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
310:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
319:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
328:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
337:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
346:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
355:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
364:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
373:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
382:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)), Tuple(Integer(0), Integer(0), Integer(2)), Tuple(Integer(0), Integer(1), Integer(0)), Tuple(Integer(0), Integer(1), Integer(1)), Tuple(Integer(0), Integer(1), Integer(2)), Tuple(Integer(0), Integer(2), Integer(0)), Tuple(Integer(0), Integer(2), Integer(1)), Tuple(Integer(0), Integer(2), Integer(2)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
391:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
400:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
409:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
418:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
427:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
436:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
445:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
454:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple()), ",
463:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), ",
472:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
481:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
490:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
499:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
508:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
517:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
526:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
535:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
544:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
553:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
562:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
571:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
580:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
589:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
598:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
607:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
616:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
625:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
634:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(0), Integer(2)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
643:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
652:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
661:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
670:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)), Tuple(Integer(0), Integer(0), Integer(2)), Tuple(Integer(0), Integer(1), Integer(0)), Tuple(Integer(0), Integer(1), Integer(1)), Tuple(Integer(0), Integer(1), Integer(2)), Tuple(Integer(0), Integer(2), Integer(0)), Tuple(Integer(0), Integer(2), Integer(1)), Tuple(Integer(0), Integer(2), Integer(2)), Tuple(Integer(0), Integer(3), Integer(0)), Tuple(Integer(0), Integer(3), Integer(1)), Tuple(Integer(0), Integer(3), Integer(2)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)))), ",
679:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
688:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
697:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
706:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
715:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
724:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
733:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
742:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple()), ",
751:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), ",
760:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
769:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
778:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
787:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
796:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
805:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
814:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
823:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
832:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
841:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
850:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
859:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
868:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
877:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
886:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
895:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
904:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
913:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
922:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(0), Integer(2)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
931:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
940:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
949:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
958:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)), Tuple(Integer(0), Integer(0), Integer(2)), Tuple(Integer(0), Integer(1), Integer(0)), Tuple(Integer(0), Integer(1), Integer(1)), Tuple(Integer(0), Integer(1), Integer(2)), Tuple(Integer(0), Integer(2), Integer(0)), Tuple(Integer(0), Integer(2), Integer(1)), Tuple(Integer(0), Integer(2), Integer(2)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
967:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
976:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
985:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
994:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1003:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1012:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1021:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1030:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple()), ",
1039:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)))), ",
1048:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1057:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1066:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1075:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1084:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1093:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1102:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1111:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1120:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1129:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1138:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)))), ",
1147:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1156:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1165:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1174:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1183:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1187:    "PY_S11CC2_REP_INVARIANCE_LAB_OPERAND_RHO4_CONSTANT": {
1192:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1196:    "PY_S11CC2_REP_INVARIANCE_MATERIAL_OPERAND_RHO4_CONSTANT": {
1201:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1205:    "PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHO4_CONSTANT": {
1210:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1214:    "PY_S11CC2_REP_INDEPENDENCE_MATERIAL_NORMAL_OPERANDS_RHO4_CONSTANT": {
1219:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)))), ",
1223:    "PY_S11CC2_REP_INDEPENDENCE_MATERIAL_OPERAND_RHO4_CONSTANT": {
1228:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1232:    "PY_S11CC2_REP_INDEPENDENCE_MATERIAL_RESIDUAL_RHO4_CONSTANT": {
1237:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1241:    "PY_S11CC2_REP_INDEPENDENCE_LAB_OPERAND_RHO4_CONSTANT": {
1246:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1250:    "PY_S11CC2_REP_INDEPENDENCE_LAB_NORMAL_OPERANDS_RHO4_CONSTANT": {
1255:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(0), Integer(0), Integer(1)))), ",
1259:    "PY_S11CC2_REP_INDEPENDENCE_LAB_CORRUPTED_OPERAND_RHO4_CONSTANT": {
1264:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1268:    "PY_S11CC2_REP_INDEPENDENCE_LAB_RESIDUAL_RHO4_CONSTANT": {
1273:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1277:    "PY_S11CC2_REP_INDEPENDENCE_MATERIAL_COMPANION_OPERAND_RHO4_CONSTANT": {
1282:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1286:    "PY_S11CC2_REP_INVARIANCE_LAB_OPERAND_RHOBR_CONSTANT": {
1291:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(0), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1295:    "PY_S11CC2_REP_INVARIANCE_MATERIAL_OPERAND_RHOBR_CONSTANT": {
1300:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(0)), Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(0)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1304:    "PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHOBR_CONSTANT": {
1309:      "multigrade": "Tuple(Str('MULTIGRADE'), Tuple(Tuple(Integer(1), Integer(0), Integer(1)), Tuple(Integer(1), Integer(1), Integer(1)))), ",
1313:    "PY_S11CC2_REP_INDEPENDENCE_MATERIAL_NORMAL_OPERANDS_RHOBR_CONSTANT": {

exec
/bin/bash -lc "rg -n 'project_s11c_c_state|S11c-c2|c2 STEP|next_prompt_c2_wolfram|8f3a017f|aa76105a' memory MEMORY.md . 2>/dev/null | head -n 260" in /var/projects/toy_physics
 succeeded in 136ms:
splay': '(1, 0, 0)', 'value_kind': 'COMPUTED_OBJECT', 'class': 'STRUCTURAL', 'step': 'S11c-c2', 'route': 'F9A_ABSENT'},
./research/pde_ledger_v3/scripts/S11c_c2_exports.py:79:'s11cc2Y3': {'value': _restore("Symbol('s11cc2Y3', real=True)"), 'display': 's11cc2Y3', 'value_kind': 'COMPUTED_OBJECT', 'class': 'COORDINATE', 'step': 'S11c-c2', 'route': 'F9A_ABSENT', 'dimension_key': 's11cc2Y3Dimension'},
./research/pde_ledger_v3/scripts/S11c_c2_exports.py:80:'s11cc2Y3Dimension': {'value': _restore('ImmutableDenseMatrix([[Integer(1)], [Integer(0)], [Integer(0)]])'), 'display': '(1, 0, 0)', 'value_kind': 'COMPUTED_OBJECT', 'class': 'STRUCTURAL', 'step': 'S11c-c2', 'route': 'F9A_ABSENT'},
./STATUS.md:5:## ⭐⭐⭐ c2 STEP A (physics review) + STEP B/C (export repair) DONE — fold physics SOUND; export 60 MB → 21.4 MB; both committed; NEXT = light spec fix → blind Wolfram engine (2026-09-06)
./STATUS.md:6:**(A) c2 PHYSICS 2-leg review DONE — fold physics SOUND, 0 confirmed defects, committed `8f3a017f`.** Real independent legs (Codex-written script → **fresh Claude agent + Grok**, identical prompt, blindness-by-absence). Both agree the fold WIRING + A/C/D1–D6 are SOUND (δp not J_s, operator-inverse, kernel bridge, V_s→face_velocity, computed jets, ε-strip; FORM ablations respond). **Grok flagged F/G/E as CONFIRMED defects — ALL FALSE POSITIVES**, refuted by my own computation (rule 13; `_measurements/S11c_c2_physics_review_adjudication.md` + `…_adjudication_verify_{F,EG}.py/.out`): **F** uniform-limit — the genuine closure-induced coupling **DECOUPLES** (Integral integrand ≡0, `.doit()`→0 all 4 cases; Grok saw `Integral(...)` without evaluating it — its recurring pattern, cf. c1 "tangential freeze"); **G** adjointness — the builder never claimed adjoint-by-construction, both blocks emitted, the self-energy is **DIRECTIONAL/one-way** (reverse block ≡0, expected since δp lives only in the θ/mech rows) — dissipative directional, not a defect; **E** N6 — leading O(ε)/O(εη) rep-invariance **HOLDS** (σ_W→0 annihilates the whole residual), σ_W-sector remnant DEFERRED. **Carry-forwards (step-record, ⛔ NOT defects):** F wording "genuine coupling decouples" NOT "increment vanishes" (→ light §5e/§3c spec fix owed); E σ_W-sector deferred; G directional self-energy; the 2 S11c-b sign conventions → the WL comparator.
./STATUS.md:7:**(B/C) EXPORT REPAIR DONE + committed `aa76105a` — 60,516,900 → 22,441,522 bytes (63% smaller), publication-only, physics untouched.** Directive gated `a5f7a06c` (2 decision legs Codex sol xhigh + Grok, one pass, the R2 guard-spec blocker folded — `_measurements/S11c_c2_export_repair_decision_{gate,codex.md,grok.txt}`). astra build (`gpt-6-astra` high, setsid-detached) **BEHAVED — no overstep** (the fixed §7 fence held) but SLOW (~50 min / 3 full runs). Repair: drop `SelfEnergyIncrement` from `EXPORT_ROOTS`+the export-key map → EMIT-only; new `publication_compact` (transparent-factors **only** the per-case VALUE leaves, keeps the cased tree, preserves Integrals + the pole set — **reverts** if factoring would change it — and raises if any temp Dummy escapes); a **separate** strict-recursive Integral-aware semantic guard (srepr roundtrip kept alongside); `display`→root-name. **2 re-review legs** (astra-written → fresh Claude + Grok): **Grok CLEAR** (guards bite — `out+1` compaction + membership re-add/drop all hard-stop). The **fresh-Claude leg WEDGED** on reaped background waiters and produced **2 FALSE POSITIVES** — both artifacts of its **cross-process** re-run (fresh Integral bound-dummies vs the stored compact: `k·d1·X − k·d2·X`). **My rule-13 verify** (same-process, `_measurements/S11c_c2_export_repair_adjudication_verify.py/.out`): FAITHFUL **5/5**, GUARD BITES **5/5** (doubling, marker-add, in-integrand) → repair SOUND; adjudication `_measurements/S11c_c2_export_repair_rereview_adjudication.md`.
./STATUS.md:9:**NEXT:** (a) light **§5e/§3c spec-wording clarification** ("genuine coupling decouples", not "increment vanishes"; review-until-clear); (b) the **blind Wolfram engine** (imports nothing, re-derives) → its **2 legs** → the **c2 T7 comparator + reconcile** → the **c2 step record** (NO per-substep card — one S11c roll-up after S11c-e). Live handoff: [[project_s11c_c_state]]. ⛔ c1 STANDS. ⚠ `MEMORY.md` ~20 KB — compact next fresh session.
./STATUS.md:15:**NEXT:** (a) run the **REAL independent legs** (Codex-written → **fresh Claude agent + Grok**, orchestrator-launched) on the **physics** (fold map, ε, kernel bridge, N6 rep-invariance, the uniform-limit `delta_p_minus` concern astra self-flagged, the 6 re-adjudications, mandatory FORM ablation + one-sided corruption); (b) **export repair** (membership + held/factored — a script re-run); (c) adjudicate → re-review until clear → commit; then the **blind Wolfram engine** → **T7 comparator + reconcile** → **c2 step record**. Live handoff: [[project_s11c_c_state]]. ⚠ `MEMORY.md` ~20 KB — compact next fresh session. **The clause below is SUPERSEDED — kept for the artifact map.**
./STATUS.md:18:Three things landed after the c1 close. **(1) c2 SHARED_PHYSICS spec (the self-energy fold) GATED + committed `16849fc6`** — 2 decision legs (Codex sol xhigh + Grok) both ⚙-computational + CONVERGED, v1 NOT sound (central **isolation claim FALSE** — the substitution increment does NOT cancel S11c-b's two unvalidated SIGN CONVENTIONS; they multiply the substituted δp_s slots, both legs computed `2·carrier·incr`), folded ONE pass → **v2** (9 findings: real θ-row has `J_s` ELIMINATED ⇒ substitute **closed δp_s+w-jets** NOT J_s; `Λ_X` traction-only; **assembled two-face** operator; **3→6 re-adjudications**; positional `load_model` + `s11c_c1_` keys; `Z→0≠Λ_A=0`). Record `_measurements/S11c_c2_spec_review.md`. **(2) c1 SPEC RETRO-REVIEW DONE** (user's *'full correct path'* — the c1 spec was folded once, never re-legged): 2 legs DISAGREED (Grok CLEAR / Codex BLOCK-reopen), adjudicated (rule 13, `_measurements/S11c_c1_spec_retro_review_adjudication.md`) — ⭐ **c1's ENGINES/EXPORTS STAND, ⛔ do NOT reopen** (PY carries density OPAQUELY, 0 derivatives; re-binding to `background_density_map` recovers the O(εη) channel exactly → c2 is **specified** to recover it, c2 not yet built; Codex-sol verify CONFIRMED the override), but record corrections OWED — **3 MUST-level** (seal-5 "harmless O(η²)"→**O(εη)-recoverable-representational**; grazing→`NOT_ESTABLISHED`; independence-scoping) **+ lower-severity carry-forward** (energy-residual orientation, `h_s`/DtN-vs-impedance terminology, density-as-mult-operator, `K_a` Hermitian, η²/ησ_W/σ_W² caveat, drain-projection, flat-`Z₀` leakage; full list in the adjudication record). **(3) CLAUDE.md RESTRUCTURED + committed `3da095c6`** (17 rules→M/E/G/S + an **artifact→review-discipline table**; the fix: **specs=review-until-clear**, not fold-once; ⚠ the preservation legs caught + removed a self-introduced **FABRICATION** in the rewrite). **NEXT = (a) the 3 c1 RECORD corrections (2 legs, review-until-clear; ⛔ NOT engine reopens), THEN (b) c2 SymPy BUILD DIRECTIVE** (own 2 decision legs BEFORE the astra build; determine `IMPORT_KEYS` root set + fold symbol map vs the real files). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = c2) is SUPERSEDED — kept for the artifact map.**
./STATUS.md:21:Both c1 engines per-engine SOUND (clause below). The **T7 comparator** was gated (`84686a54`, 2 decision legs, 12 folds) → built → reviewed baseline `7141e6ad` → scoped repair (astra) `704308af`, all re-review legs CLEAN = **SOUND** (measurement-only, prints/decides nothing). ✅ **CROSS-ENGINE RECONCILE RUN** (family-scoped, 46/50 families, 30 GB, peak ~317 MB) + **c1 STEP RECORD committed `072d0b75`** (`steps/S11c_c1_curved_bulk_closure.md`; record `_measurements/S11c_c1_comparator_reconcile.md`; reproducible `_measurements/S11c_c1_reconcile_reproduce.py`). Method = the **staged representational bridge** ([[feedback_reconcile_representational_bridge]]): comparator PRINTS raw (nonzero≠disagree), I apply identifications read from the engine sources (⛔ not tuned to zero; adversarial corruptions bite) and test collapse. ⭐ **EARNED cross-engine AGREE:** the two-momentum **DtN KERNEL** — all 4 (anchoring×face) cases EXACT ZERO off-diagonal (Stage-2 residual collapsed to the dispersion relation itself) — and the **face-response COEFFICIENTS** (LAB_HELD off-diagonal + MATERIAL on-diagonal); seals 1–4 AGREE. **2 step-record legs (Grok + Codex sol xhigh) SCOPED an earlier "closure-wide AGREE" over-claim** — folded once (`_measurements/S11c_c1_step_record_review_{grok,codex_sol}.txt`). ⛔⛔ **UNDECIDED — ⛔ NOT cross-engine-closed (c2 must NOT inherit these):** (5) **bg-density = a SURFACED rule-17 freeze** (PY bare constant vs WL live field `(ρbr/W₀)WBg(x)`; functional dependence agrees + 0 derivatives, but a scalar≠a field globally — c2 re-adjudication MANDATORY); the raw **`dtn_operator`** whole-form (kernel-AGREE does NOT extend to it); **ENERGY** (closed-form vs far-field integral); an off-diagonal **flat-resolvent leg-labeling**. ⚠ **4 giants + full symbolic residual UNMEASURED — DEFERRED ≥64 GB** (`DEFERRED_HEAVY_RUNS.md`, ⛔ not pre-adjudicated). No cross-engine DISAGREEMENT found. **The kernel is cross-engine closed; the rest is UNDECIDED/deferred.** ⚠ Codex-sol correction-verify pass was run on the fold (in flight/complete). **NEXT = c2** (self-energy fold, a CODE build → gpt-6-astra, own gated directive with 2 decision legs BEFORE the build; carry the mandatory seal-5 + `dtn_operator` re-adjudication). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = comparator) is SUPERSEDED — kept for the artifact map.**
./STATUS.md:24:The export migration + delta rebuild (clause below) is **DONE** (`f90e7630`); the c1 SymPy engine is **DONE** (`d6e16471`, 906 KB own-rows delta export); the c1 **blind Wolfram** engine is **DONE** (`c2b0f9dd`→`e139bc61`→`13f0bd2c`→ **repair-2 `dd34d564`**, fully gated, both re-review legs CLEAR). ✅ **A FULL REMEDIATION is COMPLETE** (`_measurements/S11c_c1_wl_remediation_plan.md`): the WL repair DIRECTIVE had skipped its rule-7 decision legs; the remediation re-ran repair-2 **WITH** its decision legs, which **caught that a first-draft fix targeted a non-defect** (would have corrupted the correct `PERMEABLE_PORT_HERMITIAN`) — the gate working. Both c1 `.out` are committed (`4a14100a`, annex+GIN; WL SHA-256 `d63e616786…`). ⛔ **Both c1 engines are BUILT + reviewed — do NOT rebuild.** **NEXT = the c1 T7 cross-engine comparator** — frozen contract `directives/S11c_c1_SHARED_PHYSICS.md:580-587` (N8 inherited verbatim; ⛔ point at it, don't paraphrase): it **computes and prints, deciding nothing** (join by object name with the S11c-a axis-typed keys, injective `mechanical_lower_camel`, paired residual operands + per-family accounting, repoint ablation); the three-valued AGREE/DISAGREE/UNDECIDED characterization is the POST-RUN adjudication, ⛔ not a script verdict; ⛔ no representational fold pre-registered. Reuse `scripts/S11c_a_cross_engine_comparator.py`. Orchestrator-written → Codex+Grok decision legs BEFORE the build → build → 2 legs. Known representational residuals to adjudicate AFTER the run (never pre-registered): μ_θ, the ω-assumption artifact, the two-momentum leg names, and the new proportional-4 parity NIT. ⚠ The full cross-engine residual may need the ≥64 GB box (this box is 30 GB) — build per-engine here, defer the heavy residual if it OOMs. Live handoff: [[project_s11c_c_state]]. **The two clauses below (c1 build PAUSED / NEXT = migration) are SUPERSEDED — kept for the artifact map.**

exec
/bin/bash -lc "stat -c '%n %s bytes' /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/MEMORY.md 2>&1 || true && nl -ba /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md | sed -n '1,260p' && rg -n 'project_s11c_c_state|8f3a017f|aa76105a|N6|sigma|σ|export repair' /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/MEMORY.md 2>/dev/null | sed -n '1,160p'" in /var/projects/toy_physics
 succeeded in 0ms:
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md 31326 bytes
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/MEMORY.md 21303 bytes
     1	---
     2	name: project_s11c_c_state
     3	description: "S11c-c (curved bulk closure): 2-way split c1/c2. c1 DONE: both engines per-engine SOUND + T7 comparator SOUND (7141e6ad→704308af) + ✅ CROSS-ENGINE RECONCILE + c1 STEP RECORD committed 072d0b75 → correction-verify fold b8a83206 (2 step-record legs + 1 verify leg SCOPED the verdict twice). EARNED cross-engine AGREE: two-momentum DtN KERNEL (all cases EXACT ZERO, Codex-confirmed) + δp_s/J_s response coefficients. ⛔ UNDECIDED (c2 must NOT treat as closed): seal-5 density (surfaced rule-17 freeze), t_s traction (WL 4-vec vs PY scalar), raw dtn_operator whole-form, ENERGY, flat-leg; 4 giants+full symbolic residual DEFERRED ≥64 GB. c2 SPEC (self-energy fold) GATED + folded to v2 `16849fc6`: 2 decision legs (Codex sol xhigh + Grok) both ⚙-computational + CONVERGED, v1 NOT sound — central isolation claim FALSE (increment drops bulk/kinetic base but the face-force + #90 closure-fold SIGN CONVENTIONS multiply the substituted δp_s slots and leak; both legs computed 2·carrier·incr), fold undefined vs real θ-row (J_s ALREADY eliminated #90 → substitute closed δp_s+w-jets NOT closed J_s; Λ_X is traction-only), re-adj set 3→6. c1 SPEC RETRO-REVIEW DONE (user's 'full correct path on c1' — the spec was folded once, never re-legged): 2 legs, DISAGREED (Grok CLEAR / Codex BLOCK-reopen), adjudicated (rule 13, `_measurements/S11c_c1_spec_retro_review_adjudication.md`) — ⭐ c1's ENGINES/EXPORTS STAND (do NOT reopen: PY carries density OPAQUELY, 0 derivatives; re-binding to background_density_map recovers the O(εη) channel exactly → c2 is SPECIFIED to recover it per c2-v2 §3d.1, c2 not yet built; a fresh Codex-sol verify CONFIRMED the override), but record corrections OWED — 3 MUST-level (seal-5 'harmless O(η²)'→O(εη)-recoverable-representational; grazing→NOT_ESTABLISHED; independence-scoping) + lower-severity carry-forward (energy-residual orientation; h_s/DtN-vs-impedance terminology; density-as-mult-operator; K_a Hermitian; η²/ησ_W/σ_W² caveat; drain-projection; flat-Z₀ leakage — full list in `_measurements/S11c_c1_spec_retro_review_adjudication.md`). CLAUDE.md RESTRUCTURED+committed `3da095c6` (specs=review-until-clear now explicit). ✅ STEP 1 c1 RECORD corrections DONE+committed `efae0337` (2 rounds review-until-clear, both legs CLEAR; folded grazing-inverse-misattribution [N⁻¹/Z singular, permeable resolvent REGULAR] + energy-sign; c1 STANDS 4× confirmed). ✅ STEP 2 c2 SymPy BUILD DIRECTIVE (thin, points at spec 16849fc6) GATED+committed `c1de32b0`: 2 decision legs (both ran ledger_fold), wiring SOUND (18 closure-covering roots+coupling_kernel 19th; exact-lookup IMPORT_KEYS rule; bare-vs-prefixed face_response hazard), 7 symbol-map defects folded one-pass (DELTA_P not PRESSURE; w-jets COMPUTED; ε-normalization avoid O(ε²); V_s→face_velocity not kinematic_balance; dtn_operator→kernel bridge; coupling_kernel promoted; clauses verbatim). ⚠ astra build RAN (19min EXIT0, faithful physics per builder report) but (a) OVERSTEPPED — read /build skill, ran its own INVALID self-reviews + ~40-file sprawl (quarantined to ext_logs/c2_builder_overstep; [[feedback_builder_directive_no_orchestrator_process]]), (b) OVER-EXPORTED — 60MB = 3 fully-expanded closed objects (CouplingKernel 24.7 / SlabOperator 19.9 / Increment 15.9). Directive fixed `1ae6c336` + refined post-Codex-verify (2 defects: §7 builder-bounds; §5 EMIT≠EXPORT — export BOTH closed operators S11c-d binds: CLOSED_SLAB_OPERATOR [full closed op, d needs the closure-modified DIAGONAL for its spectrum/resolvent-poles] + CLOSED_COUPLING_KERNEL [off-diagonal, for Born/mixing], per d's DECLARED scope [S11c_decisions.md:47,52,83] — ⛔ NOT coupling-kernel-only [my earlier guess, Codex-corrected]; TRANSPARENT-factored not sp.expand/UnevaluatedExpr + a casewise `canonicalize(expanded−decode(compact))==0` check; Increment+term_origins+§3d+§5 controls EMIT-only→.out). Codex compact-prep verify (`_measurements/S11c_c2_compact_prep_verify_codex.txt`): physics-legs + export-repair are SEPARABLE (physics ends script:951, publish separate :952) → physics legs first on frozen script/stdout, then PUBLICATION-ONLY export repair (rebuild only if it crosses into construction code). ✅ STEP A c2 PHYSICS 2-leg review DONE (fresh Claude agent + Grok, identical prompt, blindness-by-absence) → fold physics SOUND, 0 confirmed defects, committed `8f3a017f`. Grok flagged F/G/E — ALL FALSE POSITIVES refuted by my own computation (rule 13, `_measurements/S11c_c2_physics_review_adjudication.md` + verify_F/EG.py/.out): F uniform-limit — the genuine closure-induced coupling DECOUPLES (Integral integrand ≡0, .doit()→0 all 4 cases; Grok saw Integral(...) without evaluating it, its known pattern); G adjointness — builder never claimed adjoint-by-construction, both blocks emitted, self-energy is DIRECTIONAL/one-way (reverse block ≡0, expected since δp lives only in θ/mech rows); E N6 — leading O(ε)/O(εη) rep-invariance HOLDS (σ_W→0 kills the whole residual), σ_W-sector remnant DEFERRED. CARRY-FORWARDS (step-record, ⛔ NOT defects): F wording 'genuine coupling decouples' NOT 'increment vanishes' → light §5e/§3c spec fix owed; E σ_W deferred; G directional self-energy; the 2 S11c-b sign conventions → WL comparator. ✅ STEP B/C EXPORT REPAIR DONE + committed `aa76105a`: 60MB→21.4MB (publication-only, physics untouched). Directive gated `a5f7a06c` (2 decision legs Codex sol xhigh+Grok, R2 blocker folded 1 pass; `_measurements/S11c_c2_export_repair_decision_{gate,codex,grok}`). astra build BEHAVED (NO overstep — §7 fence held) but SLOW (~50min, 3 full runs). Repair = drop increment→EMIT-only + `publication_compact` (transparent-factor VALUE leaves only, keep cased tree, preserve Integrals/poles, REVERTS if pole-set changes, raises if a temp Dummy escapes) + strict-recursive Integral-aware semantic guard (kept the srepr roundtrip separately) + display→root-name. 2 re-review legs (astra-written→fresh Claude+Grok): Grok CLEAR (guards bite: out+1 + membership re-add/drop); the Claude leg WEDGED on reaped bg-waiters + produced 2 FALSE POSITIVES (cross-process re-run → fresh-Integral-dummy artifacts `k·d1·X−k·d2·X`). My rule-13 verify (SAME-process, `_measurements/S11c_c2_export_repair_adjudication_verify.py/.out`): FAITHFUL 5/5, GUARD BITES 5/5 (x2, +marker, in-integrand) → repair SOUND. ⚠ LESSON: Integral-bearing operators need SAME-PROCESS comparison (cross-process re-run mismatches bound dummies → spurious FAIL). ⚠ background tasks REAPED all session (setsid+marker+watcher for the astra build; harness watchers reaped ~30s each — even trivial loops) [[feedback_background_tasks_can_die_spuriously]]. ⭐ astra assessment (user asked): KEEP astra for CODE — this build's code was correct+defect-free under adversarial review; the earlier c1 overstep was MY directive-fence defect (fixed), not astra; only weakness = speed (redundant full re-runs). NEXT = light §5e/§3c spec-wording fix (review-until-clear) → blind WL engine (imports nothing, re-derives) → its 2 legs → c2 T7 comparator+reconcile → c2 step record (NO per-substep card). ⛔ c1 STANDS. ⚠ MEMORY.md ~20KB — compact next fresh session."
     4	metadata: 
     5	  node_type: memory
     6	  type: project
     7	  originSessionId: 53620ffb-59f9-482d-b804-aef04f767516
     8	  modified: 2026-09-06T06:27:13.955Z
     9	---
    10	
    11	**S11c-c = the curved-interface bulk closure** — the third sub-step of the S11c staged family. It solves the
    12	perturbed **curved two-face outgoing bulk** acoustic problem for the nonlocal **DtN/impedance** operator and
    13	folds it back into the S11c-b slab operator (whose θ-row still carries a **symbolic `δp_s`** in
    14	`closure_shape_deriv`). This is the generalization of S11b's flat-face **B0b** (impedance `Z`, three regimes,
    15	per-face inertial loading) and **B0c** (permeable face response) to the tilted faces S11c-a shape-differentiated.
    16	
    17	**⛔ 2-WAY SPLIT (user chose it 2026-09-03; N2 refinement of the S11c-c decision-list row):**
    18	- **c1 = curved bulk closure** — the nonlocal DtN/impedance operator + the permeable B0c face response +
    19	  Fredholm loci + dissipation audit. **Exports** the closed face response `(δp_s,J_s,t_s)(V_s,μ_θ)` + the DtN
    20	  operator (for c2). Tag prefix `S11CC1_`.
    21	- **c2 = self-energy fold** — fold the closed response into `S11CB_SLAB_OPERATOR`, **re-extract** the closed
    22	  off-diagonal kernel from the CLOSED full operator, → the coupled nonlocal **self-energy operator**. NOT yet
    23	  authored (author it against c1's reviewed exports; tag prefix set at c2 authoring).
    24	
    25	**STATE: c1 spec DONE + COMMITTED `db5cbf88`** (`directives/S11c_c1_SHARED_PHYSICS.md`; the committed spec is the
    26	build input — the SymPy engine's `BUILD_INPUT_DIGESTS` pins it; the c1 SymPy engine + 906 KB delta export are now
    27	committed (`d6e16471`); the first build's 145 MB accumulated products were deleted and replaced by the own-rows
    28	delta under the new architecture). The
    29	rule-7 gate worked: v1 MONOLITH spec drew **~13 verified defects** (2 legs, saved SymPy
    30	derivations) → split decision; c1 got its own 2 legs → both **"safe after folds, NOT a re-author"** → **10 folds
    31	applied one pass** (rule 7: fold and go, NOT re-legged). Records: `_measurements/S11c_c_spec_review.md` (v1 +
    32	split), `_measurements/S11c_c1_SHARED_PHYSICS.md` (c1 folds), leg derivations under `directives/_legs/`.
    33	
    34	**Locked c1 physics (the folds):**
    35	- **DtN is a two-momentum OPERATOR** `Z_s(ω;k,k′)` — composition `N₀∘M_{h_s}∘N₀` + `Div(h∇)` + `κ²h`, AND the
    36	  two-momentum kernel with **both** branch legs `q_out(k),q_out(k′)` explicit (`Z_1∝Ŵ_bg(k−k′)·[q(k)q(k′)+k·k′−ω²/c_s0²]/(q(k)q(k′))`).
    37	  ⛔ NO single-`k` multiplier, ⛔ NO left-quantized single-leg `a(x,k)=W_bg(x)σ(q(k))` (both delete the mode
    38	  mixing — the rule-17 freeze). Rigid-shift (`k=k′`, `Ŵ(0)`) cancellation is a NAMED residual.
    39	- **Dissipation = 3 distinct objects:** (1) bulk-radiation Hermitian part `H_a[Z]=(Z+Z†ₐ)/2`; (2) the two-port
    40	  power-conjugate permeable Hermitian form on `(V,μ_s)→(δp+Λ_X𝒜,J)` (S11b `:705-717`) — the bare `Z` carries NO
    41	  `Λ_I`, so it audits radiation ONLY; (3) an INDEPENDENT energy balance: face operand = true-area **traction**
    42	  pairing `½Re Σ_s∫ a_s^0 t_s·v_face,s*` (NOT the toothless typed `½Re(δp V*)`, which = the bulk flux and never
    43	  sees `t_s`), bulk operand = outgoing **far-field Poynting flux** (NOT δp at the face). Reactive part `K_a`; no
    44	  sign-definiteness on the flat evanescent nullspace at 1st order (O(η²) leakage → S11c-e).
    45	- **Controls:** `Σ_E`-advection and `μ_R,bg` are **structurally absent from c1's μ_θ-symbolic objects** →
    46	  reserved for c2 (mutate `ρ_br,bg` in `μ_s` if a c1 density test is wanted, NOT `Σ_E`); rep-invariance route 2 =
    47	  **Hanzawa/layer-potential** (⛔ NOT the global scaling `w′=[w−ζ_c]/[W_bg+δW]` — secular at infinity → wrong
    48	  radiation branch; drop N4/MATERIAL tags on the bulk DtN); zero-jet target = the **UNMODIFIED** S11b half-space
    49	  `Z` (⛔ NOT B0b re-solved at gap `W̄₀(1+η)` — that reproduces the cavity error it catches); both one-leg branch
    50	  freezes with a fixed `(k` output`,k′` input`)` convention.
    51	- Fredholm noninvertibility condition (formal) vs the algebraic locus protocol (flat/finite-dim only; profile
    52	  resolvent → S11c-d); non-uniform grazing domain (grazing = strict `v_bulk_normal_0=0`; boundary-layer
    53	  exclusion + subsonic condition away from it). N11a: rest-frame limit, ⛔ NO convective operator.
    54	
    55	✅ **c1 MIGRATION DONE + first build reviewed (2026-09-04)** ([[project-export-bind-closure-architecture]]):
    56	- **Migration committed `f90e7630`** — c1 spec §1/§7 + build directive + N1 on the fold/own-rows-delta topology;
    57	  2-leg reviewed (Codex+Grok), reconciled **44-root `IMPORT_KEYS`** (removed c2-reserved `mu_R_bg`/`m1_profile`/
    58	  `mu_R` + slab-side `rho_4D_bg_*`/`e_W_bg`; added `omega` + `epsilon_shape`). Grounded on the real base
    59	  (44 roots → 193-row closure).
    60	- **First c1 SymPy build DONE + 2-leg-reviewed; baseline `65afa1cd`.** Export = **906 KB / 44-row OWN-ROWS
    61	  DELTA** (vs the 145 MB accumulated first attempt — the redesign's payoff, validated end-to-end: check_consumer
    62	  accepts it, F9c predecessor preserved, zero fold overwrites). ⛔ The build needed the **DETACHED-launch
    63	  workaround** (`setsid`+marker+Monitor) — the harness reaped 3 `run_in_background` tasks at ~87s; NOT OOM (1.4 GB
    64	  peak, 23 GB free) [[feedback-background-tasks-can-die-spuriously]]. CORE PHYSICS 2-leg-confirmed SOUND
    65	  (two-momentum DtN kernel — Grok's "tangential freeze" adjudicated a FALSE POSITIVE, `_measurements/c1_tangential_kernel_adjudication.py`;
    66	  operator-inverse response; Λ_X-only-traction; opaque μ_θ; delta topology; lookup smoke-test bites). **5
    67	  IMPLEMENTATION defects, all EMIT-ONLY controls** (record `_measurements/S11c_c1_build_review.md`): (1)+(2) both
    68	  legs — the energy-balance route + the §5a rep-invariance control are **A−A** (fix: real far-field Poynting from
    69	  φ at |w|→∞; real Hanzawa/layer-potential route); (3) Claude — on-shell `xreplace` misses (`sp.factor` breaks
    70	  the Add-key) → rigid-shift off-shell + §5d contaminated; (4) Grok — Hermitian/port fed `Z^(1)`-only (vanish at
    71	  (η,σ_W)=0); (5) Grok — `assert_delta_is_minimal` theater (own_closure = the delta).
    72	- ✅ **REPAIR CLEARED + committed `d6e16471`** — all 5 controls now BITE (2 fresh legs, both CLEAR/SOUND,
    73	  independent derivations + /tmp ablations: R1 far-field-Poynting energy route, R2 Hanzawa layer-potential 2nd
    74	  route, R3 on-shell reduction, R4 Hermitian Z_0+Z_1, R5 intended-closure minimum-mode). Core byte-identical to
    75	  `65afa1cd`; export = 1-line self-digest diff, all 44 rows srepr-equal. Record `_measurements/S11c_c1_repair_directive.md`.
    76	  ⚠ **Carry to the c1 step record:** the `ZERO_JET_RESIDUAL` **omega-assumption artifact** (S11b plain `omega` vs
    77	  c1 `omega(real=True)` — numerically 0, NOT thickness dep; the ω-assumption difference is a **computed residual to
    78	  adjudicate AFTER the run** (⛔ NOT pre-canonicalized into the comparator — rule 5 / contract
    79	  `S11c_c1_SHARED_PHYSICS.md:585`; [[project_dual_engine_symbol_transliteration]]); a dead `route="HANZAWA"` branch (inert).
    80	- ⭐ **c1's SymPy per-engine side is DONE.** NEXT = the **blind Wolfram engine** (imports nothing, re-derives the
    81	  §§1–2 setup + T-substrate + μ_θ operand from the sibling specs — `S9_export_chain_rebuild_directive.md:16-18`
    82	  is the ONLY cross-engine control) → the **T7 comparator** (join by object name, paired residual operands,
    83	  reject a native boolean, three-valued, repoint ablation) → reconcile → c1 step record → c2. ⛔ Mathematica:
    84	  2-seat licence, `--sandbox danger-full-access`, watch RSS (the full cross-engine residual may need ≥64 GB);
    85	  ⛔ use DETACHED launch (harness reaps run_in_background).
    86	
    87	🔨 **c1 WOLFRAM engine IN PROGRESS (2026-09-04):** WL build directive `c2b0f9dd` (orchestrator-written; **2
    88	decision legs** Codex+Grok, both EXIT=0, **5 findings folded one pass** — T-a re-derivation no longer
    89	answer-matches the spec §2a exposed normal; deleted the μ_θ "cancels-out" pre-registered agreement; μ_θ is a
    90	composite row + PY per-case locals so WL builds its OWN opaque `muTheta` with the residual adjudicated after the
    91	run; copy-test rescoped to run-time import-freedom + `Get`/`Import`/abs-path scan; shared-symbol table → exact
    92	reserved names in verified `mechanical_lower_camel` spelling `WBg`/`w1Profile`/`etaBg`/`LambdaA0` (capital L)/
    93	`rhoBrBgRho4Constant`). Detached Mathematica build (~2h, 30 GB box, peak ~1.4 GB — LIGHT) → **baseline `e139bc61`**
    94	(1862 lines, 51 tags, blindness byte-identical isolated-vs-in-repo). **2 build legs SERIALIZED** (fresh Claude +
    95	Grok; both ablate kernels, 2-seat): legs DISAGREED (Claude SOUND/0-MF; Grok 3 MUST-FIX) — I confirmed all 3 in
    96	the engine code (rule 13): (1) `DTN_OPERATOR` composition freezes the input leg (both `N_0` at `momentumOutput`,
    97	inconsistent with the verified two-leg `DTN_KERNEL`); (2) energy BULK operand = face quantity × free
    98	`farFieldPhase`, not a genuine Poynting flux from φ at `|w|→∞` (routes not independent — SAME class as SymPy R1);
    99	(3) energy FACE operand = reconstructed impermeable traction, not the response `t_s` (no `Λ_X`) — can't catch a
   100	response `t_s` error. Core physics 2-leg-SOUND (two-momentum kernel both legs, rigid-shift, operator inverse, T-a
   101	re-derivation, §5 controls bite, blindness, names, μ_θ). ✅ **REPAIR CLEARED + committed `13f0bd2c`** (directive
   102	`S11c_c1_wl_repair_directive.md`: R1 composition momentum labels → both legs; R2 real far-field Poynting from φ
   103	at |w|→∞, `farFieldPhase` removed; R3 energy face operand binds the response `t_s`; R4 four NITs). +443/−103
   104	(`.wl`-scoped; the full commit `e139bc61..13f0bd2c` is 4 files +691/−103) vs `e139bc61`, core BYTE-IDENTICAL
   105	(`directBoundaryDerivation` sha `557ed758`=baseline). **2 fresh re-review legs
   106	(Claude+Grok, serialized), both "no MUST-FIX"** — each control confirmed BITING by independent one-sided
   107	corruption (refreeze rightmost → nonzero; bulk-only q flip moves only bulk, no structural-zero collapse; response
   108	t_s flip moves only face). Records: `_measurements/S11c_c1_wl_build_directive.md`,
   109	`_measurements/S11c_c1_wl_build_review.md`, `_measurements/S11c_c1_wl_repair_directive.md`.
   110	
   111	⛔⛔ **BUT THE REPAIR DIRECTIVE SKIPPED ITS RULE-7 DECISION LEGS (2026-09-04) — a real gap, user-flagged.** I
   112	authored `S11c_c1_wl_repair_directive.md` (orchestrator-written, physics-bearing) and built against it WITHOUT its
   113	2 decision legs, rationalizing via "SymPy precedent" + "repairs are lighter" — both FALSE
   114	([[feedback_directive_design_review]] now closes these by name). Run RETROACTIVELY (Codex+Grok, both in) the gate
   115	found the directive **not sound** — and the two legs caught DIFFERENT defects (one leg alone would have shipped it):
   116	(1) R1 invariant LEAKS residual=0 + re-enters at `DTN_KERNEL` (rule 5; Grok; code OK — genuinely unfroze);
   117	(2) `NONINVERTIBILITY_CONDITION` **over-protected → the input-leg freeze PERSISTS** in `fredholmFunctionSpaceOperator`
   118	(`.wl:580-597`, both `gZero` on `momentumOutput`, used at `.wl:1537-1549`; PROPAGATED, re-review structurally blind; emit-only, NOT exported) — ✅ FIXED by repair-2 R1;
   119	(3) ⚠ R4b "parity SWAP in `PERMEABLE_PORT_HERMITIAN`" — **REFUTED + RELOCATED by the repair-2 decision gate** (Codex+Grok + rule-13): `PERMEABLE_PORT_HERMITIAN` is CORRECT (the congruence `Aᵀ·diag(P₊,P₋)·A` under `A=faceToParityMatrix`, both diagonal blocks EVEN, coupling odd, correct given `V_s=s∂_tζ_s`; the naive `δW↔(ζ₊−ζ₋)` is a category error that vanishes the thickness port at equal faces). The REAL propagated parity defect is the **dead parity axis** in `PERMEABLE_DISSIPATION_VS_OMEGA_TAU` (per-face form under both parity keys) — ✅ FIXED by repair-2 R2 (both re-review legs independently re-derived the congruence + confirmed the fix bites);
   120	(4) R2 "independent data" + no φ-driver (Codex MUST-FIX / Grok NIT; code OK — bound `V_s`); (5) R1 Ŵ-over-spec,
   121	probe-scope, blanket-"no-residual" disclaimer inconsistency (NITs, not inherited). Arc-verify (Codex) also: rule-2
   122	audit trail incomplete (records point to ephemeral `/tmp` logs, not literal output; "+443/−103 one file" is
   123	engine-scoped not the 4-file commit) + **neither c1 engine's canonical `.out` committed** (WL nor SymPy) — so
   124	"per-engine DONE" is real but "c1 DONE" overstates.
   125	
   126	✅ **REMEDIATION COMPLETE (2026-09-04; plan `_measurements/S11c_c1_wl_remediation_plan.md`).** (A) **repair-2**
   127	`dd34d564` — its 2 decision legs ran BEFORE the build (this time), and **caught that the first-draft R2 targeted a
   128	NON-defect** (would have corrupted the correct `PERMEABLE_PORT_HERMITIAN`); folded once → build (+144/−45, 51 tags,
   129	blindness intact) → 2 re-review legs (fresh Claude Agent + Grok, serialized) BOTH CLEAR (R1 Fredholm two-leg
   130	re-freeze moves the probe 2→0; R2 dissipation parity-combination one-sided +face corruption moves DELTA_W/ZETA_C
   131	differently; protected core byte-identical; PORT_HERMITIAN congruence independently re-derived by both). (B) rule-2
   132	records corrected + the retroactive repair-directive decision-leg record added. (C) both `.out` committed
   133	`4a14100a` via GIN (WL SHA-256 `d63e616786…` byte-identical to the reviewed engine; SymPy run reproduced the
   134	committed export). (D) "DONE" overstatement corrected here + in the repair record. ⚠ **2 benign NITs → the c1 step
   135	record:** `LOCAL_TAG_NAMES` bookkeeping churn; the two diagonal parity Hermitian blocks are PROPORTIONAL (ratio 4,
   136	both EVEN) so the real independent parity content is the off-diagonal coupling (the R2 control still distinguishes
   137	repaired from the dead-axis bug). Core physics 2-leg-SOUND + exports UNaffected throughout (the propagated defects
   138	were emit-only; c2 consumes DtN operator/kernel/face response, all sound).
   139	
   140	✅ **c1 T7 CROSS-ENGINE COMPARATOR — BUILT + re-reviewed SOUND (2026-09-05).** `scripts/S11c_c1_cross_engine_comparator.py`
   141	+ `test_S11c_c1_cross_engine_comparator.py`. Frozen contract `S11c_c1_SHARED_PHYSICS.md:580-587` (N8, measurement-only:
   142	computes/PRINTS, decides NOTHING; three-valued residual OBJECTS printed — NOT a verdict token; join by object name on
   143	S11c-a axis-typed keys; INJECTIVE `mechanical_lower_camel`; reject native boolean; repoint ablation; ⛔ NO pre-registered fold).
   144	Reuses the S11c-a/S11c-b mechanical base. **50 shared non-LOCAL families (identical both engines)**, 1080 joins across 27.
   145	- **Directive GATED (84686a54)**: orchestrator-written → 2 decision legs (Codex+Grok) caught **12 fold areas** before the build
   146	  (the self-contradicting `qOut` spelling entry that would have UNSEALED the two-momentum seal; FACE {1,-1} vs DIRECTION {1,2,3}
   147	  Integer mistyping; wrong per-family VALUE shapes; `raw_control_case` over-reach; `join>0` manufacture-joins pressure; the
   148	  suppressed three-valued residual) + my 2 rule-13 catches (cS0 non-inverse [WRONG, see repair]; the rule-17 bg-density field seal).
   149	- **Codex build (493k tokens) → reviewed baseline `7141e6ad`**: 30/30 + 48/48 tests; peak RSS 436 MB; full symbolic residual DEFERRED.
   150	- **2 re-review legs (fresh Claude Agent + Grok) BOTH CLEAN = SOUND**: no false agreement, no hidden disagreement; every SEAL
   151	  load-bearing (proven by name-map ablations — two-momentum `qOut[...]` stays applied, μ_θ opaque, ω-assumption, bg-density field,
   152	  regime/parity all UNRECONCILED); raw whitelist exactly `{DTN_OPERATOR, NONINVERTIBILITY OPERATOR}`; no silent 0-extract.
   153	- ✅ **SCOPED REPAIR applied + re-reviewed CLEAN (`704308af`)** — Codex `gpt-6-astra` high (first astra code job, 34k tokens, +14 comparator lines; both re-review legs Claude Agent+Grok CLEAN, every claim ablated): **R1** add
   154	  `cS0←c_s0` (⚠ I mis-computed `mechanical_lower_camel("c_s0")` as `cs0`; it is **`cS0`** — a legitimate bare-symbol fold the build
   155	  wrongly excluded; latent now, MATTERS for the deferred PERMEABLE run); **R2** held-parse the WL **4-arg triple-range**
   156	  `Inactive[Integrate]` + nested `Inactive[Limit][…,outwardDistance->Infinity]` across 3 ENERGY families (display/completeness,
   157	  c1-preprocessor-only, inherited 1-range BOUND_INTEGRAL untouched). Repair directive got its own 2 decision legs (sol+Grok) —
   158	  rule 7 honored (decision legs caught my R2 mis-description + that R1 broke a test).
   159	
   160	✅ **CROSS-ENGINE RECONCILE DONE + c1 STEP RECORD committed `072d0b75`** (`steps/S11c_c1_curved_bulk_closure.md`;
   161	record `_measurements/S11c_c1_comparator_reconcile.md`; reproducible `_measurements/S11c_c1_reconcile_reproduce.py`).
   162	Ran comparator family-scoped (46/50 families, 30 GB, peak ~317 MB). Method = the **staged representational bridge**
   163	[[feedback_reconcile_representational_bridge]]: comparator PRINTS raw (nonzero≠disagree), I apply JUSTIFIED
   164	identities (opaque q_out↔live qOut; FT-of-derivative jet identity; on-shell dispersion; σ_W binding; ε-placement)
   165	and test collapse. ⭐ **EARNED cross-engine AGREE:** (3) two-momentum **DtN KERNEL** — all 4 (anchoring×face) cases
   166	EXACT ZERO off-diagonal (Stage-2 residual collapsed to the dispersion relation itself; adversarial corruptions
   167	bite); (1) μ_θ + the **δp_s (pressure) + J_s (relative-flux) response COEFFICIENTS** (`ε·A≡B` exact zero at physical
   168	kinematics — sweep AGREE=54, `_measurements/S11c_c1_reconcile_coeff_sweep.py`); (2) ω; (4) regime/parity naming
   169	(keying convention). **2 step-record legs (Grok + Codex sol) THEN a correction-verify leg (Codex sol) scoped my
   170	verdict TWICE** — folded (`072d0b75`→`b8a83206`; leg reports `_measurements/S11c_c1_step_record_review_*.txt` +
   171	`_measurements/S11c_c1_reconcile_correction_verify_codex.txt`). ⚠ **My "clean anchoring split" was a TEST ARTIFACT**
   172	(FLAT leaves are diagonal half-space objects — test on-diagonal; the flat-resolvent leg-label q_out/q_in is equal
   173	only at k=k'). [[feedback_matching_number_is_not_evidence]]
   174	⛔⛔ **UNDECIDED — c2 must NOT treat as cross-engine-closed:** (5) **bg-density = a SURFACED rule-17 freeze** (PY
   175	bare constant vs WL live field `(ρbr/W₀)WBg(x)`; functional dependence agrees + 0 derivatives, but a scalar ≠ a
   176	field globally — folding it is the rule-17 collapse to avoid; **c2 re-adjudication MANDATORY** — c2's fold sums over
   177	the face); the **t_s (traction) response leaf** (WL zero-padded 4-vec `(0,0,0,scalar)` vs PY scalar — scalar-vs-vector
   178	rep, 6 TRACTION.FLAT nonzero); the raw **`dtn_operator`** whole-form (kernel-AGREE does NOT extend to it); **ENERGY**
   179	(PY closed-form vs WL far-field integral); the off-diagonal **flat-resolvent leg-labeling**. ⚠ **4 giants + full
   180	symbolic residual UNMEASURED — DEFERRED ≥64 GB** (⛔ NOT pre-adjudicated;
   181	`DEFERRED_HEAVY_RUNS.md`). NO cross-engine DISAGREEMENT found. **The kernel is cross-engine closed; the rest is
   182	UNDECIDED/deferred.** NEXT = **c2**. [[project_dual_engine_symbol_transliteration]] [[feedback_never_freeze_a_varying_field]]
   183	
   184	**⛔ ALL of c1 is DONE (both engines + comparator + reconcile + step record `072d0b75`→`b8a83206`; do NOT rebuild/re-reconcile). c2 SPEC v2 GATED + COMMITTED `16849fc6`** (`directives/S11c_c2_SHARED_PHYSICS.md`; gate record `_measurements/S11c_c2_spec_review.md`; leg reports `_measurements/S11c_c2_spec_review_{grok,codex_sol}.txt`). ⭐ **The gate PAID OFF — v1 had a wrong CENTRAL design claim** (the substitution-increment does NOT isolate c2 from S11c-b's two unvalidated SIGN CONVENTIONS — face-generalized-force PY `+diff`/WL `−linearVirtualVariation` + #90 closure-fold — they're the COEFFICIENTS of the substituted δp_s slots ⇒ leak `2·carrier·incr`; both legs ⚙-computed it). v2 fixes: increment drops the bulk/kinetic base ONLY; SURFACE the 2 conventions + adjudicate the face-force sign via the NEW **traction-vs-slab mechanical-power pairing** (c1-assigned to c2). Other folds: real θ-row has `J_s` ALREADY eliminated (#90) ⇒ substitute **closed δp_s + w-jets** into `delta_p_±`/`d_w_delta_p_±` slots, ⛔ NOT a closed `J_s` (double-count); `Λ_X` is **traction-only** (not in `closure_shape_deriv`); emit the **assembled two-face operator per (α,ρ)** (slab sums `J₊+J₋`); **3→6 re-adjudications** (added traction-slab pairing, flat-resolvent leg-labeling, μ_R,bg form control); 3 named operators (closed-coupling-kernel/ordering-commutator/self-energy-increment); positional `load_model` + `s11c_c1_*` write-keys (v1 bound bare S11b `face_response`); `Z→0 ≠ Λ_A=0 ≠ impermeable`; density field-vs-const via `background_density_map` not a `∇ρ` jet; N6 2-route = the two anchorings. ⚠ RULE-15 ARMED: a 2nd heavy round at the build-directive gate ⇒ hand re-author to Codex. **NEXT = c2 SymPy BUILD DIRECTIVE** (own 2 decision legs BEFORE the astra build; nail the IMPORT_KEYS root set + fold symbol map vs the real 2441+44-row files) → astra build → 2 re-review legs → WL engine → T7 → reconcile → c2 step record. **c2 carries the held folds:** extract-then-close ordering (close the full
   185	operator FIRST, then re-extract — extract/eliminate DON'T commute; counterexample `R_x=x+p,p=αy` → 0 vs α), the
   186	θ-row `Λ_X`(traction)/`J_s`(mass) row-routing, and the **substitution-increment** emit (c2's self-energy =
   187	closed − open-symbolic, so S11c-b's still-deferred cross-engine residual cancels out of **c2's** residual — c1
   188	does no fold, so c1's residual never involves the S11c-b slab operator).
   189	
   190	⚠ **MEMORY/OPS:** the c1 SymPy engine MEASURED LIGHT (1.4 GB peak) — the old "S11c-c heavier than S11c-b" fear
   191	was wrong for it; but the WL/Mathematica engine and the FULL cross-engine T7 residual may still be heavy (S11c-b's
   192	cross-engine residual needed ≥64 GB). Serialize the 2 engines, watch RSS, defer
   193	heavy controls in-band→out-of-band (`DEFERRED_HEAVY_RUNS.md`); the full cross-engine residual may again need the
   194	≥64 GB box. ⛔ NEVER 2 memory-heavy CAS jobs concurrently. Imports S11c-b's `af560257` coupling (per-engine
   195	leg-verified, cross-engine residual DEFERRED — [[project_s11c_b_state]]). NO per-substep card (one S11c roll-up
   196	after S11c-e). [[project_s11c_b_state]] [[project_s11c_a_state]] [[feedback_never_freeze_a_varying_field]]
   197	[[feedback_directive_design_review]] [[feedback_handcode_comparison_never_blanket_collapse]]
   198	
   199	---
   200	
   201	✅ **c2 STEP A (physics) + STEP B/C (export repair) DONE (2026-09-06)** — full detail in the frontmatter
   202	description; commits `8f3a017f` (physics review, fold SOUND) → `a5f7a06c` (export-repair directive gated) →
   203	`aa76105a` (export repair committed, 60 MB → 21.4 MB). Gate/verify records under `_measurements/S11c_c2_*`.
   204	- **The astra export-repair build worked** — `publication_compact` + a strict-recursive Integral-aware semantic
   205	  guard; my own adversarial verify (`_measurements/S11c_c2_export_repair_adjudication_verify.py/.out`) confirmed
   206	  faithful 5/5 + guard bites 5/5. No overstep (the fixed §7 fence held). Only cost: ~50 min / 3 full runs (slow).
   207	- **Two review-leg lessons this session:** (1) the fresh-Claude re-review leg produced FALSE-POSITIVE semantic
   208	  FAILs by comparing a **re-run** emitted (fresh Integral bound-dummies) to the stored compact (astra's dummies) —
   209	  Integral-bearing operators MUST be compared **same-process**; Grok and I did, and got clean. (2) **Grok's
   210	  recurring failure mode**: it flags a nonzero-looking object as a defect **without evaluating it to zero** (F's
   211	  `Integral(0,…)`; the c1 "tangential freeze"). Both are why rule-13 self-verification is non-negotiable.
   212	- **Ops:** every `run_in_background` task was reaped this session (~30 s, even trivial `until`-loops; NOT OOM,
   213	  23 GB free). The astra build survived only via `setsid`-detached launch (own session) + a DONE-marker +
   214	  re-armed watcher polling. [[feedback_background_tasks_can_die_spuriously]]
   215	- **NEXT (unchanged):** light §5e/§3c spec-wording fix ("genuine coupling decouples", not "increment vanishes";
   216	  review-until-clear) → blind Wolfram engine → 2 legs → c2 T7 comparator + reconcile → c2 step record.
43:- [✅ S11c-c: c1 DONE; c2 physics SOUND + export 60→21 MB; NEXT=spec fix → blind WL engine](project_s11c_c_state.md) — c2 fold physics 2-leg SOUND `8f3a017f` (Grok's F/G/E all FALSE POSITIVES, my rule-13 compute); export repair `aa76105a` (increment→EMIT-only, transparent-factor, guards bite 5/5). ⭐ astra=KEEP for code. ⚠ bg-tasks REAPED (setsid+marker); Integral ops need SAME-PROCESS compare. Full arc + carry-forwards + NEXT (§5e/§3c wording → WL engine → T7 → step record) in topic file. ⛔ c1 STANDS
46:- [✅✅ S11c-b CLOSED (per-engine) — step record `bcb9f7d7`; cross-engine residual DEFERRED to ≥64 GB; NEXT=S11c-c STARTED [[project_s11c_c_state]]](project_s11c_b_state.md) — variable-coeff slab operator + off-diagonal kernel; basis 40, un-frozen op, pin B constraint-fold, #90 face+response, all per-engine 2-leg CLEAR; exports `af560257`. Full arc + deferrals + the 2 sign conventions in the topic file. ⛔ do NOT re-litigate pin B/P2a/P2b

exec
/bin/bash -lc "rg -n 'gpt-6-astra|gpt-5\\.6-sol|MODEL POLICY|CODE→|DOC/decision|2-seat|setsid-detached|DONE-marker|re-armed watcher|run_in_background' CLAUDE.md .claude research/pde_ledger_v3/directives research/pde_ledger_v3/_measurements /home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md 2>/dev/null | sed -n '1,320p' && nl -ba research/pde_ledger_v3/S9_export_chain_rebuild_directive.md | sed -n '1,30p' 2>/dev/null || true && nl -ba research/pde_ledger_v3/directives/S9_export_chain_rebuild_directive.md | sed -n '1,35p' 2>/dev/null || true && nl -ba research/pde_ledger_v3/directives/S11c_decisions.md | sed -n '35,95p'" in /var/projects/toy_physics
 succeeded in 0ms:
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:63:  workaround** (`setsid`+marker+Monitor) — the harness reaped 3 `run_in_background` tasks at ~87s; NOT OOM (1.4 GB
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:84:  2-seat licence, `--sandbox danger-full-access`, watch RSS (the full cross-engine residual may need ≥64 GB);
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:85:  ⛔ use DETACHED launch (harness reaps run_in_background).
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:95:Grok; both ablate kernels, 2-seat): legs DISAGREED (Claude SOUND/0-MF; Grok 3 MUST-FIX) — I confirmed all 3 in
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:153:- ✅ **SCOPED REPAIR applied + re-reviewed CLEAN (`704308af`)** — Codex `gpt-6-astra` high (first astra code job, 34k tokens, +14 comparator lines; both re-review legs Claude Agent+Grok CLEAN, every claim ablated): **R1** add
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:212:- **Ops:** every `run_in_background` task was reaped this session (~30 s, even trivial `until`-loops; NOT OOM,
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:213:  23 GB free). The astra build survived only via `setsid`-detached launch (own session) + a DONE-marker +
/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:214:  re-armed watcher polling. [[feedback_background_tasks_can_die_spuriously]]
.claude/skills/build/SKILL.md:25:3. Launch Codex at xhigh effort using the **Bash tool with `run_in_background: true`**. ⛔ Do NOT use a
.claude/skills/review-legs/SKILL.md:228:1. Start Grok using the **Bash tool with `run_in_background: true`**, and ⛔ **no shell `timeout`**.
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_decision_gate.md:4:(orchestrator-written) → **two decision legs = Codex `gpt-5.6-sol` xhigh + Grok** (G1), one two-leg pass, fold
research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md:4:astra/`gpt-6-astra` build). Emitted output `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out` (499 MB), navigated by
research/pde_ledger_v3/directives/S11c_c1_wl_build_directive.md:274:  wrapper). The harness reaps `run_in_background` Bash tasks at ~87 s; a naive `run_in_background` launch of a
research/pde_ledger_v3/directives/S11b_wl_engine_repair_directive.md:67:One kernel at a time (2-seat licence); orchestrator writes the committed `.out` after review; demonstration
research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:6:symbol map, the emit tags, and the three script clauses. Model: **CODE build → `gpt-6-astra` high**.
research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:238:Detached launch (`setsid` + a completion marker + `Monitor`; the harness reaps `run_in_background` ~87 s). c2's fold
research/pde_ledger_v3/directives/S11b_wl_engine_fix_round2_directive.md:41:One kernel at a time (2-seat licence); orchestrator writes the committed `.out` after review; demonstration
research/pde_ledger_v3/directives/S11c_b_p1_wl_residual_emit_directive.md:74:    MemAvailable floor of 2.5 GB, ONE kernel at a time (2-seat licence), and a concrete `timeout` set ABOVE the
research/pde_ledger_v3/directives/S11c_b_p1_wl_residual_emit_directive.md:99:`--sandbox danger-full-access`, xhigh) ⇒ 2 build legs (fresh Claude agent + Grok; ⛔ SERIALIZE — Mathematica 2-seat
research/pde_ledger_v3/directives/S11c_c2_export_repair_directive.md:3:You are the builder (`gpt-6-astra`). This is a **publication-only** repair of one already-reviewed script's
research/pde_ledger_v3/directives/S11c_a_comparator_reemit_plan.md:203:background via `run_in_background`, ⛔ never `&`; logs OUTSIDE repo, prompts UNDER project; rule-2 twin +
research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:20:(Codex `gpt-5.6-sol` xhigh + Grok) before any builder**; the build directive that follows gets its own two decision
research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:354:  (harness reaps `run_in_background` ~87 s). Mathematica: 2-seat licence, `--sandbox danger-full-access`.
research/pde_ledger_v3/directives/_measurements/S11c_c1_spec_retro_review_adjudication.md:11:codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_repair_directive.md:35:Orchestrator-written repair directive → legs = Codex(`gpt-5.6-sol` xhigh, doc-review policy) + Grok (identical
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_repair_directive.md:39:$ codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh "$(cat _legs/S11c_c1_comparator_repair_directive_review.md)" > <scratch>/codex_repair_dir_review.txt 2>&1   # detached; -m flag accepted (v0.153.4)
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_repair_directive.md:43:RESULTS (both EXIT=0; converged; `-m gpt-5.6-sol` accepted on Codex v0.153.4). Both legs: **R1 CLEAN**
research/pde_ledger_v3/directives/_measurements/S11c_c1_comparator_repair_directive.md:66:**Repair build (Codex `gpt-6-astra` high — first astra code job; new model policy):** `codex exec -m gpt-6-astra
research/pde_ledger_v3/directives/_measurements/S11c_c1_wl_build_directive.md:9:## Launch (detached — escapes the ~87 s `run_in_background` reap)
research/pde_ledger_v3/directives/_measurements/S11c_b_wl_build_review.md:5:**Legs (Codex-written deliverable → fresh Claude agent + Grok), serialized (Mathematica 2-seat), prompt:**
research/pde_ledger_v3/directives/_measurements/S11c_b_wl_repair_review.md:6:serialized (Mathematica 2-seat), prompt `_legs/S11c_b_wl_repair_review.md`. **Raw:** Grok
research/pde_ledger_v3/directives/_measurements/S11c_c2_spec_review.md:14:codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
research/pde_ledger_v3/directives/_measurements/S11c_c2_build_directive_decision_review.md:11:codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
research/pde_ledger_v3/directives/_measurements/S11c_c1_wl_repair_directive.md:33:# Re-review leg 2 — Grok, serialized SECOND (both ablate Mathematica; 2-seat)
research/pde_ledger_v3/directives/_measurements/S11c_c1_wl_build_review.md:5:(2-seat licence): fresh-Claude first (in-process Agent), then Grok (detached). Leg prompt:
research/pde_ledger_v3/directives/_measurements/S11c_c1_wl_build_review.md:15:# Leg 2 — Grok, serialized SECOND (both ablate Mathematica kernels; 2-seat licence)
research/pde_ledger_v3/directives/_measurements/S11c_c1_wl_repair2_directive.md:8:Mathematica; 2-seat licence). Re-review prompt: `directives/_legs/S11c_c1_wl_repair2_review.md`.
research/pde_ledger_v3/directives/_measurements/S11c_c1_record_corrections_review_adjudication.md:19:codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
research/pde_ledger_v3/directives/_measurements/S11c_b_wl_admissibility_build_review.md:5:(serialized for the 2-seat licence), and derived the normalization independently in SymPy. Both verdict:
research/pde_ledger_v3/directives/_measurements/claude_md_restructure_review.md:5:**Codex `gpt-5.6-sol` xhigh + Grok `grok-4.6`** reviewed (whatever-writes-doesn't-review: Astra≠reviewer).
research/pde_ledger_v3/directives/_measurements/claude_md_restructure_review.md:12:codex exec -m gpt-6-astra -c model_reasoning_effort=high --sandbox danger-full-access \
research/pde_ledger_v3/directives/_measurements/claude_md_restructure_review.md:15:codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
research/pde_ledger_v3/directives/_legs/S11b_wl_fix2_review_prompt.md:34:raise it; ⛔ one kernel at a time (2-seat licence); ⛔ kill on RSS>6 GB, `free -h` first if you do. Save
research/pde_ledger_v3/directives/_legs/S11c_c1_wl_build_directive_review.md:83:   the §6 locus protocol only for finite-dim algebraic loci)? Is the run discipline right (one kernel / 2-seat;
research/pde_ledger_v3/directives/_legs/S11c_a_conormal_adjudication_leg.md:54:- wolframscript (only if you use it): wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat
research/pde_ledger_v3/directives/_legs/S11c_a_wl_build_review_prompt.md:49:   symbolically; one kernel at a time (2-seat licence); kill criteria (600 s no-output / 6 GB) by PID;
research/pde_ledger_v3/directives/_legs/S11c_c1_wl_repair_review.md:59:⛔ NEVER raise the timeout; run AT MOST ONE kernel at a time (2-seat licence — another leg may run after you).
research/pde_ledger_v3/directives/_legs/S11_wl_wall2_strata_diagnosis_prompt.md:57:  Python/SymPy (available on this box). This is a 2-seat licence machine and memory is the
research/pde_ledger_v3/directives/_legs/S11c_b_jet_depth_reconciliation_consult.md:97:- ⛔ If you run a CAS kernel (SymPy is free; Mathematica has a 2-seat licence — assume the other advisor may hold
research/pde_ledger_v3/directives/_legs/S11c_a_projection_window_adjudication_leg.md:52:- wolframscript (only if you use it): wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat
research/pde_ledger_v3/directives/_legs/S11c_b_p1_wl_build_review.md:45:⛔ NEVER raise the timeout; ⛔ run ONE kernel at a time (2-seat licence). Kill an orphaned kernel by exact pid; if a
research/pde_ledger_v3/directives/_legs/S11c_b_89b_wl_operator_build_review.md:105:- ⛔ This is a Mathematica artifact: wrap EVERY kernel run in `timeout 600`; run ONE kernel at a time (2-seat
research/pde_ledger_v3/directives/_legs/S11c_b_89b_wl_repair_rereview.md:57:- ⛔ Mathematica: wrap EVERY kernel run in `timeout 600`; ONE kernel at a time (2-seat licence); ablate `/tmp`
research/pde_ledger_v3/directives/_legs/S11c_a_wl_bgcurrent_build_review.md:45:run in `timeout 600`; a hit is a FAILED ablation — report and move on. ⛔ ONE kernel at a time (2-seat licence).
nl: research/pde_ledger_v3/S9_export_chain_rebuild_directive.md: No such file or directory
     1	# S9 — build the export chain, starting at the beginning
     2	
     3	**Revision 2.** Revision 1 went to two independent review legs (Codex, Grok) and **both said do not build
     4	it.** Ten findings; every one is applied below. Where the legs disagreed, the disagreement is adjudicated
     5	here and the measurement is given.
     6	
     7	**Read `research/pde_ledger_v3/S9_REWRITE_PLAN.md` first.** It is the governing artifact and it is CLOSED.
     8	
     9	Repository root: `/var/projects/toy_physics`. Ledger root: `research/pde_ledger_v3`.
    10	
    11	---
    12	
    13	## What this is for
    14	
    15	Each step's SymPy engine writes a module the next step imports, so cross-step consistency is true by
    16	construction rather than reconciled against a second copy. Each step's Wolfram engine imports nothing and
    17	re-derives independently — **that is the only blindness control in this design, and nothing else in it
    18	should be built pretending to be one.**
    19	
    20	S9 is where the chain begins: it is the first step of the light sector and consumes nothing from an earlier
    21	v3 step.
    22	
    23	A registry (`reduction/`) previously held declared quantities, declared dimensions and hand-written
    24	name→name pair tables. **It has been deleted — 63 files, commit `fb29bba2`.** Do not look for it, do not
    25	restore any part of it, and do not create anything shaped like it.
    26	
    27	---
    28	
    29	## The files you may create or change — anything else is a STOP
    30	
    31	| may create | may change |
    32	|---|---|
    33	| `scripts/S9_exports.py` (generated by the run) | `scripts/S9_light_requires_shear_sympy_audit.py` |
    34	| one knob-extractor script under `scripts/` | `mathematica/S9_light_requires_shear_mathematica_audit.wl` — **emitted name strings ONLY** |
    35	
    35	`scripts/ledger_fold.py`, §D3) — ⛔ **not** the accumulated whole-model file. S11c-a and S11c-b, already built
    36	under the prior accumulate topology, are the frozen base, byte-unchanged. The Wolfram engines are **blind** and
    37	re-derive. A single roll-up ledger **card** closes the family. `BUILD_INPUT_DIGESTS` pins, per sub-step,
    38	`{that sub-step's audit, the imported exports, that sub-step's spec}` — and, for a delta-writing sub-step,
    39	**also `scripts/ledger_fold.py`** (a shared executable input, §D3). ⛔ **Do not reuse the `S11b` slug.**
    40	
    41	**N2 · The split — decided WITH the legs; ⛔ not pre-committed, ⛔ not frozen. ⭐[folded]**
    42	Both legs independently found the draft's three-seam sketch **too coarse** — its "result" seam re-merged
    43	spectrum + leakage + falsification + conclusion (the largest surface, and the one that re-invites the
    44	forbidden global spectrum, `N5`). The **adopted starting family** (Codex's five, into which Grok's
    45	kinematics/matrix-element/photon three-way maps):
    46	
    47	| sub-step | scope — named objects, ⛔ not recipes | exports |
    48	|---|---|---|
    49	| **S11c-a** background & geometry | which quantities vary (`W₀(x)`, `μ_R(x)`, `ρ_4D⁰`); the profile's **physical anchoring** (material-advected vs externally-held, and the stationarity/force that holds it, `N4`); the **(ε,η) power counting** (`N12`); the **FULL tilted-face shape-derivative** of every interface law — normal, relative flux, traction, bulk-field-at-shifted-face, DtN — ⛔ not only normal face motion (`N3`) | kinematics + boundary shape-derivatives |
    50	| **S11c-b** variable-coefficient brane operator | the divergence-form slab operator + the **off-diagonal** transverse→`{θ,e_W,u_L}` kernel; inherit S11b invariants with variable coefficients and **emit** new gradient-of-background invariants as **results** (`N15`); the **representation-invariance control** (Eulerian ≡ material after field redefinition, `N4`/`N6`) | the coupling kernel |
    51	| **S11c-c** curved-interface bulk closure | the perturbed two-face outgoing bulk **DtN/impedance** (flux, traction, permeability/memory, face parity); the `v_bulk_normal_0` carry-or-restrict decision (`N11`) | coupled nonlocal self-energy/operator |
    52	| **S11c-d** profile-conditioned spectrum/scattering | for an **explicitly named profile class** (localized→Born/scattering kernel; periodic→Bloch; slowly-varying→WKB) — ⛔ **not** a generic "full spectrum"; the profile-conditioned **leakage rates** live here; ⛔ **no global dispersion relation** (`N5`) | scattering amplitudes / resonances / local spectrum |
    53	| **S11c-e** leakage observable & falsification | the flux-normalized **dimensionless conversion FORM** for a supplied profile (⭐ magnitude needs the interior `R1` — out of scope, `V3_STEP_PLAN.md:1179`); **confinement** interpreted here (`N13`); the bench-optics bound **withheld**, diffed orchestrator-side (`N7`) | the conversion observable + leakage |
    54	
    55	⚠ This is the **starting** structure; each sub-step's SPEC refines its own boundary. ⛔ Do not treat a–e as
    56	frozen before their specs. ⚠ Sequencing/scale (how many run now) is the user's call (rule 11).
    57	
    58	---
    59	
    60	## B · The five requirements — re-validated & CORRECTED ⭐[folded]
    61	
    62	**N3 · Tilted faces — the FULL first-order boundary shape-derivative, ⛔ not only normal face motion.**
    63	⭐[folded — both legs.] With `W₀(x)` varying, an in-plane displacement across `∇W₀` gives normal face
    64	motion `∼ u·∇W₀` at the coupling's own gradient order (correct as far as it goes). ⛔ But the same tilt
    65	changes, **at the same order**, every interface law: the outward normal `n_s = (−½∇W₀, s)+O(|∇W₀|²)`, hence
    66	the normal derivative `n_s·∇₄`, the **relative flux** (`S11b_SHARED_PHYSICS.md:195`, flux is measured along
    67	the face's outward normal), the traction/pressure work, the evaluation of bulk fields at the **shifted**
    68	face, and the perturbed **DtN/impedance** operator. Keeping only the normal-motion term can give the wrong
    69	coefficient or the wrong **face parity**. ⇒ S11c-a demands the full level-set/shape derivative of every
    70	interface law.
    71	
    72	**N4 · Eulerian vs material is a REPRESENTATION, ⛔ not a free physics choice; the physics is the profile's
    73	anchoring.** ⭐[folded — both legs, verified against `S11b_SHARED_PHYSICS.md:73,320-341`.] S11b already
    74	fixed the representation: `θ` is **Eulerian** (:73) and the binding constraint is **material**
    75	(`δ_v Σ_mat = 0`; `δ_vΣ_E = 0` is explicitly forbidden, :341). Eulerian and material variables are related
    76	by the map + Jacobian (`Δρ = δρ_E + u·∇ρ⁰`) and **must agree after that field redefinition** — treating the
    77	choice as "physics the spec picks" can double-count or omit `u·∇ρ⁰`. ⭐ The **real** physical decision is
    78	whether the inhomogeneous background profile is **advected with the material** or **held fixed in lab/
    79	Eulerian space** — S11c-a states it, and requires the two derivations to map to the **same** operator (the
    80	representation-invariance control, `N6`). ⚠ The load-bearing advective term is `u·∇Σ_E⁰ = u·∇(ρ_4D⁰W₀)`,
    81	tied to the freeze (`N11b`); the bare `u·∇ρ_4D` can drop while the coupling does not.
    82	
    83	**N5 · Plane waves are not eigenmodes — and the OUTPUT must be a profile-conditioned object with an
    84	order-count, ⛔ not a "full spectrum."** ⭐[folded — both legs.] The named trap (a global dispersion
    85	relation) is avoided, but `N10`'s "full variable-coefficient slab spectrum" **re-invites it**: there is no
    86	universal spectrum for an **unspecified** coefficient function. ⇒ S11c-d must **name a profile class** and
    87	fix the output as the matching object (localized→Born/scattering kernel; periodic→Bloch; slowly-varying→
    88	WKB), ⛔ never `ω(k)` for generic `W₀(x)`. ⚠ And it must fix the **expansion**: for inhomogeneity `O(η)`,
    89	the off-diagonal coupling is `O(η)` and the leading leakage rate `O(η²)` — so an engine cannot discard the
    90	leading coupling as "second order" nor confuse `O(η²)` leakage with the excluded nonlinear-intensity program
    91	(`N12`). ⛔ "Perturbation + off-diagonal coupling" is a recipe (rule 3); the **object** is the linear mixing
    92	between the uniform transverse and thickness sectors.
    93	
    94	**N6 · The uniform-limit control stays as a REGRESSION; the real control is the independent shape/coordinate
    95	route + one-sided corruption, named NOW.** ⭐[folded — both legs.] The uniform limit cannot validate the

exec
/bin/bash -lc "nl -ba research/pde_ledger_v3/directives/S11c_decisions.md | sed -n '90,125p' && nl -ba .claude/skills/build/SKILL.md | sed -n '1,80p' && nl -ba .claude/skills/review-legs/SKILL.md | sed -n '200,255p'" in /var/projects/toy_physics
 succeeded in 0ms:
    90	leading coupling as "second order" nor confuse `O(η²)` leakage with the excluded nonlinear-intensity program
    91	(`N12`). ⛔ "Perturbation + off-diagonal coupling" is a recipe (rule 3); the **object** is the linear mixing
    92	between the uniform transverse and thickness sectors.
    93	
    94	**N6 · The uniform-limit control stays as a REGRESSION; the real control is the independent shape/coordinate
    95	route + one-sided corruption, named NOW.** ⭐[folded — both legs.] The uniform limit cannot validate the
    96	gradient coefficient/sign/parity (S11b: coupling identically zero, `steps/S11b_interface_coupling_law.md:75`)
    97	— but it is a useful **smoke test** for a forbidden gradient-**independent** term; keep it as secondary.
    98	⛔ A "gradient-order control" phrased as "coupling `∝ ∇W₀` ⇒ vanishes at `∇W₀=0`" **is the vacuous uniform
    99	limit renamed**, and `∇W₀→0` is ⛔ **not** an accepted corruption. The genuine control (rule 14): derive the
   100	off-diagonal kernel by direct level-set/graph linearization, derive it **again** after flattening faces into
   101	material coordinates, transform Eulerian↔material exactly and compare, then **corrupt one route only** —
   102	flip one face's slope term, or omit its advective-density contribution — and require a nonzero residual.
   103	⚠ There are ≥2 same-order channels (tilt, `N3`; advection, `N4`); the one-sided corruption is the
   104	**independence test between them**, ⛔ so "the gradient channel" (singular) is wrong.
   105	
   106	**N7 · Falsification: emit the dimensionless conversion FORM; ⛔ the `O(1)`/grating reductio is withheld
   107	from the builder.** ⭐[folded — both legs; `O(1)` STRIPPED.] The lab bounds a **dimensionless conversion**
   108	(photons lost at an edge), ⛔ not the bare, dimensionful `∂²U/∂u_T∂e_W` — a nonzero coupling still carries an
   109	undetermined normalization until "what couples to what" is fixed (`S11b_SHARED_PHYSICS.md:824-825`). ⇒
   110	S11c-e emits the **flux-normalized dimensionless conversion FORM** for a supplied profile (⭐ the FORM is
   111	computable now; the **magnitude** needs the throat interior `R1`, out of scope — `V3_STEP_PLAN.md:1179`,
   112	which also fixes the structure: `∝ k·a`, mixing driven by `∇μ_R ≠ 0`). ⛔⛔ The `O(1)`-fraction / diffraction-
   113	grating argument is an order-of-magnitude **hint** (rule 5: the builder iterates to any target it can see) —
   114	it belongs in the **step record**, ⛔ never in the builder-facing acceptance text; the withheld numeric bound
   115	is diffed orchestrator-side. ⚠ A slit edge is an order-unity localized gradient, ⛔ not the small-`∇W₀`
   116	regime of a Born kernel — the spec must state which, so we do not compare a weak-gradient coefficient to a
   117	non-perturbative lab situation.
   118	
   119	---
   120	
   121	## C · Requirements the five OMITTED — added ⭐[folded]
   122	
   123	**N12 · Two-parameter power counting, on every result.** ⭐[folded — Codex's sixth requirement, Grok's
   124	order-bookkeeping.] Every S11c object carries an order label in **wave amplitude `ε`** and **background
   125	inhomogeneity `η`**: the transverse↔thickness coupling is `O(εη)`; the leakage probability/rate `O(ε²η²)`.
     1	---
     2	name: build
     3	description: Launch one Codex build and its two independent review legs as a single operation for the PDE ledger. Leak-gates the directive before launch, backgrounds an xhigh Codex run with absolute paths and no shell timeout, verifies the deliverable rather than the session, then automatically starts fresh-agent and Grok review before the caller reads the result.
     4	allowed-tools: Bash, Read, Edit, Write, Agent
     5	user_invocable: true
     6	---
     7	
     8	# Build With Automatic Review
     9	
    10	Invoke as `/build <absolute-directive> <absolute-deliverable> --check "<physics to check>"`.
    11	⛔ **There is no `--do-not-read` argument.** ⚠ It was one until 2026-08-12 — a denylist means the
    12	architecture is wrong, and this file's own CUT table at §143 had already said so.
    13	
    14	This is one operation: build, verify the deliverable, and launch both review legs. Never return a
    15	successful build without running `.claude/skills/review-legs/SKILL.md` yourself.
    16	
    17	## Runbook
    18	
    19	1. Resolve the directive, deliverable, log, review prompt, and review outputs to absolute paths before
    20	   launch. Background shells retain their cwd between calls; a relative write followed by an absolute
    21	   read has already produced a silent empty prompt.
    22	2. Leak-gate the directive by grepping it with fixed-string `rg` searches for every pre-registered answer it must not
    23	   contain. If anything matches, repair and re-run the gate **before launch**. Codex snapshots the prompt
    24	   into argv, so editing the file after launch changes nothing.
    25	3. Launch Codex at xhigh effort using the **Bash tool with `run_in_background: true`**. ⛔ Do NOT use a
    26	   shell `&` — this harness detaches the job itself and notifies you when it exits; a `&` inside a
    27	   foreground call leaves the build untracked and unreported.
    28	
    29	   ```bash
    30	   codex exec -c model_reasoning_effort=xhigh "$(</absolute/directive.md)" \
    31	     > /absolute/path/OUTSIDE/the/repo/codex-build.log 2>&1
    32	   ```
    33	
    34	   ⭐ **Write the raw transcript OUTSIDE the repository — ⚠ as TREE HYGIENE, ⛔ not as a blindness claim.**
    35	   ⚠ A transcript carries the engine's complete tag values verbatim and is noise in the tree.
    36	   ⛔⛔ **It is NOT a leak to be plugged by relocation:** the same measurement (2026-08-03, several such
    37	   files sitting in `_scratch/`, reachable by none of the naming conventions) is evidence that ⭐ **hiding
    38	   cannot work**, ⛔ not that it should be done harder ⇒ the CUT table at §137 below, and `CLAUDE.md`
    39	   rule 12.
    40	
    41	   Add `--sandbox danger-full-access` when the build must run Mathematica. ⛔ Never wrap the command in a
    42	   shell `timeout` — SIGKILL has cost 300k+ tokens.
    43	4. Do not poll for completion; the harness re-invokes you when the job exits. On that notification,
    44	   verify the **deliverable**, not exit status: require the
    45	   requested artifact to exist and be non-empty, and require the transcript/token usage to be plausible
    46	   for the requested build. Exit 0 plus `hook: Stop` has accompanied an empty prompt and no work; the
    47	   measured tell was about 3k tokens instead of 37k+.
    48	5. Do not open the deliverable or read its results. Immediately read and execute
    49	   `.claude/skills/review-legs/SKILL.md` with the deliverable and `--check` from this invocation.
    50	   ⛔ **There is no `--do-not-read` argument to pass** — it was cut 2026-08-12 (rule 12).
    51	   Launching the legs is the build skill's responsibility; do not tell the caller
    52	   to invoke `/review-legs` later.
    53	   ⭐ The deliverable here is **Codex-written**, so its two legs are **a fresh Claude agent + Grok**.
    54	   ⛔ Codex does not review what Codex wrote — see that skill's authorship table.
    55	6. Return only after both independent review legs have completed. Keep their findings separately
    56	   attributed so the caller can filter them before acting.
    57	
    58	## ⭐⭐⭐ EVERY SCRIPT DIRECTIVE CARRIES THESE THREE CLAUSES — non-negotiable
    59	
    60	> **1. The script may PRINT computed objects. It may NOT state conclusions.** An `emit`/`Print` payload
    61	> must be a CAS object — an expression, a solved root, a boolean from a symbolic test. ⛔ Never prose
    62	> describing a result.
    63	> **2. PRINT the residual; do NOT assert it.** `assert residual == 0` **is the builder writing down the
    64	> expected output**, and it turns an informative value into a binary crash. Compute → emit → *then* assert.
    65	> **3. Interpretation belongs to the STEP RECORD.** ⛔ The script does not editorialise.
    66	
    67	⚠⚠ **Measured 2026-08-04: named tags at named lines in three independently-built steps are typed sentences
    68	with no CAS object behind them** (e.g. `S11bB:348-360`, and a whole transverse block whose symbol `mu_R`
    69	appears in **no expression in the file**). ⛔ **Do not quote a fraction of the corpus** — two review legs
    70	rejected that as unmeasurable. ⛔ Cross-engine agreement on such a tag is **vacuous** —
    71	both engines carry the same author's sentence — and **eight fidelity review legs did not catch it**, because
    72	*"does it say the right thing"* and *"does it depend on anything"* are different questions.
    73	⇒ `research/pde_ledger_v3/REBUILD_HANDOFF.md` · `[[feedback-scripts-print-never-assert]]`.
    74	
    75	⭐ **Copy the one shape that worked** (`S11bB` lines 421–443): type a candidate, compute the object
    76	**independently**, `emit` the symbolic difference, hard-stop if nonzero.
    77	
    78	## ⭐⭐⭐ FOUR COROLLARIES — all measured 2026-08-04, all of them defeat the three clauses as written
    79	
    80	**1. ⛔⛔ A HAND-TYPED CAS OBJECT IS STILL HAND-TYPED.** Clause 1 bans a *prose* payload. It does **not**
   200	⭐ **Put these in the prompt file for any artifact that spawns a CAS kernel:**
   201	```
   202	⛔ Wrap EVERY kernel run in `timeout 600`. A 600s hit is a FAILED ablation — report it and move on.
   203	⛔ NEVER raise the timeout, and ⛔ never run more than one kernel at a time (the licence has TWO seats).
   204	⛔ Copy the artifact to /tmp and ablate the COPY. ⛔ Never modify the working tree.
   205	⭐ Save every ablation script AND its literal stdout to named absolute paths, and report those paths.
   206	```
   207	⛔⛔ **AND AN ORPHANED KERNEL LEAKS MEMORY, NOT JUST A SEAT — the symptom looks unrelated to Mathematica.**
   208	⚠ **Measured 2026-08-04:** a kernel left over from an over-budget leg reached **17.9 GB RSS at 100% CPU**;
   209	the machine went from **22 GB free to 500 MB**, and the OOM killer then took **four consecutive unrelated
   210	background jobs** — two `grok` legs and two `codex` builds — while earlier runs of the same jobs had
   211	completed. ⛔ None of them produced a Mathematica error, or any error of their own.
   212	⇒ ⭐⭐ **When a background job is killed with a healthy-looking log, check `free -h` BEFORE anything else.**
   213	⭐ `ps -eo pid,rss,pcpu,etime,comm --sort=-rss | head` finds the orphan in one command.
   214	⚠ ⭐ **Its parent will not be a `claude` process** — so "is another session competing?" is the wrong first
   215	question, and it cost an hour of misdiagnosis.
   216	⇒ ⭐ **This is why the budget is enforced ON SIGHT:** the kernel that orphaned was one observed running
   217	**13× over budget** and allowed to finish because it looked like it was converging.
   218	
   219	⚠ **An over-budget kernel does not announce itself** — it looks identical to a busy one. ⭐ `etime ≈ cputime`
   220	on a single kernel process tells you it is **one continuous computation with no idle**, ⛔ not a batch of
   221	short ones. ⚠ **It does NOT tell you whether it will converge** — the measured case above ran **~23 minutes
   222	and then finished normally.** ⇒ ⭐ **Budget, ⛔ not convergence, is the thing to enforce**: a run past
   223	`timeout 600` is a failed ablation whether or not it would eventually return, because ⛔ you cannot tell
   224	which from outside, and waiting to find out is how a seat gets held for half an hour.
   225	
   226	## Launch in Parallel (non-Mathematica artifacts)
   227	
   228	1. Start Grok using the **Bash tool with `run_in_background: true`**, and ⛔ **no shell `timeout`**.
   229	   ⛔ Do NOT use a shell `&` — this harness detaches the job and notifies you on exit; `&` inside a
   230	   foreground call leaves it untracked.
   231	
   232	   ```bash
   233	   grok --prompt-file /absolute/review-prompt.md --cwd /var/projects/toy_physics \
   234	     --model grok-4.5 --effort high --permission-mode bypassPermissions \
   235	     --output-format plain > /absolute/grok-review.txt 2>&1
   236	   ```
   237	
   238	2. In the **same message**, launch the **second leg chosen by the authorship table above** — a Codex
   239	   review for orchestrator-written artifacts, a fresh `general-purpose` Agent for Codex-written ones.
   240	   ⛔ For the agent case use a **fresh** agent, never a fork — a fork inherits the caller's context,
   241	   including the results the reviewer must not see. Give either leg no Grok output and no prior context.
   242	3. Do not poll; the harness notifies you as each leg finishes. Preserve the fresh-agent and Grok reports
   243	   **separately attributed**, and ⛔ do not turn either finding into an edit or rebuild — filter first.

exec
/bin/bash -lc 'git cat-file -e 8f3a017f:research/pde_ledger_v3/scripts/S11c_c2_exports.py 2>&1; echo git_cat_exit=$?; find /tmp /var/projects/toy_physics_ext_logs -xdev -type f -size 60516900c -printf '"'%p %s bytes\\n' 2>/dev/null | head -n 40; git show 8f3a017f:research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md | nl -ba | sed -n '112,132p'; git show aa76105a:research/pde_ledger_v3/_measurements/S11c_c2_export_repair_builder_report.md | nl -ba | sed -n '1,32p'" in /var/projects/toy_physics
 succeeded in 8054ms:
/tmp/S11c_c2_exports.before.py 60516900 bytes
/tmp/s11cc2-review-grok-ablate/S11c_c2_exports.py 60516900 bytes
   112	- `slab_operator_term_origins`
   113	- `traction`
   114	
   115	## Export and measurement evidence
   116	
   117	The declared outgoing roots are `s11cc2ClosedSlabOperator`, `s11cc2ClosedCouplingKernel`, and `s11cc2SelfEnergyIncrement`. Publication computes their own bind-closure over the real parent fold, adds referenced new coordinate/function declarations and their dimension rows, checks exact-key collisions, runs `assert_delta_is_minimal`, and checks serialization against the live values. The generated module pins the five required input digests. The runtime evidence records the exact exported key set and guard operands.
   118	
   119	Production stdout is `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`. Each §5 control has its own literal operand and residual tags there. No residual is asserted against a physics target. The completion marker is `_measurements/S11c_c2_sympy_completion.txt`; `/usr/bin/time -v` measures the Python process in `_measurements/S11c_c2_sympy_resources.txt`; `_measurements/S11c_c2_sympy_monitor.tsv` samples that process's PID, elapsed seconds, RSS, VSZ, and CPU percentage. `_measurements/S11c_c2_sympy_progress.json` records the current case/control.
   120	
   121	`_measurements/S11c_c2_sympy_source_index.json` records the script hash, line count, Python lexical-token count, and exact function spans. The runtime tag index provides the exact computation/emit location and fresh write-key for each emitted object. Prior development attempts and the single-case smoke run are separate `/tmp/S11c_c2_attempt*` / `/tmp/s11cc2_smoke*` files; they are not the production completion evidence.
   122	
   123	## Final measured evidence
   124	
   125	The production completion marker contains `0`. The audit contains 1053 lines, 82,702 bytes, and 22,953 Python lexical tokens (a source-size plausibility measure, not a model-token billing count). Its SHA-256 is `6360721815b0b8c4a52393c388c37b20c858bca04615d8e671964af49a87c7dd`. The plain, non-symlink export contains 71 own rows and 60,516,900 bytes. Publication reloaded the serialized module and emitted its literal equality tuple before the round-trip guard.
   126	
   127	The actual lookup manifest has 34 keys, exactly matching the recorded accesses; its consumer closure has 543 keys. The three-parent fold has 2556 rows. Own delta membership equals the computed own bind-closure, imported-key collisions are `[]`, and the independent file reread records all five digest comparisons in `_measurements/S11c_c2_postbuild_verify.json`. Closed primary pressure/normal-slot references are `[]`; raw c1 whole-form/resolvent references are `[]`.
   128	
   129	The stdout index records 154 unique tags across 498,811,405 bytes, duplicate tags `[]`, and 0 untagged records. Every record carries a computed multigrade and dimension object. The complete literal metadata, branch definitions, and residual values remain in the indexed stdout; they are not converted into physics dispositions here.
   130	
   131	The streaming dimension-metadata scan found missing dimension fields on `[]` and literal `nan` entries on `[]`. This records the dimension walker's output; the independent corruption probe below tests its sensitivity.
   132	
     1	# S11c-c2 export repair — builder report
     2	
     3	Build → verify → report completed. Publication-only; physics construction and emissions are unchanged.
     4	
     5	The generated `scripts/S11c_c2_exports.py` is **22,441,522 bytes**, down from **60,516,900 bytes**: **38,075,378 bytes removed (62.92%; 2.70× smaller)**.
     6	
     7	## Exact edit sites
     8	
     9	Only `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` was edited as source:
    10	
    11	- `EXPORT_ROOTS`, current line 48: removed `s11cc2SelfEnergyIncrement`; retained both declared deliverables.
    12	- New adjacent publication helper `publication_compact`, current line 806.
    13	- `publish`, current line 852: compact only case `VALUE` trees; bounded root-name displays; retain declaration schema, F9, closure, minimality, digest structure, and structural roundtrip; add strict emitted-versus-restored semantics and reciprocal guards, literal evidence and size measurements; install only after passing.
    14	- The sole `export_key` dictionary inside `run`, current line 1051: removed increment export routing. The increment's existing emit remains.
    15	
    16	AST comparison against the pre-edit source passed after excluding only these authorized sites. All other construction, extraction, control, grading, emission, and physics-loop code is identical. No construction change was needed.
    17	
    18	## Representation and size
    19	
    20	Transform: exact `collect` followed by `factor_terms(fraction=False)` on per-case VALUE expressions, with reciprocal/calculus atoms protected locally and restored immediately; Integral boundaries and container types retained; any subexpression whose reciprocal-power set would change keeps its original form. No expansion is used to produce export values, no CSE is used, and no temporary/hold symbol survives in the delivered values. Expansion occurs only in the separate semantic guard.
    21	
    22	`s11cc2ClosedCouplingKernel` VALUE totals: 18,599,401 → 13,917,507 bytes (25.17% smaller). `s11cc2ClosedSlabOperator` VALUE totals: 15,031,549 → 8,422,000 bytes (43.97% smaller). `s11cc2SelfEnergyIncrement` VALUE totals: 11,997,906 → 0 bytes (100.00% smaller; absent from delta, still emitted).
    23	
    24	Sizes below are UTF-8 byte lengths of each case's `srepr(VALUE)`, excluding its unchanged metadata. Baseline sizes were measured from the original 60,516,900-byte module; the final build's emitted sizes match those baseline sizes exactly. Zero means absent from the delta.
    25	
    26	| Object | Case | Before bytes | After bytes | Reduction |
    27	| --- | --- | ---: | ---: | ---: |
    28	| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHO4_CONSTANT` | 4,355,993 | 3,375,460 | 22.51% |
    29	| `s11cc2ClosedCouplingKernel` | `LAB_HELD/RHOBR_CONSTANT` | 3,407,032 | 2,596,031 | 23.80% |
    30	| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHO4_CONSTANT` | 5,990,000 | 4,422,717 | 26.16% |
    31	| `s11cc2ClosedCouplingKernel` | `MATERIAL_ADVECTED/RHOBR_CONSTANT` | 4,846,376 | 3,523,299 | 27.30% |
    32	| `s11cc2ClosedSlabOperator` | `LAB_HELD/RHO4_CONSTANT` | 4,002,816 | 2,271,711 | 43.25% |

exec
/bin/bash -lc "stat -c '%n %s bytes' research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md && nl -ba research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md | sed -n '1,120p' && sha256sum research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md && git status --porcelain=v2 --branch" in /var/projects/toy_physics
 succeeded in 0ms:
research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md 3998 bytes
     1	# Verify — the S11c-c2 compact-prep state (docs accurate? next-prompt correct? clear to compact?)
     2	
     3	The orchestrator is about to compact context. Verify, against the REAL repo + git history, that (1) the state
     4	docs are accurate and not overstated, (2) the carry-forward set is right, and (3) the hand-back prompt is correct
     5	and complete. ⛔ Do not rubber-stamp; if something is wrong or overstated, say so plainly with file:line / commit
     6	evidence. This is a document + light-git verification; ⛔ do not modify the working tree.
     7	
     8	## What just happened (claimed)
     9	Two things closed this session, both committed on `ledger-v3-rebuild`:
    10	- **c2 fold physics** was independently 2-leg reviewed (fresh Claude agent + Grok) → **SOUND, 0 confirmed defects**,
    11	  committed `8f3a017f`. Grok flagged three "defects" (F uniform-limit, G adjointness, E N6) that the orchestrator
    12	  adjudicated as FALSE POSITIVES by its own computation (rule 13).
    13	- **c2 export repair** (publication-only, astra build) shrank `scripts/S11c_c2_exports.py` **60,516,900 →
    14	  22,441,522 bytes**, directive gated `a5f7a06c` (2 decision legs), committed `aa76105a`, re-reviewed clear.
    15	
    16	## 1. STATUS + adjudication accuracy (verify against git + the files)
    17	Read the new top clause of `STATUS.md` (2026-09-06, "c2 STEP A + STEP B/C DONE") and the two adjudication records
    18	`_measurements/S11c_c2_physics_review_adjudication.md` and `_measurements/S11c_c2_export_repair_rereview_adjudication.md`.
    19	Check against `git log --oneline -8`, `git show --stat 8f3a017f aa76105a`, and the artifacts:
    20	- Are commits `8f3a017f` / `a5f7a06c` / `aa76105a` / `d2befb7c` real and as described?
    21	- Is `scripts/S11c_c2_exports.py` actually ~22.4 MB, does it import, is `s11cc2SelfEnergyIncrement` absent, and are
    22	  both `s11cc2ClosedSlabOperator` + `s11cc2ClosedCouplingKernel` present (4-case trees)?
    23	- Is the diff `8f3a017f`→working of `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` genuinely publication-only
    24	  (only `EXPORT_ROOTS`, the `export_key` map in `run()`, `publish`, and the new `publication_compact` helper — NO
    25	  construction function changed)?
    26	- Is "physics SOUND, 0 defects" appropriately supported, or overstated? In particular: are the F/G/E dispositions
    27	  (F genuine coupling decouples / G directional self-energy / E leading-order rep-invariance holds, σ_W deferred)
    28	  defensible from the adjudication's own evidence, or does any read as the orchestrator rationalizing away a real
    29	  leg finding?
    30	
    31	## 2. The carry-forward set (right + complete?)
    32	The step record is claimed to owe: (F) "genuine coupling decouples", not "increment vanishes" (+ a light §5e/§3c
    33	spec-wording fix); (E) σ_W-sector rep-invariance remnant deferred; (G) directional self-energy; the two S11c-b
    34	sign conventions surfaced by the WL comparator. Against the physics adjudication + `directives/S11c_c2_SHARED_PHYSICS.md`
    35	(§5e, §3c, §7): is this set correct and complete, or is a real carry-forward missing / miscategorized (e.g. is any
    36	of these actually a build defect rather than a record note)?
    37	
    38	## 3. The hand-back prompt (correct + complete + clear to compact?)
    39	Read `/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md`.
    40	Verify: the STATE line matches reality; the NEXT sequence is right (STEP 1 light spec-wording fix = review-until-clear;
    41	STEP 2 blind WL engine with 2 decision legs BEFORE + 2 build legs SERIALIZED; STEP 3 T7 comparator + reconcile →
    42	step record); the model policy + ops notes are correct; nothing outstanding is being skipped. Is anything in the
    43	prompt wrong, missing, or an overstatement? Is there any reason we are NOT clear to compact (an uncommitted
    44	load-bearing artifact, an open review leg, a claimed-committed thing that isn't)?
    45	
    46	## Output
    47	For each of 1–3: your finding + evidence (commit / file:line). End with: **CLEAR TO COMPACT** or the exact list to
    48	fix first, and (if any) the precise wording corrections for the STATUS clause / next-prompt.
5f422eacc497deb212c2b22b5200c825229ef7ced55b3cc405770969cea9d095  research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md
5773b328ae03a54bb94dcb4c0fed0807320da0ca20ca4f85eea0e422f9ff1367  /tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md
# branch.oid d2befb7c84016dd7bd9ecac9e1b053c389e6b3ee
# branch.head ledger-v3-rebuild
? research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md

codex
One additional sequencing issue is now clear: `S11c_c2_exports.py` pins the shared-physics spec in `BUILD_INPUT_DIGESTS`. Therefore the proposed “edit the spec, commit, then go directly to Wolfram” sequence would immediately stale the committed export’s provenance digest unless the delta is regenerated (and its publication review updated) after the spec change.
exec
/bin/bash -lc "rg -n --max-columns 1000 'BUILD_INPUT_DIGESTS' research/pde_ledger_v3/scripts/S11c_c2_exports.py && sha256sum research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py && env PYTHONDONTWRITEBYTECODE=1 python -B - <<'PY'
import importlib.util
p='research/pde_ledger_v3/scripts/S11c_c2_exports.py'
s=importlib.util.spec_from_file_location('x',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m)
for k,v in sorted(m.BUILD_INPUT_DIGESTS.items()): print(k,v)
PY" in /var/projects/toy_physics
 succeeded in 10604ms:
directives/S11c_c2_SHARED_PHYSICS.md e0b8e5deb1b2cb321d128c2f3ba65b1209f68fac38d8d5839b2811222e69fe89
scripts/S11c_b_exports.py 053af8042ac0e4e00b04986ce5c19e41b8d16fcf9b99dbbb9aa1d40f258f72b9
scripts/S11c_c1_exports.py 8e29b4a4afc8d0d3cd98c6febff63b25c17f787724415e9fca78db778681533f
scripts/S11c_c2_selfenergy_fold_sympy_audit.py 24d3ecb40c999604d6dc529cd004d62a9db9cb4186d50896663d35d46499d2e3
scripts/ledger_fold.py 97da78c6631c09e5575dd14511760506f1ac08771cf07948978004e92f0b2333

exec
/bin/bash -lc "stat -c '%n %s bytes' /tmp/S11c_c2_selfenergy_fold_sympy_audit.out /tmp/c2_ablate_agent/load_tag.py 2>&1 || true && git ls-files 'research/pde_ledger_v3/scripts/out/*S11c_c2*' 'research/pde_ledger_v3/_measurements/*S11c_c2*.out' | sed -n '1,120p' && rg -n '/tmp/c2_ablate_agent|/tmp/S11c_c2_selfenergy' research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md" in /var/projects/toy_physics
 succeeded in 0ms:
/tmp/S11c_c2_selfenergy_fold_sympy_audit.out 498811405 bytes
/tmp/c2_ablate_agent/load_tag.py 1883 bytes
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.out
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out
research/pde_ledger_v3/_measurements/S11c_c2_physics_review_adjudication.md:4:astra/`gpt-6-astra` build). Emitted output `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out` (499 MB), navigated by
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py:2:import sys; sys.path.insert(0,'/tmp/c2_ablate_agent')
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py:12:import sys; sys.path.insert(0,'/tmp/c2_ablate_agent')

exec
/bin/bash -lc "nl -ba CLAUDE.md | sed -n '176,206p' && rg -n 'canonical.*\\.out|\\.out.*commit|annex.*out|GIN' research/pde_ledger_v3/STATUS.md STATUS.md research/pde_ledger_v3/_measurements/S11c_c2* research/pde_ledger_v3/directives/S11c_c2* 2>/dev/null | sed -n '1,180p'" in /var/projects/toy_physics
 succeeded in 0ms:
   176	## Repository infrastructure — the `.out` transcripts live in git-annex + GIN (set up 2026-09-01)
   177	
   178	The v3 CAS audit transcripts (`research/pde_ledger_v3/scripts/out/*.out` and
   179	`research/pde_ledger_v3/mathematica/out/*.out`, ~370 MB) are **git-annex content backed by GIN**
   180	(`gin.g-node.org/trevnorris/toy_physics`, public), NOT plain git blobs — one exceeded GitHub's 100 MB/file cap.
   181	The policy is the root `.gitattributes` (last-match-wins): everything is `annex.largefiles=nothing` (plain git)
   182	**except** those two `out/*.out` globs, which are `anything` (annex). GitHub is `annex-ignore` by design (git +
   183	tiny pointers only); GIN holds the bytes.
   184	
   185	- **Generating/updating a v3 `.out`**: just `datalad save -m "…"` — the policy annexes `out/*.out` automatically
   186	  and keeps everything else in git. Then publish with **both** `datalad push --to gin` (content → GIN) **and**
   187	  `git push origin <branch>` (git + pointers → GitHub). ⛔ Never `git add -f` a big `.out`. ⛔ Never annex an
   188	  `*_exports.py` — they are hash-chained plain-git inputs the next step imports.
   189	- **After any fresh checkout/clone the `.out` are annex symlinks with NO content** until `datalad get <path>`.
   190	  `grep`-by-line still resolves once content is present (the symlink is followed), so a script/directive/leg
   191	  that reads a committed `.out` must `datalad get` it first — otherwise it reads a dangling link.
   192	- The GIN token is stored in the datalad keyring under credential name **`gin`**; use `--credential gin` with any
   193	  `create-sibling-gin`/publish. Full record, exact commands, and open follow-ups: the
   194	  `project-datalad-gin-out-storage` auto-memory.
   195	
   196	---
   197	
   198	## Evidence ledger — why each control exists
   199	
   200	The measured incident or core rationale behind each control, quoted from the original rule — or **excerpted**
   201	where the rest of that rule's content now lives in its canonical M/E/G/S entry above (so a few entries are the
   202	war-story extract, not the rule's full text). ⛔ Nothing invented. A control whose original text carried no
   203	separate war-story has **no entry here** — its rule above is self-contained (R3, R5, R6, R9, R11, R15, R16
   204	are such). Each control that has an entry references it below as `L-R#`.
   205	
   206	- **L-R2 (E1) — the orchestrator's own typed conclusions.** Measured 2026-08-12: four export-chain designs
research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:183:⇒ S11CC2_CLOSED_SLAB_OPERATOR (assembled, per (α,ρ)) , S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS (per face) ,
research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:199:⇒ S11CC2_CLOSED_COUPLING_KERNEL , S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS ,
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:26:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:53:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:314:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:341:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:602:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:629:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:890:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:917:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
STATUS.md:14:**(3) astra SymPy build RAN** (`gpt-6-astra` high, detached, 19 min, EXIT 0) — the **builder REPORTS** implementing the intended construction (`_measurements/S11c_c2_sympy_builder_report.md`; ⚠ physics **NOT independently verified** — that report even carries `DERIVED_OR_DECLARED: FAIL` + a parser-limitation, and — now stale — claims no reviewer was launched) — **BUT two problems, both from MY directive defects (now fixed `1ae6c336`, refined post-Codex-verify):** (a) **OVERSTEP** — astra read the `/build` skill and ran its own INVALID self-reviews + ~40 orchestrator files (self-check trap; §7 had embedded the orchestrator's review-launch); **quarantined (37 files)** to `/var/projects/toy_physics_ext_logs/c2_builder_overstep/`; [[feedback_builder_directive_no_orchestrator_process]]. (b) **OVER-EXPORT** — 60 MB = 3 fully-expanded closed objects (`ClosedCouplingKernel` 24.7 / `ClosedSlabOperator` 19.9 / `SelfEnergyIncrement` 15.9); §5 didn't split EMIT (→`.out`) from EXPORT (→delta). **Directive fixed:** §7 bounds the builder to build→verify→report; §5 exports **BOTH closed operators** (`CLOSED_SLAB_OPERATOR` full closed op — d needs the closure-modified diagonal for its spectrum — + `CLOSED_COUPLING_KERNEL` off-diagonal, per d's declared scope `S11c_decisions.md:47,52,83`), in a **transparent factored** encoding with a semantic-equivalence check, and the increment (a comparison/emit artifact) + term_origins + §3d + §5 controls EMIT-only (Codex-verify corrected an earlier "coupling-kernel-only" guess). **Artifacts SALVAGED** (`scripts/S11c_c2_selfenergy_fold_sympy_audit.py`, `scripts/S11c_c2_exports.py`, the builder report — **UNcommitted**, pending review + export repair).
STATUS.md:18:Three things landed after the c1 close. **(1) c2 SHARED_PHYSICS spec (the self-energy fold) GATED + committed `16849fc6`** — 2 decision legs (Codex sol xhigh + Grok) both ⚙-computational + CONVERGED, v1 NOT sound (central **isolation claim FALSE** — the substitution increment does NOT cancel S11c-b's two unvalidated SIGN CONVENTIONS; they multiply the substituted δp_s slots, both legs computed `2·carrier·incr`), folded ONE pass → **v2** (9 findings: real θ-row has `J_s` ELIMINATED ⇒ substitute **closed δp_s+w-jets** NOT J_s; `Λ_X` traction-only; **assembled two-face** operator; **3→6 re-adjudications**; positional `load_model` + `s11c_c1_` keys; `Z→0≠Λ_A=0`). Record `_measurements/S11c_c2_spec_review.md`. **(2) c1 SPEC RETRO-REVIEW DONE** (user's *'full correct path'* — the c1 spec was folded once, never re-legged): 2 legs DISAGREED (Grok CLEAR / Codex BLOCK-reopen), adjudicated (rule 13, `_measurements/S11c_c1_spec_retro_review_adjudication.md`) — ⭐ **c1's ENGINES/EXPORTS STAND, ⛔ do NOT reopen** (PY carries density OPAQUELY, 0 derivatives; re-binding to `background_density_map` recovers the O(εη) channel exactly → c2 is **specified** to recover it, c2 not yet built; Codex-sol verify CONFIRMED the override), but record corrections OWED — **3 MUST-level** (seal-5 "harmless O(η²)"→**O(εη)-recoverable-representational**; grazing→`NOT_ESTABLISHED`; independence-scoping) **+ lower-severity carry-forward** (energy-residual orientation, `h_s`/DtN-vs-impedance terminology, density-as-mult-operator, `K_a` Hermitian, η²/ησ_W/σ_W² caveat, drain-projection, flat-`Z₀` leakage; full list in the adjudication record). **(3) CLAUDE.md RESTRUCTURED + committed `3da095c6`** (17 rules→M/E/G/S + an **artifact→review-discipline table**; the fix: **specs=review-until-clear**, not fold-once; ⚠ the preservation legs caught + removed a self-introduced **FABRICATION** in the rewrite). **NEXT = (a) the 3 c1 RECORD corrections (2 legs, review-until-clear; ⛔ NOT engine reopens), THEN (b) c2 SymPy BUILD DIRECTIVE** (own 2 decision legs BEFORE the astra build; determine `IMPORT_KEYS` root set + fold symbol map vs the real files). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = c2) is SUPERSEDED — kept for the artifact map.**
STATUS.md:20:## ⭐⭐⭐ S11c-c1 DONE — T7 comparator + CROSS-ENGINE RECONCILE + step record committed `072d0b75`; NEXT = c2 (2026-09-05) [SUPERSEDED by the clause above]
STATUS.md:21:Both c1 engines per-engine SOUND (clause below). The **T7 comparator** was gated (`84686a54`, 2 decision legs, 12 folds) → built → reviewed baseline `7141e6ad` → scoped repair (astra) `704308af`, all re-review legs CLEAN = **SOUND** (measurement-only, prints/decides nothing). ✅ **CROSS-ENGINE RECONCILE RUN** (family-scoped, 46/50 families, 30 GB, peak ~317 MB) + **c1 STEP RECORD committed `072d0b75`** (`steps/S11c_c1_curved_bulk_closure.md`; record `_measurements/S11c_c1_comparator_reconcile.md`; reproducible `_measurements/S11c_c1_reconcile_reproduce.py`). Method = the **staged representational bridge** ([[feedback_reconcile_representational_bridge]]): comparator PRINTS raw (nonzero≠disagree), I apply identifications read from the engine sources (⛔ not tuned to zero; adversarial corruptions bite) and test collapse. ⭐ **EARNED cross-engine AGREE:** the two-momentum **DtN KERNEL** — all 4 (anchoring×face) cases EXACT ZERO off-diagonal (Stage-2 residual collapsed to the dispersion relation itself) — and the **face-response COEFFICIENTS** (LAB_HELD off-diagonal + MATERIAL on-diagonal); seals 1–4 AGREE. **2 step-record legs (Grok + Codex sol xhigh) SCOPED an earlier "closure-wide AGREE" over-claim** — folded once (`_measurements/S11c_c1_step_record_review_{grok,codex_sol}.txt`). ⛔⛔ **UNDECIDED — ⛔ NOT cross-engine-closed (c2 must NOT inherit these):** (5) **bg-density = a SURFACED rule-17 freeze** (PY bare constant vs WL live field `(ρbr/W₀)WBg(x)`; functional dependence agrees + 0 derivatives, but a scalar≠a field globally — c2 re-adjudication MANDATORY); the raw **`dtn_operator`** whole-form (kernel-AGREE does NOT extend to it); **ENERGY** (closed-form vs far-field integral); an off-diagonal **flat-resolvent leg-labeling**. ⚠ **4 giants + full symbolic residual UNMEASURED — DEFERRED ≥64 GB** (`DEFERRED_HEAVY_RUNS.md`, ⛔ not pre-adjudicated). No cross-engine DISAGREEMENT found. **The kernel is cross-engine closed; the rest is UNDECIDED/deferred.** ⚠ Codex-sol correction-verify pass was run on the fold (in flight/complete). **NEXT = c2** (self-energy fold, a CODE build → gpt-6-astra, own gated directive with 2 decision legs BEFORE the build; carry the mandatory seal-5 + `dtn_operator` re-adjudication). Live handoff: [[project_s11c_c_state]]. **The clause below (NEXT = comparator) is SUPERSEDED — kept for the artifact map.**
STATUS.md:23:## ⭐⭐⭐ S11c-c1 BOTH ENGINES DONE + per-engine reviewed; NEXT = the c1 T7 CROSS-ENGINE COMPARATOR (2026-09-04) [SUPERSEDED by the clause above]
STATUS.md:24:The export migration + delta rebuild (clause below) is **DONE** (`f90e7630`); the c1 SymPy engine is **DONE** (`d6e16471`, 906 KB own-rows delta export); the c1 **blind Wolfram** engine is **DONE** (`c2b0f9dd`→`e139bc61`→`13f0bd2c`→ **repair-2 `dd34d564`**, fully gated, both re-review legs CLEAR). ✅ **A FULL REMEDIATION is COMPLETE** (`_measurements/S11c_c1_wl_remediation_plan.md`): the WL repair DIRECTIVE had skipped its rule-7 decision legs; the remediation re-ran repair-2 **WITH** its decision legs, which **caught that a first-draft fix targeted a non-defect** (would have corrupted the correct `PERMEABLE_PORT_HERMITIAN`) — the gate working. Both c1 `.out` are committed (`4a14100a`, annex+GIN; WL SHA-256 `d63e616786…`). ⛔ **Both c1 engines are BUILT + reviewed — do NOT rebuild.** **NEXT = the c1 T7 cross-engine comparator** — frozen contract `directives/S11c_c1_SHARED_PHYSICS.md:580-587` (N8 inherited verbatim; ⛔ point at it, don't paraphrase): it **computes and prints, deciding nothing** (join by object name with the S11c-a axis-typed keys, injective `mechanical_lower_camel`, paired residual operands + per-family accounting, repoint ablation); the three-valued AGREE/DISAGREE/UNDECIDED characterization is the POST-RUN adjudication, ⛔ not a script verdict; ⛔ no representational fold pre-registered. Reuse `scripts/S11c_a_cross_engine_comparator.py`. Orchestrator-written → Codex+Grok decision legs BEFORE the build → build → 2 legs. Known representational residuals to adjudicate AFTER the run (never pre-registered): μ_θ, the ω-assumption artifact, the two-momentum leg names, and the new proportional-4 parity NIT. ⚠ The full cross-engine residual may need the ≥64 GB box (this box is 30 GB) — build per-engine here, defer the heavy residual if it OOMs. Live handoff: [[project_s11c_c_state]]. **The two clauses below (c1 build PAUSED / NEXT = migration) are SUPERSEDED — kept for the artifact map.**
STATUS.md:42:✅ **S11c-b step record `steps/S11c_b_variable_coefficient_operator.md` committed `bcb9f7d7` + pushed** — the honest close: basis 40 + un-frozen operator + constraint-fold pin B + #90 face+response coupling, all PER-ENGINE leg-verified; the FULL cross-engine `row_residual` (P2a/P2b comparator builds + WL primaries `.out` with term-origins) + #88 re-adjudication + 2 control-hardenings DEFERRED to a ≥64 GB box (user choice — the lighter core-only ~8 GB path was offered, not taken). The 2 whole-row sign conventions (kinetic, face) + #90's 2 flags stay cross-engine-unvalidated. Both document legs (Codex+Grok) confirmed the per-engine-vs-deferred split is HONEST; ~16 accuracy/provenance defects folded (review `_measurements/S11c_b_step_record_review.md`). ✅ `exports.py` regen COMMITTED `af560257` (faithful — `BUILD_INPUT_DIGESTS` match the committed folded+#90 engine/spec; plain-git blob; per user direction). ⚠ its coupling content is per-engine leg-verified, NOT yet cross-engine-validated ⇒ S11c-c imports per-engine-verified coupling. Committed in-tree `.out` remain STALE (fresh ones live in `~/.s11_build/` scratch; regen is ≥64 GB work). **NEXT = S11c-c (curved-interface bulk closure / DtN)** — `directives/S11c_decisions.md`; NO per-substep card (one S11c roll-up after S11c-e). ⚠ Do NOT run 2 memory-heavy CAS jobs concurrently (OOM lesson). Anchors: step record `bcb9f7d7`, deferral `66e8d021`, P1-WL `06048d15`.
STATUS.md:44:## ⭐⭐⭐ S11c-b CROSS-ENGINE RESIDUAL DEFERRED to a ≥64 GB box (STEP 0 OVERTURNED); S11c-b CLOSES on per-engine leg-verification; earlier NEXT was the honest step record (2026-09-03)
STATUS.md:45:⛔ **STEP 0 was WRONG.** It measured only `evaluatedModel` + one `FINAL_KERNEL` (7.95 GB) and concluded the residual fit 30 GB. The full `S11CB_PRIMARIES_ONLY` production emit ALSO builds `mainKernelOrigins` (a FULL `extractCoupling` per origin, ~6/case, for `COUPLING_KERNEL_TERM_ORIGINS`) → ONE case hit **15.6 GB still growing** → OOM. ⇒ **the full cross-engine residual is DEFERRED to a ≥64 GB box** (user decision 2026-09-03: defer whole residual, do not take the lighter core-only path). S11c-b closes on the per-engine leg-verification (fold + #90 + #89/a/b, all 2-leg CLEAR) + the coarse cross-engine consistency established 2026-09-02. The P2a/P2b comparator decision lists are committed as the SPEC for the ≥64 GB run (`DEFERRED_HEAVY_RUNS.md`). ⚠ The two whole-row SIGN conventions (kinetic −K PY vs +K WL; face) + #90's two flags (closure-fold sign; uniform-Λ survivors) remain cross-engine-UNVALIDATED (deferred). NEXT = honest step record `steps/S11c_b_*.md` (2 legs) documenting the arc + the deferrals; exports handling; then S11c-c. Below is the (now superseded) in-flight plan, kept for the artifact map.
STATUS.md:48:- ✅ **P1-WL gate committed `06048d15` + pushed** (origin+GIN synced): additive +17-line `S11CB_PRIMARIES_ONLY` in the WL engine (mirrors PY's), set-mode emits primaries for all 4 cases + skips controls/tower, unset byte-identical to HEAD. **2 build legs CLEAR** (⚠ TWO fresh Claude agents — Grok substituted for the gate's 2nd leg due to an xAI capacity outage; low-risk additive change, objective checks: 14 tags/0 controls, extractors pass, payloads byte-identical set-vs-HEAD, both FORM ablations bite). v1→v2 after 2 decision legs rejected a single-case selector (branch-scoped energy families; `row_residual` raises on non-aligned keys → the fix is the symmetric PRIMARIES_ONLY, both engines emit matching 4-case primaries).
STATUS.md:49:- ✅ **Fresh PY `.out` done** (`~/.s11_build/S11c_b_step0_scope/py_primaries_fresh.out`, 183 MB, folded+#90, PRIMARIES_ONLY). 🔄 **Production WL `.out` running** (`~/.s11_build/S11c_b_production/wl_primaries.out`, task `bgxqpd4ok`, memory-watched, ~8 GB/case as STEP 0 predicted, ~2h). These are the residual inputs (scratch, not committed).
STATUS.md:51:- ⛔ **OOM LESSON (2026-09-03):** ran the WL production run (~8 GB) + the P2b Codex build (parsing the full 156 MB WL + 84 MB PY `.out` into SymPy = tens of GB) CONCURRENTLY → system OOM killed BOTH (watchdog never fired — min avail 8.6 GB; no dmesg access; memory recovered to 23 GB after). ⇒ NEVER run two memory-heavy CAS jobs concurrently on this 30 GB box; SERIALIZE. The incomplete P2b build changes were discarded (committed baseline intact); the P2b directive gained a memory-bound build note (test against a §3a-invariant-only / single-case extract, not the full pipeline).
STATUS.md:53:- **NEXT (serialized, memory-safe):** production `.out` completes → **P2b build** (Codex, bounded, alone) → 2 build legs (fresh Claude+Grok) → commit → **P2a build** (alone) → 2 build legs → commit → **P3: run `row_residual --py <fresh PY> --wl <fresh WL>`** → the cross-engine residual (fold's constraint-reduced rows + #90's closure-fold sign / uniform-Λ survivors; the whole-row sign conventions surface here for step-record adjudication) → #88 KINETIC+θ re-adjudication + 2 owed control-hardenings → honest step record `steps/S11c_b_*.md` + exports (2 legs) → close of S11c-b. ⚠ NO per-substep card (one S11c roll-up after S11c-e); S11c-c/d/e remain. Records: `_measurements/S11c_b_{p1_wl_residual_emit,p2_comparator_update,p2a_slab_row_join,p2b_gamma_bridge}_directive.md`. Post-06048d15 anchors: P1-WL `06048d15`, P2b `69bcc18d`.
STATUS.md:55:## ⭐⭐⭐ S11c-b CROSS-ENGINE STEP 0 (memory scope) DONE — the single-case residual FITS this 30 GB box; integration pass is a 3-builder sequence; NEXT = P1 residual-mode single-case emit (2026-09-03)
STATUS.md:57:- ⚠ **REFRAMING (verified):** `row_residual` and the comparator PARSE the two committed `.out` (`--py`/`--wl`, default the committed files) — they do NOT run the engines live. **Both committed `.out` are STALE:** WL at `d4adbd99` (#89a, pre-#89b → FROZEN operator) and PY at the migration checkpoint (PRIMARIES_ONLY, pre-fold, pre-#90). The stale PY `.out` even has the pre-fold `THETA_BALANCE = μ_θ` structure. ⇒ the integration pass MUST regenerate **fresh single-case `.out`** from both current engines first, then point the comparator at them.
STATUS.md:62:## ⭐⭐⭐ S11c-b #90 §3c COUPLING CONTENT BUILT + leg-gated + committed `7677aa18` (pushed origin+GIN) — the "first physics number" (face+response) is IN; NEXT = cross-engine integration pass (this box) (2026-09-02)
STATUS.md:69:## ⭐⭐⭐ S11c-b PY CONSTRAINT-FOLD BUILT + leg-gated + committed `82f53828` (pushed origin+GIN) — pin (B) implemented; NEXT = single-case cross-engine `row_residual` (this box) → #90 (2026-09-02)
STATUS.md:74:- ⭐ **CROSS-ENGINE ATTEMPT (2026-09-02) → KEY DISCOVERY: the ≥64 GB box is NOT needed for the residual.** A guarded single-case WL U-row run (EULERIAN/LAB_HELD/RHO4_CONSTANT, full basis) FINISHED at **0.9 GB peak** (not 16 GB) — the ~16 GB is only the full run's tower-depth control variants (`operatorActivated/Truncated/Extended` + kernel variants), which the residual doesn't need. ⇒ the whole cross-engine integration (single-case AND likely full-4-case primaries) is doable **on this 30 GB box**, gating only the tower controls. Coarse cross-engine consistency established (both U-rows order-3; both carry live-`W_bg` `eta_bg*w1` coupling — neither froze the constraint). ⚠ RULE 4: a prose shortcut (a clean shared `(η w1 −1)` factor) was REFUTED by the all-jets check (72/75 PY jet coeffs not cleanly divisible). A real COEFFICIENT-level residual needs the full S11c-b symbol bridge (~35 §3a basis coefficients + jet + DOF transliteration) + the `extract_slab` update to the folded PY structure — a reviewed physics-bearing cross-engine JOIN (⛔ never blanket-collapse), doable HERE, deferred to a dedicated integration pass. Evidence + record `~/.s11_build/S11c_b_constraint_fold/{cross_engine_single_case_attempt.md, wl_urow_labheld_rho4.txt, py_urow_labheld_rho4.txt}`.
STATUS.md:82:## ⭐⭐⭐ S11c-b #89b WL OPERATOR UN-FREEZE + REPAIR DONE — engine committed `a1be8d8f` + PUSHED (origin + GIN); full `.out` regen DEFERRED to a ≥64 GB box; PY sibling CLEAR (no freeze); NEXT = integration (reconcile strong-row jet depth) → #90 (2026-09-02)
STATUS.md:83:The §3b variable-coefficient slab operator is un-frozen: coefficients stay LIVE through every order-raising step (EL + `activateSpatialDivergences` + `extractCoupling` + face EL + constraint), `Inactive[Div]` split preserved, jet-retaining reduction LAST (rule 17). Tractability (naively HUNG 2h+ → ~6 min): innermost-first Div activation (no `If`-on-Association) + per-summand `Series` linearity, `PossibleZeroQ`-verified. ⚠ **MEMORY WALL — INTRINSIC:** the correct un-frozen operator must hold its full jet tower un-reduced (~16 GB/case) until the final reduction, so BOTH the 2 heavy equivalence controls (behind `S11CB_SKIP_HEAVY_CONTROLS`) AND the whole 4-case in-band `.out` regen are DEFERRED to a **≥64 GB box** (`research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md`); the committed `.out` is UNCHANGED — **the deliverable is the ENGINE.**
STATUS.md:93:- ⚠ #89a is COMMITTED + **PUSHED** (HEAD `6f5c8c38` atop basis `d4adbd99`; `origin/ledger-v3-rebuild` + GIN both at `6f5c8c38`; the 156 MB `.out` has a `gin` annex copy). The other session's `memory/` work was untouched.
STATUS.md:96:The #89 PY engine repair (un-freeze the background jet tower across the 4 frozen consumers: basis quotient, strong rows, MATERIAL pullback, coupling cascade) is BUILT (Codex, directive `directives/S11c_b_89_sympy_3a_repair_directive.md` committed `660a3055`, 2 decision legs folded) and EMITS the corrected basis **40 = 10+15+15** (confirmed against the #86 reference — established 4 ways — reduces to frozen 26). ✅ **CLEARED by 2 independent build legs** (fresh Claude agent + Grok, both VERDICT CLEAR with agreeing from-scratch derivations: 15/src live [nullity 0] → 40, 8/src → 26, the mandatory 40→26 form ablation bites, and lever-C value-exactness [0 value diffs; 17 srepr-form-only, comparator value-based]; record `directives/_measurements/S11c_b_89_sympy_buildleg_clearance.md`, evidence `~/.s11_build/S11c_b_89_buildleg_{claude,grok2}/` — ⚠ Grok/leg-B saved the full A1–C8 set; the Claude/leg-A dir saved A1–B7 and reported C8 soundness in its verdict). The engine is **COMMITTED** (pre-migration checkpoint `fce14c1a`, rewritten to `f655ea65` by the DataLad migration; `.out` now git-annex-backed, basis 40). ⚠ **Provenance caveat:** the committed `.out` was made with `S11CB_PRIMARIES_ONLY` so the HESSIAN_FREEZE / PROJECTION_EQUIVALENCE / FORM controls are SKIPPED in-band — the 40→26 reduction and lever-C identity are proven OUT-OF-BAND (both legs' ablations/numerics); the in-band control run is deferred pending the kernel-build optimization (full run OOMs).
STATUS.md:103:- **§3c CONTENT verdict — SETTLED: INCLUDE / INCLUDE ⇒ WL spec-correct, PY under-extracts.** The two engines' coupling kernels are structurally near-disjoint at the retained grades: PY = bulk stored-energy only (`gamma·profile-jet`); WL = that bulk PLUS reversible tilted-face geometry (`A_T`, ~1734 terms) PLUS irreversible frequency-dependent response coupling (`A_T·Λ(ω)`, ~660 terms). §3c mandates INCLUDING both: the "weak restriction under the stored-energy/kinetic pairing" is the EXTRACTION INSTRUMENT, not a content filter — §3c forbids the "parallel direct-variation route" (= PY's bulk-only recipe) and "filtering to a single channel", requires `TERM_ORIGINS` to classify face/flux, and its adjointness residual is expressly NOT ∂²U (rules out an energy-only object). Q2 (dissipation): §1c's "not by putting an irreversible response kernel in an ordinary action" is a construction-ROUTE rule not a deletion; §0's "memory kernel" ban targets the S11c-c curved-bulk solve (`δp=Z·v_bulk`), NOT the supplied flat-face `Λ` (verified via the S11c-a **T-i** seam, which states T-i is "not B0c's bulk-response assembly"). Λ stays SYMBOLIC (no bulk elimination). **Adversarially confirmed** by Codex+Grok across 2 rounds (each tried to break it; §3c prohibitions verified verbatim by me). User endorsed the intent. ⇒ repair PY to extract the full operator block + a one-sentence **§0 clarity pin** (Λ = supplied closure, not the S11c-c kernel).
STATUS.md:105:- **PLAN — user chose FIX THE FOUNDATION FIRST** (tasks #86–#90, sequenced): **#86** ✅ **DONE** (corrected basis = 40, settled 4 ways) → **#88** ✅ **DONE** (blast radius SETTLED, commit `05cb1ea5`, record `_measurements/S11c_b_88_blast_radius_result.md`): correcting 26→40 disturbs **every strong stored-energy EL row** at retained grade **non-absorbably** — `RANK_GAIN 8/8/8/6/4 > 0` over the const-coeff span (u-momentum×3, θ, thickness), Hessian-zeroing ⇒ 0, so the non-absorbable content **IS the background curvature** (Hessian). Quadruple-confirmed (orchestrator anchor + Codex instrument + Grok + fresh-Claude legs, all matching; termcounts 98/98/98/108/186). ⇒ the **KINETIC (u-momentum + thickness) and θ-strong verdicts are INVALIDATED** — #89 must re-adjudicate across families, **not only coupling**. Scope: PY-side witness, LAB_HELD stored-energy rows only; mass/ADVECTIVE (constraint `ε·u_t·∇ρ_br`, not an energy-EL) likely spared; WL side + MATERIAL_ADVECTED + admissibility-ε⁰ operator are #89 (a zero here ≠ clearance). OWED (fold at #89 rebuild): harden CONTROL_ENGINE (source engine row independently of the shared density) + CONTROL_JACOBIAN (assert jet-identity). → **#87** ✅ **DONE** (commit `bab2b828`, record `_measurements/S11c_b_87_wl_subspace_result.md`): WL's 8 hand-coded invariants span a STRICT 8-dim subspace of the correct 15 (rank 8/15/15; all 8 in-span; same under the quotient) ⇒ undercomplete by exactly 7, validating the WL fix = complete to 15 → **#89 (PY ✅ CLEARED — see top clause; WL side NEXT)** both-engine §3a repair — fixes DIFFER per engine (PY: retain the Hessian in the quotient — DONE + CLEARED by 2 build legs, commit `9f40c18e`; WL: complete the enumeration to 15 — NOT started); spec-confirm intent → directive → 2 decision legs → build → 2 build legs each → rebuild operator, re-run row-residual + **re-adjudicate KINETIC + θ** (#88); repaired engines must emit 40) → **#90** PY §3c content fix (face+response) + §0 pin, on the CORRECTED basis. ⛔ COMPUTE don't assert (rule 4); ⛔ a (both-)engine fix is a SPEC question FIRST; ⛔ verify basis SPANS not counts.
STATUS.md:139:## ⭐⭐⭐ CURRENT FRONT — S11c-b (variable-coefficient brane operator + off-diagonal coupling kernel; S11c's FIRST physics number). **ADJUDICATION LAYER COMPLETE (3 phases reviewed-sound) + STEP 1 DONE; STEP 2 IN PROGRESS — only ADMISSIBILITY-θ is a SETTLED verdict, advective/kinetic/coupling are LIKELY/PENDING (2026-08-29). THREE commits: `bb0bfc02` comparator time-order fix (canon_jet_name recorded time differentiation as a Boolean `has_time` → WL `D[,{time,2}]`→`u_t` collapsed while PY `u_tt` stayed = asymmetric, flattening transverse-trial ∂²ₜ in COUPLING; count time tokens, emit `_tt`; both build legs SOUND, S11c-a NOT regressed [28/28 diffs canonically zero]; ⭐ discovered BY this adjudication's round-3 decision-list legs); `2e5c6755` adjudication layer v1 `scripts/S11c_b_adjudicated_comparison.py` (Bridge A `bRho↦B_rho_3/W_0` + sort-routing + total-bijection containers + exact-multiset accounting + jet-conservation; both build legs SOUND, 38 MATCH all genuine NO false MATCH via triple ablation, ⭐ Bridge A creates ZERO match — the 38 agreements are RENAME-LEVEL); `afc276b7` adjudication layer v2 (Bridge D = the engine's own committed `PROFILE_GRADE_SUBS` background expansion + IBP/total-in-plane-divergence classifier for weak-pairing DENSITIES ONLY + strong-operator EXACT + atom-gated protection; the v2 directive took 2 decision-list leg rounds catching a naive-chain-rule Bridge D [engines use `sigma_W`/`L_W` not `W_0·eta_bg`], an over-broad divergence classifier [valid only for weak densities not strong operators], and Bridge-D-∂ non-commutation; both build legs SOUND, one reimplementing the Euler/divergence test FROM SCRATCH). ⭐ MEASURED CROSS-ENGINE RESULT (rule 6 — a disagreement IS the measurement): route counts `MATCH=38, PROTECTED_UNREDUCED=32, FLAG=12, RESIDUAL_BULK=8, STRUCTURE_INCOMPLETE=57, COVERAGE=84` (231, multiset-equal, JET 231/0). The adjudication routed 38 rename-level MATCH (ADMISSIBILITY operator(16)/support(20)/COUNT(2)); 32 REPRESENTATIONAL [protected ENERGY_BASIS quotient reps 07/10 + gamma-DivGrad — non-unique by design]; and 20 INDEPENDENTLY-NONZERO operator/coupling differences (not rename/bridge-level — ⚠ step 2 refines which are genuine-physics vs representational): ADMISSIBILITY-THETA `∇²w1` (PY carries, WL 0), SLAB_TERM_ORIGINS ADVECTIVE, KINETIC `W_bg²−W_0²`, 8 COUPLING_KERNEL bulk certified non-IBP by the v2 legs. `REPRESENTATIONAL_DIVERGENCE=0` is REAL (both legs; a from-scratch Euler operator with exact rational zero-test certified all 8 coupling residuals as genuine bulk; the classifier recovers a verified `V` for known divergences up to 5218 ops). ⚠ the original "systematic higher-background-order WL truncation" hypothesis was ⛔ FALSIFIED by STEP 1. **⭐ STEP 1 DONE + COMMITTED `e5bf4122` (2026-08-29): the `(eta_bg,sigma_W)` MULTIGRADE INSTRUMENT `scripts/S11c_b_background_multigrade.py` (both build legs SOUND — fresh Claude agent + Grok each re-derived the coefficients by a DIFFERENT route with ZERO mismatch, form-ablated every guard live; one tautological guard RECONSTRUCTION folded out; extractor `scripts/S11c_b_grade_fingerprint.py`; run `~/.s11_build/S11c_b_multigrade_run.out`; records `_measurements/S11c_b_multigrade_instrument_{build_directive,directive_review,build_review}.md`). ⭐ COMPUTED FINGERPRINT (`(a,b)=(eta_bg^a,sigma_W^b)`; A=PY, B=blind WL) — the 20 differences are a per-family MIX, NOT uniform WL under-retention: ADMISSIBILITY-THETA A={(0,1),(1,1)} B=∅ (WL lacks ∇²w1 entirely); KINETIC THICKNESS_ROW A={(0,0),(1,0),(2,0)} B={(0,0)} (u-momentum rows AGREE); ADVECTIVE A={(0,1)} B={(0,0),(0,1),(1,0)} ⇒ WL carries ∇·u_t that PY LACKS (OPPOSITE direction); COUPLING BIDIRECTIONAL A={(0,1)..(4,1)}+rem(5,1) B={(0,1),(0,2),(1,1),(1,2)}. ⭐ STEP 2 ADJUDICATION (WIP, recorded `_measurements/S11c_b_step2_adjudication.md`; user chose FULL RIGOR all families 2026-08-29): (1) ADMISSIBILITY-θ = GENUINE WL under-retention [SETTLED] — §3a mandates retaining ∇²W_bg (2nd spatial deriv still first shape order), WL operand is 0. Mechanism (leading hypothesis, CONFIRM at repair by running WL): WL's `truncateBackground` Series-to-first-background-order truncation drops it (vs a missing energy term) while PY keeps it at σ_W¹. CLEAR repair candidate. (2) SLAB ADVECTIVE = LIKELY representational provenance re-bucketing (PY = ∇W_bg·u_t, WL = full ∇·(W_bg u_t); PY assigns W_bg∇·u_t to ACCUMULATION/FLUX) — CONFIRM via the deferred ENERGY_BASIS quotient reconciliation. (3) SLAB KINETIC = LIKELY genuine WL under-retention (W_bg²→W_0²; direction PY-has-more). (4) COUPLING = certified non-IBP bulk (v2 legs) with a MEASURED bidirectional grade pattern; the §2a per-cell "which engine is spec-correct" verdict is PENDING. ⛔ NEXT (full rigor): (A) SLAB advective/kinetic — the ENERGY_BASIS quotient reconciliation via VARIATIONAL-EQUIVALENCE (⛔⛔ NOT `classify_total_divergence` — v2 legs proved it invalid for strong operators); (B) COUPLING §2a/§3a per-cell adjudication; (C) then repair ONLY the confirmed genuine engine gaps (⇒ admissibility ∇²w1, likely kinetic; NOT advective if representational) → re-run WL engine + comparator + adjudication → re-adjudicate survivors (rule 13) → honest step record + S11c card + close. ⛔ COMPUTE don't assert (rule 4); ⛔⛔ NEVER blanket-collapse (the four families adjudicate DIFFERENTLY); ⛔ a one-engine fix is a SPEC question FIRST (if §3d/§2a ambiguous for a cell → fix spec). OWED (carry to card): v2 N1/N2/N3; 57 STRUCTURE_INCOMPLETE + 12 control NAMESPACE cross-engine-owed; S11c-a control-family keying; admissibility §5 control-coverage.** **[SUPERSEDED — historical; the reconcile phase below LED HERE. Read as prior state, not the result.] RECONCILE PHASE (2026-08-28): engines review-clean + transcripts committed; T7 comparator BUILT + BOTH BUILD LEGS SOUND (fresh Claude + Grok FORM-ablated every fold) + committed `5f01f9fa`; my run of the COMMITTED comparator is byte-identical to both legs (224727199 bytes; families=28 all-join, 7 unpaired, 0 verdict tokens). One comparator finding (Grok: FACE_FLUX "silent drop") verified + downgraded — visible py_only via COUPLING_KERNEL_TERM_ORIGINS, not integrity. Hand-coded reconcile `S11c_b_handcoded_comparison.py` BUILT + BOTH build legs SOUND (map COMPLETE, no false MATCH, `--drop-rename` surgical/load-bearing) + committed `82ec3b5f`; run (`~/.s11_build/S11c_b_reconcile_run.out`): **2 MATCH** (ADMISSIBILITY_SUPPORT 20/20, ENERGY_BASIS_COUNT bare-int 26), **14 FLAG** (all core objects), **12 NAMESPACE_INCOMPLETE** (control families, WL `<|…|>` unparsed → OWED). DIMENSIONS: the digest shows the dimension VALUES agree (a container-shape diff), but its reconcile verdict is still `FLAG` — confirm in adjudication, do not pre-call it. Both legs verified the 14 FLAGs are NOT naming noise (map complete); representational-vs-finding is the pending adjudication, not yet computed. **⛔ NEXT + FINAL = the ADJUDICATION (decides the number; the never-blanket-collapse danger zone):** apply the further REVIEWED bridges + re-check zero — (i) bRho→`B_rho*W_0` scale, (ii) gamma-index selected-representatives, (iii) WL `f(xOne,xTwo,xThree,time)`→PY bare GATED on a computed 0-spatial-jet check per field (S11c-a INERT_APPLIED gate), (iv) integral linearity; ⛔⛔ NEVER collapse `W_bg`/`mu_R_bg` (real background-gradient jets = the variable-coeff physics; blanket applied→bare hid the S11c-a in-plane current freeze). Reduces-to-0 ⇒ engines AGREE on that object (coupling-kernel number); survives ⇒ genuine finding (rule 1/6). ENERGY_BASIS = non-unique quotient, never fold a representative. Should be a REVIEWED artifact (reconcile-v2 with the gated bridges → 2 legs → commit) — COMPUTE, don't assert (rule 4). Likely mirrors S11c-a: agreement up to reviewed representational identities.** Journey: spec drafted → 2 legs (8 defects) `7023420d` → both blind engines built (Codex) → 4 build legs found a CONFIRMED rule-7 shared-blind-spot (both engines independently produced a vacuous ε→0 admissibility from an ambiguous §3d — a comparator alone would have read it as agreement) + a 4-leg tautological Clairaut adjoint + per-engine kernel-extraction bugs → my spec round-2 fix was STILL insufficient (2 re-legs) → **rule-15 author change: Codex rewrote §3a/§3c/§3d/§5a, 2 re-legs CLEAN** `0c0e9a8a` → engine repairs: WL W1-W3 clean both legs `49d9fad1`; SymPy B2/B3/B4 clean, B1 took a round-2 (the round-1 legs DISAGREED — Grok caught a scalar over-promotion the Claude-agent leg missed; the sibling WL engine's independent-correct §3d construction was the tiebreaker → SymPy-only fix) → **SymPy B1 round-2 clean both legs `103cdfdb`**. Both engines now review-clean. **OWED (non-blocking, carry to the card):** the admissibility operand is verified correct but has no EMITTED §5 discriminating control (independence omits on structural absence, rep-invariance is a background-order structural zero, §5b/§5c don't name it) — a control-coverage refinement, physics verified; plus S11c-a's owed control-family keying. **[SUPERSEDED — the live NEXT is the gated-bridge adjudication in the RECONCILE-PHASE front clause above; comparator + reconcile are built/legged/committed/run. The following is the historical plan record.]** Transcripts committed (SymPy `55abd09b` 167MB, WL `73b4e639` 97MB); T7 comparator (Codex `bqj3ftjqc`, brief `directives/S11c_b_comparator_build_directive.md`) was verified + 2 build legs + run + reconciled per the front. Finding records `_measurements/S11c_b_*`. [Historical spec-build detail below.] **SPEC DONE + 2 legs folded + committed `7023420d` (2026-08-27).** The SHARED PHYSICS spec `directives/S11c_b_SHARED_PHYSICS.md` was authored (mirrors the S11c-a template; inherits S11c-a §§1–3 by reference — SymPy imports `S11c_a_exports`, blind Wolfram re-derives). Two decision-list legs (orchestrator-written → Codex + Grok) ran BEFORE any engine (rule 7) and were **very productive**: they converged on 2 serious shared-blind-spots and surfaced 6 more — all 8 verified against source (rule 13) and folded ONCE (not iterated to green). The 8: **D1** ζ_s=ζ_c+sδW/2 not ζ_±=±δW/2 (was forcing ζ_c=0); **D2 (serious)** confine "u only through gradients" to the uniform quotient — the non-uniform background admits undifferentiated-u spurion couplings, the very N15 channels the step emits (my gloss would have made both engines agree the coupling is absent); **D3** declare which quantities vary (W_bg,μ_R,bg,density) vs constant moduli (N12); **D4** enumerate ONE thickness coordinate + impose e_W,bg map before the rank test; **D5 (deep)** the total-divergence quotient does NOT lift to variable coefficients (`c∇·F≡−(∇c)·F` ⇒ first-jet terms are PHYSICS) — dropped the leaky μ_⊥ citation + removed the pre-registered representative fold from the comparator (could have masked the coupling); **D6 (serious)** admissibility = background-order (ε⁰) balance in the 𝒮_hold⁰ pairing, not ε→0 of the wave operator (vacuously 0) — new §3d; **D7** sectors via local curl/div structure not a global projector (N5), adjoint w.r.t. the supplied pairing; **D8** separate one-source form ablations for ∂W_bg and ∂μ_R,bg (independent profiles). Finding record `directives/_measurements/S11c_b_spec_review.md`; leg prompt `_legs/S11c_b_spec_review.md`; raw transcripts `~/.s11_build/S11c_b_spec_review_{grok,codex}.txt`. Rule-5 leak scan clean (no coupling grade/sign, basis count, or μ_⊥ identity stated). **[historical, DONE — the plan at spec-close was]** the two blind engine builds against the spec (SymPy chained on `S11c_a_exports`; Wolfram blind), each with 2 build legs, then this sub-step's own T7 comparator + step record. OWED (carry in from S11c-a, not blocking): bridge the control-family keying (REP_INVARIANCE/CONTROL_INDEPENDENCE PY-missing-DENSITY; CONTROL_FORM FACE-granularity) + a CONTROL_FORM re-characterization.
STATUS.md:143:⛔ **[SUPERSEDED / OVERCLAIMED — see the corrected front above] The block below reported a clean "THE TWO ENGINES AGREE, every residual representational." That conclusion depended on the blanket collapse that hid the #4 in-plane current freeze and the FACE_SHIFT non-join. Read it as the (partly wrong) prior state, not the result.**
STATUS.md:145:⭐⭐⭐ **T7 CROSS-ENGINE COMPARATOR — BUILT + 2 BUILD LEGS PASS + COMMITTED `50f43123` (2026-08-27); RESULT: ENGINES AGREE.**
STATUS.md:161:  raw. Non-blocking: PROJECTION_TERM_ORIGINS labels DYNAMIC/STATIC positionally (symmetric → can only surface a
STATUS.md:307:⭐⭐ **S11c-a SymPy ENGINE CLOSED — committed `9b6438fa`** (build `buo60i510` → reviewed baseline `488c2a65`
STATUS.md:320:⭐⭐ **S11c-a blind WL ENGINE — COMMITTED `ddecdbc2` (repair-round-2 done + verified).**
STATUS.md:330:/ 12 changed (exactly T-0+5×PROJECTION_*+6×CONTROL_FORM/UNIFORM_LIMIT). ⚠ LESSON: my "88 ∫1 in ORIGINS" alarm
STATUS.md:332:⭐⭐⭐ **T7 CROSS-ENGINE — the engines DISAGREE on CASE STRUCTURE; user chose FULL RECONCILE; fixes land on BOTH engines.**
STATUS.md:344:_RESIDUAL/_TERM_ORIGINS [COMPUTATION]); (C) virtual-work → **WL correct** (full grid §4 T-d:419; PY missing
STATUS.md:525:⭐ CERTIFIED (both legs + orchestrator, rounds 2–4): reducer object-level arithmetic; population reconciliation 0-gap; all 14 repair classes live by per-class ablation; every plant able-to-fail at the prior commit; every verdict transition traced to a repair class; both round-4 Required semantics load-bearing. ⭐⭐ GENUINE CENSUS CORES (now certified measured facts feeding the engine round): **917 decided-undecided (BOTH ENGINES UNDER-DECIDE)**; 136 WL + 35 PY = 171 spurious (XKIN STACKED_DROP/rank-drop; sampler-decorrelation catches, 4 new D3 tags); 104 membership-driven witness failures (= register `_REAL_WITNESS` defect); **omitted 72 = D3/D4 ROOT_COINCIDENCE (2 genuine families each) + 70 PY**; D2 ROOT_COINCIDENCE now certified honest `COMPLETENESS_UNDECIDED` (NOT an omission).
STATUS.md:530:⭐⭐⭐ WL FIX ROUND 2 CLOSED `a4cf6539` — THE CANONICAL WL RECORD IS 21/21 COMPLETE, SPLICED, AND LEG-REVIEWED. Build self-reported FAILED per obligation 5 (final-code D4 memguard-killed at STRATUM10, kernel RSS 24.9 GB against the 1 GiB floor; the earlier "complete" run had only ~122 MB margin — RSS sidecar profiles measured the accumulation both times, ~17.6 GB entering the strata, +1.3 GB/stratum, peak ~25 GB). ⭐ USER RULING (2026-08-16): a guard-floor margin kill is an infrastructure limit, NOT a script failure — D4 re-run isolated with the floor at 256 MB (variant runner `~/.s11_build/fix2_final/run_guarded_cell_floor256.sh`; pinned original untouched): COMPLETE guard=NONE (5,412 em, 4,579 s, min available 755 MB — it dips ~290 MB BELOW the old floor, so the old floor was the only killer). D2 followed (967 s, 1,898 em, 7 GB margin). Final-engine records verified vs pre-patch runs: tag sequences identical; all payload diffs = the time-vs-memory expiry race (semantic content equal) or the final patch replacing arithmetic-over-`Failure` with clean propagation. Splice: canonical `mathematica/out/S11_...audit.out` = 16,587 lines, 21 RUN_PAIRS, untouched regions byte-identical (verified by diff). TWO SCRIPT LEGS run serialized on `d6316978` (fresh agent + grok-4.6, identical committed prompt `S11_wl_fix2_script_review_prompt.md`): round-2 diff SURVIVED — fallbacks uniform at every call site, measured-expiry-only selection, no identity gates, form ablations moved physics (721/455-line diffs). Obligation-6 HALT (44→16/12 arity on `ROOT{2,3}_RANK_DROP_JOINT`) DISCHARGED by both legs + orchestrator: same solution content, admissible sets identical, coverage improved (20 UNDECIDED→decided), the added `{sRho->0}` is a cleared-denominator artifact EXCLUDED downstream. ⭐ FOUR NEW DEFECTS REGISTERED (entries 4–7, every one verified by orchestrator re-computation; none is a round-2 regression): omitted `Solve` branch in `ROOT_COINCIDENCE_COEFF` (confirmed `{muR==bComp, sRho==1}` residuals 0); spurious branch `{muR->bComp}` in `RANK_DROP_JOINT` failing 1/16 minors both sheets (engine's own STRATUM2 rank=CONSTANT 3 corroborates; grok's `{sRho->1}` claim REFUTED 16/16; ⭐ SymPy emits a 2-branch complex k-only object for the SAME tag — the T7 comparator will flag this family); `_REAL_WITNESS` quantifies away non-solve parameters (witness fails emitted equations at generic k); latent route gaps (empty-factor-support degradation, 0 live occurrences; missing SOLUTION_ATTEMPTS on 3 secondary paths). ⚠ OBLIGATION 4 NEVER RAN: both instruments (WL probe census, SymPy containment) parse-fail on the real record format (0 probes ever executed; probe census exits 0 while dying) — calibrated only on planted synthetic records; the legs' independent containment substituted for the HALT loci + 5 cheap loci only. The instrument repair belongs to the defects round. Archives: `~/.s11_build/fix2_final/final_engine_complete_runs.tar.gz`, `fix2_d4_rerun_memguard_death.tar.gz`, `fix2_first_complete_runs_snapshot.tar.gz`, grok report `fix2_final/grok_script_leg.txt`, agent leg scripts `/tmp/claude-1000/s11_fix2_leg_agent/`. NEXT: (a) the registered defects' round (now SEVEN entries, incl. the obligation-4 instrument), (b) the FROZEN COMPARATOR CONTRACT T7 (semantic compare for expiry payloads + additive provenance + WL-only strata + entry-5's known cross-engine divergence), (c) S11 closes.
STATUS.md:534:Prior top block: BUILD MID-FLIGHT, ALL 21 WL CELLS COMPLETE FOR THE FIRST TIME (uncertified). Diagnosis (2 analysts, verified): D4 died in the unbounded `Solve` at wl:312 (one terminal call; the object decides in <1 s by exact routes); D2 died in `assumedRank[stacked]` (wl:1076, `FullSimplify` ZeroTest) after ~1,223 s of unbounded-call ACCUMULATION; six unbounded `Solve` sites + zero `MemoryConstrained` existed engine-wide. Brief folded ONCE `a8f26909` after two BLOCKING legs (Codex 7 + grok 4 findings, 7 defective-repair constructions, all verified): hard D2/D4 completion (no measurement exception), killer classes extended to the whole live path, route-uniformity generalized to fallbacks (+ out-of-cell exercise per branch), no return coercion (spec:245/271), armed probe (live operands on undecided count-class records, no-starvation budgets, planted-record calibration, residual recompute, SymPy completeness containment), guard sha256 pinned, RSS instrumented. THIRD defect REGISTERED (additive count-payload provenance extension — deliberate, round-1 fold). Codex builder launched `b5a9cbae`; measured so far: D3 269 s / 2,501 (byte-shape = baseline, rerun reproducible), ⭐ D4 COMPLETE FIRST TIME EVER (rc=0, guard=NONE, 5,412 emissions, 4,603 s), ⭐ D2 COMPLETE (rc=0, 1,898 emissions, 964 s — faster than its death run). ⚠ WIP SAFETY SNAPSHOT `14b2d56e` (user-requested) holds the repaired engine — NOT certified, no legs yet; first-complete records archived `~/.s11_build/fix2_first_complete_runs_snapshot.tar.gz` (guard runner OVERWRITES per stem). ⛔ REMAINING: builder's acceptance outputs (probe census after planted-record calibration, manifest census, partial-record compare, 19-cell byte regression, RSS profiles) → orchestrator verifies the DELIVERABLE (not the self-report) → archive builder scratches → TWO SCRIPT LEGS (fresh agent + grok, SERIALIZED, form ablation) → verify findings → scoped legs on any post-fold delta → splice canonical .out → real commit replaces the WIP label → STATUS/memory. Then: (b) the THREE registered defects' round, (c) the FROZEN COMPARATOR CONTRACT (semantic compare for expiry payloads + additive provenance fields + WL-only strata surface: SymPy `STRATUM_ORDERING` is EMPTY on all 21 cells — measured).
STATUS.md:545:⭐ **(done — the engine was built and has had three fix rounds; see the S11 PY ENGINE block below.)**
STATUS.md:556:⭐⭐ **S11 PY ENGINE — ROUND 3 CLOSED `9fb45365` (brief folded `1b2f8cf9`). ⭐ D1, D2 AND THE PUBLISH ARE
STATUS.md:609:`~/.s11_build/wl_percell_driver.sh` concatenating into the canonical `.out` (engine output only;
STATUS.md:623:forms, values identical); commits with the successful sweep. ⛔ The `.out` is committed only from a
STATUS.md:632:⚠⚠ **THE PREVIOUS ENGINE DID NOT "MISS" §Q8b — ⛔ there was no spec.** `S11_SHARED_PHYSICS.md` first exists
STATUS.md:675:⭐⭐ **THE WL SIBLING ENGINE IS BUILT (BLIND) AND LEG-CLEAN** — directive `b8395704` (2 legs, 8 defects
STATUS.md:748:⛔⛔ **S11 STARTS WITH A SPEC REPAIR, ⛔ NOT AN ENGINE.** ⚠ Its 914-line spec is **closed and incorrect** in
STATUS.md:859:⛔⛔ **DUAL-ENGINE AGREEMENT IS VACUOUS WHERE BOTH SIDES ARE HAND-DECLARED LITERALS — recorded
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1417:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1429:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1545:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1557:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1673:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1685:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1801:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_guard_evidence.json:1813:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:176:- `S11CC2_CLOSED_SLAB_OPERATOR` (assembled two-face), `…_TERM_ORIGINS` (per face), `…_PARITY_BLOCKS` — §3a.
research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:177:- `S11CC2_CLOSED_COUPLING_KERNEL` (both off-diagonal blocks), `…_TERM_ORIGINS`; `S11CC2_SELF_ENERGY_ADJOINTNESS_RESIDUAL` (only if an independent route exists) — §3b.
research/pde_ledger_v3/directives/S11c_c2_sympy_build_directive.md:211:binds it), every `…_TERM_ORIGINS` / `…_PARITY_BLOCKS`, the six §3d re-adjudication objects, and every §5 control
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:37:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:81:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:523:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:567:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:1009:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:1053:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:1495:    "PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_object_provenance.json:1539:    "PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT": {
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:219:| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsLabHeldRho4Constant` | 4024863 / 4396607 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:222:| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHO4_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsLabHeldRho4Constant` | 16321247 / 4189623 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:251:| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsLabHeldRhobrConstant` | 115414017 / 3056958 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:254:| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_LAB_HELD_RHOBR_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsLabHeldRhobrConstant` | 124751968 / 2904370 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:283:| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsMaterialAdvectedRho4Constant` | 198222339 / 5050165 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:286:| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHO4_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsMaterialAdvectedRho4Constant` | 213313761 / 9173840 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:315:| `PY_S11CC2_CLOSED_SLAB_OPERATOR_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 550, 560 | 893 | `s11cc2ClosedSlabOperatorTermOriginsMaterialAdvectedRhobrConstant` | 338897151 / 3576380 |
research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:318:| `PY_S11CC2_CLOSED_COUPLING_KERNEL_TERM_ORIGINS_MATERIAL_ADVECTED_RHOBR_CONSTANT` | 887 | 893 | `s11cc2ClosedCouplingKernelTermOriginsMaterialAdvectedRhobrConstant` | 350634273 / 7145727 |

codex
The mechanical/git claims mostly check out, but the state is overstated in two material ways: E/N6 is unresolved, not a demonstrated false positive; and the export repair never received the required second usable review report. The proposed spec-first sequence would also stale the export’s pinned digest.

## 1. STATUS and adjudications

### Repository and export facts — verified

- The commits are real and ordered as claimed:

  - `8f3a017f`: physics-review baseline and adjudication.
  - `a5f7a06c`: export-repair directive plus two decision reports.
  - `aa76105a`: final compact export and publication changes.
  - `d2befb7c`: current `HEAD`; changes only `STATUS.md`.

- The current export is exactly `22,441,522` bytes. The pre-repair file is not in commit `8f3a017f`, but an exact `60,516,900`-byte copy still exists at `/tmp/S11c_c2_exports.before.py`, and the original builder report records that size at [S11c_c2_sympy_builder_report.md:123](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:123).

- A real import succeeded in 10.9 seconds. Runtime inspection found 70 rows, no `s11cc2SelfEnergyIncrement`, and four `(anchoring,density)` cases under each required root. The roots are serialized at [S11c_c2_exports.py:11](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11) and [S11c_c2_exports.py:12](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:12).

- The `8f3a017f`→working-tree script diff is genuinely publication-only. Zero-context hunks occur only at:

  - `EXPORT_ROOTS`: [S11c_c2_selfenergy_fold_sympy_audit.py:48](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:48)
  - new `publication_compact`: [line 806](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:806)
  - `publish`: [line 852](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:852)
  - the export-key map in `run()`: [line 1051](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1051)

  No construction function changed.

### F — defensible, with the already-acknowledged spec correction

The adjudication does refute Grok’s alleged surviving closed `Z₀μθ` coupling. All four cases reduce the closed contribution to `Integral(0,…)`, and `.doit()` returns zero: [S11c_c2_adjudication_verify_F.out:1](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out:1).

But the literal increment remains nonzero because of `−extract(open)`. Therefore current §5e is false as written—“the off-diagonal self-energy increment must vanish” at [S11c_c2_SHARED_PHYSICS.md:321](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:321). Treating this as a wording/representation defect is reasonable; calling the literal increment zero would not be.

### G — directionality is supported; “dissipative” is not

The evidence establishes that the increment is one-way while the full closed coupling kernel contains both directions: [S11c_c2_adjudication_verify_EG.out:1](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out:1). The spec requires an adjointness residual only when there is a genuine independent route: [S11c_c2_SHARED_PHYSICS.md:187](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:187). Grok attributed an adjoint-by-construction claim the builder did not make.

However, neither the directionality check nor builder line 70 establishes dissipativity or passivity. Remove “dissipative directional” from [STATUS.md:6](/var/projects/toy_physics/STATUS.md:6) and “a dissipative directional self-energy” from the next prompt. Use “directional/one-way SymPy increment” only.

### E — not defensibly adjudicated as a false positive

This is the major physics overstatement.

The computation proves only:

- the full residual is nonzero;
- setting `σ_W=0` kills it.

See [S11c_c2_adjudication_verify_EG.out:15](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out:15). The actual retained grades are `O(εσ_W)` and `O(εησ_W)`: [S11c_c2_sympy_stdout_index.json:1205](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:1205).

That does not establish N6. The governing spec says the residual must vanish after field redefinition: [S11c_c2_SHARED_PHYSICS.md:303](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:303). It also inherits `η` and `σ_W` as independent first-shape-order bookkeepers; `O(εσ_W)` is retained physics, not higher order: [S11c_c1_SHARED_PHYSICS.md:223](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:223).

The adjudication’s association with the drain-projection `O(σ_W²)` deferral is especially unsupported: a residual linear in `σ_W` is not that second-order correction. The only explicit ≥64 GB deferral in the c2 scope is the full cross-engine/four-giant-family work, not this already-emitted single-engine N6 residual: [S11c_c2_SHARED_PHYSICS.md:39](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:39).

Correct disposition: **N6 is established only after projecting to `σ_W=0`; the retained σ-sector is unresolved and may be a `representation_pullback` build defect.** Grok’s E finding is not disproven.

Therefore “physics SOUND, 0 confirmed defects” and “F/G/E ALL FALSE POSITIVES” in [STATUS.md:6](/var/projects/toy_physics/STATUS.md:6) are overstated.

### Export re-review — technical evidence good, process not clear

The artifact itself looks technically sound. But [S11c_c2_export_repair_rereview_adjudication.md:17](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md:17) explicitly calls the fresh-Claude leg incomplete. It never produced a usable full report or CLEAR verdict. The record substitutes the orchestrator’s verification and says the two-leg discipline was met “in substance” at [line 45](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md:45).

That contradicts repository policy:

- both reports must be usable: [CLAUDE.md:46](/var/projects/toy_physics/CLAUDE.md:46);
- orchestrator verification is never a leg: [CLAUDE.md:119](/var/projects/toy_physics/CLAUDE.md:119);
- no commit before both legs report: [CLAUDE.md:137](/var/projects/toy_physics/CLAUDE.md:137).

Thus “2 re-review legs” and “re-reviewed clear” are procedural overclaims. A replacement full fresh-Claude publication review is owed.

## 2. Carry-forward set

| Item | Correct classification |
|---|---|
| F | Correct carry: genuine closed coupling decouples; literal increment does not. Fix §3c/§5e wording. |
| G | Correct measured SymPy feature, pending the blind-WL comparison. Do not call it dissipative. |
| Two S11c-b signs | Correct comparator carry; §3c and §7 explicitly say surface, do not normalize: [§3c](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:203), [§7](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:394). |
| E | Misclassified. It is an unresolved retained-order N6 finding/potential build defect, not merely a deferred record note. |

If the list is intended to be exhaustive for the step record, it is also incomplete. The blind WL/comparator still owes the §3d questions: density field-vs-field, traction representation, DtN whole-form versus kernel, traction/slab power pairing, flat-symbol usage, and the `μ_R,bg` form control: [S11c_c2_SHARED_PHYSICS.md:224](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:224). The prompt’s generic “§3d representation questions” partly covers these, but should name them.

Also preserve the fact that c1 ENERGY remains UNDECIDED. The c2 spec says it enters load-bearingly at [S11c_c2_SHARED_PHYSICS.md:84](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:84), while the builder explicitly says the c2 traction pairing does not adjudicate the separate c1 far-field ENERGY audit at [S11c_c2_sympy_builder_report.md:70](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:70). At minimum the step record must retain that deferral and avoid claiming energy closure.

## 3. Hand-back prompt and compaction readiness

The post-gate sequence—blind WL directive, two decision legs before build, serialized fresh-Claude/Grok build legs, T7 comparator/reconcile, then two record legs—is otherwise correct. The model-role policy and Mathematica/setsid/same-process notes are consistent with the repo.

But the prompt is not ready:

1. Its STATE repeats the unsupported “SOUND,” “all E false-positive,” and “re-reviewed clear” claims at [next_prompt_c2_wolfram.md:3](/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md:3).

2. It skips the N6 resolution and missing fresh-Claude export-review leg.

3. Its STEP 1 would edit the shared spec and then go directly to WL. That makes the committed export provenance stale: `S11c_c2_exports.py` pins the current spec digest at [S11c_c2_exports.py:9](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:9), as required by [S11c_c2_SHARED_PHYSICS.md:386](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:386). After any spec edit, the export must be regenerated or lawfully repinned and reverified before advancing.

4. The live `project_s11c_c_state` memory repeats the same E and export-clear overclaims in its frontmatter and closeout: [project_s11c_c_state.md:3](/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:3), [line 201](/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:201).

5. The working tree is not clean. The only dirty path is the untracked [S11c_c2_wolfram_compact_prep_verify_prompt.md](/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md:1). It is this audit prompt, not an uncommitted build product, but if it is intended as review evidence it must be committed with the resulting record or intentionally discarded.

6. The 498,811,405-byte reviewed SymPy stdout exists only at `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`, and both committed adjudication scripts depend on uncommitted `/tmp/c2_ablate_agent/load_tag.py`: [verify_F.py:12](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py:12), [verify_EG.py:2](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py:2). That is non-durable load-bearing evidence under the repository’s annex policy: [CLAUDE.md:176](/var/projects/toy_physics/CLAUDE.md:176).

### Precise STATUS correction

Replace “SOUND, 0 confirmed defects / ALL FALSE POSITIVES” with:

> The c2 fold wiring and A/C/D constructions are supported. F is resolved narrowly: the genuine closed coupling decouples, while the literal §3c increment retains the open-slot term and requires a §3c/§5e wording correction. G establishes a directional SymPy increment; no dissipativity claim is made. E remains OPEN: the N6 residual is nonzero at retained `O(εσ_W)`/`O(εησ_W)` grades; `σ_W→0` proves confinement only, not full representation invariance. Resolve as a possible `representation_pullback` defect or through a separately reviewed scope/spec decision before the WL build.

Replace “2 re-review legs / re-reviewed clear” with:

> Export evidence currently supports the repair: Grok CLEAR plus same-process orchestrator verification found no technical defect. The fresh-Claude attempt was incomplete; because orchestrator verification is not a leg, one usable fresh-Claude re-review remains owed before procedural clearance.

### Precise next-prompt sequence

Prepend before the current STEP 1:

> **STEP 0 — close the open gates.** Resolve E/N6 without presuming deferral: determine whether the retained σ-sector residual is a `representation_pullback` defect or requires a physics-bearing scope/spec amendment; review any script/spec change until clear. Obtain a complete usable fresh-Claude publication-review report. Do not launch the WL builder while either gate is open.

Change STEP 1’s ending to:

> After the reviewed spec correction, regenerate or lawfully repin `S11c_c2_exports.py` so `BUILD_INPUT_DIGESTS` matches the accepted spec/script bytes; verify unchanged root membership and semantics, and complete the required publication re-review before committing the refreshed delta.

Add to the carry paragraph:

> The SymPy §3d dispositions remain pending blind-WL comparison, including density, traction representation, DtN whole-form, traction/slab pairing, flat-symbol usage, and `μ_R,bg` form control. c1 ENERGY remains explicitly UNDECIDED/deferred. Do not preload these results into the blind WL packet.

Remove both occurrences of “dissipative.”

**NOT CLEAR TO COMPACT — fix first:**

1. Correct E/N6 from “false positive/deferred” to OPEN and resolve it before WL.
2. Obtain one complete usable fresh-Claude export-repair review leg.
3. Correct `STATUS.md`, the next prompt, and `project_s11c_c_state`.
4. Add the mandatory export regeneration/repin step after the spec edit.
5. Preserve the 499 MB stdout and its loader durably, or make the committed verification scripts self-contained.
6. Commit the untracked audit prompt with its result, or explicitly discard it as non-record scratch.
209,237
The mechanical/git claims mostly check out, but the state is overstated in two material ways: E/N6 is unresolved, not a demonstrated false positive; and the export repair never received the required second usable review report. The proposed spec-first sequence would also stale the export’s pinned digest.

## 1. STATUS and adjudications

### Repository and export facts — verified

- The commits are real and ordered as claimed:

  - `8f3a017f`: physics-review baseline and adjudication.
  - `a5f7a06c`: export-repair directive plus two decision reports.
  - `aa76105a`: final compact export and publication changes.
  - `d2befb7c`: current `HEAD`; changes only `STATUS.md`.

- The current export is exactly `22,441,522` bytes. The pre-repair file is not in commit `8f3a017f`, but an exact `60,516,900`-byte copy still exists at `/tmp/S11c_c2_exports.before.py`, and the original builder report records that size at [S11c_c2_sympy_builder_report.md:123](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:123).

- A real import succeeded in 10.9 seconds. Runtime inspection found 70 rows, no `s11cc2SelfEnergyIncrement`, and four `(anchoring,density)` cases under each required root. The roots are serialized at [S11c_c2_exports.py:11](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:11) and [S11c_c2_exports.py:12](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:12).

- The `8f3a017f`→working-tree script diff is genuinely publication-only. Zero-context hunks occur only at:

  - `EXPORT_ROOTS`: [S11c_c2_selfenergy_fold_sympy_audit.py:48](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:48)
  - new `publication_compact`: [line 806](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:806)
  - `publish`: [line 852](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:852)
  - the export-key map in `run()`: [line 1051](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py:1051)

  No construction function changed.

### F — defensible, with the already-acknowledged spec correction

The adjudication does refute Grok’s alleged surviving closed `Z₀μθ` coupling. All four cases reduce the closed contribution to `Integral(0,…)`, and `.doit()` returns zero: [S11c_c2_adjudication_verify_F.out:1](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out:1).

But the literal increment remains nonzero because of `−extract(open)`. Therefore current §5e is false as written—“the off-diagonal self-energy increment must vanish” at [S11c_c2_SHARED_PHYSICS.md:321](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:321). Treating this as a wording/representation defect is reasonable; calling the literal increment zero would not be.

### G — directionality is supported; “dissipative” is not

The evidence establishes that the increment is one-way while the full closed coupling kernel contains both directions: [S11c_c2_adjudication_verify_EG.out:1](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out:1). The spec requires an adjointness residual only when there is a genuine independent route: [S11c_c2_SHARED_PHYSICS.md:187](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:187). Grok attributed an adjoint-by-construction claim the builder did not make.

However, neither the directionality check nor builder line 70 establishes dissipativity or passivity. Remove “dissipative directional” from [STATUS.md:6](/var/projects/toy_physics/STATUS.md:6) and “a dissipative directional self-energy” from the next prompt. Use “directional/one-way SymPy increment” only.

### E — not defensibly adjudicated as a false positive

This is the major physics overstatement.

The computation proves only:

- the full residual is nonzero;
- setting `σ_W=0` kills it.

See [S11c_c2_adjudication_verify_EG.out:15](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out:15). The actual retained grades are `O(εσ_W)` and `O(εησ_W)`: [S11c_c2_sympy_stdout_index.json:1205](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json:1205).

That does not establish N6. The governing spec says the residual must vanish after field redefinition: [S11c_c2_SHARED_PHYSICS.md:303](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:303). It also inherits `η` and `σ_W` as independent first-shape-order bookkeepers; `O(εσ_W)` is retained physics, not higher order: [S11c_c1_SHARED_PHYSICS.md:223](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c1_SHARED_PHYSICS.md:223).

The adjudication’s association with the drain-projection `O(σ_W²)` deferral is especially unsupported: a residual linear in `σ_W` is not that second-order correction. The only explicit ≥64 GB deferral in the c2 scope is the full cross-engine/four-giant-family work, not this already-emitted single-engine N6 residual: [S11c_c2_SHARED_PHYSICS.md:39](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:39).

Correct disposition: **N6 is established only after projecting to `σ_W=0`; the retained σ-sector is unresolved and may be a `representation_pullback` build defect.** Grok’s E finding is not disproven.

Therefore “physics SOUND, 0 confirmed defects” and “F/G/E ALL FALSE POSITIVES” in [STATUS.md:6](/var/projects/toy_physics/STATUS.md:6) are overstated.

### Export re-review — technical evidence good, process not clear

The artifact itself looks technically sound. But [S11c_c2_export_repair_rereview_adjudication.md:17](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md:17) explicitly calls the fresh-Claude leg incomplete. It never produced a usable full report or CLEAR verdict. The record substitutes the orchestrator’s verification and says the two-leg discipline was met “in substance” at [line 45](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_adjudication.md:45).

That contradicts repository policy:

- both reports must be usable: [CLAUDE.md:46](/var/projects/toy_physics/CLAUDE.md:46);
- orchestrator verification is never a leg: [CLAUDE.md:119](/var/projects/toy_physics/CLAUDE.md:119);
- no commit before both legs report: [CLAUDE.md:137](/var/projects/toy_physics/CLAUDE.md:137).

Thus “2 re-review legs” and “re-reviewed clear” are procedural overclaims. A replacement full fresh-Claude publication review is owed.

## 2. Carry-forward set

| Item | Correct classification |
|---|---|
| F | Correct carry: genuine closed coupling decouples; literal increment does not. Fix §3c/§5e wording. |
| G | Correct measured SymPy feature, pending the blind-WL comparison. Do not call it dissipative. |
| Two S11c-b signs | Correct comparator carry; §3c and §7 explicitly say surface, do not normalize: [§3c](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:203), [§7](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:394). |
| E | Misclassified. It is an unresolved retained-order N6 finding/potential build defect, not merely a deferred record note. |

If the list is intended to be exhaustive for the step record, it is also incomplete. The blind WL/comparator still owes the §3d questions: density field-vs-field, traction representation, DtN whole-form versus kernel, traction/slab power pairing, flat-symbol usage, and the `μ_R,bg` form control: [S11c_c2_SHARED_PHYSICS.md:224](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:224). The prompt’s generic “§3d representation questions” partly covers these, but should name them.

Also preserve the fact that c1 ENERGY remains UNDECIDED. The c2 spec says it enters load-bearingly at [S11c_c2_SHARED_PHYSICS.md:84](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:84), while the builder explicitly says the c2 traction pairing does not adjudicate the separate c1 far-field ENERGY audit at [S11c_c2_sympy_builder_report.md:70](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_builder_report.md:70). At minimum the step record must retain that deferral and avoid claiming energy closure.

## 3. Hand-back prompt and compaction readiness

The post-gate sequence—blind WL directive, two decision legs before build, serialized fresh-Claude/Grok build legs, T7 comparator/reconcile, then two record legs—is otherwise correct. The model-role policy and Mathematica/setsid/same-process notes are consistent with the repo.

But the prompt is not ready:

1. Its STATE repeats the unsupported “SOUND,” “all E false-positive,” and “re-reviewed clear” claims at [next_prompt_c2_wolfram.md:3](/tmp/claude-1000/-var-projects-toy-physics/53620ffb-59f9-482d-b804-aef04f767516/scratchpad/next_prompt_c2_wolfram.md:3).

2. It skips the N6 resolution and missing fresh-Claude export-review leg.

3. Its STEP 1 would edit the shared spec and then go directly to WL. That makes the committed export provenance stale: `S11c_c2_exports.py` pins the current spec digest at [S11c_c2_exports.py:9](/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_exports.py:9), as required by [S11c_c2_SHARED_PHYSICS.md:386](/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md:386). After any spec edit, the export must be regenerated or lawfully repinned and reverified before advancing.

4. The live `project_s11c_c_state` memory repeats the same E and export-clear overclaims in its frontmatter and closeout: [project_s11c_c_state.md:3](/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:3), [line 201](/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md:201).

5. The working tree is not clean. The only dirty path is the untracked [S11c_c2_wolfram_compact_prep_verify_prompt.md](/var/projects/toy_physics/research/pde_ledger_v3/directives/_legs/S11c_c2_wolfram_compact_prep_verify_prompt.md:1). It is this audit prompt, not an uncommitted build product, but if it is intended as review evidence it must be committed with the resulting record or intentionally discarded.

6. The 498,811,405-byte reviewed SymPy stdout exists only at `/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`, and both committed adjudication scripts depend on uncommitted `/tmp/c2_ablate_agent/load_tag.py`: [verify_F.py:12](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.py:12), [verify_EG.py:2](/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.py:2). That is non-durable load-bearing evidence under the repository’s annex policy: [CLAUDE.md:176](/var/projects/toy_physics/CLAUDE.md:176).

### Precise STATUS correction

Replace “SOUND, 0 confirmed defects / ALL FALSE POSITIVES” with:

> The c2 fold wiring and A/C/D constructions are supported. F is resolved narrowly: the genuine closed coupling decouples, while the literal §3c increment retains the open-slot term and requires a §3c/§5e wording correction. G establishes a directional SymPy increment; no dissipativity claim is made. E remains OPEN: the N6 residual is nonzero at retained `O(εσ_W)`/`O(εησ_W)` grades; `σ_W→0` proves confinement only, not full representation invariance. Resolve as a possible `representation_pullback` defect or through a separately reviewed scope/spec decision before the WL build.

Replace “2 re-review legs / re-reviewed clear” with:

> Export evidence currently supports the repair: Grok CLEAR plus same-process orchestrator verification found no technical defect. The fresh-Claude attempt was incomplete; because orchestrator verification is not a leg, one usable fresh-Claude re-review remains owed before procedural clearance.

### Precise next-prompt sequence

Prepend before the current STEP 1:

> **STEP 0 — close the open gates.** Resolve E/N6 without presuming deferral: determine whether the retained σ-sector residual is a `representation_pullback` defect or requires a physics-bearing scope/spec amendment; review any script/spec change until clear. Obtain a complete usable fresh-Claude publication-review report. Do not launch the WL builder while either gate is open.

Change STEP 1’s ending to:

> After the reviewed spec correction, regenerate or lawfully repin `S11c_c2_exports.py` so `BUILD_INPUT_DIGESTS` matches the accepted spec/script bytes; verify unchanged root membership and semantics, and complete the required publication re-review before committing the refreshed delta.

Add to the carry paragraph:

> The SymPy §3d dispositions remain pending blind-WL comparison, including density, traction representation, DtN whole-form, traction/slab pairing, flat-symbol usage, and `μ_R,bg` form control. c1 ENERGY remains explicitly UNDECIDED/deferred. Do not preload these results into the blind WL packet.

Remove both occurrences of “dissipative.”

**NOT CLEAR TO COMPACT — fix first:**

1. Correct E/N6 from “false positive/deferred” to OPEN and resolve it before WL.
2. Obtain one complete usable fresh-Claude export-repair review leg.
3. Correct `STATUS.md`, the next prompt, and `project_s11c_c_state`.
4. Add the mandatory export regeneration/repin step after the spec edit.
5. Preserve the 499 MB stdout and its loader durably, or make the committed verification scripts self-contained.
6. Commit the untracked audit prompt with its result, or explicitly discard it as non-record scratch.
