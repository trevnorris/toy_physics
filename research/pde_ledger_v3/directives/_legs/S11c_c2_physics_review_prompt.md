# Independent physics review — the S11c-c2 self-energy fold (SymPy audit script)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/scripts/S11c_c2_selfenergy_fold_sympy_audit.py`
(1053 lines; a SymPy engine that folds a closed curved-bulk face response into a variable-coefficient slab
operator and re-extracts an off-diagonal self-energy). Its emitted output is the 499 MB file
`/tmp/S11c_c2_selfenergy_fold_sympy_audit.out`, navigable by byte offset via the index
`/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json`
(a 154-entry `tags` dict, each `PY_S11CC2_<QUANTITY>_<ANCHORING>_<DENSITY>` → `{byte_offset, byte_count, sha256}`).

You are ONE of two independent legs. You will not see the other's output. Derive the physics yourself; the
script is the thing under test, not a source of truth.

## Artifact role
This is a **CODE build** (SymPy). Your job is to establish, by independent computation, whether the script's
**physics** is right — the fold construction, the six re-adjudications, the controls — and to catch any way the
physics could be **wrong**. A build is reviewed **until clear**; a disagreement you find is a *finding*, not a
failure to fix. ⛔ You are NOT checking whether the script matches any instruction sheet — you are checking whether
its physics is correct against the spec's governing equations, which you re-derive from first principles.

## What you are handed (read any of these; nothing is withheld that you may use)
- The script above (the deliverable).
- **The physics authority**: `/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md`
  — the governing setup (§1 inherited models + face closure laws + Λ-routing), the fold construction (§2 ordering,
  §3 outputs), the six §3d re-adjudications, and the §5 controls. This states **no** component value, sign, order,
  or parity — deriving those is your job.
- **The real imported inputs** (so you can verify the import wiring binds the correct objects, not a false-equal):
  `.../scripts/S11c_b_exports.py`, `.../scripts/S11c_c1_exports.py`, and the fold machinery `.../scripts/ledger_fold.py`.
- **Predecessor step records for physics context** (consult as needed):
  c1 `.../steps/S11c_c1_curved_bulk_closure.md` (the closed face response `(δp_s,J_s,t_s)(V_s,μ_θ)`, the permeable
  resolvent, the w-jets); S11c-b `.../steps/S11c_b_variable_coefficient_operator.md` (the slab operator rows and
  its two UNVALIDATED sign conventions).
- The emitted `.out` + its byte-offset index (above).

## Required method — this is a SCRIPT, so DERIVE INDEPENDENTLY
1. **Write your own derivation script BEFORE you rely on the artifact's output.** Save both the script and its
   literal stdout to named absolute paths and report those paths. ⛔ A prose derivation ("I worked it out and got
   X") is **discarded** — unless it is in CAS with visible output from inputs, it does not count.
2. For every claim you make about the script, name **which line computed it** (or report it as uncomputed). Two
   hand derivations agreeing is not evidence; a runnable script and its stdout is.
3. **Ablate every load-bearing check and report its literal output.** Probe for: a value verified with the same
   predicate/definition that produced it; a conclusion emitted as an unconditional literal; a hand-typed CAS
   payload with no data dependence on the derivation (delete the derivation → does the output move?); an `assert`
   that precedes the value it guards.
4. **⛔⛔ A FORM ABLATION IS MANDATORY, NOT OPTIONAL.** Change the *structure* of at least one load-bearing object
   (flip a sign AND an off-diagonal; collapse two independent symbols into one; break the closure threading) and
   report the **literal** diff. A **coefficient** rescale tests only arithmetic; only a **FORM** change tests
   physics. This is the only instrument that has ever caught this program's worst defect (a typed sentence with no
   computation behind it, byte-identical under ablation).
5. **One-sided corruption for independence.** Where the script claims two independent routes to one object and
   emits their residual, a zero residual proves nothing until you corrupt ONE route and show only that route's
   object moves. If breaking route A also moves route B, they were never independent and the residual is decoration.

## The physics to settle (each is a computation whose disposition is a FINDING — ⛔ not a value to match)

**A. The fold map (spec §3a / §1c / §1d).**
- The fold substitutes the c1 **closed face pressure `δp_s(V_s,μ_θ)`** (from `s11c_c1_face_response`) and its
  **w-jets `d_w_delta_p_s`** into the symbolic `delta_p_±` / `d_w_delta_p_±` slots of the slab operator. Verify it
  substitutes the **pressure**, and that it does **NOT** add a closed `J_s` term (there is no `J_s` slot — the
  θ-row flux is already `Λ_A𝒜_s + Λ_V V_s`; adding a `J_s` would double-count). `build_face`@477, `build_case`@532.
- The w-jets `d_w_delta_p_s` must be **computed** (differentiated), not hand-typed. Which line differentiates them?
- The response elimination is the c1 **operator inverse `[I+(Λ_A/ρ_m²)Z]⁻¹`**, ⛔ not a scalar division. Verify.
- The **kernel bridge** `dtn_operator → dtn_kernel` (`kernel_bridge`@366): the load-bearing bulk object must be the
  cross-engine-AGREE'd two-momentum **kernel** carrying both legs `q_out(k), q_out(k′)` — not the UNDECIDED raw
  whole-form, and not a single-k multiplier (a one-leg freeze deletes mode mixing — the rule-17 trap).
- **The symbol map** (how c1's `V_s`, opaque `μ_θ`, per-case resolvent symbols land on the slab-row fields): derive
  which slab-substrate field each c1 symbol *physically* is (e.g. is `V_s` the interfacial normal velocity, i.e. the
  S11c-a `face_velocity`, or something else?), and check the script's map is the physically correct one.

**B. ε-normalization (spec §4, `N12`).** The transverse↔thickness coupling is the inherited **O(εη)**, ⛔ NOT O(ε²).
c1 builds its response inputs as `ε·symbol`, so `δp_s` already carries one ε; the slab slot carrier carries an
outer ε. A naive direct substitution would give a spurious **O(ε²)**. Verify the script strips exactly one ε so the
self-energy increment lands at O(εη). Compute the multigrade of the emitted increment yourself and compare.

**C. Close-then-extract ordering (spec §2, §5a).** The construction must be `extract(close(SLAB))`, and
`extract(close)` ≠ `close(extract)`. The emitted `S11CC2_ORDERING_COMMUTATOR = extract(close(SLAB)) −
close(extract(SLAB))` must be **genuinely nonzero in general** — a zero/byte-identical commutator is the tell that
the closure was NOT threaded through the full operator (so no self-energy was induced). FORM-ablate the closure
threading and confirm the commutator responds. `extract`@326, and the run() assembly at lines 875–906.

**D. The six §3d re-adjudications** (each a c1 UNDECIDED / rule-17 carry-in the fold makes load-bearing):
  1. **Background density, field-vs-field (rule 17, c1 seal 5).** `rho_br_bg_rho4_constant` must be bound to its
     **live** relation `W_bg·ρ_br/W_0` with `W_bg=W_0(1+η w₁)`, ⛔ NOT a bare constant — else the fold silently
     repeats the c1 PY freeze while passing the manifest guards. `S11CC2_DENSITY_LIVE_MINUS_FROZEN` must be a
     **live-field difference**, ⛔ NOT a `∇ρ` jet, ⛔ NOT `ρ_m` (the bulk fluid density). FORM-ablate: freeze the
     field and confirm the increment changes (the freeze is observable), and confirm `∇ρ→0` is *not* what's tested.
  2. **`t_s` traction representation.** Carried in its native covector form `t_s = −(δp_s + Λ_X𝒜_s) n̂_s`
     (a covector along the curved outward normal), ⛔ not pre-collapsed to a scalar. `traction_pairing`@772.
  3. **DtN kernel vs whole-form.** `S11CC2_DTN_WHOLEFORM_DEPENDENCE` must separately expose which self-energy terms
     depend on the **raw whole-form `dtn_operator`** beyond the AGREE'd kernel. Ablate it: does any *load-bearing*
     term secretly consume the UNDECIDED whole-form (which would inherit an unvalidated cross-engine object)?
  4. **Traction-vs-slab mechanical-power pairing.** `S11CC2_TRACTION_SLAB_POWER_PAIRING` + residual settle the
     load-bearing **face-generalized-force sign** that the §3c increment does NOT cancel. A **one-sided `t_s`-sign
     corruption** (`flip=True`, emitted as `TRACTION_SIGN_RESIDUAL`) must move the mechanical-power residual against
     the slab kinetic term. Verify that residual is a **genuine independent route**, not a structural zero — corrupt
     one side and show it moves alone.
  5. **Flat-resolvent leg-labeling.** `S11CC2_FLAT_SYMBOL_USAGE`: the PY-output-leg vs WL-input-leg convention is
     UNDECIDED, equal on the `k=k′` diagonal, differing off-diagonal in MATERIAL_ADVECTED. If the fold uses
     `dtn_flat_symbol` **only** as the uniform-limit diagonal regression operand, that's fine; if **any** MATERIAL
     off-diagonal self-energy term consumes it, that is an unre-adjudicated UNDECIDED item — a finding. Check.
  6. **`μ_R,bg` FORM ablation.** `S11CC2_MU_R_FORM_RESIDUAL` must move under a **FORM** change of `μ_R,bg(x)` in the
     composed `μ_θ` operand (a coefficient rescale is insufficient — only a form change tests the coupling).
     `modulus_form`@959, `control(...,'MU_R_FORM',...)`@926.

**E. N6 representation-invariance (spec §5c).** The two anchorings `α∈{LAB_HELD, MATERIAL_ADVECTED}` ARE the
independent-route pair: mapped Eulerian↔material by `Δρ = δρ_E + u·∇ρ⁰`, the residual `S11CC2_REP_INVARIANCE_RESIDUAL`
is predicted to **vanish** (same operator, two representations). Settle by **computation** whether the emitted
residual actually vanishes (simplify, or numeric Schwartz–Zippel on the noncommutative-symbol-freed form) — ⛔ do
NOT assume it does. Then confirm the one-sided `FLIP_FACE_SLOPE` corruption moves ONE route only. Lines 931–951.

**F. Uniform limit (spec §5e) — settle astra's self-flagged concern by computation.** The spec predicts: in the
uniform limit (`W_bg→W̄₀`, `η→0`) the off-diagonal self-energy increment must **vanish** (S11b decoupling). The
emitted `S11CC2_UNIFORM_LIMIT_RESIDUAL` is **megabytes** in each case (not a visibly-zero payload). **Determine by
computation whether it is actually zero** (a large but structurally-zero unsimplified expression) **or genuinely
nonzero** (a decoupling failure — e.g. a minus-face `delta_p_minus` term that does not cancel in the uniform limit).
This is the load-bearing open question. `control(...,'UNIFORM_LIMIT',...)`@925.

**G. Tautology / operand-theatre (rule 2 corollary 3).** The §3c increment is an **export representation**, not a
check — verify it is not dressed as a residual. The builder emitted **no** `SELF_ENERGY_ADJOINTNESS_RESIDUAL`,
claiming the two off-diagonal blocks are adjoint by construction with no independent second route — verify that is
honest (the blocks really are structurally adjoint) rather than a suppressed check. For every emitted `*_RESIDUAL`,
ask: were its two operands produced by **independent routes**? If `q:=A/B` and it emits `A−q·B`, that is zero for
any input including a wrong one.

## Physics filter
Report a finding only if it catches a way the **physics** could be wrong. ⛔ Do not report "the script would be
wrong on a different input" or style/robustness nits. If nothing survives this filter, say so — but a leg that finds
nothing is weak evidence, so show the ablations you ran and their literal output.

## Ablation sandbox and operational limits (SymPy — heavy)
- ⛔ **Never modify the working tree.** Copy what you need into a private sandbox and ablate the COPY:
  ```bash
  mkdir -p /tmp/c2_ablate/scripts /tmp/c2_ablate/_measurements
  cp /var/projects/toy_physics/research/pde_ledger_v3/scripts/{S11c_c2_selfenergy_fold_sympy_audit.py,ledger_fold.py,S11c_b_exports.py,S11c_c1_exports.py} /tmp/c2_ablate/scripts/
  # edit the COPY under /tmp/c2_ablate/scripts/, then run (ROOT resolves to /tmp/c2_ablate — writes stay there):
  cd /tmp/c2_ablate && S11CC2_PACKAGE=TRIAGE timeout 1200 python3 scripts/S11c_c2_selfenergy_fold_sympy_audit.py > /tmp/c2_ablate/run.out 2>&1
  ```
- `S11CC2_PACKAGE=TRIAGE` runs **one** case and emits only the 3 core objects (closed slab operator, coupling
  kernel, increment) — fast; use it to FORM-ablate the **core fold** (close/extract/increment/ordering). The §3d
  re-adjudications, the controls (uniform limit, μ_R form, traction sign), N6, and the export run **only in the full
  path** (no TRIAGE env), which is ~19 min / ~1.4 GB peak — budget it, run **one heavy job at a time**, and watch
  `free -h`. A re-run past `timeout 1200` (20 min) is a **failed ablation**: report it and move on, ⛔ do not raise
  the timeout.
- Save every ablation script AND its literal stdout to named absolute paths and report them.

## Output
For each of A–G: your finding + the evidence (your derivation script path + its stdout, the ablation diff, and the
script `file:line` you're judging). Separate CONFIRMED defects (physics is wrong) from questions you could not
settle. End with a one-line disposition per item and an overall: is the fold's physics sound, or what must change?
