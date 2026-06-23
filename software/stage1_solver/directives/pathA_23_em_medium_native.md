# Directive pathA_23 — Light as BRANE in-plane SHEAR: a candidate medium-native EM re-founding (brane-elastic fork)

**Status:** DRAFT v5 (GLM tertiary fixes applied to v4). v3 = re-pointed from the scalar-bulk framing; v4 = Codex's 10 fixes;
v5 = GLM's 3 required fixes + folded concerns. For Codex **confirm-pass on v5** → (if SOUND) staged execute. Supersedes v2.
Review record: `_scratch/codex_pathA_23_v3_review.md` (10 fixes), `_scratch/codex_pathA_23_v4_confirmpass.md` (SOUND-AS-IS),
GLM tertiary (NOT-SOUND → 3 required fixes + 10 concerns, all incorporated here; prompt `_scratch/pathA_23_v4_GLM_tertiary_prompt.md`).
**Date:** 2026-06-22.
**Owner of this doc:** orchestrator (Claude). Codex confirm-passes it, then executes stage-by-stage AFTER sign-off.
**Resume/context:** conceptual foundation `decisions/15` (picture; MacCullagh template §11 incl. the specific gauge obstruction;
routing §12; λγ subtlety §13; **honest conceptual costs §14**); program state `decisions/13` §0; value map `decisions/14`;
front door `STATUS.md`.
**Engine:** **Mathematica leads** (symbolic constitutive/stress/couple-stress tensors, principal symbols, polarization
eigenvectors, dispersion relations, the MacCullagh↔Maxwell map); **SymPy cross-checks** wherever a second engine can verify
([[feedback-dual-engine-required]]; `feedback-mathematica-single-seat`: ≤2 concurrent `math -script`).

---

## Decision (user, 2026-06-22) + honest status of the claim

The EM sector is **re-pointed to a candidate medium-native fork** in which the resolution lives in the **brane**, not the
scalar bulk. The accepted near-theorem (`decisions/15` §2, GLM): a single-component **scalar** superfluid cannot carry
transverse light. Adding **bulk** shear is disfavored — it would threaten the Magnus/magnetic force that already works
(`decisions/15` §3). The working **hypothesis** (`decisions/15` §5), to be TESTED, not assumed:

- **Our 3D space = the brane = a taut ELASTIC membrane** embedded in the 4+1D bulk.
- **LIGHT = the brane's in-plane SHEAR waves** (candidate 2 transverse polarizations, non-dispersive at `c=√(μ_brane/ρ_brane)`).
- **Shear confined to the brane, not the bulk** → the bulk would remain a pure shear-free fluid, leaving gravity (bulk flow)
  and magnetism (bulk swirl/Magnus) intact **IF there is no leak** — which Stages 1 & 3 must TEST, not presume.

**THE TEMPLATE (a hypothesis to test, `decisions/15` §11).** **MacCullagh's rotationally-elastic ether (1839)** —
`U=½μ(∇×u)²` — maps onto Maxwell *in that template* (FitzGerald/Kelvin/Larmor), transverse only; Kelvin's gyrostatic ether
realizes it via spinning substructure; Cosserat/micropolar carries an EM-like gauge structure. **But curl-only *potential*
energy alone does NOT prove gauge redundancy** — and there is a SPECIFIC known obstruction (GLM C5, `decisions/15` §11): the
**kinetic** term `½ρ(∂_t u)²` is not invariant under time-dependent `u→u+∇χ`, so the longitudinal mode is a *constrained zero
mode, not gauge*; Maxwell needs a scalar potential `φ` (`E=-∂_tA-∇φ`) to compensate, which the bare MacCullagh ether lacks.
Whether the MacCullagh→Maxwell map survives for our brane action — via a scalar-potential analog OR a constraint — is **exactly
what Stages 2–5 test (it is NOT assumed).**

**HONEST CLASSIFICATION (load-bearing).** This **adds a brane-elastic sector** to the parent action — a **`NEW_PARENT_ACTION`**,
NOT a variational rewrite. Do not call it "derived." Three provenance disciplines (GLM FIX 1, FIX 2; `decisions/15` §14):
- **Constitutive form:** either **derived from an *independently-motivated* microstructure** (motivated by the medium's known
  properties, NOT reverse-engineered to give the MacCullagh form — that is circular), or **`POSTULATED`** from symmetry.
- **CONDITIONAL-VERDICT RULE (FIX 1, deepest point):** if the form is `POSTULATED`, a full stage-pass means only **"Maxwell
  structure follows from a POSTULATED rotational elasticity"** — it does **NOT** establish that *the medium* produces it. That
  result is **`CONDITIONAL`**, is explicitly weaker than a derivation, and is **NOT sufficient on its own** to reclassify `λγ`
  in `decisions/14` or to integrate into the papers without the user's explicit acceptance of the conditional status.
- **Brane↔bulk coupling (FIX 2):** the coupling term carries the **same** provenance label (derived / symmetry-allowed-postulated
  / ad hoc); a *postulated* coupling weakens the Stage-6a charge test (it could be tuned to give the desired charge).

**Honest conceptual cost to carry (GLM C1/C9, `decisions/15` §14).** A GP/NLS superfluid is a fluid (zero shear modulus), so
the brane's shear rigidity is **not** a consequence of the mean-field — it requires the **substructure** beneath it. Honest
framing: "GP/NLS as the effective medium + a deeper substructure supplying the brane elasticity," and the construction is
bulk+brane+coupling. `FAIL_NO_TRANSVERSE_STIFFNESS` is the *most likely* Stage-2 failure. The treatment is **classical**;
quantization is deferred (C10).

**Falsification stance (user, load-bearing — `feedback-falsification-is-the-goal`).** A genuine TEST. Plausible deaths, in
priority order: **(D1) leak — the brane sector perturbs the bulk Magnus/gravity** (Stages 1, 3); **(D2) the substructure gives
Cauchy / no transverse stiffness / postulate-only elasticity, or fails energy/angular-momentum closure** (Stage 2); **(D3) not
2 transverse / dispersive / non-universal `c` / a residual longitudinal mode / a massless `u_w` fifth force** (Stage 4);
**(D4) fails Maxwell structure / gauge (the C5 obstruction) / E↔B / Ward identity / preferred-frame** (Stage 5); **(D5) charge
inconsistent with the `η_Q` ontology / energy unbounded / cone-lock needs tuning** (Stage 6). Any terminal FAIL is a
**first-class result that FALSIFIES (this construction of) the concept** — reported plainly, NEVER rescued/softened/worked
around. A clean "it all works" is suspicious and re-checked HARDER. Each stage distinguishes a **concept-fatal** FAIL from a
**needs-ingredient-X** FAIL (itself classified principled-vs-ad-hoc).

---

## §0. Why this fork — the chain (compressed; full record `decisions/15` §1)

1. GR-quadrupole verdict needs `λγ=c_γ/c_s` (`BETA_GENUINE_GAP`); EM-anchor scoping → only a **speed-cone** observable isolates it.
2. `c_γ(ρ0)` derivation → current action is **Type-4 DECOUPLED** (`c_γ`⊥`ρ0`, `c_s∝ρ0²`) ⇒ cones can't lock; fingerprint of a
   unified-gauge-on-flat-metric **drift** from the single-medium concept (`reports/pathA_cgamma_of_rho_derivation.md`).
3. Re-founding medium-native exposed the scalar near-theorem ⇒ light cannot come from the scalar bulk.
4. Candidate re-pointing: light = **brane in-plane shear**, with **MacCullagh rotational elasticity** as the template to test
   and **Cosserat/spinning-substructure** as a candidate mechanism (subject to the independent-motivation rule above).

**A consequence the fork INTRODUCES (`decisions/15` §13).** `c_γ=√(μ_brane/ρ_brane)` is a **brane** property; `c_s` a **bulk**
property ⇒ **`λγ=1` is not automatic** — a no-tuning relation between two distinct elastic sectors. Because shear and compression
speeds differ generically, **`CONE_LOCK_ABSENT`/`CONE_SCALING_ONLY` is the EXPECTED outcome** (GLM C7); a clean lock is
suspicious, not a triumph.

---

## §1. Goal + acceptance bars

**Goal.** Specify a candidate **brane-elastic** EM sector (light = in-plane shear), classify it honestly (`NEW_PARENT_ACTION`;
constitutive form derived-vs-postulated → conditional-verdict rule; coupling provenance), and **TEST** whether it (a) leaves the
working bulk gravity/magnetism undisturbed, (b) yields genuine electromagnetism (2 transverse Maxwell modes + gauge + E↔B + a
gauge-invariant conserved source), and (c) admits a characterized, no-tuning cone relation feeding `λγ`.

**Acceptance bars (ALL required; a clean pass on ANY is adversarially re-checked HARDER).**
0. **Provenance + admissibility.** Constitutive form AND brane↔bulk coupling each labeled derived(independently-motivated) /
   postulated / ad hoc; the conditional-verdict rule applied. A cheap **current-admissibility pre-check** (does the defect-brane
   coupling structure admit a gauge-invariant conserved current AT ALL? GLM C8) run early — fail early if not. *(FIX 1, FIX 2, C8.)*
1. **No leak — kinematic + constitutive.** (i) constitutive-independent kinematic audit over ALL candidate forms incl. `u_w`
   (Stage 1); (ii) post-constitutive traction/couple-stress closure comparing any induced bulk transverse source to the
   Magnus/gravity terms with units restored + a named small parameter (Stage 3). The coupling must source the
   **longitudinal/static (Coulomb) sector WITHOUT sourcing transverse shear** — the explicit decomposition that lets no-leak and
   charge-sourcing coexist (GLM C2). *(D1.)*
2. **Constitutive provenance + closure.** Brane elastic energy derived(independently-motivated) or postulated, with a
   per-ingredient provenance table; plus stress AND couple-stress tensors, boundary traction/work, **total angular-momentum
   conservation incl. spin/couple**, and Hamiltonian **positivity for all Fourier modes**. Any Cosserat micro-rotation **gap
   scale is itself a potential tuning** — classify it in the provenance framework (GLM C6). *(D2.)*
3. **DOF + spectrum.** Two transverse propagating polarizations, **non-dispersive**, single universal `c`; the longitudinal
   sector shown gauge OR an honest residual; the off-brane bending `u_w` shown **stable and gapped/confined-to-`w`/constrained**
   — a **massless unsuppressed `u_w` is EXCLUDED (fifth force)** (GLM FIX 3). *(D3.)*
4. **Maxwell structure.** Homogeneous AND inhomogeneous Maxwell from the brane action with a **generic conserved `J_brane` +
   Ward identity**; explicit gauge redundancy **confronting the C5 obstruction** (test whether a scalar-potential analog `φ` or
   a constraint restores full Maxwell gauge); dual dictionary `u↔(A,E,B)` pinned by derivation. *(D4.)*
5. **Lorentz + E↔B.** Boost law of `E`,`B` from the brane variables (not by assuming `A^μ` is a four-vector), reproducing EM
   with `c→c_γ`; exact-vs-leading scope; **preferred-frame effects bounded from BOTH the `w`-anisotropy AND the bulk FLOW `v_r`**
   (GLM C3). *(D4.)*
6. **Charge.** `J_brane` derived from the defect/throat coupling; `∂_μ J_brane^μ=0`; Gauss flux → charge; sign/quantization/
   normalization matched to the topological `η_Q e_*` ontology (`pde.tex:279-312`) — NOT circulation sign — `q∝ρ0ΓA`
   reinterpreted/rejected under the charge firewall; no source double-counting; a postulated coupling weakens this test (FIX 2). *(D5.)*
7. **Energy/consistency.** Hamiltonian bounded below (confront the MacCullagh antisymmetric-stress/ghost concern); no
   double-counting of medium energy bulk↔brane; charge conservation. *(D5.)*
8. **Cone payoff — no tuning.** Both `c_γ²=μ_brane/ρ_brane` and bulk `c_s²` from the **same** microparameters (brane density
   dimension; `w`-Jacobian); classified `CONE_LOCK_DERIVED` / `CONE_SCALING_ONLY` / `CONE_LOCK_TUNED_POSTULATE` /
   `CONE_LOCK_ABSENT` (the EXPECTED outcome, C7); throat-local deviations characterized. *(D5.)*

---

## §2. Staged plan (ONE stage at a time; STOP + explicit user gate between stages — `feedback-sequential-audit-chunks`)

Codex designs each stage's derivation route (REQUIREMENTS + acceptance only here). Each stage gets a verdict-first report under
`reports/` + a post-exec **fidelity** (`feedback-transliteration-fidelity-audit`) and **adversarial** review (clean agents)
before the gate. No roll-forward without a user gate.

### Stage 0 — Action + embedding contract + microstructure contract + coupling provenance + early admissibility + DOF + classification
- **Action:** bulk scalar sector (KEEP) + **brane-elastic sector** (in-plane `u_∥` + off-brane bending `u_w`) + the
  **brane↔bulk coupling** term.
- **Embedding/symmetry contract (fix #8 + GLM C4):** brane coordinates, induced metric, normal vector, extrinsic-curvature
  assumptions; whether `u_∥` is physical displacement or material relabeling; **whether the brane is thin (δ-function in `w`) or
  finite-thickness, and how the shear modulus is distributed in `w`**; BCs at defects/throats; residual symmetries; how the bulk
  flow intersects the brane; how the old `Z(w)`/`W(w)` projection firewall is retained/replaced.
- **Microstructure contract (fix #1 + GLM FIX 1):** EITHER specify a minimal substructure/cohesion model (symmetry class,
  internal variables, coarse-graining rule, allowed-invariant list) that is **independently motivated by the medium's known
  properties** (NOT reverse-engineered to yield MacCullagh — flag circularity), from which Stage 2 can *derive* the form; OR
  explicitly declare the form will be *postulated* (→ the conditional-verdict rule applies).
- **Coupling provenance (GLM FIX 2):** label the brane↔bulk coupling derived / symmetry-allowed-postulated / ad hoc.
- **Current-admissibility pre-check (GLM C8):** does the defect-brane coupling structure admit a gauge-invariant conserved
  current at all? If structurally impossible, FAIL EARLY (before investing in Stages 2–5).
- Enumerate candidate constitutive forms: (i) Cauchy/Navier; (ii) rotational/MacCullagh (curl-only); (iii) Cosserat/micropolar.
- **Classify the action** `NEW_PARENT_ACTION`; DOF count of the brane EM content vs the fundamental-`A_M` action (`pde.tex:257-261`).
- Outcomes: `ACTION_SPECIFIED_CLASSIFIED` / `FAIL_ILLDEFINED_BRANE_EMBEDDING` / `FAIL_NO_CONSISTENT_COUPLING` /
  `FAIL_UNSPECIFIED_SUBSTRUCTURE` / `FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED` / `FAIL_NO_ADMISSIBLE_CURRENT` /
  `FAIL_PREFERRED_W_ANISOTROPY_AT_PRINCIPAL_LEVEL` / `FAIL_PROJECTION_REDUCTION_FIREWALL`.
- Acceptance: action + both contracts + coupling provenance + admissibility pre-check + DOF + classification. **STOP.**

### Stage 1 — Kinematic / coupling leak audit (cheap, pre-constitutive; first death-check, D1)
- From the Stage-0 coupling *structure* alone, test whether a bulk leak channel exists **for every candidate constitutive form**
  (a kinematic theorem) — into bulk shear, the Magnus force, or the gravity flow `v_r`. **Include `u_w`** (can bending mediate a
  kinematic bulk leak?). Test the **sector decomposition (GLM C2):** the coupling should be able to source the
  longitudinal/static (Coulomb) sector WITHOUT sourcing transverse shear; flag if it cannot. It must NOT pass by zeroing the
  coupling (Stage 6 needs bulk drains to source charge).
- Outcomes: `NO_KINEMATIC_LEAK_FOR_ALL_CANDIDATES` / `FAIL_COUPLING_LEAKS_INDEPENDENT_OF_CONSTITUTIVE` /
  `FAIL_COUPLING_SOURCES_TRANSVERSE` / `LEAK_CONDITIONS_DEFERRED` (to the throat-profile + Stage-3 closure) / `FAIL_BENDING_BULK_LEAK`.
- Acceptance: the kinematic leak audit (incl. `u_w` + sector decomposition), constitutive-dependence deferred explicitly. **STOP.**

### Stage 2 — Constitutive form from the Stage-0 microstructure (THE CRUX, D2)
- **USER DECISION (2026-06-22) — track B (stronger test):** Stage 2 MUST make a genuine, bounded attempt to **independently
  motivate** a microstructure from the medium's known properties and **derive** the constitutive form — falling back to
  `POSTULATED` (→ CONDITIONAL verdict) ONLY if it honestly cannot, with the reason recorded. Do not default to postulate to
  save effort; but do not reverse-engineer a microstructure to yield MacCullagh either (`FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED`).
- If a microstructure was specified: **derive** the brane elastic energy and classify (rotational vs Cauchy vs Cosserat),
  confirming it was **independently motivated** (else `FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED`). If postulated: state it, grade
  `ROTATIONAL_POSTULATED_NOT_DERIVED`, and mark the downstream verdict `CONDITIONAL`. **Provenance table required.**
- **Operationalize the obstructions (fix #4):** stress AND couple-stress tensors from the action, boundary work/traction,
  **total angular-momentum conservation incl. spin/couple**, **Hamiltonian positivity for all Fourier modes**. For Cosserat:
  explicit spectrum + coupling classification of **every** micro-rotation mode, and **classify the gap scale as a potential
  tuning** (GLM C6).
- Confront head-on: GP/NLS gives zero shear ⇒ `FAIL_NO_TRANSVERSE_STIFFNESS` is the *most likely* failure; or Cauchy ⇒ a
  propagating (faster) longitudinal mode.
- Outcomes: `ROTATIONAL_CURL_ONLY_CANDIDATE` (NOT a Maxwell verdict) / `COSSERAT_MICROPOLAR(spectrum,gap)` /
  `ROTATIONAL_POSTULATED_NOT_DERIVED(→CONDITIONAL)` / `MICROSTRUCTURE_DERIVES_CAUCHY` / `FAIL_NO_TRANSVERSE_STIFFNESS` /
  `FAIL_CAUCHY_STRAY_LONGITUDINAL` / `FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED` / `FAIL_ANTISYMMETRIC_STRESS_NO_SPIN_CLOSURE` /
  `FAIL_ENERGY_GHOST` / `FAIL_EXTRA_MICROPOLAR_LIGHT_MODE`.
- Acceptance: constitutive form + provenance table + stress/couple-stress + angular-momentum + positivity + gap classification. **STOP.**

### Stage 3 — Constitutive no-leak closure (the real Magnus-break gate, D1)
- With the Stage-2 form + coefficients fixed (fix #2): vary the full action w.r.t. bulk density/phase/velocity and the brane
  embedding; derive brane boundary/jump conditions; compute normal+tangential **traction** and **couple-stress**; determine
  whether the bulk Euler/vorticity equations acquire a **transverse source**; compare its order to Magnus/gravity-flow with
  **units restored + a named small parameter**. Confirm the C2 decomposition survives at the constitutive level.
- Outcomes: `NO_LEAK_CLOSED(order)` / `FAIL_CONSTITUTIVE_TRACTION_LEAK` (= `FAIL_LEAK_BREAKS_MAGNUS`, concept-fatal) /
  `LEAK_BOUNDED_CONDITIONAL(condition+price)`.
- Acceptance: the constitutive leak closure with order/condition explicit. **STOP.**

### Stage 4 — Spectrum, polarizations, dispersion, `u_w` (make-or-break, D3)
- From the Stage-2 law: derive the **principal symbol**; **count propagating modes**; require **two transverse polarizations**,
  **non-dispersive** (`ω=ck`), single universal `c`. Show the longitudinal sector is gauge OR flag a residual physical mode.
- **`u_w` audit:** stable? **gapped / confined-to-`w` / constrained?** mixing with in-plane shear? extra scalar radiation
  channel? **A massless unsuppressed `u_w` in our 3D space is a fifth force → FAIL** (GLM FIX 3).
- Outcomes: `TWO_TRANSVERSE_NONDISP_UNIVERSAL` / `FAIL_DISPERSIVE` / `FAIL_WRONG_POLARIZATION_COUNT` / `FAIL_NONUNIVERSAL_C` /
  `FAIL_RESIDUAL_LONGITUDINAL_ZERO_MODE` / `FAIL_BENDING_INSTABILITY` / `FAIL_BENDING_SHEAR_MIXING` /
  `FAIL_BENDING_MASSLESS_FIFTH_FORCE`.
- Acceptance: mode content + speeds + dispersion + `u_w` audit, honest pass/fail. **STOP.**

### Stage 5 — Maxwell structure: dictionary, gauge (the C5 obstruction), field equations, Lorentz/E↔B (D4)
- Pin the dictionary `u↔(A,E,B)` **by derivation**. Derive **homogeneous** AND **inhomogeneous** Maxwell from the brane action
  with a **generic conserved `J_brane`** (physical source derived in Stage 6) + the **Ward identity**.
- **Confront the C5 gauge obstruction explicitly:** the bare curl-only kinetic term breaks `u→u+∇χ`; test whether the brane
  theory supplies a **scalar-potential analog `φ`** (`E=-∂_tA-∇φ`, restoring full Maxwell gauge) OR a **constraint** removing
  the longitudinal mode — state which is expected and which the action actually realizes.
- Derive the **Lorentz** transformation of the brane perturbations (invariant `c_γ`) + **E↔B mixing**; state exact-vs-leading
  scope; **bound preferred-frame effects from BOTH the `w`-anisotropy AND the bulk flow `v_r`** (GLM C3).
- Outcomes: `MAXWELL_STRUCTURE_RECOVERED(via φ|via constraint)` / `FAIL_CURL_ONLY_NOT_GAUGE` / `FAIL_NO_GAUGE` /
  `FAIL_NOT_MAXWELL` / `FAIL_SOURCE_BREAKS_GAUGE` / `APPROX_EB_MIXING(scope)` / `FAIL_EB_MIXING` / `FAIL_PREFERRED_FRAME_FLOW`.
- Acceptance: dictionary + derived Maxwell (generic `J`) + Ward identity + gauge resolution + E↔B + preferred-frame bounds. **STOP.**

### Stage 6 — Charge, energy bookkeeping, cone payoff, leftovers (split sub-verdicts, fix #7; D5 + λγ)
Each sub-verdict must be **able to fail independently**; a concept-fatal sub-FAIL stops the process.
- **6a — Charge (D5).** Derive `J_brane` from the defect/throat coupling; prove `∂_μ J_brane^μ=0`; Gauss flux → charge; match
  sign/quantization/normalization to the topological `η_Q e_*` ontology (NOT circulation sign); reinterpret/reject `q∝ρ0ΓA`
  under the charge firewall; no source double-counting. Note: a *postulated* coupling (Stage 0) weakens this (FIX 2). Outcomes:
  `CHARGE_RECONCILED` / `FAIL_CHARGE_SIGN_FIREWALL` / `FAIL_SOURCE_NOT_GAUGE_INVARIANT` / `FAIL_NONQUANTIZED_CHARGE`.
- **6b — Energy bookkeeping + EM-carrier reconciliation (D5).** Hamiltonian bounded below (incl. MacCullagh concern); no
  double-counting bulk↔brane. **Reconcile the two EM carriers:** the frozen `A_M` Maxwell sector still lives in `pde.tex`
  alongside the new brane-elastic EM sector, so a gate must **remove, identify, or demote one carrier** — and if neither can be
  demoted, the construction carries TWO photons (concept-fatal, not merely "to be avoided"). Outcomes: `ENERGY_CONSISTENT` /
  `FAIL_ENERGY_UNBOUNDED` / `FAIL_ENERGY_DOUBLE_COUNT` / `FAIL_DOUBLE_COUNT_IRREDUCIBLE` (cannot demote either EM carrier).
- **6c — Cone payoff (λγ, no tuning).** Express `c_γ²=μ_brane/ρ_brane` and `c_s²` from the same microparameters (units restored;
  brane density dimension; `w`-Jacobian). **`CONE_LOCK_ABSENT`/`CONE_SCALING_ONLY` is the EXPECTED result** (C7); a clean lock
  is re-checked HARDER. Outcomes: `CONE_LOCK_DERIVED` / `CONE_SCALING_ONLY` / `CONE_LOCK_TUNED_POSTULATE` / `CONE_LOCK_ABSENT` /
  `FAIL_DIMENSIONAL_BRANE_DENSITY`.
- **6d — Leftovers (`decisions/15` §7a).** Catalog residual modes as candidate features (cf. `S_leak`) — but only if
  gapped/confined (a massless one already failed in Stage 4).
- **Do NOT reclassify `λγ` in `decisions/14`, do NOT relax the EM-anchor bookkeeping, and do NOT integrate papers until ALL
  gates pass AND (if the constitutive form is POSTULATED) the user explicitly accepts the CONDITIONAL status** (FIX 1). Only
  then: orchestrator updates `decisions/14` + `decisions/13` §0 + `STATUS.md`; then (gated, last) paper integration (Codex
  applies edits, Claude reviews — `feedback-codex-is-fix-applier`).
- Acceptance: the four sub-verdicts (6a–6d) with honest provenance + conditional-status flag. **STOP.**

---

## §3. Discipline (binding on every stage)

- **TARGET-BLIND / ABLE-TO-FAIL.** Every stage has concrete reachable FAIL outcomes. A clean success is re-checked HARDER; a
  terminal FAIL falsifies (this construction of) the concept and is never softened/rescued (`feedback-falsification-is-the-goal`).
  Guard actively against: *choosing the desired constitutive law and interpreting its null longitudinal sector as gauge without
  proving the full constrained, sourced theory* — and against *bypassing the central crux by postulation without flagging the
  verdict CONDITIONAL.*
- **NO TAUTOLOGY / NO CIRCULAR MICROSTRUCTURE.** A constitutive "derivation" from an unspecified substructure is forbidden; a
  microstructure reverse-engineered to yield MacCullagh is circular (`FAIL_MICROSTRUCTURE_REVERSE_ENGINEERED`). A gate whose
  positive result is a necessary consequence of its setup decides nothing (`feedback-decisive-test-not-tautological`).
- **HONEST PROVENANCE + CONDITIONAL-VERDICT RULE.** Action = `NEW_PARENT_ACTION`; constitutive form & coupling each labeled
  derived(independently-motivated)/postulated/ad hoc; a postulated form ⇒ CONDITIONAL verdict, insufficient alone for
  `decisions/14`/papers without explicit user acceptance. Dual dictionary pinned by derivation.
- **NO-LEAK BEFORE PAYOFF, IN TWO PASSES** (kinematic Stage 1, constitutive Stage 3); cheapest concept-fatal checks first.
- **CONFRONT THE KNOWN OBSTRUCTIONS** (not name-check): the C5 kinetic-term gauge obstruction (Stage 5); antisymmetric stress /
  spin-couple closure / energy positivity / micropolar gap-tuning (Stage 2); residual longitudinal & massless-`u_w` fifth force
  (Stage 4); preferred-frame from `w` AND bulk flow (Stages 0, 5); charge firewall (6a). The brane elasticity needs substructure
  beyond GP/NLS (carry the conceptual cost, `decisions/15` §14).
- **DIMENSIONAL-CONSISTENCY with units RESTORED** (no hidden brane-volume/`w`-Jacobian) + **DUAL-ENGINE** symbolic checks
  (Mathematica primary, SymPy cross-check; `feedback-transliteration-fidelity-audit` after every computational step).
- **No leading / conclusion-shaped language** ("forced", "by construction", "dissolves", "preserved", "reproduces Maxwell
  exactly", "resolution") in stage reasoning.
- **Codex designs the route + writes/runs all scripts; Claude reviews only** (`feedback-claude-reviews-codex-codes`). Scripts
  ≤10 min with `timeout 600` on the scripts (never on the Codex session). Classical scope; quantization deferred (C10).
- **No paper / `decisions/*` edits during Stages 0–6.** Reports → `reports/`; orchestrator owns the canonical docs.

## §4. Deliverables
- Per stage: `reports/pathA_23_stageN_*.md` (verdict-first) + Mathematica `.wl` (+ SymPy cross-check) under the solver tree.
- Final (gated, only after full pass + conditional-status acceptance if applicable): the result with honest `NEW_PARENT_ACTION`
  + constitutive/coupling-provenance + conditional flag; updated `decisions/14` / `decisions/13` §0 / `STATUS.md`; limitations
  ledger (incl. the §14 conceptual costs); (last) paper integration.

## §5. Review plan
1. **Codex confirm-pass** on this v5 (`-c model_reasoning_effort=xhigh`; confirm log shows `reasoning effort: xhigh`): are the
   3 GLM required fixes (conditional-verdict rule; coupling provenance; massless-`u_w` FAIL) + the folded concerns correctly
   incorporated WITHOUT regressing the v4 structure, and is the plan still SOUND (able-to-fail, non-tautological, ordered)?
2. **Three-way agreement reached** (Codex v4 SOUND + GLM v4 fixes applied + Codex v5 confirm) → **execute stage-by-stage**, each
   with a post-exec fidelity + adversarial review and a user gate before the next stage.
