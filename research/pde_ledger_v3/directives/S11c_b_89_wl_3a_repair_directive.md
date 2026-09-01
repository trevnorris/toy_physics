# S11c-b #89 (WL) — repair the §3a energy basis: complete the O(3)-Kronecker field-bilinear enumeration

## 0 · Role and single deliverable

Modify **in place** the Wolfram engine
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
so that its §3a new-invariant family is the **complete** set of scalar field-bilinears the S11b symmetry
group admits with the background first jets carrying their transformation — not a hand-picked subset — with
the exact local-thickness map imposed before the rank test, and re-run the engine to completion. The run
regenerates the tracked output `research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out`.
That regeneration is authorized and expected; change nothing else in the tree.

**Scope: the §3a energy BASIS only.** This engine also **freezes the background jet in its slab operator**
(see §2C) — a real, separate defect that this build does **not** fix; it is deferred to a follow-up (#89b).
This directive requires you to (a) complete the basis enumeration and impose the thickness map, and (b) emit
a **diagnostic** that documents the operator freeze — not to repair the operator. Do not claim the operator
is repaired, and do not treat an operator-level cross-engine disagreement as closed by this build.

This is a **spec-compliance repair of an existing implementation defect**, not new physics and not a spec
change. The defects, the spec obligations they violate, and the code paths are supplied below as verified
inputs. You are asked to make the enumeration complete and the map consistent — to construct **every**
symmetry-allowed field-bilinear the spec names — and let the engine's existing rank selection compute and
emit the independent count, whatever it is.

This engine is written to **re-derive** everything; it imports nothing from the SymPy engine. Build WL's
enumeration from the physics stated below, Wolfram-idiomatically. The completed count is withheld from you on
purpose (§4); the enumeration is fully determined by the physics, so there is nothing to tune toward.

## 1 · The spec obligation (SUPPLIED, verbatim — the object to satisfy)

From `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`:

- **§3a (245–256):** "enumerate the DOF fields and their first gradients `{u,∇u,θ,∇θ,e_W,∇e_W}` — one
  thickness coordinate only, **with the exact map `e_W,bg=(W_0/W_bg)e_W` imposed before any
  independence/rank test** so the thickness sector is not double-counted — together with the supplied
  background first jets `{∂W_bg, ∂μ_R,bg}` treated as symmetry-breaking spurion data … and construct
  **every** scalar bilinear in the DOF data allowed by the S11b symmetry group with those spurions carrying
  their transformation. … **Independence is judged as field bilinears with B1's constraint not applied;
  carry every independent invariant with its own free symbolic coefficient.**"
- **The S11b symmetry group (§1c ~109–122, and S11b B0):** in-plane **O(3)** — rotations **and parity** —
  with the wave/contrast reflection `w→−w` **unbroken**. The spurion is a vector (a background first jet) and
  carries that transformation; it breaks in-plane isotropy but not parity.
- **§3a (262–265):** "⛔ Do **not** obtain the variable-coefficient energy by the substitution `W₀→W_bg(y)`
  … that smuggles the answer and omits the newly admitted terms; ⛔ and do **not** forbid a symmetry-allowed
  new invariant. Whether any new invariant appears, and with what coefficient, is the computed result."
- **(deferred-scope reference) §3b (276–278):** "Retain the full spatial dependence of every background
  coefficient … and its first jet; **do not freeze a coefficient at its constant binding before
  differentiation.**" — the obligation the operator (§2C) violates; repaired in #89b, not here.

⚠ Design constant, not a target: each new-invariant coefficient carries exactly **one** background first jet
(`σ_W¹`); a background factor differentiated any number of times stays one factor, and only a **product** of
background factors (`σ_W²`+) is dropped by the retained-grade projection.

## 2 · The defects (SUPPLIED, code-verified)

**A · The new-invariant family is a hand-picked subset.** It is sourced from **four parallel, position-keyed,
length-8 hand-coded lists**: `newInvariantExpressions[spurion,u,theta,ew]` (`…mathematica_audit.wl:417–436`),
`newInvariantLabels` (L410–415), `newCoefficientSymbols` (L390–398), `newCoefficientNames` (L400–408). Both
`constructEnergyData` (L501–520) and `constructFullFieldBackgroundEnergy` (L570–576) consume them;
`basisRepresentativeIndices` (L622), from `independentRepresentativeIndices` (L599–618, incremental
`MatrixRank` L612 over `rankVariables` L590–597), drives the emitted count (L1244), the new-invariant set
(L1260), and the dimension/coefficient tables (L1750, L1831). ⇒ Completing the four lists to the full
symmetry-allowed family propagates through the existing rank selection; the count stays a computed rank.

**B · The exact thickness map is not imposed on the new invariants.** The uniform terms use the mapped
local field `localEw = (WZero/anchoredWidth) ewWave` (L453, used L463–467), but the new invariants are built
from the **raw** field: `newInvariantExpressions[spurion, uWave, thetaWave, ewWave]` (L505–506) and
`newInvariantExpressions[…, ewVariation]` (L572–574) pass `ewWave`/`ewVariation`, not `localEw`. This
violates §3a (246–247, "map imposed before any independence/rank test") and changes the retained `η¹σ_W¹`
grade of the e_W-bearing invariants.

**C · [KNOWN — DEFERRED to #89b, NOT fixed here] The slab operator freezes the background jet.** The slab
operator `evaluatedModel` (L898–911) applies `applyProfile` (L363; `profileDerivativeRules` `higherRules`
L304–307 zero every 2nd/3rd background jet and replace the first jet by a constant symbol) to the energy
records **before** the Euler–Lagrange variation `constrainedRows`→`variationalSource` (L790–795). So the EL
treats the spurion as a constant and drops the higher jet a divergence must generate — the opposite of §3b
(276–278). The live path exists (`applyBackgroundProfileWithGeneratedJets`, L360) but is used only for
admissibility (L1337–1341), not the operator. ⚠ This directive does **not** repair this; it requires the
diagnostic of §5.G so the `.out` documents the freeze. (⚠ Do not cite `constructFullFieldBackgroundEnergy`'s
`gradient[anchoredWidth]` at L568 as evidence the operator is unfrozen — L568/L576 are the admissibility
constructor, a different consumer.)

Ground truth: the committed `.out` emits `ENERGY_BASIS_COUNT … "COUNT_OPERAND" -> 26` (both anchorings),
`RANK_SELECTED_INDICES -> {1,…,26}` — the **incomplete-family** value and the sole public regression target
of §4.

## 3 · The object to realise — a systematic, complete O(3)-Kronecker enumeration

Replace the hand-picked family with a **systematic in-script enumerator** and derive the four lists from it.
Name the object precisely — do not paraphrase it into looser recipe words:

1. **Building blocks.** The DOF perturbation data `{u, ∇u (the full rank-2 shear ∂u_i/∂x_j, of which the
   trace is divU), θ, ∇θ, e_W, ∇e_W}` with the exact map `e_W,bg=(W_0/W_bg)e_W` imposed **before** the rank
   test (impose it on the new invariants exactly as the uniform terms already do — use the mapped local
   field and its mapped gradient, not the raw field); and the spurion `s` = the background first jet
   (`∂W_bg` / `∂μ_R,bg`), a vector.
2. **Enumerate EVERY distinct O(3) full contraction.** Form each candidate as one spurion factor × a
   quadratic in the DOF data, and contract **all** tensor index slots to a scalar **using the metric δ only
   (Kronecker contractions)**. A rank-2 factor such as the shear admits **several independent** full
   contractions against a spurion and a vector — enumerate **all** of them; the current family keeps only
   some. ⛔ **Forbid Levi-Civita (ε) contractions**: an ε-contraction is a parity-odd pseudoscalar, and
   S11b parity / `w→−w` are unbroken, so it is not a symmetry-allowed energy invariant. Symmetry-reduce and
   deduplicate to distinct scalars. The bare invariant need **not** be of energy dimension — the free
   coefficient supplies the remaining dimension.
3. **Let the existing rank selection compute the independent count.** Pass the full generated family through
   `independentRepresentativeIndices`; carry every independent representative with its own free symbolic
   coefficient (generate fresh coefficient symbols/labels/names for the new representatives, following the
   existing `WJET_…`/`MUJET_…` string convention — Wolfram symbol names carry no `_`, only the emitted string
   labels do). **Derive** each new invariant's dimension independently from its factor content and compare it
   to any stored value before forming the coefficient dimension (§5.E); regenerate every table keyed by the
   label list (dimensions L1750, coefficient-dimension map L1831, omissions L1275) from the completed labels.

⛔ Do not hard-code any count, any specific invariant form, or any operator term. If you type a specific
surviving-invariant count or paste a specific "missing" form, stop — that is the defect (rule 2 corollary 1:
a hand-typed CAS object is still hand-typed). The only by-hand symbol combination permitted is the
action/energy construction and the ansatz (already built); the family and its count are reached by the
enumeration + rank. The 10-term uniform family (`uniformLabels`, L382–386) carries no spurion and must keep
computing exactly what it computes now.

## 4 · The one acceptance you may check against — the incomplete-family public regression

Provide a **dedicated switch** that restricts the enumeration to **only the eight originally-coded forms per
source** (the committed hand list), leaving everything else as built. Under that switch, and only under it:

- `ENERGY_BASIS_COUNT`'s `COUNT_OPERAND` = the committed public **26** (both anchorings). Emit the
  restricted-switch count and its symbolic residual against `26` (print the residual; do not `Assert` it
  before emitting — §6).

⛔⛔ **26** is the **only** target in this directive. The completed-family count, its per-source count, and
the operator structure are **withheld on purpose** and are **not** acceptance criteria. ⛔ Do not tune the
enumeration toward any count — the O(3)-Kronecker family of §3 is fully determined by the physics, so emit it
as computed; the comparison against the reference and the sibling engine happens **outside this build**,
where a disagreement is a **finding**, not a build failure.

## 5 · Controls and independent routes (must remain able to FAIL)

- **A · Form ablation (decisive).** The §4 restrict-to-8 switch is a **form** ablation: it must **move**
  `COUNT_OPERAND` back to `26` and move the emitted new-invariant set back to the committed 8/source. A
  completed count byte-identical to the restricted count means the enumeration added no independent direction
  (rule 14). Emit both counts and their difference.
- **B · Completeness / uniqueness of the enumeration (the real protection).** The completed family must be a
  **superset** of the committed 8 and must be the **closed** O(3)-Kronecker family. Emit: (i) the incremental
  rank showing the committed 8 per source lie within the span of your generated family (add no rank on top);
  (ii) evidence that adding **one** ε (Levi-Civita) pseudoscalar candidate **would** raise the raw rank (so
  the parity exclusion is load-bearing, not vacuous) — emit that trial rank as an operand, then exclude ε
  from the family; (iii) that your generated independents are pairwise independent under `rankVariables` (no
  zero rows, no duplicate directions). ⛔ Do not describe the mere subspace check (i) alone as completeness.
- **C · Coefficient control (must NOT move the count).** Rescaling a generated candidate's coefficient by a
  nonzero constant is arithmetic; `COUNT_OPERAND` and the selected rank must be invariant. Emit enough to
  distinguish this from the form ablation.
- **D · Per-anchoring independent rank (fixes a tautology).** The engine computes `basisRepresentativeIndices`
  once from `LAB_HELD` (L621–623) and reuses it for both branches (L1208–1211), so the two anchoring counts
  agree by construction. Compute the independent rank and selected span **separately for each anchoring**,
  emit their span-equivalence residual, and reuse common indices only after that residual is emitted (do not
  gate on it — emit it).
- **E · Independent dimension derivation (fixes a self-validating check).** `newInvariantDimensions` is a
  position-keyed hand list (L1750–1764); coefficient dim = energy − invariant dim (L1766–1768) and
  homogeneity adds them back (L1828–1831), so a wrong invariant dimension passes automatically. Derive each
  new invariant's dimension from its generated factor metadata, emit the residual against any stored value,
  then form the coefficient dimension.
- **F · Thickness-map residual.** Emit the difference between the new-invariant family built with the mapped
  local field vs the raw field (a computed object showing the map is imposed and changes the retained grade),
  per source.
- **G · Operator-freeze DIAGNOSTIC (documents the deferred #89b defect — emit, do not repair, do not assert).**
  Emit the independent **rank of the completed operator's EL rows reduced two ways**: with the production
  frozen path (`applyProfile` before EL, as `evaluatedModel` does now) and with the live path (EL taken while
  `widthBase`/`modulusBase` are still functions, then reduced by `applyBackgroundProfileWithGeneratedJets`),
  and their difference — per anchoring. ⛔ Do not assert them equal; a nonzero difference is the emitted
  record that the operator still freezes the jet (to be repaired in #89b). This is the load-bearing rule-17
  object; the energy-monomial frozen-vs-live rank (below) is only a Hessian-in-energy guard.
- **H · Hessian-in-energy guard.** Emit the completed energy-monomial rank with vs without the background
  2nd-jet atoms substituted to zero. For a correct first-jet energy basis this difference is 0 (no Hessian
  belongs in the energy); it must **move** if a Hessian atom is wrongly admitted into the energy family.
- **Existing controls must still run and pass** (`REP_INVARIANCE`, `INDEPENDENCE`, `FORM`, `UNIFORM`,
  `HOMOGENEITY`, and every currently-emitted admissibility/operator tag). ⛔ Do not weaken any; keep both
  anchorings (`LAB_HELD`, `MATERIAL_ADVECTED`) and both density representatives throughout.

## 6 · Method and script obligations (non-negotiable)

1. The script **PRINTS computed objects; it does not state conclusions.** Every emit payload is a CAS object.
2. **PRINT the residual; do not `Assert` it** before emitting. Compute → emit → only then may a structural
   guard assert. The restrict-to-8 residual against `26` is emitted, then may be guarded.
3. Interpretation belongs to the step record, not the script.
4. No hard-coded results; no tautological residual (zero for any input tests nothing — corollary 3); no
   emission conditioned on a payload's value (corollary 4). Tag names name the object, never its value or the
   shape of the answer (corollary 2).
5. Match the engine's Wolfram idioms (strip `ConditionalExpression`, pole idioms `1/x==0`).
6. ⚠ Operational: this is a full engine regeneration (~97 MB `.out`). Do **not** wrap the full run in
   `timeout 600` — that is the review-ablation budget, and the full WL run is bounded by output-silence /
   memory, not elapsed wall-clock (a 600 s kill of a live regeneration wastes a round). Wrap only your own
   short ablation/diagnostic kernels in `timeout 600`, one kernel at a time (the licence has two seats, and
   other work may hold one). Launch with `--sandbox danger-full-access` and `-c model_reasoning_effort=xhigh`.
   Run the full engine to completion and verify the `.out` is regenerated and non-empty (a real run is large).

## 7 · Builder report (return these)

- The unified diff of every definition changed, with L-numbers, and the enumeration procedure you built.
- The completed `ENERGY_BASIS_COUNT` `COUNT_OPERAND` (both anchorings) and per-source new-invariant count —
  whatever they are — plus the emitted new-invariant label set.
- The restrict-to-8-switch `COUNT_OPERAND` and its literal residual against `26`.
- Controls A–H literal outputs: form ablation Δ; completeness (8-in-span, the ε-trial rank, pairwise
  independence); coefficient control; per-anchoring ranks + span residual; dimension residual; thickness-map
  residual; the operator-freeze frozen-EL vs live-EL rank + difference; the Hessian-in-energy guard Δ.
- Confirmation the existing controls still emit and pass, and that the full engine ran to completion with the
  `.out` regenerated (byte size). ⛔ `.out` is git-annex — do **not** `git add -f` it; leave saving to the
  orchestrator.

## 8 · Supplied vs computed

SUPPLIED and unfalsifiable within this build: the spec obligation (§1), the four hand-coded lists / thickness
asymmetry / operator-freeze code paths (§2), the incomplete-family public value `26` (§4). COMPUTED and
subject to external review: the completed basis count, the completed new-invariant set, every control
residual, and the operator-freeze diagnostic. ⛔ Do not present a passing restrict-to-8 regression as
confirmation of the completed physics — it confirms only that the incomplete-family limit is unchanged. ⛔ Do
not present this build as repairing the operator freeze (§2C) — that is #89b.
