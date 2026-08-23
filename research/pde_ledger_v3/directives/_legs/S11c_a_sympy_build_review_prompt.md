# Review — the S11c-a SymPy engine BUILD DIRECTIVE (wiring, not physics)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_a_sympy_build_directive.md`

This is an **orchestrator-written build directive** for the S11c-a SymPy engine. ⚠ It is **not** the physics
spec — the physics lives in the **closed** `directives/S11c_a_SHARED_PHYSICS.md` (HEAD `2926c71c`), already
reviewed and locked. This directive adds only the **chain wiring** (imports, exports, F1–F9, digests, publish
gate) plus one naming pin. Your job is to find a **wiring** defect that would corrupt the engine or its
export chain — ⛔ not to re-review the spec's physics. It gets two legs before the builder launches (rule 7).

## Read first
- `directives/S11c_a_SHARED_PHYSICS.md` (the spec it wraps — esp. §4 tasks T-0…T-i, §5 controls, §7
  names/F9/tags, §8 supplied-vs-computed, §3a's `μ_s^α`/§4 T-i's `μ_s` glyph).
- `directives/S11b_sympy_build_directive.md` (the S11b precedent this adapts).
- `scripts/S11b_exports.py` (the LEDGER this engine imports; check the bind targets exist:
  `c_s0`, `mu_R`, `rho_br`, `W_0`, `e_W`, `v_bulk_normal_0`).
- The chain rules it points at: `directives/S11_export_chain_decisions_v2.md` (F1, **F9** :171-186 totality,
  :203 prefix generalization), `directives/S11_C17_C18_spec_repair_decisions_v2.md:53-60` (T7),
  `S9_REWRITE_PLAN.md` (D1, D3 :217, class vocab :41).
- The rule-2 twin `directives/_measurements/S11c_a_sympy_build_directive.md` (spot-check the anchors).

## What to check
1. **Faithful pointing, ⛔ not restatement/paraphrase.** The directive must POINT at the spec's §§ and the
   inherited chain rules, not restate them. ⚠ The S11b analog's measured defect was a **paraphrase of F9**
   that dropped F9's "TOTAL over the imported row shapes" clause. Confirm F9 is applied **whole** and not
   re-rendered; confirm no §-restatement silently narrows the spec.
2. **Bind targets exist and are not re-declared.** `c_s0→LEDGER['c_s0']`, `mu_R→LEDGER['mu_R']`,
   `rho_br→LEDGER['rho_br']`, `W_0→LEDGER['W_0']`, `e_W→LEDGER['e_W']`, `v_dr→LEDGER['v_bulk_normal_0']` —
   do all six keys exist in `S11b_exports.py`? Are the varying profiles (`W_bg`, `mu_R_bg`, densities, …)
   correctly originating **here** with fresh names (§7), never reusing a constant key?
3. **Primary vs control split.** The directive exports **§4 T-0…T-i** and treats **§5 controls (5a/5b/5c)
   and §6 homogeneity controls as ablation evidence, NOT exports** (the S11b-`B8` analog). Is that split
   correct against the spec — is any primary object wrongly excluded, or any control wrongly exported?
4. **F9 / prefix / digests / publish.** Is the `F9c` prefix `s11c_a_` correct (are `s11_`/`s11b_` really
   taken in the imported LEDGER)? Is the "no locus protocol ⇒ inherited T5 default" claim right for S11c-a?
   Does `BUILD_INPUT_DIGESTS` pin exactly the three inputs `{this engine source, S11b_exports.py,
   S11c_a_SHARED_PHYSICS.md}`? Does F6 publish only if T-0…T-i completed?
5. **The three G8 deviations.** Comparator separate (net-new T7, not built here); D3 round-trip restored;
   `_RELATIONALS` reviver included (present even if S11c-a emits no relational — ⛔ but the directive must not
   manufacture one).
6. **The μ_θ pin.** The directive pins §3a's branched `μ_s^α ≡ μ_θ^α/ρ_br,bg^{0,α}` over T-i's
   un-superscripted glyph, treating `μ_θ` as the reserved S11c-b operand branch-anchored at `χ`. Is this a
   pure naming/anchoring convention (⛔ not new physics, ⛔ not constructing `μ_θ`)? Does it correctly foreclose
   the two engines anchoring `μ_θ` differently on MATERIAL_ADVECTED?
7. **No leaked value / no physics / no self-review machinery (rule 5, rule 2, rule 12).** Does the wiring add
   any expected value, acceptance criterion, or physics the spec withholds? Any check that would audit its
   own input, or a completeness registry the spec's task structure doesn't already carry?

## Physics filter
Report a finding only if it catches a way the **engine or its export chain** would be wrong (a bad bind, a
dropped/duplicated export, an F9 paraphrase, a wrong prefix/digest, a leaked value, added physics). ⛔ Not
spec physics (closed), not style. "Nothing survives" is weak evidence — say what you checked.

## Output
Findings most-serious first (source quote + directive line + the concrete wiring failure); then: is the
primary/control export split correct? is F9 applied whole? is the directive safe to hand to the builder?
