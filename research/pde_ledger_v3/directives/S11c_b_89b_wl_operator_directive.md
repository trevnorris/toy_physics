# S11c-b #89b (WL) — un-freeze the §3b operator: take the EL variation with the background jet tower LIVE

## 0 · Role and single deliverable

Modify **in place** the Wolfram engine
`research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl`
so that the emitted variable-coefficient §3b slab operator — and every object derived from it (the §3c
coupling kernel, the μ_θ operator, the divergence-form source, and their term-origins) — is constructed with
the background coefficients kept as **live functions of the spatial coordinates through every differentiation
that can raise a background-jet order**, and the jet-retaining reduction applied only as the **last** step.
The operator must carry the full first-amplitude-order jet tower of every background coefficient, including
the retained-grade second (and, in the kernel, third) spatial jets — the background curvature that the current
construction zeroes. Re-run the engine to completion. The run **regenerates one tracked file**: its output
`research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out`. That regeneration is
authorized and expected; change nothing else in the tree.

Also fold the one control-quality fix in §6. This is a **spec-compliance repair of an existing implementation
defect** (a freeze-before-differentiate on the operator), not new physics and not a spec change. The defect,
the spec obligation it violates, and the correct pattern already present in this file are supplied below as
verified inputs — you are not asked to re-discover them, only to make the emitted operator satisfy the spec.

⛔ You build your **own** production integration of the live path (§3). ⛔ Do **not** import or read the SymPy
engine — the two engines re-derive independently; their disagreement is the measurement (rule 1).

## 1 · The spec obligation (SUPPLIED, verbatim — the object to satisfy)

From `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`:

- **§3b (274–282):** "compute the first-order equations of motion … in **divergence form** — the variable
  coefficients sit inside the in-plane divergences, so their first jets appear explicitly. **Retain the full
  spatial dependence of every background coefficient (`μ_R,bg`, `W_bg`, `ρ_4D,bg⁰`, `ρ_br,bg⁰`, and the
  `Σ_E⁰` map) and its first jet; do not freeze a coefficient at its constant binding before
  differentiation.**"
- **§3a (248–254):** "**first-background-jet order** bounds the number of independent background-**amplitude**
  factors … **not the spatial-derivative order**. … A divergence or variation may generate a higher spatial
  derivative of a background profile, and that term is retained at the background-bookkeeper order of its
  originating factor: a second spatial derivative of `W_bg` is still first order in background amplitude,
  `O(η)`/`O(σ_W)`, and is **not dropped**."
- **§1d (163–168):** the variable-coefficient representative difference (`c∇·F ≡ −(∇c)·F`) is "physics in the
  operator/kernel, not a representational identity."
- **§3c (290–319):** the coupling kernel is a **weak variational restriction of the §3b operator itself**
  (compact-support IBP), ⛔ "not … a parallel direct-variation route."

⚠ Design constant, not a target: a background factor differentiated any number of times in space stays **one**
background-amplitude factor; only a **product** of two background factors would raise the amplitude grade and
be truncated. **This is parity with the PY side: PY #89 already un-froze both its §3a-quotient and its §3b
operator (kept the jet tower live). No spec change; a different mechanism, the same operator object.**

## 2 · The defect (SUPPLIED, code-verified — the current construction freezes before it differentiates)

The emitted production operator is built by `evaluatedModel` (located by name; post-#89a ~L1130–1195), whose
outputs form `mainModels` (~L1798, `evaluatedModel["EULERIAN", …]`) and feed `SLAB_OPERATOR` (~L1806),
`MU_THETA_OPERATOR`, `COUPLING_KERNEL` (via `extractCoupling`, ~L1376/1829), the divergence-form source, and
their term-origins. Inside `evaluatedModel`:

- the energy is built on the LIVE functions `widthBase`/`modulusBase` (`constructEnergyData`, ~L1136), but
- **`applyProfile` is applied to the energy records (~L1137–1140), the constraint (~L1141), `width`/`rhoBrValue`
  (~L1144/1146), the faces (~L1148), and the evolution terms (~L1152) BEFORE** `constrainedRows` (~L1143) and
  `variationalSource` (~L1181–1184) take the EL variation.
- `applyProfile` (~L363) reduces with `profileRules` → `profileDerivativeRules`, whose `higherRules`
  (~L304–307) map **every** jet with `2 ≤ i+j+k ≤ 3` to **0**.

⭐ **The measured mechanism (from the two decision legs — this is WHY a naive fix fails).** The retained-grade
Hessian curvature is the **product-rule content of a live divergence** (`∂(c·F) = ∂c·F + c·∂F`). Once the
profile rule replaces `W^{(n)}(x)` by a **constant** jet-symbol, that product rule is gone: a live-then-reduce
that does not differentiate *through the coefficient as a function* reproduces the **frozen** operator
byte-for-byte, and a reduction applied *before* a later differentiation (e.g. the kernel's `curl`/`div`)
drops the jet that differentiation would have raised. So the freeze is not one `applyProfile` call to move —
it is the **ordering of reduction relative to every order-raising differentiation.**

The engine already emits an honest diagnostic — `operatorFreezeRankDiagnostic` (~L1219–1258) /
`WL_S11CB_OPERATOR_BACKGROUND_FREEZE_DIAGNOSTIC`. ⭐ **Keep it unchanged.** It is an abstract per-record EL
rank (a different object from the production operator); ⛔ do **not** add a production-rank residual against
it — that would turn a withheld magnitude into a tuning target.

## 3 · The object to realise — the construction INVARIANT (name the object, not a recipe)

Realise the variable-coefficient operator whose live product-rule / retained-grade curvature is **explicit
before any reduction to jet-symbols**. The invariant that makes it spec-compliant:

> **`widthBase`/`modulusBase` and the density constraint stay LIVE FUNCTIONS through EVERY differentiation
> that can raise a background-jet order** — the bulk EL variation (`constrainedRows`/`variationalSource`), the
> in-plane divergence handling, the §3c weak-block extraction (`extractCoupling`'s IBP split and its Helmholtz
> `curl`/`div` on the remainder), and the face EL — **and `applyBackgroundProfileWithGeneratedJets` (~L360,
> backed by `profileRulesRetainingGeneratedJets` ~L339–357) is the ABSOLUTE LAST step**, once no further
> order-raising differentiation remains on that object.

Two consequences the legs proved you must honour, ⛔ or the fix is partial and still passes a weak control:

- **The §3c kernel split must be taken on the un-reduced operator with its `Inactive[Div]` structure intact.**
  `extractCoupling` → `splitDivergenceRows` (~L1340) separates `Inactive[Div]` terms (IBP via `weakPairTerm`)
  from the direct remainder (Helmholtz `curl`/`div`, ~L1356/1368). If you fully activate/expand the
  divergences first, `splitDivergenceRows` mis-classifies **everything** as remainder and the Helmholtz step
  hits the Hessian terms — a **different weak block** than §3c defines. Keep the coefficients live but the
  `Inactive[Div]` split intact into `extractCoupling`; reduce the kernel **after** the pairing.
- **⛔ `rawModel` (~L1049) + `processModel` (~L1108) + reduce is NOT a sufficient realisation.** It never
  activates the divergences (so it reproduces the frozen operator), and `processModel` (~L1120–1121) leaves
  `DIVERGENCE_FORM_SOURCE` unreduced. You may reuse or refactor any function, but the emitted objects must
  satisfy the invariant above — verified by the §5 controls, not by relocating one `applyProfile` call.

## 3bis · Tractability — the operator regen MUST complete in bounded time (SUPPLIED diagnosis + approach)

⚠ The naive live-operator construction **does not finish in bounded time** — measured 2h+ at 100% CPU with no
output, before any `Series` even runs. Two distinct blow-ups were localized by an independent timing
diagnosis (evidence `~/.s11_build/S11c_b_89b_perf_diag/`, scripts + literal stdout). Both are ordinary
Mathematica evaluation traps, **not** physics — they are SUPPLIED so you do not re-discover them:

1. **`activateSpatialDivergences` (the top-down `FixedPoint[ReplaceAll[…]]`, ~L1197) does not terminate on the
   multi-component operator.** `variationalSource`/`linearVirtualVariation` emit **nested** `Inactive[Div]`
   (an outer Div wrapping coefficients that still hold inner Divs). Top-down `ReplaceAll` activates the OUTER
   Div first, so `D[…]` differentiates a vector still containing an inner `Inactive[Div]`, and `D` of an
   inactive `Div` spawns fresh `Inactive[Div]` terms faster than the rule removes them. Measured per pass:
   `#Div 16→19→22→25→28→31→34…`, `LeafCount 32k→57k→86k→121k→…` — monotonic, unbounded. Individual rows
   converge in ~0.05 s; the **bundle** hangs. This is the 2h wall.
2. **The grade truncation's global nested `Series` (`backgroundSeries`/`truncateScalar`, ~L99–125) is
   quadratically wasteful.** A single global `Series[·,{etaBg,0,1}]` then `Series[·,{sigmaW,0,1}]` on the whole
   operator fully distributes the grade-≥2 expansion just to discard it. Measured: one 146-term row's
   `{etaBg,0,1}` Series ≈ 49 s; across all rows × 4 cases that alone is hours.

**REQUIRED PROPERTY:** the operator regen must **complete** (terminate — the failure was non-termination, not
slowness) AND be **algebraically equivalent to the naive computation on the retained-grade content**.
⛔⛔ It MUST NOT drop a jet term, truncate a retained-grade jet, or reduce a coefficient before a later
differentiation — any of those silently RE-FREEZES the operator and defeats the entire build. Tractable means
*cheaper*, never *smaller*. ⛔ **Wall-clock is NOT a pass condition** — the acceptance is the equality residual
below, and the run is silence/RSS-bounded (§7), not time-bounded. Dropping retained-grade content to make the
run faster is a **failed build even if every §5-control-8 residual is 0**; the partial-drop teeth are §5.1
(Hessian-zero) and §5.5 (per-surface atom set), which must still bite.

**APPROACH (diagnosis-verified — realise it; ⛔ you write the production code, this is the mechanism not a
paste target):**
- **Make the divergence activation terminate AND evaluate.** ⚠ The precise trap (decision-leg-verified): the
  hang is `ReplaceAll` of the rule's `If[SameQ[c, spatialCoordinates], Sum[D[…]], leftover Div]` RHS **inside an
  `Association`** — the `If` condition does not evaluate in the Association context, so the `Inactive[Div]` is
  never removed and `ReplaceAll` re-fires unboundedly (`nDiv 2→3→4→5→6→7…`). The *same* expression as a `List`
  terminates (`2→0` in two passes). ⛔⛔ Adding a Div-free `FreeQ` guard but **keeping the `If`** STILL hangs on
  an `Association` — so "innermost-first" alone is not the fix. The fix must do **one** of:
  (i) drop the `If` — use a `/; SameQ[c, spatialCoordinates]` **pattern constraint** with a fully-**evaluating**
  bare `Sum` RHS, so each generated `D` collapses to an explicit `Derivative[…]`; or
  (ii) never run the rewrite on an `Association` — `Map` the activation over the list-valued rows, or `Join`
  the rows to a list first (`elRowVector` ~L1205 already does this; there top-down terminates AND equals
  innermost, verified). Either way the divergences activate innermost-first (a Div is resolved only once its
  vector is Div-free) so `D` never differentiates an unresolved inner Div.
  **Name the object:** the fully-activated divergence `Σ_i ∂v_i/∂x_i` in **derivative-normal form** — ⛔ no
  leftover `Inactive[Div]`, ⛔ no held/unevaluated `D`/`Sum`/`If`. ⚠⚠ "no `Inactive[Div]` left" is NOT enough:
  a container can hold an unevaluated `Sum[D[…]]` (Div-count 0) whose derivative never fired, and then profiling
  substitutes the background before it fires and **silently drops the retained jets** — the exact re-freeze this
  build exists to prevent. Verified equal to top-down wherever top-down terminates; apply to every activation
  consumer (`elRowVector` etc.), reproducing the top-down result, not changing existing behaviour.
- **Exploit `Series` linearity.** `Series` is linear, so apply the grade truncation per top-level `Plus`
  summand instead of one global call; each summand's expansion is small and grade-≥2 cross-terms across
  summands never form. Verified: one row 59 s→1.3 s, a heavier row (was aborting >120 s)→14 s,
  `PossibleZeroQ[Expand[perSummand − global]] = True`. A stronger grade-truncated-product variant is
  allowed **only** if it too emits a zero equivalence residual (below).

## 4 · Blast radius (every surface the fix must reach and keep consistent)

- **Emit sites**, each for **both** anchoring branches and **both** density representatives (~L1798–1838,
  ~L1951–1987): `SLAB_OPERATOR`, `SLAB_OPERATOR_TERM_ORIGINS`, `MU_THETA_OPERATOR`, `COUPLING_KERNEL`,
  `COUPLING_KERNEL_TERM_ORIGINS`, and — nested on the `SLAB_OPERATOR` payload via `modelRecord` (~L1290–1292)
  — **`DIVERGENCE_FORM_SOURCE`** and **`FACE_SHAPE_SUBSTRATE`**. All of these must be reduced consistently with
  the operator (live-then-reduce), or be *explicitly* pre-reduction geometry stated as such. ⛔ Do not leave
  `DIVERGENCE_FORM_SOURCE` raw or frozen beside a corrected operator.
- **The density constraint / `THETA_VARIATION_RULE`** (built in `constrainedRows` from the constraint, ~L1028;
  currently `applyProfile`'d at ~L1141) is the ρ_4D,bg⁰/ρ_br,bg⁰ relation §3b requires live — put it in the
  same live-through-EL set as the energy, ⛔ do not `applyProfile` it before `constrainedRows`.
- **`extractCoupling`** must consume the corrected operator per §3 (Inactive-Div split preserved).
- **Term-origin provenance** (bulk-energy vs face/flux vs advective) must be re-derived from the corrected
  rows, not carried frozen; emit an origins-sum-minus-production residual for both slab and kernel.
- ⚠ **§3d ADMISSIBILITY** (`backgroundBalanceFromModel`, ~L1845–1859) already does expanded EL on the
  full-field energy then `applyBackgroundProfileWithGeneratedJets` — it is a **separate, already-live** route.
  ⛔ Do **not** fold it into the wave operator; **verify** its emitted operand is unchanged by this fix and
  report that.
- **Faces are not a second Hessian factory** (the Eulerian LAB_HELD virtual work has no `W` jets; the virtual
  displacement is undifferentiated) — but the face rows inherit `muTheta`, so they are corrected transitively
  by the bulk-EL + constraint fix. Keep `formOnly`/corrupted/parametric branches (`formOnly` early-return
  ~L1165; parametric `formOnly->True` ~L2064/~L2203; the `corrupted` 6th-arg call ~L1979) and the §5c
  uniform-limit (~L2203–2215) working.

## 5 · Controls and independent routes (must remain able to FAIL)

⛔ Emit operands and residuals; never assert a residual (§7). Every control re-enters the chain at the
energy/action, ⛔ never at an already-computed result (§7).

1. **FORM ablation — MANDATORY, the load-bearing one (rule 14): Hessian-zero of the corrected operator.** On a
   copy of the **corrected** operator rows, set the order-2 (and order-3, where the kernel generates it)
   background-jet atoms → 0, and emit the ablated object. It must **move back toward the frozen construction**
   (the curvature disappears). ⭐ This — not replaying the old path — is the form ablation that proves the
   un-freeze carries physics. Take the ablated atoms **from the emitted EL**, ⛔ not from a hand-inserted term.
2. **Regression witness (keep, but label as regression, not the form control).** Re-freeze
   (`applyProfile`-before-EL) on a copy and emit that operator; it must reproduce the **current committed
   frozen operator**. This is the public regression anchor.
3. **Tower-depth control.** Emit the operator reduced with the retaining depth **truncated one order below**
   the generated depth (an emitted object must MOVE) and **extended one order above** it (must NOT move). ⛔ Do
   not state a numeric depth as a pass condition (rule 5) — emit live / depth-truncated / depth-extended and
   let each residual move or not.
4. **Grade-support control (global, across `Inactive[Div]`).** The current `truncateScalar` (~L99–117) protects
   inactive operators as tokens and truncates their exterior and interior **separately**, so a σ_W outside ×
   σ_W inside a `Inactive[Div]` (reachable via `muTheta·thetaRule` ~L1031 and the face path) is **not**
   truncated though plain σ_W²·f is. Require a **global multigrade projection** that accounts σ_W/η factors
   **jointly across** `Inactive[Div]` boundaries, and emit a grade-support control proving σ_W² is absent
   (including across a Div) while the allowed ησ_W survives.
5. **Per-surface Hessian witness.** For **each** emitted surface (`SLAB_OPERATOR` rows, `MU_THETA_OPERATOR`,
   both kernel blocks, both origin payloads, `DIVERGENCE_FORM_SOURCE`) emit the live-minus-frozen background-
   jet **atom set**, computed from that surface's own EL. A corrected slab row must not coexist with a frozen
   or raw source/kernel — each surface's difference set must show the retained-grade atoms it should carry.
6. **One-sided corruption (independence of the two anchorings).** Inject a defect into **one anchoring
   branch's energy/background input** and **rerun the complete constructor**; emit both branches' operators —
   only the corrupted branch must move. ("Routes" here means the two anchoring branches.) ⛔ Do not corrupt a
   post-reduction result.
7. **Uniform-limit regression (§5c).** The corrected operator's uniform limit must still reproduce S11b's
   decoupled zero for the off-diagonal kernel; emit the limit and its residual against the supplied `N1a`
   uniform value (§1d), do not assert it.
8. **Tractability equivalence controls (§3bis — MANDATORY, and these were re-hardened by a decision leg that
   broke the weak version).** The tractable computation must be proven **equal to the naive one on the full
   operator**, ⛔ never assumed, and ⛔ never inferred from a leftover-`Div` count. Emit:
   (a) **Full-operator equality vs a per-row top-down reference.** The naive top-down activation TERMINATES
   **per row** (only the whole-`Association` bundle hangs) — so build the reference by top-down-activating each
   row separately (as a `List`/`Join`, not the `Association`) and reassembling the container, then emit
   `tractable_full_operator − reassembled_per_row_top_down` residual (must be 0). ⛔ This must cover **every**
   production row, including at least one **nested-`Div`** row (e.g. `U_INTERNAL`, whose top-down terminates in
   ~0.05 s) — a nest-free row (e.g. `EW_INTERNAL`) does NOT exercise nested activation and is not sufficient.
   (b) **Derivative-normal postcondition.** Emit the count of leftover `Inactive[Div]` **and** the count of
   held/unevaluated `D`/`Sum`/`If` on background functions after activation (both must be 0). ⛔ Div-count 0 with
   a held `Sum[D[…]]` is a re-freeze that profiling will drop (measured counterexample: controls green while the
   full residual carried `−σ_W·wJetTwo·u − 2σ_W·wJetOne·∂u`).
   (c) **Series-linearity residual.** `per_summand_series − global_series` on a bounded case with a
   `thetaRule`-like **rational** (`n/d`) coefficient present (must be 0). ⚠ If the stronger grade-truncated-
   **product** variant of the reduction is used, (c) must exercise a **product + re-truncation**
   (`Series[f]·Series[g] ≠ Series[f·g]` until re-truncated) and show its residual 0 too.
   ⛔ Emit each residual (rule 2); a nonzero residual is a finding that the speed-up dropped or altered a term.
   Apply (a)/(b) to every activation consumer (`elRowVector` included).
9. Keep every existing control that still applies (`REP_INVARIANCE`, `INDEPENDENCE`, `FORM`, `UNIFORM`,
   `HOMOGENEITY`, `operatorFreezeRankDiagnostic`) emitting.

## 6 · The folded control-quality fix — §5.E dimension residual

`…audit.wl:2269–2290`: `RESIDUAL_DERIVED_MINUS_STORED` differences two copies of the same
`dimensionGradient + Σ(factor dims)` sum (`deriveInvariantDimensionFromFactorContent`, ~L2269, vs
`blueprint["STORED_INVARIANT_DIMENSION"]`, ~L482), both drawn from the **same** `dofFactorSpecifications`
table — structurally 0, cannot catch a wrong invariant dimension (a check auditing its own input; verified:
mutating the shared table leaves the residual 0). Replace it with a **genuinely independent** second operand:
the dimension obtained by assigning **primitive atom** dimensions to the emitted invariant **expression**
(`uOne`, `D[uOne, xOne]`, `Derivative[1,0,0][widthBase]/WZero`, θ, e_W, …) and reducing the actual term. ⛔ The
second operand must **not** be computed by walking `FACTOR_NAMES` / `dofFactorSpecifications` (that reproduces
the tautology). If no independent route is achievable, **delete** the residual and emit only the derived
dimension with a one-line note that there is no independent second route here. ⛔ Do not keep a structural zero
dressed as a check.

## 7 · Method and script obligations (non-negotiable)

> **1. The script may PRINT computed objects. It may NOT state conclusions.** Every emit payload is a CAS
> object — an expression, a rank, a symbolic residual — never prose describing a result.
> **2. PRINT the residual; do NOT assert it.** Compute → emit → *then* (optionally) guard. An `assert
> residual == 0` before the emit turns an informative value into a crash and hides it from review.
> **3. Interpretation belongs to the step record**, not the script.

⭐ **The ONLY place the physical symbols may be combined by hand is in constructing the energy/action and the
ANSATZ. Every other expression — the operator, the kernel, every residual, every control's ablated object —
must be REACHED BY COMPUTATION from the energy. Every control re-enters the chain at the energy, never at a
result.**
⛔ A tag **name** names the object, never its value, rank, sign, or the shape of the answer.
⛔ No tautological residual: before emitting a difference, confirm the two operands came from **independent
routes**; where none exists, emit the object and say so.
⛔ Emission of a tag must never be conditional on its payload's value (only on which object/anchoring it is).
⛔ Do not hard-code any operator rank, termcount, or specific operator term. If you type a specific number as
a target or an expected operator form, you have inverted the method.

Mathematica hygiene: this is one full-engine run (`--sandbox danger-full-access`). ⛔ Do **not** wrap the full
regeneration in a shell `timeout` — it is silence/RSS-bounded, not budget-bounded, and the `.out` is ~156 MB.
`ClearSystemCache[]` / `Clear[…]` between heavy cases as the file already does. One kernel.

## 8 · The one public regression you may check against; everything else is WITHHELD

- **PUBLIC** (given): the spec obligation (§1); the freeze diagnosis and mechanism (§2); the construction
  invariant (§3); the energy-basis count **40** (emitted and cleared in #89a); the current **frozen** operator
  (in the committed `.out`) — your regression control (§5.2) must reproduce it.
- ⛔ **WITHHELD** (rule 5): the corrected operator's rank, its per-row term counts, and any rank-gain figure.
  ⛔ Do not tune any object to a number. Your job **ends at compute-and-print**; the corrected operator's
  content is diffed against the sibling engine and the prior-step measurements on the review side, where a
  disagreement is a **finding**, not a build failure. ⛔⛔ The acceptance is the spec-compliant **construction**
  (the §3 invariant, the §5 ablations biting), verified by cross-engine **content/spans** — never a matching
  rank.

## 9 · Builder report (return these)

1. The diff summary: which functions changed, where the coefficients are kept live, and where the single
   final reduction now enters (bulk rows, μ_θ, faces, evolution, constraint, kernel).
2. How you preserved the `Inactive[Div]` split into `extractCoupling` (§3), and how `DIVERGENCE_FORM_SOURCE` /
   `FACE_SHAPE_SUBSTRATE` are made consistent.
3. The absolute path of the regenerated `.out` and its size; confirm the `datalad`-tracked file regenerated.
4. For each control in §5: the emit tag names and a one-line statement of what its literal payload showed
   (moved / didn't move / residual value) — ⛔ not an interpretation, the literal object.
5. Confirmation that the §3d admissibility operand is unchanged, and the §5.E residual is replaced or deleted.
6. Anything you could not make compute, and why — do not paper over it.
