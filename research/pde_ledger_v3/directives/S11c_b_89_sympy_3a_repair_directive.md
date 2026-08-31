# S11c-b #89 (PY) — repair the §3a variable-coefficient Euler–Lagrange: keep the background jet tower LIVE

## 0 · Role and single deliverable

Modify **in place** the SymPy engine
`research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py`
so that every total in-plane divergence in its §3a energy-basis quotient and §3b operator construction is the
**variable-coefficient** one — the background profile and its spatial-jet tower `{W_bg, ∂W_bg, ∂²W_bg, …}`,
`{μ_R,bg, ∂μ_R,bg, ∂²μ_R,bg, …}` are differentiated, not held constant — to the jet depth each consumer
needs. Re-run the engine to completion. The run **regenerates two tracked files**: its output
`research/pde_ledger_v3/scripts/out/S11c_b_brane_operator_sympy_audit.out` **and** the generated ledger
`research/pde_ledger_v3/scripts/S11c_b_exports.py` (written at L2930). Both regenerations are authorized and
expected; change nothing else in the tree.

This is a **spec-compliance repair of an existing implementation defect**, not new physics and not a spec
change. The defect, the spec obligation it violates, and the (partial) correct pattern already in this file
are supplied below as verified inputs. You are not asked to re-discover them; you are asked to make every
frozen total-divergence compute the variable-coefficient object the spec already names, and to emit the
results.

## 1 · The spec obligation (SUPPLIED, verbatim — the object to satisfy)

From `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md`:

- **§3a (251–254):** "A divergence or variation may generate a higher spatial derivative of a background
  profile, and that term is retained at the background-bookkeeper order of its originating factor: a second
  spatial derivative of `W_bg` is still first order in background amplitude, `O(η)` / `O(σ_W)`, and is **not
  dropped**." — and (248–250): "**first-background-jet order** bounds the number of independent
  background-amplitude factors — powers of `η`/`σ_W` — **not the spatial-derivative order**."
- **§3a (255–256):** "Independence is judged as field bilinears with B1's constraint not applied; carry every
  independent invariant with its own free symbolic coefficient."
- **§1d (163–168):** "the uniform quotient does **not** lift trivially to variable coefficients … `c∇·F ≡
  −(∇c)·F` … representatives that were equivalent uniformly differ by a first-jet invariant that is **physics
  in the operator/kernel**, not a representational identity."
- **§3b (276–278):** "Retain the full spatial dependence of every background coefficient … and its first
  jet; **do not freeze a coefficient at its constant binding before differentiation**."

⚠ Consequence made explicit (this is the design constant, not a target): a background factor differentiated
any number of times in space stays **one** background-amplitude factor (`σ_W¹`). So every spatial jet of a
profile — 2nd, 3rd, … — is retained-grade and must survive; only a **product** of background factors
(`σ_W²` or higher) is dropped by the retained-grade projection. The construction rule §3a (244–265) is
otherwise binding and unchanged; ⛔ do not obtain the energy by `W₀→W_bg(y)`, ⛔ do not forbid a
symmetry-allowed invariant.

## 2 · The defect — the constant-coefficient divergence is frozen PERVASIVELY (SUPPLIED, code-verified)

The engine's total-divergence machinery treats a background spatial jet as a constant, so the **next** jet a
divergence must generate is silently zeroed. This occurs in **four** consumers that feed the emitted §3a
basis and §3b operator, not two. All four are the same defect and take the same correction.

- **A · §3a basis quotient.** `basis_euler_signatures` (L936) builds its divergence map (L940–945) only from
  the `fields` it is passed — `basis_fields` (L1025–1028), DOF-only; the spurion `bg` (L1014) is not among
  them, so `basis_dx` (L947–956) returns `0` on a spurion factor. `quotient_independent_indices` (L972) then
  selects independence modulo this frozen divergence, merging candidates that are independent once `∂(bg)`
  is retained. ⚠ Note the outer signature loop (L961–967) also iterates `fields` to take `δ/δ(field)`;
  the spurion must **not** enter that loop (§3 below).
- **B · §3b LAB_HELD strong rows.** `operator_from_density` (L1459) forms the U/θ/e_W rows with the **global**
  `dx` (L1473/1482/1489). `dx` (L616) reads `DERIVATIVE_MAP`, whose background entries (L612–613) are only
  `W_bg→grad_W[i]`, `mu_R_bg→grad_mu[i]`; there is no `grad_W[i]→(2nd jet)` entry, so `dx(grad_W[i])=0`.
  The θ operand `mu_theta_amplitude` (L1480–1483) is part of this. (This is the #88 blast radius.)
- **C · §3b MATERIAL_ADVECTED pullback.** `material_pullback` (L1348) builds `e_material = e_W +
  u·grad_W/W_bg` (L1357) and `theta_material` (L1355) and differentiates them with the **global** `dx`
  (L1363–1364); the pulled density then reaches `operator_from_density` (L1535). Because `e_material`
  carries `grad_W`, the frozen `dx` drops its `∂²W_bg` term before the MATERIAL slab-operator and the
  coupling-kernel MATERIAL representative are formed. Verified: the frozen-minus-live `e_material` derivative
  is `−σ_W·u·w1_profile_d1d·/(L_W·W_bg)` (nonzero).
- **D · the coupling differentiation cascade.** The coupling kernel differentiates the strong rows once more
  through `operator_divergence`/`operator_curl` (L1876/1881), built on `operator_dx` (L1850). Once rows A/B/C
  carry a 2nd jet, this outer differentiation must generate the **3rd** jet — still `σ_W¹`, retained — but
  `operator_dx`'s local map (L1864–1869) stops at the 2nd jet. Verified: `operator_dx` on a Hessian row keeps
  `σ_W·u_d1·w1_profile_d1d1/L_W` but **drops** `σ_W·u·w1_profile_d1d1d1/L_W²`.

**The partial correct pattern.** `operator_dx` (L1850) / `background_dx` (L2122) already clone
`DERIVATIVE_MAP` and add the **2nd-jet** entries `grad_W[j]→σ_W·w1_profile_dij/L_W`,
`grad_mu[j]→μ_R·σ_W·m1_profile_dij/(W0·L_W)`. This is the right shape but only one level deep, and it is not
used by A/B/C. The ledger already binds the `w1_profile_d{i}d{j}` 2nd-jet atoms (via `bind_additional_inherited`,
≈L353–360); the `m1_profile_d{i}d{j}` atoms are created **locally** in `operator_dx` (not ledgered), and there
is **no** abstract-space `bg` 2nd-jet table and no 3rd-jet atoms anywhere. New scaffolding is required.

Ground truth: `S11CB_ENERGY_BASIS_COUNT` in the committed `.out` is `Integer(26)` (both anchorings) — the
**frozen** value and the sole public regression target of §4.

## 3 · The object to realise — one live total-derivative, spurion as background data

Build a single **live total in-plane derivative** `D_i` that differentiates the full background-profile jet
tower — each spatial jet mapping to the next (`W_bg→∂W_bg→∂²W_bg→∂³W_bg→…`, μ likewise), every level carried
at background-amplitude order `σ_W¹` — and route every one of the four frozen consumers (A basis quotient, B
LAB_HELD strong rows, C MATERIAL pullback, D coupling `operator_dx`/divergence/curl) through it. Concretely:

1. **The spurion is background DATA the divergence differentiates, NEVER a variational field.** In the §3a
   quotient the spurion `bg` must appear in the **total-divergence / derivative map** (so `D_i bg_j = H_ij`
   generates the 2nd jet) **and must not** be added to the set of fields varied in the Euler–Lagrange
   signature (the L961–967 `δ/δ(field)` loop). Adding `bg` to `basis_fields` would do both and wrongly vary
   a fixed background — do not. Provide the spurion's 2nd jet as a **symmetric** table `H_ij = H_ji` of
   independent knobs (mirroring the DOF `basis_second` tables), in the divergence map only.
2. **Extend the jet tower to the depth each consumer needs, as shared canonical tables.** Define shared,
   ledgered `W`/`μ` jet tables so the same 2nd- and 3rd-jet atoms are used everywhere (fix the current
   split where `w1_profile_dij` is ledgered but `m1_profile_dij` is created locally, and add the 3rd-jet
   atoms the coupling cascade needs). The quotient (A) and strong rows (B, C) reach one divergence deep
   (2nd jet); the coupling cascade (D) reaches two deep (3rd jet). The tower is deep enough when **extending
   it one further level changes no emitted object** (make that a control, §5).
3. **Apply the retained-grade projection consistently.** The live `D_i` on a real invariant also produces
   `σ_W²`/`σ_W³` cross-terms (products of background factors); these are **not** retained and must be removed
   by the retained-grade projection `first_shape_series` (L713–725, keeps `η≤1 ∧ σ_W≤1`). Emit and compare
   physical densities/rows at the **grade-projected** level so a single background factor differentiated many
   times survives while products are dropped; do not leave un-projected `σ_W≥2` terms in an emitted row.

⛔ Do not hard-code any count, invariant, or operator term. The corrected count and corrected rows are
computed outputs of the live derivative; if you type a specific surviving-invariant count or a specific added
term, stop — that is the defect (rule 2 corollary 1). The only by-hand symbol combination is the
action/energy construction (already built) and the ansatz; every operator expression is reached by
computation. The UNIFORM candidate family (`UNIFORM_CANDIDATES`, L1093–1098) carries no spurion and must keep
computing exactly what it computes now.

## 4 · The one acceptance you may check against — the FROZEN-limit public regression

Provide a **dedicated Hessian-freeze switch** that turns the newly-live jet tower **off at the 2nd jet and
above while leaving the first jet live** (equivalently, restores A/B/C/D to the constant-coefficient
divergence). ⛔ This is **not** `uniformize`/the UNIFORM control: `uniformize` (L2511) zeroes the entire first
jet (`σ_W`, `grad_W`, `grad_mu`, `w1_grad`, `m1_grad`), which deletes every spurion invariant and cannot
reproduce `26`. Keep UNIFORM as the full uniform-limit check against S11b; add a separate switch for this
regression.

Under the Hessian-freeze switch, and only under it, at the **grade-projected** level:
- `S11CB_ENERGY_BASIS_COUNT` = the committed public `26` (both anchorings). Emit the frozen-switch count and
  its symbolic residual against `26` (print the residual; do not `assert` it before emitting — §6).
- The grade-projected LAB_HELD and MATERIAL strong rows = the **grade-projected** committed rows. Emit the
  per-row residual (grade-projected live-frozen vs grade-projected committed). A raw byte comparison against
  the committed rows may be emitted only as a **separate diagnostic**, never as the regression (the live
  path introduces `σ_W≥2` terms the projection removes, so raw bytes will not match and that is expected).

⛔⛔ The frozen `26` is the **only** target value in this directive. The corrected (tower-live) basis count,
its per-source count, and the corrected operator structure are **withheld on purpose** and are **not**
acceptance criteria. Do not tune the repair toward any live count — emit the live objects as computed; the
comparison against the reference happens outside this build, where a disagreement is a finding, not a build
failure. A live count that merely differs from `26` but is wrong will be caught by the build's review legs
and our diff — that safety works only if you do not aim at a hidden number.

## 5 · Controls and independent routes (must remain able to FAIL)

- **Form ablation (decisive), at the grade-projected level.** The Hessian-freeze switch of §4 is a **form**
  ablation of the new physics: it must **move** the grade-projected basis count (back to `26`) and move the
  grade-projected rows (back to frozen). A count byte-identical with the tower live vs frozen means the
  repair did nothing (rule 14 / form-ablation first catch). Because the ablation and the grade projection are
  separate operations, project first, then ablate, so the ablation isolates the Hessian **form**, not the
  truncation.
- **Coefficient control (must NOT move the physics).** Rescaling a live jet-atom magnitude by a nonzero
  constant is arithmetic; the basis **count** and **rank** must be invariant while individual coefficients
  scale. Emit enough to distinguish this from the form ablation.
- **Full-source jet-zero ablation (collapse / monotonicity).** ⚠ The existing `ablated_jets` (L1218) zeroes
  only **one direction of one source**, leaving a nonzero multi-component spurion and a nonempty enumerated
  family — it does **not** collapse the family. Add a switch that zeroes **all** components of a source's
  first jet **and** its whole generated jet tower, and check the monotonic inclusion: full-first-jet-zero ⊆
  Hessian-freeze ⊆ tower-live (the strongest ablation empties the new-invariant family for that source).
- **Tower-depth ablation.** Extending the jet tower one level deeper than chosen must change **no** emitted
  object (proves the chosen depth is sufficient); truncating it one level shallower must change an emitted
  object (proves it is necessary).
- **Existing §5 controls** (`CONTROL_TASKS` L39: REP_INVARIANCE, INDEPENDENCE, FORM, UNIFORM, HOMOGENEITY)
  must still run and pass; ⛔ do not weaken any to accommodate the repair, and do not repurpose UNIFORM as the
  frozen-Hessian regression.

## 6 · Method and script obligations (non-negotiable)

1. The script **PRINTS computed objects; it does not state conclusions.** Every emit payload is a CAS object,
   never a prose verdict.
2. **PRINT the residual; do not `assert` it** before emitting. Compute → emit → only then may a structural
   guard assert. The frozen-limit residual against `26` is emitted, then may be guarded.
3. Interpretation belongs to the step record, not the script.
4. No hard-coded results, no tautological residuals (zero for any input tests nothing — corollary 3), no
   emission conditioned on a payload's value (corollary 4). Tag names name the object, never its value or the
   shape of the answer (corollary 2).
5. Run the full engine to completion; keep **both** anchorings (LAB_HELD, MATERIAL_ADVECTED) and **both**
   density representatives throughout — do not narrow scope to converge. Verify the `.out` and
   `S11c_b_exports.py` are regenerated and non-empty (a real run is tens of thousands of tokens).
6. After building, run `reduction/derived_or_declared.py` and `reduction/engine_output_checks.py` on the
   deliverable/output and report their literal verdicts.

## 7 · Builder report (return these)

- The unified diff of every function changed, with L-numbers.
- The literal live `S11CB_ENERGY_BASIS_COUNT` (both anchorings) and per-source new-invariant count —
  whatever they are — as computed objects.
- The Hessian-freeze-switch count and its literal residual against `26`; the grade-projected per-row
  residuals for **both** LAB_HELD and MATERIAL rows.
- Form-ablation vs coefficient-control behaviour; the full-source jet-zero collapse; the tower-depth ablation
  (deeper = no change, shallower = change); the chosen tower depth per consumer.
- The `derived_or_declared.py` and `engine_output_checks.py` literal output.
- Confirmation the full engine ran to completion and both `.out` and `S11c_b_exports.py` were regenerated,
  with byte sizes.

## 8 · Supplied vs computed

SUPPLIED and unfalsifiable within this build: the spec obligation (§1), the four code-verified frozen
consumers and the partial `operator_dx` pattern (§2), the frozen-limit public value `26` (§4). COMPUTED and
subject to external review: the corrected basis count, corrected new-invariant set, corrected LAB_HELD and
MATERIAL operators, and every residual. Do not present a passing frozen-limit regression as confirmation of
the corrected physics — it confirms only that the frozen limit is unchanged.
