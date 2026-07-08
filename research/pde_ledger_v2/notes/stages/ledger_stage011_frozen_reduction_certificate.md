# ledger_stage011_frozen_reduction_certificate

## Status

**Part II — Gravity. II-G1a (build-order 011).** Reshape **plus a de-rig** of the **frozen-reduction leg** of gate
**`pathA_30`** (the gate splits into stages 011 + 012 per the finalized Part-II split). Source top-line, verbatim:
**`DN_UNITTEST_BC_DEPENDENT`** — a JOINT verdict. **This stage carries the reduction-certificate component**: it shows
that, on the postulated straight finite throat with the wall FROZEN, the reduced longitudinal operator collapses to the
constant-coefficient Helmholtz resonator operator — *produced* by a per-term reduction, not typed — together with the
reduction certificate (the validity window). The D/N determinant → DtN `−(ω/c_S)tan(L0ω/c_S)`, the half-shifted pole
ladder, the Robin counterfactual, and the `BC_DEPENDENT` landing (`bc_derivation_emitted=False`, a banked calibration
input) are **stage 012's**; the cut is drawn at the general-solution `dsolve` (never called here).

- **POSTULATED (labeled):** the straight finite throat, brane/mouth at `s=0`, regular pinch-off cap at `s=L0`
  (`R0(L0)=0`), uniform centerline density `ρ0(s)=ρ*`, `A_M0=0`, frozen wall `η=0`; `L0` (throat depth) is
  `ACTION`-geometry (like stage009's slab spacing `d`, NOT a medium constant), and the reference throat-radius taper
  `R0(s)` used to LOCATE the cap is postulated.
- **CITED (integrity-checked, no re-derivation):** the sound-speed law `c_S² = 5Kρ⁴/m_GNLS` = Part I edge R1
  (`ledger_stage005`), evaluated at the throat centerline density `ρ*`; the EOS exponent-5 closure `P=Kρ⁵` stays
  `IMPOSED`; stage004's `{L,T,M}` dictionary + the `[K]=[P]−5[ρ]` chain; `ħ,m` = stage004 ACTION primitives;
  `ξ=ħ/(m c_s)` DERIVED (R2 family).
- **EARNED (the de-rig — this stage's central burden):** the reduced operator `L_s` is ASSEMBLED from the Helmholtz
  core plus the Σ of computed intruding terms, so `operator_is_helmholtz` can genuinely FAIL; the wave speed is
  EXTRACTED from the assembled operator's ψ-coefficient (bulk, not renormalized); the domain cap endpoint is SOLVED
  from the pinch-off; `unsuppressed_operator_intrusion` is COMPUTED; the c_S² dimensional leg + the corrupt-`[K]`
  probe; the ξ≠ℓ_c firewall.

Ledger-local earned-label (NOT a source verdict token): `FROZEN_REDUCTION_HELMHOLTZ_CERTIFIED`. The 011-scoped verdict
`REDUCTION_CERTIFIED` composes into the joint `DN_UNITTEST_BC_DEPENDENT = (011: REDUCTION_CERTIFIED) ∧ (012: D/N ladder
+ BC_DEPENDENT landing)`.

## Purpose

Before the D/N boundary problem can be posed (stage 012), the interior operator on the frozen throat must be pinned.
This stage certifies that pinning: with the wall frozen (`R=R0, η=0`) and only the matter field dynamic
(`∝ e^{−iωt}`), the longitudinal fluctuation on the straight uniform centerline reduces to a **constant-coefficient
Helmholtz resonator** `ψ''(s) + (ω/c_S)²ψ(s) = 0` on `s ∈ [0,L0]`. The source pair asserted this via three tautologies
(an `X≡X` operator equality, a literal `speed_is_cs=True`, and `L0==L0`); the reshape's burden is to **produce** each of
those facts from the reduction so each can fail. The result is banked forward to stage 013 (harmonic β lift) and stage
017 (calibration input).

## The derivation (both engines, own routes)

- **The consumed speed (CITED R1 at ρ*).** `c_S² = 5Kρ*⁴/m` is the barotropic sound speed of the imposed EOS
  `P=Kρ⁵` (`c_s² = (1/m)dP/dρ = 5Kρ⁴/m`), evaluated at the throat centerline reference density `ρ*`. It is cited from
  `ledger_stage005` (edge R1), integrity-checked dual-site: site A the literal `5Kρ⁴/m`, site B the EOS route
  `∂_ρ(Kρ⁵)/m` (both `= 5Kρ⁴/m`, two independent constructions), plus an explicit frozen-export anchor
  `consumed − 5Kρ*⁴/m ≡ 0`. `k = ω/c_S`. (pathA_30's bare `m` is the same stage004 `m_GNLS` primitive.)

- **The assembled operator (the de-rig).** The reduction produces a dimensionally-consistent operator
  `L_s = ψ'' + M·ψ' + N·ψ − B·(ħ²/4m²c_S²)·ψ''''`, where:
  - `M = d/ds(log√γ0)` is the projection-measure first-derivative coefficient (`[M]=L⁻¹`);
  - `N = (ω/c_s(s))²` is the ψ-coefficient with `c_s(s)² = 5Kρ0(s)⁴/m` (`[N]=L⁻²`);
  - `B ∈ {0,1}` is the Bogoliubov inclusion flag — the `k⁴` correction is a genuine **fourth-derivative** intrusion,
    from the identity `ψ'' + (ω/c_S)²ψ − (ħ²/4m²c_S²)ψ'''' = 0` for the full BdG dispersion `ω² = c_s²k² + ħ²k⁴/4m²`
    (verified on `ψ=e^{iks}`); retaining it makes `L_s` fourth-order, NOT Helmholtz.
  `operator_is_helmholtz := expr_equal(L_s − (ψ'' + (ω/c_S)²ψ), 0)` is computed on the ASSEMBLED `L_s` (not two
  hand-typed twins), and `unsuppressed_operator_intrusion := (M≢0) ∨ (N−(ω/c_S)²≢0) ∨ (retained B) ∨ (nonzero
  {δV_conf, ∇Q} witness)` is computed (not the source's hardwired `False`) and gates the verdict.

- **The reduction certificate (the validity window, EARNED).** Under the frozen background each intruding quantity is
  computed to vanish or is deferred:
  - projection measure: `√γ0(s) = A_perp0` constant → `M = d/ds(log√γ0) = 0` (no `ψ'` term);
  - density gradient: `ρ0(s) = ρ*` constant → `ρ0' = 0`, hence `c_s' = 0` and `N = (ω/c_S)²` — the **bulk** sound
    speed, not wall/healing renormalized;
  - confinement: `δV_conf = 0` in the frozen `η=0` test (a computed intrusion witness); `ℓ_c` INERT;
  - quantum potential: `Q(ρ0) = −(ħ²/2m)·∂ₛ²√ρ0/√ρ0`; with uniform `ρ*`, `∂ₛQ = 0` (computed by differentiation, not
    asserted);
  - BdG: deferred (`B=0`) only under `kξ≪1` (`ξ=ħ/(m c_s)`), with smallness witness
    `ħ²k²/(4m²c_S²) = (kξ/2)²` — carried as a validity condition, NOT set to 0.
  So `L_s` is const-coefficient Helmholtz **conditional on** `{ρ0'/ρ0=0, √γ0 const, δV_conf=0, ∇Q=0, kξ≪1}` (edge
  R26). The equivalent ODE artifact `ψ''(s) + (ω/c_S)²ψ(s) = 0` is emitted.

- **Speed extraction (de-rig, EARNED).** The ψ-coefficient of the assembled `L_s` is extracted and shown to equal
  `(ω/c_S)²` with the bulk `c_S² = 5Kρ*⁴/m`; `csgrad ≡ 0` (no spatial renormalization). `speed_is_cs` is set from
  these computed facts (not literal `True`). A wall/healing renormalization would make the extracted speed differ from
  the bulk value.

- **Domain solve (de-rig, EARNED).** The cap endpoint is solved from the pinch-off regularity: with a postulated
  monotone reference taper `R0(s) = R_mouth·(1 − s/L0)` (`R0(0)=R_mouth>0`, `R0(L0)=0`), `cap_endpoint = solve(R0(s)=0,
  s) = L0` — independent of `R_mouth`, so `R_mouth` cancels. `domain_is_L0 := expr_equal(cap_endpoint − L0, 0)`,
  domain `= [0, L0]`. This is a labeled SELECTION on the postulated taper (not a first-principles derivation of `L0`),
  but strictly better than `L0==L0` (a wrong taper root fails it).

- **The ξ≠ℓ_c firewall (EARNED; edge R27).** The healing length `ξ = ħ/(m c_s)` (the BdG-deferral scale) and the
  confinement length `ℓ_c` (in `V_wall(Σ/ℓ_c)`) are kept as DISTINCT symbols, never substituted for one another —
  both are `[L]`, so a dimensional check alone cannot separate them (analogous to R22's `μ_R≠μ_R⁽⁴⁾`). `ℓ_c` is INERT
  here (`δV_conf=0`).

- **The c_S² dimensional leg + corrupt-`[K]` probe (EARNED).** With `[energy]=(2,1,-2)` (order `L,M,T`),
  four-volume `(4,0,0)`, `[P]=[energy]−[four-volume]=(-2,1,-2)`, `[ρ]=(-4,0,0)`, the chain
  `[K]=[P]−5[ρ]=(18,1,-2)` gives `[c_S²=5Kρ*⁴/m]=(2,0,-2)`. The corrupt-`[K]` probe adds one power of `L` to `[K]` →
  `[c_S²]` flips to `(3,0,-2) ≠ (2,0,-2)` → `DN_UNITTEST_FAIL_DIMENSIONAL` (`mutation_fires=True`), with a two-verdict
  self-ablation (with mutation → FAIL_DIMENSIONAL, without → `REDUCTION_CERTIFIED`, `fail_suppressed=True`).

- **The 011-scoped verdict.** Computed from the 011 rungs — `FAIL_DIMENSIONAL` → `FAIL_OPERATOR_INTRUSION`
  (unsuppressed) → `FAIL_OPERATOR_INTRUSION` (¬helmholtz) → `FAIL_WRONG_SPEED` → `FAIL_WRONG_DOMAIN` → else
  `REDUCTION_CERTIFIED`. The joint `DN_UNITTEST_BC_DEPENDENT` is printed as the composition with stage 012, NOT typed
  as 011-earned.

## Consumed inputs

**Cited — no file reads; genuine DUAL-SITE citation-integrity (two independently-constructed sites, `siteA − siteB ≡
0`, plus an explicit frozen-export anchor `consumed − 5Kρ*⁴/m ≡ 0`):**
- **From `ledger_stage005` (I-2, edge R1):** the sound-speed law `c_S² = 5Kρ⁴/m_GNLS` (EOS exponent-5). Site A the
  literal `5Kρ⁴/m`, site B the EOS route `∂_ρ(Kρ⁵)/m`; evaluated at `ρ*`. Every one-site exponent corruption
  (`5→4`, `5→6`) fails BOTH engines; a coordinated both-site drift escapes `siteA−siteB` but is caught by the
  frozen-export anchor. NOT re-derived. Carry the "R1 at ρ*" caveat; EOS exponent-5 stays `IMPOSED`.

## Exports

- The frozen const-coefficient Helmholtz operator `L_s`, the domain `[0,L0]` (cap `R0(L0)=0`), the speed `c_S`, and
  the reduction-certificate **validity window** `{ρ0'/ρ0=0, √γ0 const, δV_conf=0, ∇Q=0, kξ≪1}` → stage 013 (harmonic
  β lift) + stage 017 (calibration input). (Distinct from stage009/010's *bulk* slab Helmholtz mode → 024/026 — the
  **frozen-throat** longitudinal mode, not the bulk mode.)
- Register: **zero new counted knobs**. `L0` (throat depth) = POSTULATED `ACTION`-geometry (tracked, not counted, like
  stage009's `d`); `ℓ_c` = INERT here (`δV_conf=0`; tracked, not counted); `R_mouth` cancels out of the pinch-off root
  (a construction scale); `ξ=ħ/(m c_s)` is R2-family (`DERIVED`). Two **structural** edges: R26 (the frozen-reduction
  validity record) and R27 (the ξ≠ℓ_c firewall) — both discharge NOTHING (the speed is banked R1, so the stage adds
  no reduction), both able-to-fail.

## Verification

- **Reshape (blueprint §5) — bridge-severing:** stripped the `.py` scratch-YAML/`_sympy_exprs.wl` export, the
  Mathematica-YAML re-read, the `digest_mapping`/yaml helpers, the `mma_exports`/`expression_digest`, and the
  report/feed/YAML writers; and the `.wl`'s `Get[sympyExprFile]`, all `sympy*`-comparison equalities, and the
  `Export`. Both engines standalone, print-only, **zero file I/O**, float-free, ledger idioms; the `.wl` derives the
  certificate, the `L_s` assembly, the speed extraction, the domain solve, and the c_S² dim leg natively (own
  `D`/`Solve`/`Coefficient`/`Simplify`), with an arity self-check and a native trig-basis operator-nature cross-check;
  no messages silenced. **Clean 011/012 cut:** no `dsolve`, no D/N determinant, no DtN/pole/Robin/static, no
  `tan_argument`/`Z00` dim legs (all stage 012).
- **Dual-engine:** SymPy **61 PASS / 0 FAIL** · Mathematica **71 PASS / 0 FAIL** (71 = 61 shared + 8 `.wl` arity
  self-checks + 2 native trig-basis confirmations), both exit 0, CWD-independent; runner transcripts under
  `scripts/output/` + `mathematica/output/`.
- **Tri-review (fresh agents):** arbiter re-run via the runners (both engines, repo root + foreign CWD);
  **`FIDELITY_CLEAN`** (full 011-slice coverage diff, every value hand-re-derived incl. the BdG fourth-order identity,
  `[c_S²]=(2,0,-2)`, the taper root, `bdg_ratio=(kξ/2)²`; all three booleans + the intrusion flag confirmed genuinely
  de-rigged; the 011/012 boundary clean); **`ADVERSARIAL_CLEAN`** (14+ mutant matrix — the de-rig proven able-to-fail
  (operator via M/N/B, speed via renormalization, domain via taper), BdG decisively a real fourth-order term
  (coeff of `Derivative(4)` = −ħ²/4m²c_S², speed unshifted), dual-site R1 real with the anchor catching the
  coordinated-drift escape, `.wl` independent, no 012 leakage).
- **Remediation (3 teeth/label nits) → fresh-agent `REVERIFY_CLEAN`:** the ξ≠ℓ_c firewall tooth was made genuine (it
  had been a vacuous `xi≠xi` predicate → replaced by a real `ℓ_c→ξ` conflation that flips `firewall_ok` to False); a
  dedicated witness-path tooth (nonzero `δV_conf` with `M=0, N=bulk, B=0`) was added so `unsuppressed_operator_intrusion`
  is independently load-bearing (hardwiring the flag to `False` now fails this tooth in both engines — the escape the
  adversarial leg found); and the `.wl` arity self-check labels were corrected (5-arg/7-arg, not 6/8). Tallies rose
  58/68 → **61/71** (the tooth-3b asserts). The re-verify confirmed all three fixes genuine and both escapes closed.
- **Teeth (8+1, all fire):** (1) nonconstant `√γ0` → `M≠0`; (2) retained BdG → `ψ''''` term; (3) nonuniform `ρ0` →
  `N` s-dependent (each → `FAIL_OPERATOR_INTRUSION`); (3b) nonzero `δV_conf` witness → flag-driven
  `FAIL_OPERATOR_INTRUSION`; (4) wall renormalization → `FAIL_WRONG_SPEED`; (5) corrupt taper → `FAIL_WRONG_DOMAIN`;
  (6) `ℓ_c→ξ` conflation → firewall fires; (7) dual-site R1 exponent corruption → both engines (both one-site + the
  anchored coordinated both-site); (8) corrupt-`[K]` → `[c_S²]=(3,0,-2)` → `FAIL_DIMENSIONAL` + self-ablation.

## Provenance

- Source gate: `software/stage1_solver/tools/pathA_30_dn_unit_test_{sympy.py,.wl}` (011 slice; sources unchanged);
  `software/stage1_solver/reports/pathA_30_dn_unit_test.md` (`## Reduced Operator` :3–15 + `## Reduction Certificate`
  :17–24; `## Dimensional Check` `cs_squared_from_EOS` leg :70–81).
- Reshape directive + tri-review artifacts: `research/pde_ledger_v2/_scratch/ledger_stage011_*` +
  `_scratch/adv_stage011/` + `_scratch/reverify_stage011/`. Running-start source map:
  `research/pde_ledger_v2/notes/stage011_pathA30_frozen_reduction_source_map.md`.
- Split row: `research/pde_ledger_v2/notes/part2_gravity_atomic_split.md` (id 011). Opens the pathA_30 fold; stage 012
  completes it (the D/N pole ladder + Robin falsifier + the `BC_DEPENDENT` landing).
