# S11c-d clean-condition packet — round-1 review disposition

**Packet (v1, orchestrator-written):**
- A: `directives/S11c_d_clean_condition.md`
- B: `directives/S11c_d_zinvariant_operator_blocks_directive.md`

**Review prompt:** identical for both legs, `directives/_legs/S11c_d_clean_condition_review_prompt.md` (6909 bytes).

**Legs (G1: orchestrator-written → Codex + Grok):**

```
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access "$(<prompt)" < /dev/null > <log> 2>&1   # exit 0
grok --prompt-file <prompt> --cwd /var/projects/toy_physics --model grok-4.6 --effort high \
  --permission-mode bypassPermissions --output-format plain > <log> 2>&1                                                # exit 0
```

**Reports:**
- `S11c_d_clean_condition_review_codex_sol.txt` (the final message)
- `S11c_d_clean_condition_review_grok.txt`

The legs' scripts and stdout are in `S11c_d_clean_condition_review_scripts/{codex,grok}/`. The Codex leg's
`*.copy.*` engine copies are omitted, since they are verbatim copies of tracked files.

**Verdicts:**
- **Codex:** "Not cleared."
- **Grok:** "H-planar holds on this operator; the packet still has five defects that change what would be computed
  or claimed, plus two citation-use mismatches."

Both legs computed H-planar's parity bookkeeping and H-round's inversion parities independently and confirmed
them (Grok `01`–`04`; Codex `reflection_operator_audit`, `round_toroidal_audit`). ⚠ That is leg evidence about
the term structure. It does not substitute for the Part 1 build.

## Orchestrator verification (G4) — mechanical lookups, run 2026-09-30 from the repo root

```
$ sed -n 151,154p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
The three `τ_I` are independent. Equations of motion are obtained by S11b's method — balance laws, the
binding virtual-displacement rule, variational derivatives with held-fixed fields named, and prescribed
external virtual work — **not** by putting an irreversible response kernel in an ordinary action.
$ sed -n 183,187p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
fresh varying profiles on the anchor coordinate `y`,
ξ ≡ y/L_W ,   W_bg(y) ≡ W̄₀[1+η w₁(ξ)] ,   μ_R,bg(y) ≡ μ̄_R[1+η m₁(ξ)] ,   σ_W ≡ η W̄₀/L_W ,
∂_{yᵢ}W_bg = σ_W ∂_{ξᵢ}w₁ ,   ∂_{yᵢ}μ_R,bg = (μ̄_R/W̄₀) σ_W ∂_{ξᵢ}m₁ ;
$ grep -n '^DIRECTIONS =' research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
57:DIRECTIONS = range(3)
$ sed -n 145,148p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
J_s = Λ_A(ω)𝒜_s + Λ_V(ω)V_s ,   Λ_I(ω)=Λ_I⁰/(1−iωτ_I) ,  I∈{A,V,X} ,
𝒜_s = μ_s − δp_s/ρ_m ,   μ_s = μ_θ/ρ_br⁰ ,
t_s = −(δp_s + Λ_X(ω)𝒜_s)n̂_s ,   n̂_s·v_bulk,s = V_s + J_s/ρ_m ,
∂_tΣ + ∇_x·(Σ v) = −(J₊+J₋) ,   Σ ≡ Σ_E ≡ ρ_4D W ,   v ≡ ∂_t u ,
$ sed -n 350,354p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
J_s^α ≡ ρ_m (v_bulk,s − v_face,s^α)·n̂_s^α ,
n̂_s^α·v_bulk,s = V_s^α + J_s^α/ρ_m ,
𝒜_s^α ≡ μ_s^α − δp_s^α/ρ_m ,
t_s^α ≡ −(δp_s^α + Λ_X(ω) 𝒜_s^α)n̂_s^α ,
J_s^α − Λ_A𝒜_s^α − Λ_VV_s^α = 0 .
$ sed -n 315p docs/toy_model_ontology_summary.md   (excerpt)
The fields h_+ and h_- describe the interfaces only while each remains a single-valued graph over the brane
coordinates. They are useful background, far-field, and separated-interface collective coordinates, not a
complete nonlinear throat topology. …
$ sed -n 957p docs/toy_model_ontology_summary.md   (excerpt)
… The support mode is a spectrally normalizable bound state or acceptably long-lived resonance of the complete
variable-coefficient transverse operator; …
$ sed -n 2377p docs/native_light_em_and_vortex_throat_interpretation.md
### 9.8 A trapped standing wave is not automatically spinning
$ sed -n 8,12p AGENTS.md
For future S11c Python constructors, validators and export jobs, use
`scripts/s11c_guarded_run.py` around the existing supervisor. It requires a
host systemd user manager, verifies a whole-job 2 GiB memory cap, disables
job swap, limits process count and CPU affinity, lowers CPU/I/O priority,
records resource samples and fails closed. Keep native memory limits too. Never fall back to an unguarded launch if containment fails.
$ ls -la scripts/s11c_guarded_run.py
-rw-rw-r-- 1 trevnorris trevnorris 10714 Sep 30 09:51 scripts/s11c_guarded_run.py
$ sed -n 3p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_flow_calibration_assessment.md   (excerpt)
The current pilot answers a strict rest-bulk, LAB_HELD/RHO4_CONSTANT development-input question. It does not yet
answer leakage in the calibrated, draining medium. …
```

## Dispositions

Every finding below changes what is computed or what may be claimed. **All are accepted.** None is rejected.

| # | finding (legs) | verified by | disposition |
|---|---|---|---|
| 1 | Background jets exist in all three in-plane directions, and `y` is S11c's anchor-coordinate *vector*. "Perturbations independent of `z`" is not an invariant class unless every jet with a 2- or 3-index is zeroed (Codex 1; Grok 4). A's rotation claim holds per tangential Fourier component on a one-direction background only (Codex 1). | `DIRECTIONS = range(3)`; spec `:183–187` | **B:** state the restriction by engine index before construction: every background jet and every background datum carrying index 2 or 3 is zero. **A:** scope the rotation claim. |
| 2 | K1/K2 are under-specified and admit null representatives (Codex 5; Grok 1, 2). | the legs' scripts (Grok `03`; Codex `reflection_operator_audit`) | **B:** pin both terms as equations with live coefficients. K1 contracts `ê` into a field index; K2 is a single-ε term that is nonzero on the class. |
| 3 | "Governing action" would put the `Λ` kernels in an ordinary action (Grok 3). | spec `:151–154` | **B:** construct by the spec §1c method. The structural rule names the stored energy, the supplied face laws and the ansatz. Controls enter the stored energy only. |
| 4 | Part 1 is not a closed brane–bulk operator: there is no `φ`/`δp` row, and the bulk enters as face operands only (Codex 2). | spec `:145–148`; S11c-b performs no curved-bulk solve | **B:** print every bulk-facing operand and face kinematic quantity (`δp_±`, `n̂·v_bulk`, `J_s`, `V_s`, `𝒜_s`) as explicit rows/columns, with their dependence on every field. **A:** Part 1 establishes the slab/face block structure. "No bulk leakage" additionally needs the closure argument, stated as an interpretation for review, and the drain-on operator. |
| 5 | P2 is the wrong premise: the face laws carry vector operands (Codex 3; Grok citation). | S11c-a `:350–354` | **A:** P2 becomes equivariance of every scalar, polar-vector and axial-vector operand and of the support/boundary data. Spin/director fields are added to it. The citation is corrected. |
| 6 | The throat coverage is overclaimed. H-round covers a hypothetical `O(3)` background with a normalizable toroidal mode, not a represented throat (Codex 4). | ontology `:315`, `:957` | **A:** the coverage becomes conditional. New falsifiers: no normalizable toroidal mode; full throat not `O(3)`; extra parent/core/support fields. |
| 7 | B's model point omits the 40-term basis coefficients (Codex 6). | record `:35` ("40 = 10 uniform + 15 … + 15") | **B:** all 40 accepted basis coefficients stay live, plus both bookkeepers, both density representatives, both anchorings, and every allowed higher jet subject to #1. |
| 8 | Feasibility: AGENTS.md requires the guarded runner (2 GiB cap, no unguarded fallback). The 20 GB clause contradicts it. There is no blockwise cross-engine residual (Codex 7). | AGENTS.md `:8–12`; runner exists | **B:** every engine and comparator job runs under `scripts/s11c_guarded_run.py`. A contained failure is a reported result. Specialize the background before construction. Add a blockwise comparator with the whole-row sign conventions mapped explicitly. |
| 9 | F2: `m ≠ 0` alone does not carry angular momentum (a standing `cos mφ` carries none). An isotropic microrotation/director field could enlarge the field content without breaking `O(3)` (Codex 8). | native_light `:2377` heading | **A:** F2 is reworded (a circulating complex combination, not `m ≠ 0`). New falsifier for the extra-field case. |
| 10 | "Transfers to the calibrated model" is an overclaim (Codex 9; Grok 5). | assessment `:3` | **A:** the theorem transfers *if* the completed calibrated, draining operator has the same symmetry and field content. That is not yet shown. |
| 11 | §4 wording: fast bulk is "suppressed, not exact"; the flow horizon is unverified while there is no convective operator; "no clean zero" should be "no generic symmetry-protected zero" (Codex 9). | the legs' `speed_gap_scattering_audit` | **A:** reworded. |
| 12 | Citations: the `v=0` pinpoint quotes the grazing clause; P2's pointer; Molz–Beamish wording too strong; "λγ = 1" is a calibrated, uncommitted target; B's citations are missing from the lookups file (Codex; Grok). | lookups file + files | **A/B:** corrected. The lookups file is extended. |

**Next:** fold into v2 (A, B), then a fresh two-leg round. This is a physics spec/directive, so the rule is
review-until-clear (G2/G4), ⛔ not fold-once.

---

# Round 2 — review of v2

**Prompt:** `directives/_legs/S11c_d_clean_condition_review_round2_prompt.md` (7944 bytes). It is identical for
both legs and is rendered from the round-1 template. It adds a fold-check item and hands over this disposition.
The commands are the same as in round 1.

**Reports:**
- `S11c_d_clean_condition_review_r2_codex_sol.txt` (final message)
- `S11c_d_clean_condition_review_r2_grok.txt`

The scripts are in `S11c_d_clean_condition_review_scripts/r2_{codex,grok}/`.

**Verdicts:**
- **Codex:** "not cleared". The fixed-mirror H-planar, pin B, both anchorings, K1/K2 as FORM controls and the H-round
  parity table all clear.
- **Grok:** "H is right for this operator on the class A/B now name. The packet is **not clear**. Two defects still
  change what Part 1 would compute or what the citations support."

**⚠ Process slip, recorded.** I began the v3 fold by editing v2 in place before committing v2 as the reviewed
baseline (G4). v2 was recovered exactly:
- **Source:** the Codex leg's literal `nl -ba` dump of both files at review time (transcript line 11011 onward;
  the transcript is outside the repo).
- **Check:** byte-identical (`cmp`) to the copies the Grok leg saved independently (`/tmp/s11cd_clean_review_r2_grok/{A,B}.md`).

The v2 files are therefore committed as reviewed, and v3 is applied after this commit.

## Orchestrator verification — mechanical lookups, 2026-09-30

```
$ sed -n 112,114p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
  in-band-but-verified-out-of-band: the #89 PY skipped controls + the #89b WL heavy controls.) **Two whole-row SIGN
  CONVENTIONS** to adjudicate there (kinetic −K PY vs +K WL; face generalized-force PY `+diff` vs WL
  `−linearVirtualVariation`) — the comparator SURFACES them, does not normalize them (rule 1/6). **#90's two flags**
$ sed -n 35,37p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
density representatives ρ4D/ρbr, two anchorings LAB_HELD/MATERIAL_ADVECTED), S11c-b computes: (1) the **§3a energy
basis** — the O(3)-Kronecker field-bilinear invariant family, corrected to **40 = 10 uniform + 15 ∂W_bg-spurion + 15
∂μ_R,bg-spurion**; (2) the **variable-coefficient slab OPERATOR ROWS** — the equations of motion for U/θ/e_W with the
$ sed -n 2377,2381p docs/native_light_em_and_vortex_throat_interpretation.md
### 9.8 A trapped standing wave is not automatically spinning

A real linearly polarized standing wave can have zero time-averaged angular
momentum. Two degenerate modes with a relative phase can form a circularly
polarized bound pattern that carries angular momentum. Orbital winding can
$ grep -n "CENTER_FACE_GENERALIZED_ROW" research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | head -2
2237:        "CENTER_FACE_GENERALIZED_ROW": center_face,
3112:            Str("CENTER_FACE_GENERALIZED_ROW"),
$ sed -n 343,345p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md ; sed -n 365,366p (same file)
δ_vx_s^α ≡ δ_vR_s^α|_X ,              v_face,s^α ≡ ∂_tR_s^α|_X ,
V_s^α ≡ V_{n,s}^α ≡ n̂_s^α·v_face,s^α .
δ_v𝒲_bulk^α ≡ Σ_s a_s^α t_s^α·δ_vx_s^α ,
∂_tΣ^α + ∇_x·(Σ^α v) = −Σ_s a_s^α J_s^α ,       v = ∂_t u .
$ sed -n 195,196p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
𝔅⁰ ≡ {W_bg, μ_R,bg, ρ_4D,bg⁰, ρ_br,bg⁰, θ⁰, V_s⁰, J_s⁰, 𝒜_s⁰, boundary loads} ,   θ⁰ ≡ 0 ,
𝒮_hold⁰ ≡ {f_hold⁰(x), t_hold,s⁰(x)} ,   V_s⁰ = J_s⁰ = 𝒜_s⁰ = 0 .
$ sed -n 190,195p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    for i in range(1, 4)
)
w1_grad = tuple(
    inherited_symbol(f"w1_profile_d{i}", "KNOB", f"thickness-profile first jet {i}", DIM_ZERO)
    for i in range(1, 4)
)
$ grep -o 'w1_profile_d[0-9]' research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c
  43040 w1_profile_d1
  34147 w1_profile_d2
  27186 w1_profile_d3
$ ls CHARTER.md ../../CHARTER.md   (from research/pde_ledger_v3)
ls: cannot access '../../CHARTER.md': No such file or directory
CHARTER.md
```

## Dispositions (all accepted; v3 folds them)

| # | finding (leg) | verified by | v3 change |
|---|---|---|---|
| R2-1 | "Index 2 or 3" is ambiguous. SymPy loops `DIRECTIONS = range(3)` (0-based), but symbol names are 1-based (`u_1`, `w1_profile_d1`) (Grok 1). | engine `:190–195`; export symbol counts | **B:** R1 is bound to symbol-name labels, and the 0-based loop variable is excluded explicitly. |
| R2-2 | The rotation-to-every-incidence claim needs **full `O(2)` invariance** of background values about direction 1. Independence of directions 2/3 plus one mirror is not enough: a constant direction-2 vector breaks it (Codex 1). | spec `:195–196` (hold vectors exist) | **A:** an explicit `O(2)` premise, with the fallback that only incidence planes containing a mirror are protected. P2 separates law covariance from background-value invariance. **B:** R1 requires background vectors along direction 1 only, and tensors rotation-invariant about it. |
| R2-3 | The comparator's pre-residual convention mapping conflicts with the record: "the comparator SURFACES them, does not normalize them" (Codex 2). | record `:112–114` | **B:** raw operands and raw residual first; a mapped diagnostic only after, under distinct names. |
| R2-4 | The operator/dependency matrix is not uniquely defined. It omits the `ζ_c` centre row, `t_s`, `v_face`, `a_s` and `δ_v x_s`; `δp_s` is an input (Codex 3). | engine `:2237`; S11c-a `:343–354`, `:365–366` | **B:** one rectangular Fréchet map with ordered inputs (slab fields + bulk trace inputs) and ordered outputs (all rows including the centre row + every face quantity). **A:** the closure paragraph lists the full face path. |
| R2-5 | The SymPy import of `S11c_b_exports.py` bypasses "restrict before construction" (Codex 4). | B v2 `§4` vs `§2` | **B:** an import whitelist (symbol definitions + the accepted 40-term basis only); no exported operator/kernel/face row. |
| R2-6 | The F2 angular-momentum wording is wrong. Two real standing modes in quadrature carry angular momentum; "circulating" need not mean travelling (Codex 5). **New in v2.** | native_light `:2377–2381` | **A:** F2 quotes the source and names the degenerate-pair-in-quadrature case. |
| R2-7 | R-LEAK-1 mixes a linear selection rule with particle stability, so F3/F4 cannot falsify a linear theorem (Codex 6). | Codex `operator_symmetry_audit` stdout: `LINEAR_ODD_SCALAR_HESSIAN_AT_BACKGROUND 0`, `NONLINEAR_SCALAR_SOURCE lambda*odd_amp**2` | **A:** R-LEAK-1 is the particle-stability requirement. H is its linear clean condition, with falsifiers F1/F2/F2b/F5/F6. F3/F4 become nonlinear gates N1/N2. H alone is stated not to deliver stability. |
| R2-8 | "Exponential for smooth profiles" is wrong for C∞ compact bumps, which have stretched-exponential tails (Codex 7). | Codex `smooth_profile_fourier_audit` stdout | **A:** "rapid suppression (exponential for suitable analytic profiles)". |
| R2-9 | Citations: record `:35` → `:35–37`; `CHARTER.md` path; leg-evidence claims grounded only by a directory listing; §3c pairing text not reproduced (Grok 2; Codex). | lookups above | **A/B:** corrected. The grounding files reproduce the leg stdout excerpts and the §3c text. |

**Round-count note.** v2 introduced two new defects in material just changed (R2-3 conflicts with the record; R2-6
mis-states the source). Under G4, if v3's review again finds defects bred by the fold itself, the author changes
(⛔ no fourth fold). The physics core (H, its parity bookkeeping) cleared in both rounds. The open items are the
directive's mechanics and the doc's claim scope.

---

# Round 3 — review of v3

**Prompt:** `directives/_legs/S11c_d_clean_condition_review_round3_prompt.md`. It is identical for both legs and
is rendered from round 2 with the version, the fold item (R2-1 … R2-9) and the sandbox path updated. The commands
are the same as in round 2 (`codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access
… < /dev/null`; `grok --model grok-4.6 --effort high`).

**Reports:**
- `S11c_d_clean_condition_review_r3_codex_sol.txt` (final message)
- `S11c_d_clean_condition_review_r3_grok.txt`

The scripts and literal stdout are in `S11c_d_clean_condition_review_scripts/r3_{codex,grok}/`.

**Verdicts:**
- **Codex:** "Not cleared. The symmetry theorem itself is substantially correct, but Part 1's rectangular map is not
  yet the physical or complete map needed to test it."
- **Grok:** "**Not cleared.** H-planar is right for this operator's energy, pin B, both anchorings, and the
  *contracted* face laws on class `P` under R1. Three defects still change what Part 1 would be claimed to show, or
  what K1 would compute."

**Cleared by both, for the third time:**
- H-planar on R1/P (Kronecker basis, pin B, both anchorings, `Λ` laws);
- the `O(2)` rotation about direction 1;
- the H-round parities (ℓ ≥ 1);
- K2 as a FORM control;
- the import whitelist, raw-first comparator, blind Wolfram, freezes, lane, and Part 2 scoping.

## Orchestrator verification — mechanical lookups, 2026-09-30

````
$ sed -n 150,165p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py   (excerpt)
# Traced bulk perturbations and their normal jets at the flat reference faces.
dw_delta_v_bulk = { s: tuple(symbol(f"d_w_delta_v_bulk_{...}_{i}", "COORDINATE", ...) for i in range(1, 5)) ... }
dw_delta_p = { 1: symbol("d_w_delta_p_plus", ...), -1: symbol("d_w_delta_p_minus", ...) }
$ sed -n 636,646p (same file)   → pressure/velocity traces built by affine_bulk_perturbation(trace, d_w jet, face)
$ grep -o -E 'd_w_delta_p_(plus|minus)' research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c
   1477 d_w_delta_p_minus
   1477 d_w_delta_p_plus
$ grep -o -E 'd_w_delta_v_bulk_(plus|minus)_[0-9w]' research/pde_ledger_v3/scripts/S11c_b_exports.py | sort | uniq -c
    221 (each of the eight, plus/minus × 1–4)
$ sed -n 95p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
Inherited from S11c-a §1b unchanged: the rest-frame bulk fields `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`,
$ sed -n 37,38p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
- `u(x,t)` is a three-component in-plane displacement and has no `w` component. The material map is
  `x(X,t) = X + u(X,t)` with `𝒥_x = det(∂x/∂X)`.
$ sed -n 340,344p research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
For each branch, the only permitted virtual face displacement and face velocity vector are obtained from
the supplied parameterization:
δ_vx_s^α ≡ δ_vR_s^α|_X ,              v_face,s^α ≡ ∂_tR_s^α|_X ,
$ sed -n 214,217p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
u = tuple(
    inherited_symbol(f"u_{a}", "COORDINATE", f"displacement component {a}", DIM_L)
    for a in range(1, 4)
)
$ sed -n 14,17p research/pde_ledger_v3/CHARTER.md
- ⭐⭐ **v3 then TOOK a method change (2026-08-01, user decision): REQUIREMENTS-FIRST.** Each force
  sector states what it needs to survive. **Brane and bulk are defined LAST, at the knit**, by asking
  whether one medium can satisfy every requirement at once. **A no-go between requirements IS the
  falsification.**
$ sed -n 362p docs/toy_model_ontology_summary.md
- **Structural support:** energy \(E_{\rm support}\) in a trapped transverse brane-shear standing mode—the model's light or photon mode—helps hold the aperture open.
$ sed -n 100p docs/toy_model_ontology_summary.md   (excerpt)
... distributed return transfers material back into the ordered state. The bulk may store the corresponding response as compression, flow, internal energy, entropy, or other unresolved excitations of the same medium. ...
$ sed -n 151,152p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
(blank)
The three `τ_I` are independent. Equations of motion are obtained by S11b's method — balance laws, the
$ sed -n 20p research/pde_ledger_v3/steps/S11c_b_variable_coefficient_operator.md
> exceeded 30 GB (one case > 15.6 GB and growing → OOM; STEP 0's on-box claim was overturned —
$ sed -n 12p research/pde_ledger_v3/directives/_measurements/S11c_b_step0_residual_scope.md   (excerpt)
`kernelTruncated/Extended`) that dominate the full run's ~16 GB/case peak.
````

The leg stdout relied on below is quoted from the filed `r3_*` files:
- `r3_codex/face_trace_and_map_audit.stdout.txt`:
  - `V_BULK_3_TRACE 0`
  - `INDEPENDENT_FOUR_VELOCITY_COMPONENTS_ON_P False`
  - `D_DELTAV_R3_D_U3 0`
  - `D_VFACE3_D_U3 d_t`
  - `DELTA_V_X_IS_PHYSICAL_OUTPUT_OF_SAME_FRECHET_MAP False`
- `r3_codex/reflection_and_round_parity.stdout.txt`:
  - `REFLECTION_BLOCK vface3 <-u3 ALLOWED 1`
  - `REFLECTION_BLOCK V <-u3 FORBIDDEN -1` (also `J`, `affinity`)
  - `ROUND_PARITY 0 … TOROIDAL_EXISTS False`
  - `R1_LOST_O2_ALLOWED_COMPONENT B`
- `r3_grok/03_K1_K2_form_controls.out`:
  - `K1_MIXED_ON_P [('d1u3', 'd1theta', a_K1*cos(beta)), ('d2u3', 'd2theta', a_K1*cos(beta))]`
  - `K1_MIXED_IF_R1_ZEROS_E3 []`
  - `K2_MIXED_ON_P_R1 [('d2u3', 'theta', a_K2*g1)]`
- `r3_grok/05_O2_R1_tensors_and_rotation.out`:
  - `T_SO2_INVARIANT_ANSATZ … T23 …`
  - `ANTISYMMETRIC_T23_BREAKS_MIRROR True`
  - `AXIAL_MIXING_ON_P A1*d2u3*theta`

## Dispositions (all accepted)

"Fold-bred" marks a defect introduced by the v3 fold in material it changed.

| # | finding (leg) | verified by | resolution owed in v4 |
|---|---|---|---|
| R3-1 | B's input list omits the bulk-trace **normal jets** `∂_wδp_s` and `∂_w v_bulk,s` (4 components each). The engines carry them as live coordinates, entering through the shifted-trace law (Codex 1). Partly fold-bred: R2-4 was not fully resolved. | S11c-a audit `:150–165`, `:636–646`; export counts | B's inputs are the bulk-trace coordinates **as the engines carry them**, enumerated from the code. P2 classifies the jets (`∂_w` is reflection-even). |
| R3-2 | "Four independent `v_bulk,s` components" contradicts the supplied rest-frame bulk `v_bulk=∇₄φ`, `δp=−ρ_m∂_tφ`. On `P`, `v_3 = ∂_3φ = 0` and `v_2` is fixed by `δp`, so an independent odd bulk column is spurious (Codex 2). Fold-bred. | spec `:95`; Codex stdout | B emits the formal operand map **and** its pullback to the supplied potential-flow trace subspace on `P`, each labelled. |
| R3-3 | A's closure criterion "no face quantity depends on the twist-sector field" is **false** for the correct operator. `v_face,3 = ∂_t u_3` and `δ_v x_3 = δ_v u_3` are allowed odd→odd kinematic blocks; the even face scalars `V_s`, `J_s`, `𝒜_s` and `t_s` (with `n_3 = 0`) do not see `u_3`. The criterion is **parity block-diagonality**: no reflection-odd ↔ reflection-even block. F1's "any bulk-facing quantity" has the same defect (Codex 3, Grok 1). Fold-bred: v3 added `v_face`/`δ_v x` to the list. | S11c-a `:37–38`, `:340–344`; both legs' stdout | A §1 and F1 restated as a parity-block criterion. Note: on `P` the rest-frame potential bulk has no odd trace (R3-2), which ties to `native_light…:1343`. |
| R3-4 | `δ_v x_s` is a **virtual/test** quantity. It is not an output of the physical-input Jacobian (`D_{u_3}δ_vR_{s,3}=0` while `D_{u_3}v_face,3=∂_t`) (Codex 4). Fold-bred. | S11c-a `:340–344`; Codex stdout | B types it correctly: a separate virtual-kinematic map, or one weak bilinear form with explicit trial and test domains. |
| R3-5 | R1 is inconsistent. Clause 1 (zero every symbol carrying label 2/3 in **any** position) erases `O(2)`-allowed tensor components (`diag(A,B,B)` → `diag(A,0,0)`). Clause 4 (rotation-invariant only) admits an axial `T_23` that breaks the mirror and mixes `θ`–`∂_2u_3`. A's "R1 requires full `O(2)`" is therefore not what B's R1 says (Codex 5, Grok 3). Absent from this Kronecker constitutive law, so H-planar on this operator stands. Fold-bred. | both legs' stdout | R1 is stated as the **symmetry**: every background datum is invariant under the full `O(2)` (rotations **and** reflections) about direction 1, and every profile depends on `x_1` alone. Label-zeroing applies only to spatial-derivative indices of profile jets. |
| R3-6 | K1's `ê = (sin β, 0, cos β)` is a background vector, so R1 clause 3 zeroes `ê_3`. K1 then has no `u_3`–`θ` block and cannot expose an unimplemented block (Grok 2). Fold-bred: R1.3 was added in v3. K2 is unaffected. | Grok stdout `K1_MIXED_IF_R1_ZEROS_E3 []` | B states that control terms are added **after** R1 and are not subject to it. |
| R3-7 | The toroidal sector is empty at `ℓ = 0` (`T_00 = r×∇Y_00 = 0`) (Codex 6). Pre-existing since v1. | Codex stdout | A: "for `ℓ ≥ 1`". |
| R3-8 | Flow-horizon rejection: "leaked energy returns to the brane, not to the particle" is not established by ontology `:100` (material return ≠ coherent return of mode energy) (Codex 7). Pre-existing. | ontology `:100` | A: "no mechanism has been shown to return that energy to the same particle". |
| R3-9 | F5/F6 are **applicability** tests of R-LEAK-1's linear condition, not falsifiers of the conditional selection rule. F1 is the operator falsifier (Codex 8). Partly fold-bred, from the v3 restructure. | A §3 vs A §1 conditionals | A: relabel. |
| R3-10 | Citations (both legs): B `:190–195`/`:420` for `u_1` → `u_1` is at `:214–217`; CHARTER `:14` → `:14–17`; A's "performs no curved-bulk response solve" is at spec `:95–97` only (`:145–148` supports "through face quantities"); ontology `:362` says "aperture", not "throat"; spec method quote is at `:152–154`. | lookups above | Corrected, and grounding regenerated. |

**Feasibility (both legs; not a packet defect).** Under the mandatory 2 GiB guard, completing Part 1 is
unestablished. The S11c-b full Wolfram run peaked at about 16 GB per case and was OOM-killed at 15.6 GB. B already
treats a contained failure as a result. Whether to run it under the guard, on the remote ≥64 GB plan, or under a
user-changed cap is a **user decision**, surfaced at handover.

## Round-count / author decision (G4, M2)

Successive folds keep breeding defects in the material just changed:
- v2 bred R2-3 and R2-6;
- v3 bred R3-2, R3-3, R3-4, R3-5, R3-6, part of R3-1 and R3-9, and the `u_1` citation.

All three rounds concentrate on B's enumeration of the map and on A's closure wording. The physics core cleared
every time, which is the recipe-creep tell: the HOW keeps being specified from paraphrase, not from the code.
**⇒ Change the author:**
- **v4 is authored by Codex** (`gpt-5.6-sol`, xhigh), working from the committed v3 baseline, all three
  dispositions and the round-3 reports.
- **Authorship becomes mixed:** v1–v3 are orchestrator-written; v4 is Codex-revised.
- **The valid non-author pairing for v4 review is a fresh Claude agent + Grok.** Neither authored any version. ⛔
  The Codex authoring instance does not review.
- ⛔ No fourth orchestrator fold.

---

# Round 4 — review of v4 (Codex-authored)

**Authorship and pairing.** v1–v3 are orchestrator-written; v4 was revised by a Codex author (prompt
`_legs/S11c_d_clean_condition_v4_author_prompt.md`, `gpt-5.6-sol` xhigh, `--sandbox workspace-write`). That is
mixed authorship. The valid non-author pairing is a **fresh Claude agent** (Opus, general-purpose, no prior
contact) **+ Grok** (`grok-4.6`, high). Neither authored any version.

**Prompt:** `_legs/S11c_d_clean_condition_review_round4_prompt.md`. It is identical for both legs; the Claude leg
received the file text verbatim.

**Reports:**
- `S11c_d_clean_condition_review_r4_claude.txt` (the leg's final report, verbatim)
- `S11c_d_clean_condition_review_r4_grok.txt`

The scripts and literal stdout are in `S11c_d_clean_condition_review_scripts/r4_{claude,grok}/`. The Claude leg's
`engine_copy/` and `import_probe/` were not filed (about 47 MB of copies of committed sources).

**Verdicts:**
- **Claude:** "NOT CLEARED. H itself holds on the actual operator. Findings 1–5 change what Directive B would
  compute. Findings 6–9 change what Record A may claim."
- **Grok:** "**cleared** on the physics filter."

**Adjudication:** **not cleared.** Every Claude finding is verified below. Grok's clear is outweighed by verified
defects it did not catch, among them a pairing that prints a vacuous zero (R4-4).

**Strongest evidence to date (review-leg, single engine).** The Claude leg read the **actual exported S11c-b slab
operator**, applied R1/P, and split every leaf by mirror parity: `LEAVES 307 MIXED_PARITY_LEAVES 0` in all four
cases (`LAB_HELD`/`MATERIAL_ADVECTED` × `RHO4`/`RHOBR`). Its untransformed direction-3 profile-jet control lifts
this to 80/74/117/110 (`r4_claude/02_engine_reflection_parity.stdout.txt`). This is SymPy-export evidence from a
review leg, ⛔ not the blind dual-engine build.

## Orchestrator verification — mechanical lookups, 2026-10-01

````
$ sed -n 484,547p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | grep -o -E '"(delta_rho_4D_face|d_w_delta_rho_4D_face|delta_j_bulk|d_w_delta_j_bulk|trace_grad_f|d_w_trace_grad_f|delta_rho_4D_bulk_t)[^"]*"' | sort -u
"d_w_delta_j_bulk_{face_name}_{component}"
"d_w_delta_rho_4D_face_{face_name}"
"d_w_trace_grad_f_{component}"
"delta_j_bulk_{component}"
"delta_rho_4D_bulk_t"
"delta_rho_4D_face_{face_name}"
"trace_grad_f_{component}"
$ sed -n 647,652p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
    density_perturbation = affine_bulk_perturbation(
        delta_rho4_face[face], dw_delta_rho4_face[face], face,
    )
    current_perturbation = tuple(
        affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face)
        for i in range(4)
$ sed -n 584,585p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
    reference_height = sp.Rational(face, 2) * W0
    return reference_value + (w - reference_height) * reference_normal_jet
$ grep -n 'formal operand domain' research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
141:`:578–592,629–646`. These 20 trace coordinates define a **formal operand domain**, not 20 independent physical
$ sed -n 328,329p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
For the in-plane domain, take all trial and test fields to have compact support in its interior, so the
in-plane integration-by-parts boundary term is fixed to zero; the inherited face boundary conditions still
$ grep -n 'compact support' research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
257:before the weak object exists. The pairing domain is class `P` with compact support in the in-plane interior, as
$ sed -n 11,13p research/pde_ledger_v3/directives/S11c_b_p2b_gamma_bridge_directive.md
- `I_PY = W_0·I_WL` (W-family) and `I_PY = μ_R·I_WL` (μ-family): `EXACT_UNIQUE 0 / SCALED_UNIQUE 30`. WL spurion
  `∇W/W_0`, `∇μ/μ_R` (`…audit.wl:721-722`); PY raw jets `grad_W`/`grad_mu` (`sympy_audit.py:182`). ⇒ energy terms
  `γ·I` equal ⟺ `γ_WL = W_0·γ_PY` (resp. `μ_R·γ_PY`).
$ grep -n -c -i 'p2b' research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
0
$ grep -n -E '^\s*print\(' …/S11c_d_clean_condition_v4_author_scripts/potential_trace_pullback_audit.py   (excerpt)
31:print("CLASS_P_ODD_BULK_TRACE_COMPONENT", 0)
32:print("CLASS_P_ODD_BULK_NORMAL_JET_COMPONENT", 0)
44:print("D_VIRTUAL_X3_D_TEST_U3", 1)
$ git log -1 --format='%h %ad %s' --date=iso 14861016
14861016 2026-09-30 18:07:45 -0600 Resume central benchmark and scope equal-speed feasibility
$ sed -n 47p research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_equal_speed_feasibility.md   (excerpt)
The same PID4097233 was resumed with SIGCONT after checking start ticks85113744, exact command/cgroup, 2GiB/zero-swap/task containment, …
$ grep State /proc/4097233/status
State:	R (running)
$ sed -n 116,117p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
symmetry group in full (`S11b_SHARED_PHYSICS.md:280–288`): **in-plane translation invariance** (so `u` enters
only through its gradients, never undifferentiated), **in-plane `O(3)` isotropy and parity**, **reflection
$ sed -n 114,116p docs/native_light_em_and_vortex_throat_interpretation.md
- Compare candidate carriers: circulating intake, brane-tangent vortex flow,
  mixed \(a\)-\(w\) circulation, trapped chiral shear, and an independent
  microrotation of the ordered substructure.
$ grep -n -i -E 'chiral shear|a.w circulation|helic|mean flow' research/pde_ledger_v3/directives/S11c_d_clean_condition.md
(no output)
$ sed -n 8,9p AGENTS.md
For future S11c Python constructors, validators and export jobs, use
`scripts/s11c_guarded_run.py` around the existing supervisor. It requires a
````

The leg stdout relied on below:
- `r4_claude/01_slab_operator_structure.stdout.txt`:
  - `UNDECLARED_INPUT_SYMBOLS ['delta_j_bulk_1', 'delta_j_bulk_2', 'delta_j_bulk_3', 'delta_rho_4D_bulk_t']` (projection_* leaves)
  - `d_w_trace_grad_f_1..4` (conormal_deriv)
  - `d_w_delta_j_bulk_minus_1..` (face_shift)
- `r4_claude/06_face_controls_pullback_pairing.stdout.txt`:
  - `PART_C PAIRING_SAME_CLASS_P 0`
  - `PART_C PAIRING_CONJUGATE 2*I*pi*U_amp*V_amp`
- `r4_claude/09_round_sector_checks.stdout.txt`:
  - `HELICITY_PURE_TOROIDAL 0`
  - `HELICITY_TOROIDAL_PLUS_c_POLOIDAL 5*sqrt(2)*pi**(3/2)*c/4`
  - `CORIOLIS_IMAGE_BREATHING_PROJECTION … -sqrt(2)*pi**(3/2)*Omega/4`
  - `CORIOLIS_IMAGE_L2_PROJECTION … -sqrt(2)*pi**(3/2)*Omega/8`
- `r4_claude/07_author_script_form_ablation.stdout.txt`: every claim tag is byte-identical under the
  rotational-bulk corruption.

## Dispositions (all accepted)

"Fold-bred" marks a defect introduced by v4 in material it changed.

| # | finding (Claude leg) | verified by | resolution owed in v5 |
|---|---|---|---|
| R4-1 | R3-1 is only partly resolved. The codomain depends on **25 more** engine coordinates outside B's 20-coordinate domain. They include the face density traces and their jets, the bulk current `delta_j_bulk_1..4` (no face label) and its face jets, `trace_grad_f_1..4` and its jets, and `delta_rho_4D_bulk_t`; some are odd (`…_3`). The physical-map entries also carry test symbols `delta_v_*`, untyped. No density or current pullback is supplied. | engine `:484–547`; S11c-a `:647–652`; leg `01` | The domain is defined as **every free perturbation coordinate the codomain actually depends on** under R1/P, computed from the constructed objects (⛔ not a hand list). Each is typed (slab / bulk trace / bulk interior / test / background) and declared, or held with a stated reason. Density and current pullbacks are supplied as equations from a cited governing relation. The face assignment of unlabelled coordinates is stated. |
| R4-2 | The pullback's trace height is unfixed. The engine traces are at the **flat reference face** `w = sW_0/2`; a background-face reading double-shifts (the WL c2 defect class). Fold-bred. | S11c-a `:584–585`; leg `06` Part B | B fixes the evaluation height to the engine's reference face. |
| R4-3 | Controls entering only at the stored energy cannot expose an omission in the face-trace or kinematic columns (`dJ/dv_b3`, `dV/du3_t`, face-force entries are unchanged under K1/K2). A fixed, untransformed direction-3 **background** datum does bite (leg `02`, `06`). Pre-existing design gap. | leg `06` Part A; leg `02` control | Add a FORM control that is an untransformed direction-3 background datum introduced after R1, with the structural rule widened accordingly. |
| R4-4 | Class `P` (plane wave in `x_2`, constant in `x_3`) contradicts "compact support in the in-plane interior". A same-class test field gives a **vacuous zero** (`PAIRING_SAME_CLASS_P 0`). Fold-bred. | spec `:328–329`; B `:257`; leg `06` Part C | Test fields carry the conjugate exponential, compact support in `x_1` only, and the pairing is per `x_2` period and per unit `x_3` length. Or an equivalent, stated domain. |
| R4-5 | The comparator has no account of the P2b coefficient normalization (`γ_WL = W_0·γ_PY`, resp. `μ_R·γ_PY`; P2b deferred ≥64 GB). Every spurion entry would show a representational raw residual, and the builder would have to invent the bridge. The entry representation and join key are also unspecified. | P2b `:11–13`; record `:56`; B has 0 `p2b` mentions | Raw comparison stays first. The P2b map is cited as the mapped-diagnostic convention, or the spurion comparison is declared deferred with P2b. The entry representation (jet-indexed, after `∂_2→ik_2`) is fixed for both engines. |
| R4-6 | The "author evidence" scripts are **typed** (literal `print(…, 0)`, typed parity dicts). Under FORM ablation (rotational bulk) every tag is byte-identical. A and the change log present them as computed evidence. Fold-bred. | author script `:31–32,44`; leg `07` | Remove the evidence claim, or relabel the scripts as typed bookkeeping. Correct the change log, which overstated R3-1 as "resolved". |
| R4-7 | A §6 is stale. PID 4097233 was **resumed** (commit `14861016`; feasibility `:47`; `/proc` state `R`). | lookups | A states that the job was resumed and that Part 1 waits on its completion (B §0.2 already stops on "exists in any state"). |
| R4-8 | P1 (no parity-odd constitutive term) is a **supplied input** of the basis ("in-plane `O(3)` isotropy and parity", spec `:116–117`). Part 1 cannot test it. | spec `:116–117` | A: P1 holds by construction of the supplied basis and is untestable by Part 1. The medium's achirality becomes an explicit part of R-LEAK-1's adopted condition, and F2's chiral term is its applicability test. |
| R4-9 | Two of the cited spin-carrier candidates are unclassified: **trapped chiral shear** (toroidal + poloidal mix has helicity ≠ 0 and mixes parities, so it is unprotected) and **mixed `a–w` circulation** (an in-plane vector; net spin breaks `O(3)`). N2 also omits an induced **mean flow** (a Coriolis image projects onto `ℓ = 0, 2`). | native_light `:114–116`; leg `09` | A: add an applicability test that the required trapped mode is pure twist-type (zero helicity); classify `a–w` circulation; N2 names the induced mean flow. |
| R4-10 | Minor points. B cites `AGENTS.md:8–12` as requiring guarded **Wolfram** runs, but those lines cover Python constructors and validators, and WL containment under the guard is untested. A §4 "would also trap gravity changes" is ungrounded. | AGENTS `:8–9`; A `:169` | B makes guarded execution its own requirement, cites AGENTS for the Python scope, and flags WL containment as untested. A drops or grounds the gravity-trap clause. |

**Author decision (G4).** Codex keeps authorship for v5. This is its first fold, and it bred R4-2, R4-4 and R4-6.
If v5 again breeds defects in material it changed, the author changes again (⛔ no repeated folds by a
defect-breeding author). v5 review: a fresh Claude agent + Grok, as new instances.

**Stop point surfaced to the user (they own the cut).** The physics question has strong review-leg evidence
(R4 Claude `02`). B is a blind dual-engine confirmation that likely cannot complete under the 2 GiB guard. The
leg measured 1.58 GB RSS just to load the exported slab payload.

---

# Round 5 — review of v5 (Codex-authored)

**Authorship and pairing:**
- v1–v3 are orchestrator-written; v4–v5 are Codex-revised (prompt `_legs/S11c_d_clean_condition_v5_author_prompt.md`).
- Legs: a **fresh Claude agent** (Opus, new instance) and **Grok** (`grok-4.6`, high).
- Prompt `_legs/S11c_d_clean_condition_review_round5_prompt.md`, identical for both. It corrects the stale
  "suspended" wording for PID 4097233 to "running".

**Reports:**
- `S11c_d_clean_condition_review_r5_claude.txt` (verbatim)
- `S11c_d_clean_condition_review_r5_grok.txt`

The scripts and stdout are in `S11c_d_clean_condition_review_scripts/r5_{claude,grok}/`.

**Verdicts:**
- **Claude:** "**not clear.** H-planar holds on the operator we actually have … Two findings change what may be
  claimed (F-1, F-2), and four defects change what B computes or how its output can be read (F-3 to F-6)."
- **Grok:** "**not cleared.** Hypothesis H is right for this model's operator on the class A names. Two defects in
  directive B change what Part 1 would compute."

**Adjudication: not cleared.** Every finding from both legs is verified below.

**Second independent measurement of H-planar on the real operator** (Claude S-a, `r5_claude/05_*`): `MIXED 0`
and `WRONG_SIDE 0` in all four cases (e.g. `{'EVEN': 237, 'ODD': 15, 'ZERO': 41}`), and the control gives
80/74/125/118. This agrees with R4's measurement. Both are single-engine review-leg measurements on the SymPy
export.

## Orchestrator verification — mechanical lookups, 2026-10-01

````
$ grep -n -E ':2416|:3103|:4156' research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md
114:at `:4156–4183`, from the whitelisted inputs in §4. The SymPy operator's code-defined top-level paths are exactly:
126:The first four are created at engine `:2416–2433`; the face and binding paths are attached at `:3103–3127`.
$ sed -n 3095,3096p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    operator["THETA_BALANCE"] = sp.Tuple(
        sp.Tuple(Str("SOURCE_OPERAND"), mass_balance),
$ sed -n 2967,2967p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    operator["U_BODY_BALANCE"] = sp.Tuple(
$ sed -n 290,291p research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md
**The slab momentum and thickness rows are the CONSTRAINT-REDUCED equations under S11b's binding material
virtual-displacement rule, not the held-fixed variational derivatives.** First compute the constitutive operand
$ sed -n 101,105p research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
x1 = inherited("x_1", "COORDINATE", "first Eulerian in-plane coordinate")
x2 = inherited("x_2", "COORDINATE", "second Eulerian in-plane coordinate")
x3 = inherited("x_3", "COORDINATE", "third Eulerian in-plane coordinate")
t = inherited("t", "COORDINATE", "time coordinate")
w = inherited("w", "COORDINATE", "bulk normal coordinate")
$ git show bcb9f7d7:research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | sed -n 2325p
            - (epsilon * rhobr * u_tt[a] if include_kinetic else 0)
$ sed -n 2375p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
            + (u_kinetic[a] if include_kinetic else 0)
$ sed -n 3034p research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py
    face_multiplier = named_tuple_row(face_normalization, "ACTION_TO_STORED_ROW_MULTIPLIER")
$ git log -1 --format='%h %ad %s' --date=short a74da30a ; git log -1 --format='%h %ad %s' --date=short c643112a
a74da30a 2026-09-10 S11c-b inertia repair + regenerated b/c1/c2 artifacts; S11c-d WIP checkpoint
c643112a 2026-09-12 S11c mechanical-load action normalization + regenerated d current and spectra; repair checkpoint
$ for k in U_BODY_BALANCE THETA_BALANCE E_W_BALANCE ADVECTIVE_MASS_OPERAND FACE_FLUX_BOUNDARY_OPERANDS projection_term_origins evolution_term_origins; do grep -c "$k" research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl; done
0 (each of the seven)
$ sed -n 112,120p docs/native_light_em_and_vortex_throat_interpretation.md | grep -c -i director
0
$ sed -n 530p research/pde_ledger_v3/V3_STEP_PLAN.md   (excerpt)
**Register:** **C13** is adjacent but ⛔ **not** this step's job (a gravitational wave is not a brane
$ sed -n 101p research/pde_ledger_v3/directives/S11b_SHARED_PHYSICS.md   (excerpt)
There is a **steady background transfer** of material across the interface. Let **`v_dr`** be the resulting
$ sed -n 100p docs/toy_model_ontology_summary.md   (excerpt)
… localized throat drainage transfers material, momentum, and energy from ordered brane degrees of freedom into de-structured bulk degrees of freedom, while distributed return transfers material back into the ordered state. …
````

The leg stdout relied on below:
- `r5_claude/03_drain_odd_channel.stdout.txt`:
  - `ODD_CHANNEL_ODE_AT_U0 … [0]`
  - `ODD_CHANNEL_SOLUTION_U Eq(a(w), C1*exp(I*omega*w/U))`
  - `FACE_TRACE_SOLUTION_m_nonzero [u3_t]`
  - `ODD_TO_ODD_FACE_BLOCK d(jump)/d(u3_t) -m`
  - `CONTROL_m0_ZERO_ROOTS [-I*sqrt(K)/sqrt(M), I*sqrt(K)/sqrt(M)]`
- `r5_claude/04_round_support_quadrupole.stdout.txt`:
  - `A_L2_PROJECTION_M0_OVER_TRACELESS_B33 -4*pi/5`
  - `CONTROL_ISOTROPIC_INCOHERENT_L2 {… all 0}`
  - `B_COUPLING_L2 T_x_plus_iT_y {… 1: -8*I*pi*delta/5 …}`
  - `CONTROL_ROUND_DELTA0 T_x_plus_iT_y {… all 0}`
  - `PARITY_MIRROR_Y T_z -1 Y_20 1 deformation_P2 1`
- `r5_grok/06_prefold_vs_emitted.stdout.txt`: `S11CA_76_100_HAS_X1 False`, `S11CA_101_HAS_X1 True`.

## Dispositions (all accepted)

| # | finding (leg) | verified by | what it changes |
|---|---|---|---|
| R5-1 | **The drain opens an odd channel** (Claude F-1). The parity selection rule survives a live normal drain, but "no leak" needs an **empty** odd channel. The rest-frame potential bulk has none (`ODD_CHANNEL_ODE_AT_U0 … [0]`). With mass flux across the face, the tangential momentum balance forces `v3_bulk = u3_t`, and the twist motion is advected away (`exp(iωw/U)`), an **allowed odd→odd** block. In the leg's lumped model the mode damps at `m0/(2M)`; the drain-off control is undamped. The drain is located at throats, where trapped light lives. ⇒ The no-leak **zero holds only with the drain frozen**, which is the M3 / "no static background" freeze. The leg's model is minimal (typed assumptions, computed roots), and the coefficient depends on the missing drain specification. | ontology `:100`; S11b `:101`; S11c-a `:368–369`; leg `03` | A's "exact zero" and "no linear leak" hold only for the **rest-frame (J⁰ = 0) operator**. With the drain live, H is a selection rule and R-LEAK-1 is a **bound**. A needs an empty-odd-channel premise and a falsifier for the drain's odd→odd channel. "Carries over" applies to the selection rule only. B Part 2 must surface the toroidal↔toroidal block including the bulk vortical sector. |
| R5-2 | **A spin-carrying trapped mode deforms its own throat** (Claude F-2). A coherent `ℓ=1` toroidal pattern's time-averaged energy always has an `ℓ=2` part, so a throat it holds open is not round. On that throat the `m=±1` (angular-momentum-carrying) pattern couples to `Y_21`; only `m=0` stays protected (mirror-odd vs even deformation). | leg `04` (controls: isotropic → 0; `δ=0` → 0) | A: state that the deformation is unavoidable for a coherent `ℓ=1` support mode; qualify F2's quadrature-pair sentence; tie F6/N2 to it. |
| R5-3 | B cites the pre-fold creation site (`:2416–2433`) and the face attachment (`:3103–3127`). The emitted rows are **overwritten** at `:2967–3101` by pin B, the face forces and the mass-evolution `θ` row (spec `:290–291`: constraint-reduced, not held-fixed). A builder following the citations computes a different operator for the even/TM blocks (Grok 1). Fold-bred. | engine `:2967`, `:3095–3096`; spec `:290–291` | B points at the emitted object, not at construction sites. |
| R5-4 | Typing spans: S11c-a `:76–100` does not contain `x_1…w` (at `:101–105`); `:568–577` lies inside `:552–599` (double-typed `ζ_c` family); following the spans literally holds `θ`, `e_W`, `δp_±` (Claude F-3, Grok 2). Fold-bred. | S11c-a `:101–105` | Disjoint, correct typing, or typing by registry description. |
| R5-5 | The sign-convention source (record `:112–114`, 09-03) predates the SymPy engine changes `a74da30a` (U-row kinetic sign) and `c643112a` (`face_multiplier`). These are among the **unreviewed upstream S11c-d repairs** (Claude F-4). | engine `:2375`, `:3034`; bcb9f7d7 `:2325`; git log | Each engine emits its own convention anchor; the record is not the map source. |
| R5-6 | The join key uses SymPy code paths, none of which exists in the WL engine. A "blind" WL would have to transcribe the SymPy construction (M1) (Claude F-5). Fold-bred. | WL grep: 0 for all seven | Join on physics objects, with a SymPy-emitted path→object table. |
| R5-7 | Citations: "director" is not in native_light `:112–120` (it is at `:2552`); C13 is uncited (V3 `:530`) (Claude F-6). | lookups | Corrected. |

**Open question raised by R5-1** (the orchestrator's question, ⛔ not a claim; a computation is owed). Ontology
`:100` also has **distributed return**, material flowing back into the brane away from throats. If that puts a
nonzero mass flux across the faces everywhere, does passing TE-like light meet the same odd→odd channel? R5-1's
leg model treats the face jump for flux `m` without a sign restriction, but nobody has computed this case.

## Process decision

- **Defect breeding.** v5 again bred defects in the material it changed (R5-3, R5-4, R5-6), so G4's author-change
  trigger fires a second time.
- **Recipe-creep.** This is the fifth round on B's specification of engine internals while the selection rule
  cleared every time (M2 recipe-creep tell).
- **The physics moved.** R5-1 shows that what B's Part 1 computes, the drain-frozen rest-frame operator, is exactly
  the freeze that manufactures the zero. Two review legs have already measured the selection rule on the real
  export (R4, R5).

⇒ **The B iteration is paused, and the scope goes to the user** (the user owns the cut; this is a scope change).
No further fold until the user chooses. v5 is preserved as reviewed.
