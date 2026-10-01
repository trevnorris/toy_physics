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
