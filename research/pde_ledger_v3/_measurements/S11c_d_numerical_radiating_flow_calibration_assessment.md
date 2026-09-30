# Flow, calibration and anchoring: source-only decision checkpoint

The current pilot answers a strict rest-bulk, LAB_HELD/RHO4_CONSTANT development-input question. It does not yet answer leakage in the calibrated, draining medium. These limitations were explicit in the governing scope, but should have been foregrounded in the pilot decision. No result can remove them through numerical convergence alone.

On the user's request to assess before continuing and stop for go/no-go, the existing central-balance-v2 scientific child was suspended in memory, not terminated or restarted. The pause receipt records PID 4097233, its process start identity, cgroup and exact command. It had completed two matrix groups after twenty restorations and saved 12,288 batches of the next group; no finite solve had been reached. The guard/supervisor remain active, no deadline was added, and pinned sources are unchanged. Resume requires the user's decision and rechecking that process identity. Do not launch another job alongside it.

This assessment reads source, JSON and hashes only. No scientific payload was restored. The accompanying assessment JSON verifies the three consumed modal JSON files against the accepted end-map checkpoint and records the arithmetic and source hashes.

## Speed ratio

The development input sets mu_R=1, rho_br=1 and c_s0=10 in the L_ref/T_ref unit frame. Under the bare S9/R4 definition c_gamma^2=mu_R/rho_br, this gives c_gamma/c_s=0.1, not 1. The step plan explicitly distinguishes the derived ratio definition from the calibrated, uncommitted equality lambda_gamma=1; the equality has not been imposed on these inputs.

There is an additional distinction: that bare coefficient formula must not be substituted for the actual full S11c-d transverse dispersion. The saved omega=3 LEFT incoming branch has normal momentum 2.439262183530097 and tangential norm squared 1/20. Its fixed-point phase speed omega/sqrt(k_normal^2+1/20) is 1.2247448713915874, or 0.12247448713915873 of c_s0. At the full-contrast RIGHT end the corresponding ratio is 0.12186666955535794; its zero-contrast ratio returns to the LEFT value. These are arithmetic readouts of saved full-pencil modes, not a new root solve or a demonstrated identification with the calibrated light cone. The extra retained elastic structure and coefficient conventions cannot be silently discarded.

Thus neither the bare input ratio nor the actual selected modal ratio realizes equal light/bulk speeds. Identifying this development slice with the calibrated analog-light band remains OPEN. Changing that calibration would change physical inputs and the branch/phase-matching problem, not just relabel a plot.

Sources: `steps/S9_light_requires_shear.md:78`, `V3_STEP_PLAN.md:1107–1111`, `S11c_d_variable_profile_development_input.json`, accepted `S11c_d_numerical_radiating_end_maps_checkpoint.json` and the source paths/hashes in this assessment JSON.

## Omitted drain flow

The solved bulk operator is the rest-frame wave equation. The convective operator with (partial_t+v_bulk_normal_0 partial_w)^2 is not constructed. This is an explicit inherited S11c scope decision, not evidence that the physical drain vanishes. The driven-background S11bB discussion also requires a named reservoir and stated power budget before interpreting drive-fed channels; it supplies no numerical drain/response bound usable for this pilot.

The required result label is: **strict v_bulk_normal_0=0, LAB_HELD/RHO4_CONSTANT, uncalibrated development speeds, physical face permeability/memory retained**. A stationary held profile and omitted bulk flow are separate assumptions; the present calculation still contains time-dependent fluctuations and physical face dissipation.

At omega=3 and the saved tangents, the physical depth momentum satisfies

    q^2 = omega^2/c_s0^2 - |k_parallel|^2 - k_normal^2
        = 1/25 - k_normal^2.

The rest-bulk propagating interval is |k_normal|<1/5 and grazing occurs at its endpoints. Away from grazing the governing specification requires BOTH

    |q v|/3 << 1,       3 |v|/(100 |q|) << 1,

as well as the independent subsonic condition |v|<10. Equivalently, |v| must be much smaller than both 3/|q| and 100|q|/3, in reference speed units. At k_normal=0, q=1/5: the stricter smallness scale is 20/3. This is a scale, not a certified permissible bound. The first condition also becomes more restrictive at large evanescent |q|.

Crucially, the second condition fails arbitrarily close to q=0 for every fixed nonzero v. The integration crosses these grazing endpoints, so the first inequality alone cannot justify applying the whole rest-bulk result to a flowing medium. The amendment already says exactly this. No numerical value or applicable smallness bound for the model's drain has been established in the supplied pilot inputs. Moving-medium validity is therefore unresolved, not known to hold. A nonzero drain can change dispersion and energy exchange; its quantitative effect is not calculated here.

Sources: `directives/S11c_decisions.md:180–186`, `directives/S11c_c1_SHARED_PHYSICS.md:204–220`, `directives/S11c_d_SCATTERING_FORM_AMENDMENT.md:373–383`, `steps/S11bB_interface_assembly.md:45–65`.

## MATERIAL_ADVECTED feasibility

The two anchorings are physical alternatives: LAB_HELD uses Q_bg(x), while MATERIAL_ADVECTED uses Q_bg(chi(x,t)), with chi the inverse brane material map. Both currently have a supplied support-stabilized background. Selecting MATERIAL_ADVECTED does not introduce the missing bulk-normal drain or establish a driven background power budget.

There is a reusable source route: accepted remaining-case live-frequency records cover all four cases, and the end-input catalogue covers all eight physical ends. Those catalogue acceptances do not supply validated omega=3 material-advected candidate paths, current/face maps or finite-response controls. Existing material numerical rows at 1-0.01i are not omega=3 results.

A material-advected comparison is feasible as a separate case using the same numerical machinery, with its own source bindings, end continuation/sheet/current/face checks, affected matrices and zero-contrast/scaling/domain/regulator/refinement controls. Reuse is permitted only where actual operands match. It is not a configuration switch on the current LAB_HELD matrices, and it would consume additional solve budget (a second full three-frequency/ten-control suite would exceed the present 30-solve scope). Jobs would run sequentially under containment. No runtime estimate for that extension is established by this source inspection.

Sources: `directives/S11c_a_SHARED_PHYSICS.md:233–268`, `S11c_d_remaining_case_frequency_sources_report.md`, `S11c_d_remaining_case_frequency_end_inputs_report.md`, and the retained case-specific rational/row reports.

## Decision boundary

Do not add MATERIAL_ADVECTED merely to address the flow concern: it does not address it. First decide whether to finish this explicitly labelled rest-bulk development benchmark, or redirect the physics scope toward the calibrated, flowing background. The latter requires a specified background flow and calibration, and an applicable boundary/operator/power treatment; its cost is not established here and it is not authorized by this assessment. Keep the current calculation suspended and preserve all saved work until the user's go/no-go.

The independent review findings and numerical interpretation checks remain pending as recorded in `S11c_d_numerical_radiating_replacement_review_disposition.md`. This scope assessment neither clears those findings nor accepts a deficit, physical loss, analog-light calibration, full Green/FORM or A11/A12.
