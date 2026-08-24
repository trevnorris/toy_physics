# Measurements twin — S11c_a_sympy_engine_repair_directive.md

Every claim in the directive carries the command that produced it (CLAUDE.md rule 2).

## Leg findings (the two quoted legs)
- Grok leg log: `~/.s11_build/S11c_a_sympy_engine_grok.log` (bg task bui56192f, exit 0, 10583 bytes).
  Its own derivation/ablation scripts: `/tmp/s11ca_review/{derive_s11ca_independent,ablations/*}.py` + `/measurements/*.stdout.txt`.
- Agent leg (af1c189cb7f8e6890) result: recorded in session; scripts `/tmp/s11ca_{agent_deriv,route_identity,onesided}.py` + `.out`,
  ablations `/tmp/s11ca_ablate_{normal,grad}.out`, baseline `/tmp/s11ca_baseline.out`.
- Both quoted verbatim in the directive. Both converge; orchestrator code-verified below.

## Defect 1 code facts
CMD: sed -n '604,616p' research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py
  ⇒ EULERIAN: grad_h = dx(h_exact,i,scales); MATERIAL else: tangent_vertical = dx(h_exact,i,scales) [SAME jets];
     cofactor=(-face*component,...,face); norm=sqrt(dot(cofactor,cofactor)). face^2=1 ⇒ norm=sqrt(1+|grad_h|^2)=denominator
     ⇒ MATERIAL normal ≡ EULERIAN normal (algebraically identical).
CMD: sed -n '970,984p' .../S11c_a_interface_geometry_sympy_audit.py
  ⇒ virtual_constraint_route MATERIAL: flattened_mass=density*thickness*jacobian; return eps*shape(flattened_mass/background_mass).
     EULERIAN else: material_mass=density*thickness*jacobian; return eps*shape(material_mass/background_mass). IDENTICAL product.
CMD: sed -n '490,497p' .../S11c_a_interface_geometry_sympy_audit.py
  ⇒ source_jet_scales(...,reverse_upper_x1,face): if reverse_upper_x1 and face==1: scales[0]=-1. Applied in build_face_source
     (line 597) BEFORE the route split ⇒ feeds dx() for BOTH routes ⇒ shared reverse moves both operands.
CMD: sed -n '1202,1232p' .../S11c_a_interface_geometry_sympy_audit.py
  ⇒ task_rep_invariance: RELATIVE_FLUX/TRACTION/CLOSURE_SHAPE_DERIV/VIRTUAL_WORK use route="MATERIAL" (aliased) ⇒ residual 0 by construction.
     VIRTUAL_CONSTRAINT material=build_virtual_constraint_raw(route="MATERIAL") [aliased]; EVOLUTION material=build_evolution_raw(route="MATERIAL") [GENUINE].
  ⇒ task_independence: geometry quantities corrupted via reverse_upper_x1 (shared); T-g/T-h via corrupt_direct.

## T-h genuine (leave alone)
CMD: sed -n '1002,1042p' .../S11c_a_interface_geometry_sympy_audit.py
  ⇒ evolution_route EULERIAN: density_time + dilatation(sigma0*trace_grad(grad_u_t)) + advection(u_t·dx(sigma0)).
     MATERIAL: material_mass=(sigma_exact+evaluation_shift)*det_jacobian(grad_u); density_time=dt(shape(...)); dilatation=advection=0.
     STRUCTURALLY DISTINCT ⇒ genuine second route.

## Defect 2 code facts
CMD: sed -n '698,705p' .../S11c_a_interface_geometry_sympy_audit.py
  ⇒ kinematic_raw: lhs=eps*shape(n·v_bulk); residual=factor_terms(lhs - face_velocity_raw - relative_flux_raw/rho_m); return Eq(residual,0).
     relative_flux_raw/rho_m = eps*shape((v_bulk-v_face)·n); face_velocity_raw=eps*shape(n·v_face) ⇒ residual ≡ 0 identically. Bare Eq(0,0).

## Spec authority (route-2 = w', unchanged, already 2-legged in the closed spec)
CMD: sed -n '452,481p' research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
  ⇒ §5a route-2 = material-coordinate diff via w'=[w-ζ_c]/[W_bg+δW], mapped back to Eulerian; one-sided control mutates ONLY the direct route.
