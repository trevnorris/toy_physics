# _measurements — S11b SymPy fix-round-2 directive (rule 2)

Every claim carries the command that produced it. SCRIPT-GENERATED
(`scratchpad/gen_s11b_fix_round2_twin.sh`), never transcribed. Engine at round-1-repair baseline `6d57b27e`.
The behavioral probes (kernel global-transpose -> residuals [0,0,0]; homogeneity term-corruption True->False)
live in the Codex consult log `~/.s11_build/s11b_round2_disposition_consult.log` and the two round-1-repair
script-leg reports; the commands below substantiate each row's structure, its export status, and the onsager
arithmetic.

### Item 1 — parity_even_jw AND parity_odd_jw both reduce to the outer integral_coefficient x typed placeholder algebra (667-673)

```
$ sed -n '667,673p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    upper_coefficient = source_finite.coeff(upper_product)
    lower_coefficient = source_finite.coeff(lower_product)
    even_pair = sp.expand(integral_coefficient * (Omp * je + (-Omp) * je))
    odd_pair = sp.expand(integral_coefficient * (Omp * jo + (-Omp) * (-jo)))
    even_boundary = sp.expand(upper_coefficient * Om * je + lower_coefficient * Om * je)
    odd_boundary = sp.expand(upper_coefficient * Om * jo + lower_coefficient * Om * (-jo))
    emit("PARITY_EVEN_JW", sp.Tuple(even_pair, even_boundary, sp.simplify(even_pair + even_boundary)), key="parity_even_jw")
```

### Item 2 — control_no_reciprocal_traction: BOTH routes drop Lambda_X; export=False (1919-1942)

```
$ sed -n '1919,1942p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    control_f_model = derive_model(substitutions={Lambda_X0: 0})
    operator_substituted = sp.cancel(MODEL["thickness"].subs(Lambda_X0, 0))
    operator_residual = sp.cancel(control_f_model["thickness"] - operator_substituted)
    V = symbol("V_control_F", "CONTROL", "control-F face velocity amplitude")
    A = symbol("A_control_F", "CONTROL", "control-F affinity amplitude")
    p = symbol("p_control_F", "CONTROL", "control-F face pressure amplitude")
    J = symbol("J_control_F", "CONTROL", "control-F face flux amplitude")
    mu_s = A + p / rho_m
    v_bulk = V + J / rho_m
    full_power = sp.re((p + Lambda_X * A) * sp.conjugate(V) + mu_s * sp.conjugate(J)) / 2
    recomputed_power = sp.re(p * sp.conjugate(V) + mu_s * sp.conjugate(J)) / 2
    substituted_power = sp.simplify(full_power.subs(Lambda_X0, 0))
    power_residual = sp.simplify(recomputed_power - substituted_power)
    control_f_compression = compressional_response(control_f_model)
    emit(
        "CONTROL_NO_RECIPROCAL_TRACTION",
        sp.Tuple(
            control_f_model["inplane"], control_f_model["thickness"],
            control_f_compression["ratio"], control_f_model["determinant"],
            operator_substituted, operator_residual,
            recomputed_power, substituted_power, power_residual,
        ),
        key=None,
        export=False,
```

### Item 3 — onsager crosscheck = flux_res - force_res, with force_res = -flux_res by construction (1370-1385)

```
$ sed -n '1370,1385p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    epsilon = symbol("epsilon_port", "CONTROL", "formal mixed-law conversion coefficient", nonzero=True, real=True)
    a_entry, b_entry, c_entry, d_entry = Lambda_A, Lambda_V, Lambda_X, epsilon
    all_flux = sp.ImmutableMatrix([
        [a_entry - b_entry * c_entry / d_entry, b_entry / d_entry],
        [-c_entry / d_entry, 1 / d_entry],
    ])
    flux_reciprocity_residual = sp.factor(sp.together(all_flux[0, 1] - all_flux[1, 0]) * epsilon)
    all_force = sp.ImmutableMatrix([
        [1 / a_entry, -b_entry / a_entry],
        [c_entry / a_entry, d_entry - c_entry * b_entry / a_entry],
    ])
    force_reciprocity_residual = sp.factor(sp.together(all_force[0, 1] - all_force[1, 0]) * a_entry)
    reciprocity_equation = sp.Eq(flux_reciprocity_residual, 0, evaluate=False)
    reciprocity_crosscheck = sp.simplify(flux_reciprocity_residual - force_reciprocity_residual)
    emit("ONSAGER_CONDITION", sp.Tuple(all_flux, reciprocity_equation), key="onsager_condition")
    emit("ONSAGER_RECIPROCITY", sp.Tuple(all_force, flux_reciprocity_residual, force_reciprocity_residual, reciprocity_crosscheck), key="onsager_reciprocity")
```

### Item 3 — LIVE arithmetic probe on the exported onsager_reciprocity value: force = -flux, cross = 2*flux

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3/scripts && python3 -c "import S11b_exports as m, sympy as sp; v=m.LEDGER['onsager_reciprocity']['value']; print('flux+force =', sp.simplify(v[1]+v[2])); print('cross - 2*flux =', sp.simplify(v[3]-2*v[1]))"
flux+force = 0
cross - 2*flux = 0
```

### Item 4 — kernel_orientation uses the module-global Lambda_A/V/X on both sides (942-957)

```
$ sed -n '942,957p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    closure_law = MODEL["face"]["closure_template"]
    live_load_plus = live_face_loads(MODEL)[0]
    extracted = (
        sp.diff(closure_law.subs(velocity_coordinate, 0), affinity_coordinate),
        sp.diff(closure_law.subs(affinity_coordinate, 0), velocity_coordinate),
        sp.diff(live_load_plus - delta_p_plus, A_plus_affinity),
    )
    supplied = (Lambda_A, Lambda_V, Lambda_X)
    orientation_rows = []
    for name, left, right, coefficient, relaxation in zip(
        ("A", "V", "X"), extracted, supplied, (Lambda_A0, Lambda_V0, Lambda_X0), (tau_A, tau_V, tau_X)
    ):
        residual = sp.cancel(left - right)
        status, test = zero_test(sp.together(residual).as_numer_denom()[0])
        indistinguishable = sp.Or(sp.Eq(coefficient, 0, evaluate=False), sp.Eq(relaxation, 0, evaluate=False))
        orientation_rows.append(sp.Tuple(Str(name), left, right, residual, status, test, indistinguishable))
```

### Item 4 — causality_check carries the same kernel-orientation payload (1123-1127)

```
$ sed -n '1123,1127p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    causality = causality_objects()
    emit("KERNEL_ORIENTATION_IDENTITIES", causality["orientation"], key="kernel_orientation_identities")
    emit("KERNEL_PROPAGATION_RESIDUALS", causality["propagation"], key="kernel_propagation_residuals")
    emit("KERNEL_POLE_LOCATIONS", causality["poles"], key="kernel_pole_locations")
    emit("CAUSALITY_CHECK", sp.Tuple(causality["orientation"], causality["propagation"], causality["poles"]), key="causality_check")
```

### Leave — homogeneity: trace_dimension RECURSES into every nested Add (numerator/denominator) (1598-1616)

```
$ sed -n '1598,1616p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
def trace_dimension(expression: sp.Expr, symbol_dimensions: Mapping[sp.Symbol, sp.MatrixBase]) -> DimensionTrace:
    """Derive an expression's dimension and every internal additive test."""
    expression = sp.sympify(expression)
    if expression.is_Number or expression in (I, sp.oo, -sp.oo, sp.nan):
        return DimensionTrace(ZERO_DIM, ())
    if isinstance(expression, sp.Symbol):
        if expression not in symbol_dimensions:
            raise KeyError(f"no live dimension route for {expression}")
        return DimensionTrace(sp.ImmutableMatrix(symbol_dimensions[expression]), ())
    if isinstance(expression, sp.Add):
        traces = tuple(trace_dimension(argument, symbol_dimensions) for argument in expression.args)
        reference = traces[0].vector
        tests = tuple(test for trace in traces for test in trace.tests)
        tests += tuple(dimension_match(trace.vector, reference) for trace in traces[1:])
        return DimensionTrace(reference, tests)
    if isinstance(expression, sp.Mul):
        traces = tuple(trace_dimension(argument, symbol_dimensions) for argument in expression.args)
        vector = sum((trace.vector for trace in traces), ZERO_DIM)
        tests = tuple(test for trace in traces for test in trace.tests)
```

### Export status of every round-2 row (parity+onsager+kernel+causality exported; control NOT)

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3/scripts && python3 -c "import S11b_exports as m; L=m.LEDGER; [print(f'{k}: exported={k in L}') for k in ['parity_even_jw','parity_odd_jw','onsager_reciprocity','kernel_orientation_identities','causality_check','homogeneity_thickness_eom','control_no_reciprocal_traction']]"
parity_even_jw: exported=True
parity_odd_jw: exported=True
onsager_reciprocity: exported=True
kernel_orientation_identities: exported=True
causality_check: exported=True
homogeneity_thickness_eom: exported=True
control_no_reciprocal_traction: exported=False
```
