# _measurements — S11b SymPy fix-round-1 directive (rule 2)

Every claim in `S11b_sympy_engine_fix_round1_directive.md` carries the command that produced it.
SCRIPT-GENERATED (`scratchpad/gen_s11b_fix_round1_twin.sh`), never transcribed. Engine at round-1 baseline
`7dd89076`. The one-sided-corruption / FORM-ablation evidence lives in the two round-1 script-leg reports and
the two directive-review-leg reports (`~/.s11_build/s11b_*_leg*.txt` / `.log`); the commands below
substantiate every worklist row's presence, its typed/tautological element, and the export-row membership.

### Fix 1 — general impedance TYPED at module level (138); bulk solve computes the same, emits Z_IMPERMEABLE (612-620)

```
$ sed -n '138p;612,620p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
Z_GENERAL = sp.cancel(rho_m * omega / q_out)
    phi_plus = A_plus * sp.exp(I * q_out * (w - W0 / 2))
    phi_minus = A_minus * sp.exp(I * q_out * (-w - W0 / 2))
    v_plus = sp.diff(phi_plus, w).subs(w, W0 / 2)
    v_minus = -sp.diff(phi_minus, w).subs(w, -W0 / 2)
    p_plus = I * rho_m * omega * phi_plus.subs(w, W0 / 2)
    p_minus = I * rho_m * omega * phi_minus.subs(w, -W0 / 2)
    z_plus = sp.cancel(p_plus / v_plus)
    z_minus = sp.cancel(p_minus / v_minus)
    emit("Z_IMPERMEABLE", sp.Tuple(z_plus, z_minus), key="z_impermeable")
```

### Fix 1 — assembly DEFAULTS to typed Z_GENERAL (478); z_plus/z_minus consumed nowhere else

```
$ grep -n 'Z_GENERAL\|z_value: sp.Expr = \|z_plus\|z_minus\|Z_IMPERMEABLE' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
138:Z_GENERAL = sp.cancel(rho_m * omega / q_out)
478:    z_value: sp.Expr = Z_GENERAL,
618:    z_plus = sp.cancel(p_plus / v_plus)
619:    z_minus = sp.cancel(p_minus / v_minus)
620:    emit("Z_IMPERMEABLE", sp.Tuple(z_plus, z_minus), key="z_impermeable")
628:    z_prop = real_axis_simplify(Z_GENERAL.subs(q_out, q_prop))
629:    z_evan = real_axis_simplify(Z_GENERAL.subs(q_out, q_evan))
630:    grazing_sides = sp.Tuple(sp.limit(Z_GENERAL, q_out, 0, dir="+"), sp.limit(Z_GENERAL, q_out, 0, dir="-"))
644:        sp.Tuple(Str("THICKNESS"), thickness_velocities, sp.simplify(Z_GENERAL * thickness_velocities)),
645:        sp.Tuple(Str("CENTRE_SHIFT"), centre_velocities, sp.simplify(Z_GENERAL * centre_velocities)),
698:    face = face_solution(Z_GENERAL, mu_drive, V_face)
1111:    sliced_face = face_solution(Z_GENERAL, sp.Integer(0), V_face)
1150:    response = face_solution(Z_GENERAL, rho_br * mu_s, V)
1165:    z_dependence = sp.Tuple(*sorted(H_port.free_symbols.intersection(Z_GENERAL.free_symbols), key=sp.default_sort_key))
1506:    face = face_solution(Z_GENERAL, mu_theta_value, sp.Integer(0))
1593:    control_e_power_response = face_solution(Z_GENERAL, sp.Integer(0), V_face)
```

### Fix 1 — regime impedances RE-TYPED (not via the solve): 714-720 and 990-996

```
$ sed -n '714,720p;990,996p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    q_squared = omega**2 / c_s0**2 - k**2
    z_prop = rho_m * sp.Abs(omega) / sp.sqrt(q_squared)
    kappa = next(sym for sym in DECLARED_SYMBOLS if sym.name == "kappa_out")
    z_evan = -I * rho_m * omega / kappa
    coefficient_by_regime = []
    for regime_name, z_value in (("PROPAGATING", z_prop), ("EVANESCENT", z_evan), ("GRAZING", sp.oo)):
        response = face_solution(z_value, sp.Integer(0), V_face)["pressure"]
    q_squared = omega**2 / c_s0**2 - k**2
    kappa = next(sym for sym in DECLARED_SYMBOLS if sym.name == "kappa_out")
    regime_z = (
        (Str("PROPAGATING"), rho_m * sp.Abs(omega) / sp.sqrt(q_squared)),
        (Str("EVANESCENT"), -I * rho_m * omega / kappa),
        (Str("GRAZING"), symbol("Z_grazing_control", "CONTROL", "grazing face-response limit coordinate", positive=True)),
    )
```

### Fix 1 — breathing k=0 slice impedance RE-TYPED z_value=rho_m*c_s0 (1336-1339)

```
$ sed -n '1336,1339p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    slice_model = derive_model(
        substitutions={Lambda_A0: 0, Lambda_V0: 0, Lambda_X0: 0, k: 0},
        z_value=rho_m * c_s0,
    )
```

### Fix 2 — every worklist self-check row is EXPORTED as a ledger row

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3/scripts && python3 -c "import S11b_exports as m; L=m.LEDGER; rows=['convention_check_inplane','convention_check_conservative','pressure_work_sign_check','full_two_port_balance_check','energy_sinks','energy_sources','unattributed_sink_terms','unattributed_exchange_terms','kernel_orientation_identities','kernel_propagation_residuals','conservative_positivity_inequality','two_port_power_identity','parity_even_jw','parity_odd_jw','homogeneity_inplane_eom','homogeneity_ablation_demo','dim_route_kind_z','parity_interval','sheet_of_each_root','branch_degenerate_point']; [print(f'{k}: exported={k in L}') for k in rows]"
convention_check_inplane: exported=True
convention_check_conservative: exported=True
pressure_work_sign_check: exported=True
full_two_port_balance_check: exported=True
energy_sinks: exported=True
energy_sources: exported=True
unattributed_sink_terms: exported=True
unattributed_exchange_terms: exported=True
kernel_orientation_identities: exported=True
kernel_propagation_residuals: exported=True
conservative_positivity_inequality: exported=True
two_port_power_identity: exported=True
parity_even_jw: exported=True
parity_odd_jw: exported=True
homogeneity_inplane_eom: exported=True
homogeneity_ablation_demo: exported=True
dim_route_kind_z: exported=True
parity_interval: exported=True
sheet_of_each_root: exported=True
branch_degenerate_point: exported=True
```

### convention_check_inplane — fresh placeholders, not MODEL['inplane'] (937-941); MODEL['inplane'] IS emitted (923)

```
$ sed -n '923p;937,941p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    emit("INPLANE_EOM", sp.Eq(MODEL["inplane"], 0, evaluate=False), key="inplane_eom")
    sigma_coordinate = symbol("sigma_eta_coordinate", "COORDINATE", "independent explicit strain stress")
    mu_coordinate = symbol("mu_theta_coordinate", "COORDINATE", "independent chemical potential density")
    constrained_stress = sigma_coordinate - mu_coordinate
    inplane_coefficient = sp.diff(constrained_stress, mu_coordinate)
    emit("CONVENTION_CHECK_INPLANE", sp.Tuple(constrained_stress, inplane_coefficient, sp.Eq(inplane_coefficient, -1, evaluate=False)), key="convention_check_inplane")
```

### convention_check_conservative — both operands from U_LONG, never MODEL['thickness']; omega_squared ornament (943-953)

```
$ sed -n '943,953p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    reduced_energy = sp.expand(U_LONG.subs({k: 0, kappa_W: 0, theta: -e_W, eta: 0}))
    K_check = sp.diff(reduced_energy, e_W, 2)
    p_w = sp.diff(U_LONG, e_W)
    mu_theta_value = sp.diff(U_LONG, theta)
    equation_stiffness = sp.diff((p_w - mu_theta_value).subs({k: 0, kappa_W: 0, theta: -e_W, eta: 0}), e_W)
    stiffness_residual = sp.simplify(equation_stiffness - K_check)
    omega_squared = sp.solve(sp.Eq(-mu_W * W0**2 * omega**2 + equation_stiffness, 0), omega**2, dict=True)[0][omega**2]
    conservative_inequality = sp.Gt(omega_squared, 0, evaluate=False)
    positive_energy = sp.And(sp.Gt(B_rho_3, 0, evaluate=False), sp.Gt(B_rho_3 * k_W * W0**2 - (C * W0)**2, 0, evaluate=False))
    positivity_implication = sp.Implies(positive_energy, conservative_inequality, evaluate=False)
    emit("CONVENTION_CHECK_CONSERVATIVE", sp.Tuple(equation_stiffness, K_check, stiffness_residual, omega_squared, sp.diff(omega_squared, B_rho_3)), key="convention_check_conservative")
```

### MODEL['thickness'] is the live assembled equation the check should police (503)

```
$ sed -n '503p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    thickness = sp.cancel(
```

### pressure_work_sign_check / full_two_port_balance_check / energy sinks-sources / unattributed (853-968 head+tail)

```
$ sed -n '853,858p;885,889p;912,913p;962,968p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    Vm = symbol("V_minus", "COORDINATE", "lower outward face velocity amplitude")
    pp = symbol("delta_p_plus", "COORDINATE", "upper face pressure amplitude")
    pm = symbol("delta_p_minus", "COORDINATE", "lower face pressure amplitude")
    Ap = symbol("A_plus_affinity", "COORDINATE", "upper affinity amplitude")
    Am = symbol("A_minus_affinity", "COORDINATE", "lower affinity amplitude")
    Jp = next(sym for sym in DECLARED_SYMBOLS if sym.name == "J_plus")
    Zprop = symbol("Z_prop_real", "COORDINATE", "real outgoing propagating face response", positive=True)
    Vabs2 = symbol("V_abs_squared", "COORDINATE", "squared face-velocity amplitude", nonnegative=True)
    slab_pressure_route = sp.simplify(-2 * sp.re(Zprop) * Vabs2 / 2)
    bulk_flux_route = sp.simplify(2 * sp.re(Zprop) * Vabs2 / 2)
    pressure_residual = sp.simplify(slab_pressure_route + bulk_flux_route)
        "pressure_check": sp.Tuple(REAL_AXIS, slab_pressure_route, -bulk_flux_route, pressure_residual),
        "two_port_check": sp.Tuple(REAL_AXIS, *order_rows),
    ))
    emit("ENERGY_SINKS", sinks, key="energy_sinks")
    emit("ENERGY_SOURCES", sources, key="energy_sources")
    emit("UNATTRIBUTED_SINK_TERMS", sp.Tuple(), key="unattributed_sink_terms")
    emit("UNATTRIBUTED_EXCHANGE_TERMS", sp.Tuple(), key="unattributed_exchange_terms")
    emit("PRESSURE_WORK_SIGN_CHECK", accounting["pressure_check"], key="pressure_work_sign_check")
    emit("FULL_TWO_PORT_BALANCE_CHECK", accounting["two_port_check"], key="full_two_port_balance_check")
```

### kernel_orientation_identities — d(Lambda*coord)/dcoord - Lambda on a typed copy (808-826)

```
$ sed -n '808,826p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
def causality_objects() -> dict[str, object]:
    affinity_coordinate = symbol("affinity_coordinate", "COORDINATE", "independent affinity amplitude")
    velocity_coordinate = symbol("velocity_coordinate", "COORDINATE", "independent outward velocity amplitude")
    closure_law = Lambda_A * affinity_coordinate + Lambda_V * velocity_coordinate
    traction_law = p_face + Lambda_X * affinity_coordinate
    extracted = (
        sp.diff(closure_law.subs(velocity_coordinate, 0), affinity_coordinate),
        sp.diff(closure_law.subs(affinity_coordinate, 0), velocity_coordinate),
        sp.diff(traction_law - p_face, affinity_coordinate),
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

### kernel_propagation_residuals — BOTH routes call derive_model() (828-843), and MODEL=derive_model() (2035)

```
$ sed -n '828,843p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py; echo '---'; sed -n '2035p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    placeholder_model = derive_model(lambda_a=ell_A, lambda_v=ell_V, lambda_x=ell_X)
    replacements = {ell_A: Lambda_A, ell_V: Lambda_V, ell_X: Lambda_X}
    propagation_rows = []
    for name in ("face", "mass", "inplane", "thickness", "determinant"):
        if name == "face":
            placeholder_object = placeholder_model["face"]["pressure"]
            actual_object = MODEL["face"]["pressure"]
        else:
            placeholder_object = placeholder_model[name]
            actual_object = MODEL[name]
        propagated = casify(placeholder_object).subs(replacements)
        if isinstance(actual_object, sp.MatrixBase):
            residual = sp.ImmutableMatrix(propagated - actual_object)
        else:
            residual = sp.Add(propagated, -actual_object, evaluate=False)
        propagation_rows.append(sp.Tuple(Str(name.upper()), propagated, actual_object, residual))
---
    MODEL = derive_model()
```

### conservative_positivity_inequality — positive_energy is a TYPED And(Gt,Gt) (950-955)

```
$ sed -n '950,955p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    conservative_inequality = sp.Gt(omega_squared, 0, evaluate=False)
    positive_energy = sp.And(sp.Gt(B_rho_3, 0, evaluate=False), sp.Gt(B_rho_3 * k_W * W0**2 - (C * W0)**2, 0, evaluate=False))
    positivity_implication = sp.Implies(positive_energy, conservative_inequality, evaluate=False)
    emit("CONVENTION_CHECK_CONSERVATIVE", sp.Tuple(equation_stiffness, K_check, stiffness_residual, omega_squared, sp.diff(omega_squared, B_rho_3)), key="convention_check_conservative")
    emit("CONSERVATIVE_POSITIVITY_INEQUALITY", sp.Tuple(conservative_inequality, positive_energy, positivity_implication), key="conservative_positivity_inequality")

```

### two_port_power_identity — both sides from the same fresh V,mu_s,p,J ⇒ residual 0 (1141-1148)

```
$ sed -n '1141,1148p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    A = sp.expand(mu_s - p / rho_m)
    v_bulk = sp.expand(V + J / rho_m)
    left_power = sp.expand((p + Lambda_X * A) * sp.conjugate(V) + mu_s * sp.conjugate(J))
    right_power = sp.expand(p * sp.conjugate(v_bulk) + A * sp.conjugate(J) + Lambda_X * A * sp.conjugate(V))
    left_power_real = real_axis_simplify(sp.re(left_power) / 2)
    right_power_real = real_axis_simplify(sp.re(right_power) / 2)
    power_residual = sp.simplify(left_power_real - right_power_real)
    emit("TWO_PORT_POWER_IDENTITY", sp.Tuple(REAL_AXIS, left_power_real, right_power_real, power_residual), key="two_port_power_identity")
```

### parity_even_jw / parity_odd_jw — hardcoded algebra on placeholders; live source at 576-577 (576,581-590)

```
$ sed -n '576,577p;581,590p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    source_finite = sp.Integral(sp.diff(Omega, w) * jw, (w, w1, w2)) - (Omega * jw).subs(w, w2) + (Omega * jw).subs(w, w1)
    source_infinite = sp.Integral(sp.diff(Omega, w) * jw, (w, -sp.oo, sp.oo))
    Om = symbol("Omega_w", "COORDINATE", "window value at a reflected point")
    Omp = symbol("Omega_prime_w", "COORDINATE", "window derivative at a reflected point")
    je = symbol("j_even_w", "COORDINATE", "even normal-current value")
    jo = symbol("j_odd_w", "COORDINATE", "odd normal-current value")
    even_pair = sp.expand(Omp * je + (-Omp) * je)
    odd_pair = sp.expand(Omp * jo + (-Omp) * (-jo))
    even_boundary = sp.expand(-Om * je + Om * je)
    odd_boundary = sp.expand(-Om * jo + Om * (-jo))
    emit("PARITY_EVEN_JW", sp.Tuple(even_pair, even_boundary, sp.simplify(even_pair + even_boundary)), key="parity_even_jw")
    emit("PARITY_ODD_JW", sp.Tuple(odd_pair, odd_boundary, sp.simplify(odd_pair + odd_boundary)), key="parity_odd_jw")
```

### homogeneity family stamped + ablation demo corrupts the STAMP not a live term; DIM values kept (1450-1472)

```
$ sed -n '1450,1472p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
        emit(f"DIM_{name}", vector, key=key_name)
        emit(f"DIM_ROUTE_KIND_{name}", sp.Tuple(route_kind, route_operands), key=f"dim_route_kind_{name.lower()}")

    homogeneity = {
        "INPLANE_EOM": (M_DIM - 3 * L_DIM - 2 * T_DIM,) * len(sp.Add.make_args(sp.expand(MODEL["inplane"]))),
        "THICKNESS_EOM": (energy_density,) * len(sp.Add.make_args(sp.expand(MODEL["thickness"]))),
        "MASS_BALANCE": (mass_flux,) * len(sp.Add.make_args(sp.expand(MODEL["mass"]))),
        "AFFINITY": (specific_energy, pressure_4d - (M_DIM - 4 * L_DIM)),
        "CLOSURE": (mass_flux, dimensions["LAMBDA_A_0"][0] + specific_energy, dimensions["LAMBDA_V_0"][0] + L_DIM - T_DIM),
        "FACE_RESPONSE": (pressure_4d, dimensions["Z"][0] + L_DIM - T_DIM),
        "TWO_PORT_POWER_IDENTITY": (M_DIM - L_DIM - 3 * T_DIM,) * 3,
        "DISPERSION_DETERMINANT": (3 * M_DIM - 7 * L_DIM - 5 * T_DIM,),
    }
    for name, term_dimensions in homogeneity.items():
        emit(f"HOMOGENEITY_{name}", homogeneous_payload(term_dimensions), key=f"homogeneity_{name.lower()}")

    correct = homogeneous_payload(homogeneity["INPLANE_EOM"])
    corrupted_terms = list(homogeneity["INPLANE_EOM"])
    corrupted_terms[0] = corrupted_terms[0] + L_DIM
    corrupted = homogeneous_payload(corrupted_terms)
    restored = homogeneous_payload(homogeneity["INPLANE_EOM"])
    emit("HOMOGENEITY_ABLATION_DEMO", sp.Tuple(correct, corrupted, restored), key="homogeneity_ablation_demo")

```

### unattributed_* emitted as empty Tuple() (965-966)

```
$ sed -n '965,966p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S11b_interface_coupling_law_sympy_audit.py
    emit("UNATTRIBUTED_SINK_TERMS", sp.Tuple(), key="unattributed_sink_terms")
    emit("UNATTRIBUTED_EXCHANGE_TERMS", sp.Tuple(), key="unattributed_exchange_terms")
```

### Out of scope — the G13 criterion is an exported computed object (value must NOT change under Fix 1)

```
$ cd /var/projects/toy_physics/research/pde_ledger_v3/scripts && python3 -c "import S11b_exports as m; print('zperm_slice_map exported =', 'zperm_slice_map' in m.LEDGER)"
zperm_slice_map exported = True
```
