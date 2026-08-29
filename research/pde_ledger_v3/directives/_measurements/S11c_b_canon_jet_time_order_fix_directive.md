# Measurements backing `S11c_b_canon_jet_time_order_fix_directive.md`
Regenerated from the commands below (rule 2). Each block is a claim in the directive and the literal
stdout that verifies it. The comparator is now fixed, so the PRE-fix defect is shown via `git show HEAD:`.

## Claim: pre-fix canon_jet_name collapsed time-order via a Boolean `has_time`
$ git show HEAD:research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py | sed -n "798,812p"
    has_time = False
    directions: list[int] = []
    for token in derivatives:
        if token == "t":
            has_time = True
        elif token == "dw":
            directions.append(99)
        else:
            directions.extend(int(item) for item in re.findall(r"d(\d+)", token))
    directions.sort()
    suffix = "_t" if has_time else ""
    if directions:
        suffix += "_" + "".join("dw" if item == 99 else f"d{item}" for item in directions)
    return "_".join(base_parts) + suffix


## Claim: jet_suffix_from already PRESERVES the time count (so the collapse was solely in canon_jet_name)
$ sed -n "638,641p" research/pde_ledger_v3/scripts/S11c_a_cross_engine_comparator.py
            codes.extend([code] * order)
    time_codes = [code for code in codes if code == "t"]
    other_codes = [code for code in codes if code != "t"]
    return "_".join((*time_codes, *other_codes))

## Claim: PY spells order-1 time as _t and order-2 as a single tt token
$ grep -nE "_t\"|_tt\"" research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py | grep -iE "symbol" | head -6
214:    inherited_symbol(f"u_{a}_t", "COORDINATE", f"displacement velocity {a}", DIM_VELOCITY)
217:e_t = inherited_symbol("e_W_t", "COORDINATE", "thickness-fraction time derivative", -DIM_T)
249:    symbol(f"u_{a + 1}_tt", "COORDINATE", f"displacement acceleration {a + 1}", dim_div(DIM_L, 2 * DIM_T))
252:e_tt = symbol("e_W_tt", "COORDINATE", "thickness-fraction second time derivative", -2 * DIM_T)
270:uT_t = tuple(symbol(f"u_T_{a + 1}_t", "COORDINATE", f"curl-sector velocity probe {a + 1}", DIM_VELOCITY) for a in DIRECTIONS)
271:uL_t = tuple(symbol(f"u_L_{a + 1}_t", "COORDINATE", f"divergence-sector velocity probe {a + 1}", DIM_VELOCITY) for a in DIRECTIONS)

## Claim: the WL COUPLING carries 864x second-time (4th-slot order 2) transverse-trial derivatives -> COUPLING is NOT tt-free
$ grep -oE "Derivative\[[0-9], [0-9], [0-9], 2\]\[transverseTrialPotential[A-Za-z]+\]" research/pde_ledger_v3/mathematica/out/S11c_b_brane_operator_mathematica_audit.out | sort | uniq -c
    864 Derivative[0, 1, 1, 2][transverseTrialPotentialOne]
    864 Derivative[1, 0, 1, 2][transverseTrialPotentialTwo]
    864 Derivative[1, 1, 0, 2][transverseTrialPotentialThree]

## Claim: the pre-fix reconcile carried e_W_tt (inertia) residuals in SLAB_OPERATOR/THICKNESS_ROW
$ grep -n "e_W_tt" /home/trevnorris/.s11_build/S11c_b_reconcile_run.out | sed -E "s/(.{120}).*/\\1/" | head
68:A_minus_B SLAB_OPERATOR (OBJECT=THICKNESS_ROW, BRANCH=LAB_HELD, DENSITY=RHO4_CONSTANT, DOF=E_W) = Add(Mul(Symbol('C')
69:A_minus_B SLAB_OPERATOR (OBJECT=THICKNESS_ROW, BRANCH=LAB_HELD, DENSITY=RHOBR_CONSTANT, DOF=E_W) = Add(Mul(Integer(-1
70:A_minus_B SLAB_OPERATOR (OBJECT=THICKNESS_ROW, BRANCH=MATERIAL_ADVECTED, DENSITY=RHO4_CONSTANT, DOF=E_W) = Add(Mul(I,
71:A_minus_B SLAB_OPERATOR (OBJECT=THICKNESS_ROW, BRANCH=MATERIAL_ADVECTED, DENSITY=RHOBR_CONSTANT, DOF=E_W) = Add(Mul(I
