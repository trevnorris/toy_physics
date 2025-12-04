import sympy

def find_scaling_solution_robust():
    # Variables
    a = sympy.Symbol('a', positive=True)
    M = sympy.Symbol('M', positive=True)

    print("Searching for scaling combinations that yield M ~ a...")

    # 1. Define Physics Options

    # Flow Geometry: Controls Velocity v = Flux / Area
    flow_geometries = {
        "4D_Throat (v ~ 1/a^3)": a**3,
        "3D_Hole (v ~ 1/a^2)": a**2
    }

    # Interaction Area: The surface the pressure pushes against to open the hole
    interaction_areas = {
        "Volume_Boundary (S ~ a^3)": a**3,
        "Surface_Area (S ~ a^2)": a**2
    }

    # Restoring Force: The 'Impedance' trying to close the hole
    restoring_forces = {
        "Vacuum_Pressure (F ~ a^2)": a**2,   # Force ~ Pressure * Area
        "Surface_Tension (F ~ a)": a,        # Force ~ Tension * Length (or Energy/Area * Length?)
                                             # Actually Surface Tension Force F ~ gamma * Length ~ a
        "Line_Tension (F ~ 1)": 1            # Force ~ Constant
    }

    matches = []

    for flow_name, area_flow in flow_geometries.items():
        for interact_name, area_interact in interaction_areas.items():
            for restore_name, f_close_scale in restoring_forces.items():

                # Logic:
                # F_open = P_dynamic * Area_interact
                # P_dynamic ~ v^2 ~ (M / area_flow)^2
                # F_open ~ (M^2 / area_flow^2) * area_interact

                # Equilibrium: F_open = F_close
                # (M^2 / area_flow^2) * area_interact = f_close_scale

                # Solve for M
                # M^2 = f_close_scale * area_flow^2 / area_interact
                m_squared = f_close_scale * area_flow**2 / area_interact
                m_scaling = sympy.sqrt(m_squared)

                # Extract Exponent k where M ~ a^k
                # k = d(log(M))/da * a
                exponent = sympy.simplify(sympy.diff(sympy.log(m_scaling), a) * a)

                if exponent == 1:
                    matches.append({
                        "Flow": flow_name,
                        "Area": interact_name,
                        "Restore": restore_name,
                        "Exponent": exponent
                    })
                elif abs(float(exponent) - 1.0) < 0.1: # Near match
                     print(f"NEAR MATCH (k={exponent}): {flow_name} + {interact_name} + {restore_name}")

    print(f"\nFound {len(matches)} Perfect Linear Matches (M ~ a):")
    for match in matches:
        print(f"  - Flow: {match['Flow']}")
        print(f"    Push: {match['Area']}")
        print(f"    Pull: {match['Restore']}")
        print("    -----------------------------")

find_scaling_solution_robust()
