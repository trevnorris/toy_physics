import sympy

def analyze_throat_scaling():
    # Variables
    a = sympy.Symbol('a', positive=True) # Throat radius
    n = sympy.Symbol('n', positive=True) # Polytropic index (P ~ rho^n)

    # 1. Geometric Scaling (4D Bulk)
    # The throat cross-section is a 3D ball of radius a.
    # Volume/Area for flux ~ a^3
    Area_4D = a**3

    # 2. Soliton Density Scaling
    # The "healing length" xi scales such that gradient energy ~ potential energy.
    # 1/a^2 ~ rho^(n-1)  =>  rho_core ~ a^(-2/(n-1))
    rho_core = a**(-2/(n - 1))

    # 3. Velocity Scaling (Choked Flow)
    # Velocity is limited by sound speed c_s at the throat.
    # c_s^2 ~ rho^(n-1)  =>  c_s ~ a^-1
    v_sonic = a**-1

    # 4. Mass Flux Scaling
    # M ~ rho * Area * v
    M_scaling = rho_core * Area_4D * v_sonic

    # Calculate Exponent: M ~ a^k
    # k = -2/(n-1) + 3 - 1
    exponent_M = sympy.simplify(-2/(n - 1) + 2)

    print(f"General Mass Scaling: M ~ a^({exponent_M})")

    # Check specific Equations of State
    print(f"For n=5 (Stiff Vacuum): M ~ a^{exponent_M.subs(n, 5)}")
    print(f"For n=3 (Intermediate): M ~ a^{exponent_M.subs(n, 3)}")
    print(f"For n=2 (BEC):          M ~ a^{exponent_M.subs(n, 2)}")

analyze_throat_scaling()
