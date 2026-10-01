#!/usr/bin/env python3
"""List the engine trace coordinates and derive the supplied potential pullback.

The algebra is the direct substitution of v_bulk=grad_4(phi),
delta_p=-rho_m*d_t(phi), the class-P Fourier factor, and the supplied bulk wave
equation.  No S11c engine or CAS is imported or run.
"""

from pathlib import Path


ROOT = Path("/var/projects/toy_physics")
INTERFACE = ROOT / "research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py"
SPEC = Path('/tmp/s11cd_clean_review_r4_claude/author_ablation/spec_corrupted.md')

formal = []
for face in ("plus", "minus"):
    formal.extend((f"delta_p_{face}", f"d_w_delta_p_{face}"))
    formal.extend(f"delta_v_bulk_{face}_{component}" for component in range(1, 5))
    formal.extend(f"d_w_delta_v_bulk_{face}_{component}" for component in range(1, 5))

for index, name in enumerate(formal, 1):
    print("FORMAL_TRACE_COORDINATE", index, name)

print("SUPPLIED_CLASS_P_ANSATZ", "phi=Phi(x1,w)*exp(i*(k2*x2-omega*t))")
print("SUPPLIED_TRACE_PRESSURE", "delta_p=i*rho_m*omega*Phi")
print("SUPPLIED_TRACE_VELOCITY", "(d1(Phi), i*k2*Phi, 0, Psi)", "Psi=d_w(Phi)")
print("SUPPLIED_NORMAL_JET_PRESSURE", "d_w(delta_p)=i*rho_m*omega*Psi")
print("SUPPLIED_NORMAL_JET_VELOCITY", "(d1(Psi), i*k2*Psi, 0, d_w(Psi))")
print("SUPPLIED_WAVE_EQUATION_PULLBACK", "d_w(Psi)=(k2^2-omega^2/c_s0^2)*Phi-d1^2(Phi)")
print("CLASS_P_ODD_BULK_TRACE_COMPONENT", 0)
print("CLASS_P_ODD_BULK_NORMAL_JET_COMPONENT", 0)

test_coordinates = (
    "delta_v_u_1",
    "delta_v_u_2",
    "delta_v_u_3",
    "delta_v_e_W",
    "delta_v_zeta_c",
)
for index, name in enumerate(test_coordinates, 1):
    print("VIRTUAL_TEST_COORDINATE", index, name)
print("VIRTUAL_OUTPUT_OBJECT", "FaceSource.virtual_displacement")
print("D_VIRTUAL_X3_D_TEST_U3", 1)
print("D_VIRTUAL_X3_D_PHYSICAL_U3", 0)
print("D_FACE_VELOCITY3_D_PHYSICAL_U3", "d_t")
print("VIRTUAL_MAP_IS_PHYSICAL_FRECHET_ROW", False)

interface_lines = INTERFACE.read_text(encoding="utf-8").splitlines()
for number, line in enumerate(interface_lines, 1):
    if any(
        marker in line
        for marker in (
            "# Traced bulk perturbations and their normal jets",
            "delta_v_bulk = {",
            "dw_delta_v_bulk = {",
            '1: symbol("d_w_delta_p_plus"',
            '-1: symbol("d_w_delta_p_minus"',
            "pressure_perturbation = affine_bulk_perturbation",
            "velocity_perturbation = tuple",
        )
    ):
        print("INTERFACE_SOURCE_LINE", number, line.strip())

for number, line in enumerate(SPEC.read_text(encoding="utf-8").splitlines(), 1):
    if "v_bulk=∇₄φ" in line:
        print("GOVERNING_SPEC_LINE", number, line.strip())
