#!/usr/bin/env python3
"""Audit the actual S11c-a trace coordinates against directive B's map.

Only the Python standard library is used.  The algebra below is the direct
first-order composition of the face trace and of v_bulk=grad(phi).
"""

from pathlib import Path
import hashlib
import re


ROOT = Path("/var/projects/toy_physics")
SOURCE = ROOT / "research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py"
EXPORT = ROOT / "research/pde_ledger_v3/scripts/S11c_a_exports.py"
SPEC_B = ROOT / "research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md"
DIRECTIVE = ROOT / "research/pde_ledger_v3/directives/S11c_d_zinvariant_operator_blocks_directive.md"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def occurrences(path, needles):
    counts = {needle: 0 for needle in needles}
    overlap = max(map(len, needles)) - 1
    tail = ""
    with path.open("r", encoding="utf-8") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            data = tail + chunk
            safe = len(data) - overlap
            prefix = data[:safe] if safe > 0 else ""
            for needle in needles:
                counts[needle] += prefix.count(needle)
            tail = data[safe:] if safe > 0 else data
        for needle in needles:
            counts[needle] += tail.count(needle)
    return counts


print("SOURCE_SHA256", digest(SOURCE))
print("EXPORT_SHA256", digest(EXPORT))
print("DIRECTIVE_SHA256", digest(DIRECTIVE))

normal_trace_coordinates = ["d_w_delta_p_plus", "d_w_delta_p_minus"]
normal_trace_coordinates += [
    f"d_w_delta_v_bulk_{face}_{component}"
    for face in ("plus", "minus")
    for component in range(1, 5)
]
counts = occurrences(EXPORT, normal_trace_coordinates)
for name in normal_trace_coordinates:
    print("EXPORTED_NORMAL_TRACE_COORDINATE", name, "OCCURRENCES", counts[name])

directive_text = DIRECTIVE.read_text(encoding="utf-8")
listed_inputs = ["delta_p_s", "four components of v_bulk,s"]
print("DIRECTIVE_LISTED_BULK_INPUT_CLASSES", listed_inputs)
for name in normal_trace_coordinates:
    print("DIRECTIVE_NAMES_NORMAL_TRACE_INPUT", name, name in directive_text)

# The S11c-a ansatz evaluates a perturbation p + (w-s W0/2) p_w at the
# background face h0=s Wbg/2.  Delta_h0 is background order, so its product
# with the wave coordinate p_w belongs to the retained linear operator.
for face in (1, -1):
    print(
        "SHIFTED_PRESSURE_TRACE",
        "face", face,
        "= p_s +", f"({face}/2)*(W_bg-W_0)*d_w_p_s",
    )
    print(
        "D_SHIFTED_PRESSURE_TRACE_D_DW_P",
        "face", face,
        "=", f"({face}/2)*(W_bg-W_0)",
    )
    for component in range(1, 5):
        print(
            "D_SHIFTED_VELOCITY_TRACE_D_DW_V",
            "face", face,
            "component", component,
            "=", f"({face}/2)*(W_bg-W_0)",
        )

# Apply the supplied potential-flow identities on P:
# phi=Phi(x1,w) exp(i(k2*x2-omega*t)), d3 phi=0.
print("POTENTIAL_TRACE_ANSATZ", "phi=Phi(x1,w)*exp(i*(k2*x2-omega*t))")
print("DELTA_P_TRACE", "i*rho_m*omega*phi")
print("V_BULK_1_TRACE", "d1(phi)")
print("V_BULK_2_TRACE", "i*k2*phi = k2*delta_p/(rho_m*omega)")
print("V_BULK_3_TRACE", 0)
print("V_BULK_4_TRACE", "d_w(phi)")
print("INDEPENDENT_FOUR_VELOCITY_COMPONENTS_ON_P", False)

# The face map is R=(X+u,h).  Its virtual displacement is a test variation,
# not a physical output coordinate.  Its physical Frechet derivative is zero,
# whereas the actual face velocity has the expected physical derivative.
print("FACE_MAP", "R_s=(X1+u1,X2+u2,X3+u3,h_s)")
print("VIRTUAL_FACE_DISPLACEMENT", "delta_v_R_s=(delta_v_u1,delta_v_u2,delta_v_u3,delta_v_h_s)")
print("D_DELTAV_R3_D_U3", 0)
print("D_VFACE3_D_U3", "d_t")
print("DELTA_V_X_IS_PHYSICAL_OUTPUT_OF_SAME_FRECHET_MAP", False)

# Show the exact source declaration lines for the normal trace coordinates.
source_lines = SOURCE.read_text(encoding="utf-8").splitlines()
patterns = (
    "# Traced bulk perturbations and their normal jets",
    "dw_delta_v_bulk =",
    "d_w_delta_p_plus",
    "pressure_perturbation = affine_bulk_perturbation",
    "velocity_perturbation = tuple",
    "def traction_raw",
    "def closure_raw",
)
for number, line in enumerate(source_lines, 1):
    if any(pattern in line for pattern in patterns):
        print("SOURCE_LINE", number, line.strip())

# Show the supplied bulk identity directly from the governing spec.
for number, line in enumerate(SPEC_B.read_text(encoding="utf-8").splitlines(), 1):
    if "v_bulk=∇₄φ" in line:
        print("GOVERNING_SPEC_LINE", number, line.strip())
