"""Challenge the archive selection and force and secant identities on disposable scripts."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    tools = Path("software/stage1_solver/tools").resolve()
    env = {**os.environ, "PYTHONPATH": str(tools)}
    force = "mouth_stress_response.py"
    calibration = "mouth_density_calibration_check.py"
    cases = (
        ("face one plane early", force, "abs(field[-1])**2", "abs(field[-2])**2",
         "Density is not from the z=1 mouth plane"),
        ("edge contribution removed", force, "boundary = weight*r**2*n0*j1", "boundary = 0",
         "Angular differentiation and radial reduction disagree"),
        ("receiving jump reversed", force, "jump1 = 2*kappa/c", "jump1 = -2*kappa/c",
         "assert stress_force != 0"),
        ("bath read from rate", calibration, 'bath_density = float(packet["bath_density"])',
         'bath_density = float(packet["rate"])', "Receiver preparation or geometry differs"),
        ("inlet rate read from bath", calibration, 'rate = float(packet["rate"])',
         'rate = float(packet["bath_density"])', "Receiver preparation or geometry differs"),
        ("far face used as join", calibration, 'int(packet["sleeve_cells"]) * width',
         '(packet["w_index"].max()+1) * width', "Receiver preparation or geometry differs"),
        ("face used instead of first center", calibration, '(packet["w_index"].min() + 0.5) * width',
         'packet["w_index"].min() * width', "Receiver preparation or geometry differs"),
        ("drop one preparation", calibration, 'for _, name in files]', 'for _, name in files[:1]]',
         "Expected the draining and zero-rate receiver preparations"),
        ("secant uses center twice", calibration, "second.subs(q,-h)", "second.subs(q,h)",
         "and secant == slope"),
    )
    for label, filename, original, replacement, expected in cases:
        text = (tools / filename).read_text()
        if text.count(original) != 1:
            raise ValueError(f"Could not locate the one-sided change for {label}")
        with tempfile.TemporaryDirectory(prefix="mouth-stress-review-") as directory:
            copy = Path(directory) / filename
            copy.write_text(text.replace(original, replacement))
            result = subprocess.run([sys.executable, str(copy)], env=env, capture_output=True, text=True)
        if result.returncode != 1 or expected not in result.stderr:
            raise ArithmeticError(f"{label}: exit {result.returncode}, stderr {result.stderr[-500:]}")
        print(f"{label}: rejected ({expected})", flush=True)


if __name__ == "__main__":
    main()
