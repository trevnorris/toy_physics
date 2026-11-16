"""Apply Grok's five source changes to disposable copies and separate its two stencil controls."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


def main():
    root = Path("software/stage1_solver/tools").resolve()
    sys.path.insert(0, str(root))
    from held_crossflow_budget import archived_packet
    from mouth_feed_weak_check import dx4

    cases = [
        ("mouth_confinement_projection", "work - r**2 * density", "work + r**2 * density", "ArithmeticError: Inconsistent angular"),
        ("held_crossflow_budget", "(pos - neg) / (2 * speed)", "(pos - neg) / (2 * speed**2)", "AssertionError:"),
        ("mouth_transmission_work", "translated != sp.diff(energy, y)", "translated != drive", "ArithmeticError: Interface momentum"),
        ("unequal_mouth_boundary_work", "force = -np.sum", "force = +np.sum", "AssertionError:"),
        ("held_crossflow_budget", "weight = np.maximum(1 - x*x - y*y - z*z, 0)**3", "weight = np.maximum(1 - x*x - y*y - z*z, 0)**2", "AssertionError:"),
    ]
    for number, (module, old, new, error) in enumerate(cases, 1):
        source = (root/(module+".py")).read_text()
        assert source.count(old) == 1
        with tempfile.TemporaryDirectory() as directory:
            copy = Path(directory)/(module+".py")
            copy.write_text(source.replace(old, new))
            run = subprocess.run([sys.executable, str(copy)], env=dict(os.environ, PYTHONPATH=str(root)),
                                 capture_output=True, text=True)
        assert run.returncode == 1 and error in run.stderr, (number, run.stderr[-300:])
        print(f"M{number} rejected by {error} (exit {run.returncode})")
        if number == 2:
            print(run.stdout.splitlines()[-1])

    path = archived_packet("mouth_intake_redistribution.json", "n32_restoring_slope_u0.01.npz")
    with np.load(path, allow_pickle=False) as data:
        n, mouth, dx = abs(data["field"][-1])**2, data["mouth"], float(data["dx"])
    coords = (np.arange(n.shape[0])+.5-n.shape[0]/2)*dx
    x,y,z = np.meshgrid(coords,coords,coords,indexing="ij",sparse=True)
    slope = -2*np.maximum(1-x*x-y*y-z*z,0)**3*(1-1/np.sqrt(n))
    for label, normal in [("compact", slope), ("extended", slope+.01*(~mouth))]:
        weight = dx**3/4
        full_parts = weight*np.sum(normal*dx4(n,dx))
        full_direct = -weight*np.sum(n*dx4(normal,dx))
        mask_parts = weight*np.sum(normal*dx4(n*mouth,dx))
        mask_direct = -weight*np.sum(n*mouth*dx4(normal,dx))
        alternate = -weight*np.sum(n*np.gradient(normal,dx,axis=0))
        print(f"M6 {label}: periodic full adjoint={full_parts:.15g} direct={full_direct:.15g}; "
              f"aperture adjoint={mask_parts:.15g} direct={mask_direct:.15g}; "
              f"fixed support, alternate stencil={alternate:.15g}")
        np.testing.assert_allclose([full_parts,mask_parts], [full_direct,mask_direct], rtol=0, atol=1e-14)
        assert abs(full_parts-alternate)>1e-11
    print("Changing slope support alone preserves full-plane periodic summation by parts.")


if __name__ == "__main__":
    main()
