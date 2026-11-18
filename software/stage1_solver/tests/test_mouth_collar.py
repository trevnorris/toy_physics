"""Checks for independent terminal calibration and signed material-packet playback."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from stage1_solver.mouth_collar import fit_terminal_load


def test_end_load_needs_two_observable_motions_and_a_passing_holdout():
    v = np.eye(2)
    f = np.array([[1.,0.2],[0.2,0.5]])
    check = np.array([[1.,1.],[-1.,1.]])
    with pytest.raises(ValueError,match="independent motions"):
        fit_terminal_load(v[:1],f[:1],check,check@f,1e-6)
    with pytest.raises(ValueError,match="fails calibration"):
        fit_terminal_load(v,f,check,2*check@f,1e-6)
    with pytest.raises(ValueError,match="fails calibration"):
        fit_terminal_load(v,-f,check,-check@f,1e-6)


def test_signed_material_packet_with_distributed_load(tmp_path):
    root = Path(__file__).resolve().parents[3]
    field = root/"_scratch/s11c/mouth_feed_depth01/fields/n32_restoring_slope_u0.npz"
    if not field.exists():
        pytest.skip("Requires the archived radius-one drive bundle")
    z = np.array([0,0.06,0.25,0.4,0.5])
    # Manufactured uniform accelerations of both fields. Applied force is
    # tabulated independently for all drives. The cut traction must vanish.
    histories = {}
    for name,a in (("zero",0),("negative",-0.02),("positive",0.02),
                   ("negative_half",-0.01),("positive_half",0.01)):
        force = np.broadcast_to(np.array([a,0.4*a])[:,None],(2,len(z)))
        histories[name] = dict(times=[0,0.001],values=[force.tolist()]*2)
    packet = dict(baseline_sha256=hashlib.sha256(field.read_bytes()).hexdigest(),inward_rate=0.25,
                  z=z.tolist(),cut_index=2,wall_line_inertia=1,wall_tension=1,
                  standing_amplitude=[1]*len(z),transverse_gradient_norm=1,
                  shear_inertia=2,shear_modulus=1,coupling_per_length=0,
                  terminal_readings=dict(train_velocity=[[1,0],[0,1]],train_force=[[0,0],[0,0]],
                                         check_velocity=[[1,1],[1,-1]],check_force=[[0,0],[0,0]],
                                         relative_tolerance=1e-6),distributed_force=histories)
    material = tmp_path/"fixture.json"
    material.write_text(json.dumps(packet))
    out = tmp_path/"out"
    result = subprocess.run([sys.executable,str(root/"software/stage1_solver/tools/signed_mouth_collar.py"),
                             "--material",str(material),"--output",str(out)],capture_output=True,text=True)
    assert result.returncode == 0, result.stderr
    for file in out.glob("*.npz"):
        with np.load(file) as row:
            a = float(row["acceleration"])
            expected = np.array([1,0.4])[:,None]*0.25*a*row["times"]
            np.testing.assert_allclose(row["stored"], expected,atol=1e-17)
            np.testing.assert_allclose(row["cut_flux"],0,atol=1e-17)
    assert len(list(out.glob("*.npz"))) == 5
