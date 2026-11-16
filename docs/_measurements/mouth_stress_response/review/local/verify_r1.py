"""Scoped R1 review of changed archive readers and the frozen numerical output."""
import contextlib
import hashlib
import json
import re
import sys
import traceback
from pathlib import Path

import numpy as np


root = Path.cwd()
local = Path(__file__).resolve().parent
review_home = local.parent
candidate = review_home/"r1/candidate"
scripts = candidate/"software/stage1_solver/tools"
measurements = candidate/"docs/_measurements/mouth_stress_response"
sys.path.insert(0,str(root/"software/stage1_solver/tools"))
for row in (candidate.parent/"candidate.sha256").read_text().splitlines():
    digest, name = row.split()
    assert hashlib.sha256((candidate/name).read_bytes()).hexdigest()==digest
    print("hash OK",name,digest)


def archived_face(manifest, run, speed):
    path = Path("_scratch/s11c")/run/"fields"/f"n32_restoring_slope_u{speed:g}.npz"
    records = json.loads((root/"docs/_measurements"/manifest).read_text())["runs"]
    expect, = [digest for record in records for nm,digest in record.get("fields",{}).items()
               if Path(nm)==path]
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expect
    outcome = json.loads((path.parent.parent/"outcome.json").read_text())
    assert outcome["limitsVerified"] and outcome["exitCode"]==outcome["childOutcome"]["exitCode"]==0
    with np.load(path) as p:
        count = round(1./p["dz"])
        assert count==64 and p["field"].shape==(count+1,32,32,32)
        face = p["field"][count]
        return face.real**2+face.imag**2, float(p["dx"])


n0, dx = archived_face("mouth_feed_weak_run.json","mouth_feed_depth01",0)
x = (np.arange(32)-15.5)*dx
axis = x[:,None,None]
envelope = np.maximum(1-axis**2-x[None,:,None]**2-x[None,None,:]**2,0)
# Fourier symbol for the knot derivative of the uniform periodic cubic spline.
modes = 2*np.pi*np.fft.fftfreq(32)
symbol = (3j/dx)*np.sin(modes)/(2+np.cos(modes))
sx = np.fft.ifft(np.fft.fft(n0,axis=0)*symbol[:,None,None],axis=0).real
readouts = (measurements/"force.stdout").read_text().splitlines()
for speed in (.01,.02):
    plus, spacing = archived_face("mouth_intake_redistribution.json","mouth_redistribution_depth01",speed)
    assert spacing==dx
    # Evaluate the weight derivative at x and at its parity-partner gridpoint;
    # series uses no virtual translate of either wavefunction or spline.
    linear = (plus-plus[::-1])/(2*speed)
    even = .5*(plus+plus[::-1])-n0
    value = -dx**3*np.sum(envelope**3*sx*linear + 6*axis*envelope**2*even*linear)
    row, = [line for line in readouts if line.startswith(f"|U|={speed:g} ")]
    reported = float(re.search(r"odd_force/U=(\S+)",row).group(1))
    print(f"independent Fourier spline on z=1 U={speed:g} value={value:.15g} frozen={reported:.12g}")
    assert abs(value-reported)<5e-17


def exercise(name,stem,entry,change=None,expect_error=None,nominal_output=None):
    src = (scripts/(stem+".py")).read_text()
    if change:
        before,after = change
        assert src.count(before)==1, name
        src = src.replace(before,after)
    directory = local/"probes"
    directory.mkdir(exist_ok=True)
    path = directory/(name+".py")
    path.write_text(src)
    out = directory/(name+".stdout")
    scope = {"__name__":"r1_probe","__file__":str(path)}
    err = None
    with out.open("w") as stream, contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
        try:
            exec(compile(src,str(path),"exec"),scope)
            scope[entry]()
        except Exception as exc:
            err = exc
            traceback.print_exc()
    if expect_error:
        assert err is not None and expect_error in str(err), (name,err)
        print(name,type(err).__name__,expect_error)
    else:
        assert err is None, (name,err)
        assert out.read_text().splitlines()==nominal_output
        print(name,"PASS; literal output matches frozen archive lines")


force = "mouth_stress_response"
cal = "mouth_density_calibration_check"
exercise("force_nominal",force,"archived_traces",nominal_output=readouts[5:])
cal_rows = (measurements/"calibration.stdout").read_text().splitlines()
exercise("receiver_nominal",cal,"receiver_packets",nominal_output=[s for s in cal_rows if s.startswith("open sleeve:")])
exercise("force_neighbor",force,"archived_traces",
         ("density = abs(field[-1])**2","density = abs(field[-2])**2"),"z=1 mouth plane")
exercise("force_reference_depth",force,"archived_traces",
         ("depths == 1.0","depths == 1.0-float(packet['dz'])"),"z=1 mouth plane")
exercise("force_state_drop_phase_part",force,"archived_traces",
         ("density = abs(field[-1])**2","density = field[-1].real**2"),"z=1 mouth plane")
for tag,old,new in [
        ("receiver_use_rate_for_bath",'bath_density = float(packet["bath_density"])','bath_density = float(packet["rate"])'),
        ("receiver_use_bath_for_rate",'rate = float(packet["rate"])','rate = float(packet["bath_density"])'),
        ("receiver_far_as_join",'int(packet["sleeve_cells"]) * width','int(1+packet["w_index"].max()) * width'),
        ("receiver_face_as_center",'packet["w_index"].min() + 0.5','packet["w_index"].min() + 0.0')]:
    exercise(tag,cal,"receiver_packets",(old,new),"Receiver preparation or geometry differs")
exercise("receiver_duplicate_prep",cal,"receiver_packets",
         ("expected_rates = {","files = [files[0],files[0]]\n    expected_rates = {"),
         "draining and zero-rate receiver preparations")
