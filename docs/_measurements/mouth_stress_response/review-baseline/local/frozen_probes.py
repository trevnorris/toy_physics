"""Independent archived-field postprocessing and focused mutations of frozen copies."""
from contextlib import redirect_stdout, redirect_stderr
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback

import numpy as np
import sympy as sp

base = Path("_scratch/s11c/mouth_stress_review01/baseline")
case = base/"candidate"
local = Path(__file__).resolve().parent
frozen_tools = case/"software/stage1_solver/tools"
sys.path.insert(0, str(Path("software/stage1_solver/tools").resolve()))
for line in (base/"candidate.sha256").read_text().splitlines():
    digest, name = line.split()
    assert sha256((case/name).read_bytes()).hexdigest() == digest


def read_face(run, manifest, speed):
    filename = f"n32_restoring_slope_u{speed:g}.npz"
    path = Path("_scratch/s11c")/run/"fields"/filename
    md = json.loads((Path("docs/_measurements")/manifest).read_text())
    digest, = [v for row in md["runs"] for k,v in row.get("fields",{}).items()
               if Path(k) == path]
    assert sha256(path.read_bytes()).hexdigest() == digest
    outcome = json.loads((path.parent.parent/"outcome.json").read_text())
    assert outcome["limitsVerified"] and outcome["exitCode"] == outcome["childOutcome"]["exitCode"] == 0
    with np.load(path) as pack:
        fld = pack["field"]
        assert fld.shape == (65,32,32,32)
        trace = np.square(np.abs(fld[64]))
        feed = np.square(np.abs(fld[0]))
        z_end = 64*pack["dz"]
        rate = pack["dx"]**3 * np.sum(pack["current"])
        print(f"feed U={speed:g}: range=({feed.min():g},{feed.max():g}); last_z={z_end:g};"
              f" mouth n range=({trace.min():.12g},{trace.max():.12g}); inward rate={rate:.12g}")
        assert z_end == 1 and feed.min() == feed.max() == 1 and abs(rate-.25)<1e-12
    return trace


ref = read_face("mouth_feed_depth01","mouth_feed_weak_run.json",0)
N = len(ref)
dx = 4/N
x = (np.arange(N)+.5)*dx-2
X,Y,Z = np.meshgrid(x,x,x,indexing="ij",sparse=True)
p = np.maximum(1-(X**2+Y**2+Z**2),0)
# At uniform periodic knots, cubic spline derivatives satisfy this cyclic system.
identity = np.eye(N)
matrix = 4*identity + np.roll(identity,1,axis=1) + np.roll(identity,-1,axis=1)
rhs = (3/dx)*(np.roll(ref,-1,axis=0)-np.roll(ref,1,axis=0))
grad = np.linalg.solve(matrix,rhs.reshape(N,-1)).reshape(ref.shape)
for speed in [.01,.02]:
    plus = read_face("mouth_redistribution_depth01","mouth_intake_redistribution.json",speed)
    # Split into even and odd densities. The F'-term from translating h involves
    # the product of the odd density and the finite-speed even departure.
    symmetry = (plus+plus[::-1])/2
    odd = (plus-plus[::-1])/(2*speed)
    a = -dx**3 * np.sum(p**3*odd*grad)
    b = -6*dx**3 * np.sum(X*p**2*odd*(symmetry-ref))
    print(f"archive by cyclic spline/even-odd density, U={speed:g}, ref={a:.12g}, weight={b:.12g}, total={a+b:.12g}")
    assert abs(a+b-({.01:-1.69608024596e-5,.02:-1.69606371567e-5}[speed]))<1e-16

for target in ["open_sleeve_current01","open_sleeve_zero01"]:
    pth = Path("_scratch/s11c")/target/"fields/level_32.npz"
    assert sha256(pth.read_bytes()).hexdigest()+"  "+str(pth) in Path("docs/_measurements/open_sleeve_current.txt").read_text()
    with np.load(pth) as pck:
        w = (pck["w_index"]+.5)*pck["dw"]
        print("sleeve",target,"bath",pck["bath_density"],"rate",pck["rate"],"wmin/max",w.min(),w.max(),
              "join",pck["dw"]*pck["sleeve_cells"],"K",pck["eos_k"])

# An explicit smooth surface-energy law tested at varying rest control f.
r,f,n = sp.symbols("r f n", positive=True)
loc = 1+r*r + (1+r)*f + f*f
law = r*n**2
control_ratio = (sp.diff(law.subs(n,loc),f)/sp.diff(loc,f)).subs(f,0)
print("rest feed derivative / local density derivative:",sp.factor(control_ratio))
print("rest feed derivative alone:",sp.factor(sp.diff(law.subs(n,loc),f).subs(f,0)))
assert sp.simplify(control_ratio-sp.diff(law,n).subs(n,loc.subs(f,0)))==0


def ablate(filename, tag, entry, original, replacement, expect_pass):
    src = (frozen_tools/filename).read_text()
    assert src.count(original)==1, (tag,original)
    changed = src.replace(original,replacement)
    folder = local/"mutants"
    folder.mkdir(exist_ok=True)
    script = folder/(tag+".py")
    script.write_text(changed)
    scope = {"__name__":"disposable_probe","__file__":str(script)}
    failure = None
    out = folder/(tag+".stdout")
    with out.open("w") as data, redirect_stdout(data), redirect_stderr(data):
        try:
            exec(compile(changed,str(script),"exec"),scope)
            scope[entry]()
        except Exception as exc:
            failure = type(exc).__name__
            traceback.print_exc()
    status = failure or "PASS"
    print(f"{tag} {status}")
    assert expect_pass == (failure is None), (tag, status)

force = "mouth_stress_response.py"
for tag,old,new in [
        ("angular_tangent_removed","+ (1-mu**2)/r*sp.diff(jump, mu)","+ 0"),
        ("edge_deleted","boundary = weight*r**2*n0*j1","boundary = 0"),
        ("edge_control_monopole","n1:0, j0:0, j1:r","n1:0, j0:0, j1:r**2"),
        ("jump_orientation","jump1 = 2*kappa/c","jump1 = -2*kappa/c"),
        ("anchor_frozen","ref = 1+b*radius2","ref = 1+b*(x**2+t)/a**2")]:
    ablate(force,tag,"angular_projection",old,new,False)
for tag,old,new,p in [
        ("weight_translation_omitted",".5*hx*delta**2 - h*delta*ref_x","-h*delta*ref_x", False),
        ("reflection_wrong_axis","neg = pos[::-1]","neg = pos[:,::-1]",False),
        ("field_plane_off_by_one",'packet["field"][-1]','packet["field"][-2]',True)]:
    ablate(force,tag,"archived_traces",old,new,p)

cal = "mouth_density_calibration_check.py"
for tag,entry,old,new,expect in [
        ("center_anchor_removed","finite_density_difference","gamma*q*(q-h)*(q+h)","gamma*(q-h)*(q+h)", False),
        ("finite_secant_one_sided","finite_density_difference","second.subs(q,-h)","second.subs(q,0)",False),
        ("size_of_feed_domain","feed_packets",'(len(field)-1)*float(packet["dz"])','len(field)*float(packet["dz"])',False),
        ("feed_read_on_mouth","feed_packets","abs(field[0])","abs(field[-1])",False),
        ("sleeve_bath_to_rate","receiver_packets","float(packet['bath_density'])","float(packet['rate'])",True),
        ("sleeve_shoulder_to_far_end","receiver_packets","int(packet['sleeve_cells'])*width","(packet['w_index'].max()+1)*width",True)]:
    ablate(cal,tag,entry,old,new,expect)
