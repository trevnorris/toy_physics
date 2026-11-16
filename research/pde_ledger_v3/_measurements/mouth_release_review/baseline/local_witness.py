"""Independent algebra and one-sided controls for the r0 two-mouth review.

Run with the repository as cwd under scripts/s11c_guarded_run.py. The first
argument is a frozen tree (or the checkout of the candidate). The second is
the directory for disposable copies. The baseline stdout is read as recorded.
"""
from hashlib import sha256
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import sympy as s
import yaml


tool = Path("software/stage1_solver/tools/mouth_slip_traction.py")
dynamics = Path("software/em_charge_attribute/u1_body_dynamics_inputs.yaml")
mechanics = Path("software/em_charge_attribute/u1_body_mechanics_inputs.yaml")
measurement = Path("research/pde_ledger_v3/_measurements/mouth_slip_traction.txt")
frozen = Path(sys.argv[1]).resolve()
copies = Path(sys.argv[2]).resolve()
copies.mkdir(parents=True, exist_ok=True)
texts = {}
for path, expected in (
    (tool, "0fbea23c3426338f71b39830fd3d759b73f25d6eefae6074717fe491617a31a9"),
    (dynamics, "eef5cdc2178ec76f10f3e5e361cdd1461eb2eee3980ffa99f814df02475866ee"),
    (mechanics, "6bb2ebe8ecf72ec59efba0f4d1929f97684b2c90f9e56857bc36f25d8d89e826"),
    (measurement, "e8f4fa14a12a35d55db2542877dfea725d000fe154f0e1b575fd642ffbde62fa"),
):
    source = frozen / path if (frozen / path).exists() else path
    data = source.read_bytes()
    assert sha256(data).hexdigest() == expected, path
    texts[path] = data.decode()

K, gamma = s.symbols("K gamma", positive=True)
v, V, U, T, T_U, T_V = s.symbols("v V U T T_U T_V", real=True)
names = {str(x): x for x in [K, gamma, v, V, U, T, T_U, T_V]}
R = gamma*(v - V)**2/2
exterior = K*v**2/2 - T*v
# Stationarity of the total quadratic functional. Derive the force from
# the exterior momentum balance, and its slopes by implicit differentiation.
equation = s.diff(exterior + R, v)
vstar, = s.solve(equation, v)
wall = s.expand(-s.diff(exterior, v))
load = s.factor(wall.subs(v, vstar))
rest = {v: 0, V: 0, U: 0, T: 0}
direct = s.factor(T_U - K*(T_U - s.diff(equation, U))/s.diff(equation, v))
drift = s.factor(T_V - K*(T_V - s.diff(equation, V))/s.diff(equation, v))
fv, fV = -s.diff(R, v), -s.diff(R, V)
power = s.factor(fv*v+fV*V)
print("stationarity =", equation, "| wall from exterior =", wall, "| solved =", load)
print("force on material =", fv, "| force on wall =", fV, "| power =", power)
print("implicit responses =", direct, drift, "| sum =", s.factor(direct+drift))

body = yaml.safe_load(texts[dynamics])
e5 = body["endpoints"]["E5"]["trace_system"]
ynames = {"gammaSigma": gamma, "tangentDtN": K}
matrix = s.Matrix([[s.sympify(i, locals=ynames) for i in row] for row in e5["matrix"]])
rhs = s.Matrix([s.sympify(i, locals=ynames) for i in e5["rhs"]])
normal, archive_trace, shear = map(s.factor, matrix.inv()*rhs)
print("separate archived 3x3 normal tangent shear =", normal, archive_trace, shear)
assert s.factor(vstar.subs(T, 0)/V-archive_trace) == 0
assert s.simplify(load - fV.subs(v, vstar)) == 0
assert s.factor(power+2*R) == 0

def residuals(output):
    def parse(pattern):
        m = re.search(pattern, output, re.MULTILINE)
        assert m, pattern
        return [s.sympify(a.strip(), locals=names) for a in m.groups()]
    op, = parse(r"^tangential operator = (.+?)  forcing = .+$")
    rayleigh, pow_report = parse(r"^Rayleigh density = (.+?)  summed force power = (.+)$")
    tr, _ = parse(r"^solved trace = (.+?)  equation residual = (.+)$")
    wall_report, solved = parse(r"^force on surface = (.+?)  solved force = (.+)$")
    tr_zero, force_zero = parse(r"^zero forcing: trace per V = (.+?)  force per V = (.+)$")
    a, b, sum_out = parse(r"^surface force slopes: A = (.+?)  B = (.+?)  A\+B = (.+)$")
    # Cross the independently supplied tangential equation and Rayleigh force.
    # Separate controls also detect common-mode sign reversals or wrong exterior.
    tests = {
        "surface/normal": int("surface = Sigma  normal condition = v.normal=V.normal\n" not in output),
        "input matching": op - K*v - s.diff(rayleigh, v),
        "rayleigh": rayleigh-R, "operator": op-(equation+T),
        "trace": tr-vstar, "surface": wall_report-fV, "load": solved-load,
        "zero trace": tr_zero-archive_trace, "zero load": force_zero-load.subs(T,0)/V,
        "power": pow_report-power, "A": a-direct, "B": b-drift,
        "sum": sum_out-direct-drift,
    }
    return {k:s.factor(x) for k,x in tests.items() if s.factor(x)!=0}

recorded = texts[measurement].split("Literal stdout\n")[1].split("\n\n$")[0] + "\n"
checked = residuals(recorded)
print("r0 archived stdout discrepancies =", checked)
assert not checked

tests = [
    ("tangent drops wall velocity", dynamics, "gammaSigma*(v.tangent-V.tangent)=0", "gammaSigma*v.tangent=0", 0, "input matching"),
    ("Rayleigh deleted", mechanics, 'rayleigh_integrand: "gammaSigma*(v_tangent-V_tangent)^2/2"', 'rayleigh_integrand: "0"', 0, "input matching"),
    ("cubic slip term", mechanics, '(v_tangent-V_tangent)^2/2"', '(v_tangent-V_tangent)^2*(1+v_tangent-V_tangent)/2"', 1, "power"),
    ("Rayleigh depends on velocity sum", mechanics, '(v_tangent-V_tangent)^2/2"', '(v_tangent+V_tangent)^2/2"', 1, "surface"),
    ("independent normal unbound", dynamics, "v.normal=V.normal; tangentDtN", "v.normal=0; tangentDtN", 0, "surface/normal"),
    ("assigned to aperture", mechanics, "rayleigh_surface: Sigma", "rayleigh_surface: aperture", 0, "surface/normal"),
    ("bad forced trace", tool, "surface_force = -sp.diff(loss, V)", "trace += T*V\n    surface_force = -sp.diff(loss, V)", 1, "trace"),
    ("force at T=0 used for all T", tool, "surface_force.subs(v, trace)", "surface_force.subs(v, trace.subs(T, 0))", 0, "load"),
    ("interchange incident and drift slopes", tool, "for speed in (U, V)", "for speed in (V, U)", 0, "A"),
    ("wrong sum only", tool, 'sp.factor(A+B))', 'sp.factor(A-B))', 0, "sum"),
    ("force_per_V printed as zero", tool, "sp.factor(force.subs(T, 0) / V)", "sp.factor(force.subs({T: 0, V: 0}) / V)", 0, "zero load"),
]
for title, file, old, new, expect, witness in tests:
    with tempfile.TemporaryDirectory(dir=copies) as directory:
        tree = Path(directory)
        for path in (tool, dynamics, mechanics):
            data = texts[path]
            if path == file:
                assert data.count(old) == 1, (title,path,old,data.count(old))
                data = data.replace(old, new)
            (tree / path).parent.mkdir(parents=True, exist_ok=True)
            (tree / path).write_text(data)
        result = subprocess.run([sys.executable, str(tree/tool)], capture_output=True, text=True)
    checked = residuals(result.stdout)
    print("\nMUTATION:", title, "| candidate exit =", result.returncode)
    print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.splitlines()[-1])
    print("independent discrepancies =", checked)
    assert result.returncode == expect and witness in checked

# Independent 0PN balance, kinematics and release checks. Positive mass rate
# crosses both traces in the common feeding-to-receiving direction.
m, vm, vp, sigm, sigp, rim, W = s.symbols("m vm vp sigm sigp rim W")
def reaction(in_velocity, out_velocity):
    return m*in_velocity+sigm-m*out_velocity-sigp+rim
jump = reaction(vm, vp)
print("\nphysical opening signed momentum jump + rim =", jump)
print("observer boost difference, same crossing mass rate =", s.expand(reaction(vm-W,vp-W)-jump))
print("common phase and equal normal rate -> carried jump =", (m*(vm-vp)).subs(vp,vm))
print("equal stresses, no rim, unequal tangential velocities -> force =", jump.subs({sigm:sigp,rim:0}))
body_speed, incoming, reservoir = s.symbols("body_speed incoming reservoir")
print("boost invariants (U,V) =", tuple(s.expand(x) for x in ((incoming-W)-(reservoir-W),(body_speed-W)-(reservoir-W))))
A, B = s.symbols("A B")
relative_law = A*(U-V)
print("relative-only slopes and isolated mechanical power =", s.diff(relative_law,U), s.diff(relative_law,V), (relative_law*V).subs(U,0))

H1,H2 = s.symbols("H1 H2")
H=s.Matrix([H1,H2])
f1,f2=s.symbols("f1 f2")
f=s.Matrix([f1,f2])
# Positive definite unequal inertia including a nonzero cross term.
inertia=s.Matrix([[3,1],[1,2]])
alpha=inertia.inv()*(f-H)
chi=alpha.jacobian(f)
print("coupled release chi =", chi, "| removal acceleration =", alpha.subs({f1:0,f2:0}).T)
print("held acceleration =", alpha.subs({f1:H1,f2:H2}).T, "| release + chi H =", (alpha.subs({f1:0,f2:0})+chi*H).T)
added=inertia-s.eye(2)
print("incremental constraint load at release for H=(-1,+1), f=0 =", (-added*alpha).subs({H1:-1,H2:1,f1:0,f2:0}).T)
d,Q1,Q2,M1,M2,A1,A2=s.symbols("d Q1 Q2 M1 M2 A1 A2",positive=True)
x=s.symbols("x",positive=True)
velocity=-s.diff(1/x,x)  # magnitude per Q/(4*pi)
G1=d**2/M2*(A1/M1)*Q2/(4*s.pi)*velocity.subs(x,d)
G2=d**2/M1*(A2/M2)*Q1/(4*s.pi)*velocity.subs(x,d)
print("two Newton comparisons =", s.factor(G1), s.factor(G2))
