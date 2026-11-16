"""Review r1; from repo root: python3 witness.py CANDIDATE_ROOT SCRATCH_DIR.

Use the project guard. The candidate can be the frozen tree or a checkout of
the exact reviewed files. The YAMLs can come from cwd, checked by digest.
Only mutant subprocesses execute the candidate. Reference stdout is archived.
"""
from hashlib import sha256
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import sympy as s


candidate = Path(sys.argv[1]).resolve()
scratch = Path(sys.argv[2]).resolve()
scratch.mkdir(parents=True, exist_ok=True)
tool = Path("software/stage1_solver/tools/mouth_slip_traction.py")
body = Path("software/em_charge_attribute/u1_body_dynamics_inputs.yaml")
mechanics = Path("software/em_charge_attribute/u1_body_mechanics_inputs.yaml")
meas = Path("research/pde_ledger_v3/_measurements/mouth_slip_traction.txt")
record = Path("research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md")
sources = Path("research/pde_ledger_v3/_measurements/mouth_flow_and_motion_sources.txt")
pins = {
    tool: "f388b98c58aac58eecdab301f58d8610fc6ccc53d084922cde3e5365d798b84a",
    body: "eef5cdc2178ec76f10f3e5e361cdd1461eb2eee3980ffa99f814df02475866ee",
    mechanics: "6bb2ebe8ecf72ec59efba0f4d1929f97684b2c90f9e56857bc36f25d8d89e826",
    meas: "d034b06a58247111541cf9f848718b2b1db189dd4ee46f120f149383d63c0b61",
    record: "6059f3855ab6cf33b708341fe09a8800ed7c86f6b9550168ca1edcfc7cbeae0e",
    sources: "a47743bd342673e11f7856152bf7a9723c4f8c9d4c10c79ad60ecd55f1b0d227",
}
texts = {}
for path, digest in pins.items():
    actual = candidate / path
    if not actual.is_file() and path in (body, mechanics):
        actual = path
    data = actual.read_bytes()
    assert sha256(data).hexdigest() == digest, path
    texts[path] = data.decode()
print("reviewed candidate and its two independent YAML inputs: all six hashes match")
stdout = texts[meas].split("Literal stdout\n",1)[1].split("\n\n$",1)[0] + "\n"
assert sha256(stdout.encode()).hexdigest() == "9b2dc85787dcd27acc434a98cb5bc757561aaf78f1a03d6ee18887cac63c23d7"

K, gamma = s.symbols("K gamma", positive=True)
U,V,v,T,T_U,T_V = s.symbols("U V v T T_U T_V",real=True)
h = s.symbols("h", positive=True)
names = {str(x): x for x in [K,gamma,U,V,v,T,T_U,T_V]}
wall_force = gamma*(v-V)
fluid_force = -wall_force
# First row balances exterior and exchanged tractions. Second row imposes
# the relative-velocity law. Unknowns are the trace and wall force.
system = s.Matrix([[K,1],[-gamma,1]])
trace,force = map(s.factor, system.inv()*s.Matrix([T,-gamma*V]))
def odd(u, w):
    return s.factor((force.subs({T:(T_U*u+T_V*w)*h,V:w*h}) -
                     force.subs({T:-(T_U*u+T_V*w)*h,V:-w*h})) / (2*h))
print("two-unknown trace, wall force =", trace, force)
print("odd slopes [incident, translation, comoving] =", [odd(1,0),odd(0,1),odd(1,1)])
print("virtual power =", s.factor(wall_force*V+fluid_force*v))
first = re.match(r"Rayleigh density = (.+?)  applied traction = (.+?)  solved trace = (.+)\n",stdout)
assert first
for name,encoded,expected in zip(("Rayleigh density","traction","trace"),first.groups(),(gamma*(v-V)**2/2,T,trace)):
    print("archive and independent",name,"=",encoded,expected)
    assert s.factor(s.sympify(encoded,locals=names)-expected)==0
expected = {"tangential operator": K*v+wall_force,"forced trace in traction equation":T,
            "wall force before solving":wall_force,"solved wall force":force,
            "wall and material power":wall_force*V+fluid_force*v,
            "zero forcing trace per V":trace.subs(T,0)/V,
            "zero forcing wall force per V":force.subs(T,0)/V,
            "wall incident slope":odd(1,0),"wall translation slope":odd(0,1),"wall slope sum":odd(1,1)}
for label, value in expected.items():
    line, = [x for x in stdout.splitlines() if x.startswith(label+":")]
    m = re.fullmatch(r".+?: result=(.+?) reference=(.+?) residual=(.+)",line)
    actual, reference, residual = [s.sympify(x,locals=names) for x in m.groups()]
    print("archive and independent",label,"=",actual,s.factor(value))
    assert all(s.factor(x-value)==0 for x in (actual,reference)) and residual == 0

controls = [
    ("tangent omits wall motion", body, 'gammaSigma*(v.tangent-V.tangent)=0','gammaSigma*v.tangent=0','tangential operator'),
    ("Rayleigh action deleted", mechanics, 'rayleigh_integrand: "gammaSigma*(v_tangent-V_tangent)^2/2"', 'rayleigh_integrand: "0"', 'tangential operator'),
    ("nonquadratic slip", mechanics, '(v_tangent-V_tangent)^2/2"', '(v_tangent-V_tangent)^2*(1+v_tangent-V_tangent)/2"', 'wall and material power'),
    ("slip depends on absolute velocity sum", mechanics, '(v_tangent-V_tangent)^2/2"', '(v_tangent+V_tangent)^2/2"', 'wall force before solving'),
    ("wrong exterior operator in independent input", mechanics, 'tangent_DtN: "tangentDtN"', 'tangent_DtN: "tangentDtN*(1+V)"','tangential operator'),
    ("normal condition no longer follows wall", body, 'v.normal=V.normal; tangentDtN','v.normal=0; tangentDtN','normal motion on wall'),
    ("no essential normal trace", mechanics, 'essential: {normal: 1}\n    natural: [tangent, shear]\n    primitive_id: E5','essential: {}\n    natural: [tangent, shear]\n    primitive_id: E5','normal trace imposed'),
    ("aperture replaces wall domain", mechanics, 'rayleigh_surface: Sigma','rayleigh_surface: Aperture','Rayleigh surface'),
    ("corrupt forced trace, retaining zero traction branch",tool, 'surface_force = -sp.diff(loss, V)','trace += T*V\n    surface_force = -sp.diff(loss, V)','forced trace in traction equation'),
    ("wrong wall force sign",tool,'surface_force = -sp.diff(loss, V)','surface_force = sp.diff(loss, V)','wall force before solving'),
    ("wrong fluid force sign",tool,'force_on_fluid = -sp.diff(loss, v)','force_on_fluid = sp.diff(loss, v)','wall force before solving'),
    ("surface work omitted",tool,'surface_force * V + force_on_fluid * v','force_on_fluid * v','wall and material power'),
    ("solve force on archived branch only",tool,'surface_force.subs(v, trace)','surface_force.subs(v, trace.subs(T,0))','solved wall force'),
    ("swapped drive and drift slopes",tool,'for speed in (U, V)','for speed in (V, U)','wall incident slope'),
    ("zero reported translation slope",tool,'slope_sum = sp.factor(A + B)','B = sp.Integer(0)\n    slope_sum = sp.factor(A + B)','wall translation slope'),
    ("wrong slope sum alone",tool,'slope_sum = sp.factor(A + B)','slope_sum = sp.factor(A - B)','wall slope sum'),
    ("wrong zero-force drift readout alone",tool,'zero_force = sp.factor(force.subs(T, 0) / V)','zero_force = sp.Integer(0)','zero forcing wall force per V'),
    ("wrong archived matrix trace",body,'rhs: [2, "1+gammaSigma", "gammaSigma/(tangentDtN+gammaSigma)"]','rhs: [2, "1", "gammaSigma/(tangentDtN+gammaSigma)"]','zero forcing trace per V'),
]

for name, file, old, new, expected_failure in controls:
    with tempfile.TemporaryDirectory(dir=scratch) as temporary:
        for path in (tool,body,mechanics):
            text=texts[path]
            if path==file:
                assert text.count(old)==1, (name,old,text.count(old))
                text=text.replace(old,new)
            copy=Path(temporary)/path
            copy.parent.mkdir(parents=True,exist_ok=True)
            copy.write_text(text)
        result=subprocess.run([sys.executable,str(Path(temporary)/tool)],text=True,capture_output=True)
    print("\nMUTATION",name,"exit",result.returncode)
    print(result.stdout.rstrip())
    print(result.stderr.splitlines()[-1] if result.stderr else "stderr empty")
    # A parsing exception is not an effective physics control. Demand the
    # particular reported inconsistency and its emitted nonmatching operands.
    assert result.returncode==1 and "ArithmeticError: Inconsistent E5 wall calculation: " in result.stderr
    failed_labels = result.stderr.splitlines()[-1].split("calculation: ")[1].split(", ")
    assert expected_failure in failed_labels, (name,failed_labels)
    line, = [x for x in result.stdout.splitlines() if x.startswith(expected_failure+":")]
    assert not line.endswith(" residual=0"), (name,line)
    if "supplied=" in line:
        domain=re.search(r"supplied=(.+) required=(.+)",line)
        assert domain.group(1)!=domain.group(2)

print("\nRELEASE AND MOMENTUM")
m,u1,u2,t1,t2,rim,w=s.symbols("m u1 u2 t1 t2 rim w")
jump=m*(u1-u2)+t1-t2+rim
print("equal inflow rates, common-direction total momentum jump plus rim =", jump)
print("common tangential trace =",jump.subs(u2,u1),"| equal stress alone =",jump.subs(t2,t1))
print("change of observer for matched rates =",s.expand(jump.subs({u1:u1-w,u2:u2-w},simultaneous=True)-jump))
# Determine the first-order support-removal relation for a map with a genuine
# nonlinear term. The finite increments of its two geometric accelerations
# differ at second order even with exactly identical prepared independent data.
h1,h2,p,q,e=s.symbols("h1 h2 p q e")
F=s.Matrix([p,q]); H=e*s.Matrix([h1,h2]); net=F-H
acc=s.Matrix([[2,1],[1,3]])*net+s.Matrix([net[0]**2,net[0]*net[1]])
held={p:H[0],q:H[1]}; free={p:0,q:0}
chi=acc.jacobian(F).subs(free)
difference=(acc.subs(free)+chi*H).applyfunc(s.expand)
print("held acceleration =",acc.subs(held).T,"| free plus chi*H =",difference.T)
print("support-removal residual at order one =",difference.applyfunc(lambda x:x.coeff(e,1)).T)
# Positive volume rates, potentials oriented along laboratory +x.
d,Q1,Q2,M1,M2,A1,A2=s.symbols("d Q1 Q2 M1 M2 A1 A2",positive=True)
x=s.symbols("x",real=True)
u_at_left=s.diff(Q2/(4*s.pi*(d-x)),x).subs(x,0)
u_at_right=s.diff(Q1/(4*s.pi*x),x).subs(x,d)
print("left and right incident flows =",u_at_left,u_at_right)
print("conditional two Newton coefficients =",s.factor(d**2*A1*u_at_left/(M1*M2)),s.factor(-d**2*A2*u_at_right/(M2*M1)))
assert all(element.coeff(e,1)==0 for element in difference)
print(len(controls),"one-sided controls rejected for the expected reason")
