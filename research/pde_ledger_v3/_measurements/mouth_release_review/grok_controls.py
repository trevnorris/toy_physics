"""Exercise the two traction-preparation controls specified by the external referee."""
from pathlib import Path
import re
import sys
import sympy as sp

U, V, T_U, T_V = sp.symbols("U V T_U T_V", real=True)
K, gamma = sp.symbols("K gamma", positive=True)
symbols = {str(item):item for item in (U,V,T_U,T_V,K,gamma)}
report = Path(sys.argv[1]).read_text()

def reported(label):
    pattern = rf"^{label}: result=(.*?) reference=.*? residual=0$"
    return sp.sympify(re.search(pattern, report, re.M)[1], locals=symbols)

incident = reported("wall incident slope")
drift = reported("wall translation slope")
total = reported("wall slope sum")
for name, applied in (("relative-flow control", {T_U:K, T_V:0}),
                      ("changed incident forcing", {T_U:2*K, T_V:0})):
    # Insert the forcing into the referee's wall traction before differentiating.
    traction = gamma*(T_U*U + T_V*V - K*V)/(K+gamma)
    expected = [sp.diff(traction,U), sp.diff(traction,V)]
    values = [sp.factor(item.subs(applied)) for item in (incident,drift,total)]
    differences = [sp.factor(value-target.subs(applied))
                   for value,target in zip(values, [*expected, sum(expected)])]
    false_relation = sp.factor(values[1]-(-values[0]))
    print(name, "incident, drift, sum =", values, "wall residuals =", differences,
          "drift minus imposed opposite =", false_relation)
    assert differences == [0,0,0]
    assert (false_relation == 0) == (name == "relative-flow control")
