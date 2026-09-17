#!/usr/bin/env python3
"""T3: compact native Abel operand/normalization check; no production driver."""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE/'scripts/S11c_d_mixing_scattering_sympy_audit.py'
INPUT = BASE/'_measurements/S11c_d_variable_profile_development_input.json'
REPORT = BASE/'_measurements/S11_lean_analytic_error_source_checks.json'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    sources = [SOURCE, INPUT, BASE/'directives/S11c_d_SHARED_PHYSICS.md']
    before = {str(p.relative_to(BASE)): sha(p) for p in sources}
    tree = ast.parse(SOURCE.read_text())
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'EdgeReduction')
    init = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == '__init__')

    def assigns(node, field):
        return isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Attribute) and isinstance(t.value, ast.Name)
            and t.value.id == 'self' and t.attr == field for t in node.targets)

    start = next(i for i,n in enumerate(init.body) if assigns(n, 'regulator'))
    stop = next(i for i,n in enumerate(init.body) if assigns(n, 'distribution_record'))
    statements = init.body[start:stop+1]
    norm_start = next(i for i,n in enumerate(init.body) if assigns(n, 'regulated_ansatz'))-2
    norm_stop = next(i for i,n in enumerate(init.body) if assigns(n, 'normalization_operands'))
    normalization = init.body[norm_start:norm_stop+1]
    # Execute only the original half-line transform construction, with the
    # same coordinate/profile declarations. No source import or constructor.
    declaration_fields = ['xi', 'profiles']
    declarations = [next(n for n in init.body if assigns(n, f)) for f in declaration_fields]
    obj = SimpleNamespace()
    code = ast.Module(body=declarations+normalization+statements, type_ignores=[])
    ns = {'self': obj, 'sp': sp}
    exec(compile(ast.fix_missing_locations(code), str(SOURCE), 'exec'), ns)
    a,s = obj.regulator,obj.transfer
    # ell is positive; momentum transfer q is unrestricted real.
    ell = sp.Symbol('ell', positive=True)
    q = sp.Symbol('q', real=True)
    negative,positive = obj.abel_halves
    equalities = {
        'native_Fourier_mass': sp.simplify(obj.fourier_mass-2*sp.pi),
        'negative_half': sp.cancel(negative-1/(a-sp.I*s)),
        'positive_half': sp.cancel(positive-1/(a+sp.I*s)),
        'constant': sp.cancel(obj.abel_constant-2*a/(a*a+s*s)),
        'step_even': sp.cancel(obj.abel_even-a/(a*a+s*s)),
        'step_odd': sp.cancel(obj.abel_odd+sp.I*s/(a*a+s*s)),
    }
    h = a/ell
    poisson = h/(sp.pi*(h*h+q*q))
    conjugate = q/(sp.pi*(h*h+q*q))
    physical_constant = ell*obj.abel_constant.subs(s,ell*q)/obj.fourier_mass
    physical_step = ell*positive.subs(s,ell*q)/obj.fourier_mass
    equalities['constant_measure_and_width'] = sp.cancel(physical_constant-poisson)
    equalities['step_measure_and_PV_sign'] = sp.cancel(physical_step-(poisson-sp.I*conjugate)/2)
    assert all(v == 0 for v in equalities.values()),equalities
    supplied = json.loads(INPUT.read_text())
    assert sp.sympify(supplied['parameters']['L_W']) == 10
    controls = {
        'wrong_width_a_times_ell': sp.cancel(physical_constant-(a*ell)/(sp.pi*((a*ell)**2+q*q))) != 0,
        'wrong_PV_sign': sp.cancel(physical_step-(poisson+sp.I*conjugate)/2) != 0,
        'missing_step_half': sp.cancel(physical_step-(poisson-sp.I*conjugate)) != 0,
        'discard_constant': physical_constant != 0,
    }
    assert all(controls.values()),controls
    assert before == {str(p.relative_to(BASE)): sha(p) for p in sources}
    result = {'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
      'instrument_sha256':sha(Path(__file__)),'source_sha256':before,
      'native_AST_selection':{'class':'EdgeReduction','method':'__init__',
        'declarations':declaration_fields,'start_line':statements[0].lineno,'end_line':statements[-1].end_lineno,
        'Fourier_normalization_lines':[normalization[0].lineno,normalization[-1].end_lineno],
        'selection_sha256':hashlib.sha256(ast.dump(code,include_attributes=False).encode()).hexdigest()},
      'identities':{k:str(v) for k,v in equalities.items()},'controls':controls,
      'native_operands':{'negative_half':str(negative),'positive_half':str(positive),
        'constant':str(obj.abel_constant),'step_even':str(obj.abel_even),'step_odd':str(obj.abel_odd),
        'delta_mass':str(obj.distribution_record['DELTA_MASS'])},
      'parameter_map':{'s':'L_W * (k-k_prime)','physical_width':'a/L_W','approved_L_W':'10','approved_width':'a/10','origin':'xi=0'},
      'limits':['Exact selected native symbolic operands, not a kernel-certified CAS bridge.',
        'No production driver, closed kernel reconstruction or running exports touched.',
        'No proof of full folded-amplitude moments, kernel tails, outgoing inverse margin, quadrature errors or scattering convergence.',
        'Constant/step Fourier pairing with actual test/trial spaces remains an application obligation.']}
    REPORT.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: native Abel operands, measure/width/PV identity and four controls.')


if __name__ == '__main__':
    main()
