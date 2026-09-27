#!/usr/bin/env python3
"""Selected finite-solve AST operations on synthetic inputs only; host guard required.
No production module imports, constructors, saved operands or physical solves.
"""
import ast
import copy
import hashlib
import json
import platform
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
M = BASE / '_measurements'
PATHS = [M/'S11c_d_finite_scattering.py', M/'S11c_d_finite_scattering_resolution.py',
         M/'S11c_d_finite_scattering_domain.py', BASE/'directives/S11c_d_SHARED_PHYSICS.md',
         BASE/'directives/S11c_d_SCATTERING_FORM_AMENDMENT.md']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    before = {str(p.relative_to(BASE)): sha(p) for p in PATHS}
    tree = ast.parse(PATHS[0].read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'construct')
    def begins(node, name):
        return isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name) and node.targets[0].id == name
    def block(start, end):
        lo = next(i for i, n in enumerate(function.body) if begins(n, start))
        hi = next(i for i, n in enumerate(function.body) if begins(n, end))
        assert lo <= hi
        return copy.deepcopy(function.body[lo:hi+1])
    operations = {'balance_solve_unscale': block('row_scale', 'scaled_residual'),
                  'observe_current': block('values', 'incoming_flux'),
                  'residual': block('residual', 'scaled_residual')}
    # Only selected statements, never the enclosing constructor or its imports.
    for nodes in operations.values():
        assert not any(isinstance(n, (ast.Import, ast.ImportFrom)) for stmt in nodes for n in ast.walk(stmt))
        assert not any(isinstance(n, ast.Name) and n.id in ['atomic_pickle', 'load', 'open']
                       for stmt in nodes for n in ast.walk(stmt))
    ast_records = {k: {'sha256': hashlib.sha256(ast.dump(ast.Module(body=v, type_ignores=[])).encode()).hexdigest(),
                       'lines': [v[0].lineno, v[-1].end_lineno], 'statements': len(v)} for k,v in operations.items()}
    checks = []
    def check(name, actual, expected):
        a, b = np.asarray(actual), np.asarray(expected)
        assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), name
        error = float(np.max(abs(a-b))) if a.size else 0.
        scale = max(1., float(np.max(abs(b))) if b.size else 0.)
        assert error <= 1e-12*scale, (name, error, scale)
        checks.append({'name': name, 'kind': 'identity/example', 'passed': True, 'max_abs_error': error,
                       'tolerance': 1e-12*scale})
    def reject(name, actual, wrong):
        a,b = np.asarray(actual),np.asarray(wrong)
        delta = float(np.max(abs(a-b)))
        assert np.isfinite(delta) and delta > 1e-8, (name, delta)
        checks.append({'name': name, 'kind': 'wrong-formula control', 'passed': True, 'separation': delta})
    def bounded(name, value, bound):
        assert np.isfinite(value) and np.isfinite(bound) and value <= bound+1e-12, (name,value,bound)
        checks.append({'name': name, 'kind': 'conditional bound example', 'passed': True,
                       'actual': float(value), 'bound': float(bound)})
    def require(condition, message):
        if not condition: raise AssertionError(message)
    def execute(key, env, nodes=None):
        module = ast.fix_missing_locations(ast.Module(body=copy.deepcopy(operations[key] if nodes is None else nodes),type_ignores=[]))
        exec(compile(module, '<selected finite '+key+'>', 'exec'), env)

    size = 2
    matrix = np.kron(np.eye(5), np.array([[2,1],[0,4]],complex))
    inverse = np.kron(np.eye(5), np.array([[.5,-.125],[0,.25]],complex))
    exact = (np.arange(20).reshape(10,2)-8)/8 + .25j
    rhs = matrix@exact
    env = {'np': np, 'require': require, 'matrix': matrix, 'rhs': rhs}
    execute('balance_solve_unscale', env)
    R = np.diag(1/env['row_scale']); D = np.diag(env['column_scale'])
    check('synthetic_exact_inverse',inverse@matrix,np.eye(10))
    check('actual_balancing',env['balanced'],R@matrix@np.linalg.inv(D))
    check('actual_rhs_and_unscale',env['coefficients'],exact)
    check('actual_scaled_residual',env['scaled_residual'],env['balanced']@(D@env['coefficients'])-R@rhs)
    reject('omitted_row_scaling',env['balanced'],matrix@np.linalg.inv(D))
    reject('omitted_unknown_scaling',env['balanced'],R@matrix)
    reject('omitted_unscaling',exact,D@exact)

    delta = np.zeros_like(exact);delta[0]=[.125,.25j];delta[3]=[.0625j,-.125]
    env['coefficients'] = exact+delta
    execute('residual',env)
    check('nonzero_residual_identity',inverse@env['residual'],delta)
    check('nonzero_scaled_residual_identity',env['scaled_residual'],R@env['residual'])
    reject('raw_residual_as_scaled',env['scaled_residual'],env['residual'])
    wrong = copy.deepcopy(operations['residual'])
    changed = 0
    for node in ast.walk(wrong[0]):
        if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Sub): node.op=ast.Add();changed+=1
    assert changed == 1
    wrong_env=dict(env);execute('residual',wrong_env,wrong)
    reject('original_residual_sign_mutation',env['residual'],wrong_env['residual'])
    for column in range(2):
        bounded('residual_inverse_bound_column_'+str(column), np.max(abs(delta[:,column])),
                .625*np.max(abs(env['residual'][:,column])))
    reject('small_residual_is_error',1.,1/1024)

    derivative={0:np.array([[1,-1],[1,1]],complex)}
    channels={}
    for index,end in enumerate(['LEFT','RIGHT']):
        V=np.eye(5,dtype=complex);V[0,0]=2;V[0,1]=.5j
        current=np.eye(6,dtype=complex);current[0,1]=.5j;current[1,0]=-.5j;current[5,5]=-2-index
        incoming=np.zeros((5,1),complex);incoming[0,0]=1+index
        channels[end]={'right':V,'incomingValues':incoming,'incoming':[{'MATRIX_COLUMN':5}],
                       'outgoing':[{'kind':'open' if j<2 else 'evanescent','MATRIX_COLUMN':j} for j in range(5)],
                       'current':current}
    env.update(size=size,derivative=derivative,channels=channels)
    execute('observe_current',env)
    observed=env['scattering'];J=env['flux']
    C=[];offset=[]
    for index,end in enumerate(['LEFT','RIGHT']):
        T=np.kron(np.eye(5),derivative[0][index:index+1])
        Vinv=np.linalg.inv(channels[end]['right'])
        t=np.zeros((5,2),complex);t[:,index:index+1]=channels[end]['incomingValues']
        C.append((Vinv@T)[:2]);offset.append((-Vinv@t)[:2])
    C=np.vstack(C);offset=np.vstack(offset)
    check('affine_modal_observation',observed,C@env['coefficients']+offset)
    reject('missing_incoming_subtraction',observed,C@env['coefficients'])
    exact_env={**env,'coefficients':exact};execute('observe_current',exact_env)
    a=exact_env['scattering'];e=observed-a
    check('fixed_affine_offset_cancels',e,C@delta)
    check('native_full_current',env['outgoing_flux'],np.array([np.vdot(observed[:,i],J@observed[:,i]).real for i in range(2)]))
    reject('diagonal_only_current',env['outgoing_flux'],np.diag(observed.conj().T@np.diag(np.diag(J))@observed).real)
    mass=lambda z:float(np.sum(abs(z)))
    Cbound=float(np.sum(abs(C)))  # sum of component dual-l1 norms for coefficient sup norm
    beta=float(np.max(abs(J)))
    for column in range(2):
        bounded('affine_amplitude_bound_column_'+str(column),mass(e[:,column]),Cbound*np.max(abs(delta[:,column])))
        dq=abs(env['outgoing_flux'][column]-exact_env['outgoing_flux'][column])
        budget=beta*(2*mass(a[:,column])*mass(e[:,column])+mass(e[:,column])**2)
        bounded('full_current_bound_column_'+str(column),dq,budget)
        denominator=env['incoming_flux'][column]
        bounded('fixed_denominator_bound_column_'+str(column),dq/denominator,budget/denominator)
    check('cross_and_quadratic_scalar',3**2-2**2,5)
    reject('omitted_cross_scalar',5,1**2)
    reject('omitted_quadratic_scalar',5,2*2*1)
    reject('affine_offset_omitted_in_current',4**2-1**2,3**2)
    check('current_change_scalar',2*3**2-1*3**2,9)
    reject('omitted_current_change',9,0)
    bounded('perturbed_fraction_signed',abs(3/(-2)-2/(-4)),abs(3-2)/2+abs(2)*abs(-2+4)/4)
    reject('omitted_denominator_change',abs(2/2-2/4),0)
    reject('omitted_denominator_scale',abs(1/.25-0/.25),1)
    check('positive_denominator_margin',abs(2)-abs(1-2),1)
    check('zero_margin_does_not_exclude_zero',abs(0-1),abs(1))
    assert all(sha(BASE/p)==h for p,h in before.items())
    report={'status':'PASS','source_sha256':before,'instrument_sha256':sha(Path(__file__)),
            'selected_ast':ast_records,'runtime':{'python':platform.python_version(),'numpy':np.__version__},
            'checks':checks,'scope':'Synthetic actual balancing, least-squares/unscaling, residual, affine modal extraction and full current operations; conditional sensitivity witnesses.',
            'norms':'Coefficient sup norm, amplitude l1 mass, current maximum entry, observation sum of component dual-l1 row norms.',
            'limits':['No production module imports, full constructor, saved operands or physical solve.',
                      'NumPy tiny synthetic solves use a declared 1e-12 scaled absolute tolerance, not kernel-certified floating point.',
                      'No physical inverse bound or error certificate follows from measured condition numbers.',
                      'Original complete source hashes and selected AST ranges are recorded; unselected boundary/physical construction is a supplied premise.']}
    (M/'S11_lean_sensitivity_source_checks.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
