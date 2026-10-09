#!/usr/bin/env python3
"""S9b single-site copy mutations from repair build directive item 14.

Run this harness through s11c_guarded_run.py --pool s9b --memory-gib 16.
Workers are serial children of that same cgroup. No worker deadline or retry.
Baseline is the unchanged canonical engine. Mutants write only in their own
scratch replica. Every baseline/mutant tag, including metadata, is differenced.
--baseline-output permits reuse of a completed unchanged baseline, provided its
source digest sidecar matches. It never silently replays a completed worker.
"""
import argparse
import ast
import hashlib
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True
import sympy as sp
from sympy.core.symbol import Str
from sympy.core.function import AppliedUndef
from sympy.printing.repr import ReprPrinter
from sympy.functions.elementary.piecewise import ExprCondPair

HERE = Path(__file__).resolve()
ENGINE = HERE.with_name('S9b_light_bending_sympy_audit.py')
ROOT = HERE.parent.parent
REPOSITORY = ROOT.parent.parent

# Exactly one replacement in exactly one construction per copy. The mutation
# manifest contains no expected payload or bite criterion.
KNIVES = {
    'K1': ('kinetic = (omega-U*kr)**2', 'kinetic = omega**2'),
    'K2': ('inverse_metric = sp.diag(1/A, 1/r**2)', 'inverse_metric = sp.diag(1, 1/r**2)'),
    'K3': ('local_speed_squared = chi', 'local_speed_squared = c0**2'),
    'K4a': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 0'),
    'K4b': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 1'),
    'K4c': ('slope = slope  # K4c', 'slope = sp.S.Zero  # K4c'),
    'K5': ('return_orientation = -1', 'return_orientation = 1'),
    'K6': ('log_source = roundtrip_time', 'log_source = forward_time'),
    'K7': ('fixed_impact = None', "fixed_impact = sp.Symbol('s9b_fixed_b',positive=True)"),
    'K8': ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)"),
    'K9a': ('fixed_ratio_response = sp.sqrt(cs_ratio)', 'fixed_ratio_response = sp.S.One'),
    'K9b': ('power_response = responses[2]', 'power_response = responses[2].subs(f_symbol,0)'),
    'K10': ('observable_gate = ray_domain', 'observable_gate = sp.Not(ray_domain)'),
    'K11': ('azimuthal_velocity = sp.zeros(3, 1)',
            "azimuthal_velocity = sp.Symbol('s9b_azimuthal',real=True)*sp.Matrix([-xyz[1],xyz[0],0])/radius**2"),
}

RELATIONALS = {name: (lambda a,b, f=f: f(a,b,evaluate=False)) for name,f in
    (('Equality',sp.Eq),('Unequality',sp.Ne),('StrictGreaterThan',sp.Gt),
     ('StrictLessThan',sp.Lt),('GreaterThan',sp.Ge),('LessThan',sp.Le))}

FORWARD_CONDITION_TAGS = frozenset(
    'PY_S9B_FORWARD_NO_FAR_ZONE_LOSS_'+stage+'_'+quantity+'_CONDITION'
    for stage in ('DELTA_XI_ZERO','DELTA_XI_LIVE',
                  'DELTA_XI_LIVE_C_CONSTANT','DELTA_XI_LIVE_C_FIXED_RATIO',
                  'DELTA_XI_LIVE_C_POWER')
    for quantity in ('DEFLECTION','RADAR'))

# Amendment 2: finite exact samples, not a universal path-independence claim.
# These dimensionless coordinate/unit choices are evaluation inputs only.
# In particular both copies receive identical nonzero profile amplitudes.
COMPACT_PARAMETERS = {'d':sp.Rational(1,10), 'v':sp.Rational(1,8),
                      'w':sp.Rational(1,6), 'p':sp.Integer(2),
                      'q':sp.Integer(1), 'h':sp.Integer(1),
                      'L':sp.Integer(1), 'c0':sp.Integer(1)}
COMPACT_AZIMUTHAL = sp.Rational(1,12)
COMPACT_POINTS = ((3,4,0),(4,0,3),(0,6,8))
COMPACT_TARGETS = {'PY_S9B_NONRECIPROCITY_ONEFORM':'oneform',
                   'PY_S9B_NONRECIPROCITY_EXTERIOR_DERIVATIVE':'curl',
                   'PY_S9B_NONRECIPROCITY_PATH_DEPENDENCE':None}


def assigned_names(statement):
    return {node.id for target in statement.targets for node in ast.walk(target)
            if isinstance(node,ast.Name) and isinstance(node.ctx,ast.Store)}


def compact_copy(path):
    """Execute an assignment dependency slice, verbatim, from THIS source.

    No expression is replaced by a handwritten physical formula. Configuration
    is injected at symbol declarations; coordinates remain symbolic through
    the engine's differentiation. Only the three target emissions are reached.
    """
    source = path.read_text()
    module = ast.parse(source,filename=str(path))
    build = next(node for node in module.body if isinstance(node,ast.FunctionDef)
                 and node.name == 'build')
    stop = next(i for i,node in enumerate(build.body)
                if isinstance(node,ast.Assign) and 'curl' in assigned_names(node))
    wanted, selected = {'oneform','curl'}, []
    for node in reversed(build.body[:stop+1]):
        if not isinstance(node,ast.Assign):
            continue
        names = assigned_names(node)
        if wanted & names:
            selected.append(node)
            wanted = (wanted-names)|{a.id for a in ast.walk(node)
                                      if isinstance(a,ast.Name) and isinstance(a.ctx,ast.Load)}
    selected.reverse()
    # Loading the module defines its own Jet and construction helpers, but
    # does not call build(), import the fold's rows, or publish an export.
    namespace = {'__file__':str(path),'__name__':'s9b_compact_copy'}
    exec(compile(module,str(path),'exec'),namespace)
    parameter_map = {}
    for node in selected:
        exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),namespace)
        for name in assigned_names(node):
            if name in COMPACT_PARAMETERS:
                parameter_map[namespace[name]] = COMPACT_PARAMETERS[name]
                namespace[name] = COMPACT_PARAMETERS[name]
        if 'oneform' in assigned_names(node):
            # The engine obtains its radial comparator by replacing the live
            # azimuthal symbol with zero. Keep that symbol through this source
            # construction; specialize only its completed one-form, before
            # differentiating with respect to the still-symbolic coordinates.
            azimuthal_map = {symbol:COMPACT_AZIMUTHAL
                             for symbol in namespace['azimuthal_velocity'].free_symbols
                             if symbol.name == 's9b_azimuthal'}
            parameter_map.update(azimuthal_map)
            namespace['oneform'] = [component.subs(azimuthal_map)
                                    for component in namespace['oneform']]
    # Extract the actual predicate emission rather than choosing its polarity.
    predicate_emit = next(node for node in build.body if isinstance(node,ast.Expr)
        and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name)
        and node.value.func.id == 'emit' and node.value.args
        and isinstance(node.value.args[0],ast.Constant)
        and node.value.args[0].value == 'NONRECIPROCITY_PATH_DEPENDENCE')
    reached = {}
    for point in COMPACT_POINTS:
        point_map = dict(zip(namespace['xyz'],map(sp.Integer,point)))
        sampled = {tag:namespace['cas'](namespace[name]).subs(point_map)
                   for tag,name in COMPACT_TARGETS.items() if name is not None}
        sampled = {tag:sp.simplify(value) for tag,value in sampled.items()}
        local = dict(namespace, curl=sampled['PY_S9B_NONRECIPROCITY_EXTERIOR_DERIVATIVE'])
        predicate = eval(compile(ast.Expression(predicate_emit.value.args[1]),str(path),'eval'),local)
        sampled['PY_S9B_NONRECIPROCITY_PATH_DEPENDENCE'] = predicate
        for tag,value in sampled.items():
            reached.setdefault(tag,[]).append(value)
    lines = tuple((node.lineno,node.end_lineno) for node in selected)
    return ({tag:sp.Tuple(*values) for tag,values in reached.items()},parameter_map,
            namespace['xyz'],namespace['GRADES'],lines)


def compact_reduce(value, parameters, coordinates):
    """The identical parameter/coordinate reduction of a full baseline tag."""
    return sp.Tuple(*(sp.simplify(value.subs(parameters).subs(
        dict(zip(coordinates,map(sp.Integer,point))))) for point in COMPACT_POINTS))


def compact_difference(a,b):
    # Arithmetic subtraction, including the 0/1 indicators of exact booleans.
    # An unresolved Boolean is not silently converted into a truth value.
    if a in (sp.true,sp.false) and b in (sp.true,sp.false):
        return sp.Integer(int(bool(b))-int(bool(a)))
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
        return sp.Tuple(*(compact_difference(x,y) for x,y in zip(a,b)))
    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
        return (b-a).applyfunc(sp.factor)
    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
        return sp.factor(b-a)
    raise ValueError(('compact subtraction undefined',type(a),type(b)))


def compact_k11(baseline, source, mutant, folder):
    unmutated = mutant.with_name('S9b_compact_unmutated.py')
    unmutated.write_text(source)
    original,parameters,coordinates,grades,lines = compact_copy(unmutated)
    corrupted,other_parameters,other_coordinates,other_grades,other_lines = compact_copy(mutant)
    coverage = sp.Tuple(
        sp.Tuple(Str('method'),Str('verbatim AST assignment slice; exact parameter specialization; derivatives before coordinate sampling')),
        sp.Tuple(Str('profile_parameters'),sp.Tuple(*(sp.Tuple(Str(k),v) for k,v in COMPACT_PARAMETERS.items()))),
        sp.Tuple(Str('azimuthal_parameter'),COMPACT_AZIMUTHAL),
        sp.Tuple(Str('sample_points'),sp.Tuple(*(sp.Tuple(*point) for point in COMPACT_POINTS))),
        sp.Tuple(Str('radial_retained_grades'),sp.Tuple(*(sp.Tuple(*grade) for grade in grades))),
        sp.Tuple(Str('additional_truncation'),Str('NONE; nonradial correction is the exact engine expression')),
        sp.Tuple(Str('arithmetic'),Str('EXACT; no seed, rounding, or numerical error')),
        sp.Tuple(Str('boolean_difference'),Str('corrupted minus unmutated truth indicators; pointwise predicate only')),
        sp.Tuple(Str('extracted_line_ranges'),sp.Tuple(*(sp.Tuple(*pair) for pair in lines))))
    print('PY_LOCAL_S9B_ABLATION_K11_CONFIGURATION: '+render(coverage),flush=True)
    if coordinates != other_coordinates or grades != other_grades or lines != other_lines:
        raise ValueError('compact extraction/configuration mismatch')
    if any(other_parameters.get(key)!=value for key,value in parameters.items()):
        raise ValueError('compact parameter mismatch')
    reached,not_evaluated = [],[]
    for tag,full in sorted(baseline.items()):
        if tag in original and tag in corrupted:
            reduced = compact_reduce(full,parameters,coordinates)
            delta = compact_difference(original[tag],corrupted[tag])
            residual = compact_difference(reduced,original[tag])
            payload = sp.Tuple(full,original[tag],corrupted[tag],
                sp.Tuple(Str('DIFFERENCE'),delta,coverage),
                sp.Tuple(Str('METHOD_RESIDUAL'),residual,coverage),
                sp.Tuple(Str('REDUCED_FULL_BASELINE'),reduced))
            reached.append(tag)
        else:
            missing = Str('NOT_EVALUATED')
            payload = sp.Tuple(full,missing,missing,
                sp.Tuple(Str('DIFFERENCE'),missing,coverage),
                sp.Tuple(Str('METHOD_RESIDUAL'),missing,coverage))
            not_evaluated.append(tag)
        print('PY_LOCAL_S9B_ABLATION_K11_'+tag+': '+render(payload),flush=True)
    # Receipt is run state, not an engine output or an acceptance verdict.
    import json
    (folder/'compact-receipt.json').write_text(json.dumps({
        'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'mutant_sha256':hashlib.sha256(mutant.read_bytes()).hexdigest(),
        'harness_sha256':hashlib.sha256(HERE.read_bytes()).hexdigest(),
        'evaluated_tags':reached,'not_evaluated_tags':not_evaluated,
        'configuration':render(coverage)},indent=2)+'\n')


class PayloadPrinter(ReprPrinter):
    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef):
            return 'Function(%r, **%r)' % (expr.__name__, dict(sorted(expr._kwargs.items())))
        return super()._print_FunctionClass(expr)


def render(value):
    return PayloadPrinter().doprint(value)


def parse(path):
    namespace = {'__builtins__': {}, **vars(sp), 'Str':Str,
                 'ExprCondPair':ExprCondPair, **RELATIONALS}
    result = {}
    for line in path.read_text().splitlines():
        tag, sep, payload = line.partition(': ')
        if not sep or not (tag.startswith('PY_S9B_') or tag.startswith('PY_LOCAL_S9B_')):
            raise ValueError(('untagged line',line[:100]))
        if tag in result:
            raise ValueError(('duplicate tag',tag))
        result[tag] = eval(payload,namespace)
    if not result:
        raise ValueError(('empty stream',str(path)))
    missing = FORWARD_CONDITION_TAGS-set(result)
    if missing:
        raise ValueError(('missing amended forward conditions', sorted(missing)))
    return result


def difference(a,b):
    if a == b:
        return sp.S.Zero
    if isinstance(a,sp.Tuple) and isinstance(b,sp.Tuple) and len(a)==len(b):
        return sp.Tuple(*(difference(x,y) for x,y in zip(a,b)))
    if isinstance(a,sp.MatrixBase) and isinstance(b,sp.MatrixBase) and a.shape==b.shape:
        return b-a
    if isinstance(a,sp.Expr) and isinstance(b,sp.Expr):
        return b-a
    if isinstance(a,sp.logic.boolalg.Boolean) and isinstance(b,sp.logic.boolalg.Boolean):
        return sp.Xor(a,b)
    return sp.sympify(a!=b)


def run(script, folder):
    folder.mkdir(parents=True,exist_ok=False)
    output = folder/'stdout'
    with output.open('w') as out, (folder/'stderr').open('w') as err:
        completed = subprocess.run([sys.executable,'-B',str(script)],stdout=out,stderr=err,
                                   stdin=subprocess.DEVNULL,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    (folder/'exit-code').write_text(str(completed.returncode)+'\n')
    if completed.returncode:
        raise RuntimeError(('worker stopped',str(folder),completed.returncode))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch',type=Path,required=True)
    parser.add_argument('--baseline-output',type=Path)
    parser.add_argument('--knife',choices=tuple(KNIVES),action='append')
    args = parser.parse_args()
    scratch = args.scratch.resolve()
    scratch.relative_to(REPOSITORY/'_scratch')
    scratch.mkdir(parents=True,exist_ok=False)
    source = ENGINE.read_text()
    digest = hashlib.sha256(ENGINE.read_bytes()).hexdigest()
    if args.baseline_output:
        output = args.baseline_output.resolve()
        if output.with_name(output.name+'.source-sha256').read_text().strip()!=digest:
            raise ValueError('baseline source pin mismatch')
    else:
        output = run(ENGINE,scratch/'baseline')
        output.with_name(output.name+'.source-sha256').write_text(digest+'\n')
    baseline = parse(output)
    for knife in args.knife or KNIVES:
        before,after = KNIVES[knife]
        count = source.count(before)
        print('PY_LOCAL_S9B_ABLATION_'+knife+'_CONSTRUCTION: '+render(sp.Tuple(Str(before),Str(after),sp.Integer(count))),flush=True)
        if count != 1:
            raise ValueError(('construction-site count',knife,count))
        replica = scratch/knife/'replica'
        scripts = replica/'scripts'
        scripts.mkdir(parents=True)
        (replica/'directives').symlink_to(ROOT/'directives',target_is_directory=True)
        for name in ('ledger_fold.py','S11c_b_exports.py','S11c_c1_exports.py','S11c_c2_exports.py'):
            (scripts/name).symlink_to(ROOT/'scripts'/name)
        mutant = scripts/ENGINE.name
        mutant.write_text(source.replace(before,after,1))
        if knife == 'K11':
            compact_k11(baseline,source,mutant,scratch/knife)
            continue
        corrupted = parse(run(mutant,scratch/knife/'run'))
        for tag in sorted(set(baseline)|set(corrupted)):
            a = baseline.get(tag,Str('MISSING_TAG'))
            b = corrupted.get(tag,Str('MISSING_TAG'))
            print('PY_LOCAL_S9B_ABLATION_'+knife+'_'+tag+': '+render(sp.Tuple(a,b,difference(a,b))),flush=True)
        if set(baseline)!=set(corrupted):
            raise ValueError(('tag inventory changed',knife))


if __name__=='__main__':
    main()
