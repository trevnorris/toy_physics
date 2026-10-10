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
    'K1': ('kinetic = (omega-advected_velocity.dot(momenta))**2', 'kinetic = omega**2'),
    'K2': ('inverse_metric = sp.diag(1/A,1,1)', 'inverse_metric = sp.eye(3)'),
    'K3': ('local_speed_squared = chi', 'local_speed_squared = c0**2'),
    'K4a': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 0'),
    'K4b': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 1'),
    'K4c': ('RADIAL_FREEZE = None', 'RADIAL_FREEZE = 2'),
    'K5': ('return_orientation = -1', 'return_orientation = 1'),
    'K6': ('log_source = roundtrip_time', 'log_source = forward_time'),
    'K7': ('fixed_impact = None', "fixed_impact = sp.Symbol('s9b_fixed_b',positive=True)"),
    'K8': ('mass_density = rho', "mass_density = sp.Symbol('s9b_constant_density',positive=True)"),
    'K9a': ('fixed_ratio_response = sp.sqrt(cs_ratio)', 'fixed_ratio_response = sp.S.One'),
    'K9b': ('power_response = (1+f_symbol)**s', 'power_response = sp.S.One'),
    'K10': ('observable_gate = real_domain', 'observable_gate = sp.Not(real_domain)'),
    'K11': ('azimuthal_component = sp.S.Zero', "azimuthal_component = sp.Symbol('s9b_azimuthal',real=True)"),
}

RELATIONALS = {name: (lambda a,b, f=f: f(a,b,evaluate=False)) for name,f in
    (('Equality',sp.Eq),('Unequality',sp.Ne),('StrictGreaterThan',sp.Gt),
     ('StrictLessThan',sp.Lt),('GreaterThan',sp.Ge),('LessThan',sp.Le))}

FORWARD_CONDITION_TAGS = frozenset(
    ['PY_S9B_F_'+stage+'_B_CONDITION_'+quantity
     for stage in ('FLOW','LIVE') for quantity in ('DEFLECTION','RADAR')]+
    ['PY_S9B_F_LIVE_C_'+response+'_CONDITION_'+quantity
     for response in ('CONSTANT','FIXED_RATIO','POWER') for quantity in ('DEFLECTION','RADAR')])


def progress(object_name):
    print('S9B_HARNESS_WORK '+str(object_name), file=sys.stderr, flush=True)

# Amendment 2: finite exact samples, not a universal path-independence claim.
# These dimensionless coordinate/unit choices are evaluation inputs only.
# In particular both copies receive identical nonzero profile amplitudes.
COMPACT_PARAMETERS = {'d':sp.Rational(1,10), 'v':sp.Rational(1,8),
                      'w':sp.Rational(1,6), 'p':sp.Integer(2),
                      'q':sp.Integer(1), 'h':sp.Integer(1),
                      'L':sp.Integer(1), 'c0':sp.Integer(1)}
COMPACT_AZIMUTHAL = sp.Rational(1,12)
COMPACT_POINTS = ((3,4,0),(4,0,3),(0,6,8))
COMPACT_TARGETS = {'PY_LOCAL_S9B_NONRECIPROCITY_ONEFORM':'oneform',
                   'PY_LOCAL_S9B_NONRECIPROCITY_EXTERIOR_DERIVATIVE':'curl',
                   'PY_LOCAL_S9B_NONRECIPROCITY_PATH_PREDICATE':'predicate'}


def compact_copy(path):
    """Call the evaluated copy's own vector-dispersion/one-form construction."""
    source=path.read_text()
    namespace={'__file__':str(path),'__name__':'s9b_compact_copy'}
    exec(compile(source,str(path),'exec'),namespace)
    r=sp.Symbol('s9b_r',positive=True)
    c0=sp.Symbol('c_0',positive=True)
    config=COMPACT_PARAMETERS
    templates={namespace['Delta']:config['d']*(config['L']/r)**config['p'],
               namespace['Velocity']:config['c0']*config['v']*(config['L']/r)**config['q'],
               namespace['Embedding']:config['w']*sp.log(r/config['L'])}
    # The xi template differentiates to the configured h=1 slope. It is a
    # compact-evaluation input only, never the engine's shared-profile ansatz.
    replacements={function:(lambda x,value=value:value.subs(r,x)) for function,value in templates.items()}
    result=namespace['nonreciprocity'](r,config['c0'],templates[namespace['Delta']],
        templates[namespace['Velocity']],templates[namespace['Embedding']])
    parameters={c0:config['c0'],sp.Symbol('s9b_azimuthal',real=True):COMPACT_AZIMUTHAL}
    reached={tag:sp.Tuple(*(sp.simplify(namespace['cas'](result[name]).subs(parameters).subs(
            dict(zip(result['xyz'],map(sp.Integer,point))))) for point in COMPACT_POINTS))
             for tag,name in COMPACT_TARGETS.items()}
    functions=[n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)
               and n.name in ('vector_dispersion','nonreciprocity')]
    lines=tuple((node.lineno,node.end_lineno) for node in functions)
    reducer=(namespace['profile_substitute'],replacements,parameters)
    return reached,reducer,result['xyz'],result['grades'],lines


def compact_reduce(value, reducer, coordinates):
    apply,replacements,parameters=reducer
    reduced=apply(value,replacements).subs(parameters)
    return sp.Tuple(*(sp.simplify(reduced.subs(dict(zip(coordinates,map(sp.Integer,point)))))
                      for point in COMPACT_POINTS))


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
    progress('K11 compact unmutated dispersion/nonreciprocity')
    original,parameters,coordinates,grades,lines = compact_copy(unmutated)
    progress('K11 compact mutated dispersion/nonreciprocity')
    corrupted,other_parameters,other_coordinates,other_grades,other_lines = compact_copy(mutant)
    coverage = sp.Tuple(
        sp.Tuple(Str('method'),Str('engine vector_dispersion and nonreciprocity calls; exact profile specialization; derivatives before coordinate sampling')),
        sp.Tuple(Str('profile_parameters'),sp.Tuple(*(sp.Tuple(Str(k),v) for k,v in COMPACT_PARAMETERS.items()))),
        sp.Tuple(Str('azimuthal_parameter'),COMPACT_AZIMUTHAL),
        sp.Tuple(Str('sample_points'),sp.Tuple(*(sp.Tuple(*point) for point in COMPACT_POINTS))),
        sp.Tuple(Str('radial_retained_grades'),sp.Tuple(*(sp.Tuple(*grade) for grade in grades))),
        sp.Tuple(Str('additional_truncation'),Str('NONE beyond the engine retained grades; all vector velocity components enter the same graded dispersion')),
        sp.Tuple(Str('arithmetic'),Str('EXACT; no seed, rounding, or numerical error')),
        sp.Tuple(Str('boolean_difference'),Str('corrupted minus unmutated truth indicators; pointwise predicate only')),
        sp.Tuple(Str('extracted_line_ranges'),sp.Tuple(*(sp.Tuple(*pair) for pair in lines))))
    print('PY_LOCAL_S9B_ABLATION_K11_CONFIGURATION: '+render(coverage),flush=True)
    if coordinates != other_coordinates or grades != other_grades or lines != other_lines:
        raise ValueError('compact extraction/configuration mismatch')
    reached,not_evaluated = [],[]
    for tag,full in sorted(baseline.items()):
        progress('K11 '+tag+': compact evaluation, difference and method residual')
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
        progress(str(folder)+': engine start')
        child = subprocess.Popen([sys.executable,'-B',str(script)],stdout=out,stderr=subprocess.PIPE,
                                 text=True,stdin=subprocess.DEVNULL,
                                 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        # Relay operational progress immediately, keeping stdout purely the
        # payload grammar and preserving a separate worker stderr transcript.
        for line in child.stderr:
            err.write(line)
            err.flush()
            print(str(folder)+': '+line.rstrip('\n'),file=sys.stderr,flush=True)
        code = child.wait()
        child.stderr.close()
    (folder/'exit-code').write_text(str(code)+'\n')
    progress(str(folder)+': engine exit '+str(code))
    if code:
        raise RuntimeError(('worker stopped',str(folder),code))
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
    progress('baseline: parse '+str(output))
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
        worker_output = run(mutant,scratch/knife/'run')
        progress(knife+': parse '+str(worker_output))
        corrupted = parse(worker_output)
        for tag in sorted(set(baseline)|set(corrupted)):
            progress(knife+' '+tag+': baseline, corrupted and difference')
            a = baseline.get(tag,Str('MISSING_TAG'))
            b = corrupted.get(tag,Str('MISSING_TAG'))
            print('PY_LOCAL_S9B_ABLATION_'+knife+'_'+tag+': '+render(sp.Tuple(a,b,difference(a,b))),flush=True)
        if set(baseline)!=set(corrupted):
            raise ValueError(('tag inventory changed',knife))


if __name__=='__main__':
    main()
