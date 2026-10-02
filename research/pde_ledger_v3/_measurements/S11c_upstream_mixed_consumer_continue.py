#!/usr/bin/env python3
"""Candidate selected direct-term source and slab-consumer diagnostic.

Scientific imports occur only after pooled containment and immutable input pins.
This prints operands, identities and controls. Independent review and a separate
record must decide their meaning; successful execution is not physics acceptance.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import re
import sys
import textwrap
import time
import traceback
from types import SimpleNamespace

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(value, message):
    if value is not True:
        raise ValueError(message)


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def extract_function(source, name):
    tree = ast.parse(source)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
    require(len(nodes) == 1, 'unique native function ' + name)
    node = nodes[0]
    return ast.get_source_segment(source, node)


def literal_record(source, key):
    found = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Dict):
            continue
        for k, v in zip(node.keys, node.values):
            if isinstance(k, ast.Constant) and k.value == key:
                fields = {ast.literal_eval(a): b for a, b in zip(v.keys, v.values)}
                call = fields['value']
                require(isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                        and call.func.id == '_restore' and len(call.args) == 1,
                        'literal restore source ' + key)
                found.append(ast.literal_eval(call.args[0]))
    require(len(found) == 1, 'unique literal record ' + key)
    return found[0]


def containment():
    memory = 4 * 1024**3
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines()
                 if x.startswith('0::'))
    base = Path('/sys/fs/cgroup') / group.lstrip('/')
    actual = {k: (base/k).read_text().strip() for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    actual.update(cgroup=str(base), affinity=sorted(os.sched_getaffinity(0)),
                  threads={k: os.environ.get(k) for k in THREADS})
    require(actual['memory.max'] == str(memory) and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == '32' and len(actual['affinity']) == 1
            and all(v == '1' for v in actual['threads'].values()), 'pooled resource limits')
    require('S11C_POOLED_GUARD_MANIFEST' in os.environ, 'pooled guard required')
    require(resource.getrlimit(resource.RLIMIT_CPU) == (resource.RLIM_INFINITY,) * 2,
            'no inherited CPU deadline')
    resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    return actual | dict(nativeAddressSpace=memory, durationLimits=None)


class Journal:
    def __init__(self, out):
        self.out = out
        self.completed = []
        self.active = None

    def encode(self, x):
        if isinstance(x, (sp.Basic, sp.MatrixBase)):
            return dict(text=str(x), srepr=sp.srepr(x))
        if isinstance(x, dict):
            require(all(isinstance(k, str) for k in x), 'literal string journal keys')
            return {k: self.encode(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [self.encode(v) for v in x]
        return x

    def emit(self, name, value):
        self.active = name
        path = self.out/(name + '.json')
        save(path, self.encode(value))
        return dict(path=path.name, sha256=sha(path), bytes=path.stat().st_size)

    def stage(self, name, inputs, function):
        self.active = name
        operand = self.emit(name + '-input', inputs)
        result = function()
        receipt = self.emit(name + '-return', result)
        self.completed.append(dict(name=name, input=operand, result=receipt))
        self.active = None
        return result

    def zero(self, name, left, right):
        residual = sp.factor(sp.cancel(sp.simplify(left-right)))
        self.emit(name, dict(left=left, right=right, residual=residual))
        require(residual == 0, name)



def assignment(source, function, target, occurrence=0):
    """Verbatim native assignment; call sites declare any diagnostic replacement."""
    body = ast.parse(extract_function(source, function)).body[0].body
    nodes = [n for n in body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == target for t in n.targets)]
    require(len(nodes) > occurrence, 'native assignment ' + target)
    return ast.get_source_segment(extract_function(source, function), nodes[occurrence])


def selected_constructor(source, key, case, outer=None):
    raw = literal_record(source, key)
    node = ast.parse(raw, mode='eval').body
    def label(x):
        require(isinstance(x, ast.Call) and len(x.args) == 1, 'literal case label')
        return ast.literal_eval(x.args[0])
    if outer:
        node = next(x.args[1] for x in node.args if label(x.args[0]) == outer)
    matches = [x.args[1] for x in node.args if [label(k) for k in x.args[0].args] == list(case)]
    require(len(matches) == 1, 'unique selected native case')
    value = next(x.args[1] for x in matches[0].args if label(x.args[0]) == 'VALUE')
    return ast.get_source_segment(raw, value)



def selected_increment(row, eta, sigma, D, slots, ref_factor, jet_factor):
    """The same selected-slot substitution for the physical and ablated rows."""
    pplus,pminus,jplus,jminus=slots
    increment=row.subs({pplus:eta*sigma*D*ref_factor,jplus:eta*sigma*D*jet_factor,
                        pminus:0,jminus:0},simultaneous=True)
    mixed=sp.cancel(sp.diff(increment,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
    return increment,mixed


def exact_nonzero(value):
    """Same exact zero predicate; retain unknown and never use a numeric tolerance."""
    direct = value.is_zero
    evidence = dict(value=value, directZero=direct, finite=value.is_finite,
                    decision=None, route='direct-symbolic-flag')
    if value.is_finite is False:
        return evidence
    if direct is not None:
        evidence['decision'] = direct is False
        return evidence
    if value.free_symbols or value.is_number is not True:
        return evidence
    expanded = sp.expand_complex(value)
    real, imaginary = (sp.cancel(x) for x in expanded.as_real_imag())
    evidence.update(route='exact-constant-real-imaginary', expanded=expanded,
                    real=real, imaginary=imaginary,
                    componentZero=[real.is_zero, imaginary.is_zero],
                    componentFinite=[real.is_finite, imaginary.is_finite])
    if real.is_finite is True and imaginary.is_finite is True:
        if real.is_zero is False or imaginary.is_zero is False:
            evidence['decision'] = True
        elif real.is_zero is True and imaginary.is_zero is True:
            evidence['decision'] = False
    return evidence


def forbidden_profile_keys(profiles):
    return sorted(k.name for k in profiles if k.name in ('eta_bg','sigma_W')
                  or k.name.startswith(('w1_profile','m1_profile','gamma_')))


def zero_grade_leftovers(value, profiles):
    names = {k.name for k in profiles} | {'eta_bg','sigma_W'}
    return sorted(a.name for a in value.free_symbols if a.name in names
                  or a.name.startswith(('w1_profile','m1_profile')))


def load_saved_json(path):
    """Decode saved constructor pairs only inside the already-contained worker."""
    raw = json.loads(Path(path).read_text())
    def decode(value):
        if isinstance(value, dict):
            if set(value) == {'text','srepr'}:
                return sp.sympify(value['srepr'], locals={'Str':Str})
            return {k:decode(v) for k,v in value.items()}
        if isinstance(value,list):
            return [decode(v) for v in value]
        return value
    return raw,decode(raw)


def nested_fragment(source, function, name):
    outer=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name==function)
    child=next(n for n in outer.body if isinstance(n,ast.FunctionDef) and n.name==name)
    return textwrap.dedent(ast.get_source_segment(source,child))


def resumed_operator_form(expr, J, serial):
    from sympy.core.function import AppliedUndef
    name='resumed-operator-%03d' % serial
    J.emit(name+'-input',dict(expression=expr,algorithm='Original independent-jet reconstruction; exact cancel of its residual before literal zero decision.'))
    jets=sorted(expr.atoms(sp.Derivative,AppliedUndef),key=sp.default_sort_key)
    dummies={j:sp.Dummy('independentJet'+str(i)) for i,j in enumerate(jets)}
    polynomial=expr.xreplace(dummies)
    coefficients={str(j):sp.cancel(sp.diff(polynomial,d)) for j,d in dummies.items()}
    raw_remainder=sp.expand(polynomial-sum(coefficients[str(j)]*d for j,d in dummies.items()))
    dependent_coefficients={name:sorted(v.free_symbols & set(dummies.values()),key=sp.default_sort_key)
                            for name,v in coefficients.items()}
    J.emit(name+'-raw-reconstruction',dict(expression=expr,jets=jets,
        substitutions=[[j,d] for j,d in dummies.items()],polynomial=polynomial,
        coefficients=coefficients,rawRemainder=raw_remainder,
        coefficientJetDependencies=dependent_coefficients,
        originalLiteralZero=raw_remainder==0))
    remainder=sp.cancel(sp.together(raw_remainder))
    J.emit(name+'-exact-residual',dict(rawRemainder=raw_remainder,
        cancelledRemainder=remainder,literalZero=remainder==0,
        noNumericalTolerance=True))
    require(not any(dependent_coefficients.values()),'coefficients independent of trial jets')
    require(remainder==0,'linear arbitrary-trial jet reconstruction after exact cancellation')
    vals=list(coefficients.values())
    decisions={name:exact_nonzero(v) for name,v in coefficients.items()}
    state='zero' if all(v==0 for v in vals) else (
        'nonzero' if any(v['decision'] is True for v in decisions.values()) else 'unresolved')
    result=dict(expression=expr,independentJetCoefficients=coefficients,
        coefficientNonzeroEvidence=decisions,rawRemainder=raw_remainder,
        remainder=remainder,state=state)
    J.emit(name+'-return',result)
    return result


def persist_control_list(J, name):
    class SavedList(list):
        def append(self, value):
            J.emit(name+'-%03d-return' % len(self), value)
            super().append(value)
    return SavedList()


def resumed_tail(base_source):
    tree=ast.parse(base_source);function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='science')
    first=next(n for n in function.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='downstream' for t in n.targets))
    tail='\n'.join(base_source.splitlines()[first.lineno-1:function.end_lineno])
    tail=tail.replace('downstream=[]',"downstream=persist_control_list(J,'source-consumer-control')",1)
    tail=tail.replace('consumer_controls=[]',"consumer_controls=persist_control_list(J,'pressure-consumer-control')",1)
    tail=tail.replace("            selected_amplitude=", "            J.emit('source-consumer-control-%03d-input' % len(downstream),dict(row=name,channel=channel,factor=factor,sourcePlane=source_plane,sourceChange=source_change,epsilon=epsilon))\n            selected_amplitude=",1)
    tail=tail.replace("            damaged_source=flats[1]-piece", "            J.emit('source-consumer-control-%03d-input' % len(downstream),dict(row=name,channel='divergence-omission',factor=factor,source=flats[1],jet=jet,piece=piece,movement=movement,epsilon=epsilon))\n            damaged_source=flats[1]-piece",1)
    tail=tail.replace("        raw=row_source[name]", "        J.emit('pressure-consumer-control-%03d-input' % len(consumer_controls),dict(row=name,raw=row_source[name],pplus=pplus,slots=slots,probeSource=probe_source,epsilon=epsilon))\n        raw=row_source[name]",1)
    return tail


def science(manifest,J):
    import shutil
    prior=Path(manifest['priorDirectory'])
    base_source=Path(manifest['originalWorker']).read_text()
    copies={}
    copy_root=J.out/'prior-copy';copy_root.mkdir()
    require(sorted(p.name for p in prior.iterdir() if p.is_file())==sorted(manifest['priorFiles']), 'exact prior result census')
    for name,pin in manifest['priorFiles'].items():
        require(sha(prior/name)==pin,'prior result pin '+name)
        shutil.copyfile(prior/name,copy_root/name)
        require(sha(copy_root/name)==pin,'prior copy '+name)
        copies[name]=dict(path='prior-copy/'+name,sha256=pin,bytes=(copy_root/name).stat().st_size)
    J.emit('prior-copy-index',copies)
    # Reuse every complete published observation by exact byte copy. Decode only
    # operands needed for unfinished controls. No old scientific function is called.
    raw,restored={},{}
    for name in manifest['restoreFiles']:
        raw[name],restored[name]=load_saved_json(copy_root/name)
    J.emit('restore-index',dict(files=[dict(name=n,sha256=copies[n]['sha256']) for n in restored],
        originalScienceCalled=False,completeSourceGradeRestrictionConsumerWorkReplayed=False,
        oldTopLevelCompletedStages=0,sourceOfOperands='Published JSON receipts, not a successful top-level return.'))
    end=restored['end-to-end-controls-input.json'];binding=restored['binding-context.json']
    controls=restored['source-controls.json'];fourier=restored['selected-fourier-contraction.json']
    plus=restored['plus-source-restrictions.json'];weak=restored['weak-directions.json']
    expected={
        'source00':raw['end-to-end-controls-input.json']['upperSource00']==raw['plus-source-restrictions.json']['source00'],
        'sourcePlane':raw['end-to-end-controls-input.json']['sourcePlane']==raw['selected-fourier-contraction.json']['sourceAmplitude'],
        'rowFactors':raw['end-to-end-controls-input.json']['rowFactors']==raw['selected-fourier-contraction.json']['rowFactorsPerCombinedSource'],
        'chemical':raw['end-to-end-controls-input.json']['chemicalChannel']==raw['source-controls.json']['removedChemicalChannel'],
        'velocityPiece':raw['end-to-end-controls-input.json']['velocityPiece']==raw['source-controls.json']['removedVelocityPiece'],
        'sourcePins':raw['source-consumer-input.json']['sourcePins']==json.loads(Path(manifest['originalManifest']).read_text())['sourcePins']}
    expected['scalarRows']=all(raw['end-to-end-controls-input.json']['scalarConsumers'][n]==raw[n+'-consumer-input.json']['raw'] for n in ('THETA_BALANCE','E_W_BALANCE'))
    expected['divergenceChoices']=all(a['jet']==b['jet'] and a['piece']==b['piece'] and a['restrictedPiece']==b['curlSourceMovement'] for a,b in zip(raw['end-to-end-controls-input.json']['divergenceChoices'],raw['source-controls.json']['divergenceOmissions'])) and len(end['divergenceChoices'])==len(controls['divergenceOmissions'])
    J.emit('prior-argument-joins',expected)
    require(all(expected.values()),'actual saved argument joins')
    eta,sigma=end['grades'];epsilon=end['epsilon']
    row_source={name:restored[name+'-consumer-input.json']['raw'] for name in ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')}
    all_objects=[*row_source.values(),end['upperSource00'],end['sourcePlane'],binding['density'],*end['rowFactors'].values()]
    symbols=set().union(*(v.atoms(sp.Symbol) for v in all_objects))
    def atom(name):
        matches=[v for v in symbols if v.name==name]
        require(len(matches)==1,'unique saved symbol '+name)
        return matches[0]
    profiles={eq.lhs:eq.rhs for eq in binding['density'][0].atoms(sp.Equality)
              if isinstance(eq.lhs,sp.Symbol) and eq.lhs!=sigma}
    require({str(k):v for k,v in profiles.items()}==binding['profileEqualities'],'restored complete profile equality map')
    density_map={atom('rho_br_bg_rho4_constant'):binding['density'][1]}
    require({str(k):v for k,v in density_map.items()}==binding['densityMap'],'restored density argument')
    numeric={a:binding['numeric'][a.name] for a in symbols if a.name in binding['numeric']}
    pplus=atom('delta_p_plus');pminus=atom('delta_p_minus');jplus=atom('d_w_delta_p_plus');jminus=atom('d_w_delta_p_minus')
    slots=(pplus,pminus,jplus,jminus)
    require(pplus==end['upperPressureSlot'],'saved native pressure-slot symbol')
    # D and H are formal labels. No kernel/integral is reconstructed.
    old_mixed=restored['THETA_BALANCE-consumer-input.json']['mixedPerCombinedSource']
    D=next(a for a in old_mixed.free_symbols if a.name=='inherited_whole_bare_mixed_kernel')
    amplitudes={name:next((a for a in end['sourcePlane'].free_symbols if a.name=='amplitude_'+name),sp.Symbol('amplitude_'+name)) for name in ('u_1','u_2','u_3','theta','e_W')}
    X=tuple(sp.Symbol('s11cc2X'+str(i),real=True) for i in (1,2,3));time_symbol=sp.Symbol('s11cc2Time',real=True)
    c2=Path(manifest['nativeTrialSource']).read_text()
    ns=dict(sp=sp,re=re,X=X,TIME=time_symbol,NEW_DIMENSIONS={})
    for name in ('field','wave_jet'):
        exec(compile(extract_function(c2,name),'<native-'+name+'>','exec'),ns)
    env=dict(globals(),J=J,profiles=profiles,density_map=density_map,numeric=numeric,
        eta=eta,sigma=sigma,epsilon=epsilon,row_source=row_source,row_factors=end['rowFactors'],
        flats={1:end['upperSource00']},source_plane=end['sourcePlane'],
        velocity_movement=controls['velocityMovement'],theta_movement=controls['thetaMovement'],
        choices=[(a['jet'],a['piece'],a['restrictedPiece']) for a in end['divergenceChoices']],
        restrictions={1:plus['restrictions']},scalar_rows=('THETA_BALANCE','E_W_BALANCE'),
        amplitudes=amplitudes,pplus=pplus,pminus=pminus,jplus=jplus,jminus=jminus,
        slots=slots,D=D,ref_factor=end['referenceFactor'],jet_factor=end['jetFactor'],
        H=sp.Function('selectedKernelAppliedSource')(*X,time_symbol),X=X,c2=c2,
        vector=weak['vectorPerSourceImage'],curl=weak['reverseCurl'],
        records=json.loads((copy_root/'native-source-census-join.json').read_text())['saved'],ns=ns,
        pattern=re.compile(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)'))
    for name in ('bind','wave_atoms','restrict'):
        fragment=nested_fragment(base_source,'science',name)
        exec(compile(fragment,'<unchanged-'+name+'>','exec'),env)
    serial=[0]
    def operator_form(expr):
        index=serial[0];serial[0]+=1
        return resumed_operator_form(expr,J,index)
    env['operator_form']=operator_form
    J.emit('resumed-context',dict(grades=[eta,sigma],epsilon=epsilon,slots=slots,
        rowFactors=end['rowFactors'],source00=end['upperSource00'],sourcePlane=end['sourcePlane'],
        referenceFactor=end['referenceFactor'],jetFactor=end['jetFactor'],D=D,
        numeric={str(k):v for k,v in numeric.items()},profiles={str(k):v for k,v in profiles.items()},
        originalNumericKeysUnusedHere=sorted(set(binding['numeric'])-{a.name for a in numeric}),
        reconstructedUnsavedContext=['formal H','missing zero-channel amplitude symbol','native trial helper callables','coordinate symbols','selected bind mappings from saved density/equalities/numerics'],
        firstUnfinishedSection='downstream source-omission controls; prior in-memory unsaved prefix is recomputed from published operands, not labelled restored'))
    # This tail is source-extracted verbatim except per-control persistence and
    # the exact reconstruction predicate above. Never call original science().
    tail=resumed_tail(base_source)
    J.emit('resumed-source-tail',dict(source=tail,originalWorker=manifest['originalWorker'],
        originalSha256=sha(manifest['originalWorker']),equationsChanged=False))
    exec(compile('def resume_selected_controls():\n'+tail,'<saved-consumer-unfinished-tail>','exec'),env)
    return env['resume_selected_controls']()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads(args.inputs.read_text())
    for path,pin in manifest['sourcePins'].items():
        require(sha(path)==pin,'source pin '+path)
    require(manifest['sourcePins'][str(Path(__file__).resolve())]==sha(__file__),'worker pin')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c')
    args.out.mkdir(exist_ok=False)
    started=time.monotonic(); J=None; result={}; code=1
    try:
        save(args.out/'containment.json',containment())
        global sp,Str
        import sympy as sp
        from sympy.core.symbol import Str
        J=Journal(args.out)
        result=J.stage('source-consumer-continue',dict(sourcePins=manifest['sourcePins'],
            priorDirectory=manifest['priorDirectory'],restoreFiles=manifest['restoreFiles']),
            lambda:science(manifest,J))
        result['executionStatus']='COMPLETED_CANDIDATE_OPERANDS'
        code=0
    except BaseException:
        result=dict(executionStatus='FAILED_PRESERVED',traceback=traceback.format_exc(),
                    incompleteOperation=None if J is None else J.active,automaticRetry=False)
        save(args.out/'failure.json',result)
    finally:
        post={p:dict(expected=h,actual=sha(p)) for p,h in manifest['sourcePins'].items()}
        save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):
            result['integrityFailure']=True; code=1
        save(args.out/'operation-index.json',[] if J is None else J.completed)
        save(args.out/'artifact-index.json',[
            dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size)
            for p in sorted(args.out.iterdir()) if p.is_file()])
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result)
        sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':
    sys.exit(main())
