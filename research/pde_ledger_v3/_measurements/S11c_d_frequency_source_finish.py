#!/usr/bin/env python3
"""Finish actual complex threshold diagnostics from completed source operands."""
import argparse,ast,copy,hashlib,inspect,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_frequency_source as q

f=q.f;engine=q.engine;boundary=q.boundary
PLAN=f.M/'S11c_d_frequency_source_threshold_recovery_plan.md'
REPAIR=f.M/'S11c_d_frequency_source_threshold_repair.json'
PREVIOUS=f.STORE/'s11c-frequency-source-20260919/recovery-01/complete'
DIAGNOSTIC=f.STORE/'s11c-frequency-source-20260919/threshold-repair/elimination.pickle'


def real_axis_analysis(poly,coordinate):
    """Keep the complex polynomial and isolate its exact real common roots."""
    normalized=sp.Poly(poly.monic().as_expr(),coordinate,domain=sp.QQ_I)
    expression=normalized.as_expr();real=sp.Poly(sp.re(expression).expand(),coordinate,domain=sp.QQ);imag=sp.Poly(sp.im(expression).expand(),coordinate,domain=sp.QQ)
    decomposition=sp.expand(expression-real.as_expr()-sp.I*imag.as_expr())
    common=sp.gcd(real,imag);squarefree=common.sqf_part()
    real_quotient,real_remainder=sp.div(real,common);imag_quotient,imag_remainder=sp.div(imag,common)
    if imag.is_zero:
        bezout_real=sp.Poly(1/real.LC(),coordinate,domain=sp.QQ);bezout_imag=sp.Poly(0,coordinate,domain=sp.QQ);gcd=real.monic()
    else:
        bezout_real,bezout_imag,gcd=sp.gcdex(real,imag)
    bezout_residual=sp.expand(bezout_real.as_expr()*real.as_expr()+bezout_imag.as_expr()*imag.as_expr()-common.as_expr())
    intervals=squarefree.intervals(eps=sp.Rational(1,10**12))
    factorization=sp.factor_list(expression,coordinate,extension=sp.I)
    reconstruction=sp.expand(factorization[0]*sp.prod(v**n for v,n in factorization[1])-expression)
    result={'normalized':normalized,'complexSquarefree':normalized.sqf_part(),'realPart':real.as_expr(),'imaginaryPart':imag.as_expr(),
      'decompositionResidual':decomposition,'realCommonFactor':common.as_expr(),'realSquarefree':squarefree.as_expr(),
      'realQuotient':real_quotient.as_expr(),'imaginaryQuotient':imag_quotient.as_expr(),'divisionResiduals':(real_remainder.as_expr(),imag_remainder.as_expr()),
      'bezoutReal':bezout_real.as_expr(),'bezoutImaginary':bezout_imag.as_expr(),'bezoutResidual':bezout_residual,
      'intervals':intervals,'factorization':factorization,'factorReconstructionResidual':reconstruction,
      'scope':'For real frequency coefficients, both actual real and imaginary polynomial parts must vanish. The complex polynomial, denominator and sheet artifacts are retained.'}
    f.require(decomposition==reconstruction==bezout_residual==0 and real_remainder.is_zero and imag_remainder.is_zero and gcd==common,'exact real-axis polynomial decomposition/divisibility/Bezout joins')
    f.require(sum(n for _,n in intervals)==squarefree.count_roots(-sp.oo,sp.oo),'complete real isolation for computed common polynomial')
    for (left,right),multiplicity in intervals:
        f.require(multiplicity==1 and squarefree.count_roots(left,right)==1,'individual real common-root interval')
    return result


def end_adapter():
    """Reversible saved-pencil hooks and complex-to-real candidate analysis."""
    original=ast.parse(inspect.getsource(q.end_sources));tree=copy.deepcopy(original);fn=tree.body[0]
    fn.args.args.extend([ast.arg(arg='completed'),ast.arg(arg='eliminations')])
    loop=next(v for v in fn.body if isinstance(v,ast.For))
    first_guard=next(i for i,v in enumerate(loop.body) if isinstance(v,ast.Expr) and isinstance(v.value,ast.Call) and isinstance(v.value.func,ast.Attribute) and v.value.func.attr=='require')
    construction=copy.deepcopy(loop.body[:first_guard])
    reuse=ast.parse("""
if label in completed:
    packet=completed[label]
    live=packet['livePencil'];wave=packet['wave'];branch=packet['branchResiduals']
    matrix_residual=(live.subs(frequency,inp.parameters['omega'])-data['acceptedEnds'][label]['freshPencil']).applyfunc(sp.cancel)
    wave_residual=sp.cancel(wave.subs(frequency,inp.parameters['omega'])-data['acceptedEnds'][label]['curve'])
    tangency=sp.cancel(wave.diff(frequency)+wave.diff(modes.q)*packet['radicalFrequencyTransport'])
    f.require(packet['momentum']==modes.k and packet['radical']==modes.q and packet['frequency']==frequency,'saved end coordinates')
    f.require(live==packet['originalAlgebraic'].xreplace(packet['mapping']).subs(packet['origin']).xreplace({r.omega:frequency}),'saved complete end source expression')
else:
    pass
""").body[0];reuse.orelse=construction;loop.body[:first_guard]=[reuse]
    def assigns(node,name):return isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==name for t in node.targets)
    det_index=next(i for i,v in enumerate(loop.body) if isinstance(v,ast.Assign) and isinstance(v.value,ast.Call) and isinstance(v.value.func,ast.Attribute) and v.value.func.attr=='rational_determinant')
    det_statements=copy.deepcopy(loop.body[det_index:det_index+3])
    det_reuse=ast.parse("""
if label in eliminations:
    saved=eliminations[label]
    f.require(saved['coefficientMatrix']==zero and saved['wave']==curve and saved['coordinate']==coordinate and saved['radicalCoordinate']==radical_coordinate,'actual saved elimination operand joins')
    numerator=saved['numerator'];denominator=saved['denominator'];cleared=saved['clearedMatrix'];row_denominators=saved['rowDenominators'];elimination=saved['elimination'];poly=sp.Poly(elimination,coordinate)
    for i in range(5):
        for j in range(5):f.require(sp.cancel(zero[i,j]*row_denominators[i]-cleared[i,j])==0,'saved actual row denominator clearing')
else:
    pass
""").body[0];det_reuse.orelse=det_statements;loop.body[det_index:det_index+3]=[det_reuse]
    checkpoint=ast.parse("""
f.atomic_pickle(base/(label.lower()+'-elimination.pickle'),{'coefficientMatrix':zero,'wave':curve,'numerator':numerator,'denominator':denominator,'clearedMatrix':cleared,'rowDenominators':row_denominators,'elimination':elimination,'coordinate':coordinate,'radicalCoordinate':radical_coordinate})
""").body
    loop.body[det_index+1:det_index+1]=checkpoint
    start=next(i for i,v in enumerate(loop.body) if assigns(v,'normalized'));stop=next(i for i,v in enumerate(loop.body) if assigns(v,'branches'))
    old_analysis=copy.deepcopy(loop.body[start:stop])
    analysis=ast.parse("""
analysis=real_axis_analysis(poly,coordinate)
normalized=analysis['normalized'];squarefree=analysis['complexSquarefree'];intervals=analysis['intervals'];factorization=analysis['factorization'];reconstruction=analysis['factorReconstructionResidual']
""").body;loop.body[start:stop]=analysis
    ti=next(i for i,v in enumerate(loop.body) if assigns(v,'threshold'));addition=ast.parse("threshold['realAxisAnalysis']=analysis").body[0];loop.body.insert(ti+1,addition)
    # Undo every explicit change and compare the complete function AST.
    reverse=copy.deepcopy(tree);rfn=reverse.body[0];rfn.args.args=rfn.args.args[:-2];rloop=next(v for v in rfn.body if isinstance(v,ast.For))
    rloop.body=[v for v in rloop.body if ast.dump(v)!=ast.dump(addition)]
    i=next(i for i,v in enumerate(rloop.body) if assigns(v,'analysis'));rloop.body[i:i+len(analysis)]=old_analysis
    i=next(i for i,v in enumerate(rloop.body) if isinstance(v,ast.If) and ast.dump(v.test)==ast.dump(det_reuse.test));rloop.body[i:i+1+len(checkpoint)]=det_statements
    rloop.body[:1]=construction
    f.require(ast.dump(reverse)==ast.dump(original),'whole end constructor reverse AST join')
    env=dict(q.__dict__,real_axis_analysis=real_axis_analysis);ast.fix_missing_locations(tree);exec(compile(tree,'<saved-frequency-end-adapter>','exec'),env)
    return env['end_sources']


def emit_result(result,r):
    q.emit_result(result,r)
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W']
    for label,value in result['ends'].items():
        a=value['thresholds']['realAxisAnalysis']
        for name in ('realPart','imaginaryPart','decompositionResidual','realCommonFactor','realSquarefree','realQuotient','imaginaryQuotient','divisionResiduals','bezoutReal','bezoutImaginary','bezoutResidual'):
            tag=q.PREFIX+'_'+label+'_REAL_AXIS_'+name;body=engine.cas(a[name]);engine.emit(tag,engine.carrier_fingerprint(body));engine.emit('METADATA_'+tag,modes.numeric_metadata(body,lambda path:(0,0,0)))
        boundary.structural_flags(q.PREFIX+'_'+label+'_REAL_AXIS_DOMAIN',{'scope':a['scope'],'coordinateMap':value['thresholds']['coordinateMap']})


def validation_tail():
    """Use the full original result/emission/validation tail without edits."""
    tree=ast.parse(inspect.getsource(q.main));body=tree.body[0].body
    index=next(i for i,v in enumerate(body) if isinstance(v,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='result' for t in v.targets))
    tail=ast.FunctionDef(name='finish',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=v) for v in ('base','data','src','ends','started')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body[index:],decorator_list=[])
    module=ast.Module(body=[tail],type_ignores=[]);ast.fix_missing_locations(module);env=dict(q.__dict__,emit_result=emit_result);exec(compile(module,'<original-frequency-validation-tail>','exec'),env)
    return env['finish']


def reuse(base,data):
    old=json.loads((PREVIOUS/'inputs.json').read_text());repair=json.loads(REPAIR.read_text())
    f.require(f.digest(Path(q.__file__))==repair['originalSourceSha256'],'unchanged completed source checker')
    for name,sha in old['sourceFiles'].items():
        f.require(f.digest(PREVIOUS/'source'/name)==sha==f.digest(f.ROOT/name),'current/frozen completed source join')
        if name in data['pins']:f.require(data['pins'][name]==sha,'shared consumed helper hash')
        else:data['pins'][name]=sha;dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
    for name,sha in old['inputPackets'].items():f.require(f.digest(Path(name))==sha,'original input packet unchanged')
    f.require(all(old['inputPackets'][name]==sha for name,sha in data['operands'].items()),'fresh/old consumed input join')
    data['operands'].update(old['inputPackets'])
    f.require(old['input']==data['manifest']['input'] and old['settings']==json.loads(json.dumps(data['manifest']['settings'])),'same physical input and finite settings')
    copies={}
    for name,sha in repair['originalArtifacts'].items():
        src=PREVIOUS/name;dest=base/name;dest.parent.mkdir(parents=True,exist_ok=True);f.require(f.digest(src)==sha,'completed operand at failure')
        shutil.copyfile(src,dest);f.require(f.digest(dest)==sha,'byte-identical saved operand copy');copies[name]=sha;data['operands'][str(src)]=sha
    src=f.unpickle(base/'frequency-sources.pickle');inventory=json.loads((PREVIOUS/'record-inventory.json').read_text());f.require(len(src['records'])==len(inventory)==375 and len(src['rowCensus'])==80 and len(src['sourceFrequencies'])==35,'complete finished source census')
    proof_count=0;cert_count=0
    for key,item in src['records'].items():
        path=base/inventory[key]['path'];f.require(f.digest(path)==inventory[key]['sha256'],'completed source record hash');record=f.unpickle(path)
        original=data['packet']['records'][key];f.require(record==item and record['original']==original['record']['ORIGINAL'] and record['address']==tuple(original['address']) and record['unit']==original['record']['UNIT'],'complete original record and source/field/unit join')
        for label,comparison in record['bindingComparisons'].items():
            f.require(comparison['normalizedResidual']==0 and all(v==0 for v in comparison['proofResiduals']),'completed exact binding proofs');proof_count+=len(comparison['proofResiduals'])
            op=comparison['operands'];f.require(tuple(hashlib.sha256(v.encode()).hexdigest() for v in op['representationStrings'])==op['representationSha256'],'saved live representation hashes')
            if comparison['certificate'] is not None:
                cert_count+=1;cert=comparison['certificate'];op=comparison['operands'];f.require(cert['LEFT']==op['left'] and cert['RIGHT']==op['right'] and cert['RESIDUAL']==0 and comparison['mutation']['RESIDUAL']!=0,'completed certificate actual-pair/mutation join')
            expected_keys=('boundAtReference','acceptedBinding') if label=='seed' else ('controlActual','controlExpected')
            f.require(comparison['operands']['left']==record[expected_keys[0]] and comparison['operands']['right']==record[expected_keys[1]],'all certificate-to-record source joins')
    f.require(all(v==0 for v in src['baselineChecks']['coefficientResiduals']) and boundary.norm(src['localBindingResidual'])==0,'accepted baseline source and local residuals')
    pencil=f.unpickle(base/'reference-frequency-pencil.pickle');saved=f.unpickle(DIAGNOSTIC);f.require(f.digest(DIAGNOSTIC)==repair['diagnosticEliminationSha256'] and saved['sourceSha256']==f.digest(base/'reference-frequency-pencil.pickle'),'saved actual elimination source');data['operands'][str(DIAGNOSTIC)]=f.digest(DIAGNOSTIC)
    validation={'copiedPackets':len(copies),'sourceRecords':375,'bindingPairs':750,'certificates':cert_count,'proofScalars':proof_count,'sourcePacketSha256':f.digest(base/'frequency-sources.pickle'),'copies':copies,'scope':'Reuse completed source calculations and certificates; their original live representation distinction remains explicit.'}
    f.save(base/'operand-validation.json',validation)
    for path in (Path(__file__),PLAN,REPAIR):
        name=str(path.resolve().relative_to(f.ROOT));data['pins'][name]=f.digest(path);dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
    data['manifest']['completedSourceReuse']=validation;data['manifest']['thresholdRepair']='Complex polynomial retained; real candidate isolation uses the exact real/imaginary GCD and Bezout proof.';f.save(base/'inputs.json',data['manifest'])
    return src,{'REFERENCE':pencil},{'REFERENCE':saved}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('frequency-source finish budget; retain all completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic();data=q.load(base);src,pencils,eliminations=reuse(base,data);ends=end_adapter()(base,data,pencils,eliminations);validation_tail()(base,data,src,ends,started)


if __name__=='__main__':main()
