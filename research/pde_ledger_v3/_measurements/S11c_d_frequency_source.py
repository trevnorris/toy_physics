#!/usr/bin/env python3
"""Keep the actual reduced frequency dependence and end threshold operands."""
import argparse,contextlib,copy,json,resource,shutil,signal,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import sympy as sp
import S11c_d_profile_form as p

f=p.f;engine=p.engine;boundary=p.boundary;grades=p.grades;interior=p.interior
PLAN=f.M/'S11c_d_frequency_source_plan.md';PREFIX='FREQUENCY_SOURCE_LAB_HELD_RHO4_CONSTANT'


def load(base):
    values=p.load(base);r,coefficients,ends,system,reference,baseline,packet,old,adapter,_,pins,operands=values
    uniform,ucp,up=f.accepted_packet(f.M/'S11c_d_uniform_source_checkpoint.json','uniform-source.pickle')
    for n,h in ucp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n)==h,('unchanged uniform source',n))
        if n in pins:f.require(pins[n]==h,'common physical source')
        pins[n]=h
    operands[str(up)]=f.digest(up)
    uniform_cp=json.loads((f.M/'S11c_d_uniform_response_checkpoint.json').read_text());uniform_path=Path(uniform_cp['runDirectory']);accepted_ends={}
    f.require(uniform_cp['status']=='PUBLISHED_ANNEX_VERIFIED','accepted actual uniform modes/current source')
    for label in ('REFERENCE','LEFT','RIGHT'):
        path=uniform_path/(label.lower()+'-symbolic.pickle');f.require(f.digest(path)==uniform_cp['artifacts'][path.name]['sha256'],'accepted end pencil packet');accepted_ends[label]=f.unpickle(path);operands[str(path)]=f.digest(path)
    control=f.M/'S11c_d_first_jet_checkpoint.json';cp=json.loads(control.read_text());f.require(cp['status']=='PUBLISHED_ANNEX_VERIFIED','preceding controls complete')
    f.require(f.digest(f.ROOT/cp['publication']['path'])==cp['publication']['sha256'],'actual accepted control output')
    for q in (Path(__file__),PLAN,control,f.M/'S11c_d_uniform_source_checkpoint.json',f.M/'S11c_d_uniform_response_checkpoint.json',f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'):
        pins[str(q.resolve().relative_to(f.ROOT))]=f.digest(q)
    for n in pins:
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest)
    frequency=sp.Symbol('s11cdProfilePoleFrequency',complex=True);engine.PHYSICAL_METADATA.dimensions.known[frequency]=(0,-1,0)
    live=engine.NumericalReducedAction(SimpleNamespace(r=r),{},adapter.input.specification)
    live.input.parameters=dict(live.input.parameters,omega=frequency)
    live.input.origin={v:v for v in adapter.input.origin}
    control_spec=copy.deepcopy(adapter.input.specification);control_spec['parameters']['omega']='11/10'
    control=engine.NumericalReducedAction(SimpleNamespace(r=r),{},control_spec)
    # This uses the original positive-real input only as the reference seed;
    # it does not extend ChannelInput's physical-chart declaration silently.
    manifest={'sourceFiles':pins,'inputPackets':operands,'settings':system['settings'],'input':adapter.input.specification,
      'frequencySymbol':str(frequency),'frequencyUnit':[0,-1,0],'referenceFrequency':str(adapter.input.parameters['omega']),
      'frequencyBinding':'Separate formal continuation variable; original input and positive-real seed unchanged.',
      'arithmeticControlFrequency':'11/10; independent rebinding only, not a new physical scattering instance',
      'scope':'Reduced source dependence and end-threshold candidates only; no profile pole solve or global analytic chart.'}
    f.save(base/'inputs.json',manifest)
    return dict(r=r,coefficients=coefficients,ends=ends,system=system,reference=reference,baseline=baseline,packet=packet,old=old,adapter=adapter,live=live,control=control,frequency=frequency,uniform=uniform,acceptedEnds=accepted_ends,pins=pins,operands=operands,manifest=manifest)


def census(expression,frequency):
    nonanalytic=tuple(sorted({v for v in sp.preorder_traversal(expression)
        if isinstance(v,sp.Basic) and v.has(frequency) and (v.func in (sp.Abs,sp.sign,sp.re,sp.im,sp.conjugate,sp.Heaviside,sp.DiracDelta) or isinstance(v,sp.Piecewise))},key=sp.default_sort_key))
    denominators=tuple(sorted({v.base for v in expression.atoms(sp.Pow) if v.exp.is_negative and v.base.has(frequency)},key=sp.default_sort_key))
    fractional=tuple(sorted({v for v in expression.atoms(sp.Pow) if v.base.has(frequency) and v.exp.is_Rational and v.exp.q!=1},key=sp.default_sort_key))
    return {'frequencyDependent':expression.has(frequency),'nonanalyticNodes':nonanalytic,'literalDenominatorBases':denominators,'fractionalPowers':fractional,
            'integralLimits':tuple(sorted({tuple(v.limits) for v in expression.atoms(sp.Integral)},key=sp.default_sort_key))}


def sources(base,data):
    r=data['r'];w=data['frequency'];origin=data['adapter'].input.origin;w0=data['adapter'].input.parameters['omega'];bound=data['old']['domain-binding.pickle']['bound'];cuts=bound['cutoffBindings']
    native=p.rebind(base/'baseline-binding.pickle',r,bound,data['packet'],data['adapter'])
    replay=p.baseline_check(native,data['old']['source-binding.pickle']);size=len(data['system']['nodes'])
    local=np.zeros_like(data['system']['localMatrix'])
    for n,m in native['local'].items():
        for i in range(5):
            for j in range(5):local[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(m[i,j],r,data['system'])[:,None]*data['system']['derivativeMatrices'][n]
    local_residual=local-data['system']['localMatrix'];f.require(boundary.norm(local_residual)==0,'actual accepted baseline local matrix')
    records={};inventory={};directory=base/'records';directory.mkdir();derivative_counts={};addresses={};frequency_controls=0
    for key,item in data['packet']['records'].items():
        original=item['record']['ORIGINAL'];live=engine.memo_xreplace(data['live'].bind(original),cuts)
        expected=engine.memo_xreplace(data['adapter'].bind(original),cuts);actual=live.subs({w:w0,**origin},simultaneous=True);residual=actual-expected
        control_frequency=data['control'].input.parameters['omega'];control_expected=engine.memo_xreplace(data['control'].bind(original),cuts);control_actual=live.subs({w:control_frequency,**origin},simultaneous=True)
        control_residual=control_actual-control_expected;frozen_difference=control_actual-expected;frequency_controls+=int(frozen_difference!=0)
        address=tuple(item['address']);first=sp.diff(live,w);second=sp.diff(first,w)
        record={'address':address,'original':original,'liveFrequencyAndGrades':live,'retainedFrequency':live.subs(origin),
          'firstFrequencyDerivative':first,'secondFrequencyDerivative':second,'boundAtReference':actual,'acceptedBinding':expected,'bindingResidual':residual,
          'controlFrequency':control_frequency,'controlExpected':control_expected,'controlActual':control_actual,'controlResidual':control_residual,'frozenFrequencyDifference':frozen_difference,
          'unit':item['record']['UNIT'],'census':census(live,w),'firstDerivativeCensus':census(first,w),'secondDerivativeCensus':census(second,w)}
        path=directory/(key+'.pickle');f.atomic_pickle(path,record);inventory[key]={'path':str(path.relative_to(base)),'sha256':f.digest(path),'bytes':path.stat().st_size};f.save(base/'record-inventory.json',inventory)
        f.require(residual==0 and control_residual==0,('exact reference and independent frequency source binding',key))
        old_limits=census(expected,w)['integralLimits'];new_limits=tuple(sorted({tuple(sp.Tuple(*(v.subs({w:w0,**origin},simultaneous=True) if isinstance(v,sp.Basic) else v for v in limit)) for limit in limits) for limits in record['census']['integralLimits']},key=sp.default_sort_key))
        f.require(new_limits==old_limits,'complete original ordered integration limits')
        records[key]=record;addresses[address]=record;derivative_counts[address[0]]=derivative_counts.get(address[0],0)+int(first!=0)
    source_frequencies={};row_census=[]
    for si in range(35):
        old=bound['sources'][0,si];value=engine.memo_xreplace(data['live'].bind(old['symbolicFrequency']),cuts)
        residual=value.subs({w:w0,**origin},simultaneous=True)-old['frequency'];f.require(residual==0,'actual source Fourier character seed')
        source_frequencies[si]={'frequency':value,'original':old['symbolicFrequency'],'accepted':old['frequency'],'residual':residual,'dependsOnFrequency':value.has(w)}
    for row in bound['rows']:
        factors=[addresses['factor',row['index'],i] for i in range(len(row['factors']))];used=sorted({v['sourceIndex'] for v in row['factors']})
        dependent=any(v['census']['frequencyDependent'] for v in factors) or any(addresses['source',i]['census']['frequencyDependent'] or source_frequencies[i]['dependsOnFrequency'] for i in used)
        row_census.append({'index':row['index'],'sources':used,'frequencyDependent':dependent,'limits':row['limits'],'sourceLimit':row['sourceLimit'],
          'sourceColumns':tuple(sorted({native['jets'][i]['column'] for i in used}))})
        f.require(len(row_census[-1]['sourceColumns'])==1,'full frequency row input-field join')
    f.require(frequency_controls>0,'actual frequency-freezing control responds')
    result={'records':records,'respondingFrequencyControls':frequency_controls,'sourceFrequencies':source_frequencies,'rowCensus':row_census,'baselineChecks':replay,'localBindingResidual':local_residual,'derivativeCounts':derivative_counts,'frequency':w,'referenceFrequency':w0,'origin':origin}
    f.atomic_pickle(base/'frequency-sources.pickle',result);return result


def end_sources(base,data):
    r=data['r'];frequency=data['frequency'];inp=data['adapter'].input;uniform=data['uniform'];ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r=r;ends.kn=sp.Symbol('s11cdSpectralNormalMomentum',real=True)
    modes=engine.FullPencilModes(ends,uniform['curl'],uniform['units']['weak']);results={}
    coordinate=sp.Symbol('s11cdFrequencyCoefficientCoordinate',real=True);radical_coordinate=sp.Symbol('s11cdRadicalCoefficientCoordinate',complex=True)
    engine.PHYSICAL_METADATA.dimensions.known[coordinate]=(0,0,0);engine.PHYSICAL_METADATA.dimensions.known[radical_coordinate]=(0,0,0)
    for label in ('REFERENCE','LEFT','RIGHT'):
        algebraic,relation,branch=modes.analytic(uniform['records'][label]['strong'])
        mapping=inp.mapping(algebraic,relation,(r.omega,modes.k,modes.q,*inp.origin))
        origin={s:sp.S.Zero if label=='REFERENCE' else v for s,v in inp.origin.items()}
        live=algebraic.xreplace(mapping).subs(origin).xreplace({r.omega:frequency});wave=relation.xreplace(mapping).subs(origin).xreplace({r.omega:frequency})
        old=data['acceptedEnds'][label];matrix_residual=(live.subs(frequency,inp.parameters['omega'])-old['freshPencil']).applyfunc(sp.cancel)
        wave_residual=sp.cancel(wave.subs(frequency,inp.parameters['omega'])-old['curve'])
        transport=-sp.diff(wave,frequency)/sp.diff(wave,modes.q);derivative=live.diff(frequency)+live.diff(modes.q)*transport
        tangency=sp.cancel(wave.diff(frequency)+wave.diff(modes.q)*transport)
        packet={'livePencil':live,'wave':wave,'originalAlgebraic':algebraic,'originalRelation':relation,'mapping':mapping,'origin':origin,'branchResiduals':branch,
          'referencePencilResidual':matrix_residual,'referenceWaveResidual':wave_residual,'radicalFrequencyTransport':transport,'frequencyDerivative':derivative,'waveTangencyResidual':tangency,
          'frequency':frequency,'momentum':modes.k,'radical':modes.q,'unitConvention':'Entry units inherit physical row/field ratios; no continuum inverse expansion is used.'}
        path=base/(label.lower()+'-frequency-pencil.pickle');f.atomic_pickle(path,packet)
        f.require(all(v==0 for v in branch) and all(v==0 for v in matrix_residual) and wave_residual==tangency==0,'actual end branch/frequency seed/tangency joins')
        f.require(not(live.free_symbols|wave.free_symbols)-{frequency,modes.k,modes.q},'complete fixed material/tangential binding')
        # Scalar elimination is a numerical coefficient-frame diagnostic. The
        # physical frequency and radical unit maps are retained separately.
        zero=live.subs(modes.k,0).xreplace({frequency:coordinate,modes.q:radical_coordinate});curve=wave.subs(modes.k,0).xreplace({frequency:coordinate,modes.q:radical_coordinate})
        (numerator,denominator),cleared,row_denominators=modes.rational_determinant(zero)
        elimination=sp.resultant(numerator,curve,radical_coordinate);poly=sp.Poly(elimination,coordinate)
        f.require(not poly.is_zero,'nonzero actual zero-normal-momentum elimination')
        normalized=poly.monic();f.require(all(v.is_real is True for v in normalized.all_coeffs()),'real coefficient-frame threshold polynomial')
        normalized=sp.Poly(normalized.as_expr(),coordinate,domain=sp.QQ)
        squarefree=normalized.sqf_part();intervals=squarefree.intervals(eps=sp.Rational(1,10**12));factorization=sp.factor_list(normalized.as_expr(),coordinate)
        reconstruction=sp.expand(factorization[0]*sp.prod(v**n for v,n in factorization[1])-normalized.as_expr())
        f.require(reconstruction==0,'full diagnostic factor reconstruction')
        branches=sp.solve(curve.subs(radical_coordinate,0),coordinate)
        threshold={'coefficientMatrix':zero,'clearedMatrix':cleared,'rowDenominators':row_denominators,'numerator':numerator,'denominator':denominator,'wave':curve,
          'elimination':elimination,'normalizedElimination':normalized.as_expr(),'squarefreeElimination':squarefree.as_expr(),'factorization':factorization,'factorReconstructionResidual':reconstruction,
          'realRootIntervals':intervals,'bulkBranchFrequencies':branches,'frequencyCoordinate':coordinate,'radicalCoordinate':radical_coordinate,
          'coordinateMap':{'frequencyUnit':(0,-1,0),'radicalUnit':(0,-1,0),'normalMomentumUnit':(-1,0,0),'normalMomentumCoefficient':0,'unitFrame':inp.frame,'fieldReferenceUnits':data['ends']['fieldUnits'],'equationReferenceUnits':data['ends']['rowUnits'],'matrixConvention':'Numerical coefficient matrix between equation and field reference-unit bases; entries and elimination coordinates are dimensionless.'},
          'scope':'Real roots are zero-normal-momentum algebraic candidates; retain opposite radical sheets and denominator artifacts. Not a profile-frequency pole set or completed physical threshold classification.'}
        f.atomic_pickle(base/(label.lower()+'-threshold-candidates.pickle'),threshold);results[label]={'pencil':packet,'thresholds':threshold}
        f.save(base/'end-inventory.json',{name:{'pencilSha256':f.digest(base/(name.lower()+'-frequency-pencil.pickle')),'thresholdSha256':f.digest(base/(name.lower()+'-threshold-candidates.pickle'))} for name in results})
    return results


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W'];zero=(0,0,0)
    def put(name,value,unit=zero):
        body=engine.cas(value);engine.emit(PREFIX+'_'+name,engine.carrier_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,unit if callable(unit) else lambda path:unit))
    for key,v in result['sources']['records'].items():
        put(key+'_SOURCE',sp.Tuple(v['original'],v['liveFrequencyAndGrades']),v['unit'])
        for order,name in ((1,'firstFrequencyDerivative'),(2,'secondFrequencyDerivative')):put(key+'_DERIVATIVE_'+str(order),v[name],tuple(a+(order if i==1 else 0) for i,a in enumerate(v['unit'])))
        put(key+'_SEED',sp.Tuple(v['boundAtReference'],v['acceptedBinding'],v['bindingResidual']),v['unit'])
        put(key+'_FREQUENCY_CONTROL',sp.Tuple(v['controlActual'],v['controlExpected'],v['controlResidual'],v['frozenFrequencyDifference']),v['unit'])
        # Domains are operands, not a single generic certificate. Each node's
        # original frequency-live expression is retained in the source packet.
        boundary.structural_flags(PREFIX+'_'+key+'_DOMAIN',{'frequencyDependent':v['census']['frequencyDependent'],'nonanalyticNodes':tuple(map(str,v['census']['nonanalyticNodes'])),'literalDenominatorBases':tuple(map(str,v['census']['literalDenominatorBases'])),'fractionalPowers':tuple(map(str,v['census']['fractionalPowers'])),'scope':'Literal expression census; branches, zeros and global operator domains remain to be resolved.'})
    put('FREQUENCY',result['sources']['frequency'],(0,-1,0));put('REFERENCE_FREQUENCY',result['sources']['referenceFrequency'],(0,-1,0))
    for label,v in result['ends'].items():
        end=v['pencil'];fields=result['fieldUnits'];rows=result['rowUnits']
        def unit(path,extra=0):
            i,j=divmod(path[0],5);return tuple(a-b+(extra if axis==1 else 0) for axis,(a,b) in enumerate(zip(rows[i],fields[j])))
        for name in ('livePencil','referencePencilResidual'):put(label+'_'+name,end[name],unit)
        put(label+'_frequencyDerivative',end['frequencyDerivative'],lambda path:unit(path,1))
        put(label+'_WAVE',sp.Tuple(end['wave'],end['referenceWaveResidual']),(0,-2,0));put(label+'_WAVE_TANGENCY',end['waveTangencyResidual'],(0,-1,0));put(label+'_RADICAL_TRANSPORT',end['radicalFrequencyTransport'],zero)
        threshold=v['thresholds']
        for name in ('coefficientMatrix','clearedMatrix','rowDenominators','numerator','denominator','wave','elimination','normalizedElimination','squarefreeElimination','factorization','factorReconstructionResidual','realRootIntervals','bulkBranchFrequencies'):
            put(label+'_COEFFICIENT_FRAME_'+name,threshold[name],zero)
        boundary.structural_flags(PREFIX+'_'+label+'_DOMAIN',{'coordinateMap':threshold['coordinateMap'],'scope':threshold['scope'],'frequencyCoordinate':str(threshold['frequencyCoordinate']),'radicalCoordinate':str(threshold['radicalCoordinate'])})
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'scope':result['scope'],'rowFrequencyCensus':result['sources']['rowCensus'],
      'retainedOperatorStatus':'Frequency-source grades remain live. End pencils use the full retained contrast at the approved input; no pole-inverse re-expansion. Numerical seed and elimination data are coefficients in the declared unit frame.',
      'fixedSpacePlan':'Finite Chebyshev coefficient and equation spaces with inherited field/row units. Analytic outgoing endpoint identification and profile-frequency inverse remain unconstructed.'})


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('frequency-source budget; preserve completed records and end packets')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);data=load(base);src=sources(base,data);ends=end_sources(base,data)
    result={'sources':src,'ends':ends,'fieldUnits':data['ends']['fieldUnits'],'rowUnits':data['ends']['rowUnits'],'sourceFiles':data['pins'],'inputPackets':data['operands'],
      'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),'scope':'Frequency-live reduced sources and algebraic end-threshold candidates; no profile-pole solve, analytic outgoing boundary chart, full inverse or physical-bound classification.'}
    f.atomic_pickle(base/'frequency-source.pickle',result);before=f.digest(base/'frequency-source.pickle');engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,data['r']);keys={tag:'s11cdFrequencySource'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique frequency-source tags');entries[tag]=grades._restore(body)
    old_emit=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('frequency-source full output replay',tag));seen.add(tag)
    engine.emit=replay
    try:emit_result(result,data['r']);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=old_emit
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'complete frequency-source keys/payloads')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES','_DOMAIN'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)};f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'frequency source units');f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'frequency source grades/lambda');paths+=1 if structural else len(fields['PATHS'])
    last='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[last]},list(entries)[:list(entries).index(last)])
    f.require(before==f.digest(base/'frequency-source.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'frequency packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'frequency sources/inputs post hashes')
    checks={'runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'records':len(src['records']),'rows':len(src['rowCensus']),'sources':len(src['sourceFrequencies']),
      'respondingFrequencyControls':src['respondingFrequencyControls'],'frequencyDependentRows':[v['index'] for v in src['rowCensus'] if v['frequencyDependent']],'derivativeCounts':src['derivativeCounts'],
      'nonanalyticRecordCount':sum(bool(v['census']['nonanalyticNodes']) for v in src['records'].values()),'endThresholdDegrees':{k:int(sp.degree(v['thresholds']['squarefreeElimination'],v['thresholds']['frequencyCoordinate'])) for k,v in ends.items()},
      'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'frequency-source.pickle'),
      'artifacts':{str(q.relative_to(base)):{'bytes':q.stat().st_size,'sha256':f.digest(q)} for q in base.rglob('*') if q.suffix in ('.pickle','.out') and 'source' not in q.relative_to(base).parts},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
