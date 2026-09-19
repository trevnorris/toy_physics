#!/usr/bin/env python3
"""Local analytic scalar-kernel chart derived from every accepted operand."""
import argparse,contextlib,hashlib,json,resource,shutil,signal,time
from functools import lru_cache
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_frequency_source as q

f=q.f;engine=q.engine;boundary=q.boundary;grades=q.grades
PLAN=f.M/'S11c_d_frequency_chart_plan.md';PREFIX='FREQUENCY_CHART_LAB_HELD_RHO4_CONSTANT'


def load(base):
    data=q.load(base);packet,checkpoint,path=f.accepted_packet(f.M/'S11c_d_frequency_source_checkpoint.json','frequency-source.pickle')
    f.require(checkpoint['status']=='PUBLISHED_ANNEX_VERIFIED','accepted full frequency sources')
    for name,sha in packet['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==sha,'unchanged frequency-source helper')
        if name in data['pins']:f.require(data['pins'][name]==sha,'shared physical source')
        data['pins'][name]=sha
    data['operands'][str(path)]=f.digest(path)
    for path in (Path(__file__),PLAN,f.M/'S11c_d_frequency_source_checkpoint.json',f.ROOT/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        data['pins'][str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    for name in data['pins']:
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
    data['frequencyPacket']=packet;data['manifest']['scope']='Local scalar-kernel analytic continuation and end rational tables only; outgoing subspaces and frequency inverse remain unconstructed.'
    f.save(base/'inputs.json',data['manifest']);return data


def root_chart(base,data):
    src=data['frequencyPacket']['sources'];w=src['frequency'];center=src['referenceFrequency'];radius=sp.Rational(1,4);maximum=abs(center)+radius
    bases=sorted({v.base for record in src['records'].values() for v in record['census']['fractionalPowers']},key=sp.default_sort_key);roots={}
    for index,radicand in enumerate(bases):
        momenta=tuple(radicand.free_symbols-{w});f.require(len(momenta)==1 and momenta[0].is_real is True,'actual one-real-momentum radical')
        p=momenta[0];poly=sp.Poly(radicand,w,p);cw=poly.coeff_monomial(w**2);cp=poly.coeff_monomial(p**2);constant=poly.coeff_monomial(1)
        residual=sp.expand(radicand-(cw*w**2+cp*p**2+constant));f.require(residual==0 and cw>0 and cp<0 and constant<0,'actual quadratic acoustic radical')
        margin=-constant-cw*maximum**2;f.require(margin>0,'strict positive real part of continued inner radicand')
        carrier=sp.Symbol('s11cdFrequencyChartRadical'+str(index),complex=True);engine.PHYSICAL_METADATA.dimensions.known[carrier]=(0,-1,0)
        analytic=sp.I*sp.sqrt(-radicand);seed=sp.sqrt(radicand.subs(w,center));seed_residual=sp.simplify(analytic.subs(w,center)-seed)
        f.require(seed_residual==0 and sp.expand(analytic**2-radicand)==0,'actual positive-frequency seed and radical relation')
        roots[radicand]={'index':index,'momentum':p,'carrier':carrier,'radicand':radicand,'analytic':analytic,'seed':seed,'seedResidual':seed_residual,
          'relationResidual':sp.expand(analytic**2-radicand),'quadraticResidual':residual,'frequencySquareCoefficient':cw,'momentumSquareCoefficient':cp,'constant':constant,
          'innerRealPartLowerBound':margin,'modulusLowerBound':sp.sqrt(margin),'unit':(0,-1,0)}
    f.require(len(roots)==3,'all three actual Fourier radical arguments')
    result={'frequency':w,'center':center,'radius':radius,'frequencyModulusUpperBound':maximum,'roots':roots,
      'coordinateConvention':'Bounds use numerical frequency/momentum coefficients in the inherited reference-unit frame. Physical radical and frequency have unit T^-1.',
      'scope':'Real momentum coordinates and this closed frequency disk only; no outgoing-subspace or inverse assertion.'}
    f.atomic_pickle(base/'root-chart.pickle',result);return result


def lift(expression,w,roots,proofs):
    if expression==sp.sign(w):return sp.S.One
    if expression in roots:return roots[expression]['carrier']**2
    if isinstance(expression,sp.Pow) and expression.base in roots and expression.exp.is_Rational:
        f.require((2*expression.exp).is_Integer,'actual half-integer radical power');return roots[expression.base]['carrier']**(2*expression.exp)
    if isinstance(expression,sp.Piecewise):
        values=[lift(v,w,roots,proofs) for v,_ in expression.args];first=values[-1]
        for value in values:
            raw=value-first;certificate=None;normalized=raw;residuals=()
            if raw!=0:
                certificate=engine.BoundedSourceFourierAssembly.reconstruction_certificate(value,first,shared=False);normalized=certificate['RESIDUAL'];residuals=tuple(certificate['REPLAY_RESIDUALS'])+tuple(v[1] for v in certificate['PHASE_SPLITS'].values())+tuple(v[2] for v in certificate['RADICAL_POWERS'].values())
            proofs.append({'originalPiecewise':expression,'left':value,'right':first,'rawResidual':raw,'normalizedResidual':normalized,'proofResiduals':residuals,'certificate':certificate})
            f.require(normalized==0 and all(v==0 for v in residuals),'actual positive-seed Piecewise branch identity')
        return first
    if not expression.args:return expression
    return expression.func(*(lift(v,w,roots,proofs) for v in expression.args))


def denominator_chart(base,data,chart):
    src=data['frequencyPacket']['sources'];w=chart['frequency'];roots=chart['roots'];maximum=chart['frequencyModulusUpperBound']
    original=sorted({v for record in src['records'].values() for v in record['census']['literalDenominatorBases']},key=sp.default_sort_key)
    # Derive the normalized relaxation factor from the actual affine bases.
    affine=[]
    for value in original:
        if value.free_symbols<={w} and value.is_polynomial(w) and sp.degree(value,w)==1:
            polynomial=sp.Poly(value,w);affine.append(sp.cancel(value/polynomial.coeff_monomial(1)))
    f.require(affine and all(sp.cancel(v-affine[0])==0 for v in affine),'unique actual normalized relaxation factor')
    relaxation=affine[0];slope=sp.diff(relaxation,w);f.require(sp.diff(slope,w)==0 and relaxation.subs(w,0)==1,'actual relaxation normalization')
    relaxation_bound=1-abs(slope)*maximum;f.require(relaxation_bound>0,'relaxation denominator excluded on disk')
    families=[('constant',sp.S.One,sp.S.One),('relaxation',relaxation,relaxation_bound),('relaxationSquared',relaxation**2,relaxation_bound**2)];ratios={}
    for root in roots.values():
        c=root['carrier'];ratio=maximum/(relaxation_bound*root['modulusLowerBound']);f.require(ratio<1,'coupled denominator strict triangle bound')
        ratios[root['index']]=ratio;families.extend([(f"radical{root['index']}",c,root['modulusLowerBound']),(f"radicalSquared{root['index']}",c**2,root['innerRealPartLowerBound']),(f"coupled{root['index']}",1+w/(relaxation*c),1-ratio)])
    records=[]
    for index,value in enumerate(original):
        proofs=[];algebraic=lift(value,w,roots,proofs);matches=[]
        for name,family,bound in families:
            quotient=sp.cancel(algebraic/family)
            if not quotient.free_symbols and quotient!=0:
                residual=sp.cancel(algebraic-quotient*family);f.require(residual==0,'literal denominator factor identity');matches.append((name,quotient,family,bound,residual))
        f.require(matches,('unclassified actual frequency denominator',value));name,constant,family,bound,residual=matches[0]
        record={'original':value,'algebraic':algebraic,'family':name,'constant':constant,'familyExpression':family,'factorResidual':residual,'coefficientFrameModulusLowerBound':abs(constant)*bound,'branchProofs':proofs}
        f.atomic_pickle(base/f'denominator-{index}.pickle',record);records.append(record)
    result={'records':records,'relaxation':relaxation,'relaxationSlope':slope,'relaxationModulusLowerBound':relaxation_bound,'coupledRatios':ratios,'scope':'Exact classification of every original frequency-dependent denominator; local scalar bounds in the declared coefficient frame.'}
    f.atomic_pickle(base/'denominator-chart.pickle',result);return result


@lru_cache(maxsize=None)
def seed_root_identity(radicand,w,momentum,frequency):
    value=sp.expand(radicand.subs(w,frequency));poly=sp.Poly(value,momentum)
    quadratic=poly.coeff_monomial(momentum**2);constant=poly.coeff_monomial(1)
    reconstruction=sp.expand(value-quadratic*momentum**2-constant)
    f.require(momentum.is_real is True and reconstruction==0 and quadratic<0 and constant<0,'strictly negative real seed radicand for all real momenta')
    root=sp.I*sp.sqrt(-value);residual=sp.simplify(sp.sqrt(value)-root)
    f.require(residual==0 and sp.expand(root**2-value)==0,'principal negative-real seed root identity')
    return {'radicand':value,'momentum':momentum,'frequency':frequency,'quadraticCoefficient':quadratic,'maximumOnRealAxis':constant,'root':root,'principalRoot':sp.sqrt(value),'residual':residual,'reconstructionResidual':reconstruction,'scope':'This real seed frequency and real momentum only; strict negative radicand fixes the principal root sign.'}


def seed_comparison(base,actual,expected,unit,chart,frequency):
    """Join the analytic seed using its explicit negative-real root proofs."""
    base.mkdir(parents=True,exist_ok=False);forms=tuple(sp.srepr(v) for v in (actual,expected,actual-expected))
    operands={'left':actual,'right':expected,'rawResidual':actual-expected,'unit':unit,'representationStrings':forms,'representationSha256':tuple(hashlib.sha256(v.encode()).hexdigest() for v in forms),'exactEqualLive':actual==expected}
    f.atomic_pickle(base/'operands.pickle',operands);mapping={};identities=[]
    for power in expected.atoms(sp.Pow):
        if not(power.exp.is_Rational and power.exp.q==2):continue
        for root in chart['roots'].values():
            identity=seed_root_identity(root['radicand'],chart['frequency'],root['momentum'],frequency)
            scale=sp.cancel(power.base/identity['radicand'])
            if not scale.free_symbols and scale.is_positive is True:
                replacement=(sp.sqrt(scale)*identity['root'])**(2*power.exp)
                mapping[power]=replacement;identities.append({'originalPower':power,'replacement':replacement,'positiveScale':scale,'scaleResidual':sp.expand(power.base-scale*identity['radicand']),'rootIdentity':identity});break
    normalized_expected=engine.memo_xreplace(expected,mapping)
    f.atomic_pickle(base/'root-join.pickle',{'identities':identities,'expected':expected,'normalizedExpected':normalized_expected,'mapping':mapping})
    joined=q.binding_comparison(base/'normalized',actual,normalized_expected,unit)
    result=dict(joined,operands=operands,rawResidual=actual-expected,normalizedOperands=joined['operands'],rootIdentities=identities,
      proofResiduals=tuple(joined['proofResiduals'])+tuple(v for item in identities for v in (item['scaleResidual'],item['rootIdentity']['residual'],item['rootIdentity']['reconstructionResidual'])))
    f.require(all(v==0 for v in result['proofResiduals']),'all actual seed-root proof residuals')
    f.atomic_pickle(base/'comparison.pickle',result);return result


def source_chart(base,data,chart):
    src=data['frequencyPacket']['sources'];w=chart['frequency'];roots=chart['roots'];mapping={v['carrier']:v['analytic'] for v in roots.values()};records={};inventory={};directory=base/'records';directory.mkdir()
    for key,old in src['records'].items():
        proofs=[];algebraic=lift(old['liveFrequencyAndGrades'],w,roots,proofs);analytic=engine.memo_xreplace(algebraic,mapping)
        f.require(not any(isinstance(v,sp.Piecewise) or v.func in (sp.sign,sp.Abs,sp.conjugate,sp.re,sp.im) for v in sp.preorder_traversal(analytic) if isinstance(v,sp.Basic) and v.has(w)),'complete analytic constructor census')
        record={'address':old['address'],'originalLive':old['liveFrequencyAndGrades'],'algebraic':algebraic,'analytic':analytic,'unit':old['unit'],'branchProofs':proofs,'firstDerivative':sp.diff(analytic,w),'secondDerivative':sp.diff(analytic,w,2)}
        raw=base/'raw-records'/(key+'.pickle');raw.parent.mkdir(exist_ok=True);f.atomic_pickle(raw,record)
        pairs={}
        for name,frequency,expected in (('seed',src['referenceFrequency'],old['acceptedBinding']),('control',old['controlFrequency'],old['controlExpected'])):
            actual=analytic.subs({w:frequency,**src['origin']},simultaneous=True)
            pairs[name]=seed_comparison(base/'comparisons'/key/name,actual,expected,old['unit'],chart,frequency)
        limits=tuple(sorted({tuple(v.limits) for v in analytic.atoms(sp.Integral)},key=sp.default_sort_key))
        f.require(limits==old['census']['integralLimits'],'all ordered source/profile/momentum limits unchanged')
        record.update(bindingComparisons=pairs,limits=limits);path=directory/(key+'.pickle');f.atomic_pickle(path,record);inventory[key]={'path':str(path.relative_to(base)),'sha256':f.digest(path),'bytes':path.stat().st_size};f.save(base/'record-inventory.json',inventory);records[key]=record
    f.atomic_pickle(base/'analytic-sources.pickle',records);return records


def path_checks(base,chart):
    w=chart['frequency'];results=[]
    for root in chart['roots'].values():
        p=root['momentum'];carrier=root['carrier'];relation=carrier**2-root['radicand'];native=engine.JointBulkSheetPath(relation,w,p,carrier,sp.sqrt(root['radicand']))
        evaluate=sp.lambdify((w,p),root['analytic'],'numpy');wrong=sp.lambdify((w,p),sp.sqrt(root['radicand']),[{'sqrt':np.lib.scimath.sqrt},'numpy'])
        for momentum in (-4.,0.,4.):
            for frequency in (1.+.1j,1.-.1j,1.1-.1j):
                trace=native.trace(((1.,momentum),(frequency,momentum)));value=complex(evaluate(frequency,momentum));residual=value-trace['END_Q'];mutation=complex(wrong(frequency,momentum))-value
                record={'rootIndex':root['index'],'momentum':momentum,'frequency':frequency,'trace':trace,'chartValue':value,'chartTransportResidual':residual,'naivePrincipalRootDifference':mutation}
                f.atomic_pickle(base/f'path-{len(results)}.pickle',record);results.append(record)
                f.require(trace['STATUS']=='TRANSPORTED' and trace['PATH_DEFINED'] and abs(residual)<1e-10*(1+abs(value)) and abs(trace['ODE_DIFFERENCE'])<1e-9*(1+abs(value)),'native joint-path/root chart join')
                if frequency.imag<0:f.require(abs(mutation)>abs(value),'actual naive-principal-root sign control')
    f.atomic_pickle(base/'path-checks.pickle',results);return results


def end_tables(base,data):
    results={}
    for label in ('LEFT','RIGHT'):
        source=data['frequencyPacket']['ends'][label]['pencil'];k,z=source['momentum'],source['radical'];rows=[]
        for i in range(5):
            for j in range(5):
                original=source['livePencil'][i,j];num,den=sp.fraction(sp.cancel(original));a=sp.Poly(num,k,z);b=sp.Poly(den,k,z)
                reconstructed=sum(v*k**p*z**t for (p,t),v in a.terms())/sum(v*k**p*z**t for (p,t),v in b.terms());residual=sp.cancel(reconstructed-original)
                f.require(residual==0,'full rational end entry reconstruction');rows.append({'row':i,'column':j,'original':original,'numeratorTerms':a.terms(),'denominatorTerms':b.terms(),'reconstructionResidual':residual})
        result={'source':source,'entries':rows,'scope':'Actual full end rational coefficients for later complete invariant-pair continuation; no end mode has been continued here.'};f.atomic_pickle(base/(label.lower()+'-rational-end.pickle'),result);results[label]=result
    return results


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W']
    def put(name,value,unit=(0,0,0)):
        body=engine.cas(value);engine.emit(PREFIX+'_'+name,engine.carrier_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,unit if callable(unit) else lambda path:unit))
    for key,record in result['records'].items():
        put(key+'_SOURCE',sp.Tuple(record['originalLive'],record['algebraic'],record['analytic']),record['unit'])
        for order,name in ((1,'firstDerivative'),(2,'secondDerivative')):put(key+'_'+name,record[name],tuple(v+(order if i==1 else 0) for i,v in enumerate(record['unit'])))
        for name,pair in record['bindingComparisons'].items():put(key+'_'+name,sp.Tuple(pair['operands']['left'],pair['operands']['right'],pair['rawResidual'],pair['normalizedResidual']),record['unit'])
        boundary.structural_flags(PREFIX+'_'+key+'_DOMAIN',{'address':record['address'],'branchProofs':[v['normalizedResidual'] for v in record['branchProofs']],'pairProofs':{k:v['proofResiduals'] for k,v in record['bindingComparisons'].items()},'limits':record['limits']})
    chart=result['chart'];w=chart['frequency'];coordinate=sp.Symbol('s11cdFrequencyChartCoefficient',complex=True);engine.PHYSICAL_METADATA.dimensions.known[coordinate]=(0,0,0)
    put('FREQUENCY',w,(0,-1,0));put('COEFFICIENT_FRAME_DISK',sp.Tuple(chart['center'],chart['radius'],chart['frequencyModulusUpperBound']))
    for root in chart['roots'].values():
        put('ROOT_'+str(root['index']),sp.Tuple(root['carrier'],root['analytic'],root['seed'],root['seedResidual']),(0,-1,0))
        put('ROOT_RELATION_'+str(root['index']),sp.Tuple(root['radicand'],root['relationResidual'],root['quadraticResidual']),(0,-2,0))
        put('ROOT_COEFFICIENT_BOUND_'+str(root['index']),sp.Tuple(root['innerRealPartLowerBound'],root['modulusLowerBound']))
    for i,v in enumerate(result['denominators']['records']):
        # Bound coefficient expressions have their physical reference-unit map
        # explicitly separated, as in the accepted end candidate diagnostic.
        mapping={w:coordinate,**{root['carrier']:sp.Symbol('s11cdRootCoefficient'+str(root['index']),complex=True) for root in chart['roots'].values()}}
        for symbol in mapping.values():engine.PHYSICAL_METADATA.dimensions.known[symbol]=(0,0,0)
        body=sp.Tuple(v['algebraic'],v['familyExpression'],v['constant'],v['factorResidual'],v['coefficientFrameModulusLowerBound']).xreplace(mapping)
        put('DENOMINATOR_COEFFICIENT_FRAME_'+str(i),body)
    for label,record in result['ends'].items():
        rows=result['rowUnits'];fields=result['fieldUnits']
        for v in record['entries']:put(label+'_ENTRY_'+str(v['row'])+'_'+str(v['column']),sp.Tuple(v['original'],v['reconstructionResidual']),tuple(a-b for a,b in zip(rows[v['row']],fields[v['column']])))
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'scope':result['scope'],'chartScope':chart['scope'],'coordinateConvention':chart['coordinateConvention'],'denominatorScope':result['denominators']['scope'],'pathChecks':[{'root':v['rootIndex'],'momentum':v['momentum'],'frequency':str(v['frequency']),'status':v['trace']['STATUS'],'relativeDifference':abs(v['chartTransportResidual'])/(1+abs(v['chartValue'])),'wrongRootResponds':abs(v['naivePrincipalRootDifference'])>abs(v['chartValue'])} for v in result['paths']]})


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('frequency chart budget; preserve completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic();data=load(base);chart=root_chart(base,data);denominators=denominator_chart(base,data,chart);records=source_chart(base,data,chart);paths=path_checks(base,chart);ends=end_tables(base,data)
    result={'chart':chart,'denominators':denominators,'records':records,'paths':paths,'ends':ends,'fieldUnits':data['ends']['fieldUnits'],'rowUnits':data['ends']['rowUnits'],'sourceFiles':data['pins'],'inputPackets':data['operands'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),'scope':data['manifest']['scope']}
    f.atomic_pickle(base/'frequency-chart.pickle',result);before=f.digest(base/'frequency-chart.pickle');engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as out,contextlib.redirect_stdout(out):
        emit_result(result,data['r']);keys={tag:'s11cdFrequencyChart'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique chart tags');entries[tag]=grades._restore(body)
    old=engine.emit;seen=set()
    def replay(tag,value):
        key='PY_S11CD_'+tag;f.require(key not in seen and entries[key]==engine.cas(value),('full frequency chart emission replay',key));seen.add(key)
    engine.emit=replay
    try:emit_result(result,data['r']);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=old
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'chart payload/key coverage')
    count=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES','_DOMAIN'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)};f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']) and 'MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'complete chart metadata');count+=1 if structural else len(fields['PATHS'])
    last='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[last]},list(entries)[:list(entries).index(last)])
    f.require(before==f.digest(base/'frequency-chart.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'chart packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'all chart source/input post hashes')
    checks={'runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'records':len(records),'radicals':len(chart['roots']),'denominators':len(denominators['records']),'pathChecks':len(paths),'maximumPathDifference':max(abs(v['chartTransportResidual']) for v in paths),'lowerHalfWrongRootControls':sum(v['frequency'].imag<0 and abs(v['naivePrincipalRootDifference'])>abs(v['chartValue']) for v in paths),'endEntries':sum(len(v['entries']) for v in ends.values()),'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':count,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'frequency-chart.pickle'),'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
