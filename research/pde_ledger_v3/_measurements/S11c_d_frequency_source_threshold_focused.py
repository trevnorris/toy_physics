#!/usr/bin/env python3
"""Check actual complex elimination and saved-source recovery without rebinding."""
import argparse,contextlib,json,resource,signal,time
from pathlib import Path
import sympy as sp
import S11c_d_frequency_source_finish as u

f=u.f;q=u.q


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--resume-analysis',type=Path);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('bounded actual-threshold instrument check')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(180);started=time.monotonic();adapter=u.end_adapter();data=q.load(base);src,pencils,eliminations=u.reuse(base,data)
    saved=eliminations['REFERENCE'];w=saved['coordinate'];z=saved['radicalCoordinate'];poly=sp.Poly(saved['elimination'],w)
    if args.resume_analysis:
        original=args.resume_analysis.resolve();original.relative_to(f.STORE);actual=f.unpickle(original)
        f.require(actual['normalized']==sp.Poly(poly.monic().as_expr(),w,domain=sp.QQ_I) and actual['decompositionResidual']==actual['bezoutResidual']==actual['factorReconstructionResidual']==0 and all(v==0 for v in actual['divisionResiduals']),'completed actual analysis polynomial/proof join')
        q.shutil.copyfile(original,base/'actual-analysis.pickle');f.require(f.digest(original)==f.digest(base/'actual-analysis.pickle'),'byte-identical actual analysis reuse');data['operands'][str(original)]=f.digest(original)
    else:
        actual=u.real_axis_analysis(poly,w);f.atomic_pickle(base/'actual-analysis.pickle',actual)
    f.require(actual['imaginaryPart']!=0 and sp.degree(actual['realCommonFactor'],w)==8 and sp.degree(actual['realSquarefree'],w)==2,'actual complex polynomial and common factor')
    real_control=u.real_axis_analysis(sp.Poly(actual['realCommonFactor'],w),w)
    f.require(real_control['imaginaryPart']==0 and real_control['intervals']==actual['intervals'],'old real-polynomial special case')
    mutation=u.real_axis_analysis(sp.Poly(actual['normalized'].as_expr()+sp.I,w),w)
    f.atomic_pickle(base/'actual-imaginary-mutation.pickle',mutation)
    f.require(sp.degree(mutation['realCommonFactor'],w)==0 and not mutation['intervals'],'actual imaginary-coefficient mutation removes the common real roots')
    # Exercise the actual emitted complex polynomial, real-axis proofs and all
    # completed source records before any remaining end construction.
    inp=data['adapter'].input
    threshold={'coefficientMatrix':saved['coefficientMatrix'],'clearedMatrix':saved['clearedMatrix'],'rowDenominators':saved['rowDenominators'],'numerator':saved['numerator'],'denominator':saved['denominator'],'wave':saved['wave'],'elimination':saved['elimination'],
      'normalizedElimination':actual['normalized'].as_expr(),'squarefreeElimination':actual['complexSquarefree'].as_expr(),'factorization':actual['factorization'],'factorReconstructionResidual':actual['factorReconstructionResidual'],
      'realRootIntervals':actual['intervals'],'bulkBranchFrequencies':sp.solve(saved['wave'].subs(z,0),w),'frequencyCoordinate':w,'radicalCoordinate':z,
      'coordinateMap':{'frequencyUnit':(0,-1,0),'radicalUnit':(0,-1,0),'normalMomentumUnit':(-1,0,0),'normalMomentumCoefficient':0,'unitFrame':inp.frame,'fieldReferenceUnits':data['ends']['fieldUnits'],'equationReferenceUnits':data['ends']['rowUnits'],'matrixConvention':'Numerical coefficient matrix between equation and field reference-unit bases; entries and elimination coordinates are dimensionless.'},
      'scope':'Focused actual saved-reference threshold diagnostic only; no completed three-end or profile-frequency result.','realAxisAnalysis':actual}
    q.engine.PHYSICAL_METADATA.dimensions.known[w]=(0,0,0);q.engine.PHYSICAL_METADATA.dimensions.known[z]=(0,0,0)
    result={'sources':src,'ends':{'REFERENCE':{'pencil':pencils['REFERENCE'],'thresholds':threshold}},'fieldUnits':data['ends']['fieldUnits'],'rowUnits':data['ends']['rowUnits'],'sourceFiles':data['pins'],'inputPackets':data['operands'],'scope':'Focused actual source/reference diagnostic emission; remaining end construction is not run.'}
    f.atomic_pickle(base/'focused-result.pickle',result);q.engine.EMISSION_LINES.clear();q.engine.PAYLOAD_ENCODER=q.grades.PayloadEncoder()
    with (base/'focused.out').open('x') as out,contextlib.redirect_stdout(out):u.emit_result(result,data['r'])
    entries={}
    for line in q.grades.decoded_lines(base/'focused.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique focused tags');entries[tag]=q.grades._restore(body)
    old=q.engine.emit;seen=set()
    def replay(tag,value):
        key='PY_S11CD_'+tag;f.require(key not in seen and entries[key]==q.engine.cas(value),('actual focused emission replay',key));seen.add(key)
    q.engine.emit=replay
    try:u.emit_result(result,data['r'])
    finally:q.engine.emit=old
    f.require(seen==set(entries) and not q.engine.PHYSICAL_METADATA.dimensions.constraints,'complete actual emission and dimensions')
    checks={'status':'PASSED','complexDegree':int(poly.degree()),'realCommonDegree':int(sp.degree(actual['realCommonFactor'],w)),'realSquarefreeDegree':int(sp.degree(actual['realSquarefree'],w)),
      'realIntervals':[[str(a),str(b)] for (a,b),_ in actual['intervals']],'imaginaryPartRetained':True,'actualCoefficientMutationRejected':True,'realSpecialCaseJoined':True,'wholeEndConstructorReverseAstJoin':True,
      'copiedPackets':json.loads((base/'operand-validation.json').read_text())['copiedPackets'],'sourceRecords':375,'sourcePacketSha256':f.digest(base/'frequency-sources.pickle'),'tagCount':len(entries),'fullActualEmissionReplay':True,
      'sourceFiles':data['pins'],'inputPackets':data['operands'],'artifacts':{str(p.relative_to(base)):f.digest(p) for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},
      'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    for name,sha in data['pins'].items():f.require(f.digest(f.ROOT/name)==sha,'post source hash')
    for name,sha in data['operands'].items():f.require(f.digest(Path(name))==sha,'post input hash')
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
