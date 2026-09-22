#!/usr/bin/env python3
"""Bounded saved-input/receipt review of the one-row numerical pilot."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
import S11c_d_remaining_case_frequency_row_pilot_finish as finish

p=finish.prior
io=p.io
M,F,require,same=p.M,p.F,p.require,p.same
ROOT=F/'row-pilot-recovery-02'
CHECKS_SHA='5800400c5fd181d7bbeb08b97dc05f111b5623c1d8c30ef92229de15394b3a1c'
FINISH_SHA='8d53210471fce736e399eac96fa955eb386e413f8c8be7d14a57a93296ffa03d'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();reader=io.Reader()
    inspection=reader.json(ROOT/'completion-inspection.json');checks=reader.json(ROOT/'complete/checks.json',CHECKS_SHA)
    require(inspection['checksSha256']==CHECKS_SHA and inspection['actualExits']==[0,0,0] and inspection['emptyStrictStderr'] and inspection['checksStdoutIdentity'] and inspection['zeroCapOOMSwap'],'final clean producer inspection')
    for record in inspection['evidence'].values():reader.retain(record['path'],record['sha256'])
    g=reader.json(ROOT/'resource-guard/outcome.json');c=reader.json(ROOT/'resource-guard/child-outcome.json');n=reader.json(ROOT/'frequency_row_pilot_finish.invocation.json')
    require(g['exitCode']==c['exitCode']==n['exitCode']==0 and c['guardReason'] is None,'actual final exits')
    for path,h in ((Path(finish.__file__),FINISH_SHA),(Path(finish.recovery.__file__),finish.RECOVERY_SHA),(Path(p.__file__),finish.recovery.HELPER_SHA)):
        reader.retain(path,h)
    for path in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_pilot_review_plan.md',M/'S11c_d_remaining_case_frequency_row_pilot_finish_plan.md',M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    # Compile/reverse metadata only; no producer main/restore/evaluator is called.
    adapter,join=finish.adapted(Path(p.__file__).read_text());compile(adapter,'<saved-row-review-compile-only>','exec')
    def forbidden(*args,**kwargs):raise RuntimeError('scientific computation forbidden during saved row review')
    p.evaluate=p.source_matrix=p.independent_source=p.integrate.quad_vec=forbidden
    p.storage.Journal.write=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    def packet(record):return reader.packet(record.get('path',record.get('logical')),record['sha256'])
    def artifact(name):
        record=checks['artifacts'][name]
        return reader.json(record['path'],record['sha256']) if name.endswith('.json') else packet(record)
    for name,record in checks['artifacts'].items():require(reader.retain(record['path'],record['sha256'])['bytes']==record['bytes'],'actual artifact byte identity')
    inputs=artifact('inputs.json')
    for logical,record in inputs['consumedRoutes'].items():require(reader.retain(logical,record['sha256'])==record,'unchanged consumed logical/canonical input')
    require(artifact('finish-whole-helper-joins.json')==join,'whole saved final continuation joins')
    callers=artifact('native-and-new-numerical-callers.json')
    for name,record in callers['native'].items():
        path=reader.retain(record['file']['logical'],record['file']['sha256'])
        module=ast.parse(Path(path['canonical']).read_text())
        for fn,source in record['bodies'].items():
            current=next(node for node in module.body if getattr(node,'name',None)==fn)
            require(ast.dump(current)==ast.dump(ast.parse(source).body[0]),'whole actual numerical source/caller identity')
    native=artifact('recovery-native-power-join.json');require(native['literalPower']=='lambda x, p: np.asarray(x, dtype=complex) ** p' and not native['unsupportedNodes'] and not native['nativeCompilerCalled'],'accepted source expression branch route')
    for key in ('current','frozen'):reader.retain(native[key]['logical'],native[key]['sha256'])
    row=artifact('row-input.pickle');scope=artifact('pilot-scope.json');selected=scope['physicalInput'];source=artifact('source-action/input.pickle')
    raw=packet(selected['sourceInputs']['packet']);scalar=packet(selected['scalarInputs']);jet_route=selected['jets'][str(scope['sourceIndex'])];jet=packet(jet_route['packet'])
    for key in jet_route['keys']:jet=jet[key]
    ri,si=scope['rowIndex'],scope['sourceIndex']
    require(scope['case']=='LAB_HELD__RHOBR_CONSTANT' and ri==28 and si==22 and selected['status']=='UNMATCHED_FULL_ROW_INPUT','one actual missing row owner')
    require(same(row['row'],raw['bound']['rows'][ri]) and same(row['coefficient'],scalar['actual']['factor',ri,0]) and same(row['source'],raw['bound']['sources'][0,si]),'actual full native row/coefficient/source input')
    require(same(row['jet'],jet) and same(jet['originalBoundAmplitude'],scalar['actual']['source',si]) and same((jet['amplitudeUnit'],jet['integralUnit']),(row['source']['amplitudeUnit'],row['source']['integralUnit'])),'actual new source coefficient and complete units')
    require(same(row['fieldUnits'],raw['fieldUnits']) and same(row['equationUnits'],raw['equationUnits']) and same(row['settings'],raw['settings']),'complete field-equation/settings inputs')
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256']);require(same(physical,row['physicalRoutes']),'complete own physical context routes')
    common=packet(physical['context']);require(same(row['context'],common['contextPair'][0]),'actual own saved context')
    require(same(source['context'],row['context']) and same(source['jet'],jet) and source['size']==129 and source['bound']==row['settings']['sourceBound']==64,'full new source-action caller')
    # Compare stored arrays and their sizes; never regenerate the source rule,
    # basis/derivatives, Fourier transform, point product or quadrature.
    source_packet=packet(source['nodeRoute']['packet'])
    require(same(source['nodes'],source_packet['original']['source_nodes']) and same(source['weights'],source_packet['original']['source_weights']),'exact original source-rule nodes and measure')
    del source_packet
    matrix=artifact('source-action/value.pickle');other=artifact('source-action/independent-value.pickle');comparison=artifact('source-action/comparison.json')
    require(matrix.shape==other.shape==(1024,129) and matrix.dtype==other.dtype==np.dtype(complex) and np.isfinite(matrix).all() and np.isfinite(other).all(),'full saved source action results')
    require(comparison['maximumScaledDifference']==checks['sourceBasisScaledDifference']<comparison['tolerance']==2e-10,'saved independent source-action comparison')
    values=artifact('source-action/coefficient-values.pickle');require(len(values)==len(jet['coefficients'])==2,'actual source order')
    for order,a in enumerate(jet['coefficients']):
        prefix='source-action/coefficients/'+str(order);inp=artifact(prefix+'/input.pickle');value=artifact(prefix+'/value.pickle');receipt=artifact(prefix+'/completed.json')
        require(same(inp['expression'],a) and same(inp['variable'],row['context']['zp']) and same(value,values[order]) and value.shape==(1024,),'full actual saved coefficient evaluation')
        require(receipt['value']['sha256']==checks['artifacts'][prefix+'/value.pickle']['sha256'],'actual saved coefficient completion')
    io.save(base,'validated-source-and-native-input.json',{'case':scope['case'],'rowIndex':ri,'sourceIndex':si,'sourceActionShape':list(matrix.shape),'savedComparison':comparison,'nativeBranch':native['literalPower'],'allFullInputsAndUnitsJoined':True,'newScientificCalls':0})
    del raw,scalar,matrix,other,values
    point_indices=sorted(int(name.split('/')[1]) for name in checks['artifacts'] if name.startswith('integrand/') and name.endswith('/completed.json'))
    require(point_indices==list(range(1,6337)) and checks['newIntegrandEvaluations']==6336 and checks['newIntegrandsThisContinuation']==2240 and checks['savedPointLoadsThisContinuation']==1555,'actual complete point census')
    second_input=artifact('quadrature/1/input.pickle');variable=row['row']['limits'][0][0];keys=set();point_routes=[]
    for index in point_indices:
        prefix='integrand/'+str(index);inp=artifact(prefix+'/input.pickle');receipt=artifact(prefix+'/completed.json')
        require(receipt['input']==checks['artifacts'][prefix+'/input.pickle'] and receipt['value']==checks['artifacts'][prefix+'/value.pickle'],'immediate full point receipt')
        expected={'momentum':inp['momentum'],'variable':variable,'positions':selected['basis'],'rowInput':second_input['rowInput'],'sourceAction':second_input['sourceAction']}
        require(same(inp,expected) and type(inp['momentum']) is float and -4<=inp['momentum']<=4,'whole point caller including exact source/position routes')
        key=inp['momentum'].hex();require(key not in keys,'actual complete first point owner');keys.add(key)
        coefficient=artifact(prefix+'/coefficient-value.pickle');frequency=artifact(prefix+'/frequency-value.pickle');factors=artifact(prefix+'/factors.pickle');value=artifact(prefix+'/value.pickle')
        require(same(coefficient,factors['coefficient']) and same(frequency,factors['sourceFrequency']) and abs(frequency.imag)<1e-14,'saved literal coefficient/frequency factors')
        require(coefficient.shape==factors['sourceFourier'].shape==(129,) and value.shape==(129,129) and coefficient.dtype==factors['sourceFourier'].dtype==value.dtype==np.dtype(complex),'complete typed factor and full row sizes')
        require(np.isfinite(coefficient).all() and np.isfinite(factors['sourceFourier']).all() and np.isfinite(value).all(),'finite saved point values')
        point_routes.append({'index':index,'momentumHex':key,'input':receipt['input'],'value':receipt['value']})
        if index%256==0:io.save(base,'validated-points/'+str(index)+'.json',{'firstIndex':index-255,'lastIndex':index,'fullInputsAndTypedReturnsJoined':True})
    io.save(base,'validated-point-catalogue.json',point_routes)
    # Read successful full integrations and their diagnostic interval cache.
    # SciPy may mark evicted interval telemetry with NaN; never fill it in.
    rules=[];full_rows=[]
    for index,rule in enumerate(('gk21','gk15')):
        prefix='quadrature/'+str(index);inp=artifact(prefix+'/input.pickle');value=artifact(prefix+'/value.pickle');receipt=artifact(prefix+'/completed.json')
        require(receipt['input']==checks['artifacts'][prefix+'/input.pickle'] and receipt['value']==checks['artifacts'][prefix+'/value.pickle'],'actual full quadrature completion')
        require(inp['rule']==rule and inp['interval']==(-4.,4.) and inp['norm']=='max' and inp['workers']==1 and inp['limit']==256 and inp['cacheSize']==32*1024**2,'actual declared finite quadrature route')
        require(inp['epsabs']==(1e-8,5e-9)[index] and inp['epsrel']==(1e-6,5e-7)[index],'unchanged rule tolerance')
        require(value['success'] and value['status']==0 and value['row'].shape==(129,129) and np.isfinite(value['row']).all() and np.isfinite(value['errorEstimate']),'successful full numerical row return')
        intervals=value['intervals'];errors=value['intervalErrors'];cached=value['intervalValues']
        require(intervals.shape==(len(errors),2) and cached.shape==(len(errors),129,129) and np.isfinite(intervals).all() and np.isfinite(errors).all() and np.all(intervals[:,0]<intervals[:,1]),'actual interval telemetry shapes/order')
        finite=np.isfinite(cached).all(axis=(1,2));evicted=np.isnan(cached).all(axis=(1,2));require(np.all(finite|evicted),'saved complete-or-evicted interval cache diagnostics')
        summary={'rule':rule,'evaluations':value['evaluations'],'intervalCount':len(intervals),'evictedDiagnosticIntervals':int(np.count_nonzero(evicted)),'errorEstimate':value['errorEstimate'],'wallSeconds':value['wallSeconds'],'fullRowMaxNorm':float(np.max(abs(value['row']))),'input':receipt['input'],'value':receipt['value']}
        rules.append(summary);full_rows.append(value['row']);io.save(base,'validated-rule-'+rule+'.json',summary)
    diff=artifact('row-comparison.pickle')
    require(diff['first']==checks['artifacts']['quadrature/0/value.pickle'] and diff['second']==checks['artifacts']['quadrature/1/value.pickle'],'full saved comparison input addresses')
    require(diff['difference'].shape==(129,129) and np.isfinite(diff['difference']).all() and diff['maximumAbsoluteDifference']==checks['rowAbsoluteDifference'] and diff['maximumScaledDifference']==checks['rowScaledDifference']<1e-6 and diff['rowMaxNorm']==checks['rowMaxNorm'],'saved complete rule comparison')
    ci=artifact('controls/input.pickle');cv=artifact('controls/value.pickle');cs=artifact('controls/checks.json')
    require(ci['wrongPhaseSign']==1 and ci['weightScale']==1.001 and cs==checks['controlResponses'] and all(v>1e-12 for v in cs.values()),'saved actual phase and measure responses')
    actual=packet(ci['actualIntegrand']);require(all(v.shape==actual.shape==(129,129) and np.isfinite(v).all() and not np.array_equal(v,actual) for v in cv.values()),'actual changed control arrays and complete sizes')
    io.save(base,'validated-rule-comparison-and-controls.json',{'absoluteRuleDifference':checks['rowAbsoluteDifference'],'scaledRuleDifference':checks['rowScaledDifference'],'rowMaxNorm':checks['rowMaxNorm'],'controls':cs,'comparisonArithmeticRepeated':False})
    provenance=artifact('completed-point-and-failure-provenance.json');require(provenance['failure']['actualExits']==[1,1,1] and not provenance['completedFirstRuleRecomputed'] and not provenance['controlsRecomputed'] and not provenance['sourceActionRecomputed'],'preserved actual failure and completed-prefix provenance')
    require(checks['completedGK21Reused'] and checks['completedControlsReused'] and checks['newSourceActions']==0 and checks['priorUnfinishedRuleReuseCountNotReconstructed'],'explicit actual operation scope')
    io.save(base,'validated-continuation-provenance.json',{'joins':join,'preservedFailures':True,'completedSourceGK21ControlsReused':True,'newPointsInFinish':2240,'savedPointLoadsInFinish':1555,'oldPartialAdaptiveStateNotClaimedRestored':True})
    reader.postcheck();io.save(base,'validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_FREQUENCY_ROW_PILOT_REVIEW','producerChecksSha256':CHECKS_SHA,'case':scope['case'],'rowIndex':ri,'sourceIndex':si,'fullPointReceipts':len(point_indices),'rules':rules,'rowMaxNorm':checks['rowMaxNorm'],'rowAbsoluteDifference':checks['rowAbsoluteDifference'],'rowScaledDifference':checks['rowScaledDifference'],'sourceBasisScaledDifference':checks['sourceBasisScaledDifference'],'controlResponses':cs,'newScientificCalls':0,'allConsumedHashesUnchanged':True,'consumedLogicalPaths':len(reader.routes),'wallSeconds':time.monotonic()-started,'scope':'Saved full input/value/receipt/native-source/unit review of one fixed-setting129x129 row at1-0.01i. No evaluator, source, branch, basis, quadrature, map or response recomputation. Numerical spread is not a scattering error bound or physical tiny-effect claim.'}
    io.save(base,'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
