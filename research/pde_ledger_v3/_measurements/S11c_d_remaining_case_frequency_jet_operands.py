#!/usr/bin/env python3
"""Inspect only unmatched saved jet operands; no extraction or arithmetic."""
import argparse,json,resource,signal,time
from pathlib import Path
import sympy as sp
import S11c_d_remaining_case_frequency_end_continuation_inputs as io
M,REPO=io.M,io.REPO
CP=M/'S11c_d_remaining_case_frequency_scalar_bindings_checkpoint.json'
SHA='d1961c25f725959e896e62e49c45d389fdd503575d27803a89055175f6308b50'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path)
    base=ap.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    reader=io.Reader();cp=reader.json(CP,SHA);origin=Path(cp['runDirectory']);review=cp['independentSavedReview'];vr=Path(review['runDirectory'])
    io.require(cp['status']=='ACCEPTED_CASE_FREQUENCY_SCALAR_BINDINGS','accepted scalar operands')
    reader.retain(origin/'checks.json',cp['checksSha256']);reader.retain(vr/'checks.json',review['checksSha256'])
    reader.retain(Path(__file__).resolve());reader.retain(M/'S11c_d_remaining_case_frequency_jet_operands_plan.md')
    name='source-jet-input-catalogue.json';catalogue=reader.json(vr/name,review['artifacts'][name]['sha256'])
    cache={}
    def packet(record):
        p=record.get('logical',record.get('path'));r=reader.retain(p,record['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(p,record['sha256'])
        return cache[r['canonical']]
    # Read complete saved coefficient values from the exact catalogued packets.
    atlas=[]
    for request in catalogue['requests']:
        old=packet(request['nativeSource']['packet'])
        for si,jet in old['jets'].items():
            owner=(request['nativeSource']['packet']['canonical'],si)
            if not any(v[0]==owner for v in atlas):atlas.append((owner,jet))
        for route in request['matches']:
            value=packet(route['packet']);key=route['keys']
            jet=value[key[0]][key[1]]
            owner=(route['packet']['canonical'],key[1])
            if not any(v[0]==owner for v in atlas):atlas.append((owner,jet))
    unique=[];summaries=[]
    for request in catalogue['requests']:
        if request['matches']:continue
        case,si=request['case'],request['sourceIndex'];old=packet(request['nativeSource']['packet']);binding=packet(request['boundAmplitude']['packet'])
        expr=binding['actual']['source',si];prior=old['jets'][si];source=old['bound']['sources'][0,si];probe=prior['probe'];coordinate=probe.args[0]
        signature=(expr,source['amplitudeUnit'],source['integralUnit'],probe)
        matches=[i for i,(value,_) in enumerate(unique) if io.same(signature,value)]
        index=matches[0] if matches else len(unique)
        if not matches:unique.append((signature,request))
        derivatives=list(expr.atoms(sp.Derivative))
        io.require(all(v.expr==probe and all(x==coordinate for x,_ in v.variable_count) for v in derivatives),'actual derivative node coordinates')
        degree=max([0]+[sum(n for _,n in v.variable_count) for v in derivatives])
        derivative_routes=[]
        for order in range(degree+1):
            nodes=[]
            if order==0:nodes.append({'kind':'actual-saved-probe','packet':request['ownProbeAndCoordinate']})
            for node in derivatives:
                if sum(n for _,n in node.variable_count)==order:nodes.append({'kind':'actual-derivative-node-in-bound-amplitude','value':repr(node),'variableCount':repr(node.variable_count)})
            # Read the original actual source amplitude too; no diff is called.
            for owner,jet in atlas:
                if not io.same(jet['probe'],probe):continue
                for node in jet['originalBoundAmplitude'].atoms(sp.Derivative):
                    if node.expr==probe and sum(n for _,n in node.variable_count)==order:
                        nodes.append({'kind':'actual-derivative-node-in-saved-original-amplitude','owner':[owner[0],owner[1]],'value':repr(node)})
            derivative_routes.append({'order':order,'availableNodes':nodes})
        coeff_matches=[];stack=[((),expr)]
        while stack:
            address,node=stack.pop()
            for owner,jet in atlas:
                for order,coefficient in enumerate(jet['coefficients']):
                    if type(node) is type(coefficient) and node==coefficient:
                        coeff_matches.append({'argsAddress':list(address),'savedOwner':[owner[0],owner[1]],'coefficientIndex':order,'value':repr(node)})
            stack.extend((address+(i,),v) for i,v in enumerate(node.args))
        out={'case':case,'sourceIndex':si,'inputFamily':index,'firstOwner':[unique[index][1]['case'],unique[index][1]['sourceIndex']],
             'fullInputRoutes':request,'actualExpression':repr(expr),'expressionTree':sp.srepr(expr),
             'probe':repr(probe),'sourceCoordinate':repr(coordinate),'degree':degree,'actualDerivativeNodeRoutes':derivative_routes,
             'actualCoefficientSubtreeMatches':coeff_matches,'previousOwnCoefficients':[repr(v) for v in prior['coefficients']],
             'noDerivativeOrCoefficientConstructed':True,
             'scope':'Saved raw operands only; a node match is not a completed source-jet or numerical row acceptance.'}
        io.save(base,'cases/'+case+'/source-'+str(si)+'.json',out)
        summaries.append({k:out[k] for k in ('case','sourceIndex','inputFamily','firstOwner','degree')})
    reader.postcheck();io.save(base,'inputs.json',{'acceptedScalarCheckpoint':reader.retain(CP,SHA),'consumedRoutes':reader.routes})
    result={'status':'COMPLETED_SAVED_UNMATCHED_SOURCE_JET_OPERAND_INSPECTION','unmatchedUses':len(summaries),'literalInputFamilies':len(unique),
            'cases':summaries,'existingJetValuesInspected':len(atlas),'newScientificCalls':0,'allConsumedInputsUnchanged':True,'wallSeconds':time.monotonic()-start,
            'scope':'Only saved expression/probe/derivative/coefficient subtree routes. Inspect actual operands before choosing genuinely missing extraction or numerical evaluation; no science acceptance.'}
    io.save(base,'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
