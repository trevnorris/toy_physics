#!/usr/bin/env python3
"""Validate and atomically publish a computed four-case coordinate export check."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import pickle
import shutil
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'_measurements')]
import sympy as sp
from ledger_fold import _restore
from S11c_d_modal_current_check import engine,digest


def assoc(value):return {str(k):v for k,v in value}


def nonfinite_census(body):
    """Distinguish explicit improper-integration domains from nonfinite values.

    c2 exports retain definite integrals on the whole real line. Their exact
    +/-oo endpoint nodes describe domains; they are not evaluated coefficients.
    No other occurrence of infinity, NaN or complex infinity is exempted.
    """
    visited=set();bad=set();endpoints=0
    def visit(node):
        nonlocal endpoints
        if node in visited:return
        visited.add(node)
        if node in (sp.nan,sp.zoo,sp.oo,-sp.oo):
            bad.add(sp.srepr(node));return
        if isinstance(node,sp.Integral):
            visit(node.function)
            for limit in node.limits:
                visit(limit[0])
                for endpoint in limit[1:]:
                    if endpoint in (sp.oo,-sp.oo):endpoints+=1
                    else:visit(endpoint)
        else:
            for arg in node.args:visit(arg)
    visit(body)
    return {'nonfiniteValueNodes':sorted(bad),'improperIntegralEndpointNodes':endpoints}


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--transcript',type=Path,required=True)
    parser.add_argument('--publish',action='store_true')
    parser.add_argument('--stage',choices=('b','c2'),default='b')
    args=parser.parse_args();base=args.run_directory
    summary=json.loads((base/'checks.json').read_text())
    if summary['stage']!=args.stage:raise ValueError('export stage mismatch')
    key_prefix='s11cThicknessCoordinate'+{'b':'B','c2':'C2'}[args.stage]
    if summary['sourcePins']!=summary['sourcePinsAfter']:raise ValueError('source changed during trace')
    for name,sha in summary['sourcePins'].items():
        if digest(ROOT/name)!=sha or digest(base/summary['sourceSnapshots'][name])!=sha:raise ValueError(('source pin',name))
    if digest(base/'objects.pickle')!=summary['objectsSha256']:raise ValueError('trace packet pin')
    objects=pickle.loads((base/'objects.pickle').read_bytes())
    lines=[json.loads(line) for line in args.transcript.read_text().splitlines()]
    if len(lines)!=len(objects) or len(lines)!=summary['objects']:raise ValueError('trace census')
    paths=0;keys=[];nonfinite=[];computed_checks={};integration_domains={}
    for item,line in zip(objects,lines):
        actual={k:_restore(v) for k,v in line.items()}
        expected={k:engine.cas(v) for k,v in item['record'].items()}
        if actual!=expected:raise ValueError(('trace payload',item['key']))
        body=item['body'];keys.append(item['key'])
        value=engine.carrier_fingerprint(body) if str(actual['representation'])=='CARRIER_PIT_SHA' else body
        if value!=actual['value']:raise ValueError(('trace fingerprint',item['key']))
        domain=nonfinite_census(body)
        if domain['nonfiniteValueNodes']:nonfinite.append(item['key'])
        if domain['improperIntegralEndpointNodes']:
            integration_domains[item['key']]=domain['improperIntegralEndpointNodes']
        metadata=[assoc(v) for v in actual['metadata']]
        source_paths=[p for p,v in engine.leaves(body) if not isinstance(v,engine.Str)]
        if Counter(tuple(int(v) if not isinstance(v,engine.Str) else str(v) for v in m['path']) for m in metadata)!=Counter(source_paths):
            raise ValueError(('metadata paths',item['key']))
        for m in metadata:
            if len(m['dimensionLTM'])!=3 or any(not v.is_Integer for v in m['dimensionLTM']):raise ValueError('dimension')
            if not all(k in m for k in ('multigrade','epsilonLambdaSupport')):raise ValueError('grade data')
        paths+=len(metadata)
        if item['key'].endswith('Residual'):
            values=[v for _,v in engine.leaves(body) if not isinstance(v,engine.Str)]
            name=item['key'].removeprefix(key_prefix)
            computed_checks[name]={'scalars':len(values),'nonzero':sum(v!=0 for v in values)}
    if len(keys)!=len(set(keys)) or nonfinite:raise ValueError('duplicate/nonfinite trace')
    if computed_checks!=summary['checks']:raise ValueError('computed residual census')
    residual_scalars=sum(v['scalars'] for v in summary['checks'].values())
    nonzero=sum(v['nonzero'] for v in summary['checks'].values())
    if nonzero!=summary['nonzeroResidualScalars']:raise ValueError('nonzero residual total')
    if set(map(tuple,summary['cases']))!={(a,r) for a in ('LAB_HELD','MATERIAL_ADVECTED') for r in ('RHO4_CONSTANT','RHOBR_CONSTANT')}:raise ValueError('case coverage')
    inventory={**summary,'runDirectory':str(base),'validatorSha256':digest(Path(__file__)),
        'metadataPaths':paths,'residualScalars':residual_scalars,'nonzeroResidualScalars':nonzero,
        'nonfiniteObjects':nonfinite,'transcript':{'bytes':args.transcript.stat().st_size,'sha256':digest(args.transcript)}}
    inventory['improperIntegralEndpointNodesByObject']=integration_domains
    print(json.dumps({k:inventory[k] for k in ('objects','metadataPaths','residualScalars','nonzeroResidualScalars','transcript')},indent=2))
    if args.publish:
        stem='S11c_thickness_coordinate_'+args.stage+'_export'
        target=ROOT/'scripts/out'/(stem+'.out')
        if target.exists() or target.is_symlink():raise ValueError('publication exists')
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.thickness-coordinate-',delete=False) as f:
            tmp=Path(f.name)
            with args.transcript.open('rb') as source:shutil.copyfileobj(source,f)
            f.flush();os.fsync(f.fileno())
        if digest(tmp)!=inventory['transcript']['sha256']:raise ValueError('publication digest')
        os.replace(tmp,target)
        inventory['publication']=str(target.relative_to(ROOT))
        (ROOT/'_measurements'/(stem+'_checkpoint.json')).write_text(json.dumps(inventory,indent=2)+'\n')
    (base/'validation.json').write_text(json.dumps(inventory,indent=2)+'\n')
    if nonzero:raise ValueError('trace residual; inspect saved operands')


if __name__=='__main__':run()
