#!/usr/bin/env python3
"""Validate saved source quadrature and repair lossless structural emission."""
import argparse
import ast
import hashlib
import faulthandler
import multiprocessing
import os
import resource
import inspect
import json
from pathlib import Path
import shutil
import time

import numpy as np
import sympy as sp
import S11c_d_source_fourier_quadrature_check as original
from S11c_d_source_fourier_quadrature_check import (
    ROOT, STORE, engine, digest, save, atomic_pickle, unpickle, decoded_lines, _restore)
from S11c_d_output_codec import restore_emission_index

PLAN = ROOT/'_measurements/S11c_d_source_fourier_quadrature_emission_plan.md'
REPAIR = ROOT/'_measurements/S11c_d_source_fourier_quadrature_emission_repair.json'
REPLAY_REPAIR = ROOT/'_measurements/S11c_d_source_fourier_quadrature_replay_repair.json'
SOURCES = (*original.SOURCES, Path(__file__).resolve(), PLAN, REPAIR, REPLAY_REPAIR)


def checker_join(previous_text, current_text):
    """Undo only the finalize extraction and the two structural emit changes."""
    before, after = ast.parse(previous_text), ast.parse(current_text)
    members = {getattr(n, 'name', None): n for n in after.body}
    final = members['finalize']; main = members['main']
    # finalize has its docstring, two saved-result aliases, original tail,
    # and a returned summary. main calls it and prints the returned summary.
    if (not isinstance(final.body[-1], ast.Return) or
            ast.dump(final.body[-1].value) != ast.dump(ast.Name(id='summary', ctx=ast.Load())) or
            not isinstance(main.body[-2], ast.Assign) or
            not isinstance(main.body[-2].value, ast.Call) or
            getattr(main.body[-2].value.func, 'id', None) != 'finalize'):
        raise ValueError('unexpected finalization extraction')
    main.body = main.body[:-2]+final.body[3:-1]+main.body[-1:]
    after.body.remove(final); after.body.remove(members['emit_manifest'])
    class UndoManifest(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'emit_manifest':
                node.func = ast.Attribute(value=ast.Name(id='engine',ctx=ast.Load()),
                                          attr='physical',ctx=ast.Load())
            return node
    after = UndoManifest().visit(after)
    if ast.dump(before) != ast.dump(after):
        raise ValueError('checker changed beyond structural emission and tail extraction')
    return True


def equal(left, right, path=()):
    """Exact saved-packet join, including every numeric array element."""
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        ok = isinstance(left,np.ndarray) and isinstance(right,np.ndarray) and np.array_equal(left,right)
    elif isinstance(left, dict) and isinstance(right,dict):
        ok = list(left)==list(right)
        if ok:
            for key in left: equal(left[key],right[key],path+(key,))
    elif isinstance(left,(list,tuple)) and isinstance(right,(list,tuple)):
        ok = len(left)==len(right)
        if ok:
            for i,(a,b) in enumerate(zip(left,right)): equal(a,b,path+(i,))
    else:
        ok = left == right
    if not ok: raise ValueError(('saved packet join',path))


def entries(path):
    result = {}
    for line in decoded_lines(path):
        tag,_,body = line.rstrip('\n').partition(': ')
        if tag in result: raise ValueError(('duplicate tag',tag))
        result[tag] = body
    return result


def certificate_worker(directory, left, right):
    """Fork preserves the actual live operands; pickle can change their trees."""
    with (directory/'stdout').open('x') as out, (directory/'stderr').open('x') as err:
        os.dup2(out.fileno(),1); os.dup2(err.fileno(),2)
        resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
        faulthandler.dump_traceback_later(150,file=err)
        begin=time.monotonic()
        certificate=engine.BoundedSourceFourierAssembly.reconstruction_certificate(left,right,shared=False)
        atomic_pickle(directory/'certificate.pickle',certificate)
        proofs=tuple(certificate['REPLAY_RESIDUALS'])+tuple(v[1] for v in certificate['PHASE_SPLITS'].values())+tuple(
            v[2] for v in certificate['RADICAL_POWERS'].values())
        result={'liveExactEqual':left==right,'normalizedResidualZero':certificate['RESIDUAL']==0,
            'proofResidualCount':len(proofs),'nonzeroProofResiduals':sum(v!=0 for v in proofs),
            'representationSha256':[hashlib.sha256(sp.srepr(v).encode()).hexdigest() for v in (left,right)]}
        save(directory/'certificate-checks.json',result)
        mutation=engine.BoundedSourceFourierAssembly.reconstruction_certificate(left,2*right,shared=False)
        atomic_pickle(directory/'mutation.pickle',mutation)
        result.update(coefficientMutationNonzero=mutation['RESIDUAL']!=0,
            wallSeconds=time.monotonic()-begin,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        save(directory/'checks.json',result); print(json.dumps(result,indent=2),flush=True)
        faulthandler.cancel_dump_traceback_later()
        if result['liveExactEqual'] or not result['normalizedResidualZero'] or result['nonzeroProofResiduals'] or not result['coefficientMutationNonzero']:
            raise ValueError('live source replay certificate requires inspection')


def certify_replay(base,left,right,key,unit,inventory):
    """Save every actual replay pair before applying its exact equality guard."""
    directory=base/'source-binding-replay'/('-'.join(map(str,key)))
    directory.mkdir(parents=True,exist_ok=False)
    representations=tuple(sp.srepr(v) for v in (left,right))
    atomic_pickle(directory/'operands.pickle',{'left':left,'right':right,'rawResidual':left-right,
        'representations':representations,'dimensionLTM':unit,
        'sourceBoundPacketSha256':digest(base/'bound-sources.pickle')})
    item={'key':list(key),'path':str(directory.relative_to(base)),'liveExactEqual':left==right,
        'operandSha256':digest(directory/'operands.pickle')}
    inventory.append(item);save(base/'source-binding-replay-inventory.json',inventory)
    if left!=right:
        # A spawn/unpickle boundary may itself canonicalize the operands. Fork
        # keeps the live expression trees and performs no concurrent parent CAS.
        worker=multiprocessing.get_context('fork').Process(target=certificate_worker,args=(directory,left,right))
        worker.start(); worker.join(180)
        if worker.is_alive():
            worker.terminate(); worker.join(5)
            if worker.is_alive(): worker.kill(); worker.join()
            item['status']='timeout'
        else: item['status']='exited'
        item['exitCode']=worker.exitcode
        item['stderrBytes']=(directory/'stderr').stat().st_size if (directory/'stderr').exists() else None
        save(base/'source-binding-replay-inventory.json',inventory)
        if item['status']!='exited' or item['exitCode']!=0 or item['stderrBytes']!=0:
            raise ValueError(('bounded live source replay certificate',key,item['status']))
        checks=json.loads((directory/'checks.json').read_text())
        equal(checks['representationSha256'],[hashlib.sha256(v.encode()).hexdigest() for v in representations],(key,'live-worker-pair'))
        if checks['liveExactEqual'] or not checks['normalizedResidualZero'] or checks['nonzeroProofResiduals'] or not checks['coefficientMutationNonzero']:
            raise ValueError(('source replay worker guards',key))
        certificate=unpickle(directory/'certificate.pickle');mutation=unpickle(directory/'mutation.pickle')
        equal(certificate['LEFT'],unpickle(directory/'operands.pickle')['left'],(key,'certificate-left'))
        equal(certificate['RIGHT'],unpickle(directory/'operands.pickle')['right'],(key,'certificate-right'))
        equal(mutation['LEFT'],certificate['LEFT'],(key,'mutation-left'))
        equal(mutation['RIGHT'],2*certificate['RIGHT'],(key,'mutation-right'))
        item['checks']=checks
    else:
        item['normalizedResidualZero']=(left-right)==0
    item['artifacts']={p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in directory.iterdir() if p.is_file()}
    save(base/'source-binding-replay-inventory.json',inventory)


def validate_saved(previous, base, bound_packet, quadrature_packet, pencil, assembly, factors, tests, numerical, provenance):
    bound, result = bound_packet['result'], quadrature_packet['result']
    for packet in (bound_packet,quadrature_packet): equal(packet['provenance'],provenance)
    equal(quadrature_packet['boundPacketSha256'],digest(previous/'bound-sources.pickle'))
    adapter = engine.NumericalReducedAction(pencil,assembly,json.loads(original.native.INPUT.read_text()))
    original_rows = {row['INDEX']:row for row in factors['ROWS']}
    comparisons = bound['nativeTestIntegralComparisons']
    expected_census = {(ti,row['COLUMN'],row['ROW'],term) for ti in range(len(tests))
        for row in assembly['ROWS'] for term in range(len(row['NONLOCAL']))}
    actual_census = [(r['test'],r['column'],r['row'],r['term']) for r in comparisons]
    if len(actual_census)!=len(set(actual_census)) or set(actual_census)!=expected_census:
        raise ValueError('native field/row/term coverage')
    representation_records = []
    for i,item in enumerate(comparisons):
        c=item['comparison']; ti,j,row,term=(item[k] for k in ('test','column','row','term'))
        native_row=next(r for r in tests[ti]['assembled'] if (r['COLUMN'],r['ROW'])==(j,row))
        source_row=next(r for r in assembly['ROWS'] if (r['COLUMN'],r['ROW'])==(j,row))
        equal(item['original'],source_row['NONLOCAL'][term][0],(i,'original'))
        equal(item['original'],original_rows[item['sourceIndex']]['ORIGINAL'],(i,'factor'))
        equal(c['expected'],native_row['NONLOCAL'][term][1],(i,'expected'))
        if c['current'].limits!=c['expected'].limits or not c['limitsEqual']:
            raise ValueError(('original native limits',i))
        restored=tuple(_restore(text) for text in c['representationStrings'])
        for name,value in zip(('current','expected','rawResidual'),restored):
            equal(value,c[name],(i,'live-representation',name))
        if c['normalizedResidual']!=0 or any(v!=0 for v in (*c['replayResiduals'],*c['exponentResiduals'])):
            raise ValueError(('saved exact binding proof',i))
        certificate=c['certificate']
        if certificate is not None:
            equal(certificate['LEFT'],c['current'].function,(i,'certificate-left'))
            equal(certificate['RIGHT'],c['expected'].function,(i,'certificate-right'))
            equal(certificate['RESIDUAL'],c['normalizedResidual'],(i,'certificate-residual'))
            equal(certificate['REPLAY_RESIDUALS'],c['replayResiduals'],(i,'certificate-replay'))
            equal(tuple(v[1] for v in certificate['PHASE_SPLITS'].values())+tuple(
                v[2] for v in certificate['RADICAL_POWERS'].values()),c['exponentResiduals'],(i,'certificate-exponents'))
        elif c['current']!=c['expected'] or c['rawResidual']!=0:
            raise ValueError(('missing binding certificate',i))
        representation_records.append({'occurrence':i,'liveExactEqual':bool(c['exactEqual']),
            'restoredExactEqual':bool(c['current']==c['expected']),
            'representationSha256':[hashlib.sha256(s.encode()).hexdigest() for s in c['representationStrings']]})
    if len(bound['nativeTestIntegralJoins'])!=len(comparisons) or not all(bound['nativeTestIntegralJoins']):
        raise ValueError('binding join/proof census')
    record_map={(r['test'],r['sourceIndex']):r for r in bound['records']}
    if (len(record_map)!=len(bound['records']) or set(record_map)!=
            {(ti,si) for ti in range(len(tests)) for si in range(len(factors['SOURCE_INTEGRALS']))}):
        raise ValueError('distinct source/test coverage')
    equal(bound['testWidthsMomenta'],[(t['width'],t['momentum']) for t in tests])
    expected_momenta=tuple(pencil.r.normal_map[g[2]] for g in pencil.r.momentum_groups)
    equal(bound['momenta'],expected_momenta)
    settings=numerical['results'][0]['settings']; momentum_bound=sp.Rational(str(settings['momentumBound']))
    equal(bound['momentumDomains'],{k:(-momentum_bound,momentum_bound) for k in expected_momenta})
    equal(bound['sourceBounds'],(float(settings['sourceBound']),1.5*float(settings['sourceBound'])))
    equal(tuple(bound['orders']),(32,64,128))
    replay_inventory=[]
    for key,r in record_map.items():
        ti,si=key; source=factors['SOURCE_INTEGRALS'][si]
        uses=[(row['INDEX'],fi,f) for row in factors['ROWS'] for fi,f in enumerate(row['FACTORS']) if f['SOURCE_INTEGRAL']==source]
        equal(r['uses'],[(i,j) for i,j,_ in uses],(key,'uses'))
        equal(r['originalSourceIntegral'],source,(key,'integral'))
        field=lambda z:engine.memo_xreplace(tests[ti]['field'],{pencil.r.z:z})
        mapping={p:field for p in pencil.probes}
        for symbolic,saved in (('AMPLITUDE','boundAmplitude'),('SOURCE','boundSource')):
            certify_replay(base,r[saved],adapter.bind(engine.dag_substitute(uses[0][2][symbolic],mapping)),
                (*key,saved),engine.PHYSICAL_METADATA.dimensions.measure(uses[0][2][symbolic]),replay_inventory)
        equal(r['symbolicAmplitude'],uses[0][2]['AMPLITUDE'],(key,'symbolicAmplitude'))
        equal(r['symbolicFrequency'],uses[0][2]['FREQUENCY'],(key,'symbolicFrequency'))
        equal(r['frequency'],adapter.bind(uses[0][2]['FREQUENCY']),(key,'frequency'))
        equal(r['boundCharacter'],adapter.bind(uses[0][2]['CHARACTER']),(key,'character'))
        equal(r['range'],engine.BoundedSourceFourierQuadrature.affine_range(r['frequency'],bound['momentumDomains']),(key,'range'))
        fn=sp.lambdify(bound['momenta'],r['frequency'],'numpy')
        equal(r['assignmentResidual'],np.asarray([float(fn(*p))-nu for p,nu in zip(r['assignments'],r['frequencies'])]),(key,'assignments'))
        lo,hi=map(float,r['range']['bounds'])
        samples=list(np.linspace(lo,hi,65)); samples.extend(x for x in (0.,float(tests[ti]['momentum'])) if lo<=x<=hi)
        equal(r['frequencies'],np.asarray(sorted(set(samples)),dtype=float),(key,'frequencySamples'))
        for name,operand in (('amplitudeUnit',r['symbolicAmplitude']),('integralUnit',source)):
            equal(r[name],engine.PHYSICAL_METADATA.dimensions.measure(operand),(key,name))
    inventory=json.loads((previous/'integral-inventory.json').read_text())
    bound_inventory=json.loads((previous/'bound-record-inventory.json').read_text())
    for collection in (inventory,bound_inventory):
        for item in collection:
            src=previous/item['path']
            if src.stat().st_size!=item['bytes'] or digest(src)!=item['sha256']:
                raise ValueError(('saved per-record hash',item['path']))
            dest=base/item['path']; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(src,dest)
    expected_bound_paths={f'bound-records/joins-{ti}.pickle' for ti in range(len(tests))}|{
        f'bound-records/source-{ti}-{si}.pickle' for ti,si in record_map}
    if len(bound_inventory)!=len(expected_bound_paths) or {r['path'] for r in bound_inventory}!=expected_bound_paths:
        raise ValueError('bound artifact census')
    for item in bound_inventory:
        saved=unpickle(previous/item['path']); equal(saved['provenance'],provenance)
        name=Path(item['path']).stem
        expected=([r for r in comparisons if r['test']==int(name.split('-')[1])] if name.startswith('joins-')
                  else record_map[tuple(map(int,name.split('-')[1:]))])
        equal(saved['value'],expected,(name,))
    eval_map={tuple(v[k] for k in ('test','sourceIndex','limitIndex')):v for v in result['evaluations']}
    expected_evals={(ti,si,bi) for ti,si in record_map for bi in range(len(bound['sourceBounds']))}
    if len(eval_map)!=len(result['evaluations']) or set(eval_map)!=expected_evals:
        raise ValueError('integrated source/test/interval census')
    if len(inventory)!=len(expected_evals) or {v['path'] for v in inventory}!={
            'integrals/'+'-'.join(map(str,key))+'.pickle' for key in expected_evals}:
        raise ValueError('integral artifact census')
    for item in inventory:
        saved=unpickle(previous/item['path']); key=tuple(saved[k] for k in ('test','sourceIndex','limitIndex'))
        value=eval_map[key]; equal(saved,value,(key,'integral-packet'))
        r=record_map[key[:2]]; g=value['gaussValues']
        equal(value['frequencies'],r['frequencies'],(key,'frequencies'))
        equal(tuple(value['orders']),tuple(bound['orders']),(key,'orders'))
        limit=bound['sourceBounds'][key[2]]; width=float(tests[key[0]]['width']); ell=bound['profileWidth']
        equal(value['points'],sorted({-limit,limit,0.,*(x for x in (-ell,ell,-width,width) if -limit<x<limit)}),(key,'panels'))
        for name,want in (('sourceResidual',g[-1]-value['originalSourceValues']),
                ('adaptiveResidual',g[-1]-value['adaptiveValues']),('refinements',np.diff(g,axis=0)),
                ('measureMutationResidual',value['measureMutationValues']-g[-1])):
            equal(value[name],want,(key,name))
        for name,values in value.items():
            if isinstance(values,np.ndarray) and not np.all(np.isfinite(values)):
                raise ValueError(('nonfinite saved quadrature',key,name))
        if float(np.max(np.abs(value['measureMutationResidual'])))<=1e-12:
            raise ValueError(('insensitive measure control',key))
        if value['peakWorkspaceBytes']>bound['workspaceBytes']:
            raise ValueError(('phase workspace estimate exceeded',key))
    if len(result['domainChanges'])!=len(record_map): raise ValueError('source interval change census')
    for r,change in zip(bound['records'],result['domainChanges']):
        key=(r['test'],r['sourceIndex']); equal(key,(change['test'],change['sourceIndex']))
        equal(change['difference'],eval_map[(*key,1)]['gaussValues'][-1]-eval_map[(*key,0)]['gaussValues'][-1],(key,'intervalChange'))
    for name,items in (('integral-inventory.json',inventory),('bound-record-inventory.json',bound_inventory)):
        save(base/name,items)
    save(base/'representation-joins.json',representation_records)
    return inventory,bound_inventory,representation_records,replay_inventory


def zero_metadata_transition(old_body, new_body, item):
    """Prove the coefficient-support change of a certified raw zero on reload."""
    comparison=item['comparison']; text=comparison['representationStrings'][2]
    tree=ast.parse(text,mode='eval')
    # coefficients() has an explicit base case for expressions without its
    # three grade generators. Read symbol identities from the saved live AST,
    # before SymPy can canonicalize this particular raw expression to zero.
    symbols={n.args[0].value for n in ast.walk(tree) if isinstance(n,ast.Call)
        and isinstance(n.func,ast.Name) and n.func.id in ('Symbol','Dummy')
        and n.args and isinstance(n.args[0],ast.Constant)}
    generators=engine.PHYSICAL_METADATA.generators
    if symbols & {str(g) for g in generators}:
        raise ValueError('raw live metadata has unresolved grade generators')
    if (text==sp.srepr(sp.S.Zero) or comparison['rawResidual']!=0 or _restore(text)!=0 or
            comparison['normalizedResidual']!=0 or not comparison['limitsEqual'] or
            comparison['certificate'] is None or comparison['certificate']['RESIDUAL']!=0 or
            any(v!=0 for v in (*comparison['replayResiduals'],*comparison['exponentResiduals']))):
        raise ValueError('metadata transition lacks its saved live-to-zero certificate')
    if len(old_body)!=1 or len(new_body)!=1: raise ValueError('raw metadata leaf census')
    old={str(k):v for k,v in old_body[0]}; new={str(k):v for k,v in new_body[0]}
    if set(old)!=set(new) or set(old)!={'PATHS','DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}:
        raise ValueError('raw metadata schema differs')
    equal(old['PATHS'],sp.Tuple(sp.Tuple()),('raw-metadata-scalar-path',))
    equal(old['PATHS'],new['PATHS']);equal(old['DIMENSION_L_T_M'],new['DIMENSION_L_T_M'])
    equal(old['DIMENSION_L_T_M'],engine.cas(item['integrandUnit']))
    # Derive supports from the unchanged coefficients() base case and its
    # homotopy projection. No physical value or unit is replaced here.
    degree=tuple(0 for _ in generators)
    equal(old['MULTIGRADE'],engine.cas((degree,)))
    equal(old['EPSILON_LAMBDA_SUPPORT'],engine.cas(((degree[0],sum(degree[1:])),)))
    equal(new['MULTIGRADE'],engine.cas(tuple(engine.PHYSICAL_METADATA.coefficients(comparison['rawResidual']))))
    equal(new['EPSILON_LAMBDA_SUPPORT'],sp.Tuple())
    return {'liveRepresentationSha256':hashlib.sha256(text.encode()).hexdigest(),
        'liveAstHasGradeGenerators':False,'serializedRawIsZero':True,'certifiedResidualIsZero':True,
        'pathsUnchanged':True,'dimensionLTM':list(map(str,item['integrandUnit'])),
        'liveMultigrade':[list(degree)],'serializedMultigrade':[],
        'liveEpsilonLambdaSupport':[[degree[0],sum(degree[1:])]],'serializedEpsilonLambdaSupport':[]}


def compare_transcripts(previous,base,bound,representation_records):
    old,new=entries(previous/'full.out'),entries(base/'full.out')
    if list(old)!=list(new): raise ValueError('recovery tag/order census changed')
    prefix='PY_S11CD_'+original.PREFIX; final=prefix+'_EMISSION_LINES'; key_tag=prefix+'_WRITE_KEYS'
    indexed_tags=list(new)[:list(new).index(final)]
    index={str(k):v for k,v in _restore(new[final])}
    source_lines=restore_emission_index(index,indexed_tags)
    # The frozen emitter used engine.physical for both write-key tags. All
    # other data-emitter functions and their line locations are unchanged.
    physical_lines,first=inspect.getsourcelines(engine.physical)
    emit_lines=[first+i for i,line in enumerate(physical_lines) if line.lstrip().startswith('emit(')]
    old_lines=dict(source_lines)
    old_lines[key_tag]=emit_lines[0]
    old_lines['PY_S11CD_METADATA_'+original.PREFIX+'_WRITE_KEYS']=emit_lines[1]
    old_index=engine.emission_index(old_lines)
    equal(_restore(old[final]),engine.cas(engine.carrier_fingerprint(engine.cas(old_index))),('old-index-fingerprint',))
    keys={str(k):str(v) for k,v in _restore(new[key_tag])}
    equal(_restore(old[key_tag]),engine.cas(engine.carrier_fingerprint(engine.cas(keys))),('old-keys-fingerprint',))
    allowed={final,key_tag,'PY_S11CD_METADATA_'+original.PREFIX+'_EMISSION_LINES'}
    metadata_transitions=[]
    for item,join in zip(bound['nativeTestIntegralComparisons'],representation_records):
        i=join['occurrence']; tag=prefix+'_BINDING_REPRESENTATION_SHA256_'+str(i)
        equal(_restore(old[tag]),engine.cas(join['representationSha256']),('original-live-representation-hash',i))
        # Any representation-sensitive carrier difference is separately recorded;
        # the original stream remains preserved. Exact restored operand joins,
        # literal native limits and the saved zero certificates are mandatory.
        if not join['liveExactEqual']:
            allowed.update(prefix+'_'+kind+'_'+str(i) for kind in ('BOUND_INTEGRAL_PAIR','RAW_BINDING_RESIDUAL'))
            metadata_tag='PY_S11CD_METADATA_'+original.PREFIX+'_RAW_BINDING_RESIDUAL_'+str(i)
            if old[metadata_tag]!=new[metadata_tag]:
                proof=zero_metadata_transition(_restore(old[metadata_tag]),_restore(new[metadata_tag]),item)
                metadata_transitions.append({'occurrence':i,'tag':metadata_tag,**proof})
                allowed.add(metadata_tag)
    differences=[]
    for tag in old:
        if old[tag]!=new[tag]:
            differences.append({'tag':tag,'sourceSha256':hashlib.sha256(old[tag].encode()).hexdigest(),
                'recoverySha256':hashlib.sha256(new[tag].encode()).hexdigest(),
                'kind':'structural-manifest' if tag in (final,key_tag) or tag.endswith('_EMISSION_LINES')
                       else 'saved-live-to-pickle-representation'})
    save(base/'emission-differences.json',differences)
    unexpected=[v['tag'] for v in differences if v['tag'] not in allowed]
    if unexpected: raise ValueError(('unexpected recovery emission differences',unexpected))
    # Independently test structural decoding rejection of a wrong count and
    # wrong order, using this actual computed index rather than a toy fixture.
    for malformed,tags in ((index|{'s11cdIndexedTagCount':len(indexed_tags)+1},indexed_tags),
                           (index,list(reversed(indexed_tags)))):
        try: restore_emission_index(malformed,tags)
        except ValueError: pass
        else: raise ValueError('emission index mutation was accepted')
    return {'identicalDecodedPayloads':len(old)-len(differences),'changedPayloads':differences,
        'originalIndexFingerprintJoin':True,'originalWriteKeyFingerprintJoin':True,
        'wrongCountAndOrderRejected':True,'certifiedRawZeroMetadataTransitions':metadata_transitions}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--source-run-directory',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True); args=parser.parse_args()
    previous=args.source_run_directory.resolve(); previous.relative_to(STORE)
    base=args.run_directory.resolve(); base.relative_to(STORE); base.mkdir(parents=True,exist_ok=False)
    started=time.monotonic(); repair=json.loads(REPAIR.read_text())
    if str(previous)!=repair['sourceRunDirectory']: raise ValueError('wrong source run')
    original_preflight=json.loads((previous/'preflight.json').read_text())
    checker_name=str(Path(original.__file__).relative_to(ROOT))
    for name,sha in original_preflight['sourceFiles'].items():
        equal(digest(previous/'source'/name),sha,('frozen-source',name))
        if name!=checker_name: equal(digest(ROOT/name),sha,('unchanged-source',name))
    checker_join((previous/'source'/checker_name).read_text(),Path(original.__file__).read_text())
    for name,item in repair['sourceArtifacts'].items():
        equal(digest(previous/name),item['sha256'],('pinned-source-artifact',name))
    pins={str(p.relative_to(ROOT)):digest(p) for p in dict.fromkeys((*SOURCES,*(ROOT/n for n in original_preflight['sourceFiles'])))}
    for name in pins:
        target=base/'source'/name; target.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,target)
    pencil,assembly,factors,tests,numerical,provenance,joins=original.load()
    equal(provenance,original_preflight['provenance'])
    bp=unpickle(previous/'bound-sources.pickle'); qp=unpickle(previous/'quadrature.pickle')
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(qp['dimensionState'])
    for name in ('bound-sources.pickle','quadrature.pickle'):
        shutil.copyfile(previous/name,base/name)
    bound_hash=digest(base/'bound-sources.pickle'); packet_hash=digest(base/'quadrature.pickle')
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'sourceRunDirectory':str(previous),
        'sourcePreflightSha256':digest(previous/'preflight.json'),'wholeCheckerAstJoin':True,
        'sourceArtifacts':repair['sourceArtifacts'],'boundPacketSha256':bound_hash,'quadraturePacketSha256':packet_hash})
    inventory,bound_inventory,representations,replay_inventory=validate_saved(previous,base,bp,qp,pencil,assembly,factors,tests,numerical,provenance)
    original.progress(base,'saved_operands_validated',sources=len(bp['result']['records']),integrals=len(inventory))
    summary=original.finalize(base,started,pins,pencil,factors,provenance,joins,bp['result'],qp['result'],
                              bound_hash,packet_hash,inventory,bound_inventory)
    comparison=compare_transcripts(previous,base,bp['result'],representations)
    for name,item in repair['sourceArtifacts'].items(): equal(digest(previous/name),item['sha256'],('unchanged-original',name))
    equal(digest(base/'bound-sources.pickle'),bound_hash); equal(digest(base/'quadrature.pickle'),packet_hash)
    if engine.PHYSICAL_METADATA.dimensions.constraints: raise ValueError('unresolved dimension constraint')
    summary.update({'sourceRunDirectory':str(previous),'sourcePreflightSha256':digest(previous/'preflight.json'),
        'sourceTranscriptSha256':digest(previous/'full.out'),'wholeCheckerAstJoin':True,
        'savedOperandAndInventoryJoins':True,'quadratureRecomputed':False,'emissionComparison':comparison,
        'restoredExactBindingOccurrences':sum(v['restoredExactEqual'] for v in representations),
        'representationJoinsSha256':digest(base/'representation-joins.json'),
        'emissionDifferencesSha256':digest(base/'emission-differences.json'),
        'sourceBindingReplayArtifacts':replay_inventory,
        'sourceBindingReplayInventorySha256':digest(base/'source-binding-replay-inventory.json')})
    save(base/'checks.json',summary); save(base/'recovery.json',summary)
    equal(pins,{name:digest(ROOT/name) for name in pins},('post-recovery-sources',))
    print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
