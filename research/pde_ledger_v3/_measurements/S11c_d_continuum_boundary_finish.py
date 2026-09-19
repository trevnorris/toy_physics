#!/usr/bin/env python3
"""Validate and emit saved continuum boundary/current operands."""
import argparse
import ast
import copy
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import numpy as np
import sympy as sp
import S11c_d_continuum_boundary as c
import S11c_d_continuum_boundary_recover as recovery

f=c.f
SOURCE=f.STORE/'s11c-continuum-boundary-20260919/recovery-01/complete'
PLAN=f.M/'S11c_d_continuum_boundary_finish_plan.md'
VALIDATED=f.STORE/'s11c-continuum-boundary-20260919/finish-01/complete'


def emission_join():
    name=str(Path(c.__file__).resolve().relative_to(f.ROOT))
    before=ast.parse((SOURCE/'source'/name).read_text());after=ast.parse(Path(c.__file__).read_text())
    emitter=next(n for n in after.body if getattr(n,'name',None)=='emit_result')
    initialization=ast.dump(ast.parse('modes.eta,modes.sigma=eta,sigma').body[0])
    matches=[n for n in emitter.body if ast.dump(n)==initialization]
    f.require(len(matches)==1,'exact fingerprint parameter initialization')
    emitter.body.remove(matches[0])
    helpers=[n for n in after.body if getattr(n,'name',None)=='structural_flags'];f.require(len(helpers)==1,'single Boolean metadata adapter')
    after.body.remove(helpers[0]);calls=[]
    class UndoFlags(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n)
            if isinstance(n.func,ast.Name) and n.func.id=='structural_flags':
                calls.append(ast.dump(n));n.func=ast.parse('grades.structural',mode='eval').body
            return n
    UndoFlags().visit(emitter)
    f.require(len(calls)==1 and ast.dump(before)==ast.dump(after),'whole checker join: fingerprint references and census metadata only')
    return {'wholeCheckerReverseAstJoin':True,'initialization':ast.unparse(matches[0]),
            'booleanMetadataAdapter':ast.unparse(helpers[0]),'censusCalls':len(calls),
            'originalSha256':f.digest(SOURCE/'source'/name),'currentSha256':f.digest(Path(c.__file__))}


def load(base):
    joins=emission_join();inputs=json.loads((SOURCE/'inputs.json').read_text())
    result=f.unpickle(SOURCE/'continuum-boundary.pickle')
    f.require(result['sourceFiles']==inputs['sourceFiles'] and result['inputPackets']==inputs['inputPackets'],'completed source manifest identity')
    f.require(f.digest(SOURCE/'continuum-boundary.pickle')=='8e3e53c3b7b06e3e81f961f013b77b4b87dfede081db75b4de458503b062eed0','saved complete boundary packet')
    name=str(Path(c.__file__).resolve().relative_to(f.ROOT))
    for source,sha in inputs['sourceFiles'].items():
        f.require(f.digest(SOURCE/'source'/source)==sha,('frozen source',source))
        if source!=name:f.require(f.digest(f.ROOT/source)==sha,('unchanged current source',source))
    for source,sha in inputs['inputPackets'].items():f.require(f.digest(Path(source))==sha,('completed input',source))
    tables=json.loads((SOURCE/'current-table-inventory.json').read_text());inventory=json.loads((SOURCE/'record-inventory.json').read_text())
    f.require(len(tables)==52 and len(inventory)==14,'complete derivative/end construction inventory')
    old=inputs['savedOperandRecovery']
    f.require(old['wholeCheckerJoin'] and all(inputs['instrumentJoins'][k] for k in ('constructEndReverseAstJoin','currentTablesReverseAstJoin')),'saved exact recovery joins')
    original=Path(old['originalDirectory'])
    for n,item in old['artifacts'].items():
        f.require(f.digest(original/n)==f.digest(SOURCE/n)==item['sha256'],'six original recovered operands unchanged')
    for n,item in {**tables,**inventory}.items():
        f.require(f.digest(SOURCE/n)==item['sha256'] and (SOURCE/n).stat().st_size==item['bytes'],('complete saved artifact',n))
    copied={}
    for p in [*(SOURCE.glob('*.pickle')),*(SOURCE/'current-tables').glob('*.pickle')]:
        relative=p.relative_to(SOURCE);dest=base/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        f.require(f.digest(dest)==f.digest(p),('byte-identical finish packet',str(relative)));copied[str(relative)]={'bytes':p.stat().st_size,'sha256':f.digest(p)}
    pins=dict(result['sourceFiles']);pins[name]=f.digest(Path(c.__file__))
    for p in (Path(__file__),PLAN):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    operands=dict(result['inputPackets']);operands.update({str(SOURCE/n):v['sha256'] for n,v in copied.items()})
    operands[str(SOURCE/'inputs.json')]=f.digest(SOURCE/'inputs.json');operands[str(SOURCE/'full.out')]=f.digest(SOURCE/'full.out')
    for source in pins:
        target=base/'source'/source;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/source,target)
    rp=next(Path(p) for p in operands if Path(p).name=='reduced-action.pickle')
    r,dimensions=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(rp));dimensions.__dict__.update(result['dimensionState'])
    f.require(not dimensions.constraints,'no unresolved source dimensions')
    manifest={'sourceFiles':pins,'inputPackets':operands,'originalInputs':inputs,'emissionRepairJoin':joins,'copiedArtifacts':copied}
    f.save(base/'inputs.json',manifest)
    return r,result,pins,operands,manifest


def validate_operands(result,base):
    tables=json.loads((SOURCE/'current-table-inventory.json').read_text());proofs={}
    for end,v in result['ends'].items():
        saved=f.unpickle(base/(end.lower()+'-boundary.pickle'))
        f.require(recovery.equal(v,saved),'complete combined/end packet identity')
        f.require(len(v['census'])==18 and len(v['clusters'])==5 and sum(x['R'][(0,0)].shape[1] for x in v['clusters'])==7,'complete retained subspaces')
        f.require({(x['ROOT_DISK_INDEX'],x['NORMAL_LIFT_SIGN']) for x in v['census']}=={(i,s) for i in range(9) for s in (-1,1)},'isolated-root/lift candidate census')
        for cluster in v['clusters']:
            i=cluster['info']['INDEX'];mode=f.unpickle(base/(end.lower()+f'-mode-{i}.pickle'))
            f.require(all(recovery.equal(a,cluster[k]) for k,a in mode.items()),'every actual saved mode field')
            coeff=mode['coefficients'];R=mode['R'];K=mode['K'];shift={g:a for g,a in K.items() if g!=(0,0)}
            equations=c.J.equation(coeff,R,shift)
            f.require(recovery.equal(equations[(0,0)],mode['diagnostics']['BASE_EQUATION']),'base equation replay')
            for g in c.J.grades:
                d=mode['diagnostics']['COEFFICIENT_RESIDUALS'][g]
                f.require(recovery.equal(equations[g],d['EQUATION']) and recovery.equal(R[(0,0)].conj().T@R[g],d['GAUGE']),'every coefficient and gauge residual replay')
            f.require(recovery.equal(cluster['D'],{g:1j*a for g,a in c.J.multiply(R,K).items()}),'complete mode normal derivative replay')
        out=[a for a in v['clusters'] if a['info']['direction']=='outgoing'];inc=[a for a in v['clusters'] if a['info']['direction']=='incoming']
        for key,clusters,field in [('outgoing',out,'R'),('outgoingDerivative',out,'D'),('incoming',inc,'R'),('incomingDerivative',inc,'D')]:
            f.require(recovery.equal(v[key],c.concatenate(clusters,field)),('full trace column join',key))
        inverse,res=c.inverse(v['outgoing']);trace=c.J.multiply(v['outgoingDerivative'],inverse)
        f.require(recovery.equal(inverse,v['outgoingInverse']) and recovery.equal(trace,v['trace']) and recovery.equal(res,v['residuals']['inverse']),'inverse and trace coefficients')
        f.require(recovery.equal(c.subtract(v['incomingDerivative'],c.J.multiply(trace,v['incoming'])),v['insertion']),'complete boundary insertion coefficients')
        variables=tuple(next(x for x in v['currentTables']['bulk'][(0,0,0,0)].free_symbols if str(x)==name) for name in
            ('s11cdCurrentLeftMomentum','s11cdCurrentRightMomentum','s11cdAcousticLeftNormalMomentum','s11cdAcousticRightNormalMomentum'))
        offsets=v['offsets'];pairs=0;maxdiff=0.
        for name,table in v['currentTables'].items():
            for index,value in table.items():
                n='current-tables/'+end.lower()+'-'+name+'-'+''.join(map(str,index))+'.pickle'
                f.require(n in tables and f.unpickle(base/n)=={'index':index,'value':value},'exact current derivative entry join')
            functions={g:sp.lambdify(variables,a,'numpy',cse=True) for g,a in table.items()}
            for i,left in enumerate(v['clusters']):
                for j,right in enumerate(v['clusters']):
                    f.require(left['q'].imag>0 and right['q'].imag>0,'current pair convergent depth domain')
                    computed=c.current_pair(table,functions,left,right)
                    for g,a in computed.items():
                        expected=v['currents'][name][g][offsets[i]:offsets[i+1],offsets[j]:offsets[j+1]]
                        maxdiff=max(maxdiff,c.norm(a-expected));f.require(np.array_equal(a,expected),'complete saved current contraction replay')
                    pairs+=1
        for g in c.G:f.require(np.array_equal(v['currents']['total'][g],v['currents']['slab'][g]+v['currents']['bulk'][g]),'complete slab plus bulk current')
        f.require(c.norm(v['residuals']['invariantPair'])<1e-8 and c.norm(v['residuals']['boundary'])<1e-8,'original full subspace and boundary guards')
        f.require(all(c.norm(v['residuals']['currentHermitian'][k])/(1+c.norm(a))<1e-9 for k,a in v['currents'].items()),'original full current Hermitian guards')
        f.require(c.norm(v['residuals']['referenceSignedCurrent'])<1e-8,'original signed current guards')
        f.require(end!='RIGHT' or c.norm(v['etaForcingMutation'])>1e-10,'actual eta mutation response')
        proofs[end]={'clusters':len(v['clusters']),'currentPairsByPart':pairs,'currentReplayMaximum':maxdiff,
                     'residualMaxima':{k:c.norm(a) for k,a in v['residuals'].items()}}
    f.save(base/'operand-validation.json',proofs);return proofs


def finish_tail():
    main=next(n for n in ast.parse(Path(c.__file__).read_text()).body if getattr(n,'name',None)=='main')
    start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and ast.unparse(n.value)=='engine.EMISSION_LINES.clear()')
    body=copy.deepcopy(main.body[start:])
    # Defer stdout until the original partial-transcript and finish guards pass.
    f.require(isinstance(body[-1],ast.Expr) and ast.unparse(body[-1].value)=='print(json.dumps(summary, indent=2))','original summary print tail')
    body[-1]=ast.Return(ast.Name('summary',ast.Load()))
    function=ast.FunctionDef(name='tail',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in
        ('base','result','r','pins','operands','before','started','ends')],kwonlyargs=[],kw_defaults=[],defaults=[]),body=body,decorator_list=[])
    namespace=dict(vars(c));exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),str(Path(__file__)),'exec'),namespace)
    return namespace['tail']


def reuse_validation(result,base,pins,operands):
    old_inputs=json.loads((VALIDATED/'inputs.json').read_text())
    helper=str(Path(__file__).resolve().relative_to(f.ROOT))
    before=next(n for n in ast.parse((VALIDATED/'source'/helper).read_text()).body if getattr(n,'name',None)=='validate_operands')
    after=next(n for n in ast.parse(Path(__file__).read_text()).body if getattr(n,'name',None)=='validate_operands')
    f.require(ast.dump(before)==ast.dump(after),'unchanged completed operand validator')
    f.require(f.digest(VALIDATED/'source'/helper)==old_inputs['sourceFiles'][helper],'frozen completed validator source')
    for name,item in old_inputs['copiedArtifacts'].items():
        f.require(f.digest(VALIDATED/name)==f.digest(base/name)==item['sha256'],'exact completed validation input packets')
    old_proofs=json.loads((VALIDATED/'operand-validation.json').read_text())
    for end,v in result['ends'].items():
        expected={'clusters':len(v['clusters']),'currentPairsByPart':2*len(v['clusters'])**2,'currentReplayMaximum':0.0,
                  'residualMaxima':{k:c.norm(a) for k,a in v['residuals'].items()}}
        f.require(old_proofs[end]==expected,'complete saved validation outcome')
    shutil.copyfile(VALIDATED/'operand-validation.json',base/'operand-validation.json')
    for p in (VALIDATED/'operand-validation.json',VALIDATED/'inputs.json',VALIDATED/'source'/helper):operands[str(p)]=f.digest(p)
    inputs=json.loads((base/'inputs.json').read_text());inputs['inputPackets']=operands
    inputs['completedValidationReuse']={'directory':str(VALIDATED),'validatorAstJoin':True,'packetIdentity':True,
        'validationSha256':f.digest(VALIDATED/'operand-validation.json')};f.save(base/'inputs.json',inputs)
    return old_proofs


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('saved boundary emission/validation budget')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic()
    r,result,pins,operands,manifest=load(base);proofs=reuse_validation(result,base,pins,operands)
    before=f.digest(base/'continuum-boundary.pickle')
    summary=finish_tail()(base,result,r,pins,operands,before,started,result['ends'])
    original={}
    for line in c.grades.decoded_lines(SOURCE/'full.out'):
        tag,_,payload=line.rstrip('\n').partition(': ');original[tag]=c.grades._restore(payload)
    count=0
    for line in c.grades.decoded_lines(base/'full.out'):
        tag,_,payload=line.rstrip('\n').partition(': ')
        if tag in original:f.require(original[tag]==c.grades._restore(payload),'original decoded emission prefix identity');count+=1
    f.require(count==len(original),'all original completed emission tags retained')
    f.require(before==f.digest(SOURCE/'continuum-boundary.pickle')==f.digest(base/'continuum-boundary.pickle'),'unchanged original/finish construction packet')
    for name,item in manifest['copiedArtifacts'].items():f.require(f.digest(base/name)==f.digest(SOURCE/name)==item['sha256'],'all copied post hashes')
    f.save(base/'recovery.json',{'sourceDirectory':str(SOURCE),'emissionRepairJoin':manifest['emissionRepairJoin'],
        'originalPrefixTags':count,'operandValidation':proofs,'copiedArtifacts':manifest['copiedArtifacts'],
        'constructionPacketSha256':before,'scope':'Saved-operand validation and emission only; no mode solves or symbolic derivative construction.'})
    f.require(json.loads((base/'checks.json').read_text())==summary,'final persisted checks identity')
    signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
