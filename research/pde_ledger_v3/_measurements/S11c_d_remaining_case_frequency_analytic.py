#!/usr/bin/env python3
"""Only missing native analytic images and denominator classifications."""
import argparse
import ast
import collections
import copy
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time

import S11c_d_remaining_case_frequency_analytic_inputs as inputs

f,n,sp,engine,q,c = inputs.f,inputs.native,inputs.sp,inputs.engine,inputs.q,inputs.chart
source=inputs.source
CP=f.M/'S11c_d_remaining_case_frequency_analytic_inputs_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_frequency_analytic_plan.md'
SCOPE=('Only missing native scalar analytic images, their derivatives and new denominator identities. '
       'Saved baseline chart/images/root/seed proofs remain inputs. Real Fourier momenta on the '
       'declared scalar disk only; no outgoing-end domain, numerical pencil, pole or physical output.')


def body(fn):
    return hashlib.sha256(ast.dump(ast.parse(textwrap.dedent(inspect.getsource(fn)))).encode()).hexdigest()


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);vr=Path(cp['validation']['runDirectory'])
    f.require(cp['status']=='ACCEPTED_CASE_FREQUENCY_ANALYTIC_INPUTS'
              and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted analytic inputs')
    source.receipts.inspect_guard(vr,'validate')
    f.require(f.digest(vr/'checks.json')==cp['validation']['checksSha256']
              and (vr/'checks.json').read_bytes()==(vr/'validate.stdout').read_bytes(),'accepted final input validation')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],'scope':SCOPE,
              'acceptedAnalyticInputs':{'checkpoint':str(CP),'runDirectory':str(origin),'checksSha256':cp['checksSha256'],
                                       'validatorChecksSha256':cp['validation']['checksSha256']}}
    for name,item in cp['artifacts'].items():source.reference(base,manifest,origin/name,name,item['sha256'])
    for path,name in ((CP,'accepted-analytic-input-checkpoint.json'),(origin/'checks.json','accepted-analytic-input-checks.json'),
                      (origin/'inputs.json','accepted-analytic-input-manifest.json'),(vr/'checks.json','accepted-analytic-input-validation.json')):
        source.reference(base,manifest,path,name,f.digest(path))
    for path in (Path(__file__).resolve(),PLAN,CP):
        name=str(path.relative_to(f.ROOT));value=f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'source pin identity')
        manifest['sourceFiles'][name]=value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,('current source',name))
        if name in cp['sourceFiles']:f.require(f.digest(origin/'source'/name)==value,'accepted frozen source')
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
        f.require(f.digest(dest)==value,'new frozen helper')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(cp['cases'])


class RemoveObservations(ast.NodeTransformer):
    def visit_Expr(self,node):
        if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id in (
                'begin_image','begin_denominator','observe','observe_branch','observe_match'):
            return None
        return self.generic_visit(node)

    def visit_Call(self,node):
        if isinstance(node.func,ast.Name) and node.func.id=='observed_value':
            return self.visit(node.args[-1])
        return self.generic_visit(node)


def source_constructor(recorder,lift):
    original=ast.parse(inspect.getsource(c.source_chart));changed=copy.deepcopy(original)
    loop=next(x for x in changed.body[0].body if isinstance(x,ast.For))
    loop.body.insert(0,ast.parse('begin_image(base,key,old,src,chart)').body[0])
    class Observe(ast.NodeTransformer):
        def visit_Assign(self,node):
            node=self.generic_visit(node)
            if len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
                name=node.targets[0].id
                if name in ('algebraic','analytic'):
                    return [node,ast.parse(f"observe({name!r},{name})").body[0]]
                if name=='actual':return [node,ast.parse("observe('actual-'+name,actual)").body[0]]
            return node
        def visit_Dict(self,node):
            node=self.generic_visit(node)
            for i,key in enumerate(node.keys):
                if isinstance(key,ast.Constant) and key.value in ('firstDerivative','secondDerivative'):
                    node.values[i]=ast.Call(func=ast.Name(id='observed_value',ctx=ast.Load()),
                        args=[ast.Constant(key.value),node.values[i]],keywords=[])
            return node
    changed=Observe().visit(changed)
    reverse=RemoveObservations().visit(copy.deepcopy(changed))
    f.require(ast.dump(reverse)==ast.dump(original),'whole source_chart reverse AST: observations only')
    ns=dict(vars(c),lift=recorder.lift,sp=recorder,begin_image=recorder.image,observe=recorder.observe,observed_value=recorder.value)
    exec(compile(ast.fix_missing_locations(changed),'<native-analytic-records-observed>','exec'),ns)
    return ns['source_chart'],{'wholeSourceChartReverseAST':True,'originalSha256':body(c.source_chart),
                              'changes':'Input, algebraic/analytic, individual derivative and seed/control value persistence only.'}


def lift_constructor(recorder):
    original=ast.parse(inspect.getsource(c.lift));changed=copy.deepcopy(original)
    class Observe(ast.NodeTransformer):
        def visit_Expr(self,node):
            if (isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute)
                    and ast.unparse(node.value.func)=='proofs.append'):
                return [node,ast.parse('observe_branch(proofs[-1])').body[0]]
            return node
    changed=Observe().visit(changed);reverse=RemoveObservations().visit(copy.deepcopy(changed))
    f.require(ast.dump(reverse)==ast.dump(original),'whole native lift reverse AST: branch persistence only')
    ns=dict(vars(c),observe_branch=recorder.branch)
    exec(compile(ast.fix_missing_locations(changed),'<native-analytic-lift-observed>','exec'),ns)
    return ns['lift'],{'wholeLiftReverseAST':True,'originalSha256':body(c.lift),'branchObservers':1}


def denominator_constructor(recorder,lift):
    original=ast.parse(inspect.getsource(c.denominator_chart)).body[0]
    loop=next(x for x in original.body if isinstance(x,ast.For) and ast.unparse(x.target)=='(index, value)')
    changed=copy.deepcopy(loop);changed.body.insert(0,ast.parse('begin_denominator(base,index,value,roots,families)').body[0])
    class Observe(ast.NodeTransformer):
        def visit_Assign(self,node):
            if len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
                name=node.targets[0].id
                if name=='algebraic':return [node,ast.parse("observe('algebraic',algebraic)").body[0]]
                if name in ('quotient','residual'):
                    return [node,ast.parse(f"observe_match(name,{name!r},{name})").body[0]]
            return node
        def visit_Expr(self,node):
            if (isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='f.require'
                    and isinstance(node.value.args[0],ast.Name) and node.value.args[0].id=='matches'):
                return [ast.parse("observe('matches',matches)").body[0],node]
            return node
    changed=Observe().visit(changed);reverse=RemoveObservations().visit(copy.deepcopy(changed))
    f.require(ast.dump(reverse)==ast.dump(loop),'whole native denominator record loop reverse AST')
    args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=x) for x in ('base','original','w','roots','families')],
                       kwonlyargs=[],kw_defaults=[],defaults=[])
    fn=ast.FunctionDef(name='classify',args=args,body=ast.parse('records=[]').body+[changed]+ast.parse('return records').body,decorator_list=[])
    ns=dict(vars(c),lift=recorder.lift,begin_denominator=recorder.denominator,observe=recorder.observe,observe_match=recorder.match)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<native-new-denominator-loop>','exec'),ns)
    family=next(x for x in original.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='families' for t in x.targets))
    rootloop=next(x for x in original.body if isinstance(x,ast.For) and isinstance(x.target,ast.Name) and x.target.id=='root')
    newroot=copy.deepcopy(rootloop)
    at=next(i for i,x in enumerate(newroot.body) if isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='ratio')
    saved=copy.deepcopy(newroot.body[at]);newroot.body[at]=ast.parse("ratio=accepted['coupledRatios'][root['index']]").body[0]
    back=copy.deepcopy(newroot);back.body[at]=saved
    f.require(ast.dump(back)==ast.dump(rootloop),'native family loop with exact saved ratio routing')
    args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=x) for x in ('w','roots','relaxation','relaxation_bound','accepted')],
                       kwonlyargs=[],kw_defaults=[],defaults=[])
    fn=ast.FunctionDef(name='families_from_saved',args=args,body=[copy.deepcopy(family)]+ast.parse('ratios={}').body+[newroot]+ast.parse('return families,ratios').body,decorator_list=[])
    family_ns=dict(vars(c));exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<native-saved-domain-families>','exec'),family_ns)
    return ns['classify'],family_ns['families_from_saved'],{'wholeDenominatorLoopReverseAST':True,
        'originalFunctionSha256':body(c.denominator_chart),'originalLoopSha256':hashlib.sha256(ast.dump(loop).encode()).hexdigest(),
        'familyLoopReverseAST':True,'savedRatioSubstitutions':1,'baselineRootOrDenominatorConstruction':False}


class Recorder:
    def __init__(self):self.active=None;self.branches=0;self.router=None;self.unit=None
    def __getattr__(self,name):return getattr(sp,name)
    def image(self,base,key,old,src,chart):
        self.active=base/'operand-checkpoints'/key;self.active.mkdir(parents=True,exist_ok=False);self.branches=0
        self.unit=tuple(old['unit'])
        f.atomic_pickle(self.active/'input.pickle',{'key':key,'record':old,'variables':{k:src[k] for k in ('frequency','referenceFrequency','origin')},'chart':chart})
    def denominator(self,base,index,value,roots,families):
        self.active=base/'operand-checkpoints'/str(index);self.active.mkdir(parents=True,exist_ok=False);self.branches=0
        self.unit=None
        f.atomic_pickle(self.active/'input.pickle',{'index':index,'original':value,'roots':roots,'families':families})
    def observe(self,slot,value):
        f.require(self.active is not None,'active native operation');f.atomic_pickle(self.active/(slot+'.pickle'),value)
    def value(self,slot,value):self.observe(slot,value);return value
    def branch(self,value):self.observe('branch-'+str(self.branches),value);self.branches+=1
    def match(self,name,slot,value):self.observe(name+'-'+slot,value)
    def diff(self,expression,*variables,**kwargs):return self.router.derivative(expression,*variables,**kwargs)
    def lift(self,expression,w,roots,proofs):return self.router.lift(expression,w,roots,proofs)


class ProofRouter:
    def __init__(self,base,roots,baseline,recorder,native_lift):
        self.base=base;self.roots=roots;self.recorder=recorder;self.seed={};self.certificates={};self.seed_uses=[];self.calls=[]
        self.derivatives={};self.derivative_calls=[];self.lifts={};self.lift_calls=[];self.native_lift=native_lift
        self.original_certificate=engine.BoundedSourceFourierAssembly.reconstruction_certificate
        self.native_certificate_hash=body(self.original_certificate);self.native_seed_hash=body(c.seed_root_identity)
        raw=[]
        for key,record in baseline.items():
            self.add_lift(record['originalLive'],record['unit'],record['algebraic'],record['branchProofs'],
                          {'kind':'accepted-analytic-lift','key':key})
            for order,field in ((1,'firstDerivative'),(2,'secondDerivative')):
                self.add_derivative(record['analytic'],record['unit'],order,record[field],{'kind':'accepted-analytic-derivative','key':key,'field':field})
            for point,pair in record['bindingComparisons'].items():
                for entry in pair['rootIdentities']:
                    raw_path=base/'seed-root-inputs'/str(len(raw));raw_path.mkdir(parents=True,exist_ok=False)
                    f.atomic_pickle(raw_path/'input.pickle',{'sourceKey':key,'point':point,'join':entry,'roots':roots})
                    identity=entry['rootIdentity'];matches=[r for r in roots['roots'].values() if n.same(r['momentum'],identity['momentum'])]
                    f.require(len(matches)==1,'unique saved root momentum address');root=matches[0]
                    args=(root['radicand'],roots['frequency'],root['momentum'],identity['frequency'])
                    raw.append({'sourceKey':key,'point':point,'root':root,'join':entry,'arguments':args})
                    old=self.seed.get(args)
                    if old is not None:f.require(n.same(old['value'],identity),'literal completed seed identity reuse')
                    else:self.seed[args]={'value':identity,'owner':{'key':key,'point':point},'join':entry}
                    f.require(n.same(root['radicand'].subs(roots['frequency'],identity['frequency']),identity['radicand']),
                              'actual native seed-root argument and saved substituted input')
                self.add_pair_certificates(pair,{'kind':'accepted-analytic-comparison','key':key,'point':point})
            for index,proof in enumerate(record['branchProofs']):
                self.add_certificate(proof['certificate'],{'kind':'accepted-analytic-branch','key':key,'index':index})
        f.atomic_pickle(base/'seed-root-input-pairs.pickle',raw)
        frequencies={args[3] for args in self.seed}
        f.require(len(frequencies)==2 and len(self.seed)==len(roots['roots'])*len(frequencies),'complete saved seed/control root identities')
        for root in roots['roots'].values():
            for frequency in frequencies:
                f.require((root['radicand'],roots['frequency'],root['momentum'],frequency) in self.seed,'no missing seed proof')
        f.atomic_pickle(base/'saved-seed-root-atlas.pickle',self.seed)

    def add_derivative(self,expression,unit,order,value,owner):
        key=(expression,tuple(unit),order)
        existing=self.derivatives.get(key)
        if existing is not None:
            f.require(n.same(existing['value'],value),'same native derivative arguments have exact saved values')
        else:self.derivatives[key]={'expression':expression,'unit':unit,'order':order,'value':value,'owner':owner}

    def add_lift(self,expression,unit,value,proofs,owner):
        key=(expression,None if unit is None else tuple(unit))
        if key not in self.lifts:self.lifts[key]={'expression':expression,'unit':unit,'value':value,'proofs':proofs,'owner':owner}

    def derivative(self,expression,*variables,**kwargs):
        f.require(not kwargs and variables in ((self.roots['frequency'],),(self.roots['frequency'],2)),
                  'only exact native first/second frequency derivative call shapes')
        order=1 if len(variables)==1 else variables[1];unit=self.recorder.unit
        target=self.base/'derivative-calls'/str(len(self.derivative_calls));target.mkdir(parents=True,exist_ok=False)
        args={'expression':expression,'variables':variables,'unit':unit,'active':str(self.recorder.active)}
        f.atomic_pickle(target/'input.pickle',args)
        entry=self.derivatives.get((expression,unit,order))
        if entry is None:
            value=sp.diff(expression,*variables)
            entry={'expression':expression,'unit':unit,'order':order,'value':value,
                   'owner':{'kind':'new-native-analytic-derivative','path':str(target)}}
            f.atomic_pickle(target/'value.pickle',entry);self.add_derivative(expression,unit,order,value,entry['owner']);reused=False
        else:
            f.require(n.same((expression,unit,order),(entry['expression'],tuple(entry['unit']),entry['order'])),
                      'exact completed expression/unit/frequency/order derivative route')
            f.atomic_pickle(target/'value.pickle',entry);reused=True
        self.derivative_calls.append({'path':str(target),'reused':reused,'owner':entry['owner']})
        return entry['value']

    def lift(self,expression,w,roots,proofs):
        target=self.base/'lift-calls'/str(len(self.lift_calls));target.mkdir(parents=True,exist_ok=False)
        unit=self.recorder.unit
        f.atomic_pickle(target/'input.pickle',{'expression':expression,'frequency':w,'roots':roots,'unit':unit,
                       'priorProofs':list(proofs),'active':str(self.recorder.active)})
        f.require(not proofs and n.same((w,roots),(self.roots['frequency'],self.roots['roots'])),
                  'exact top-level native lift chart and empty accumulator')
        entry=self.lifts.get((expression,unit))
        if entry is None:
            value=self.native_lift(expression,w,roots,proofs)
            entry={'expression':expression,'unit':unit,'value':value,'proofs':list(proofs),
                   'owner':{'kind':'new-native-lift','path':str(target)}}
            f.atomic_pickle(target/'value.pickle',entry);self.add_lift(expression,unit,value,list(proofs),entry['owner']);reused=False
        else:
            f.require(n.same((expression,unit),(entry['expression'],None if entry['unit'] is None else tuple(entry['unit']))),
                      'exact saved native lift expression/unit/chart route')
            f.atomic_pickle(target/'value.pickle',entry);proofs.extend(entry['proofs']);reused=True
        self.lift_calls.append({'path':str(target),'reused':reused,'owner':entry['owner']})
        return entry['value']

    def add_pair_certificates(self,pair,owner):
        for name in ('certificate','mutation'):
            self.add_certificate(pair.get(name),dict(owner,member=name))

    def add_certificate(self,certificate,owner):
        if certificate is None:return
        f.require(certificate['SHARED_DEFINITIONS']==(),'accepted native shared=False certificate')
        key=(certificate['LEFT'],certificate['RIGHT'],False)
        bucket=self.certificates.setdefault(key,[])
        # Different native dummy identifiers may occur in independent historical
        # certificates. Keep their complete packets; select the first exact input.
        if not bucket:bucket.append({'value':certificate,'owner':owner})

    def root(self,radicand,w,momentum,frequency):
        args=(radicand,w,momentum,frequency);entry=self.seed.get(args)
        target=self.base/'seed-root-calls'/str(len(self.seed_uses));target.mkdir(parents=True,exist_ok=False)
        f.atomic_pickle(target/'input.pickle',{'arguments':args,'active':str(self.recorder.active),'nativeSha256':self.native_seed_hash})
        f.require(entry is not None,'only exact completed native seed identities may be used')
        f.require(n.same(args,(self.roots['roots'][radicand]['radicand'],self.roots['frequency'],entry['value']['momentum'],entry['value']['frequency'])),
                  'full actual seed-root call argument join')
        f.atomic_pickle(target/'value.pickle',entry)
        self.seed_uses.append({'arguments':args,'owner':entry['owner'],'path':str(target),'recomputed':False})
        return entry['value']

    def certificate(self,left,right,*,shared=True):
        target=self.base/'certificate-calls'/str(len(self.calls));target.mkdir(parents=True,exist_ok=False)
        args={'left':left,'right':right,'shared':shared}
        f.atomic_pickle(target/'input.pickle',dict(args,active=str(self.recorder.active),nativeSha256=self.native_certificate_hash))
        f.require(shared is False,'only native chart/binding shared=False certificate calls')
        entry=next((x for x in self.certificates.get((left,right,shared),[]) if n.same((left,right),(x['value']['LEFT'],x['value']['RIGHT']))),None)
        if entry is None:
            value=self.original_certificate(left,right,shared=shared)
            entry={'value':value,'owner':{'kind':'new-native-certificate','path':str(target)}}
            f.atomic_pickle(target/'value.pickle',entry);self.add_certificate(value,entry['owner']);reused=False
        else:f.atomic_pickle(target/'value.pickle',entry);reused=True
        self.calls.append({'path':str(target),'reused':reused,'owner':entry['owner']})
        f.require(n.same((left,right),(entry['value']['LEFT'],entry['value']['RIGHT'])),'exact native certificate input/result route')
        return entry['value']

    def install(self):
        self.recorder.router=self
        c.seed_root_identity=self.root
        engine.BoundedSourceFourierAssembly.reconstruction_certificate=staticmethod(self.certificate)

    def finish(self):
        f.atomic_pickle(self.base/'seed-root-call-routes.pickle',self.seed_uses)
        f.save(self.base/'certificate-call-routes.json',self.calls)
        f.save(self.base/'derivative-call-routes.json',self.derivative_calls)
        f.save(self.base/'lift-call-routes.json',self.lift_calls)
        return {'savedSeedIdentities':len(self.seed),'seedCallsReused':len(self.seed_uses),'newSeedIdentities':0,
                'certificateCalls':len(self.calls),'reusedCertificates':sum(v['reused'] for v in self.calls),
                'newCertificates':sum(not v['reused'] for v in self.calls),
                'newDerivatives':sum(not v['reused'] for v in self.derivative_calls),
                'reusedDerivatives':sum(v['reused'] for v in self.derivative_calls),
                'newLifts':sum(not v['reused'] for v in self.lift_calls),'reusedLifts':sum(v['reused'] for v in self.lift_calls)}


def prohibit():
    def forbidden(*args,**kwargs):raise RuntimeError('completed scientific construction is disabled in missing analytic stage')
    # Keep only the observed native analytic loop/lift and original comparison
    # functions captured in their compiled namespaces, and native symbolic diff.
    for module,names in ((q,('load','sources','resume_sources','end_sources','census','emit_result','main')),
                         (c,('load','root_chart','denominator_chart','source_chart','lift','path_checks','end_tables','emit_result','main')),
                         (source,('load','restore','adapter','construct')),
                         (inputs,('load','prepare')),
                         (n,('load','bind_case','operator_grades')),
                         (n.grades,('load','specs','split','check','term_joins'))):
        for name in names:
            if hasattr(module,name):setattr(module,name,forbidden)
    engine.NumericalReducedAction.__init__=engine.NumericalReducedAction.bind=forbidden
    engine.ChannelInput.__init__=engine.ReducedPencil.__init__=forbidden
    engine.BoundedSourceFourierQuadrature.__init__=engine.BoundedSourceFourierAssembly.construct=forbidden
    n.factors.context=n.factors.cases.source.restore_context=forbidden
    for name in ('solve','inv','pinv','lstsq','svd','eig','eigh','eigvals','eigvalsh','matrix_rank'):setattr(n.np.linalg,name,forbidden)
    n.np.polynomial.legendre.leggauss=forbidden;sp.integrate=sp.lambdify=forbidden


def construct(base,manifest,labels,run,domain_run,family_run,recorder,lift):
    accepted=base/'accepted';roots=f.unpickle(accepted/'accepted-chart/root-chart.pickle')
    den=f.unpickle(accepted/'accepted-chart/denominator-chart.pickle')
    baseline=f.unpickle(accepted/'accepted-chart/analytic-sources.pickle')
    pending=f.unpickle(base/'uncomputed-analytic-image-inputs.pickle')
    missing=f.unpickle(base/'uncomputed-denominator-inputs.pickle')
    f.require(not f.unpickle(base/'uncomputed-radical-inputs.pickle'),'no unproved radical substitution')
    proof_base=base/'analytic-proof-reuse';proof_base.mkdir()
    router=ProofRouter(proof_base,roots,baseline,recorder,lift)
    # Source binding certificates are already completed operands, including
    # their coefficient controls. Retain exact LEFT/RIGHT/native caller joins.
    for label in labels:
        packet=f.unpickle(accepted/'frequency-cases'/label/'frequency-source.pickle')
        for key,record in packet['records'].items():
            f.require(packet['frequency']==roots['frequency'],'actual source derivative variable')
            for order,field in ((1,'firstFrequencyDerivative'),(2,'secondFrequencyDerivative')):
                router.add_derivative(record['liveFrequencyAndGrades'],record['unit'],order,record[field],
                                      {'kind':'accepted-live-source-derivative','case':label,'key':key,'field':field})
            for point,pair in record['bindingComparisons'].items():
                router.add_pair_certificates(pair,{'kind':'accepted-frequency-binding','case':label,'key':key,'point':point})
    for index,record in enumerate(den['records']):
        router.add_lift(record['original'],None,record['algebraic'],record['branchProofs'],{'kind':'accepted-denominator-lift','index':index})
        for i,pair in enumerate(record['branchProofs']):router.add_certificate(pair['certificate'],{'kind':'accepted-denominator-branch','index':index,'branch':i})
    f.atomic_pickle(proof_base/'saved-certificate-atlas.pickle',router.certificates)
    f.atomic_pickle(proof_base/'saved-derivative-atlas.pickle',router.derivatives)
    f.atomic_pickle(proof_base/'saved-lift-atlas.pickle',router.lifts)
    router.install()
    families,ratios=family_run(roots['frequency'],roots['roots'],den['relaxation'],den['relaxationModulusLowerBound'],den)
    f.atomic_pickle(base/'saved-denominator-family-inputs.pickle',{'accepted':den,'roots':roots,'families':families,'ratios':ratios})
    f.require(n.same(ratios,den['coupledRatios']),'exact accepted scalar ratio proofs reused')
    folder=base/'new-denominators';folder.mkdir()
    ordered=sorted(missing,key=sp.default_sort_key)
    f.atomic_pickle(folder/'requested-inputs.pickle',{'ordered':ordered,'physicalUses':missing,'families':families,'chart':roots})
    new_den=domain_run(folder,ordered,roots['frequency'],roots['roots'],families)
    domain_controls=[]
    for index,record in enumerate(new_den):
        changed=record['algebraic']+1;residual=sp.cancel(changed-record['constant']*record['familyExpression'])
        op={'original':record,'changedCoefficientOperand':changed,'actualResidual':residual,'physicalUses':missing[record['original']]}
        f.atomic_pickle(folder/('coefficient-control-'+str(index)+'.pickle'),op)
        f.require(residual!=0 and record['coefficientFrameModulusLowerBound']>0,'new exact denominator coefficient/bound control')
        domain_controls.append({'index':index,'responding':True,'positiveBound':True})
    complete_domains=dict(den,records=den['records']+new_den)
    f.atomic_pickle(base/'remaining-case-denominator-chart.pickle',complete_domains)
    image_owners={('accepted-baseline-analytic',n.BASELINE,key):value for key,value in baseline.items()}
    counts={};total_new=0
    for label in labels:
        packet=f.unpickle(accepted/'frequency-cases'/label/'frequency-source.pickle')
        routes=f.unpickle(base/'analytic-input-cases'/label/'analytic-routes.pickle')
        own=[v for v in pending if v['owner']['case']==label]
        request={v['owner']['key']:v['record'] for v in own}
        folder=base/'new-analytic-records'/label;folder.mkdir(parents=True)
        f.atomic_pickle(folder/'requested-inputs.pickle',own)
        for value in own:
            key=value['owner']['key'];f.require(n.same(value['record'],packet['records'][key]),'accepted exact missing native input')
            f.require(n.same(value['input'],inputs.image_input(packet['records'][key],{k:packet[k] for k in ('frequency','referenceFrequency','origin')})),
                      'full live/seed/control/unit/limit/grade variable routing')
            f.require(value['chartSha256']==f.digest(accepted/'accepted-chart/root-chart.pickle'),'actual accepted chart hash')
        data={'frequencyPacket':{'sources':{**{k:packet[k] for k in ('frequency','referenceFrequency','origin')},'records':request}}}
        new=run(folder,data,roots) if request else {}
        f.require(set(new)==set(request),'all and only missing images')
        for key,value in new.items():image_owners['uncomputed-analytic-image',label,key]=value
        total_new+=len(new);records={};aliases={}
        for key,route in routes.items():
            own=route['analyticOwner'];original=image_owners[own['kind'],own['case'],own['key']];source_record=packet['records'][key]
            f.require(n.same((original['originalLive'],original['unit'],original['limits']),
                             (source_record['liveFrequencyAndGrades'],source_record['unit'],source_record['census']['integralLimits'])),
                      'actual analytic image/source/unit/ordered-limit alias')
            records[key]=dict(original,address=route['address'],analyticOwner=own,analyticOwnerAddress=original['address'],
                              sourceOwner=route['sourceOwner'],sourceOwnerAddress=route['sourceOwnerAddress'])
            aliases[key]=route
        target=base/'analytic-cases'/label;target.mkdir(parents=True)
        result={'records':records,'aliases':aliases,'chart':roots,'denominators':complete_domains,
                'variables':{k:packet[k] for k in ('frequency','referenceFrequency','origin')},
                'generators':packet['generators'],'fieldUnits':packet['fieldUnits'],'equationUnits':packet['equationUnits'],
                'dimensionState':packet['dimensionState'],'acceptedChartDimensionStatePath':str(accepted/'accepted-chart/frequency-chart.pickle'),
                'sourcePacketPath':str(accepted/'frequency-cases'/label/'frequency-source.pickle'),
                'characterPacketPath':str(base/'character-cases'/label/'frequency-characters.pickle'),
                'physicalInputPath':str(base/'analytic-input-cases'/label/'physical-row-term-character-inputs.pickle'),
                'scope':SCOPE,'endFamilyContinuationPending':True,'numericalRowReuseAccepted':False}
        f.atomic_pickle(target/'frequency-analytic.pickle',result)
        f.atomic_pickle(target/'record-aliases.pickle',aliases)
        summary=packet['summary'];counts[label]={'records':len(records),'newAnalyticImages':len(new),
            **{k:summary[k] for k in ('rows','terms','sources')},'aliasedAcceptedOrNewImages':len(records)-len(new)}
        f.save(base/'analytic-case-inventory.json',counts)
        for record in records.values():
            for pair in record['bindingComparisons'].values():
                f.require(pair['normalizedResidual']==0 and all(v==0 for v in pair['proofResiduals']),'all native saved/new analytic seed proofs')
            f.require(all(v['normalizedResidual']==0 and all(z==0 for z in v['proofResiduals']) for v in record['branchProofs']),
                      'all native saved/new branch proofs')
        del packet,routes,request,data,new,records,result
        gc.collect()
    proof_counts=router.finish()
    f.require(total_new==len(pending) and len(new_den)==len(missing),'exact missing image/domain counts')
    f.require(tuple(sum(v[k] for v in counts.values()) for k in ('records','rows','terms','sources'))==(1467,300,647,120),
              'all actual case records/rows/terms/source identities')
    result={'cases':counts,'newAnalyticImages':total_new,'acceptedBaselineAnalyticRecords':len(baseline),
            'newDenominators':len(new_den),'acceptedDenominators':len(den['records']),'newRadicals':0,
            'proofRouting':proof_counts,'denominatorControls':domain_controls,'scope':SCOPE,
            'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets']}
    f.atomic_pickle(base/'remaining-case-frequency-analytic.pickle',result)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    inputs.protect_references(base);manifest,labels=load(base)
    recorder=Recorder();lift,lj=lift_constructor(recorder);run,sj=source_constructor(recorder,lift)
    domain_run,family_run,dj=denominator_constructor(recorder,lift)
    joins={'source':sj,'lift':lj,'denominator':dj,'nativeSeedComparisonSha256':body(c.seed_comparison),
           'nativeBindingComparisonSha256':body(q.binding_comparison),
           'nativeCertificateSha256':body(engine.BoundedSourceFourierAssembly.reconstruction_certificate),
           'nativeSeedRootSha256':body(c.seed_root_identity),
           'exactSavedOperandRouting':{'lift':body(ProofRouter.lift),'derivative':body(ProofRouter.derivative),
               'seedRoot':body(ProofRouter.root),'certificate':body(ProofRouter.certificate),
               'sourceNamespace':body(Recorder),'nativeDerivativeCalls':['sp.diff(analytic, w)','sp.diff(analytic, w, 2)'],
               'scope':'Exact actual expression/unit/chart/frequency/order and full native LEFT/RIGHT/shared=False inputs; every owner packet remains pinned.'}}
    f.save(base/'native-analytic-constructor-joins.json',joins);prohibit()
    result=construct(base,manifest,labels,run,domain_run,family_run,recorder,lift)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'source/current/frozen posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'original input posthash')
    for name,item in manifest['referencedInputs'].items():
        path=base/name;f.require(path.is_symlink() and str(path.readlink())==item['original'] and str(path.resolve())==item['resolvedOriginal']
             and path.stat().st_size==item['bytes'] and f.digest(path)==item['sha256'],'complete reference posthash/address/size')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**{k:v for k,v in result.items() if k not in ('sourceFiles','inputPackets')},
            'status':'COMPLETED_CASE_FREQUENCY_ANALYTIC_SOURCES','nativeJoins':joins,'artifacts':artifacts,
            'newNumericalWork':0,'endFamilyContinuationPending':True,'numericalRowReuseAccepted':False,
            'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
