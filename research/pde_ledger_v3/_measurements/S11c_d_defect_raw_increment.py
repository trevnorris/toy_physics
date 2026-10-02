#!/usr/bin/env python3
"""Bounded both-face raw direct increment; science only after a pinned pooled gate.

No producer, completed constructor, integral or finite solve is replayed. The
output is an unevaluated nongrazing increment, not a physical loss result.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def require(value, message):
    if value is not True:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def replace_json(path, value):
    tmp = Path(str(path)+'.next')
    save(tmp, value)
    tmp.replace(path)


def function_source(source, name):
    nodes = [n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == name]
    require(len(nodes) == 1, 'unique native function '+name)
    return ast.get_source_segment(source, nodes[0])


def assignment_source(source, function, target):
    fragment = function_source(source, function)
    nodes = [n for n in ast.parse(fragment).body[0].body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == target for t in n.targets)]
    require(len(nodes) == 1, 'unique native assignment '+target)
    return ast.get_source_segment(fragment, nodes[0])


def posthash_records(pins):
    records={}
    for path,expected in pins.items():
        try: records[path]={'expected':expected,'actual':sha(path),'error':None}
        except OSError as error:
            records[path]={'expected':expected,'actual':None,'error':str(error)}
    return records


def expanded_sinh_arguments(expression):
    """Exact argument expansion for one saved-expression equality, not output."""
    mapping={atom:sp.sinh(sp.expand(atom.args[0])) for atom in expression.atoms(sp.sinh)}
    return expression.xreplace(mapping),mapping


def pressure_census(text):
    """Inspect every native constructor; never evaluate full-row science here."""
    tree=ast.parse(text,mode='eval').body
    names={'delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus'}
    def hits(node):
        found=[]
        for n in ast.walk(node):
            if (isinstance(n,ast.Call) and isinstance(n.func,ast.Name)
                and n.func.id in ('Symbol','Function') and n.args
                and isinstance(n.args[0],ast.Constant) and isinstance(n.args[0].value,str)
                and 'delta_p' in n.args[0].value):
                require(n.func.id=='Symbol' and n.args[0].value in names,'unknown native pressure slot')
                found.append(n.args[0].value)
        return found
    def degree(node):
        if not hits(node):return 0
        require(isinstance(node,ast.Call) and isinstance(node.func,ast.Name),'pressure constructor')
        if node.func.id=='Symbol':return 1
        if node.func.id=='Add':return max(map(degree,node.args))
        if node.func.id=='Mul':return sum(map(degree,node.args))
        if node.func.id=='Pow':
            power=node.args[1]
            require(isinstance(power,ast.Call) and isinstance(power.func,ast.Name)
                    and power.func.id=='Integer' and len(power.args)==1,'pressure integer power')
            n=ast.literal_eval(power.args[0]);require(type(n) is int and n>=0,'pressure polynomial')
            return n*degree(node.args[0])
        raise ValueError('nonaffine/unknown pressure operation')
    children=tree.args if isinstance(tree,ast.Call) and isinstance(tree.func,ast.Name) and tree.func.id=='Add' else [tree]
    selected=[];covered=[]
    for index,node in enumerate(children):
        occurrences=hits(node);covered.extend(occurrences)
        if occurrences:
            value=ast.get_source_segment(text,node)
            selected.append({'childIndex':index,'constructorText':value,
                'constructorSha256':hashlib.sha256(value.encode()).hexdigest(),
                'occurrences':occurrences,'pressureDegree':degree(node)})
    all_hits=hits(tree)
    require(sorted(covered)==sorted(all_hits),'complete native slot coverage')
    return {'totalChildren':len(children),'selected':selected,'occurrenceCounts':{n:all_hits.count(n) for n in sorted(names)},
            'maximumPressureDegree':max([0]+[r['pressureDegree'] for r in selected])}


def fourier_contract(source):
    """Exact executable AST contract, not a convention inferred from comments."""
    def same(node,text):return ast.dump(node,include_attributes=False)==ast.dump(ast.parse(text,mode='eval').body,include_attributes=False)
    def assignment(function,name):
        tree=ast.parse(function_source(source,function))
        found=[n.value for n in ast.walk(tree) if isinstance(n,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id==name for t in n.targets)]
        require(len(found)==1,'unique measure assignment '+name);return found[0]
    definition=assignment('profile_bindings','definition')
    require(isinstance(definition,ast.Call) and same(definition.func,'sp.Integral'),'native profile integral')
    require(same(definition.args[0],'phase*local_field/(2*sp.pi)**3'),'native normalized profile forward')
    require(len(definition.args)==2 and isinstance(definition.args[1],ast.Starred),'profile starred limit shape')
    require(ast.unparse(definition.args[1])=='*((y, -sp.oo, sp.oo) for y in Y)','native profile Y limits')
    records={}
    for name,numerator in [('p0','phase0*diagonal*local_source'),('p1','phase1*off_diagonal*local_source'),('p2','phase1*second*local_source')]:
        node=assignment('kernel_apply',name)
        require(isinstance(node,ast.Call) and same(node.func,'integral'),'native source apply integral')
        require(len(node.args)==(3 if name=='p2' else 2) and not node.keywords,'native integration arity')
        require(same(node.args[0],numerator+'/(2*sp.pi)**3'),'native inverse normalization '+name)
        records[name]=ast.unparse(node)
    require(ast.unparse(assignment('kernel_apply','p2').args[1])=='*limits1','native second source limits')
    require(ast.unparse(assignment('kernel_apply','p2').args[2])=='*((v, -sp.oo, sp.oo) for v in MIDDLE)','plain second middle measure')
    require(same(assignment('kernel_apply','local_source'),'inputs.at_source(source)'),'native unnormalized source forward')
    require(same(assignment('kernel_apply','phase1'),'sp.exp(sp.I*(sum(k*x for k,x in zip(kout,X))-sum(k*y for k,y in zip(kin,Y))))'),'native source phase sign')
    require(same(assignment('profile_bindings','phase'),'sp.exp(-sp.I*sum((a-b)*y for a,b,y in zip(ko,ki,Y)))'),'native profile phase sign')
    globals_=ast.parse(source).body
    declarations={name:next(n.value for n in globals_ if isinstance(n,ast.Assign) and any(isinstance(a,ast.Name) and a.id==name for a in n.targets)) for name in ('Y','MIDDLE')}
    counts={}
    for name,node in declarations.items():
        generator=node.args[0];require(isinstance(generator,ast.GeneratorExp),'native coordinate generator')
        iterator=generator.generators[0].iter
        require(same(iterator,'range(1,4)'),'native three-dimensional coordinate count')
        start,stop=map(ast.literal_eval,iterator.args);counts[name]=stop-start
    require(counts['Y']==counts['MIDDLE'],'same source/middle dimensions')
    dimensions=[n.value for n in ast.walk(ast.parse(function_source(source,'fourier_profiles'))) if isinstance(n,ast.Assign)
                and any(isinstance(a,ast.Subscript) and isinstance(a.value,ast.Name) and a.value.id=='NEW_DIMENSIONS' for a in n.targets)]
    require(len(dimensions)==1 and ast.literal_eval(dimensions[0])==(3,0,0),'native profile Fourier dimensions')
    profile_tree=ast.parse(function_source(source,'fourier_profiles'))
    spectral=[n for n in ast.walk(profile_tree) if isinstance(n,ast.Assign)
              and any(isinstance(a,ast.Subscript) and isinstance(a.value,ast.Name)
                      and a.value.id=='spectral' for a in n.targets)]
    require(len(spectral)==1 and same(spectral[0].value,'function(*(a-b for a,b in zip(kout,kin)))'),'native profile transfer arguments')
    returned=[n.value for n in ast.walk(profile_tree) if isinstance(n,ast.Return)]
    require(len(returned)==1 and same(returned[0],'expression.xreplace(spectral)'),'native profile transfer insertion')
    return {'profileForwardPower':-3,'sourceInversePower':-3,'sourceForwardPower':0,
        'coordinates':counts['Y'],'invariantEdgeCoordinates':counts['Y']-1,'nativeProfileDimension':list(ast.literal_eval(dimensions[0])),
        'profileDefinition':ast.unparse(definition),'nativeSourceIntegrals':records,
        'sourcePhase':ast.unparse(assignment('kernel_apply','phase1')),
        'profilePhase':ast.unparse(assignment('profile_bindings','phase')),
        'middleMeasure':'plain three-dimensional; one-dimensional after both edge deltas',
        'scope':'Syntactic executable-factor contract; no native integral called or evaluated.'}


def containment():
    memory = 4*1024**3
    group = next(v[3:] for v in Path('/proc/self/cgroup').read_text().splitlines() if v.startswith('0::'))
    base = Path('/sys/fs/cgroup')/group.lstrip('/')
    actual = {k:(base/k).read_text().strip() for k in ('memory.max','memory.swap.max','pids.max')}
    actual.update(cgroup=str(base), affinity=sorted(os.sched_getaffinity(0)),
                  threads={k:os.environ.get(k) for k in THREADS})
    require(actual['memory.max']==str(memory) and actual['memory.swap.max']=='0'
            and actual['pids.max']=='32' and len(actual['affinity'])==1
            and all(v=='1' for v in actual['threads'].values()), 'actual pooled limits')
    require('S11C_POOLED_GUARD_MANIFEST' in os.environ, 'pooled guard required')
    require(resource.getrlimit(resource.RLIMIT_CPU)==(resource.RLIM_INFINITY,)*2,'no CPU deadline')
    resource.setrlimit(resource.RLIMIT_AS,(memory,memory))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    return actual | dict(nativeAddressSpace=memory,durationLimits=None)


def verify_gate(path, manifest_path, manifest):
    g=json.loads(Path(path).read_text())
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'gate worker/input')
    for p,h in g['sourcePins'].items(): require(sha(p)==h,'gate pin '+p)
    r=json.loads(Path(g['buildReviewRecord']).read_text())
    require(sha(g['buildReviewRecord'])==g['buildReviewRecordSha256'],'actual build record')
    if r['independentBuildClearance'] is True:
        for k in ('workerSha256','manifestSha256'):
            require(r[k]==g[k],'review/gate '+k)
    else:
        require(g['independentBuildClearance'] is False and g['localToolingRepairAccepted'] is True,'literal review status')
        require(sha(g['repairRecord'])==g['repairRecordSha256'],'local repair pin')
        repair=json.loads(Path(g['repairRecord']).read_text())
        require(repair['toolingOnly'] is True and repair['testsPassed'] is True,'tested representation repair')
        require(repair['workerSha256']==g['workerSha256'] and repair['reviewedWorkerSha256']==r['workerSha256'],'reviewed/repaired worker join')
        require(repair['reviewRecordSha256']==g['buildReviewRecordSha256'] and r['allChecksPassed'] is True,'actual assessed record')
        require(r['reports']['claude']['literalVerdict']=='CLEAR FOR THIS BOUNDED RAW-INCREMENT BUILD'
                and r['reports']['grok']['literalVerdict']=='NEEDS REVISION','preserved literal reviews')
        require(sha(g['executionAuthority'])==g['executionAuthoritySha256'],'execution authority pin')
        authority=json.loads(Path(g['executionAuthority']).read_text())
        require(authority['scienceExecutionsAuthorized']==1 and authority['localToolingRepairAllowed'] is True,'standing scope and tooling authority')
    for k in ('guardSha256','supervisorSha256'):
        require(r[k]==g[k],'review/gate '+k)
    require(sha(g['sharedGuard'])==g['guardSha256'] and sha(g['supervisor'])==g['supervisorSha256'],'actual helpers')
    require(g['methodRecordSha256']==sha(manifest['methodRecord']),'method record pin')
    require(json.loads(Path(manifest['methodRecord']).read_text())['jointIndependentMethodClearance'] is True,'method assessed')
    require(g['scope']==manifest['scope'] and g['scientificRunsAuthorized']==1,'bounded gate scope')
    return g


class Journal:
    def __init__(self,out):
        self.out=out; self.active=None; self.completed=[]; self.artifacts={}

    def encode(self,value):
        if isinstance(value,(sp.Basic,sp.MatrixBase)):
            return {'text':str(value),'srepr':sp.srepr(value)}
        if isinstance(value,dict):
            require(all(isinstance(k,str) for k in value),'string evidence keys')
            return {k:self.encode(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)): return [self.encode(v) for v in value]
        return value

    def emit(self,name,value):
        p=self.out/(name+'.json');save(p,self.encode(value))
        r={'path':p.name,'sha256':sha(p),'bytes':p.stat().st_size}
        self.artifacts[p.name]=r
        replace_json(self.out/'artifact-index.json',self.artifacts)
        return r

    def stage(self,name,inputs,fn):
        previous=self.active;self.active=name;operand=self.emit(name+'-input',inputs)
        value=fn();result=self.emit(name+'-return',value)
        self.completed.append({'name':name,'input':operand,'result':result})
        replace_json(self.out/'operation-index.json',self.completed)
        self.active=previous;return value

    def zero(self,name,left,right):
        previous=self.active;self.active=name
        self.emit(name+'-input',{'left':left,'right':right})
        raw=left-right;self.emit(name+'-raw',{'residual':raw})
        residual=sp.cancel(sp.together(raw));self.emit(name+'-return',{'raw':raw,'cancelled':residual})
        require(residual==0,name)
        self.active=previous

    def nonzero(self,name,baseline,corrupt):
        previous=self.active;self.active=name
        self.emit(name+'-input',{'baseline':baseline,'corrupt':corrupt})
        movement=sp.simplify(corrupt-baseline)
        self.emit(name+'-return',{'movement':movement,'zero':movement.is_zero,'finite':movement.is_finite})
        require(movement.is_zero is False and movement.is_finite is True,name)
        self.active=previous

    def sinh_zero(self,name,left,right):
        previous=self.active;self.active=name
        self.emit(name+'-original-input',{'left':left,'right':right})
        self.emit(name+'-original-raw',{'residual':left-right})
        new_left,left_map=expanded_sinh_arguments(left)
        new_right,right_map=expanded_sinh_arguments(right)
        self.emit(name+'-argument-expansion',{'leftReplacements':[[a,b] for a,b in left_map.items()],
            'rightReplacements':[[a,b] for a,b in right_map.items()],
            'left':new_left,'right':new_right,'identity':'exact expansion inside sinh; intrinsic odd symmetry; no numeric tolerance'})
        self.zero(name+'-canonical',new_left,new_right)
        self.active=previous


def decode(value):
    if isinstance(value,dict):
        if set(value)=={'text','srepr'}:
            return sp.sympify(value['srepr'],locals={'Str':Str,
                'Equality':lambda *a:sp.Eq(*a,evaluate=False),
                'StrictGreaterThan':lambda *a:sp.StrictGreaterThan(*a,evaluate=False)})
        return {k:decode(v) for k,v in value.items()}
    if isinstance(value,list):return [decode(v) for v in value]
    return value


def named(value,name):
    found=[v for k,v in value if str(k)==name]
    require(len(found)==1,'unique named operand '+name)
    return found[0]


def one_symbol(objects,name):
    found=set().union(*(v.atoms(sp.Symbol) for v in objects))
    selected=[s for s in found if s.name==name]
    require(len(selected)==1,'unique source symbol '+name)
    return selected[0]


def scientific_work(manifest,J):
    # Byte copies, not reprinting, preserve all consumed observations.
    copied={};saved={};prior=J.out/'saved-operands';prior.mkdir()
    for alias,record in manifest['savedOperands'].items():
        p=Path(record['path']);require(sha(p)==record['sha256'],'saved operand '+alias)
        target=prior/(alias+'.json');shutil.copyfile(p,target)
        require(sha(target)==record['sha256'],'byte copy '+alias)
        copied[alias]={'source':str(p),'copy':str(target.relative_to(J.out)),'sha256':record['sha256'],'functionCalled':False}
        saved[alias]=decode(json.loads(target.read_text()))
    J.emit('saved-return-reuse',copied)
    native=json.loads(Path(manifest['nativeRecords']).read_text())
    for p,r in native['sourcePins'].items():require(sha(p)==r['sha256'],'native source '+p)
    c1=Path(manifest['c1Source']).read_text();c2=Path(manifest['c2Source']).read_text()
    geom=Path(manifest['geometrySource']).read_text()
    physical=json.loads(Path(manifest['physicalInput']).read_text())
    require(saved['binding']['physicalInput']==physical,'saved physical-input join')
    require(saved['point']['originalParameters']==physical['parameters'],'saved bare parameters')
    values={k:sp.sympify(v) for k,v in physical['parameters'].items()};values['omega']=sp.Integer(3)
    B=saved['bare']; mixed=B['mixed']; S0=saved['source_plus']['source00'];SM=saved['source_minus']['source00']
    qnames=('q_i','q_h','q_s','q_o');qi,qh,qs,qo=(one_symbol([mixed],n) for n in qnames)
    k=one_symbol([mixed],'k');H=one_symbol([mixed],'H')
    original_modes=saved['modes']['modes']
    S=one_symbol([row[2] for row in original_modes],'S')
    om=one_symbol([mixed],'omega');rho=one_symbol([mixed],'rho_m')
    h=saved['modes']['h'];s=saved['modes']['s']
    eta,sigma=saved['binding']['independentGrades'];eps=saved['binding']['epsilon']
    t,Q=sp.symbols('increment_transfer increment_difference',real=True)
    cs=sp.Symbol('increment_effective_bulk_speed',positive=True)
    edge=tuple(saved['point']['edge']);freq=values['omega'];mass=values['rho_m']
    J.zero('frequency-join',freq,saved['point']['omega'])
    require(edge==(sp.Rational(1,5),sp.Rational(1,10)),'saved edge')
    numeric={om:freq,rho:mass}
    # cs does not enter the restored chemical/source/consumer operands.
    independent_records=[native['chemicalSource']['valueConstructorText']]
    independent_records += [r['valueConstructorText'] for r in native['geometry']['face_velocity']['cases']]
    independent_records += [r['valueConstructorText'] for r in native['faceResponseSources']['cases']]
    restored_raw=[saved['input_'+label]['raw'] for label in ('plus','minus')]
    restored_raw += [saved['consumer_'+row]['raw'] for row in manifest['rowNames']]
    J.emit('effective-speed-reuse-domain',{'sourceSymbolsCs':any('c_s0' in x for x in independent_records),
        'restoredRawSymbols':sorted({a.name for v in restored_raw for a in v.free_symbols}),
        'savedFrequency':saved['binding']['frequency'],'newFrequency':freq,
        'scope':'These source/consumer/response pieces must be independent of cs; new depths carry effective cs.'})
    require(not any('c_s0' in x for x in independent_records),'source cs independence')
    require(not any(a.name=='c_s0' for v in restored_raw for a in v.free_symbols),'consumer cs independence')
    require(not any(a.name=='c_s0' for v in (S0,SM) for a in v.free_symbols),'restored zero source cs independence')
    def q(p):
        rad=freq**2/cs**2-sum(x*x for x in edge)-p*p
        return sp.Piecewise((sp.sqrt(rad),rad>0),(sp.I*sp.sqrt(-rad),rad<0),(sp.S.NaN,True))
    physical_depth_bindings={qi:q(k),qh:q(k+t),qs:q(k+Q-t),qo:q(k+Q)}
    J.emit('physical-depth-bindings',{'bindings':{str(a):v for a,v in physical_depth_bindings.items()},'outgoingPrescription':'same positive root; zeros are excluded from raw evaluation'})
    J.emit('physical-sheet',{'frequency':freq,'cs':cs,'edge':edge,'input':q(k),'heightRoute':q(k+t),
        'slopeRoute':q(k+Q-t),'output':q(k+Q),'zeroRule':'not pointwise evaluated; limiting omega+i0+ prescription remains uncomputed',
        'externalDomain':[sp.Ne(physical_depth_bindings[x],0,evaluate=False) for x in (qi,qo)],'reflectedRoute':'k+k_out-m; separate from native MIDDLE_Q'})
    qpoint={qi:saved['point']['inputDepth'],qo:saved['point']['outputDepth']}
    for name,p,v in [('input',0,qpoint[qi]),('output',sp.Rational(1,10),qpoint[qo])]:
        J.zero('saved-sheet-'+name,q(sp.sympify(p)).subs(cs,values['c_s0']),v)
    # Reuse the upper boundary result; derive only the missing lower boundary.
    modes=original_modes
    lower_normal_source=assignment_source(geom,'build_eulerian_face_source','normal_exact') if 'def build_eulerian_face_source(' in geom else None
    if lower_normal_source is None:
        nodes=[n for n in ast.walk(ast.parse(geom)) if isinstance(n,ast.Assign)
               and any(isinstance(x,ast.Name) and x.id=='normal_exact' for x in n.targets)
               and 'grad_h' in ast.get_source_segment(geom,n)]
        require(len(nodes)==1,'native graph normal address');lower_normal_source=ast.get_source_segment(geom,nodes[0])
    # Join actual native lower geometry before constructing its new boundary.
    native_normals={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str})
                    for r in native['geometry']['face_normal']['cases']}
    native_traces={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str})
                   for r in native['geometry']['face_shift']['cases']}
    lower_trace=native_traces[-1][0]/eps
    lower_jet=one_symbol([lower_trace],'d_w_delta_p_minus')
    width=one_symbol([lower_trace],'W_0');profile=one_symbol([lower_trace],'w1_profile')
    lab_height=sp.diff(lower_trace,lower_jet)
    J.zero('native-lower-lab-height',lab_height,-width*eta*profile/2)
    normal0=native_normals[-1][0]
    profile_jet=one_symbol([normal0],'w1_profile_d1')
    nc={'__builtins__':{},'tuple':tuple,'face':-1,'grad_h':(-s,sp.S.Zero,sp.S.Zero),'denominator':sp.sqrt(1+s*s)}
    exec(lower_normal_source,nc)
    normal_derivative=sp.diff(nc['normal_exact'][0],s).subs(s,0)
    J.emit('native-lower-geometry-operands',{'trace':lower_trace,'labHeight':lab_height,
        'nativeNormal':normal0,'graphNormal':nc['normal_exact'],'outwardSlopeDefinition':sigma*profile_jet/2})
    J.zero('native-lower-normal-slope',sp.diff(normal0[0],sigma),normal_derivative*profile_jet/2)
    J.zero('native-lower-normal-orientation',normal0[3],sp.Integer(-1))
    def lower_boundary(label,height_face=-1,slope_multiplier=1):
        face=-1;normal_context={'grad_h':(face*s,sp.S.Zero,sp.S.Zero),'face':face,'denominator':sp.sqrt(1+s*s)}
        normal_context.update(__builtins__={},tuple=tuple)
        exec(lower_normal_source,normal_context)
        normal=sp.Matrix(normal_context['normal_exact']);velocity=sp.S.Zero;pressure=sp.S.Zero
        for a,b,momentum,depth,amplitude in modes:
            for n in range(2-a):
                shifted=h**(a+n)*s**b*amplitude*(sp.I*face*depth*height_face)**n/sp.factorial(n)
                velocity+=face*sp.I*face*depth*shifted
                pressure+=sp.I*om*rho*shifted
                if b==0:velocity+=slope_multiplier*sp.diff(normal[0],s).subs(s,0)*s*sp.I*momentum*shifted
        equations=[sp.expand(velocity-1).coeff(h,a).coeff(s,b) for a,b,_,_,_ in modes]
        J.emit(label+'-boundary-operands',{'normalSource':lower_normal_source,'normal':normal,'reference':face*values['W_0']/2,
            'extensionSign':face,'labHeightSign':height_face,'modes':modes,'velocity':velocity,'pressure':pressure,'equations':equations})
        solution={}
        for i,(eq,row) in enumerate(zip(equations,modes)):
            amp=row[-1];sub=eq.subs(solution);linear=sp.diff(sub,amp);constant=sub.subs(amp,0)
            J.emit(label+'-amplitude-%d-input'%i,{'equation':sub,'amplitude':amp,'linear':linear,'constant':constant})
            require(sp.diff(linear,amp)==0 and linear!=0,'triangular boundary coefficient')
            value=sp.cancel(-constant/linear);solution[amp]=value
            J.emit(label+'-amplitude-%d-return'%i,{'value':value,'denominator':linear})
            J.zero(label+'-equation-%d'%i,eq.subs(solution),sp.S.Zero)
        coefficients=[sp.cancel(sp.expand(pressure.subs(solution)).coeff(h,a).coeff(s,b)) for a,b,_,_,_ in modes]
        J.emit(label+'-coefficients',{'coefficients':coefficients,'solution':{str(a):v for a,v in solution.items()}})
        return coefficients
    lower=J.stage('lower-boundary',{'upperSaved':mixed,'nativeNormal':lower_normal_source,'face':-1},lambda:lower_boundary('lower'))
    J.zero('lower-upper-outward-mirror',lower[3],mixed)
    # New lower linear joins against native c1, not re-running the upper derivation.
    kin=sp.symbols('increment_native_kin1:4',real=True)
    ctx=SimpleNamespace(omega=om,rho_m=rho,c_s0=cs)
    ns={'sp':sp,'Mapping':dict,'Inputs':object,'k_in':kin}
    exec(function_source(c1,'dtn_first_kernel'),ns)
    linear_native=ns['dtn_first_kernel']
    nh=linear_native(ctx,{'height_hat':1,'tilt_hat':(0,0,0)},qh,qi).subs(dict(zip(kin,(k,*edge))))
    # Same physical dispersion; no independent-depth equality is assumed.
    nh=sp.expand(nh).subs(cs**(-2),(k*k+sum(e*e for e in edge)+qi*qi)/om**2)
    J.zero('lower-native-height',lower[1],nh)
    nsl=linear_native(ctx,{'height_hat':0,'tilt_hat':(1,0,0)},qs,qi).subs(dict(zip(kin,(k,*edge))))
    J.zero('lower-native-slope',lower[2],nsl)
    J.zero('lower-native-flat',lower[0],rho*om/qi)
    # Complete numerator factorization modulo the two declared dispersion identities.
    def factorization():
        Br=-sp.I*om*rho*(qi/(qh*qo)+k*(H+2*k)/(qh*qo*(qh+qi))
              +k*(2*k+2*Q-H)/(qi*qo*(qs+qo)))
        delta=sp.together(mixed-H*Br);num,den=sp.fraction(delta)
        rules={qh**2:qi**2-H*(H+2*k),qs**2:qo**2+H*(2*k+2*Q-H)}
        J.emit('factorization-before-reduction',{'original':mixed,'candidateB':Br,'rawDifference':delta,'numerator':num,'denominator':den,
            'identities':{str(a):v for a,v in rules.items()},'domain':'qi qo qh (qh+qi) (qs+qo) nonzero; common positive physical sheet'})
        expanded=sp.expand(num)
        for term in sp.Add.make_args(expanded):
            powers=term.as_powers_dict()
            require(all(powers.get(x,0) in (0,1,2) for x in (qh,qs)),'quadratic dispersion reduction domain')
        reduced=sp.expand(expanded.subs(rules,simultaneous=True))
        J.zero('factorization-physical-numerator',reduced,sp.S.Zero)
        contact=mixed.subs({H:0,qh:qi,qs:qo},simultaneous=True)
        J.zero('upper-complete-contact',contact,sp.S.Zero)
        J.zero('lower-complete-contact',lower[3].subs({H:0,qh:qi,qs:qo},simultaneous=True),sp.S.Zero)
        return {'B':Br,'contact':contact,'lowerContact':0,'exactResidual':reduced,'domainDenominator':den}
    fact=J.stage('physical-factorization',{'coefficient':mixed,'momenta':[k,k+H,k+Q-H,k+Q]},factorization)
    L=values['L_W'];W=values['W_0'];A=lambda x:L*x/(4*sp.sinh(sp.pi*L*x/2))
    # Restore the completed transform; establish symbol/factor joins, not a new transform integral.
    transform=saved['transform'];old_t=transform['t'];old_L=transform['L']
    J.zero('saved-transform-A',transform['A'].subs({old_t:t,old_L:L},simultaneous=True),A(t))
    J.zero('saved-transform-jet',transform['jetTransform'].subs({old_t:Q-t,old_L:L},simultaneous=True),L*A(Q-t))
    J.emit('raw-ordered-before-cancel',{'heightScale':W/2,'B':fact['B'].subs(H,t),'heightNumerator':A(t),'jet':L*A(Q-t)/2,'phaseDivisor':sp.I,'numeric':{str(a):v for a,v in numeric.items()}})
    raw=sp.cancel((W/2)*fact['B'].subs(H,t)*A(t)*(L*A(Q-t)/2)/sp.I).subs(numeric)
    raw_lower=raw # only after saved lower==upper exact mirror identity above
    J.emit('raw-ordered-kernel',{'heightTransfer':t,'slopeTransfer':Q-t,'externalInput':k,'externalOutput':k+Q,
        'rawPlus':raw,'rawMinus':raw_lower,'physicalRaw':raw.xreplace(physical_depth_bindings),'contactPlus':fact['contact'],'contactMinus':fact['lowerContact'],
        'uncancelledHeight':'(W0/2)*(delta(t)/2+PV[A(t)/(i*t)])','jet':L*A(Q-t)/2,
        'PVScope':'ordinary factored density only on nonzero external-depth domain; collision distribution not computed',
        'integralEvaluated':False,'nativeSlot':'raw z_three[0,2] coefficient; eta*sigma applied once',
        'sourceInverseMeasure':'(2*pi)^-3 unchanged','middleJacobian':1,'edgeDeltaCountAfterReduction':2})
    old_action=saved['action'];old_depth=one_symbol([old_action['integrand']],'q_t')
    old_transfer=one_symbol([old_action['integrand']],'t')
    selected=raw.subs({k:0,Q:sp.Rational(1,10),qi:qpoint[qi],qo:qpoint[qo],qh:old_depth,t:old_transfer},simultaneous=True)
    J.sinh_zero('saved-selected-integrand',selected,old_action['integrand'])
    dims=saved['dimensions'];J.emit('inherited-kernel-dimensions',dims)
    # L,T,M integer dimensions from original declarations, not a numerical bound.
    reduced=list(dims['reducedOneDimensionalKernel']);full=[list(x) for x in dims['nativeKernel']]
    schema_node=next(n for n in ast.parse(c2).body if isinstance(n,ast.Assign) and any(isinstance(a,ast.Name) and a.id=='DIMENSION_SCHEMA' for a in n.targets))
    schema=ast.literal_eval(schema_node.value)
    symbol_dims={rho:tuple(schema['rho_m']),om:tuple(schema['omega']),k:(-1,0,0),H:(-1,0,0),qi:(-1,0,0),qo:(-1,0,0),qh:(-1,0,0),qs:(-1,0,0)}
    def dim(v):
        if v.is_number:return (0,0,0)
        if v in symbol_dims:return symbol_dims[v]
        if v.is_Add:
            ds=[dim(x) for x in v.args];require(all(x==ds[0] for x in ds),'homogeneous new coefficient');return ds[0]
        if v.is_Mul:return tuple(sum(dim(x)[i] for x in v.args) for i in range(3))
        if v.is_Pow and v.exp.is_Rational:return tuple(x*v.exp for x in dim(v.base))
        raise ValueError('unjoined dimension '+str(v))
    coefficient_dimension=dim(mixed)
    measure=fourier_contract(c2);J.emit('native-fourier-contract',measure)
    edge_count=measure['invariantEdgeCoordinates']
    reduced_hat_dimension=tuple(a-b for a,b in zip(measure['nativeProfileDimension'],(edge_count,0,0)))
    height_hat_dimension=tuple(a+b for a,b in zip(schema['W_0'],reduced_hat_dimension))
    jet_hat_dimension=reduced_hat_dimension;transfer_measure_dimension=(-(measure['coordinates']-edge_count),0,0)
    J.zero('edge-reduced-profile-factor',(2*sp.pi)**measure['profileForwardPower']*(2*sp.pi)**edge_count,(2*sp.pi)**-1)
    J.emit('edge-delta-reduction',{'sourceDimension':measure['coordinates'],'invariantCoordinates':edge_count,
        'profileDependsOnlyOnCoordinate':1,'deltaFactors':[sp.DiracDelta(sp.Symbol('increment_edge_transfer_'+str(i),real=True)) for i in range(2,measure['coordinates']+1)],
        'normalizedOneDimensionalFactor':(2*sp.pi)**-1,'plainReducedMiddleMeasure':True,
        'identity':'integral exp(-i*d*y) dy = 2*pi*delta(d) in each invariant edge coordinate',
        'nativeFirstAndSecondSourceInverseSame':measure['sourceInversePower']})
    computed=[sum(v[i] for v in (coefficient_dimension,height_hat_dimension,jet_hat_dimension,transfer_measure_dimension)) for i in range(3)]
    J.emit('dimension-and-measure-join',{'C':coefficient_dimension,'heightHat':height_hat_dimension,'jetHat':jet_hat_dimension,'dt':transfer_measure_dimension,'nativeRhoDimension':schema['rho_m'],'nativeOmegaDimension':schema['omega'],
        'computedReduced':computed,'savedReduced':reduced,'savedNative':full,'factoredEdges':edge_count,
        'profiles':function_source(c2,'profile_bindings'),'kernelApply':function_source(c2,'kernel_apply'),
        'fourierProfiles':function_source(c2,'fourier_profiles'),'executableContract':measure})
    require(computed==reduced and all(x==[computed[0]+2,*computed[1:]] for x in full),'native reduced/full dimensions')
    # Restore actual both-face source and full consumer observations; cs absence above allows reuse.
    responses={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str}) for r in native['faceResponseSources']['cases']}
    traces={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str}) for r in native['geometry']['face_shift']['cases']}
    native_velocities={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str}) for r in native['geometry']['face_velocity']['cases']}
    native_chemical=sp.sympify(native['chemicalSource']['valueConstructorText'],locals={'Str':Str})[1]/eps
    objects=[*responses.values(),*traces.values(),*native_velocities.values(),native_chemical,S0,SM,saved['binding']['density'],mixed,*restored_raw]
    objects += [saved['input_'+label]['combined'] for label in ('plus','minus')]
    symbols=set().union(*(v.atoms(sp.Symbol) for v in objects));atom=lambda name:one_symbol(objects,name)
    numbind={a:values[a.name] for a in symbols if a.name in values and a.name not in ('eta_bg','sigma_W','epsilon_shape','c_s0')}
    densitymap={atom(name):v for name,v in saved['binding']['densityMap'].items()}
    profiles={atom(name):v for name,v in saved['binding']['profileEqualities'].items()}
    def bind(expr):
        expr=expr.subs(densitymap,simultaneous=True)
        for _ in range(len(profiles)+2):
            new=expr.xreplace(profiles)
            if new==expr:break
            expr=new
        else:raise ValueError('profile substitution cycle')
        return expr.subs(numbind,simultaneous=True)
    factors={};trace_info={}
    closed_formula=assignment_source(c2,'kernel_bridge','three_inverse')
    for face,label in [(1,'plus'),(-1,'minus')]:
        resp=responses[face];z=atom('s11cc1_dtn_operator_lab_held_'+label);identity=atom('s11cc1_identity_operator')
        definition=named(resp,'RESOLVENT_DEFINITION')[1];a=sp.expand(definition).coeff(z)
        a0=sp.cancel(bind(a).subs({eta:0,sigma:0},simultaneous=True))
        J.zero(label+'-native-resolvent-affine',definition,identity+a*z)
        source=named(resp,'DELTA_P').subs({named(resp,'RESOLVENT'):1,z:1},simultaneous=True)/eps
        source_input=saved['input_'+label]
        J.zero(label+'-saved-source-argument',source,source_input['raw'])
        native_velocity=bind(native_velocities[face]/eps);native_mu=bind(native_chemical)
        J.zero(label+'-saved-velocity-amplitude',native_velocity,source_input['velocityAmplitude'])
        J.zero(label+'-saved-chemical-amplitude',native_mu,source_input['chemicalAmplitude'])
        V=atom('s11cc1_V_lab_held_'+label);M=atom('s11cc1_mu_theta_lab_held_'+label)
        reconstructed=bind(source.subs({V:native_velocity,M:native_mu},simultaneous=True))
        J.zero(label+'-saved-combined-source',reconstructed,source_input['combined'])
        reduced_source=sp.cancel(source_input['combined']);den0=sp.fraction(reduced_source)[1].subs({eta:0,sigma:0},simultaneous=True)
        J.emit(label+'-source-zero-grade-operands',{'combined':source_input['combined'],'reduced':reduced_source,'denominatorAtZero':den0})
        require(den0.is_finite is True and den0.is_zero is False,'regular saved source grade')
        source0=reduced_source.subs({eta:0,sigma:0},simultaneous=True)
        J.zero(label+'-saved-source-zero-grade',source0,saved['source_'+label]['source00'])
        # New general external-leg identity keeps existing first-shape symbols once.
        z01,z12,D=sp.symbols(label+'_z01 '+label+'_z12 '+label+'_raw_direct')
        Z=sp.Matrix([[mass*freq/qo,z01,D],[0,mass*freq/qh,z12],[0,0,mass*freq/qi]])
        namespace={'sp':sp,'coefficient':a0,'z_three':Z};exec(closed_formula,namespace)
        closed=namespace['three_inverse']*Z
        Rprod=qi*qo/((qi+a0*mass*freq)*(qo+a0*mass*freq))
        J.emit(label+'-closure-operands',{'nativeDefinition':definition,'coefficient':a0,'nativeConstructor':closed_formula,'matrix':Z,
            'closed':closed,'factor':Rprod,'domain':[qi+a0*mass*freq,qo+a0*mass*freq]})
        J.zero(label+'-direct-once',sp.diff(closed[0,2],D),Rprod)
        J.zero(label+'-iteration-unchanged',closed[0,2]-closed[0,2].subs(D,0),D*Rprod)
        trace=bind(traces[face][0]/eps);pp=atom('delta_p_'+label);jj=atom('d_w_delta_p_'+label)
        vc=sp.diff(trace,pp);height=sp.diff(trace,jj);constant=trace.subs({pp:0,jj:0},simultaneous=True)
        J.zero(label+'-trace-reconstruction',trace,vc*pp+height*jj+constant)
        J.zero(label+'-bound-native-height',height,bind(face*atom('W_0')*eta*atom('w1_profile')/2))
        vc0=vc.subs({eta:0,sigma:0},simultaneous=True);ht0=height.subs({eta:0,sigma:0},simultaneous=True)
        trace0=vc0+sp.I*face*qo*ht0
        J.emit(label+'-trace-domain',{'trace':trace,'valueCoefficient':vc,'height':height,'zeroTrace':trace0,'normalJet':sp.I*face*qo})
        require(trace0!=0,'nonzero native reference trace')
        reference=sp.cancel(Rprod/trace0);jet=sp.I*face*qo*reference
        factors[label]=(reference,jet);trace_info[label]={'a':a0,'trace':trace0}
        J.emit(label+'-closed-before-cancel',{'Rprod':Rprod,'mixed':mixed.subs(numeric),'raw':raw,'reference':reference,'jet':jet})
        Cclosed=sp.cancel(Rprod*mixed.subs(numeric))
        expected=sp.I*mass*freq*(-H*qi**2+k*qh*qi+k*qh*qo-k*qh*qs-k*qi**2)/(qh*(qi+a0*mass*freq)*(qo+a0*mass*freq))
        J.zero(label+'-external-depth-cancellation',Cclosed,expected)
        J.emit(label+'-closed-raw-increment',{'physical':sp.cancel(raw*Rprod),'reference':sp.cancel(raw*reference),'normalJet':sp.cancel(raw*jet),
            'explicitExternalPolesCancelled':True,'uniformGrazingLimitEstablished':False,'domain':'nonzero external depths and resolvent/trace denominators; internal branch interpretation separate'})
    fp=factors['plus'][0].subs(qpoint)
    J.zero('saved-upper-external-factors',fp,saved['trace_factors']['referencePressurePerUnitV']/saved['trace_factors']['velocityCoefficient'])
    J.zero('saved-upper-normal-jet',factors['plus'][1].subs(qpoint),saved['trace_factors']['normalJetPerUnitV']/saved['trace_factors']['velocityCoefficient'])
    # Raw difference through full published row operands, with grade restriction AFTER substitution.
    Dplus,Dminus=sp.symbols('increment_raw_plus increment_raw_minus')
    replacements={};sources={'plus':S0,'minus':SM}
    for label,D in [('plus',Dplus),('minus',Dminus)]:
        replacements[atom('delta_p_'+label)]=eta*sigma*D*factors[label][0]*sources[label]
        replacements[atom('d_w_delta_p_'+label)]=eta*sigma*D*factors[label][1]*sources[label]
    increments={};consumer_bases={}
    for row in manifest['rowNames']:
        old=saved['consumer_'+row]
        native_row=next(r for r in native['slabConsumers']['rows'] if (('U'+str(r['address'][-1])) if r['address'][-2]=='EXPANDED' else r['address'][-2])==row)
        native_expr=sp.Add(*(sp.sympify(c['constructorText'],locals={'Str':Str}) for c in native_row['selectedPressureJetChildren']))
        J.zero(row+'-native-raw-consumer-join',old['raw'],native_expr)
        record=manifest['nativeRows'][row];full_file=Path(record['path']).read_text()
        require(full_file.endswith('\n') and not full_file.endswith('\n\n'),'one row-file terminator')
        full_text=full_file[:-1]
        require(sha(record['path'])==record['sha256'],'actual expanded row pin')
        require(hashlib.sha256(full_text.encode()).hexdigest()==native_row['expandedConstructorSha256'],'native expanded row identity')
        census=pressure_census(full_text);J.emit(row+'-executed-native-census',census)
        require(census['totalChildren']==native_row['totalChildren'],'native row child count')
        require(census['occurrenceCounts']==native_row['pressureSymbolOccurrenceCounts'],'all native slot occurrences')
        expected=[(r['childIndex'],r['constructorSha256']) for r in native_row['selectedPressureJetChildren']]
        actual=[(r['childIndex'],r['constructorSha256']) for r in census['selected']]
        require(actual==expected and census['maximumPressureDegree']<=1,'complete affine native slot sum')
        base=old['bound'];J.zero(row+'-saved-bound-join',bind(old['raw']),base)
        consumer_bases[row]=base;new=base.subs(replacements,simultaneous=True)
        zero=base.subs({x:0 for x in replacements},simultaneous=True)
        J.emit(row+'-full-increment-input',{'sourceAddress':native_row['address'],'completeCensus':native_row['allPressureOccurrencesCovered'],
            'completePublishedBoundPressureSlotSum':base,'slotInsertion':{str(a):v for a,v in replacements.items()},'ungradedDifference':new-zero})
        require(native_row['allPressureOccurrencesCovered'] is True,'complete row slot census')
        increment=sp.cancel(sp.diff(new-zero,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
        raw_density=increment.xreplace({Dplus:raw,Dminus:raw_lower})
        physical_density=raw_density.xreplace(physical_depth_bindings)
        physical_replacements={slot:value.xreplace({Dplus:raw.xreplace(physical_depth_bindings),Dminus:raw_lower.xreplace(physical_depth_bindings)}).xreplace(physical_depth_bindings) for slot,value in replacements.items()}
        physical_inserted=base.subs(physical_replacements,simultaneous=True)-zero
        physical_retained=sp.diff(physical_inserted,eta,sigma).subs({eta:0,sigma:0},simultaneous=True)
        J.zero(row+'-physical-row-sheet-join',physical_density,physical_retained)
        J.emit(row+'-retained-increment',{'mixedCoefficient':increment,'grades':[1,1],'epsilon':eps,
            'rawKernelPlus':raw,'rawKernelMinus':raw_lower,'rawRowDensity':raw_density,'physicalRowDensity':physical_density,'sourcePlus':S0,'sourceMinus':SM,
            'wholeConvolutionInserted':False,'middleIntegralEvaluated':False})
        if row in ('U0','U1','U2'):require(increment==0,'native empty U pressure increment')
        increments[row]=increment
    # Controls use the actual new lower equation and actual row insertions.
    wrongheight=J.stage('lower-wrong-height-control',{'labHeightSign':1,'required':-1},lambda:lower_boundary('lower-wrong-height',height_face=1))
    omitted=J.stage('lower-slope-omission-control',{'slopeMultiplier':0},lambda:lower_boundary('lower-slope-omitted',slope_multiplier=0))
    control_k=sp.Rational(1,20);control_H=sp.Rational(1,25);control_Q=sp.Rational(1,10)
    control_depths={qi:q(control_k),qh:q(control_k+control_H),qs:q(control_k+control_Q-control_H),qo:q(control_k+control_Q)}
    sample={k:control_k,H:control_H,Q:control_Q,om:freq,rho:mass,**{a:v.subs(cs,values['c_s0']) for a,v in control_depths.items()}}
    J.emit('nonzero-input-control-point',{'sample':{str(a):v for a,v in sample.items()},'cs':values['c_s0'],'depths':control_depths})
    J.zero('control-height-dispersion',(qh**2-qi**2+H*(H+2*k)).subs(sample,simultaneous=True),sp.S.Zero)
    J.zero('control-slope-dispersion',(qo**2-qs**2+H*(2*k+2*Q-H)).subs(sample,simultaneous=True),sp.S.Zero)
    baseline=mixed.subs(sample,simultaneous=True)
    J.nonzero('wrong-lab-face-response',baseline,wrongheight[3].subs(sample,simultaneous=True))
    J.nonzero('native-slope-omission-response',baseline,omitted[3].subs(sample,simultaneous=True))
    J.nonzero('wrong-sheet-response',baseline,mixed.subs(qh,-qh).subs(sample,simultaneous=True))
    for row in ('THETA_BALANCE','E_W_BALANCE'):
        expr=increments[row];probe={a:0 for a in expr.free_symbols if a.name in ('theta','e_W_t') or a.name.startswith(('u_','e_W_d','theta_d'))}
        probe.update({a:1 for a in expr.free_symbols if a.name=='e_W'})
        probe.update({Dplus:1,Dminus:1,eps:1,qi:sample[qi],qo:sample[qo]})
        good=sp.cancel(expr.subs(probe,simultaneous=True))
        base=consumer_bases[row]
        def changed_insertion(multiplier):
            changed={slot:value.xreplace({Dminus:multiplier*Dminus}) for slot,value in replacements.items()}
            difference=base.subs(changed,simultaneous=True)-base.subs({slot:0 for slot in changed},simultaneous=True)
            J.emit(row+'-changed-insertion-'+str(multiplier),{'replacements':{str(a):v for a,v in changed.items()},'ungradedDifference':difference})
            return sp.cancel(sp.diff(difference,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
        without=sp.cancel(changed_insertion(0).subs(probe,simultaneous=True))
        J.nonzero(row+'-actual-lower-row-control',good,without)
        doubled=sp.cancel(changed_insertion(2).subs(probe,simultaneous=True))
        J.nonzero(row+'-actual-double-row-control',good,doubled)
    J.nonzero('lower-reference-jet-sign',factors['minus'][1].subs(qpoint),(-factors['minus'][1]).subs(qpoint))
    # No endpoint integration is smuggled into a raw-kernel job.
    kap=sp.sqrt(freq**2/cs**2-sum(e*e for e in edge))
    J.emit('applicability',{'rawNongrazingOnly':True,'branchPoints':[-k-kap,-k+kap,k+Q-kap,k+Q+kap],
        'realBranchCondition':freq**2/cs**2>sum(e*e for e in edge),'contact':0,
        'externalConditions':['qi!=0','qo!=0','native resolvent and reference denominators!=0'],
        'constantEndRule':'Native slope amplitude/jet is zero, hence added direct term vanishes on either constant end; baseline retained.',
        'exactMatch':'UNRESOLVED: physical two-sided closed collision and integrable local bound are not computed by this nongrazing instrument',
        'grazingScale':'Use exact branch roots; at k=0 the root scale is qi, not qi^2/(2|k|).',
        'uniformObjectsRevalidated':False,'finiteMatrixChanged':False,'lossInterpretation':'NOT SUPPLIED'})
    for label,C in [('plus',mixed),('minus',lower[3])]:J.zero(label+'-constant-end-increment',(h*s*C).subs(s,0),sp.S.Zero)
    return {'candidateRawIncrementConstructed':True,'bothFaces':True,'nongrazingDomainOnly':True,
            'exactMatchApplicability':'UNRESOLVED','computedIntegral':False,'finiteSolves':0,
            'U012IncrementStatus':'literal zero required by empty native pressure census',
            'productionChanges':False,'physicalLossSupplied':False,'replayedCompletedFunctions':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.inputs.read_text())
    for path,h in manifest['sourcePins'].items():require(sha(path)==h,'input pin '+path)
    gate=verify_gate(args.gate,args.inputs,manifest)
    inspection_pins={**manifest['sourcePins'],**gate['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate),
        gate['buildReviewRecord']:gate['buildReviewRecordSha256']}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;start=time.monotonic()
    try:
        save(args.out/'containment.json',containment())
        global sp,Str
        import sympy as sp
        from sympy.core.symbol import Str
        J=Journal(args.out)
        result=J.stage('bounded-raw-construction',{'manifestSha256':sha(args.inputs),'savedOperands':manifest['savedOperands']},lambda:scientific_work(manifest,J));result['executionStatus']='COMPLETED_CANDIDATE_RAW_INCREMENT';code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False}
        save(args.out/'failure.json',result)
    finally:
        for alias,record in manifest['savedOperands'].items():
            copy=args.out/'saved-operands'/(alias+'.json')
            if copy.exists():inspection_pins[str(copy)]=record['sha256']
        post=posthash_records(inspection_pins);save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):result['integrityFailure']=True;code=1
        if J is not None:replace_json(args.out/'operation-index.json',J.completed)
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
