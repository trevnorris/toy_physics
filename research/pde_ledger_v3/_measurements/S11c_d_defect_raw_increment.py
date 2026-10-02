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
    require(r['independentBuildClearance'] is True,'independent build assessment')
    for k in ('workerSha256','manifestSha256','guardSha256','supervisorSha256'):
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
        'externalDomain':[sp.Ne(qi,0),sp.Ne(qo,0)],'reflectedRoute':'k+k_out-m; separate from native MIDDLE_Q'})
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
        reduced=sp.expand(sp.expand(num).subs(rules,simultaneous=True))
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
    J.zero('saved-selected-integrand',selected,old_action['integrand'])
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
    height_hat_dimension=tuple(a+b for a,b in zip(schema['W_0'],(1,0,0)))
    jet_hat_dimension=(1,0,0);transfer_measure_dimension=(-1,0,0)
    computed=[sum(v[i] for v in (coefficient_dimension,height_hat_dimension,jet_hat_dimension,transfer_measure_dimension)) for i in range(3)]
    J.emit('dimension-and-measure-join',{'C':coefficient_dimension,'heightHat':height_hat_dimension,'jetHat':jet_hat_dimension,'dt':transfer_measure_dimension,'nativeRhoDimension':schema['rho_m'],'nativeOmegaDimension':schema['omega'],
        'computedReduced':computed,'savedReduced':reduced,'savedNative':full,'factoredEdges':2,
        'profiles':function_source(c2,'profile_bindings'),'kernelApply':function_source(c2,'kernel_apply'),
        'fourierProfiles':function_source(c2,'fourier_profiles'),'sourceForwardNormalized':False,'profileForwardNormalized':True})
    require(computed==reduced and all(x==[computed[0]+2,*computed[1:]] for x in full),'native reduced/full dimensions')
    # Restore actual both-face source and full consumer observations; cs absence above allows reuse.
    responses={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str}) for r in native['faceResponseSources']['cases']}
    traces={int(r['case'][1]):sp.sympify(r['valueConstructorText'],locals={'Str':Str}) for r in native['geometry']['face_shift']['cases']}
    objects=[*responses.values(),*traces.values(),S0,SM,saved['binding']['density'],mixed]
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
        J.zero(label+'-saved-source-argument',source,saved['input_'+label]['raw'])
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
        vc0=vc.subs({eta:0,sigma:0},simultaneous=True);ht0=height.subs({eta:0,sigma:0},simultaneous=True)
        trace0=vc0+sp.I*face*qo*ht0
        J.emit(label+'-trace-domain',{'trace':trace,'valueCoefficient':vc,'height':height,'zeroTrace':trace0,'normalJet':sp.I*face*qo})
        require(trace0!=0,'nonzero native reference trace')
        reference=sp.cancel(Rprod/trace0);jet=sp.I*face*qo*reference
        factors[label]=(reference,jet);trace_info[label]={'a':a0,'trace':trace0}
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
    increments={}
    for row in manifest['rowNames']:
        old=saved['consumer_'+row]
        native_row=next(r for r in native['slabConsumers']['rows'] if (('U'+str(r['address'][-1])) if r['address'][-2]=='EXPANDED' else r['address'][-2])==row)
        native_expr=sp.Add(*(sp.sympify(c['constructorText'],locals={'Str':Str}) for c in native_row['selectedPressureJetChildren']))
        J.zero(row+'-native-raw-consumer-join',old['raw'],native_expr)
        base=old['bound'];new=base.subs(replacements,simultaneous=True)
        zero=base.subs({x:0 for x in replacements},simultaneous=True)
        J.emit(row+'-full-increment-input',{'sourceAddress':native_row['address'],'completeCensus':native_row['allPressureOccurrencesCovered'],
            'fullPublishedBoundRow':base,'slotInsertion':{str(a):v for a,v in replacements.items()},'ungradedDifference':new-zero})
        require(native_row['allPressureOccurrencesCovered'] is True,'complete row slot census')
        increment=sp.cancel(sp.diff(new-zero,eta,sigma).subs({eta:0,sigma:0},simultaneous=True))
        J.emit(row+'-retained-increment',{'mixedCoefficient':increment,'grades':[1,1],'epsilon':eps,
            'rawKernelPlus':raw,'rawKernelMinus':raw_lower,'rawRowDensity':increment.xreplace({Dplus:raw,Dminus:raw_lower}),'sourcePlus':S0,'sourceMinus':SM,
            'wholeConvolutionInserted':False,'middleIntegralEvaluated':False})
        increments[row]=increment
    # Controls use the actual new lower equation and actual row insertions.
    wrongheight=J.stage('lower-wrong-height-control',{'labHeightSign':1,'required':-1},lambda:lower_boundary('lower-wrong-height',height_face=1))
    omitted=J.stage('lower-slope-omission-control',{'slopeMultiplier':0},lambda:lower_boundary('lower-slope-omitted',slope_multiplier=0))
    sample={k:0,H:sp.Rational(1,20),qi:sp.Rational(1,5),qo:sp.sqrt(3)/10,qh:sp.sqrt(15)/20,qs:sp.sqrt(15)/20,om:freq,rho:mass}
    baseline=mixed.subs(sample,simultaneous=True)
    J.nonzero('wrong-lab-face-response',baseline,wrongheight[3].subs(sample,simultaneous=True))
    J.nonzero('native-slope-omission-response',baseline,omitted[3].subs(sample,simultaneous=True))
    J.nonzero('wrong-sheet-response',baseline,mixed.subs(qh,-qh).subs(sample,simultaneous=True))
    J.nonzero('double-insertion-response',baseline,2*baseline)
    for row in ('THETA_BALANCE','E_W_BALANCE'):
        expr=increments[row];probe={a:0 for a in expr.free_symbols if a.name in ('theta','e_W_t') or a.name.startswith(('u_','e_W_d','theta_d'))}
        probe.update({a:1 for a in expr.free_symbols if a.name=='e_W'})
        probe.update({Dplus:1,Dminus:1,eps:1,**qpoint})
        good=sp.cancel(expr.subs(probe,simultaneous=True))
        without=sp.cancel(expr.subs(Dminus,0).subs(probe,simultaneous=True))
        J.nonzero(row+'-actual-lower-row-control',good,without)
        doubled=sp.cancel(expr.subs(Dminus,2*Dminus).subs(probe,simultaneous=True))
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
            'productionChanges':False,'physicalLossSupplied':False,'replayedCompletedFunctions':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.inputs.read_text())
    for path,h in manifest['sourcePins'].items():require(sha(path)==h,'input pin '+path)
    verify_gate(args.gate,args.inputs,manifest)
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
        post={p:{'expected':h,'actual':sha(p)} for p,h in manifest['sourcePins'].items()};save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):result['integrityFailure']=True;code=1
        if J is not None:replace_json(args.out/'operation-index.json',J.completed)
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
