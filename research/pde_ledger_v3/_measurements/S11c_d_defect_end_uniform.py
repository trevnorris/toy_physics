#!/usr/bin/env python3
"""Bounded source-joined end/uniform equation comparison; no producer or roots."""
import argparse,ast,base64,builtins,hashlib,importlib,io,json,math,os,pickle,re,resource,shutil,sqlite3,sys,time,traceback
from collections import OrderedDict
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
G=((0,0),(1,0),(0,1),(1,1))
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(value,message):
    if value is not True:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
    with path.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())


def definitions(text,names):
    nodes=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in nodes}==set(names),'inert definition census')
    return ast.Module(body=nodes,type_ignores=[])


def source_method(text,cls,name):
    c=next(n for n in ast.parse(text).body if isinstance(n,ast.ClassDef) and n.name==cls)
    return ast.unparse(next(n for n in c.body if isinstance(n,ast.FunctionDef) and n.name==name))


class ExpectedUnavailable(ValueError):
    """A bounded representation/domain classification, never an integrity waiver."""


class SourceMapUnresolved(ValueError):pass


def source_fragment(source,scope,expected):
    node=ast.parse(source)
    for name in scope:
        found=[n for n in ast.walk(node) if n is not node and isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name==name]
        require(len(found)==1,'unique native source scope '+name);node=found[0]
    want=ast.parse(expected).body[0];dump=ast.dump(want,include_attributes=False)
    found=[n for n in ast.walk(node) if isinstance(n,type(want)) and ast.dump(n,include_attributes=False)==dump]
    require(len(found)==1,'native source contract '+expected)
    return {'scope':list(scope),'statement':ast.unparse(found[0]),'sourceSha256':hashlib.sha256(source.encode()).hexdigest()}


def source_chart_and_scale(m,J,D,raw,physical,context,depth,U,L,profile):
    """New convention joins only. No old source, lift, limit or symbol call."""
    engine=Path(m['engineSource']).read_text();uniform=Path(m['uniformWorker']).read_text()
    ends=Path(m['endsWorker']).read_text();full=Path(m['fullWeakWorker']).read_text()
    checks=[]
    def contract(text,scope,statement):
        r=source_fragment(text,scope,statement);checks.append(r);J.emit('native-contract-%03d'%len(checks),r)
    contract(engine,('EdgeReduction','__init__'),"self.normal_map = {self.x[2]: self.z, self.y[2]: self.zp}")
    contract(engine,('EdgeReduction','__init__'),"self.tangent_momenta = tuple(k for group in self.momentum_groups for k in group[:2])")
    contract(engine,('EdgeReduction','wave_phase'),"return sp.exp(sp.I * (sum(k*x for k,x in zip(self.tangents,point[:2]))-self.omega*self.t))")
    contract(engine,('ReducedPencil','trial_ansatz','displacement','derivative'),"return self.tangent_derivatives[j]*value if j < 2 else sp.diff(value,z)")
    contract(uniform,('selected_lift',),"phase = sp.exp(sp.I*(t1*x1+t2*x2+k*z))")
    contract(full,('run_science','derive'),"multiplier = (-sp.I*3)**spec['timeOrder']*(sp.I/5)**spec['spatialOrders'][1]*(sp.I/10)**spec['spatialOrders'][2]")
    contract(ends,('run_science',),"H = {'minus':sp.S.Zero,'plus':W/2}")
    contract(ends,('run_science',),"for side,endkey in [('minus','leftEndpoint'),('plus','rightEndpoint')]:\n    val=D(c[endkey]);require(not val.free_symbols,'saved local endpoint constant')\n    term=val*(sp.I*p)**c['xOrder'];rec['newTerms'][side]=term\n    group=groups[(side,c['row'],c['field'],tuple(c['grade']))]\n    group['local'].append(term);group['localAncestry'].append({'row':c['row'],'field':c['field'],'xOrder':c['xOrder'],'grade':c['grade']})")
    contract(uniform,('science',),"for label,end in (('LEFT',-sp.oo),('RIGHT',sp.oo)):\n    operands=restored['uniformSource']['records'][label]['background']['profileLimitOperands']\n    for original,limit,endpoint in operands:\n        name=next(name for name in profiles if original.func.__name__=='s11cd'+name.upper()+'Profile')\n        evaluated=sp.limit(profiles[name],xi,end)\n        jets=tuple(sp.limit(sp.diff(profiles[name],xi,order),xi,end) for order in (1,2))\n        supplied=restored['uniformSource']['profileEndpoints'][endpoint]\n        profile_joins.append(dict(end=label,source=(original,limit,endpoint),definition=profiles[name],actualEnd=evaluated,supplied=supplied,residual=evaluated-supplied,jets=jets))")
    phase=D(raw['ends/phase-arguments.json']);bindings=phase['bindings'];p=depth['p'];edge=tuple(sp.Rational(physical['parameters']['s11cdTangentialMomentum'+str(i)]) for i in (1,2))
    J.join('phase-source-contract',phase['sourceExpression'],raw['native/weak-duality.json']['inheritedConvention']['sourcePhase'])
    J.join('phase-profile-contract',phase['profileExpression'],raw['native/weak-duality.json']['inheritedConvention']['profilePhase'])
    kout,kin,ko,ki,X,Y=(bindings[k] for k in ('kout','kin','ko','ki','X','Y'))
    J.join('actual-weak-output-tangents',tuple(kout[1:]),edge);J.join('actual-weak-input-tangents',tuple(kin[1:]),edge)
    J.join('actual-profile-output',ko,kout);J.join('actual-profile-input',ki,kin)
    require(kout[0].name=='weak_end_l' and kin[0].name=='weak_end_k' and len(set(X+Y))==6,'actual phase coordinate names')
    J.zero('source-phase-argument',phase['sourceExponent'],sp.I*(sum(k*x for k,x in zip(kout,X))-sum(k*y for k,y in zip(kin,Y))))
    J.zero('profile-phase-argument',phase['profileExponent'],-sp.I*sum((a-b)*y for a,b,y in zip(ko,ki,Y)))
    profile_rates=tuple(sp.diff(phase['profileExponent'],y) for y in Y)
    J.join('profile-normal-direction',profile_rates,(-sp.I*(ko[0]-ki[0]),sp.S.Zero,sp.S.Zero))
    # Tangent order and sign are distinguished by unequal actual tangents.
    S=sp.ImmutableMatrix([[0,0,1,0,0],[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    weak_mom=sp.ImmutableMatrix((p,*edge));old_mom=sp.ImmutableMatrix((*edge,p));mapping=S[:3,:3]
    J.join('actual-carrier-momentum-map',mapping*old_mom,weak_mom)
    J.join('actual-carrier-inverse-map',mapping.T*weak_mom,old_mom)
    require(mapping.det()==1 and mapping.T*mapping==sp.eye(3),'oriented coordinate chart')
    J.join('restored-lift-native-method',L['source']['sourceMethod'],source_method(engine,'ReducedPencil','trial_ansatz'))
    J.join('restored-lift-native-field-method',L['source']['sourceFieldLift'],source_method(engine,'ConstantEndPencil','field_lift'))
    contract(engine,('ClosedCurrentPairing','construct'),"phases = tuple(sp.exp(sign*sp.I*(sum(k*x for k,x in zip(r.tangents,r.x[:2]))+momentum*r.z-frequency*r.t+c.phase_coordinate)) for sign,momentum,frequency in ((1,kright,omega_right),(-1,kleft,omega_left)))")
    require(context['context']['numeric']['omega']==3 and raw['native/wave-profile.json']['pressureTimeRule']=="(-sp.I * 3) ** spec['timeOrder']",'same negative temporal rate')
    J.emit('actual-source-chart',{'matrixWeakFromUniform':S,'weakCarrier':weak_mom,'uniformCarrier':old_mom,'timeRate':-3*sp.I,'profileRates':profile_rates,'sourceContracts':checks,'noResidualUsedToChooseMap':True})
    heights=D(raw['controls/conclusion.json']['heights']);W=sp.Rational(physical['parameters']['W_0']);length=sp.Rational(physical['parameters']['L_W'])
    J.join('actual-half-heights',heights,{'minus':sp.S.Zero,'plus':W/2});require(length==10 and length.is_positive is True,'orientation preserving profile length')
    side_map={};side_evidence=[]
    for end in ('LEFT','RIGHT'):
        records=[r for r in profile['profileJoins'] if r['end']==end];require(len(records)==2,'two actual profile records per end')
        expected_limit=-sp.oo if end=='LEFT' else sp.oo
        for r in records:
            original,limit,endpoint=r['source'];J.emit('side-'+end+'-'+original.func.__name__,r)
            require((original,limit,endpoint) in U['records'][end]['background']['profileLimitOperands'],'actual old endpoint argument')
            require(isinstance(limit,sp.Limit),'native endpoint is Limit')
            require(limit.args[2]==expected_limit and limit.args[0]==original,'native profile orientation')
            name=next(n for n in ('w','m') if original.func.__name__=='s11cd'+n.upper()+'Profile')
            xi=next(iter(r['definition'].free_symbols));definition=sp.sympify(physical['profiles'][name],locals={'xi':xi,'tanh':sp.tanh})
            J.join('side-'+end+'-'+name+'-profile',r['definition'],definition)
            require(r['actualEnd']==r['supplied']==U['profileEndpoints'][endpoint] and r['residual']==0 and all(z==0 for z in r['jets']),'inherited exact endpoint facts')
            require(r['actualEnd']==(0 if name=='m' or end=='LEFT' else 1),'actual tanh endpoint value')
            if name=='w':
                candidates=[side for side,h in heights.items() if h==W*r['actualEnd']/2];require(len(candidates)==1,'unique native weak side')
                side_map[end]=candidates[0];side_evidence.append({'end':end,'nativeLimit':limit,'profileEnd':r['actualEnd'],'weakSide':candidates[0],'halfHeight':heights[candidates[0]],'profileLength':length})
    J.emit('actual-end-side-map',side_evidence);require(set(side_map.values())=={'minus','plus'},'bijection of native sides')
    # Join the extraction paths, not a guessed factor inferred from A or R.
    contract(engine,('ConstantEndPencil','strong_matrix','strip_phase'),"terms = self.phase_terms(sp.diff(expression,epsilon),self.r.z)")
    contract(engine,('ConstantEndPencil','strong_matrix','strip_phase'),"return terms.get(self.kn,sp.S.Zero)")
    contract(engine,('ConstantEndPencil','strong_matrix'),"return sp.ImmutableMatrix.hstack(*(sp.ImmutableMatrix(c) for c in columns))")
    contract(engine,('ClosedCurrentPairing','construct'),"algebraic,relation,joins = self.modes.analytic(a.strong)")
    contract(engine,('ClosedCurrentPairing','construct'),"physical_q = algebraic.subs(self.modes.q,a.qlegs[1]/acoustic['ACOUSTIC_RADICAL_SCALE'])")
    contract(engine,('ClosedCurrentPairing','construct'),"pencils = tuple(physical_q.xreplace(mapping) for mapping in leg_maps)")
    contract(full,('run_science','derive'),"raw_coefficient = sp.cancel(original/(eps*wave))")
    contract(full,('run_science','derive'),"zero_record(ev,'native-child-reconstruction',original,eps*wave*raw_coefficient)")
    contract(ends,('run_science',),"rec['symbol'] = sp.cancel(rec['localSum']+rec['pressureSum'])")
    contract(ends,('run_science',),"rec['sumIdentity'] = zero('new-end-symbol-sum',rec['symbol'],sum(v['local']+v['pressure'],sp.S.Zero))")
    eps=context['context']['epsilon']
    dual=D(raw['native/weak-duality.json']);J.join('inherited-source-order',dual['source'],'X(k)=hat[b D_j u](k)');J.join('inherited-consumer-order',dual['consumer'],'Y(l)=2pi hat[c v](-l)')
    require(dual['inheritedConvention']['profileForwardPower']==-3 and dual['inheritedConvention']['sourceInversePower']==-3 and dual['inheritedConvention']['sourceForwardPower']==0 and dual['inheritedConvention']['invariantEdgeCoordinates']==2,'inherited Fourier factors')
    # Two invariant edge integrations cancel two of the three inverse factors.
    edge_reduced=(2*sp.pi)**dual['inheritedConvention']['invariantEdgeCoordinates']*(2*sp.pi)**dual['inheritedConvention']['profileForwardPower']
    J.zero('one-dimensional-forward-normalization',edge_reduced,1/(2*sp.pi))
    J.zero('bilinear-Fourier-dual-normalization',(2*sp.pi)*edge_reduced,1)
    require(U['constantFourierMass']==2*sp.pi,'actual original uniform Fourier mass')
    J.emit('actual-scalar-normalization',{'epsilon':eps,'coefficientFactor':sp.S.One,'weakPairingExternalFactor':2*sp.pi,'uniformConstantFourierMass':U['constantFourierMass'],'edgeReducedForwardFactor':edge_reduced,'duality':dual,'sourceContracts':checks,'savedCoefficientNormalization':'REQUIRES saved_normalization_joins before comparison','noFieldOrTransformEvaluated':True})
    return S,side_map


class SourceStage:
    """Classify source correspondence failures consistently; never soften them."""
    def __init__(self,J,name):self.J,self.name=J,name
    def __enter__(self):return self
    def __exit__(self,kind,value,tb):
        if kind and issubclass(kind,(ValueError,KeyError,StopIteration,TypeError)):
            self.J.emit(self.name+'-unresolved',{'status':'SOURCE_MAP_UNRESOLVED','reason':str(value),'exceptionType':kind.__name__,'activeOperation':self.J.active,'traceback':traceback.format_exc()})
            raise SourceMapUnresolved(self.name+': '+str(value)) from value
        return False


def attribution_usable(status,new_remainder_zero):
    if status=='UNRESOLVED_UNSUPPORTED_DEPTH_DEPENDENCE':return False
    require(status=='ZERO','required on-wave attribution reconstruction failed: '+status)
    return new_remainder_zero


def local_key(cell):return (cell['row'],cell['field'],cell['xOrder'],tuple(cell['grade']))


def normalization_routes(raw):
    """Metadata join only: exact saved cell identities, full coverage, native child receipts."""
    local=raw['normalization/local-cells.json'];ends=raw['normalization/local-end-terms.json']
    lm={local_key(c):c for c in local};em={local_key(c['savedCell']):c for c in ends}
    expected={(r,f,n,g) for r in ROWS for f in FIELDS for n in range(4) for g in G}
    require(len(local)==len(ends)==len(lm)==len(em)==400 and set(lm)==set(em)==expected,'normalization full local coverage')
    for k,c in lm.items():require(c==em[k]['savedCell'],'actual saved cell to endpoint arguments')
    inp=raw['normalization/batch-input.json'];ret=raw['normalization/batch-return.json']
    im={r['childIndex']:r for r in inp['children']};rm={r['childIndex']:r for r in ret}
    native={c['childIndex']:c for c in raw['native/U0.json']['children']}
    require(len(im)==len(rm)==16 and set(im)==set(rm),'complete saved local batch')
    for i,c in im.items():
        require(c==native[i] and rm[i]['sourceConstructor']==c['constructorText'] and rm[i]['sourceSha256']==c['sha256'],'actual child input/return source')
        require(hashlib.sha256(c['constructorText'].encode()).hexdigest()==c['sha256'] and rm[i]['completed'] is True,'completed native child receipt')
    return lm,em,rm


def saved_normalization_joins(m,J,D,raw,cells,context,p,params):
    """New cross-stage joins from published operands; no local producer or endpoint replay."""
    lm,em,rm=normalization_routes(raw)
    receipt=next(r for r in raw['normalization/operation-index.json'] if r['name']=='U0-local-batch-007')
    for slot,alias in [('input','batch-input'),('result','batch-return')]:
        pin=m['savedInputs']['normalization/'+alias+'.json'];r=receipt[slot]
        require(r['path']==Path(pin['path']).name and r['sha256']==pin['sha256'] and r['bytes']==pin['bytes'],'actual completed normalization batch receipt')
    J.emit('normalization-inherited-batch-receipt',receipt)
    decoded={};terms={}
    for index,(k,c) in enumerate(lm.items()):
        n='normalization-local-%03d'%index;v=D(c);e=D(em[k]['newTerms']);decoded[k]=v;terms[k]=e
        proof=next(x for x in v['identities'] if x['name']=='cell-source-sum')
        J.join(n+'-sum-left',proof['left'],v['coefficient'])
        J.zero(n+'-saved-sum-argument',proof['right'],sum(v['summands'],sp.S.Zero))
        require(proof['cancelled']==0 and len(v['sourceChildren'])==len(v['summands']),'inherited local sum proof and ancestry')
        for side,endpoint in [('minus','leftEndpoint'),('plus','rightEndpoint')]:
            J.zero(n+'-'+side+'-new-term-join',e[side],v[endpoint]*(sp.I*p)**v['xOrder'])
        J.emit(n+'-inherited-proof',{'key':k,'coefficient':v['coefficient'],'sourceChildren':v['sourceChildren'],'sumProof':proof,'endpoints':(v['leftEndpoint'],v['rightEndpoint']),'newTerms':e,'oldFunctionCalled':False})
    for index,c in enumerate(cells):
        n='normalization-end-cell-%03d'%index
        keys=[(c['row'],c['field'],order,tuple(c['grade'])) for order in range(4)]
        expected=[{'row':k[0],'field':k[1],'xOrder':k[2],'grade':list(k[3])} for k in keys]
        J.join(n+'-ancestry',c['localAncestry'],expected)
        J.join(n+'-local-terms',c['local'],[terms[k][c['side']] for k in keys])
        J.zero(n+'-local-sum',c['localSum'],sum(c['local'],sp.S.Zero))
        J.zero(n+'-pressure-sum',c['pressureSum'],sum(c['pressure'],sp.S.Zero))
        J.zero(n+'-sum-proof-argument',c['sumIdentity']['right'],c['localSum']+c['pressureSum'])
        J.join(n+'-symbol-proof-left',c['sumIdentity']['left'],c['symbol'])
        require(c['sumIdentity']['cancelled']==0,'inherited actual symbol-sum return')
    eps=context['context']['epsilon'];witnesses=[]
    # Fixed temporal and two unequal tangential native children, chosen before residuals.
    for child in (122,113,114):
        v=D(rm[child]);n='normalization-native-child-'+str(child);spec=v['jet']
        original=v['original'];atoms={s.name:s for s in original.free_symbols}
        require(atoms['epsilon_shape']==eps and v['row']=='U0' and spec['channel']=='u_1' and v['fieldColumn']==0 and v['xOrder']==spec['spatialOrders'][0]==0,'actual normalization witness route')
        wave=atoms[spec['name']];proof=next(x for x in v['identities'] if x['name']=='native-child-reconstruction')
        J.join(n+'-original-proof-left',proof['left'],original)
        J.zero(n+'-original-proof-right-join',proof['right'],eps*wave*v['nativeCoefficient'])
        require(proof['cancelled']==0 and v['profileMaps']==[],'inherited native coefficient proof, constant witness')
        material={s:params[s.name] for s in original.free_symbols-{eps,wave}}
        J.join(n+'-material-bound',v['nativeCoefficient'].xreplace(material),v['boundCoefficient'])
        multiplier=(-sp.I*3)**spec['timeOrder']*(sp.I/5)**spec['spatialOrders'][1]*(sp.I/10)**spec['spatialOrders'][2]
        J.join(n+'-native-wave-factor',multiplier,v['waveMultiplier'])
        # Extract epsilon from the actual saved native child after the actual
        # material/plane-wave jet map, per unit carrier. No fresh trial scalar.
        bound_child=original.xreplace({**material,wave:multiplier})
        extracted=sp.diff(bound_child,eps)
        mapped=next(g for g in v['mappedGrades'] if g['grade']==[0,0])
        grade=next(g for g in v['grade']['table'] if g['grade']==[0,0])
        require(v['grade']['excludedRemainder']==0 and not v['boundCoefficient'].free_symbols and mapped['profileMap']==[],'witness has only constant zero grade')
        J.zero(n+'-grade-bound-join',grade['coefficient'],v['boundCoefficient'])
        J.zero(n+'-actual-epsilon-normalization',extracted,mapped['value'])
        cell=decoded[('U0','u_1',0,(0,0))];position=cell['sourceChildren'].index(child)
        J.join(n+'-mapped-to-actual-summand',mapped['value'],cell['summands'][position])
        rec={'childIndex':child,'nativeSource':original,'materialMap':material,'jet':spec,'waveMultiplier':multiplier,'boundNativeChild':bound_child,'epsilon':eps,'singleDerivative':extracted,'savedMappedValue':mapped['value'],'cellKey':('U0','u_1',0,(0,0)),'summandPosition':position,'inheritedNativeProof':proof,'endTerms':terms[('U0','u_1',0,(0,0))]}
        J.emit(n+'-evidence',rec);witnesses.append(rec)
    result={'actualNativeWitnesses':witnesses,'localCellsJoined':len(decoded),'endCellsJoined':len(cells),'inheritedEndpointValuesRecomputed':False,'inheritedLocalOrPressureProducerCalled':False,'scope':'Actual temporal/spatial native epsilon witnesses and all local-cell/end-symbol ancestry joins; pressure normalization remains inherited from the saved source/consumer proofs, not an independent new pressure derivation.'}
    J.emit('saved-normalization-conclusion',result);return result

def merge_native_dimensions(base,atoms,registries,J):
    """Restore omitted generated units; never infer, guess or exempt a symbol."""
    dimensions={name:tuple(unit) for name,unit in base.items()};added={}
    for atom in sorted(atoms,key=lambda s:s.name):
        if atom.name in dimensions:continue
        matches=[{'registry':name,'symbol':next(k for k in values if k==atom),'unit':values[atom]}
                 for name,values in registries.items() if atom in values]
        J.emit('unit-registry-'+atom.name,{'nativeSymbol':atom,'savedMatches':matches,'staticTableMissing':True})
        require(atom.name.startswith('gamma_'),'unaccounted non-generated native unit '+atom.name)
        require(bool(matches),'missing actual saved generated unit '+atom.name)
        units=[]
        for match in matches:
            unit=match['unit']
            require(isinstance(unit,(tuple,list)) and len(unit)==3,'saved generated unit vector')
            require(all((type(v) is int) or (getattr(v,'is_Rational',False) is True) for v in unit),'saved generated exact rational unit')
            units.append(tuple(unit))
        require(all(v==units[0] for v in units),'conflicting actual saved generated unit '+atom.name)
        dimensions[atom.name]=units[0];added[atom.name]=matches
    require(all(a.name in dimensions for a in atoms),'complete actual native unit coverage')
    J.emit('unit-registry-merged',{'staticEntries':len(base),'generatedEntries':added,'completeRequiredSymbols':[a.name for a in sorted(atoms,key=lambda s:s.name)],'newDimensionInference':False})
    return dimensions


def restored_unit_registries(m,J,raw,restored):
    registries={}
    for name,origin in m['unitRegistryOrigins'].items():
        checks=raw[origin['checksAlias']];pin=m['originalPackets'][name]
        require(checks['objectsSha256']==pin['sha256'],'actual unit registry object receipt '+name)
        require(checks['provenance']['producerSources']['scripts/S11c_b_exports.py']==sha(m['bSource']),'actual unit registry native slab source '+name)
        expected=checks.get('sourceFiles',checks['provenance'].get('sourceFiles',{}))
        require(expected[origin['producerRelativePath']]==sha(origin['producerSource']),'actual unit registry producer snapshot '+name)
        body=Path(origin['producerSource']).read_text()
        contract=source_fragment(body,tuple(origin['scope']),origin['statement'])
        if name.endswith('Pairing'):
            packet,registry=restored[name]
            J.join('unit-registry-'+name+'-source-summary',packet['summary']['sourceFiles'],checks['sourceFiles'])
        else:registry=restored[name]['knownDimensions']
        require(isinstance(registry,dict),'actual saved unit registry dictionary '+name)
        registries[name]=registry
        J.emit('unit-registry-'+name+'-origin',{'objectReceipt':pin,'sourceChecks':checks,'producerContract':contract,'registryEntries':len(registry),'sourceReexecuted':False})
    return registries


def exact_structure(a,b):
    if type(a) is not type(b):return False
    if isinstance(a,dict):return a.keys()==b.keys() and all(exact_structure(a[k],b[k]) for k in a)
    if isinstance(a,(list,tuple)):return len(a)==len(b) and all(exact_structure(x,y) for x,y in zip(a,b))
    if isinstance(a,np.ndarray):return a.shape==b.shape and a.dtype==b.dtype and a.tobytes()==b.tobytes()
    r=a==b
    return r is True or r is sp.S.true


def convolution(a,b,rectangle=True):
    out={}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if not rectangle or (i+k<=1 and j+l<=1):out[i+k,j+l]=out.get((i+k,j+l),0)+x*y
    return out


def rotated_name(name,scalar_names):
    """Proper cyclic basis map: weak axes 1,2,3 are uniform axes 3,1,2."""
    cycle={'1':'3','2':'1','3':'2'}
    m=re.fullmatch(r'(u_([123])|theta|e_W|w1_profile|m1_profile)((?:_t|t)*)(.*)',name)
    if m:
        head,component,time_suffix,derivatives=m.groups()
        require(re.fullmatch(r'(?:_?d[123])*',derivatives) is not None,'native jet suffix')
        if component:head='u_'+cycle[component]
        directions=re.findall(r'd([123])',derivatives)
        suffix='' if not directions else '_'+''.join('d'+i for i in sorted(cycle[d] for d in directions))
        return head+time_suffix+suffix
    m=re.fullmatch(r'grad_theta_([123])',name)
    if m:return 'grad_theta_'+cycle[m[1]]
    require(name in scalar_names,'unclassified native tensor/scalar '+name)
    return name


class Evidence:
    def __init__(self,out):self.out=out;self.active=None;self.count=0;self.previous=None;self.operations=[]
    def encode(self,v):
        if isinstance(v,(sp.Basic,sp.MatrixBase,FunctionClass)):return {'text':str(v),'srepr':sp.srepr(v)}
        if isinstance(v,np.ndarray):return {'numpyDtype':v.dtype.str,'shape':list(v.shape),'bytesBase64':base64.b64encode(v.tobytes()).decode()}
        if isinstance(v,dict):
            if all(isinstance(k,str) for k in v):return {k:self.encode(x) for k,x in v.items()}
            return {'mappingPairs':[[self.encode(k),self.encode(x)] for k,x in v.items()]}
        if isinstance(v,(list,tuple)):return [self.encode(x) for x in v]
        if isinstance(v,(set,frozenset)):return {'setValues':[self.encode(x) for x in sorted(v,key=repr)]}
        if isinstance(v,(complex,np.complexfloating)):return {'pythonComplex':[float(v.real),float(v.imag)]}
        if isinstance(v,(np.integer,np.bool_,np.floating)):return v.item()
        require(v is None or type(v) in (str,bool,int,float),'supported evidence type '+str(type(v)))
        return v
    def emit(self,name,value):
        p=self.out/(name+'.json');save(p,self.encode(value))
        ref={'path':str(p.relative_to(self.out)),'sha256':sha(p),'bytes':p.stat().st_size}
        body={'sequence':self.count,'previousSha256':self.previous,'artifact':ref}
        h=hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()
        with (self.out/'evidence-chain.jsonl').open('a') as f:
            f.write(json.dumps({'payload':body,'sha256':h})+'\n');f.flush();os.fsync(f.fileno())
        self.count+=1;self.previous=h;return ref
    def stage(self,name,arguments,callback):
        previous=self.active;self.active=name;inp=self.emit(name+'-input',arguments)
        value=callback();out=self.emit(name+'-return',value)
        record={'name':name,'status':'COMPLETE','input':inp,'result':out}
        with (self.out/'operation-index.jsonl').open('a') as f:f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
        self.operations.append(record);self.active=previous;return value
    def join(self,name,a,b):
        self.emit(name,{'actual':a,'expected':b,'join':'exact typed structure','passed':exact_structure(a,b)})
        require(exact_structure(a,b),name)
    def zero(self,name,a,b):
        previous=self.active;self.active=name;self.emit(name+'-input',{'left':a,'right':b})
        raw=a-b;v=sp.cancel(sp.together(raw));self.emit(name+'-return',{'raw':raw,'cancelled':v})
        require(v==0,name);self.active=previous;return v


def verify_gate(path,manifest_path,m):
    gate=json.loads(Path(path).read_text())
    require(gate['status']=='READY_FOR_ONE_END_UNIFORM_COMPARISON','actual readiness gate')
    require(gate['workerSha256']==sha(__file__) and gate['manifestSha256']==sha(manifest_path),'worker/manifest')
    require(gate['sourcePins']==m['sourcePins'],'source pins')
    for p,h in m['sourcePins'].items():require(sha(p)==h,'pin '+p)
    require(gate['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and gate['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'actual helper routes')
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate '+key)
    require(gate['launcher']==m['launcher'] and gate['buildReviewRecord']==m['reviewRecordWillBe'] and gate['authority']==m['executionAuthority'],'actual manifest routes')
    review=json.loads(Path(gate['buildReviewRecord']).read_text())
    require(review['independentBuildClearance'] is False and gate['independentBuildClearance'] is False and review['allChecksPassed'] is True,'literal review disposition')
    require(review['reports']['claude']['literalVerdict']=='NEEDS REVISION' and review['reports']['grok']['literalVerdict']=='CLEAR FOR THIS BOUNDED END-UNIFORM BUILD','preserved build verdicts')
    repair=json.loads(Path(gate['toolingRepairRecord']).read_text())
    require(gate['toolingRepairRecord']==m['toolingRepairRecord'] and sha(gate['toolingRepairRecord'])==m['sourcePins'][gate['toolingRepairRecord']]==gate['toolingRepairRecordSha256'],'exact tooling repair record')
    require(gate['localToolingExecutionAuthority'] is True and repair['localToolingExecutionAuthority'] is True and repair['independentBuildClearance'] is False,'standing tooling authority without author CLEAR')
    require(repair['reviewRecordSha256']==gate['buildReviewRecordSha256'] and repair['reviewedPacketSha256']==review['packetSha256'],'actual preserved review relation')
    require(repair['workerSha256']==gate['workerSha256'] and repair['launcherSha256']==gate['launcherSha256'] and repair['dimensionPredicateASTUnchanged'] is True and repair['scientificComparisonASTUnchanged'] is True,'tested exact metadata-loading correction')
    require(repair['testsPassed']>0 and repair['equationsChanged'] is False and repair['scienceRunsSoFar']==0,'tooling-only preparation')
    method=json.loads(Path(m['methodRecord']).read_text())
    require(method['jointIndependentMethodClearance'] is True and method['methodSha256']==sha(m['methodPath'])==review['methodSha256'],'actual cleared method')
    authority=json.loads(Path(gate['authority']).read_text())
    require(m['sourcePins'][gate['authority']]==gate['authoritySha256'],'pinned exact authority')
    require(authority['boundedInstrumentAuthorized'] is True and authority['scienceExecutionsAuthorized']==1 and authority['automaticScientificRetry'] is False and authority['noDeadline'] is True and authority['scope']==m['scope']==gate['scope'],'bounded standing authority')
    require(gate['durationLimits'] is None and gate['scientificRunsAuthorized']==1,'no deadline or retry')
    return gate


def polynomial(expr,x,y):
    if not expr.has(x,y):return {(0,0):expr}
    if expr==x:return {(1,0):sp.S.One}
    if expr==y:return {(0,1):sp.S.One}
    if expr.is_Add:
        out={}
        for term in expr.args:
            for key,value in polynomial(term,x,y).items():out[key]=out.get(key,0)+value
        return out
    if expr.is_Mul:
        out={(0,0):sp.S.One}
        for term in expr.args:out=convolution(out,polynomial(term,x,y),False)
        return out
    if expr.is_Pow and expr.exp.is_Integer is True and 0<=expr.exp<=128:
        base=polynomial(expr.base,x,y);out={(0,0):sp.S.One}
        for _ in range(int(expr.exp)):out=convolution(out,base,False)
        return out
    raise ExpectedUnavailable('unsupported polynomial dependence; no EX domain inference')


def coefficients(value,eta,sigma,J,name):
    J.emit(name+'-input',{'value':value,'eta':eta,'sigma':sigma})
    numerator,denominator=sp.fraction(sp.together(value))
    try:n=polynomial(numerator,eta,sigma);d=polynomial(denominator,eta,sigma)
    except ExpectedUnavailable as exc:
        out={'status':'ATTRIBUTION_UNAVAILABLE','reason':str(exc),'numerator':numerator,'denominator':denominator};J.emit(name+'-return',out);return out
    d0=sp.cancel(d.get((0,0),0))
    J.emit(name+'-domain',{'numerator':numerator,'denominator':denominator,'originDenominator':d0,'regularityCondition':sp.Ne(d0,0,evaluate=False)})
    if d0==0:
        out={'status':'REGULARITY_UNAVAILABLE','originDenominator':d0};J.emit(name+'-return',out);return out
    c={}
    for i,j in G:
        accum=sum(v*c[i-a,j-b] for (a,b),v in d.items() if (a,b)!=(0,0) and a<=i and b<=j and (i-a,j-b) in c)
        c[i,j]=sp.cancel((n.get((i,j),0)-accum)/d0)
    retained=sum(eta**i*sigma**j*c[i,j] for i,j in G);H=sp.cancel(value-retained)
    J.zero(name+'-reconstruction',value,retained+H)
    # Exact coefficient certificate for denominator*(value-retained), not a claim H=0.
    hc=convolution(d,c,False)
    for g in G:J.zero(name+'-excluded-grade-%d%d'%g,n.get(g,0)-hc.get(g,0),sp.S.Zero)
    out={'status':'EXTRACTED_ON_DECLARED_REGULAR_DOMAIN','coefficients':c,'retained':retained,'remainder':H,'originDenominator':d0}
    J.emit(name+'-return',out);return out


def finite_constant(expr):
    value=sp.simplify(expr)
    if value.free_symbols:return {'status':'UNRESOLVED_FREE_SYMBOLS','value':value}
    real,imag=sp.expand_complex(value).as_real_imag();real=sp.simplify(real);imag=sp.simplify(imag)
    reconstruction=sp.simplify(value-real-sp.I*imag)
    finite=real.is_finite is True and imag.is_finite is True and reconstruction==0
    nonzero=finite and (real.is_zero is False or imag.is_zero is False)
    return {'status':'NONZERO_CERTIFIED' if nonzero else ('ZERO' if finite and real==imag==0 else 'UNRESOLVED'),
            'value':value,'real':real,'imaginary':imag,'finite':finite,'reconstruction':reconstruction}


def wave_test(value,p,q,cs,J,name):
    J.emit(name+'-input',{'value':value,'wave':q**2-(9/cs**2-p**2-sp.Rational(1,20))})
    n,d=sp.fraction(sp.together(value))
    try:terms=polynomial(n,q,sp.Dummy('unused_q_coefficient'))
    except ExpectedUnavailable as exc:
        out={'status':'UNRESOLVED_UNSUPPORTED_DEPTH_DEPENDENCE','reason':str(exc),'value':value,'numerator':n,'denominator':d};J.emit(name+'-return',out);return out
    require(all(b==0 for a,b in terms),'only physical depth reduced')
    radical=9/cs**2-p**2-sp.Rational(1,20);a=sp.S.Zero;b=sp.S.Zero
    for (order,_),v in terms.items():
        if order%2:b+=v*radical**(order//2)
        else:a+=v*radical**(order//2)
    a=sp.cancel(a);b=sp.cancel(b);remainder=a+b*q
    # Synthetic division builds an explicit quotient; no full-domain Poly.
    coeff={i:sp.cancel(v) for (i,_),v in terms.items()};quo={}
    for i in range(max(coeff,default=0),1,-1):
        leading=coeff.get(i,0);quo[i-2]=leading;coeff[i]=0;coeff[i-2]=coeff.get(i-2,0)+leading*radical
    quotient=sum(v*q**i for i,v in quo.items());J.zero(name+'-division',n,quotient*(q*q-radical)+remainder)
    domain=sp.Ne(d,0,evaluate=False);point={p:sp.S.One,q:sp.Integer(2),cs:sp.sqrt(sp.Rational(180,101))}
    witness={'point':point,'denominator':finite_constant(d.subs(point)),'value':finite_constant(value.subs(point))}
    zero=a==0 and b==0
    state='ZERO' if zero else ('NONZERO_CERTIFIED' if witness['denominator']['status']=='NONZERO_CERTIFIED' and witness['value']['status']=='NONZERO_CERTIFIED' else 'UNRESOLVED_NONZERO_SYMBOLIC_REMAINDER')
    out={'status':state,'numerator':n,'denominator':d,'domain':domain,'quotient':quotient,'remainder':remainder,
         'coefficient0':a,'coefficient1':b,'positiveSheetWitness':witness,'symbolicCoverage':'on wave, outgoing sheet, original denominator nonzero; no sample-to-interval inference'}
    J.emit(name+'-return',out);return out


def inspect_matrix(matrix,p,q,cs,J,name):
    return [[wave_test(matrix[i,j],p,q,cs,J,name+'-%d-%d'%(i,j)) for j in range(matrix.cols)] for i in range(matrix.rows)]


def run_science(m,J,ns,codec):
    D=ns['decode'];copied={};raw={}
    for alias,pin in m['savedInputs'].items():
        source=Path(pin['path']);require(sha(source)==pin['sha256'],'saved pin '+alias)
        target=J.out/'saved'/alias;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        require(sha(target)==pin['sha256'],'saved copy '+alias)
        copied[alias]={'source':str(source),'path':str(target.relative_to(J.out)),'sha256':pin['sha256'],'bytes':target.stat().st_size}
        raw[alias]=json.loads(target.read_text())
    save(J.out/'saved-copy-index.json',copied)
    # Restore only complete old operation returns and their exact operands.
    receipts=m['uniformReceipts'];database=Path(m['uniformDatabase']);conn=sqlite3.connect('file:'+str(database)+'?mode=ro',uri=True)
    blobs={};opaque_copies={}
    try:
        for rec in receipts['blobReceipts']:
            payload,h,size=conn.execute('select payload,sha256,bytes from blobs where name=?',(rec['member'],)).fetchone()
            require(hashlib.sha256(payload).hexdigest()==h==rec['sha256'] and len(payload)==size==rec['bytes'],'opaque old blob '+rec['member'])
            target=J.out/'uniform-raw'/rec['member'];target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as f:f.write(payload);f.flush();os.fsync(f.fileno())
            opaque_copies[rec['member']]={'path':str(target.relative_to(J.out)),'sha256':h,'bytes':size}
            value=codec(payload);blobs[rec['member']]=value
            J.emit('restore-opaque-%03d'%len(blobs),{'receipt':rec,'copy':str(target.relative_to(J.out)),'decodedType':str(type(value)),'oldFunctionCalled':False})
    finally:
        conn.close();save(J.out/'opaque-copy-index.json',opaque_copies)
    operations={r['name']:r for r in receipts['completedOperations']}
    require(len(operations)==len(receipts['completedOperations'])==18,'exact completed operation census')
    def operation(name):
        r=operations[name];return blobs[r['input']['member']],blobs[r['result']['member']]
    restored={}
    for name,pin in m['originalPackets'].items():
        inp,value=operation('restore/'+name);J.join('restore-pin-'+name,inp['pin'],pin)
        require(inp['originalBytes']['sha256']==pin['sha256'] and inp['originalBytes']['bytes']==pin['bytes'],'actual original bytes '+name)
        require(sha(pin['path'])==pin['sha256'],'original source stays pinned '+name);restored[name]=value
    physical=json.loads(Path(m['physicalInput']).read_text());context=D(raw['ends/local-context.json']);cells=D(raw['ends/symbols.json']);depth=D(raw['ends/depth-domain.json'])
    require(context['physical']==physical and context['fullAssembly']['nativeEpsilonOnce'] is True and len(cells)==200,'inherited physical/epsilon/cell context')
    params={k:sp.Rational(v) for k,v in physical['parameters'].items()};params['omega']=sp.Integer(3)
    p,q,cs=depth['p'],depth['q'],depth['cs'];eta,sigma=context['context']['independentGrades']
    require((eta.name,sigma.name)==('eta_bg','sigma_W'),'actual independent grades')
    require(params['eta_bg']==sp.Rational(1,100) and params['W_0']==1 and params['L_W']==10,'actual finite origin and profile')
    require(raw['controls/conclusion.json']['status']=='SOURCE_JOINED_TRANSLATED_RETAINED_WEAK_ENDS' and raw['controls/conclusion.json']['endSymbolCells']==200,'saved concluded cells')
    J.zero('physical-wave-form',depth['radicand'],9/cs**2-p**2-sp.Rational(1,20))
    beta=depth['beta'];br,bi=sp.expand_complex(beta).as_real_imag()
    J.emit('actual-closed-depth-domain',{'source':depth,'betaReal':br,'betaImaginary':bi,'argument':'For outgoing real q>=0, Re(q+beta)>=Re(beta)>0; for q=i*y with y>=0, Re(q+beta)=Re(beta)>0.','newRootOrLimit':False})
    require(br==sp.Rational(30,109) and bi==sp.Rational(9,109) and depth['range']==[1,2],'actual beta/domain arguments')
    lift_input,L=operation('selected-lift');U=restored['uniformSource'];common=restored['uniformCommon'];response=restored['uniformResponse']
    profile=blobs['profiles/constant-end-and-zero-jets']
    with SourceStage(J,'native-source-chart-units'):
        S,side_map=source_chart_and_scale(m,J,D,raw,physical,context,depth,U,L,profile)
        # Native coordinate covariance is new missing evidence, never a producer replay.
        tex={'ast':ast,'hashlib':hashlib,'Path':Path}
        exec(compile(definitions(Path(m['textHelperSource']).read_text(),('require','text_sha','literal_record','tuple_arguments','literal_key','named','selected_case')),'inert-source-text-only','exec'),tex)
        original_text,provenance=tex['literal_record'](Path(m['bSource']),'slab_operator')
        _,case_text=tex['selected_case'](original_text,('LAB_HELD','RHO4_CONSTANT'));value_text=tex['named'](case_text,'VALUE')
        actual_rows={}
        for key in ('U_BODY_BALANCE','THETA_BALANCE','E_W_BALANCE'):
            part=tex['named'](tex['named'](value_text,key),'EXPANDED')
            if key=='U_BODY_BALANCE':actual_rows.update({'U'+str(i):v for i,v in enumerate(tex['tuple_arguments'](part))})
            else:actual_rows[key]=part
        for row in ROWS:
            require(raw['native/'+row+'.json']['source']==provenance and raw['native/'+row+'.json']['fullConstructor']==actual_rows[row],'actual native raw source row '+row)
        J.emit('native-source-text-joins',{'provenance':provenance,'rowHashes':{k:hashlib.sha256(v.encode()).hexdigest() for k,v in actual_rows.items()},'sourceConstructorCalled':False})
        contracts=D(raw['native/wave-profile.json'])
        c2_text=Path(m['c2Source']).read_text();c2_tree=ast.parse(c2_text)
        assignments={t.id:n.value for n in c2_tree.body if isinstance(n,ast.Assign) for t in n.targets if isinstance(t,ast.Name)}
        require(contracts['dimensions']==ast.literal_eval(assignments['DIMENSION_SCHEMA']),'actual native unit schema')
        wave_node=next(n for n in c2_tree.body if isinstance(n,ast.FunctionDef) and n.name=='wave_jet')
        require(contracts['waveSource']==ast.get_source_segment(c2_text,wave_node),'actual native wave/derivative convention')
        native={row:D({'text':'native source operand','srepr':raw['native/'+row+'.json']['fullConstructor']}) for row in ROWS}
        atoms=set().union(*(v.free_symbols for v in native.values()));byname={s.name:s for s in atoms};require(len(byname)==len(atoms),'unique native assumptions')
        scalar_names=set(physical['parameters'])|{'epsilon_shape','sigma_W','delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus'}
        rotation={}
        for atom in atoms:
            name=rotated_name(atom.name,scalar_names);require(name in byname,'mapped native atom present '+name)
            target=byname[name];require(atom.assumptions0==target.assumptions0,'rotation assumptions');rotation[atom]=target
        rowmap={'U0':'U2','U1':'U0','U2':'U1','THETA_BALANCE':'THETA_BALANCE','E_W_BALANCE':'E_W_BALANCE'}
        J.emit('native-constructor-footprint',{'constructorBytes':{r:len(actual_rows[r].encode()) for r in ROWS},'totalConstructorBytes':sum(len(v.encode()) for v in actual_rows.values()),'peakRssKiBBeforeCovariance':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'prospectiveMemoryGuarantee':False})
        J.emit('source-chart-operands',{'weakMomentum':(p,sp.Rational(1,5),sp.Rational(1,10)),'uniformMomentum':(sp.Rational(1,5),sp.Rational(1,10),p),'rotation':rotation,'rowMap':rowmap,'sourceRowCopies':{row:copied['native/'+row+'.json'] for row in ROWS},'nativeContract':contracts,'uniformChart':source_method(Path(m['engineSource']).read_text(),'EdgeReduction','__init__')})
        for row in ROWS:
            try:J.zero('source-cyclic-covariance-'+row,native[row].xreplace(rotation),native[rowmap[row]])
            except ValueError as exc:raise SourceMapUnresolved('native covariance '+row) from exc
        S=sp.ImmutableMatrix([[0,0,1,0,0],[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],[0,0,0,0,1]])
        require(S.det()==1 and S.T*S==sp.eye(5),'proper orthonormal chart map')
        # Dimensions are checked on unbound native operands, not guessed after binding.
        registries=restored_unit_registries(m,J,raw,restored)
        dimensions=merge_native_dimensions(contracts['dimensions'],atoms,registries,J);zero_dim=(sp.S.Zero,)*3
        def dimension(expr):
            if expr.is_Number or expr==sp.I:return zero_dim
            if expr.is_Symbol:
                require(expr.name in dimensions,'native dimensions '+expr.name);return tuple(map(sp.Rational,dimensions[expr.name]))
            if expr.is_Add:
                ds=[dimension(v) for v in expr.args if v!=0];require(not ds or all(d==ds[0] for d in ds),'native additive dimensions');return ds[0] if ds else zero_dim
            if expr.is_Mul:
                ds=[dimension(v) for v in expr.args];return tuple(map(sum,zip(*ds)))
            if expr.is_Pow and expr.exp.is_Rational:return tuple(expr.exp*d for d in dimension(expr.base))
            raise ValueError('unsupported native dimension node '+str(type(expr)))
        rowunits=[dimension(native[r]) for r in ROWS];fieldunits=[tuple(dimensions[f]) for f in FIELDS]
        J.emit('source-units',{'nativeRows':rowunits,'nativeFields':fieldunits,'frame':physical['unit_frame'],'unitMap':'proper permutation only, scalar factor one from strong epsilon coefficient; no power map'})
        lift_input,L=operation('selected-lift');U=restored['uniformSource'];common=restored['uniformCommon'];response=restored['uniformResponse']
        native_hashes=[h for path,h in U['sourceFiles'].items() if path.endswith('/S11c_b_exports.py')]
        require(native_hashes==[sha(m['bSource'])],'same original physical slab export across old/new sources')
        J.join('lift-input-curl',lift_input['uniform'],U['curl']);J.join('lift-input-common',lift_input['common'],common);J.join('lift-input-parameters',lift_input['parameters'],params)
        require(lift_input['source']==sha(m['engineSource']),'actual original lift source version')
        L={**L,'fieldUnits':response['fieldUnits'],'currentUnit':response['currentUnit']}
        J.join('lift-source-evidence',L['source'],blobs['lift/source-joins'])
        require(not any(s.name in ('eta_bg','sigma_W') for s in L['lift'].free_symbols),'grade-independent lift')
        require(tuple(L['source']['fieldOrder'])==('u1','u2','u3','theta','eW'),'source field order')
        expected_units=[fieldunits[i] for i in (1,2,0,3,4)];require(list(map(tuple,L['fieldUnits']))==expected_units,'native field unit permutation')
        J.emit('source-scalar-normalization',{'strongMatrixSource':source_method(Path(m['engineSource']).read_text(),'ConstantEndPencil','strong_matrix'),'pairingSource':source_method(Path(m['engineSource']).read_text(),'ClosedCurrentPairing','construct'),'epsilonCount':1,'factor':1,'excludedPowerMap':'PLUS_ROW_POWER_MAP','weakDuality':raw['native/weak-duality.json']})
        require(context['fullAssembly']['sourceRows']==list(ROWS) and context['fullAssembly']['fields']==list(FIELDS),'native weak row field identities')
    # Restore old proof operands and source profile values rather than calculating limits.
    profile=blobs['profiles/constant-end-and-zero-jets']
    for item in profile['profileJoins']:
        require(item['residual']==0 and all(v==0 for v in item['jets']),'inherited profile returns')
        require(item['supplied']==U['profileEndpoints'][item['source'][2]],'actual profile endpoint operand')
    with SourceStage(J,'saved-scalar-normalization'):
        saved_normalization_joins(m,J,D,raw,cells,context,p,params)
    groups={side:{g:sp.zeros(5) for g in G} for side in ('minus','plus')};closed={side:{g:sp.zeros(5) for g in G} for side in groups}
    coverage=[];seen=set()
    for index,c in enumerate(cells):
        side,row,col,g=c['side'],ROWS.index(c['row']),FIELDS.index(c['field']),tuple(c['grade']);key=(side,row,col,g)
        require(key not in seen and side in groups and g in G,'unique end cell');seen.add(key)
        require(c['sumIdentity']['cancelled']==0 and c['sumIdentity']['left']==c['symbol'],'saved sum proof actual argument')
        require(not any(s.name=='epsilon_shape' for s in c['symbol'].free_symbols),'weak symbol already coefficient of single epsilon')
        require(c['weakPairingFactor']==2*sp.pi and c['epsilonPower']==(0 if c['symbol']==0 else 1),'epsilon and Fourier contract')
        groups[side][g][row,col]=c['symbol'];closed[side][g][row,col]=c['closedGrazingValue']
    require(len(seen)==200,'all end cells')
    results={}
    for end,side in side_map.items():
        inp,R=operation(end+'/restriction');freq=restored[end+'Frequency'];pairing,dims=restored[end+'Pairing'];r=pairing['result']
        J.join(end+'-restriction-source',inp['actualSource'],r['CLOSED_PENCIL_LEGS']);J.join(end+'-restriction-frequency',inp['frequency'],freq)
        J.join(end+'-restriction-lift',inp['lift'],L);J.join(end+'-restriction-parameters',inp['parameters'],params)
        inv=blobs[end+'/restriction/full-five-row-invariance'];algebra=blobs[end+'/restriction/physical-algebraic-joins']
        for key in ('P','D','lift','gram'):J.join(end+'-invariant-'+key,R[key],inv[key])
        J.join(end+'-invariant-raw',R['invariant'],inv['rawResidual']);J.join(end+'-invariant-onwave',R['onWaveInvariant'],inv['onWaveResidual'])
        require(all(v==0 for v in R['onWaveInvariant']),'inherited on-wave zero only')
        for key in ('P','algebraic','relation','wave','conversion'):J.join(end+'-algebraic-'+key,R[key],algebra[key])
        atoms_R=set().union(*(x.free_symbols for x in (R['P'],R['D'],R['wave'],R['lift'],R['conversion'])))
        def symbol(name):
            found=[s for s in atoms_R if s.name==name];require(len(found)==1,'actual old symbol '+name);return found[0]
        w,k,cold,qold=symbol('uniformFrequency'),symbol('uniformNormal'),symbol('uniformSoundSpeed'),symbol('uniformPhysicalDepth')
        common_map={w:sp.Integer(3),k:p,cold:cs,qold:q}
        old=lambda v:v.xreplace(common_map)
        Pold,Dold,lift,Iold=map(old,(R['P'],R['D'],R['lift'],R['invariant']))
        Egrades={g:S.T*groups[side][g]*S for g in G};Ezero={g:S.T*closed[side][g]*S for g in G}
        origin=freq['origin'];require({s.name for s in origin}=={'eta_bg','sigma_W'},'old two independent grades')
        origins={s.name:v for s,v in origin.items()};require(origins=={'eta_bg':sp.Rational(1,100),'sigma_W':sp.Rational(1,1000)},'exact physical origin')
        finite_origin={eta:origins['eta_bg'],sigma:origins['sigma_W']}
        E=sum((eta**a*sigma**b*Egrades[a,b] for a,b in G),sp.zeros(5));Ephys=E.subs(finite_origin)
        Eg0=sum((eta**a*sigma**b*Ezero[a,b] for a,b in G),sp.zeros(5)).subs(finite_origin)
        # Raw binding is a missing new comparison; no old EndBinding method called.
        with SourceStage(J,end+'-native-source-binding-units'):
            source=r['CLOSED_PENCIL_LEGS'][0].xreplace(U['profileEndpoints'])
            wl,wr=r['FREQUENCY_LEGS'];kl,kr=r['NORMAL_LEGS'];ql,qr=r['BULK_LEGS']
            mapping={wl:sp.Integer(3),wr:sp.Integer(3),kr:p,qr:q}
            for s in source.free_symbols-set(mapping):
                if s.name=='eta_bg':mapping[s]=eta
                elif s.name=='sigma_W':mapping[s]=sigma
                elif s.name=='c_s0':mapping[s]=cs
                elif s.name in ('omega','s11cdFrequency'):mapping[s]=sp.Integer(3)
                elif s.name in params:mapping[s]=params[s.name]
                else:raise ValueError('raw source binding unavailable '+s.name)
            Praw=source.xreplace(mapping);J.emit(end+'-raw-source-binding',{'nativeSource':r['CLOSED_PENCIL_LEGS'][0],'endpoints':U['profileEndpoints'],'source':source,'map':mapping,'raw':Praw,'origin':finite_origin,'oldPhysical':Pold,'halfHeights':raw['controls/conclusion.json']['heights']})
            require(not any(s.name=='epsilon_shape' for s in Praw.free_symbols|r['CLOSED_PENCIL_LEGS'][0].free_symbols),'positive equation pencil already coefficient of single epsilon')
            require(not Praw.atoms(sp.Integral,sp.Derivative,sp.Limit,sp.Subs),'no unjoined original source operator')
            require(Praw.free_symbols<=set((p,q,cs,eta,sigma)),'complete raw binding')
            originjoin=inspect_matrix(Praw.subs(finite_origin)-Pold,p,q,cs,J,end+'/raw-origin')
            require(all(x['status']=='ZERO' for row in originjoin for x in row),'raw finite origin/source join')
            # Entry units derive from the unbound native rows and the inherited field frame.
            roworder=(1,2,0,3,4);native_entry_units=[[tuple(a-b for a,b in zip(rowunits[i],fieldunits[j])) for j in roworder] for i in roworder]
            require([[tuple(U['units']['strong'][(5*i+j,)]) for j in range(5)] for i in range(5)]==native_entry_units,'strong native row/field units')
            J.emit(end+'-unit-join',{'native':native_entry_units,'old':U['units']['strong'],'restoredRestrictionUnits':R['units'],'factor':1})
        action=Ephys*lift;A=(Ephys-Pold)*lift;Rnew=action-lift*Dold;Delta=Ephys-Pold
        for i in range(5):
            for j in range(2):J.zero(end+'-raw-AR-consistency-%d-%d'%(i,j),Rnew[i,j]-A[i,j],Iold[i,j])
        selected=inspect_matrix(A,p,q,cs,J,end+'/selected-A');residual=inspect_matrix(Rnew,p,q,cs,J,end+'/selected-R');full=inspect_matrix(Delta,p,q,cs,J,end+'/full-Delta')
        raw_action=Praw*lift;new_action=E*lift;grade=[]
        for i in range(5):
            for j in range(2):
                nm=end+'/grade-%d-%d'%(i,j);oldc=coefficients(raw_action[i,j],eta,sigma,J,nm+'-old');newc=coefficients(new_action[i,j],eta,sigma,J,nm+'-new')
                if not all(c['status']=='EXTRACTED_ON_DECLARED_REGULAR_DOMAIN' for c in (oldc,newc)):
                    grade.append({'row':i,'column':j,'status':'ATTRIBUTION_UNAVAILABLE','old':oldc,'new':newc});continue
                joins={g:newc['coefficients'][g]-oldc['coefficients'][g] for g in G}
                jt={g:wave_test(v,p,q,cs,J,nm+'-J-%d%d'%g) for g,v in joins.items()}
                attrib=sum(eta**a*sigma**b*joins[a,b] for a,b in G)-oldc['remainder']
                attribution=wave_test(A[i,j]-attrib.subs(finite_origin),p,q,cs,J,nm+'-attribution')
                if not attribution_usable(attribution['status'],newc['remainder']==0):
                    grade.append({'row':i,'column':j,'status':'ATTRIBUTION_UNAVAILABLE','finite':selected[i][j]['status'],'waveAttribution':attribution,'newLawRemainder':newc['remainder'],'oldRemainder':oldc['remainder']});continue
                retained='AGREEMENT' if all(v['status']=='ZERO' for v in jt.values()) else ('RETAINED_MISMATCH' if any(v['status']=='NONZERO_CERTIFIED' for v in jt.values()) else 'UNRESOLVED')
                finite=selected[i][j]['status'];classification='FINITE_TRUNCATION_DIFFERENCE' if retained=='AGREEMENT' and finite=='NONZERO_CERTIFIED' else retained
                grade.append({'row':i,'column':j,'retained':retained,'finite':finite,'classification':classification,'J':jt,'H':oldc['remainder'],'regularity':oldc['originDenominator']})
        for index,c in enumerate(cells):
            if c['side']!=side:continue
            wi=ROWS.index(c['row']);wj=FIELDS.index(c['field']);native_lift=S*lift
            weighted=[c['symbol']*native_lift[wj,j] for j in range(2)]
            coverage.append({'end':end,'cellIndex':index,'row':c['row'],'field':c['field'],'grade':c['grade'],'localAncestry':c['localAncestry'],'pressureAddresses':c['addressIds'],'selectedWeights':[native_lift[wj,j] for j in range(2)],'weightedSymbol':weighted,'zeroWeight':all(native_lift[wj,j]==0 for j in range(2)),'pressureSum':c['pressureSum'],'localSum':c['localSum']})
        # Closed grazing uses saved p and speed targets. Do not call any limit.
        grazing=[]
        for sign in (-1,1):
            args,limits=operation(end+'/exact-limits/'+str(sign));J.emit(end+'-limit-%d-actual-input'%sign,args)
            J.join(end+'-limit-%d-restriction'%sign,args['restriction'],R)
            native_args,native_return=operation(end+'/native-reconstruction')
            J.join(end+'-limit-%d-native'%sign,args['native'],native_return)
            J.join(end+'-native-%d-pairing'%sign,native_args['pairing'],pairing)
            J.join(end+'-native-%d-original'%sign,native_args['native'],restored[end+'NativeSource'])
            require(args['sign']==sign and native_args['frequency']==3,'actual limit/native arguments')
            require(limits['status']=='COMPUTED','saved limits completed')
            normal=old(limits['normal']);speed=old(R['cT']);mp={p:normal,cs:speed,q:sp.S.Zero}
            J.zero(end+'-grazing-wave-%d'%sign,(9/cs**2-p**2-sp.Rational(1,20)).subs(mp),sp.S.Zero)
            require(limits['normal']==sign*sp.sqrt(R['k23']),'actual inherited normal target')
            C=Eg0.subs(mp);LL=lift.subs(mp);DD=Dold.subs(mp);g=C*LL-LL*DD
            entries=[[finite_constant(g[i,j]) for j in range(2)] for i in range(5)]
            regular=all(finite_constant(x).get('finite',False) for obj in (C,LL,DD) for x in obj)
            paths={name:v['selected']['pencil'] for name,v in limits['paths'].items()}
            require(set(paths)=={'RADIATING','EVANESCENT'},'both actual limit paths')
            for route in paths:
                ev=blobs[end+'/limits/sign-'+str(sign)+'/'+route+'/wave-path']
                require(all(x==0 for x in ev['waveResiduals']),'saved actual wave path zero')
                t=ev['parameter'];require(t.is_positive is True,'actual approach parameter')
                require(ev['mapping'][qold]==(t if route=='RADIATING' else sp.I*t),'actual outgoing approach')
                require(ev['mapping'][k]==limits['normal'] and ev['mapping'][w]==3,'actual frequency/normal approach')
                J.emit(end+'-grazing-%d-'%sign+route+'-path',ev)
            require(all(v==0 for v in limits['selectedJoin']['pencil']),'inherited two-sided selected pencil join')
            J.emit(end+'-grazing-%d'%sign,{'point':mp,'closedE':C,'lift':LL,'D':DD,'residual':g,'entryCertificates':entries,'regular':regular,'savedSelectedPencilLimits':paths,'savedSelectedJoin':limits['selectedJoin'],'gradeAttribution':'UNAVAILABLE_NO_NEW_REMAINDER_LIMIT_COMPUTED','oldFunctionsCalled':False})
            grazing.append({'sign':sign,'point':mp,'status':'COMPARED' if regular else 'UNRESOLVED_OPERAND_FINITE_DOMAIN','R0':entries,'grazingGradeAttribution':'UNAVAILABLE'})
        # Three controls use full changed-minus-baseline row movements.
        point={p:sp.S.One,q:sp.Integer(2),cs:sp.sqrt(sp.Rational(180,101))};controls=[]
        choices=[(g,i,j,(eta**g[0]*sigma**g[1]*Egrades[g][i,j]).subs(finite_origin)) for g in G for i in range(5) for j in range(5) if Egrades[g][i,j]!=0 and any(lift[j,c]!=0 for c in range(2))]
        def record_control(name,baseline,changed,ancestry):
            movement=(changed-baseline).subs(point);certs=[finite_constant(x) for x in movement]
            baseline_certs=[finite_constant(x.subs(point)) for x in baseline];changed_certs=[finite_constant(x.subs(point)) for x in changed]
            responsive=any(c['status']=='NONZERO_CERTIFIED' for c in certs) and all(c.get('finite',False) for c in certs+baseline_certs+changed_certs)
            rec={'name':name,'point':point,'baseline':baseline,'changed':changed,'movement':movement,'certificates':certs,'baselineCertificates':baseline_certs,'changedCertificates':changed_certs,'responsive':responsive,'ancestry':ancestry,'meaning':'sensitivity only'};J.emit(end+'/control-'+name,rec);controls.append(rec);return responsive
        chosen=None
        for g,i,j,term in choices:
            move=sp.zeros(5,2)
            for c in range(2):move[i,c]=-term*lift[j,c]
            if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in move):chosen=(g,i,j,term,move);break
        if chosen:
            g,i,j,term,move=chosen
            original=next(c for c in cells if c['side']==side and c['row']==ROWS[roworder[i]] and c['field']==FIELDS[roworder[j]] and tuple(c['grade'])==g)
            record_control('omit-cell',Rnew,Rnew+move,{'row':i,'field':j,'grade':g,'removedTerm':term,'originalSavedCell':original})
        else:controls.append({'name':'omit-cell','responsive':False,'status':'INAPPLICABLE'})
        bad=lift.xreplace({p:-p});responsive=record_control('reverse-lift-normal',Rnew,Ephys*bad-bad*Dold,{'normalSignHeldInEquation':True})
        if not responsive:
            fallback=False
            for i,j in ((i,j) for i in range(5) for j in range(2) if lift[i,j]!=0):
                bad=sp.MutableDenseMatrix(lift);bad[i,j]=0;movement=Ephys*bad-bad*Dold-Rnew
                if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in movement):
                    record_control('curl-entry-ablation',Rnew,Ephys*bad-bad*Dold,{'row':i,'column':j,'directionProof':False});fallback=True;break
            if not fallback:
                rec={'name':'curl-entry-ablation','responsive':False,'status':'INAPPLICABLE_NO_CERTIFIED_MOVEMENT','searchedEntries':[(i,j) for i in range(5) for j in range(2) if lift[i,j]!=0]};J.emit(end+'/control-curl-entry-ablation',rec);controls.append(rec)
        sheet_choice=None
        for j in (3,4):
            baseline=Ephys[:,j];changed=Ephys.xreplace({q:-q})[:,j]
            applicable=all(not Egrades[g][i,j].has(q) or any(c['side']==side and c['field']==FIELDS[roworder[j]] and c['pressureSum']!=0 for c in cells) for g in G for i in range(5))
            if not applicable:
                J.emit(end+'/control-opposite-sheet-inapplicable-column-'+str(j),{'column':j,'status':'UNRESOLVED_PRESSURE_ANCESTRY','baseline':baseline,'changed':changed});continue
            if any(finite_constant(x.subs(point))['status']=='NONZERO_CERTIFIED' for x in changed-baseline):sheet_choice=(j,baseline,changed);break
        if sheet_choice:
            j,baseline,changed=sheet_choice;record_control('opposite-sheet',baseline,changed,{'unitColumn':j,'isMode':False,'pressureAncestry':[c for c in cells if c['side']==side and c['field']==FIELDS[roworder[j]] and c['pressureSum']!=0]})
        else:controls.append({'name':'opposite-sheet','responsive':False,'status':'INAPPLICABLE'})
        result={'fullDifference':full,'selectedFinite':selected,'fiveRowResidual':residual,'grades':grade,'grazing':grazing,'controls':controls,'currentTransfer':False,'comparisonDoesNotDischargeZeroWeightedPressureCells':True}
        J.emit(end+'-comparison',result);results[end]={'retainedStatuses':sorted(set(g.get('retained',g.get('status','UNRESOLVED')) for g in grade)),'finiteStatuses':sorted(set(v['status'] for row in selected for v in row)),'controlCoverage':{c['name']:c['responsive'] for c in controls}}
    J.emit('selected-cell-coverage',coverage)
    return {'executionStatus':'COMPLETE_COMPARISON_PENDING_INSPECTION','ends':results,'cells':len(cells),'restoredOperations':len(operations),'newRoots':False,'fieldOrLoss':False,'frequencyDerivativeTransfer':False,'oldFunctionsCalled':False}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    m=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,m)
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==expected and g['command'][-len(expected):]==expected and args.out.resolve()==Path(g['outputDirectory']).resolve(),'actual argv/output')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;code=1;result={};start=time.monotonic()
    pins={**m['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(m['helperSource']).read_text(),('require','containment','decode')),'inert-contained-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp,np,FunctionClass
        import sympy as sp
        import numpy as np
        from sympy.core.symbol import Str
        from sympy.core.function import FunctionClass
        ns.update(sp=sp,Str=Str)
        codec_ns={'pickle':pickle,'OrderedDict':OrderedDict,'builtins':builtins,'np':np,'sp':sp,'importlib':importlib,'io':io,'require':require}
        exec(compile(definitions(Path(m['uniformWorker']).read_text(),('SavedCodec','decode')),'saved-codec-only-no-old-guard','exec'),codec_ns)
        J=Evidence(args.out);result=J.stage('end-uniform-comparison',{'manifest':g['manifestSha256'],'build':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns,codec_ns['decode']));code=0
    except BaseException as exc:
        result={'executionStatus':'SOURCE_MAP_UNRESOLVED' if isinstance(exc,SourceMapUnresolved) else 'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        for index_name in ('saved-copy-index.json','opaque-copy-index.json'):
            if (args.out/index_name).exists():
                for row in json.loads((args.out/index_name).read_text()).values():pins[str(args.out/row['path'])]=row['sha256']
        posts={}
        for p,h in pins.items():
            try:posts[p]={'expected':h,'actual':sha(p),'error':None}
            except OSError as e:posts[p]={'expected':h,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',posts)
        if any(r['expected']!=r['actual'] for r in posts.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        if J is not None:result=J.encode(result)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
        if J is not None:save(args.out/'evidence-final-receipt.json',{'records':J.count,'lastSha256':J.previous,'chainSha256':sha(args.out/'evidence-chain.jsonl'),'completeProcess':code==0})
    return code

if __name__=='__main__':sys.exit(main())
