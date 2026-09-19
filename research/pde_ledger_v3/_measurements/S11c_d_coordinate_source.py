#!/usr/bin/env python3
"""Reversible reduced-source chart operands and literal shape mutations."""
import argparse,ast,contextlib,json,resource,signal,time,tokenize
from pathlib import Path
import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_continuum_grades as grades
import S11c_d_continuum_boundary as boundary
f=grades.f;engine=f.engine
PLAN=f.M/'S11c_d_coordinate_source_plan.md';PREFIX='COORDINATE_SOURCE_LAB_HELD_RHO4_CONSTANT'


def selected_exports(path,wanted):
    """Stream the actual top-level ledger keys and evaluate only selected rows."""
    header=[]
    with path.open() as stream:
        for line in stream:
            if line.startswith('_LEDGER'):
                break
            header.append(line)
    tree=ast.parse(''.join(header));allowed=[]
    for node in tree.body:
        if isinstance(node,(ast.Import,ast.ImportFrom)) or (isinstance(node,ast.FunctionDef) and node.name=='_restore') or (isinstance(node,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='_RELATIONALS' for x in node.targets)):
            allowed.append(node)
    namespace={};exec(compile(ast.Module(body=allowed,type_ignores=[]),str(path),'exec'),namespace)
    keys=[];selected={};active=False;await_brace=False;depth=0;expect_key=False;key=None;collect=None
    with path.open() as stream:
        for tok in tokenize.generate_tokens(stream.readline):
            text=tok.string
            if not active:
                if tok.type==tokenize.NAME and text=='_LEDGER':await_brace=True
                elif await_brace and text=='{':active=True;depth=1;expect_key=True
                continue
            if depth==1 and expect_key and tok.type==tokenize.STRING:
                key=ast.literal_eval(text);f.require(key not in keys,'injective actual export keys');keys.append(key);expect_key=False;continue
            if depth==1 and text==':' and key is not None:
                collect=[] if key in wanted else None;continue
            if depth==1 and text in (',','}'):
                if collect is not None:
                    code=tokenize.untokenize(collect);selected[key]=eval(compile(ast.parse(code.strip(),mode='eval'),str(path),'eval'),namespace)
                collect=None;key=None;expect_key=True
                if text=='}':break
                continue
            if collect is not None:collect.append((tok.type,text))
            if tok.type==tokenize.OP:
                if text in ('{','[','('):depth+=1
                elif text in ('}',']',')'):depth-=1
    f.require(active and depth==1 and keys,'complete source ledger dictionary')
    return selected,keys


def load(base):
    packet,cp,path=f.accepted_packet(f.M/'S11c_d_continuum_grade_checkpoint.json','continuum-grades.pickle')
    rp=next(Path(n) for n in packet['inputPackets'] if n.endswith('reduced-action.pickle'))
    f.require(f.digest(rp)==packet['inputPackets'][str(rp)],'native reduced source')
    r,dimensions=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(rp));dimensions.__dict__.update(packet['dimensionState'])
    fold={};audit=[]
    for name in ('S11c_b_exports.py','S11c_c1_exports.py','S11c_c2_exports.py'):
        selected,census=selected_exports(f.ROOT/'scripts'/name,{'background_density_map','L_W','omega'})
        fold.update(selected);audit.append({'path':'scripts/'+name,'keys':census})
    density=fold['background_density_map']['value']
    f.require(fold['L_W']['value']==r.ell and fold['omega']['value']==r.omega,'selective scalar/native input joins')
    pins=dict(cp['sourceFiles'])
    for p in (Path(__file__),PLAN,f.M/'S11c_d_continuum_grade_checkpoint.json',f.M/'S11c_d_uniform_response_checkpoint.json',Path(boundary.__file__),f.ACCEPTANCE):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    operands={str(path):f.digest(path),str(rp):f.digest(rp)}
    for n,h in pins.items():
        f.require(f.digest(f.ROOT/n)==h,('current consumed source',n));dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((f.ROOT/n).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'extraConsumedRow':'background_density_map','case':['LAB_HELD','RHO4_CONSTANT'],'parentFoldAudit':audit,'scope':'Affine kinematic chart and source mutation preparation; no second scattering route yet.'})
    return r,packet,density,pins,operands


def chart(r):
    x=sp.symbols('s11cdMaterialX1 s11cdMaterialX2 s11cdMaterialX3',real=True)
    ansatz=sp.Matrix([sp.Rational(6,5)*x[0]+x[2]/7,sp.Rational(4,5)*x[1],sp.Rational(5,4)*x[2]])
    F=ansatz.jacobian(x);inverse=F.inv();J=F.det();A=F[:2,:2].det();d=F[2,2];c=F[0,2]
    physical_k=sp.Matrix([*r.tangents,sp.Symbol('s11cdCoordinateNormalMomentum',real=True)])
    phase=(physical_k.T*ansatz)[0];material_k=sp.Matrix([sp.diff(phase,v) for v in x])
    B=sp.diag(F,sp.eye(2));piola=J*inverse
    # Tangential inverse Fourier measure: dk_E=dk_M/A; hence psi_M=exp(i phase_shear) B^-1 psi_E/A.
    shear=sp.expand(phase-material_k[0]*x[0]-material_k[1]*x[1]-physical_k[2]*d*x[2])
    T=A*sp.exp(-sp.I*shear)*B;U=sp.exp(sp.I*shear)*B.inv()/A
    dimensions=engine.PHYSICAL_METADATA.dimensions
    for v in x:dimensions.known[v]=(1,0,0)
    dimensions.known[physical_k[2]]=(-1,0,0)
    Xp=sp.Symbol('s11cdMaterialSourcePosition',real=True);Xi=sp.Symbol('s11cdMaterialProfileCoordinate',real=True)
    K=sp.symbols('s11cdMaterialTangentialMomentum1 s11cdMaterialTangentialMomentum2',real=True)
    definitions={r.z:(x[2],d*x[2]),r.zp:(Xp,d*Xp),r.xi:(Xi,d*Xi),r.tangents[0]:(K[0],K[0]/F[0,0]),r.tangents[1]:(K[1],K[1]/F[1,1])}
    for group,name in zip(r.momentum_groups,('Output','Input','Middle')):
        original=r.normal_map[group[2]];new=sp.Symbol('s11cdMaterial'+name+'NormalMomentum',real=True)
        definitions[original]=(new,(new-c*K[0]/F[0,0])/d)
    for old,(new,value) in definitions.items():dimensions.known[new]=dimensions.measure(old)
    forward={old:value for old,(new,value) in definitions.items()};backward={new:sp.solve(sp.Eq(value,old),new)[0] for old,(new,value) in definitions.items()}
    backward={new:value.xreplace({K[0]:F[0,0]*r.tangents[0],K[1]:F[1,1]*r.tangents[1]}) for new,value in backward.items()}
    proof={'inverse':F*inverse-sp.eye(3),'covector':material_k-F.T*physical_k,'field':(T*U-sp.eye(5)).applyfunc(sp.simplify),'normalCurrent':piola[2,:]-sp.Matrix([[0,0,A]]),'normalVolume':J-A*d}
    for old,(new,value) in definitions.items():proof['coordinate_'+str(old)]=sp.expand(value.xreplace(backward)-old)
    return {'coordinates':x,'ansatz':ansatz,'displacement':ansatz-sp.Matrix(x),'F':F,'inverse':inverse,'volumeJacobian':J,'tangentialJacobian':A,'normalJacobian':d,'shear':shear,'physicalCovector':physical_k,'materialCovector':material_k,'fieldComponents':B,'eulerianFromMaterialFourier':T,'materialFromEulerianFourier':U,'piola':piola,'proof':proof,'definitions':definitions,'forward':forward,'backward':backward,'sourcePosition':Xp,'profileCoordinate':Xi,'tangentialCoordinates':K}


def coordinate_change(expression,definitions,mapping,reverse=False):
    """Change actual bound integration/limit variables, including their Jacobians."""
    by_variable={new:(old,sp.diff(mapping[new],old)) for old,(new,value) in definitions.items()} if reverse else {old:(new,sp.diff(value,new)) for old,(new,value) in definitions.items()}
    def visit(node):
        if node in mapping:return mapping[node]
        if isinstance(node,sp.Integral):
            value=visit(node.function);limits=[]
            for limit in node.limits:
                old=limit[0]
                if old not in by_variable:limits.append(tuple(visit(v) for v in limit));continue
                new,jacobian=by_variable[old];value*=jacobian
                inverse=sp.solve(sp.Eq(mapping[old],old),new)[0]
                limits.append((new,*(sp.expand(inverse.subs(old,visit(v))) for v in limit[1:])))
            return sp.Integral(value,*limits)
        if isinstance(node,sp.Limit) and node.args[1] in by_variable:
            old=node.args[1];new,jacobian=by_variable[old];inverse=sp.solve(sp.Eq(mapping[old],old),new)[0]
            return sp.Limit(visit(node.args[0]),new,inverse.subs(old,visit(node.args[2])),dir=str(node.args[3]))
        if not node.args:return node
        return node.func(*(visit(v) for v in node.args))
    return visit(expression)


def jets(r,geometry,maximum):
    X=geometry['coordinates'][2];d=geometry['normalJacobian'];A=geometry['tangentialJacobian'];B=geometry['fieldComponents'];phase=geometry['shear'];slope=sp.diff(phase,X)
    forward={0:{0:A*B}};inverse={0:{0:B.inv()/A}}
    for n in range(1,maximum+1):
        forward[n]={};inverse[n]={}
        for j in range(n+1):
            a=forward[n-1].get(j,sp.zeros(5));prev=forward[n-1].get(j-1,sp.zeros(5))
            forward[n][j]=(a.diff(X)-sp.I*slope*a+prev)/d
            a=inverse[n-1].get(j,sp.zeros(5));prev=inverse[n-1].get(j-1,sp.zeros(5))
            inverse[n][j]=a.diff(X)+sp.I*slope*a+d*prev
    proof={}
    for n in range(maximum+1):
        for j in range(n+1):
            composed=sum((forward[n][k]*inverse[k][j] for k in range(j,n+1)),sp.zeros(5))
            proof[n,j]=(composed-(sp.eye(5) if n==j else sp.zeros(5))).applyfunc(sp.expand)
    # Independent differentiated polynomial test, with the physical polynomial evaluated at y=dX.
    y=sp.Symbol('s11cdCoordinatePolynomialPosition',real=True);engine.PHYSICAL_METADATA.dimensions.known[y]=(1,0,0)
    polynomial={}
    for degree in range(maximum+1):
        material=geometry['materialFromEulerianFourier']*(d*X)**degree
        for n in range(maximum+1):
            actual=sum((sp.exp(-sp.I*phase)*forward[n][j]*material.diff(X,j) for j in range(n+1)),sp.zeros(5))
            expected=sp.eye(5)*sp.diff(y**degree,y,n).subs(y,d*X)
            polynomial[degree,n]=(actual-expected).applyfunc(sp.simplify)
    return {'forwardWithoutPhase':forward,'inverseWithoutPhase':inverse,'phase':phase,'compositionResiduals':proof,'polynomialResiduals':polynomial,'maximumOrder':maximum}


def linear_source(expression,r,index,units):
    probes=sorted([v for v in expression.atoms(AppliedUndef) if v.func.__name__.startswith('s11cdPencilProbe')],key=str)
    f.require(len(probes)==1 and probes[0].args==(r.zp,),'one original source field')
    probe=probes[0];column=int(probe.func.__name__.removeprefix('s11cdPencilProbe'));degree=max([0]+[sum(n for _,n in v.variable_count) for v in expression.atoms(sp.Derivative) if v.expr==probe])
    fields=[sp.diff(probe,r.zp,n) for n in range(degree+1)];symbols=sp.symbols(f's11cdCoordinateSource{index}Jet0:{degree+1}')
    for n,s in enumerate(symbols):engine.PHYSICAL_METADATA.dimensions.known[s]=tuple(v-n*w for v,w in zip(units[column],(1,0,0)))
    encoded=expression.xreplace(dict(zip(fields,symbols)));carriers={};restore={}
    def compress(node):
        if not node.has(*symbols):
            if node.is_Number:return node
            if node not in carriers:
                atom=sp.Dummy(f'coordinateCoefficient{len(carriers)}');carriers[node]=atom;restore[atom]=node
            return carriers[node]
        if not node.args:return node
        return node.func(*(compress(v) for v in node.args))
    small=compress(encoded);poly=sp.Poly(small,*symbols)
    f.require(all(sum(monomial)==1 for monomial,_ in poly.terms()),'complete linear native source')
    compact=[poly.coeff_monomial(s) for s in symbols];coefficients=[v.xreplace(restore) for v in compact]
    residual=sp.expand(small-sum(c*s for c,s in zip(compact,symbols)));round_trip=small.xreplace(restore)-encoded
    f.require(round_trip==0,'exact jet-free carrier definitions')
    return {'probe':probe,'column':column,'fields':fields,'symbols':symbols,'encoded':encoded,'coefficients':coefficients,'reconstructionResidual':residual,'rawReconstructionResidual':encoded-sum(c*s for c,s in zip(coefficients,symbols)),'carrierPolynomial':small,'carrierDefinitions':restore,'carrierRoundTripResidual':round_trip,'original':expression}


def first_jet_mutation(expression,r):
    mapping={}
    for node in expression.atoms(sp.Subs):
        derivative=node.expr
        if isinstance(derivative,sp.Derivative) and derivative.expr.func==r.profiles['w'] and sum(n for _,n in derivative.variable_count)==1:mapping[node]=-node
    altered=expression.xreplace(mapping)
    return {'base':expression,'mutated':altered,'residual':altered-expression,'occurrences':tuple(mapping),'count':sum(v in mapping for v in sp.preorder_traversal(expression)),'distinctAtomCount':len(mapping)}


def construct(base,r,packet,density,progress):
    geometry=chart(r);f.atomic_pickle(base/'chart.pickle',geometry)
    f.require(all(all(x==0 for x in value) if isinstance(value,sp.MatrixBase) else value==0 for value in geometry['proof'].values()),'derived chart identities')
    maximum=max(v['address'][1] for v in packet['records'].values() if v['address'][0]=='local');source={}
    for item in packet['records'].values():
        if item['address'][0]=='source':
            index=item['address'][1];source[index]=linear_source(item['record']['ORIGINAL'],r,index,packet['fieldUnits']);maximum=max(maximum,len(source[index]['fields'])-1)
    jet=jets(r,geometry,maximum);f.atomic_pickle(base/'field-jets.pickle',jet)
    f.require(all(x==0 for values in (jet['compositionResiduals'],jet['polynomialResiduals']) for matrix in values.values() for x in matrix),'full derivative maps and polynomial controls')
    material_jets={}
    for column in range(5):
        for order in range(maximum+1):
            atom=sp.Symbol(f's11cdMaterialField{column}Jet{order}');material_jets[column,order]=atom
            engine.PHYSICAL_METADATA.dimensions.known[atom]=tuple(v-order*w for v,w in zip(packet['fieldUnits'][column],(1,0,0)))
    records={};inventory={};target=base/'records';target.mkdir()
    for key,item in packet['records'].items():
        original=item['record']['ORIGINAL'];address=item['address'];encoded=source[address[1]]['encoded'] if address[0]=='source' else original
        altered=coordinate_change(encoded,geometry['definitions'],geometry['forward']);restored=coordinate_change(altered,geometry['definitions'],geometry['backward'],True)
        residual=restored-encoded;residual=sp.expand(residual) if residual!=0 else residual
        mutation=first_jet_mutation(original,r)
        record={'address':address,'original':original,'encoded':encoded,'coordinateImage':altered,'coordinateReplay':restored,'coordinateResidual':residual,'sourceJets':source.get(address[1]) if address[0]=='source' else None,'shape':mutation,'unit':item['record']['UNIT']}
        if address[0]=='source':
            native=source[address[1]];column=native['column'];coefficient_images=[coordinate_change(a,geometry['definitions'],geometry['forward']) for a in native['coefficients']]
            phase=jet['phase'].subs(geometry['coordinates'][2],geometry['sourcePosition']).xreplace(geometry['forward'])
            material=sum(coefficient_images[n]*sp.exp(-sp.I*phase)*jet['forwardWithoutPhase'][n][j][column,k].xreplace(geometry['forward'])*material_jets[k,j] for n in range(len(coefficient_images)) for j in range(n+1) for k in range(5))
            record['materialSourceAmplitude']=material;record['materialSourceCoefficients']=coefficient_images
        path=target/(key+'.pickle');f.atomic_pickle(path,record);inventory[key]={'path':str(path.relative_to(base)),'sha256':f.digest(path),'bytes':path.stat().st_size};f.save(base/'record-inventory.json',inventory)
        f.require(residual==0 and (record['sourceJets'] is None or record['sourceJets']['reconstructionResidual']==0),'exact source-coordinate/linear reconstruction')
        records[key]=record;progress('record_'+key)
    entries=dict(density);case=next(v for k,v in entries.items() if tuple(map(str,k))==('RHO4_CONSTANT',));value=dict(case)[sp.core.symbol.Str('VALUE')]
    equalities=list(value[0]);equation=next(v for v in equalities if isinstance(v,sp.Equality) and str(v.lhs)=='rho_4D_bg_rho4_constant')
    gradient=sp.Matrix([sp.diff(equation.rhs,v) for v in geometry['coordinates']]);factor=(geometry['displacement'].T*gradient)[0]/equation.rhs
    absence={'sourceCase':case,'densityEquation':equation,'gradient':gradient,'displacement':geometry['displacement'],'factor':factor,'structuralNonzeroEntries':sum(x!=0 for x in gradient)}
    f.atomic_pickle(base/'density-advection.pickle',absence)
    limits=[]
    for term in packet['termJoins']:
        old=term['remainingLimits'];new=[]
        for limit in old:
            var=limit[0];target,value=geometry['definitions'][var];inverse=sp.solve(sp.Eq(value,var),target)[0];new.append((target,*(sp.expand(inverse.subs(var,v.xreplace(geometry['forward']))) for v in limit[1:])))
        source_limit=term['sourceLimit'];sv=source_limit[0];sm,simage=geometry['definitions'][sv];sinverse=sp.solve(sp.Eq(simage,sv),sm)[0]
        changed_source=(sm,*(sp.expand(sinverse.subs(sv,v.xreplace(geometry['forward']))) for v in source_limit[1:]))
        limits.append({'originalSourceLimit':source_limit,'materialSourceLimit':changed_source,'sourceJacobian':sp.diff(simage,sm),'row':term['row'],'column':term['column'],'term':term['term'],'integralIndex':term['integralIndex'],'originalLimits':old,'materialLimits':tuple(new)})
    f.atomic_pickle(base/'ordered-limits.pickle',limits)
    return {'chart':geometry,'fieldJets':jet,'records':records,'densityAdvection':absence,'orderedLimits':limits,'fieldUnits':packet['fieldUnits'],'equationUnits':packet['equationUnits'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))}



def emit_density_equation(name,equation):
    """Keep the equation carrier; derive metadata on its two physical sides."""
    f.require(isinstance(equation,sp.Equality),'actual density source equality')
    dimensions=engine.PHYSICAL_METADATA.dimensions
    unit=dimensions.measure(equation.rhs)
    f.require(unit is not None and all(not v.free_symbols for v in unit),'resolved source density unit')
    if equation.lhs in dimensions.known:
        f.require(dimensions.known[equation.lhs]==unit,'density equality side dimensions')
    else:
        dimensions.known[equation.lhs]=unit
    operands=sp.Tuple(equation.lhs,equation.rhs)
    engine.emit(PREFIX+'_'+name,engine.carrier_fingerprint(equation))
    engine.emit('METADATA_'+PREFIX+'_'+name,engine.PHYSICAL_METADATA.record(operands,{(0,):unit,(1,):unit}))


def emit_result(result):
    def put(name,value,unit):engine.fingerprinted(PREFIX+'_'+name,value,{p:unit(p) if callable(unit) else unit for p,_ in engine.leaves(engine.cas(value))})
    g=result['chart'];zero=(0,0,0)
    for name in ('F','inverse','volumeJacobian','tangentialJacobian','normalJacobian','fieldComponents','piola'):put('CHART_'+name,g[name],zero)
    for name in ('ansatz','displacement'):put('CHART_'+name,g[name],(1,0,0))
    for name in ('physicalCovector','materialCovector'):put('CHART_'+name,g[name],(-1,0,0))
    for name in ('eulerianFromMaterialFourier','materialFromEulerianFourier'):put('CHART_'+name,g[name],zero)
    for name,value in g['proof'].items():put('CHART_PROOF_'+name,value,(-1,0,0) if name=='covector' or 'Momentum' in name else (1,0,0) if 'Position' in name else zero)
    jet=result['fieldJets']
    for n,table in jet['forwardWithoutPhase'].items():
        for j,value in table.items():put(f'FIELD_JET_{n}_{j}',value,(j-n,0,0));put(f'FIELD_JET_INVERSE_{n}_{j}',jet['inverseWithoutPhase'][n][j],(j-n,0,0));put(f'FIELD_JET_PROOF_{n}_{j}',jet['compositionResiduals'][n,j],(j-n,0,0))
    for (degree,n),value in jet['polynomialResiduals'].items():put(f'POLYNOMIAL_PROOF_{degree}_{n}',value,(degree-n,0,0))
    for key,v in result['records'].items():
        for name in ('original','coordinateImage','coordinateReplay','coordinateResidual'):put(key+'_'+name,v[name],v['unit'])
        for name in ('base','mutated','residual'):put(key+'_SHAPE_'+name,v['shape'][name],v['unit'])
        if 'materialSourceAmplitude' in v:
            put(key+'_MATERIAL_SOURCE_AMPLITUDE',v['materialSourceAmplitude'],v['unit'])
            for name in ('rawReconstructionResidual','reconstructionResidual','carrierRoundTripResidual'):put(key+'_SOURCE_JET_'+name,v['sourceJets'][name],v['unit'])
    for i,rec in enumerate(result['orderedLimits']):
        for name in ('originalLimits','materialLimits'):put(f'ORDERED_LIMITS_{i}_'+name,rec[name],(-1,0,0))
        for name in ('originalSourceLimit','materialSourceLimit'):put(f'SOURCE_LIMITS_{i}_'+name,rec[name],(1,0,0))
        put(f'SOURCE_JACOBIAN_{i}',rec['sourceJacobian'],zero)
    d=result['densityAdvection'];unit=engine.PHYSICAL_METADATA.dimensions.measure(d['densityEquation'].rhs);emit_density_equation('RHO4_DENSITY_EQUATION',d['densityEquation']);put('RHO4_DENSITY_GRADIENT',d['gradient'],tuple(a-b for a,b in zip(unit,(1,0,0))));put('RHO4_ADVECTION_FACTOR',d['factor'],zero)
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'records':{k:{'address':v['address'],'shapeOccurrences':v['shape']['count']} for k,v in result['records'].items()},'scope':result['scope'],'densityExtraInput':'background_density_map','independentPhysicalGrades':tuple(map(str,engine.PHYSICAL_METADATA.generators))})


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('coordinate source budget; preserve completed records and chart proofs')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    def progress(stage):
        with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'stage':stage,'wallSeconds':time.monotonic()-start})+'\n')
    r,packet,density,pins,operands=load(base);progress('sources_joined');result=construct(base,r,packet,density,progress)
    result.update(sourceFiles=pins,inputPackets=operands,scope='One nonorthogonal affine material/Eulerian chart, native source-coordinate and derivative maps, literal first-jet mutation and computed constant-rho4 advection operand. No second scattering route or kernel N3/N4/N6 closure yet.')
    f.atomic_pickle(base/'coordinate-source.pickle',result);before=f.digest(base/'coordinate-source.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result);keys={tag:'s11cdCoordinateSource'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique coordinate source tag');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('full coordinate emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:emit_result(result);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'coordinate key/payload census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        for path,fields in body:
            v={str(k):x for k,x in fields};d=v['DIMENSION_L_T_M'];f.require(len(d)==3 and all(not x.free_symbols for x in d) and 'MULTIGRADE' in v and 'EPSILON_LAMBDA_SUPPORT' in v,'coordinate units and independent grades');paths+=1
    last='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[last]},list(entries)[:list(entries).index(last)])
    f.require(before==f.digest(base/'coordinate-source.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'coordinate packet/dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'unchanged coordinate sources')
    checks={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'records':len(result['records']),'nativeTerms':len(result['orderedLimits']),'sourceAmplitudes':sum(v['address'][0]=='source' for v in result['records'].values()),'shapeChangedRecords':{kind:sum(v['address'][0]==kind and v['shape']['residual']!=0 for v in result['records'].values()) for kind in ('local','cell','factor','source')},'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'coordinate-source.pickle'),'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
