#!/usr/bin/env python3
"""Finish coordinate-source validation/emission from completed source packets."""
import argparse,ast,copy,hashlib,json,resource,shutil,signal,time
from pathlib import Path
import sympy as sp
import S11c_d_coordinate_source as c
f=c.f;engine=c.engine;grades=c.grades
SOURCE=f.STORE/'s11c-coordinate-control-20260919/production/complete'
PLAN=f.M/'S11c_d_coordinate_source_recovery_plan.md'
PACKET_SHA='b822833b60cff10d6901a16a29005a5185362f9559487c947184341949b0611c'


def source_join():
    name=str(Path(c.__file__).resolve().relative_to(f.ROOT));old=ast.parse((SOURCE/'source'/name).read_text());new=ast.parse(Path(c.__file__).read_text())
    helpers=[n for n in new.body if isinstance(n,ast.FunctionDef) and n.name=='emit_density_equation'];f.require(len(helpers)==1,'one equation metadata adapter');new.body.remove(helpers[0])
    count=[]
    class Restore(ast.NodeTransformer):
        def visit_Call(self,node):
            self.generic_visit(node)
            if isinstance(node.func,ast.Name) and node.func.id=='emit_density_equation':
                f.require(ast.unparse(node)=="emit_density_equation('RHO4_DENSITY_EQUATION', d['densityEquation'])",'only actual density equation call')
                count.append(1);return ast.parse("put('RHO4_DENSITY_EQUATION',d['densityEquation'],unit)",mode='eval').body
            return node
    Restore().visit(new);f.require(len(count)==1 and ast.dump(old)==ast.dump(new),'whole checker AST: equation emission only')
    return {'wholeCheckerReverseAstJoin':True,'originalSha256':f.digest(SOURCE/'source'/name),'currentSha256':f.digest(Path(c.__file__)),'adapterAstSha256':hashlib.sha256(ast.dump(helpers[0]).encode()).hexdigest()}


def load(base):
    join=source_join();inputs=json.loads((SOURCE/'inputs.json').read_text());outcome=json.loads((SOURCE.parent/'coordinate_construct.invocation.json').read_text())
    f.require(outcome==json.loads((SOURCE.parent/'active.json').read_text()) and outcome['exitCode']==1 and outcome['status']=='failed','original output failure outcome')
    error=(SOURCE.parent/'coordinate_construct.stderr').read_text();f.require("dimension node" in error and "sympy.core.relational.Equality" in error,'recorded equality dimension failure')
    f.require(f.digest(SOURCE/'coordinate-source.pickle')==PACKET_SHA,'complete source packet identity');result=f.unpickle(SOURCE/'coordinate-source.pickle')
    f.require(result['sourceFiles']==inputs['sourceFiles'] and result['inputPackets']==inputs['inputPackets'],'original construction manifest')
    name=str(Path(c.__file__).resolve().relative_to(f.ROOT))
    for n,h in inputs['sourceFiles'].items():
        f.require(f.digest(SOURCE/'source'/n)==h,('frozen source',n))
        if n!=name:f.require(f.digest(f.ROOT/n)==h,('unchanged current source',n))
    for n,h in inputs['inputPackets'].items():f.require(f.digest(Path(n))==h,('unchanged input',n))
    inventory=json.loads((SOURCE/'record-inventory.json').read_text());f.require(len(inventory)==len(result['records'])==375,'complete source records')
    for key,v in inventory.items():f.require(f.digest(SOURCE/v['path'])==v['sha256'] and (SOURCE/v['path']).stat().st_size==v['bytes'],('native record hash',key))
    copied={}
    for path in [*SOURCE.glob('*.pickle'),*(SOURCE/'records').glob('*.pickle')]:
        relative=path.relative_to(SOURCE);dest=base/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest);f.require(f.digest(path)==f.digest(dest),'byte-identical packet copy');copied[str(relative)]={'sha256':f.digest(path),'bytes':path.stat().st_size}
    shutil.copyfile(SOURCE/'record-inventory.json',base/'record-inventory.json')
    pins=dict(inputs['sourceFiles']);pins[name]=f.digest(Path(c.__file__))
    for p in (Path(__file__),PLAN):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    operands=dict(inputs['inputPackets']);operands.update({str(SOURCE/n):v['sha256'] for n,v in copied.items()})
    for p in (SOURCE/'inputs.json',SOURCE/'record-inventory.json',SOURCE/'full.out',SOURCE.parent/'coordinate_construct.invocation.json',SOURCE.parent/'coordinate_construct.stderr'):operands[str(p)]=f.digest(p)
    for n in pins:
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
    rp=next(Path(n) for n in inputs['inputPackets'] if n.endswith('reduced-action.pickle'));r,dimensions=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(rp));dimensions.__dict__.update(copy.deepcopy(result['dimensionState']))
    manifest={'sourceFiles':pins,'inputPackets':operands,'originalInputs':inputs,'originalOutcome':outcome,'emissionRepairJoin':join,'copiedArtifacts':copied,'densityEquationMetadata':'The original equation carrier is unchanged; physical metadata paths0/1 are its exact lhs/rhs operands.'};f.save(base/'inputs.json',manifest)
    return result,pins,operands,manifest


def validate(base,result):
    grade_path=next(Path(n) for n in result['inputPackets'] if n.endswith('continuum-grades.pickle'));native=f.unpickle(grade_path);g=result['chart'];jet=result['fieldJets'];zero_count=0;raw_count=0;counts={}
    f.require(set(result['records'])==set(native['records']),'full native coefficient address census')
    f.require(f.unpickle(base/'chart.pickle')==g and f.unpickle(base/'field-jets.pickle')==jet and f.unpickle(base/'ordered-limits.pickle')==result['orderedLimits'] and f.unpickle(base/'density-advection.pickle')==result['densityAdvection'],'full chart/jet/limit/absence packet joins')
    for name,value in g['proof'].items():
        leaves=list(value) if isinstance(value,sp.MatrixBase) else [value];f.require(all(v==0 for v in leaves),'saved chart proof');zero_count+=len(leaves)
    for group in ('compositionResiduals','polynomialResiduals'):
        for value in jet[group].values():f.require(all(v==0 for v in value),'saved complete jet proof');zero_count+=len(value)
    for key,rec in result['records'].items():
        original=native['records'][key];saved=f.unpickle(base/'records'/(key+'.pickle'))
        f.require(saved==rec and rec['original']==original['record']['ORIGINAL'] and rec['unit']==original['record']['UNIT'] and rec['address']==original['address'],'full original/record/coefficient/unit join')
        kind=rec['address'][0];counts[kind]=counts.get(kind,0)+1
        f.require(rec['coordinateResidual']==0 and sp.expand(rec['coordinateReplay']-rec['encoded'])==0,'saved coordinate certificate')
        f.require(c.coordinate_change(rec['encoded'],g['definitions'],g['forward'])==rec['coordinateImage'],'actual forward coordinate image replay')
        if rec['sourceJets'] is not None:
            v=rec['sourceJets'];small=v['carrierPolynomial'];defs=v['carrierDefinitions'];symbols=v['symbols'];poly=sp.Poly(small,*symbols)
            f.require(small.xreplace(defs)==v['encoded'] and v['encoded'].xreplace(dict(zip(symbols,v['fields'])))==v['original'],'source field/carrier reconstruction')
            f.require(all(poly.coeff_monomial(s).xreplace(defs)==a for s,a in zip(symbols,v['coefficients'])),'source coefficient extraction replay')
            f.require(v['reconstructionResidual']==sp.expand(small-sum(poly.coeff_monomial(s)*s for s in symbols))==0 and v['carrierRoundTripResidual']==0,'exact source polynomial proofs')
            raw_count+=v['rawReconstructionResidual']!=0
        mut=c.first_jet_mutation(rec['original'],type('Profiles',(),{'profiles':{'w':sp.Function('s11cdWProfile')}})())
        f.require(mut==rec['shape'],'actual first-derivative occurrence/baseline/mutation/residual replay')
    f.require(counts=={'local':100,'cell':160,'factor':80,'source':35},'complete native operator/source census')
    f.require(len(result['orderedLimits'])==len(native['termJoins'])==160 and {v['integralIndex'] for v in result['orderedLimits']}==set(range(80)),'all native integral rows')
    for rec,term in zip(result['orderedLimits'],native['termJoins']):
        f.require(all(rec[k]==term[k] for k in ('row','column','term','integralIndex')) and rec['originalLimits']==term['remainingLimits'] and rec['originalSourceLimit']==term['sourceLimit'],'actual term/source/ordered-limit join')
        for old,new in zip((*rec['originalLimits'],rec['originalSourceLimit']),(*rec['materialLimits'],rec['materialSourceLimit'])):
            variable=old[0];m,value=g['definitions'][variable];inverse=sp.solve(sp.Eq(value,variable),m)[0];expected=(m,*(sp.expand(inverse.subs(variable,v.xreplace(g['forward']))) for v in old[1:]));f.require(new==expected,'computed affine native limit maps')
    d=result['densityAdvection'];f.require(d['gradient']==sp.Matrix([sp.diff(d['densityEquation'].rhs,v) for v in g['coordinates']]) and d['factor']==(g['displacement'].T*d['gradient'])[0]/d['densityEquation'].rhs,'actual source density/advection derivation')
    proofs={'recordCensus':counts,'nativeTerms':160,'nativeIntegralRows':80,'exactChartJetProofScalars':zero_count,'rawSourceReconstructionFormsNonzero':int(raw_count),'densityEquation':str(d['densityEquation']),'densityGradient':list(map(str,d['gradient'])),'densityAdvectionFactor':str(d['factor']),'shapeChangedRecords':{kind:sum(v['address'][0]==kind and v['shape']['residual']!=0 for v in result['records'].values()) for kind in counts}}
    f.save(base/'operand-validation.json',proofs);return proofs


def focused_equation(result):
    equation=result['densityAdvection']['densityEquation'];dimensions=engine.PHYSICAL_METADATA.dimensions;seen={};old=engine.emit
    def capture(name,value):seen[name]=engine.cas(value)
    engine.emit=capture
    try:c.emit_density_equation('FOCUSED_EQUATION',equation)
    finally:engine.emit=old
    f.require(seen[c.PREFIX+'_FOCUSED_EQUATION']==engine.carrier_fingerprint(equation),'unchanged original equation fingerprint')
    unit=dimensions.measure(equation.rhs);metadata=seen['METADATA_'+c.PREFIX+'_FOCUSED_EQUATION'];f.require(tuple(tuple(p) for p,_ in metadata)==((0,),(1,)),'physical equation operand paths')
    for _,v in metadata:f.require(tuple(dict((str(k),x) for k,x in v)['DIMENSION_L_T_M'])==unit,'both source equation units')
    dimensions.known[equation.lhs]=(unit[0]+1,*unit[1:]);rejected=False
    try:c.emit_density_equation('FOCUSED_BAD_UNIT',equation)
    except ValueError:rejected=True
    finally:dimensions.known[equation.lhs]=unit
    f.require(rejected,'conflicting known equality unit rejected')
    return {'carrierIdentity':True,'sidePaths':[[0],[1]],'unit':list(map(str,unit)),'changedUnitRejected':rejected}


def emission_tail():
    main=next(n for n in ast.parse(Path(c.__file__).read_text()).body if getattr(n,'name',None)=='main');start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and ast.unparse(n.value)=='engine.EMISSION_LINES.clear()');body=copy.deepcopy(main.body[start:])
    f.require(ast.unparse(body[-1])=='print(json.dumps(checks, indent=2))','unaltered final summary tail');body[-1]=ast.Return(ast.Name('checks',ast.Load()))
    function=ast.FunctionDef(name='tail',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','before','pins','operands','start')],kwonlyargs=[],kw_defaults=[],defaults=[]),body=body,decorator_list=[]);namespace=dict(vars(c));exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),str(Path(__file__)),'exec'),namespace);return namespace['tail']


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('saved coordinate output recovery budget')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    result,pins,operands,manifest=load(base);proofs=validate(base,result);focused=focused_equation(result);f.save(base/'equation-metadata-controls.json',focused)
    before=f.digest(base/'coordinate-source.pickle');checks=emission_tail()(base,result,before,pins,operands,start)
    original={line.partition(': ')[0]:grades._restore(line.rstrip('\n').partition(': ')[2]) for line in grades.decoded_lines(SOURCE/'full.out')};current={line.partition(': ')[0]:grades._restore(line.rstrip('\n').partition(': ')[2]) for line in grades.decoded_lines(base/'full.out')}
    differences=[k for k,v in original.items() if current.get(k)!=v];f.save(base/'emission-differences.json',differences);f.require(not differences and set(original)<=set(current),'every original decoded prefix payload identical')
    f.require(f.digest(SOURCE/'coordinate-source.pickle')==before==PACKET_SHA,'original/copied source packet unchanged')
    for n,v in manifest['copiedArtifacts'].items():f.require(f.digest(SOURCE/n)==f.digest(base/n)==v['sha256'],'original/copied post hashes')
    f.save(base/'recovery.json',{'sourceDirectory':str(SOURCE),'emissionRepairJoin':manifest['emissionRepairJoin'],'originalPrefixTags':len(original),'copiedArtifacts':manifest['copiedArtifacts'],'operandValidation':proofs,'equationMetadataControls':focused,'packetSha256':before,'scope':'Saved operand validation and output only; no coordinate construction or physical recalculation.'})
    f.require(json.loads((base/'checks.json').read_text())==checks,'final persisted checks/stdout identity');signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
