#!/usr/bin/env python3
"""Actual case boundary maps, with exact reuse of accepted full source families."""
import argparse,ast,copy,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_modes as modes
import S11c_d_continuum_boundary as boundary

f,engine=modes.f,modes.engine
CP=f.M/'S11c_d_remaining_case_modes_checkpoint.json'
BCP=f.M/'S11c_d_continuum_boundary_checkpoint.json'
CCP=f.M/'S11c_d_matching_channels_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_boundary_plan.md'
TABLES={}


def array_units(finite,continuum,fields,rows,current_unit):
    """Units of every map entry, derived from its actual domain and range."""
    zero=(0,0,0);minus_length=(-1,0,0);maps={}
    amplitude=[c['amplitudeUnit'] for c in continuum['clusters'] for _ in range(c['R'][(0,0)].shape[1])]
    out_count=continuum['outgoing'][(0,0)].shape[1]
    out_units,inc_units=amplitude[:out_count],amplitude[out_count:]
    def put(key,array,range_units,domain_units,extra=zero,grade=None):
        value=np.asarray(array)
        f.require(value.shape==(len(range_units),len(domain_units)) and np.isfinite(value).all(),'actual map dimension and finite array')
        units=tuple(tuple(tuple(a-b+c for a,b,c in zip(u,v,extra)) for v in domain_units) for u in range_units)
        maps[key]={'shape':value.shape,'units':units,'grade':grade}
    for key,ru,du,extra in [('trace',fields,fields,minus_length),('outgoingInverse',out_units,fields,zero),('outgoing',fields,out_units,zero),('outgoingDerivative',fields,out_units,minus_length),('incoming',fields,inc_units,zero),('incomingDerivative',fields,inc_units,minus_length),('insertion',fields,inc_units,minus_length),('incomingOriginPhase',inc_units,inc_units,zero)]:
        for g,array in continuum[key].items():put(('continuum',key,g),array,ru,du,extra,g)
    open_out=[c['amplitudeUnit'] for c in continuum['clusters'] if c['info']['direction']=='outgoing' and c['info']['kind']=='open' for _ in range(c['R'][(0,0)].shape[1])]
    for g,array in continuum['outgoingOriginPhase'].items():put(('continuum','outgoingOriginPhase',g),array,open_out,open_out,grade=g)
    for name,series in continuum['currents'].items():
        # A current is a sesquilinear form: subtract both amplitude units.
        left=[tuple(a-b for a,b in zip(current_unit,u)) for u in amplitude]
        for g,array in series.items():put(('continuum','current',name,g),array,left,amplitude,grade=g)
    for c in continuum['clusters']:
        amps=[c['amplitudeUnit']]*c['R'][(0,0)].shape[1]
        for key,ru,du,extra in [('R',fields,amps,zero),('D',fields,amps,minus_length),('K',amps,amps,minus_length)]:
            for g,array in c[key].items():put(('continuum','cluster',c['info']['INDEX'],key,g),array,ru,du,extra,g)
        put(('continuum','cluster',c['info']['INDEX'],'baseEquation'),c['diagnostics']['BASE_EQUATION'],rows,amps)
    finite_out=[tuple(v/2 for v in current_unit) if c['kind']=='open' else zero for c in finite['outgoing']]
    finite_inc=[tuple(v/2 for v in current_unit) for _ in finite['incoming']]
    for key,ru,du,extra in [('right',fields,finite_out,zero),('derivative',fields,finite_out,minus_length),('traceMap',fields,fields,minus_length),('residual',fields,finite_out,minus_length),('incomingValues',fields,finite_inc,zero),('incomingDerivative',fields,finite_inc,minus_length),('incomingBoundaryData',fields,finite_inc,minus_length)]:put(('finite',key),finite[key],ru,du,extra)
    put(('finite','openCurrent'),finite['current'],[zero]*len(finite['current']),[zero]*len(finite['current']))
    return maps


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_FOUR_CASE_MODE_SUBSPACES' and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted actual case modes')
    for n,v in cp['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(origin/'source'/n)==v,'accepted current/frozen source')
    for n,v in cp['inputPackets'].items():f.require(f.digest(Path(n))==v,'accepted mode input')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':{str(origin/'checks.json'):cp['checksSha256']},'copiedInputs':{},'input':json.loads((origin/'inputs.json').read_text())['input'],'scope':'Full case finite and retained-continuum end maps before response solves.'}
    for n,v in cp['artifacts'].items():
        modes.retain(origin/n,base/n,manifest,v['sha256'])
    bcp=json.loads(BCP.read_text());bc=Path(bcp['runDirectory'])
    f.require(bcp['status']=='PUBLISHED_ANNEX_VERIFIED' and f.digest(bc/'checks.json')==bcp['checksSha256'],'accepted baseline continuum boundary')
    for n,v in bcp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n)==v,'unchanged accepted boundary helper');manifest['sourceFiles'][n]=v
    modes.retain(bc/'continuum-boundary.pickle',base/'accepted-continuum-boundary.pickle',manifest,bcp['artifacts']['continuum-boundary.pickle']['sha256'])
    ccp=json.loads(CCP.read_text());cc=Path(ccp['runDirectory'])
    for n,v in ccp['sourceFiles'].items():f.require(f.digest(cc/'source'/n)==v,'original channel snapshot')
    old_engine=cc/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    f.require(modes.node(old_engine,'TwoEndedMatchingChannels')==modes.node(engine.HERE,'TwoEndedMatchingChannels'),'unchanged whole native matching-channel constructor')
    manifest['inputPackets'][str(old_engine)]=f.digest(old_engine)
    for end in ('LEFT','RIGHT'):modes.retain(cc/(end.lower()+'.pickle'),base/('accepted-'+end.lower()+'-channels.pickle'),manifest,ccp['artifacts'][end.lower()+'.pickle']['sha256'])
    for path in (Path(__file__).resolve(),PLAN,CP,BCP,CCP):manifest['sourceFiles'][str(path.relative_to(f.ROOT))]=f.digest(path)
    for n,v in manifest['sourceFiles'].items():
        dst=base/'source'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dst);f.require(f.digest(dst)==v,'boundary frozen source')
    f.save(base/'inputs.json',manifest)
    return f.unpickle(base/'remaining-case-end-sources.pickle'),f.unpickle(base/'remaining-case-currents.pickle'),f.unpickle(base/'accepted-continuum-boundary.pickle'),manifest


def restore_builder(pair,value,modal,inputs):
    builder=engine.ModalCurrentSubspaces(pair,value['pairing'],inputs['binding'])
    f.require(modes.same(modal['SYMBOLIC_OPERANDS']['PENCIL_PLUS'],value['pairing']['CLOSED_PENCIL_LEGS'][0]),'actual restored pencil')
    for name,key in [('CURRENT_SLAB','SLAB_CURRENT_MATRIX'),('CURRENT_BULK','BULK_NORMAL_CURRENT_DENSITY_MATRIX')]:
        original=value['pairing'][key]
        source=original.applyfunc(lambda v:dict(engine.polynomial_terms(v,(builder.epsilon,))).get((2,),sp.S.Zero))
        f.require(modes.same(source,modal['SYMBOLIC_OPERANDS'][name]),'actual source current preparation')
    builder.symbolic_operands=modal['SYMBOLIC_OPERANDS'];builder.scalar_operands=modal['SCALAR_OPERANDS'];builder.coefficient_residuals=modal['COEFFICIENT_RESIDUALS']
    f.require(all(v==0 for v in modes.currents.ends_source.scalars(builder.coefficient_residuals)),'accepted full current coefficient extraction')
    expressions={**builder.symbolic_operands,**builder.scalar_operands}
    bound={k:v.xreplace(inputs['binding']) for k,v in expressions.items()}
    f.require(not set().union(*(v.free_symbols for v in bound.values()))-set(builder.variables),'actual complete modal scalar binding')
    builder.evaluate={k:sp.lambdify(builder.variables,v,'numpy',cse=True) for k,v in bound.items()}
    def saved_prepare():
        f.require(modes.same(pair.construct(pair.anchoring,pair.end),value['pairing']),'exact completed preparation source')
        return builder.coefficient_residuals
    builder.prepare=saved_prepare
    return builder


def save_entry(base,end,name,index,value):
    folder=base/'current-tables';folder.mkdir(exist_ok=True);path=folder/(end.lower()+'-'+name+'-'+''.join(map(str,index))+'.pickle')
    f.atomic_pickle(path,{'index':index,'value':value});TABLES[str(path)]={'sha256':f.digest(path),'bytes':path.stat().st_size}
    f.save(base/'current-table-inventory.json',TABLES)


def constructor():
    tree=ast.parse(Path(boundary.__file__).read_text());original=next(v for v in tree.body if getattr(v,'name',None)=='construct_end');table=next(v for v in tree.body if getattr(v,'name',None)=='current_tables')
    new=copy.deepcopy(original);hits=[]
    class Hook(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n)
            if isinstance(n.func,ast.Name) and n.func.id=='current_tables':
                hits.append(True);n.func.id='tracked_tables';n.args=[ast.Name(v,ast.Load()) for v in ('base','end','name')]+n.args
            return n
    Hook().visit(new);f.require(len(hits)==1,'one current-table checkpoint dispatch')
    restored=copy.deepcopy(new)
    class Undo(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n)
            if isinstance(n.func,ast.Name) and n.func.id=='tracked_tables':n.func.id='current_tables';n.args=n.args[3:]
            return n
    Undo().visit(restored);f.require(ast.dump(restored)==ast.dump(original),'whole original boundary constructor joined')
    tracked=copy.deepcopy(table);tracked.name='tracked_tables';tracked.args.args=[ast.arg(v) for v in ('base','end','name')]+tracked.args.args;insertions=[]
    class Save(ast.NodeTransformer):
        def visit_Assign(self,n):
            if any(isinstance(v,ast.Subscript) and isinstance(v.value,ast.Name) and v.value.id=='table' for v in n.targets):
                insertions.append(True);index=copy.deepcopy(n.targets[0].slice)
                hook=ast.Expr(ast.Call(ast.Name('save_entry',ast.Load()),[ast.Name(v,ast.Load()) for v in ('base','end','name')]+[index,ast.Subscript(ast.Name('table',ast.Load()),copy.deepcopy(index),ast.Load())],[]))
                return [n,hook]
            return n
    Save().visit(tracked);f.require(len(insertions)==1,'one durable per-entry checkpoint')
    restored=copy.deepcopy(tracked);restored.name='current_tables';restored.args.args=restored.args.args[3:]
    class Unsave(ast.NodeTransformer):
        def visit_Expr(self,n):
            if isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='save_entry':return None
            return n
    Unsave().visit(restored);f.require(ast.dump(restored)==ast.dump(table),'whole actual derivative constructor joined')
    namespace=dict(vars(boundary),save_entry=save_entry)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[tracked,new],type_ignores=[])),__file__,'exec'),namespace)
    return namespace['construct_end'],{'wholeConstructEndJoin':True,'wholeCurrentTablesJoin':True,'checkpointHooks':2}


def verify_channels(packet,modal):
    f.require(len(packet['CANDIDATES'])==len(modal['RECORDS']),'complete finite candidate census')
    for c,m in zip(packet['CANDIDATES'],modal['RECORDS']):
        f.require(all(c['INFO'][k]==m[k] for k in ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY','K','Q')),'actual full candidate identity')
        for k,v in c['BASES'].items():f.require(np.array_equal(v,m['FORMS'][k]),'actual complete mode bases')
    f.require(boundary.norm(packet['RESIDUALS'])<1e-8,'physical channel current identities')
    f.require(np.isfinite(packet['OUTWARD_CURRENT']).all(),'finite full open current')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    sources,values,accepted,manifest=load(base);build,join=constructor();outputs={};inventory={};new=None
    manifest['constructorJoin']=join;f.save(base/'inputs.json',manifest)
    for label in values['cases']:
        outputs[label]={}
        for end in ('LEFT','RIGHT'):
            target=base/'boundary-cases'/label/end.lower();target.mkdir(parents=True,exist_ok=True)
            with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'address':label+'__'+end,'wallSeconds':time.monotonic()-started})+'\n')
            origin=base/'cases'/label/end.lower();modal,known=f.unpickle(origin/'modal.pickle');adjoint,_=f.unpickle(origin/'adjoint.pickle');inputs=f.unpickle(origin/'mode-inputs.pickle')
            reference,_=f.unpickle(base/'cases'/label/'reference/modal.pickle');baseline,_=f.unpickle(base/'accepted-modes'/end.lower()/'modal.pickle')
            baseline_adjoint,_=f.unpickle(base/'accepted-modes'/end.lower()/'adjoint.pickle')
            value=values['cases'][label][end];baseline_value=values['cases'][modes.BASELINE][end]
            old_pair=f.unpickle(base/'accepted-modes'/end.lower()/'pairing.pickle')[0]['result']
            f.require(modes.same(reference,f.unpickle(base/'accepted-modes/reference/modal.pickle')[0]),'actual common full reference basis')
            f.require(tuple(inputs['fieldUnits'])==tuple(accepted['fieldUnits']),'inherited full field units')
            identical=modes.same(modal,baseline) and modes.same(adjoint,baseline_adjoint) and modes.same(value['pairing'],old_pair) and modes.same(value['acoustic'],baseline_value['acoustic'])
            if identical:
                channel=f.unpickle(base/('accepted-'+end.lower()+'-channels.pickle'));continuum=accepted['ends'][end]
                modes.retain(base/('accepted-'+end.lower()+'-channels.pickle'),target/'channels.pickle',manifest)
                f.atomic_pickle(target/'continuum-boundary.pickle',continuum);f.require(modes.same(f.unpickle(target/'continuum-boundary.pickle'),accepted['ends'][end]),'full original continuum map payload identity');disposition='accepted-full-source-map-reuse'
            elif new is not None:
                old_inputs,channel,continuum,owner=new
                for k in ('binding','physical','relation','rootPacket','pairing','acoustic','fieldUnits','frequency','depth'):f.require(modes.same(inputs[k],old_inputs[k]),('shared actual new end-map source',k))
                for kind in ('channels','continuum-boundary'):modes.retain(owner/(kind+'.pickle'),target/(kind+'.pickle'),manifest)
                disposition='same-complete-new-source-family'
            else:
                f.require(label=='LAB_HELD__RHOBR_CONSTANT' and end=='RIGHT','actual unmatched end source')
                pair,inp,_=modes.context(base,label,end,sources,values,manifest)
                prepared=f.unpickle(origin/'prepared-modal.pickle')
                f.require(modes.same(prepared['symbolic'],modal['SYMBOLIC_OPERANDS']) and modes.same(prepared['scalars'],modal['SCALAR_OPERANDS']) and modes.same(prepared['bindings'],inputs['binding']) and modes.same(prepared['residuals'],modal['COEFFICIENT_RESIDUALS']),'entire actual native preparation reuse')
                builder=restore_builder(pair,value,modal,inputs)
                channel=engine.TwoEndedMatchingChannels(builder,modal,adjoint,end).construct()
                f.atomic_pickle(target/'channels.pickle',channel);verify_channels(channel,modal)
                # Mutate an actual consumed current coefficient, keeping its units.
                pair_operand=channel['PAIR_OPERANDS'][0];record=modal['RECORDS'][pair_operand['LEFT_RECORD_INDEX']];other=modal['RECORDS'][pair_operand['RIGHT_RECORD_INDEX']]
                mutation=record['FORMS']['FLUX_RIGHT'].conj().T@pair_operand['CURRENT_SLAB']@other['FORMS']['FLUX_RIGHT']
                f.atomic_pickle(target/'current-mutation.pickle',{'original':pair_operand,'extraSlabCoefficient':pair_operand['CURRENT_SLAB'],'difference':mutation})
                f.require(boundary.norm(mutation)>1e-10,'actual current coefficient mutation response')
                r=pair.r;dims=engine.PHYSICAL_METADATA.dimensions
                cu={tuple(dims.measure(v)[d]+accepted['fieldUnits'][i][d]+accepted['fieldUnits'][j][d] for d in range(3)) for i in range(5) for j in range(5) if (v:=modal['SYMBOLIC_OPERANDS']['CURRENT_SLAB'][i,j])!=0}
                f.require(cu=={accepted['currentUnit']},'actual physical current unit')
                row_units=[]
                for i in range(5):
                    units={tuple(dims.measure(v)[d]+accepted['fieldUnits'][j][d] for d in range(3)) for j in range(5) if (v:=modal['SYMBOLIC_OPERANDS']['PENCIL_PLUS'][i,j])!=0}
                    f.require(len(units)==1,'actual homogeneous equation row');row_units.append(next(iter(units)))
                f.require(row_units==accepted['rowUnits'],'actual new equation/field unit join')
                continuum=build(end,target,r,inp,{'REFERENCE':reference,end:modal},{end:value['pairing']},channel,accepted['currentUnit'])
                f.atomic_pickle(target/'continuum-boundary.pickle',continuum);new=(inputs,channel,continuum,target);disposition='new-actual-finite-and-continuum-end'
            verify_channels(channel,modal)
            finite=f.boundary_map(channel,{'LEFT':-1,'RIGHT':1}[end]);f.atomic_pickle(target/'finite-boundary.pickle',finite)
            f.require(boundary.norm(finite['residual'])<1e-8,'actual finite boundary reconstruction')
            f.require(len(continuum['census'])==len(reference['RECORDS']),'complete independent-grade end census')
            f.require(boundary.norm(continuum['residuals'])<1e-8,'actual complete invariant-pair/current/boundary residuals')
            units=array_units(finite,continuum,accepted['fieldUnits'],accepted['rowUnits'],accepted['currentUnit']);f.atomic_pickle(target/'array-units.pickle',units)
            outputs[label][end]={'finite':finite,'continuum':continuum,'channels':channel,'sourceInputSha256':f.digest(origin/'mode-inputs.pickle'),'disposition':disposition}
            inventory[label+'__'+end]={'disposition':disposition,'finiteCounts':channel['COUNTS'],'continuumCandidates':len(continuum['census']),'continuumClusters':len(continuum['clusters']),'continuumDirections':int(continuum['offsets'][-1]),'unitArrays':len(units),'maximumResidual':boundary.norm(continuum['residuals']),'finiteTraceTruncationDifference':boundary.norm(continuum['finiteTraceTruncationDifference']),'finiteTraceCondition':finite['condition']}
            f.save(base/'boundary-inventory.json',inventory)
        f.atomic_pickle(base/'boundary-cases'/label/'case-boundary.pickle',{'ends':{end:v['continuum'] for end,v in outputs[label].items()},'finite':{end:v['finite'] for end,v in outputs[label].items()},'fieldUnits':accepted['fieldUnits'],'rowUnits':accepted['rowUnits'],'currentUnit':accepted['currentUnit'],'case':label,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':'Actual finite and retained independent-grade boundary/channel maps, not a scattering solve.'})
    f.require(new is not None and len(inventory)==8,'one complete new family and every case/end')
    f.atomic_pickle(base/'remaining-case-boundary.pickle',{'inventory':inventory,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'currentTableInventory':TABLES})
    f.save(base/'inputs.json',manifest)
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'boundary source pre/post identity')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'boundary input pre/post identity')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'boundary copied packet identity')
    artifacts={str(p.relative_to(base)):modes.currents.ends_source.binding.factors.artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    checks={**manifest,'status':'COMPLETED_FOUR_CASE_BOUNDARY_MAPS','inventory':inventory,'newFiniteChannelConstructions':1,'newContinuumEndConstructions':1,'savedCurrentTableEntries':len(TABLES),'artifacts':artifacts,'wallSeconds':time.monotonic()-started}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
