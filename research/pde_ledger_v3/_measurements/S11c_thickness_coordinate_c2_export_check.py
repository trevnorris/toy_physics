#!/usr/bin/env python3
"""Track the measured b row change through c2's field and weak restrictions."""
import argparse
import ast
import copy
import faulthandler
import hashlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time

STARTED=time.monotonic()
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'_measurements')]
import sympy as sp
import S11c_c2_selfenergy_fold_sympy_audit as c2
import S11c_d_mixing_scattering_sympy_audit as d
from ledger_fold import _restore,load_model
from S11c_inertia_artifact_audit import export_data,locate
from S11c_d_end_pairing_emit import small_literal
from S11c_thickness_coordinate_origin_components import origin_components


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def leaves(value,units,path=()):
    if isinstance(units,sp.MatrixBase):
        yield path,value,tuple(units)
        return
    association=all(isinstance(item,sp.Tuple) and len(item)==2
                    and isinstance(item[0],c2.Str) for item in value)
    if len(value)!=len(units):raise ValueError('dimension tree arity')
    for index,(item,unit) in enumerate(zip(value,units)):
        if association:
            if item[0]!=unit[0]:raise ValueError('dimension tree key')
            yield from leaves(item[1],unit[1],path+(str(item[0]),))
        else:yield from leaves(item,unit,path+(index,))


def closure_slots(inputs):
    """Read the replacement-key assignments from the actual native closure."""
    source=ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py'
    function=next(n for n in ast.parse(source.read_text()).body
                  if isinstance(n,ast.FunctionDef) and n.name=='build_face')
    statements={n.targets[0].id:n for n in function.body if isinstance(n,ast.Assign)
                and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)}
    names=[key.id for key in statements['replacements'].value.keys]
    statements_to_run=[statements['label'],*(statements[name] for name in names)]
    code=compile(ast.Module(body=statements_to_run,type_ignores=[]),str(source),'exec')
    result=[]
    for face in c2.FACES:
        context={'inputs':inputs,'face':face}
        exec(code,context)
        result.extend(context[name] for name in names)
    return tuple(result)


def kinetic_action_from_source(source):
    """Compile the actual kinetic action assignments, before row normalization."""
    function=next(node for node in ast.parse(source.read_text()).body
        if isinstance(node,ast.FunctionDef) and node.name=='conservative_power_variation')
    names={'fields','rates','rho','kinetic','kinetic_action'}
    assignments=[copy.deepcopy(node) for node in function.body if isinstance(node,ast.Assign)
        and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in names]
    if {node.targets[0].id for node in assignments}!=names:raise ValueError('kinetic source assignment coverage')
    code=compile(ast.fix_missing_locations(ast.Module(body=assignments,type_ignores=[])),str(source),'exec')
    def compute(inputs,case):
        scope={**vars(c2),'inputs':inputs,'case':case,'density':case[1]}
        exec(code,scope)
        return scope['kinetic'],scope['kinetic_action']
    return compute


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--b-checkpoint',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory;base.mkdir(parents=True,exist_ok=False)
    stack_stream=(base/'stacks.log').open('w')
    faulthandler.dump_traceback_later(300,repeat=True,file=stack_stream)
    def progress(operation,**values):
        with (base/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'operation':operation,'wallSeconds':time.monotonic()-STARTED,**values})+'\n')
    progress('source_snapshot')
    paths=[ROOT/'scripts'/('S11c_'+stage+'_exports.py') for stage in ('b','c1','c2')]
    paths += [args.baseline/'scripts'/('S11c_'+stage+'_exports.py') for stage in ('b','c1','c2')]
    paths += [Path(__file__),args.b_checkpoint,ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
        ROOT/'scripts/ledger_fold.py',ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py',
        ROOT/'scripts/S11c_d_output_codec.py',ROOT/'directives/S11c_c2_SHARED_PHYSICS.md',
        *(ROOT/'_measurements'/name for name in ('S11c_inertia_artifact_audit.py','S11c_d_end_pairing_emit.py',
          'S11c_d_end_pairing_check.py','S11c_d_modal_current_check.py','S11c_d_joint_sheet_check.py'))]
    paths.append(ROOT/'_measurements/S11c_thickness_coordinate_origin_components.py')
    old_control_source=args.baseline/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py'
    paths.append(old_control_source)
    pins={str(p):sha(p) for p in paths};snapshots={}
    for i,path in enumerate(paths):
        relative=Path('source')/str(i)/path.name;target=base/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target);snapshots[str(path)]=str(relative)
    proof=json.loads(args.b_checkpoint.read_text())
    if proof['stage']!='b' or proof['nonzeroResidualScalars'] or proof['newExportSha256']!=sha(paths[0]):
        raise ValueError('completed current b action/export checkpoint required')
    old={};new={};export_pins={}
    for i,stage in enumerate(('b','c1','c2')):
        new[stage],export_pins[stage],_=export_data(paths[i])
        old[stage],_,_=export_data(paths[i+3])
    old_origins=origin_components(old['b']['slab_operator_term_origins'])
    new_origins=origin_components(new['b']['slab_operator_term_origins'])
    fold,_=load_model(str(paths[0]),str(paths[1]));inputs=c2.bind_inputs(fold)
    slot_keys=closure_slots(inputs)
    lam=sp.Symbol('s11cThicknessCoordinateC2ExportLambda');unused=sp.Dummy('unusedGrade')
    records=[];checks={};keys=set()
    def output(name,value,unit=(0,0,0)):
        progress('emit',object=name)
        body=d.cas(value);key='s11cThicknessCoordinateC2'+name
        if key in keys or key in fold or key in new['c2']:raise ValueError('write-key collision')
        keys.add(key);metadata=[]
        for path,leaf in d.leaves(body):
            if isinstance(leaf,d.Str):continue
            expression=leaf.xreplace(inputs.profiles)
            support=c2.grades(expression,inputs.eps,inputs.eta,inputs.sigma)
            homotopy=expression.subs({inputs.eta:lam,inputs.sigma:lam*inputs.values['W_0']/inputs.values['L_W']})
            lambda_support=c2.grades(homotopy,inputs.eps,lam,unused)
            metadata.append({'path':path,'dimensionLTM':unit(path) if callable(unit) else unit,
                'gradeConvention':'NATIVE_STRUCTURAL_RETAINED_SUPPORT',
                'multigrade':sorted(support),'epsilonLambdaSupport':sorted({g[:2] for g in lambda_support})})
        heavy=not small_literal(body)
        record={'writeKey':key,'value':d.carrier_fingerprint(body) if heavy else body,
                'representation':'CARRIER_PIT_SHA' if heavy else 'LITERAL','metadata':metadata}
        print(json.dumps({k:sp.srepr(d.cas(v)) for k,v in record.items()}),flush=True)
        records.append({'key':key,'body':body,'record':record})
        if name.endswith('Residual'):
            values=[v for _,v in d.leaves(body) if not isinstance(v,d.Str)]
            checks[name]={'scalars':len(values),'nonzero':sum(v!=0 for v in values)}
    output('SourcePins',pins)
    output('ExportPinResidual',{stage:{name:sp.Integer(sha(locate(name))!=digest) for name,digest in values.items()}
                               for stage,values in export_pins.items()})
    output('BActionCheckpointSourceJoinResidual',{
        'baselineBExport':sp.Integer(proof['sourcePins'][str(paths[3])]!=sha(paths[3]))})
    changed={stage:sorted(k for k in new[stage].keys()&old[stage].keys() if new[stage][k]!=old[stage][k])
             for stage in new}
    output('ChangedValueSerializations',changed)
    output('KeySetDifference',{stage:(sorted(new[stage].keys()-old[stage].keys()),
                                    sorted(old[stage].keys()-new[stage].keys())) for stage in new})
    # KINETIC provenance changes together with slab_operator and was separately
    # joined to the action by the current b proof. Every other provenance
    # component (including face work and stored energy) must stay identical.
    output('OriginComponentKeySetDifference',(sorted(new_origins.keys()-old_origins.keys()),
                                             sorted(old_origins.keys()-new_origins.keys())))
    if old_origins.keys()!=new_origins.keys():raise ValueError('origin key set changed; emitted')
    origin_checks={}
    origin_operands=[]
    for path in sorted(new_origins):
        if path[1:]==('VALUE','KINETIC'):
            label=''.join(part.title().replace('_','') for part in path[0])
            for suffix in ('OriginActionResidual','OriginBaselineSourceResidual','OriginDeltaAccountingResidual'):
                evidence=proof['checks'][label+suffix]
                if evidence['scalars']!=4 or evidence['nonzero']:
                    raise ValueError(('b kinetic origin evidence incomplete',path,suffix))
        else:
            before,after=old_origins[path],new_origins[path]
            origin_checks[str(path)]=sp.Integer(before!=after)
            origin_operands.append((path,hashlib.sha256(before.encode()).hexdigest(),hashlib.sha256(after.encode()).hexdigest()))
    output('OtherOriginIdentityOperands',origin_operands)
    output('OtherOriginPreservationResidual',origin_checks)
    # All other input rows are exact preservation prerequisites for transporting
    # only the measured slab delta through the native field and weak maps.
    dependency_checks={}
    for stage in ('b','c1'):
        dependency_checks[stage]={key:sp.Integer(old[stage][key]!=new[stage][key])
            for key in old[stage].keys()&new[stage].keys()
            if key not in ('slab_operator','slab_operator_term_origins')}
    output('OtherInputPreservationResidual',dependency_checks)
    if any(old[stage].keys()!=new[stage].keys() for stage in new):raise ValueError('export key set changed; emitted')
    if any(v for values in dependency_checks.values() for v in values):raise ValueError('closure dependency changed; emitted')
    if any(origin_checks.values()):raise ValueError('nonkinetic origin changed; emitted')
    historical_action=kinetic_action_from_source(old_control_source)
    native_action=kinetic_action_from_source(ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py')
    for case in inputs.slab:
        label=''.join(part.title().replace('_','') for part in case)+'KineticControl'
        current=c2.conservative_power_variation(inputs,case)
        old_energy,old_action=historical_action(inputs,case)
        new_energy,new_action=native_action(inputs,case)
        normalized=c2.tree(current['ACTION_TO_ROW_MULTIPLIER']*new_action,lambda v:c2.retained_shape(v,inputs))
        historical=c2.tree(current['ACTION_TO_ROW_MULTIPLIER']*old_action,lambda v:c2.retained_shape(v,inputs))
        output(label+'BaselineActionEnergy',old_energy,(-1,-2,1))
        output(label+'NativeActionEnergy',new_energy,(-1,-2,1))
        output(label+'ActionEnergyDelta',c2.difference(new_energy,old_energy),(-1,-2,1))
        units=lambda path:(-2,-2,1) if path[0]<3 else (-1,-2,1)
        for name,value in (
            ('BaselineActionRows',historical),('NativeActionRows',normalized),
            ('ImportedOriginRows',current['KINETIC_SOURCE_ROWS']),
            ('NormalizationResidual',current['KINETIC_NORMALIZATION_RESIDUAL']),
            ('SourceAssemblyResidual',c2.difference(new_action,current['KINETIC_ACTION_ROWS'])),
            ('BaselineToCorrectedDifference',c2.difference(historical,normalized))):
            output(label+name,value,units)
    old_b=c2.cases(_restore(old['b']['slab_operator']));new_b=c2.cases(_restore(new['b']['slab_operator']))
    expanded={};component_census=[]
    for case in new_b:
        source_delta=c2.difference(c2.expanded_rows(new_b[case]),c2.expanded_rows(old_b[case]))
        retained=c2.tree(source_delta,lambda v:c2.retained_shape(v,inputs))
        physical=c2.tree(retained,inputs.physical_fields)
        kernel=c2.tree(c2.extract(physical,inputs),lambda v:c2.retained_shape(v,inputs))
        expanded[case]={'s11cc2ClosedSlabOperator':d.cas(physical),'s11cc2ClosedCouplingKernel':d.cas(kernel)}
        label=''.join(part.title().replace('_','') for part in case)
        output(label+'ClosureReplacementKeys',tuple(str(k) for k in slot_keys))
        output(label+'ClosureSlotDependenceResidual',tuple(sp.Integer(v.has(k))
            for _,v in d.leaves(d.cas(source_delta)) for k in slot_keys))
    for root in sorted(c2.EXPORT_ROOTS):
        before_cases=dict(_restore(old['c2'][root]));after_cases=dict(_restore(new['c2'][root]))
        if before_cases.keys()!=after_cases.keys():raise ValueError('c2 case keys differ')
        for axes,payload in after_cases.items():
            case=tuple(map(str,axes));label=root.removeprefix('s11cc2')+''.join(v.title().replace('_','') for v in case)
            prior=before_cases[axes]
            units=c2.named(payload,'DIMENSION_L_T_M')
            prior_units=c2.named(prior,'DIMENSION_L_T_M')
            if prior_units!=units:raise ValueError('c2 output units changed')
            before_leaves={path:v for path,v,_ in leaves(c2.named(prior,'VALUE'),prior_units)}
            expected_leaves={path:v for path,v,_ in leaves(expanded[case][root],units)}
            for index,(path,after,unit) in enumerate(leaves(c2.named(payload,'VALUE'),units)):
                before=before_leaves[path];expected=expected_leaves[path]
                delta=c2.difference(after,before)
                atoms=before.atoms(sp.Derivative)|after.atoms(sp.Derivative)
                accelerations={v:sp.S.Zero for v in atoms if sum(n for x,n in v.variable_count if x==c2.TIME)>=2}
                before_static,after_static=before.xreplace(accelerations),after.xreplace(accelerations)
                for name,value in [('BeforeOperand',before),('AfterOperand',after),('RowDelta',delta),
                    ('ImportedDeltaAfterNativeMaps',expected),('FullDeltaAccountingResidual',c2.difference(delta,expected)),
                    ('ZeroAccelerationBeforeOperand',before_static),('ZeroAccelerationAfterOperand',after_static),
                    ('NonkineticPreservationResidual',c2.difference(after_static,before_static))]:
                    output(label+'Component'+str(index)+name,value,unit)
                component_census.append({'root':root,'case':list(case),'path':list(path),'deltaIsZero':delta==0})
            output(label+'OtherPayloadSlotPreservationResidual',{
                str(key):sp.Integer(value!=c2.named(prior,str(key))) for key,value in payload if str(key)!='VALUE'})
    output('OtherC2ExportPreservationResidual',{k:sp.Integer(new['c2'][k]!=old['c2'][k])
        for k in new['c2'] if k not in c2.EXPORT_ROOTS})
    output('SourcePinStabilityResidual',{name:sp.Integer(sha(Path(name))!=digest) for name,digest in pins.items()})
    output('ComponentCensus',component_census)
    elapsed=time.monotonic()-STARTED
    output('WallSeconds',sp.Float(elapsed),(0,1,0))
    output('PeakRssKiB',sp.Integer(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    target=base/'objects.pickle';target.write_bytes(pickle.dumps(records,protocol=5))
    summary={'stage':'c2','newExportSha256':sha(paths[2]),'sourcePins':pins,
        'sourcePinsAfter':{p:sha(Path(p)) for p in pins},'sourceSnapshots':snapshots,'checks':checks,
        'objects':len(records),'objectsSha256':sha(target),'cases':[list(case) for case in new_b],
        'components':component_census,'changedValueSerializations':changed['c2'],
        'nonzeroResidualScalars':sum(v['nonzero'] for v in checks.values()),
        'wallSeconds':elapsed,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    progress('complete',nonzeroResidualScalars=summary['nonzeroResidualScalars'])
    faulthandler.cancel_dump_traceback_later();stack_stream.close()
    if summary['nonzeroResidualScalars']:raise ValueError('c2 export residual; inspect emitted operands')


if __name__=='__main__':run()
