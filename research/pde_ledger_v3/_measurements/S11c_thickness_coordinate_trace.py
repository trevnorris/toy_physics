#!/usr/bin/env python3
"""Trace the saved unequal-frequency residual to the native kinetic coordinate."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'_measurements')]
import sympy as sp
import S11c_b_brane_operator_sympy_audit as b
from S11c_d_modal_current_check import engine,digest
from S11c_d_end_pairing_emit import small_literal
from ledger_fold import _restore


def native_density():
    """Return the actual density assignment from the unmodified native helper."""
    tree=ast.parse(inspect.getsource(b.kinetic_balance_from_energy))
    function=tree.body[0]
    index=next(i for i,node in enumerate(function.body) if isinstance(node,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id=='kinetic_density' for t in node.targets))
    function=copy.deepcopy(function);function.name='trace_native_density'
    function.body=function.body[:index+1]+[ast.Return(ast.Name('kinetic_density',ast.Load()))]
    module=ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    namespace=dict(vars(b));exec(compile(module,str(b.SCRIPT_PATH),'exec'),namespace)
    return namespace[function.name]


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--pairing-checkpoint',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory;base.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    checkpoint=json.loads(args.pairing_checkpoint.read_text());source=Path(checkpoint['runDirectory'])
    if checkpoint['end']!='RIGHT':raise ValueError('this trace requires the saved RIGHT calculation')
    for name,record in checkpoint['artifacts'].items():
        if digest(source/name)!=record['sha256']:raise ValueError(('pairing artifact pin',name))
    for name,sha in checkpoint['sourceFiles'].items():
        if digest(ROOT/name)!=sha:raise ValueError(('pairing source pin',name))
    packet,known=pickle.loads((source/'complete.pickle').read_bytes());pairing=packet['result']
    emissions=pickle.loads((source/'emissions.pickle').read_bytes())
    files=['scripts/S11c_b_brane_operator_sympy_audit.py','scripts/S11c_a_exports.py',
        'scripts/S11c_b_exports.py','scripts/S11c_c1_exports.py','scripts/S11c_c2_exports.py',
        'directives/S11b_SHARED_PHYSICS.md','directives/S11c_b_SHARED_PHYSICS.md',
        'directives/S11c_d_SHARED_PHYSICS.md','_measurements/S11c_b_inertia_action_check.py',
        '_measurements/S11c_thickness_coordinate_trace.py',*checkpoint['sourceFiles']]
    pins={name:digest(ROOT/name) for name in dict.fromkeys(files)}
    for name in pins:
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dest)
    dims=dict(b.SYMBOL_DIMENSIONS);dims.update({k:sp.ImmutableMatrix(v) for k,v in known.items()})
    d_symbols={str(a):a for value in (pairing['KINETIC_ENERGY_MATRIX'],*pairing['CLOSED_PENCIL_LEGS'])
               for a in engine.dag_free_symbols(value)}
    rename={a:d_symbols[str(a)] for a in b.SYMBOL_DIMENSIONS if str(a) in d_symbols}
    reverse={v:k for k,v in rename.items()}
    eps,eta,sigma=b.epsilon,b.eta_bg,b.sigma_W
    lam=sp.Symbol('s11cThicknessCoordinateLambda');dims[lam]=b.DIM_ZERO
    objects=[];checks={};keys=set()
    def output(name,value,units=None):
        body=engine.cas(value);key='s11cThicknessCoordinate'+name
        if key in keys or key in b.INCOMING_LEDGER:raise ValueError('write-key collision')
        keys.add(key);metadata=[]
        for path,leaf in engine.leaves(body):
            if isinstance(leaf,engine.Str):continue
            unit=units[path] if isinstance(units,dict) else units
            if unit is None:unit=b.dimension_of(leaf,dims)
            expression=leaf.xreplace(reverse)
            if not expression.has(eps):
                expression=expression.xreplace({s:eps*s for s in b.WAVE_SYMBOLS})
            expression=expression.xreplace(b.PROFILE_GRADE_SUBS)
            try:
                sp.Poly(expression,eps,eta,sigma)
                parts={'polynomial':expression}
            except sp.PolynomialError:
                numerator,denominator=sp.together(expression).as_numer_denom()
                parts={'exactNumerator':numerator,'exactDenominator':denominator}
            grade_data={}
            for part,operand in parts.items():
                terms=sp.Poly(sp.expand(operand),eps,eta,sigma)
                grades=sorted(g for g,c in terms.terms() if c!=0)
                homotopy=sp.Poly(sp.expand(operand.subs({eta:lam,sigma:lam*b.W0/b.L_W})),eps,lam)
                support=sorted(g for g,c in homotopy.terms() if c!=0)
                grade_data[part]={'multigrade':grades,'epsilonLambdaSupport':support}
            metadata.append({'path':path,'dimensionLTM':tuple(unit),'gradeData':grade_data})
        heavy=not small_literal(body)
        record={'writeKey':key,'value':engine.carrier_fingerprint(body) if heavy else body,
                'representation':'CARRIER_PIT_SHA' if heavy else 'LITERAL','metadata':metadata}
        print(json.dumps({k:sp.srepr(engine.cas(v)) for k,v in record.items()}),flush=True)
        objects.append({'key':key,'body':body,'record':record})
        if name.endswith('Residual'):
            values=[v for _,v in engine.leaves(body)]
            checks[name]={'scalars':len(values),'nonzero':sum(v!=0 for v in values)}
    output('Provenance',{'sourcePins':pins,'pairingCheckpointSha256':digest(args.pairing_checkpoint),
        'pairingPacketSha256':checkpoint['objectsSha256']})
    local=b.local_thickness_map()[0]
    declared=next(eq.rhs for eq in b.profile_definitions() if isinstance(eq,sp.Equality) and eq.lhs==b.e_W_bg)
    physical=sp.expand(b.W_bg*local)
    output('LocalThicknessMap',local)
    output('DeclaredLocalMap',declared)
    output('LocalMapJoinResidual',sp.cancel(local-declared),b.DIM_ZERO)
    output('PhysicalThickness',physical)
    native=native_density();action_rows={};native_rows={};energies={}
    for representative in b.DENSITY_REPS:
        density=b.density_pair(representative)[1]
        old=native(density)
        velocity=sp.diff(physical,b.e_W)*b.e_t
        # Supplied kinetic action, after the existing physical-coordinate map.
        action=b.epsilon**2*(density*b.dot(b.u_t,b.u_t)+b.mu_W*velocity**2)/2
        velocities=(*b.u_t,b.e_t);accelerations=(*b.u_tt,b.e_tt)
        def inertia(T):return sp.ImmutableMatrix([sum(sp.diff(T,v,w)*a for w,a in zip(velocities,accelerations))/b.epsilon for v in velocities])
        before=inertia(old);after=inertia(action)
        assembled=sp.ImmutableMatrix([*b.kinetic_balance_from_energy(density)[0],b.kinetic_balance_from_energy(density)[1]])
        label=''.join(v.title() for v in representative.lower().split('_'))
        row_units={ (i,):b.DIM_BODY_U if i<3 else b.DIM_ENERGY for i in range(4)}
        for name,obj in [('NativeEnergy',old),('MappedActionEnergy',action)]:output(label+name,obj)
        for name,obj in [('NativeInertia',before),('MappedActionInertia',after),
                         ('NativeAssemblyJoinResidual',assembled-before),('CoordinateDifference',before-after)]:
            output(label+name,obj,row_units)
        output(label+'UniformCoordinateJoinResidual',(before-after).subs(b.W_bg,b.W0),row_units)
        native_rows[representative]=before;action_rows[representative]=after
        energies[representative]=(old,action)
    # Extract the actual inertia coefficient before and after the coordinate map.
    representative='RHO4_CONSTANT'
    old_mass=sp.diff(native_rows[representative][3],b.e_tt)/b.epsilon
    action_mass=sp.diff(action_rows[representative][3],b.e_tt)/b.epsilon
    output('NativeThicknessMass',old_mass);output('MappedThicknessMass',action_mass)
    # Bind the saved RIGHT endpoint's source profile, not a generic random witness.
    source_checkpoint_path=ROOT/'_measurements/S11c_d_end_current_source_right_trace_repair_checkpoint.json'
    if digest(source_checkpoint_path)!=checkpoint['provenance']['sourceCheckpointSha256']:
        raise ValueError('endpoint source checkpoint pin')
    source_checkpoint=json.loads(source_checkpoint_path.read_text())
    source_packet=pickle.loads((Path(source_checkpoint['runDirectory'])/'objects.pickle').read_bytes())
    if digest(Path(source_checkpoint['runDirectory'])/'objects.pickle')!=source_checkpoint['artifacts']['objects.pickle']['sha256']:
        raise ValueError('endpoint source packet pin')
    w_end=next(v for k,v in source_packet['profileBindings'].items() if isinstance(k,sp.Limit)
               and k.args[0].func.__name__=='s11cdWProfile' and k.args[2]==sp.oo)
    output('RightProfileLimit',w_end)
    def project(expr):return sp.expand(sum(sp.diff(expr,eta,i,sigma,j).subs({eta:0,sigma:0})*eta**i*sigma**j for i in range(2) for j in range(2)))
    velocities=(*b.u_t,b.theta_t,b.e_t)
    mass_matrices=[sp.hessian(T,velocities).diff(b.mu_W)*b.mu_W/b.epsilon**2
                   for T in energies[representative]]
    masses=[m.applyfunc(lambda v:project(v.xreplace(b.PROFILE_GRADE_SUBS).subs(b.w1_profile,w_end)).xreplace(rename))
            for m in mass_matrices]
    left,right=pairing['FREQUENCY_LEGS'];mu=d_symbols['mu_W']
    # Derive the harmonic inertial row from the acceleration of the field ansatz.
    t=sp.Symbol('s11cThicknessCoordinateTime',real=True)
    planes=[sp.exp(-sp.I*w*t) for w in (right,-left)]
    deltas=[]
    for leg,(pencil,plane) in enumerate(zip(pairing['CLOSED_PENCIL_LEGS'],planes)):
        actual=pencil.diff(mu)*mu
        acceleration=sp.cancel(sp.diff(plane,t,2)/plane)
        row=masses[0]*acceleration
        mapped=masses[1]*acceleration
        metadata=next(p for p in emissions if p['tag'].endswith('_CONSTRUCTION_CLOSED_PENCIL_LEGS'))['metadata']
        units={}
        for item in metadata:
            fields={str(k):v for k,v in item};path=tuple(int(v) for v in fields['OBJECT_PATH'])
            if path[0]==leg:units[path[1:]]=tuple(fields['VALUE_DIMENSION_L_T_M'])
        output('Leg'+str(leg)+'ImportedInertia',actual,units)
        output('Leg'+str(leg)+'NativeSourceInertia',row,units)
        output('Leg'+str(leg)+'NativeSourceJoinResidual',(actual-row).applyfunc(sp.expand),units)
        output('Leg'+str(leg)+'MappedActionInertia',mapped,units)
        deltas.append(actual-mapped)
    source_difference=pairing['PLUS_ROW_POWER_MAP']*deltas[0]+deltas[1].T*pairing['MINUS_ROW_POWER_MAP']
    predicted=-source_difference.xreplace(packet['bindings'])
    # Restore each field-pair unit from the saved balance metadata.
    entry=next(p for p in emissions if p['tag'].endswith('_RETAINED_SLAB_BALANCE_RESIDUAL'))
    units={tuple(int(v) for v in dict((str(k),v) for k,v in item)['OBJECT_PATH']):
           tuple(dict((str(k),v) for k,v in item)['VALUE_DIMENSION_L_T_M']) for item in entry['metadata']}
    output('PredictedRightBalanceTerm',predicted,units)
    for key in ('SLAB_BALANCE_RESIDUAL','FINITE_COMPOSED_BALANCE_RESIDUAL','EQUAL_DEPTH_COMPOSED_BALANCE_RESIDUAL'):
        label=''.join(v.title() for v in key.lower().split('_'))
        output(label+'Saved',packet['retained'][key],units)
        output(label+'AccountingResidual',(packet['retained'][key]-predicted).applyfunc(sp.expand),units)
    output('EqualFrequencyTerm',predicted.subs(left,right).applyfunc(sp.expand),units)
    elapsed=time.monotonic()-started
    output('WallSeconds',sp.Float(elapsed),b.DIM_T)
    output('PeakRssKiB',sp.Integer(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),b.DIM_ZERO)
    payload=base/'objects.pickle';payload.write_bytes(pickle.dumps(objects,protocol=5))
    summary={'sourcePins':pins,'sourcePinsAfter':{n:digest(ROOT/n) for n in pins},'checks':checks,
        'objects':len(objects),'objectsSha256':digest(payload),'wallSeconds':elapsed,
        'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'pairingCheckpointSha256':digest(args.pairing_checkpoint)}
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if any(v['nonzero'] for v in checks.values()):raise ValueError('trace reconstruction residual; inspect emitted objects')


if __name__=='__main__':run()
