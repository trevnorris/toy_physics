#!/usr/bin/env python3
"""Only missing scalar frequency bindings; retain exact saved whole end-map routes."""
import argparse
import ast
import copy
import inspect
import json
from pathlib import Path
import resource
import signal
import time
from types import SimpleNamespace

import sympy as sp
import S11c_d_frequency_matrix as native
import S11c_d_remaining_case_frequency_end_pilot as storage

saved = storage.saved
M, REPO = saved.M, saved.REPO
require, same = saved.require, saved.same
CP = M / 'S11c_d_remaining_case_frequency_analytic_checkpoint.json'
CP_SHA = 'cdf2472a5882b8e999bfce73c6ffefcd1fdf996b1faf0761d9b49e388dac61bb'
PILOT = M / 'S11c_d_remaining_case_frequency_end_pilot_checkpoint.json'
PILOT_SHA = 'df9ebba15ab7c4e02bd438a70675fac1bfb13924914384dea727e560b541cf2d'
MATRIX = M / 'S11c_d_frequency_matrix_checkpoint.json'
PLAN = M / 'S11c_d_remaining_case_frequency_scalar_bindings_plan.md'
BASELINE = 'LAB_HELD__RHO4_CONSTANT'


def scalar_loop():
    original = ast.parse(inspect.getsource(native.bind)).body[0]
    changed = copy.deepcopy(original)
    count = 0
    class Route(ast.NodeTransformer):
        def visit_Call(self, node):
            nonlocal count
            if ast.unparse(node) == "record['analytic'].subs({w: frequency, **origin}, simultaneous=True)":
                count += 1
                return ast.copy_location(ast.parse('router.bind(key,record,old,w,frequency,origin)', mode='eval').body, node)
            return self.generic_visit(node)
    changed = Route().visit(changed)
    require(count == 1, 'literal native scalar frequency call')
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            if ast.unparse(node.func) == 'router.bind':
                return ast.copy_location(ast.parse("record['analytic'].subs({w:frequency,**origin},simultaneous=True)", mode='eval').body,node)
            return self.generic_visit(node)
    require(ast.dump(Reverse().visit(copy.deepcopy(changed))) == ast.dump(original), 'whole native bind reverse AST')
    init = next(n for n in changed.body if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'actual')
    first = changed.body.index(init)
    loop = next(n for n in changed.body if isinstance(n, ast.For) and ast.unparse(n.iter) == "chart['records'].items()")
    # Restore captured caller objects, then execute the unchanged scalar loop
    # only. The completed baseline jet/differentiation/matrix suffix is not run.
    body = ast.parse("origin=variables['origin']; w=variables['frequency']; chart={'records':records}; oldsrc={'records':sources}").body
    body += copy.deepcopy(changed.body[first:changed.body.index(loop)+1])
    body += ast.parse("return {'actual':actual,'mapping':mapping,'joins':joins}").body
    fn = ast.FunctionDef(name='run_scalar_loop', args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=k) for k in
             ('records','sources','variables','frequency','router')], kwonlyargs=[],kw_defaults=[],defaults=[]), body=body,decorator_list=[])
    module = ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    env = {'f':SimpleNamespace(require=require)}
    exec(compile(module,'<native scalar binding loop only>','exec'),env)
    return env['run_scalar_loop'], {'wholeNativeBindReverseAST':True,'nativeBody':ast.unparse(original),
            'executedScalarLoop':ast.unparse(module),'unexecutedSuffix':'source jets, basis derivatives, source amplitudes, local matrices, row/cell assembly'}


def call_input(record, w, frequency, origin):
    return {'expression':record['analytic'], 'mapping':{w:frequency,**origin},'kwargs':{'simultaneous':True},
            'liveSource':record['originalLive'],'unit':record['unit'],'limits':record['limits']}


class Router:
    def __init__(self, journal, atlas):
        self.journal, self.atlas = journal, atlas
        self.new = 0
        self.routes = []

    def bind(self, key, record, old, w, frequency, origin):
        requested = call_input(record,w,frequency,origin)
        require(same(record['address'],old['address']) and same(record['originalLive'],old['liveFrequencyAndGrades'])
                and same(record['unit'],old['unit']), 'own live source/address/unit join')
        bucket = self.atlas.setdefault(hash(requested['expression']), [])
        matches = [v for v in bucket if same(v['input'],requested)]
        consumer = {'case':self.case,'key':key,'address':record['address'],'analyticOwner':record['analyticOwner'],
                    'sourceOwner':record['sourceOwner'],'context':self.context,'analyticPacket':self.analytic_route,'sourcePacket':self.source_route}
        if matches:
            chosen = matches[0]
            require(all(same(chosen['value'],v['value']) for v in matches), 'all saved complete scalar calls agree')
            value, owner = chosen['value'],chosen['owner']
            route = {'consumer':consumer,'owner':owner,'kind':'REUSED_COMPLETE_CALL','allMatchingOwners':[v['owner'] for v in matches]}
        else:
            folder = 'operations/' + str(self.new)
            self.journal.write(folder+'/input.pickle',{'call':requested,'consumer':consumer})
            start = time.monotonic()
            # Exactly the accepted native call and keyword; no normalization.
            value = record['analytic'].subs({w:frequency,**origin},simultaneous=True)
            result = self.journal.write(folder+'/value.pickle',value)
            self.journal.json(folder+'/completed.json',{'value':result,'wallSeconds':time.monotonic()-start,
                               'native':'S11c_d_frequency_matrix.bind scalar assignment'})
            owner = {'kind':'new-scalar-frequency-binding','case':self.case,'key':key,'input':self.journal.artifacts[folder+'/input.pickle'],'value':result}
            bucket.append({'input':requested,'value':value,'owner':owner})
            self.new += 1
            route = {'consumer':consumer,'owner':owner,'kind':'NEW_COMPLETE_CALL'}
        self.routes.append(route)
        return value


def prohibit():
    def forbidden(*args,**kwargs):
        raise RuntimeError('completed scientific operation is outside missing scalar binding stage')
    for name in ('diff','lambdify','cancel','factor','simplify','solve','resultant','gcd','gcdex','integrate'):
        setattr(sp,name,forbidden)
    for name in ('load','bind','seed_prefix','end_maps','assemble','system_and_solve','seed_case','complex_case','main'):
        setattr(native,name,forbidden)
    for name in ('load','seeds','maps','continue_pair','focused','main'):
        setattr(native.end,name,forbidden)
    native.end.Pair.__init__ = forbidden


def main():
    parser = argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    base = parser.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    reader,journal=saved.Reader(),storage.Journal(base)
    cp=reader.json(CP,CP_SHA); origin=Path(cp['runDirectory'])
    require(cp['status']=='ACCEPTED_CASE_FREQUENCY_ANALYTIC_SOURCES','accepted complete analytic sources')
    reader.retain(origin/'checks.json',cp['checksSha256'])
    pilot=reader.json(PILOT,PILOT_SHA);require(pilot['status']=='ACCEPTED_BOUNDED_CASE_FREQUENCY_END_CONTINUATION','accepted own bounded end map')
    reader.retain(Path(pilot['runDirectory'])/'checks.json',pilot['checksSha256'])
    matrix=reader.json(MATRIX,'0c350c5f9c554c2e585920745157be826daf51775e05676a0554d8670a750663');mr=Path(matrix['runDirectory']);reader.retain(mr/'checks.json',matrix['checksSha256'])
    require(matrix['status']=='VALIDATED_FIRST_FINITE_FREQUENCY_MATRIX','accepted baseline complete frequency binding/map')
    for p in (Path(__file__).resolve(),PLAN,Path(saved.__file__),Path(storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(p)
    native_source=M/'S11c_d_frequency_matrix.py'
    h=matrix['sourceFiles']['_measurements/S11c_d_frequency_matrix.py'];reader.retain(native_source,h);reader.retain(mr/'source/_measurements/S11c_d_frequency_matrix.py',h)
    run,join=scalar_loop();journal.json('native-scalar-binding-join.json',join)
    prohibit()
    def packet(name):return reader.packet(origin/name,cp['artifacts'][name]['sha256'])
    def route(name):return reader.retain(origin/name,cp['artifacts'][name]['sha256'])
    baseline=packet('analytic-cases/'+BASELINE+'/frequency-analytic.pickle')
    baseline_source=packet('accepted/frequency-cases/'+BASELINE+'/frequency-source.pickle')
    bname='complex/frequency-binding.pickle';bpath=mr/bname
    bound=reader.packet(bpath,matrix['artifacts'][bname]['sha256']);bound_route=reader.retain(bpath)
    frequency=bound['frequency']
    require(frequency==sp.S.One-sp.I/100,'actual accepted one-frequency caller')
    baseline_by_key={row['key']:row for row in bound['sourceJoins']}
    require(set(baseline_by_key)==set(baseline['records']),'all saved baseline bound records')
    atlas={};atlas_routes=[]
    variables=baseline['variables']
    for key, record in baseline['records'].items():
        actual=baseline_by_key[key];old=baseline_source['records'][key]
        require(same((actual['address'],actual['original'],actual['analytic'],actual['unit'],actual['limits']),
                     (record['address'],old['original'],record['analytic'],record['unit'],record['limits'])), 'actual full baseline scalar caller/return')
        call=call_input(record,variables['frequency'],frequency,variables['origin'])
        owner={'kind':'accepted-baseline-frequency-binding','key':key,'packet':bound_route,'keys':['sourceJoins',list(baseline_by_key).index(key),'bound']}
        atlas.setdefault(hash(call['expression']),[]).append({'input':call,'value':actual['bound'],'owner':owner})
        atlas_routes.append(owner)
    journal.json('saved-scalar-binding-atlas.json',{'baselinePacket':bound_route,'owners':atlas_routes,'noIntermediateReconstruction':True})
    control = next(iter(atlas.values()))[0]['input']
    changed_coefficient = dict(control,expression=sp.Add(control['expression'],sp.S.One,evaluate=False))
    units = list(control['unit']);units[0] = units[0] + 1
    changed_unit = dict(control,unit=tuple(units) if isinstance(control['unit'],tuple) else units)
    mapping = dict(control['mapping']);mapping[variables['frequency']] = frequency + sp.Rational(1,100)
    changed_frequency = dict(control,mapping=mapping)
    mutations = {'original':control,'coefficient':changed_coefficient,'unit':changed_unit,'frequency':changed_frequency}
    journal.write('scalar-route-mutation-operands.pickle',mutations)
    responses = {name:not same(control,value) for name,value in mutations.items() if name!='original'}
    journal.json('scalar-route-controls.json',responses)
    require(all(responses.values()),'actual full scalar routing coefficient/unit/frequency controls')
    router=Router(journal,atlas);counts={} 
    inventory=reader.json(storage.READY/'completed-input-artifact-inventory.json',storage.READY_INVENTORY_SHA)
    maps_name='complex/frequency-end-maps.pickle';baseline_maps=reader.retain(mr/maps_name,matrix['artifacts'][maps_name]['sha256'])
    own_maps_name='maps/fine-value.pickle';own_maps=reader.retain(Path(pilot['runDirectory'])/own_maps_name,pilot['artifacts'][own_maps_name]['sha256'])
    for label in cp['cases']:
        ap='analytic-cases/'+label+'/frequency-analytic.pickle';spth='accepted/frequency-cases/'+label+'/frequency-source.pickle'
        chart=baseline if label==BASELINE else packet(ap)
        sources=baseline_source if label==BASELINE else packet(spth)
        common=packet('analytic-input-cases/'+label+'/common-input-pairs.pickle')
        for name in ('contextPair','basisPair','variablesPair','characterVariablesPair'):
            require(same(*common[name]),'saved complete shared native context/field/grade/variable inputs')
        journal.write('cases/'+label+'/unit-and-variable-input-pair.pickle',{
                      'own':{k:chart[k] for k in ('variables','fieldUnits','equationUnits','generators')},
                      'baseline':{k:baseline[k] for k in ('variables','fieldUnits','equationUnits','generators')},
                      'nativeBoundFieldUnits':bound['fieldUnits'],'nativeBoundRowUnits':bound['rowUnits']})
        require(all(same(chart[k],baseline[k]) for k in ('variables','fieldUnits','equationUnits','generators')), 'full actual shared source coordinate/unit/grade frame')
        # Native boundary units were serialized as SymPy Rational dimensions;
        # source units may contain Python ints. Retain both raw typed packets
        # and compare dimensions literally, without rewriting either packet.
        for own_units,bound_units in ((chart['fieldUnits'],bound['fieldUnits']),(chart['equationUnits'],bound['rowUnits'])):
            require(len(own_units)==len(bound_units)==5 and all(len(a)==len(b)==3 and all(x==y for x,y in zip(a,b))
                    for a,b in zip(own_units,bound_units)), 'literal native field/equation dimensions')
        physical=route('analytic-input-cases/'+label+'/physical-row-term-character-inputs.pickle')
        native_rows=route('accepted/frequency-cases/'+label+'/original-row-inputs.pickle')
        router.case=label;router.context=route('analytic-input-cases/'+label+'/common-input-pairs.pickle')
        router.analytic_route=route(ap);router.source_route=route(spth)
        journal.json('cases/'+label+'/input-routes.json',{'analytic':router.analytic_route,'source':router.source_route,'context':router.context,
                     'physicalRowsTermsCharacters':physical,'originalNativeRowsSourcesProfilesSettings':native_rows,
                     'actualFrequency':repr(frequency),'upstreamIndependentGradesPreserved':True})
        first=len(router.routes);new_before=router.new
        values=run(chart['records'],sources['records'],variables,frequency,router)
        values.update(frequency=frequency,variables=variables,fieldUnits=chart['fieldUnits'],equationUnits=chart['equationUnits'],
                      generators=chart['generators'],inputRoutes={'analytic':router.analytic_route,'source':router.source_route,'physical':physical,'nativeRows':native_rows},
                      scope='Scalar native frequency bindings only. No source-jet derivative, quadrature or row/matrix reuse accepted.')
        journal.write('cases/'+label+'/scalar-bindings.pickle',values)
        journal.write('cases/'+label+'/binding-routes.pickle',router.routes[first:])
        end_routes={}
        for end in ('LEFT','RIGHT'):
            name='cases/'+label+'/'+end.lower()+'/summary.json'
            summary=reader.json(storage.READY/'complete'/name,inventory['artifacts'][name]['sha256'])
            rows=summary['clusters'];owners={tuple(v['firstOwner'][:2]) for v in rows}
            require(len(rows)==5 and sum(v['nullity'] for v in rows)==7 and len(owners)==1,'whole saved end family caller')
            for row in rows:
                for field in ('fullOwnInput','context','rationalTable','acceptedSourceRoute'):
                    route_value=row[field]
                    require(reader.retain(route_value['logical'],route_value['sha256'])==route_value,'consumed full own end reference')
                selection_route=row['fullNativeSelection']['file']
                require(reader.retain(selection_route['logical'],selection_route['sha256'])==selection_route,'consumed whole end selection')
            owner=next(iter(owners))
            if owner==(BASELINE,end):
                require(all(any(v['completeSavedCallerMatch'] for v in row['baselineCandidates']) for row in rows),'full actual baseline caller matches')
                chosen={'file':baseline_maps,'keys':[end]}
            else:
                require(owner==(storage.OWNER,'RIGHT') and end=='RIGHT' and [v['index'] for v in rows]==[0,4,6,16,17], 'exact own pilot owner')
                chosen={'file':own_maps,'keys':[]}
            end_routes[end]={'map':chosen,'frequency':repr(frequency),'wholeInputJoin':reader.retain(storage.READY/'complete'/name),
                             'physicalSources':[v['fullOwnInput'] for v in rows],'selection':rows[0]['fullNativeSelection']['file'],
                             'scope':'Exact completed input family at this one frequency only; no source substitution.'}
        journal.json('cases/'+label+'/end-map-routes.json',end_routes)
        counts[label]={'records':len(values['joins']),'newScalarCalls':router.new-new_before,'reusedScalarCalls':len(values['joins'])-(router.new-new_before),
                       'rowsPending':cp['cases'][label]['rows'],'termsPending':cp['cases'][label]['terms'],'sourceJetsPending':cp['cases'][label]['sources'],
                       'endMaps':{k:v['map'] for k,v in end_routes.items()}}
        journal.json('cases/'+label+'/summary.json',counts[label])
        journal.json('case-inventory-'+str(len(counts))+'.json',counts)
    require(sum(v['records'] for v in counts.values())==1467,'complete own source scalar census')
    reader.postcheck()
    journal.json('inputs.json',{'acceptedAnalyticCheckpoint':reader.retain(CP),'acceptedEndPilot':reader.retain(PILOT),
                 'baselineFrequencyMatrix':reader.retain(MATRIX),'consumedRoutes':reader.routes,'allConsumedInputsUnchanged':True})
    result={'status':'COMPLETED_CASE_FREQUENCY_SCALAR_BINDINGS','cases':counts,'newScalarCalls':router.new,'totalPhysicalRecords':1467,
            'savedBaselineScalarRecords':len(baseline_by_key),'boundFrequency':repr(frequency),'endMapRoutes':8,'newEndMapCalls':0,
            'sourceAndInputsUnchanged':True,'artifacts':journal.artifacts,'wallSeconds':time.monotonic()-start,
            'scope':'One-frequency scalar bindings and exact saved end-map routes under toy-model scope. Source jets, full row contractions, interior operators, response and scoped searches remain.'}
    saved.save(base,'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
