#!/usr/bin/env python3
"""Saved-operator and own-case endpoint preparation for derivative responses."""
import argparse, copy, gc, json, resource, shutil, signal, time
from pathlib import Path
import numpy as np
import S11c_d_remaining_case_first_jet_matrices as matrices
import S11c_d_remaining_case_response as response

f, m, b = matrices.f, matrices.modes, response.b
provenance = matrices.binding.source.provenance
BASELINE = matrices.BASELINE
ICP = f.M/'S11c_d_remaining_case_first_jet_matrices_checkpoint.json'
RCP = f.M/'S11c_d_remaining_case_response_checkpoint.json'
ECP = f.M/'S11c_d_remaining_case_boundary_checkpoint.json'
JCP = matrices.JCP
PLAN = f.M/'S11c_d_remaining_case_first_jet_response_plan.md'
SCOPE = ('Complete saved first-w-derivative operators and actual own-case end/current/phase inputs. '
         'The accepted baseline control is continuum only. Four finite control solves and three missing-case continuum solves remain new work. '
         'No new integration, source assembly, mode/current construction or response solve in this preparation.')


def load(base):
    accepted = [provenance.accepted(path, status) for path, status in (
        (ICP, 'ACCEPTED_CASE_FIRST_JET_INTERIOR_MATRICES'),
        (RCP, 'PUBLISHED_ANNEX_VERIFIED'),
        (ECP, 'ACCEPTED_FOUR_CASE_BOUNDARY_MAPS'),
        (JCP, 'PUBLISHED_ANNEX_VERIFIED'))]
    (ir, ic, _), (rr, rc, _), (er, ec, _), (jr, jc, _) = accepted
    pins = {}; inputs = {}
    for root, checks, cp in accepted:
        for name, value in checks['sourceFiles'].items():
            f.require(name not in pins or pins[name] == value, 'same consumed physical/helper source'); pins[name] = value
        inputs[str(root/'checks.json')] = cp['checksSha256']
        for name, value in checks['inputPackets'].items():
            f.require(name not in inputs or inputs[name] == value, 'same original consumed operand'); inputs[name] = value
    for path in (Path(__file__), PLAN, ICP, RCP, ECP, JCP, Path(response.__file__), Path(matrices.__file__)):
        pins[str(path.resolve().relative_to(f.ROOT))] = f.digest(path)
    manifest = {'runDirectory': str(base), 'sourceFiles': pins, 'inputPackets': inputs, 'copiedInputs': {},
                'input': ic['input'], 'settings': ic['settings'], 'scope': SCOPE}
    labels = tuple(ic['preflight']['cases'])
    def retain(root, checks, name, target):
        f.require(name in checks['artifacts'], ('actual accepted artifact address', name))
        m.retain(root/name, base/target, manifest, checks['artifacts'][name]['sha256'])
    retain(ir, ic, 'accepted-finite-system.pickle', 'interiors/accepted-finite-system.pickle')
    retain(ir, ic, 'preflight.json', 'accepted-matrix-preflight.json')
    retain(ir, ic, 'selected-cell-input-joins.pickle', 'selected-cell-input-joins.pickle')
    for name in ('first-jet-binding.pickle', 'first-jet-rows.pickle', 'first-jet-matrices.pickle',
                 'first-jet-systems.pickle', 'first-jet-solutions.pickle', 'first-jet-channels.pickle',
                 'first-jet-comparisons.pickle', 'first-jet-response.pickle', 'full.out'):
        retain(jr, jc, name, 'accepted-first-jet/'+name)
    for label in labels:
        for kind in ('reduced-action', 'actions', 'assembly', 'factorization'):
            name = 'accepted-cases/'+label+'/'+kind+'.pickle'
            retain(ir, ic, name, 'interiors/'+name)
        for part, names in (('finite', ('finite-system', 'finite-solution', 'observable')),
                            ('continuum', ('coefficient-systems', 'coefficient-solutions', 'channel-response', 'continuum-response', 'formal-remainders'))):
            for name in names:
                address = 'cases/'+label+'/'+part+'/'+name+'.pickle'
                retain(rr, rc, address, 'unchanged-response/'+address)
        address = 'boundary-cases/'+label+'/case-boundary.pickle'
        retain(rr, rc, address, 'unchanged-response/'+address)
        retain(er, ec, address, address)
        original = 'binding-inputs/accepted-bindings/'+label+'/case-binding.pickle'
        retain(ir, ic, original, 'original-bindings/'+label+'/case-binding.pickle')
        if label == BASELINE:
            retain(ir, ic, 'binding-inputs/baseline-first-jet-binding-view.pickle', 'baseline-binding-view.pickle')
            retain(rr, rc, 'interiors/cases/'+label+'/interior-matrices.pickle', 'original-baseline-interior.pickle')
        else:
            for name in ('interior-matrices', 'coefficient-matrices', 'direct-native-cells', 'direct-unsplit', 'row-matrices', 'comparisons'):
                address = 'cases/'+label+'/'+name+'.pickle'; retain(ir, ic, address, 'interiors/'+address)
            retain(ir, ic, 'binding-inputs/bindings/'+label+'/case-binding.pickle', 'interiors/accepted-bindings/'+label+'/case-binding.pickle')
            for name in ('first-jet-sources', 'term-inputs', 'endpoint-operands'):
                retain(ir, ic, 'binding-inputs/cases/'+label+'/'+name+'.pickle', 'selected-sources/'+label+'/'+name+'.pickle')
    retain(er, ec, 'accepted-modes/reference/modal.pickle', 'accepted-modes/reference/modal.pickle')
    for name, value in pins.items():
        dest = base/'source'/name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(f.ROOT/name, dest)
        f.require(f.digest(dest) == value, 'frozen response-input source')
    basis = f.unpickle(base/'interiors/accepted-finite-system.pickle')
    manifest['settings'] = provenance.restore_settings(manifest['settings'], basis['settings'])
    f.save(base/'inputs.json', manifest)
    return manifest, labels


def baseline_views(base, manifest):
    saved = f.unpickle(base/'accepted-first-jet/first-jet-response.pickle')
    original = f.unpickle(base/'original-bindings'/BASELINE/'case-binding.pickle')
    binding = f.unpickle(base/'accepted-first-jet/first-jet-binding.pickle')
    packet = f.unpickle(base/'accepted-first-jet/first-jet-matrices.pickle')
    prior = f.unpickle(base/'original-baseline-interior.pickle')
    rows = f.unpickle(base/'accepted-first-jet/first-jet-rows.pickle')
    # These are storage views of accepted arrays and source records. No assembler
    # is called, and no direct native-cell proof is invented for the old control.
    grade = dict(original['grades'], records=binding['records'], termJoins=binding['termJoins'])
    case = {'binding': f.unpickle(base/'baseline-binding-view.pickle'), 'grades': grade,
            'newRows': [], 'reusedRows': [], 'endpoints': binding['endpointChecks']}
    for key in ('fieldUnits', 'equationUnits'):
        f.require(m.same(case['binding'][key], prior[key]), 'baseline physical unit input')
    f.require(m.same(packet['matrices'], saved['matrices']) and m.same(rows['rows'], saved['rows']), 'entire accepted baseline operator/row arrays')
    f.require(m.same(saved['generators'], grade['generators']), 'accepted independent generator order')
    interior = {k: prior[k] for k in ('size', 'gradeOrigin', 'blockUnits', 'fieldUnits', 'equationUnits', 'settings')}
    interior.update(matrices=packet['matrices'], generators=grade['generators'], dimensionState=saved['dimensionState'],
                    comparisons={'approved_total': packet['recombination']}, scope='View of accepted baseline first-derivative arrays; no new assembly or finite response.')
    target = base/'interiors/cases'/BASELINE; target.mkdir(parents=True)
    for name, value in (('interior-matrices', interior), ('direct-unsplit', packet['unsplit']),
                        ('row-matrices', dict(rows, settings=manifest['settings']))):
        f.atomic_pickle(target/(name+'.pickle'), value)
    p = base/'interiors/accepted-bindings'/BASELINE; p.mkdir(parents=True)
    f.atomic_pickle(p/'case-binding.pickle', case)
    record = {'parents': {str(p.relative_to(base)): f.digest(p) for p in (
        base/'accepted-first-jet/first-jet-matrices.pickle', base/'accepted-first-jet/first-jet-rows.pickle',
        base/'accepted-first-jet/first-jet-response.pickle', base/'accepted-first-jet/first-jet-binding.pickle',
        base/'original-baseline-interior.pickle', base/'baseline-binding-view.pickle')},
        'newArrayEntries': 0, 'newAssemblies': 0, 'baselineContinuumReused': True, 'baselineFiniteSolveExists': False}
    f.save(base/'baseline-storage-views.json', record)


def prepare(base, manifest, labels):
    baseline_views(base, manifest); gc.collect()
    basis = f.unpickle(base/'interiors/accepted-finite-system.pickle')
    finish, finite_prefix, finite_join = response.finite_tail()
    summaries = {}; total_rows = total_terms = 0
    for label in labels:
        target = base/'preparation'/label; target.mkdir(parents=True)
        source = base/'interiors/cases'/label
        case = f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle')
        original = f.unpickle(base/'original-bindings'/label/'case-binding.pickle')
        coefficients = f.unpickle(source/'interior-matrices.pickle')
        direct = response.unsplit_arrays(f.unpickle(source/'direct-unsplit.pickle'))
        rows = f.unpickle(source/'row-matrices.pickle')
        ends = f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle')
        old_ends = f.unpickle(base/'unchanged-response/boundary-cases'/label/'case-boundary.pickle')
        old_system = f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-system.pickle')
        old_solution = f.unpickle(base/'unchanged-response/cases'/label/'finite/finite-solution.pickle')
        f.require(m.same(ends, old_ends) and m.same(old_system['channels'], ends['finite']), 'complete own-case finite/continuum end/current/phase coordinates')
        for key in ('nodes', 'derivativeMatrices', 'settings'):
            f.require(m.same(basis[key], old_system[key]), 'actual common finite basis and settings')
        f.require(case['binding']['settings'] == rows['settings'] == coefficients['settings'] == basis['settings'], 'all actual numerical settings')
        for key, ek in (('fieldUnits', 'fieldUnits'), ('equationUnits', 'rowUnits')):
            f.require(m.same(case['binding'][key], original['binding'][key]) and m.same(case['binding'][key], coefficients[key])
                      and m.same(coefficients[key], ends[ek]), 'actual inherited field/equation units')
        endpoints = case['endpoints']
        f.require(len(endpoints) == 2 and all(v['derivativeLimit'] == v['reversedDerivativeLimit'] == 0 for v in endpoints), 'actual original and reversed endpoint derivative operands')
        f.require({v['profileLimit'] for v in endpoints} == {0, 1}, 'actual profile endpoints stay unchanged')
        field_rows = case['binding']['bound']['rows']; old_rows = original['binding']['bound']['rows']
        f.require(len(field_rows) == len(old_rows) and set(rows['rows']) == set(range(len(field_rows))), 'complete actual source rows')
        for row, old in zip(field_rows, old_rows):
            for key in ('original', 'symbolicLimits', 'unit'):
                f.require(m.same(row[key], old[key]), 'original physical row identity and ordered limits')
        f.require(len(case['grades']['termJoins']) == len(original['grades']['termJoins']), 'complete native terms retained')
        total_rows += len(field_rows); total_terms += len(case['grades']['termJoins'])
        a, rhs, retained = finite_prefix(target, direct['total'].copy(), ends['finite'], len(basis['nodes']), basis['derivativeMatrices'],
            basis['nodes'], direct['local'], basis['settings'], manifest['sourceFiles'], manifest['inputPackets'], old_solution['polynomialDerivativeResiduals'], [rows])
        f.require(np.array_equal(retained, direct['total']) and a.shape == (645,645) and rhs.shape == (645,4)
                  and np.isfinite(a).all() and np.isfinite(rhs).all(), 'full finite boundary-replaced system without solve')
        wrong = copy.deepcopy(ends['finite']); wrong['LEFT']['incomingBoundaryData'] *= -1
        _, wrong_rhs, _ = finite_prefix(target, direct['total'].copy(), wrong, len(basis['nodes']), basis['derivativeMatrices'],
            basis['nodes'], direct['local'], basis['settings'], {}, {}, old_solution['polynomialDerivativeResiduals'], [rows])
        mutation = wrong_rhs-rhs; f.require(b.norm(mutation) > 0, 'actual one-sided finite incident forcing control')
        system = {'matrix': a, 'rhs': rhs, 'unreplacedOperator': retained, 'incidentSignMutation': mutation}
        f.atomic_pickle(target/'finite-prepared.pickle', system)
        continuum, forcing = response.response.systems(coefficients, ends, basis)
        f.require(set(continuum) == set(forcing) == set(b.G) and all(x.shape == (645,645) and np.isfinite(x).all() for x in continuum.values())
                  and all(x.shape == (645,4) and np.isfinite(x).all() for x in forcing.values()), 'complete independent-grade matrices and actual forcing')
        f.atomic_pickle(target/'continuum-prepared.pickle', {'matrices': continuum, 'rhs': forcing})
        if label == BASELINE:
            old = f.unpickle(base/'accepted-first-jet/first-jet-systems.pickle')
            solved = f.unpickle(base/'accepted-first-jet/first-jet-solutions.pickle')
            fresh = f.unpickle(base/'accepted-first-jet/first-jet-response.pickle')
            unchanged = f.unpickle(base/'unchanged-response/cases'/label/'continuum/continuum-response.pickle')
            f.require(m.same(fresh['baseline']['response'], unchanged['response']), 'accepted baseline control uses exact original response coordinates')
            delta = b.subtract(continuum, old['matrices']); fd = b.subtract(forcing, old['rhs'])
            residual = b.subtract(b.J.multiply(continuum, solved['coefficients']), forcing)
            scaled = {g: v/solved['rowScale'][:,None] for g,v in residual.items()}
            f.require(b.norm(delta)==0 and b.norm(fd)==0 and b.norm(scaled)<1e-9, 'exact baseline continuum system and saved solve reuse')
            f.atomic_pickle(target/'baseline-continuum-replay.pickle', {'matrixDifferences':delta,'forcingDifferences':fd,'savedSolutionScaledResidual':scaled})
            del old,solved,fresh,unchanged,delta,fd,residual,scaled
        response.emitter('FIRST_JET_'+label)  # Compile the exact native continuum replay tail; no emission occurs here.
        summaries[label] = {'rows':len(field_rows),'terms':len(case['grades']['termJoins']),'sources':len(case['binding']['jets']),
            'unknowns':645,'incomingColumns':4,'sameCompleteEndInputs':True,'endpointDerivativeZeros':True,
            'incidentSignMutation':b.norm(mutation),'baselineContinuumReuse':label==BASELINE,
            'missingFiniteControlSolve':True,'missingContinuumControlSolve':label!=BASELINE}
        f.save(base/'case-inventory.json', summaries)
        del case,original,coefficients,direct,rows,ends,old_ends,old_system,old_solution,field_rows,old_rows,a,rhs,retained,system,continuum,forcing,wrong,wrong_rhs,mutation;gc.collect()
    f.require(total_rows==300 and total_terms==647, 'complete four-case selected source census')
    result={'cases':summaries,'nativeFiniteTail':finite_join,'continuumFunctions':{name:matrices.binding.native_body(getattr(response.response,name)) for name in ('systems','solve','channels','open_flux')},
            'nativeContinuumReplayTails':4,'rowAddresses':total_rows,'nativeTerms':total_terms,'newQuadratureNodes':0,'newInteriorAssemblies':0,'newModes':0,'newCurrentClosures':0,'newResponseSolves':0,
            'remainingFiniteControlSolves':4,'remainingContinuumControlSolves':3,'scope':SCOPE}
    f.save(base/'preflight.json', result);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',required=True,type=Path);args=p.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    manifest,labels=load(base);prepared=prepare(base,manifest,labels)
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen source pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'byte-identical accepted operand copies')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    result={**manifest,'status':'COMPLETED_FIRST_JET_RESPONSE_INPUTS','preflight':prepared,'cases':prepared['cases'],'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
