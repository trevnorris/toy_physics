#!/usr/bin/env python3
"""Complete saved boundary validation with an explicit inherited-unit join.

The producer and failed validator are immutable. No scientific constructor or
current contractor is called. Original numerical packets are read unchanged.
"""
import ast
import copy
import json
from pathlib import Path
import shutil

import S11c_d_remaining_case_coordinate_boundary_loader_recover as recovery

h, f, modes = recovery.original, recovery.f, recovery.modes
ROOT = f.STORE/'s11c-remaining-case-coordinate-20260921/boundary'
OLD = ROOT/'production-acceptance'
BASE = ROOT/'production-recovery-02/complete'
OLD_SHA = '2a4068a6cd2033c509479801142b0dde259fcede9901d267d663095c69615a35'
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_units_finish_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_boundary_unit_repair.json'
LABEL = 'LAB_HELD__RHOBR_CONSTANT'
OWNER = LABEL+'__RIGHT'


def adapt():
    """Reverse the entire validator after only the documented continuation edits."""
    assert f.digest(OLD/'validate.py') == OLD_SHA
    old = ast.parse((OLD/'validate.py').read_text()); new = copy.deepcopy(old)
    fn = next(n for n in new.body if isinstance(n, ast.FunctionDef) and n.name == 'continuum_validate')
    original_fn = copy.deepcopy(fn)
    first = next(i for i, n in enumerate(fn.body) if isinstance(n, ast.Assign)
                 and ast.unparse(n.targets[0]) == 'dims')
    saved_prefix = copy.deepcopy(fn.body[:first])
    replacement = ast.parse('saved, material, common, proofs, X, U, T, pending, prepared, whole, inp, ordered, offsets = restore_continuum_prefix(folder, end, source, maps, signature)').body
    fn.body[:first] = replacement
    dims = fn.body[1]; assert ast.unparse(dims.value) == "f.unpickle(folder / 'dimension-state.pickle')"
    dims.value = ast.parse('joined_dimensions(folder, source, signature)', mode='eval').body
    counts = {'finiteCalls': 0, 'familyWrites': 0}
    class Calls(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'finite_validate':
                node.func.id = 'reuse_completed_finite'; counts['finiteCalls'] += 1
            elif ast.unparse(node.func) == 'f.save' and node.args and 'validated-family-' in ast.unparse(node.args[0]):
                node.func = ast.Name(id='retain_family_summary', ctx=ast.Load()); counts['familyWrites'] += 1
            return node
    new = Calls().visit(new)
    assert counts == {'finiteCalls': 1, 'familyWrites': 1}
    before = next(i for i, n in enumerate(new.body) if isinstance(n, ast.For)
                  and ast.unparse(n.target) == '(n, sha)' and ast.unparse(n.iter) == 'expected.items()')
    extra = ast.parse('register_continuation(expected, HERE)').body[0]; new.body.insert(before, extra)
    final = next(i for i, n in enumerate(new.body) if isinstance(n, ast.Expr)
                 and isinstance(n.value, ast.Call) and ast.unparse(n.value.func) == 'f.save'
                 and "HERE / 'checks.json'" == ast.unparse(n.value.args[0]))
    extension = ast.parse('finish_provenance(answer, HERE)').body[0]; new.body.insert(final, extension)
    reverse = copy.deepcopy(new)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name):
                if node.func.id == 'reuse_completed_finite': node.func.id = 'finite_validate'
                elif node.func.id == 'retain_family_summary': node.func = ast.parse('f.save', mode='eval').body
            return node
    reverse = Reverse().visit(reverse)
    reverse.body = [n for n in reverse.body if not (isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                    and isinstance(n.value.func, ast.Name) and n.value.func.id in ('register_continuation', 'finish_provenance'))]
    rf = next(n for n in reverse.body if isinstance(n, ast.FunctionDef) and n.name == 'continuum_validate')
    rf.body[1].value = ast.parse("f.unpickle(folder/'dimension-state.pickle')", mode='eval').body
    rf.body[:1] = saved_prefix
    assert ast.dump(rf) == ast.dump(original_fn)
    assert ast.dump(reverse) == ast.dump(old), 'whole validator reverse AST'
    proof = {'originalValidatorSha256': OLD_SHA, 'originalAST': h.ast_sha(old),
             'continuationAST': h.ast_sha(new), 'reverseWholeValidatorAST': True,
             'finiteValidationCallsReplacedByCompletedEvidence': 1,
             'continuumPrefixStatementsRestoredFromCompletedOperands': first,
             'unitStateReaderJoins': 1, 'familySummaryByteReuseCalls': 1,
             'provenanceObservations': 2, 'originalUnitAssertionUnchanged': True}
    return ast.fix_missing_locations(new), proof


def run(here):
    here = Path(here); tree, join = adapt(); spec = json.loads(REPAIR.read_text())
    assert join == spec['reverseAST']
    assert f.digest(Path(__file__)) == spec['helperSha256']
    ns = {'__file__': str(here/'validate.py'), '__name__': '__saved_boundary_validator__'}
    original_hashes = spec['failedValidatorFiles']; unit_outputs = {}; reuse = []
    native_material_ast = h.matrices.body(h.c.material_ends)
    native_current_ast = h.matrices.body(h.engine.UniformSlabCurrent)
    def same(a, b): assert modes.same(a, b)

    def register(expected, where):
        assert where == here
        old_inv = json.loads((OLD/'validate.invocation.json').read_text())
        old_guard = json.loads((OLD/'resource-guard/outcome.json').read_text())
        assert old_inv['exitCode'] == old_guard['exitCode'] == old_guard['childOutcome']['exitCode'] == 1
        assert old_guard['limitsVerified'] and old_guard['childOutcome']['guardReason'] is None
        assert not (OLD/'checks.json').exists()
        assert 'line 182, in continuum_validate' in (OLD/'validate.stderr').read_text()
        original = here/'original-validation'; original.mkdir()
        for name, sha in original_hashes.items():
            src = OLD/name; assert f.digest(src) == sha
            dst = original/name; dst.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(src, dst)
            assert f.digest(dst) == sha; expected[str(src)] = expected[str(dst)] = sha
        for path in (Path(__file__), PLAN, REPAIR):
            target = here/'source'/path.name; target.parent.mkdir(exist_ok=True); shutil.copyfile(path, target)
            expected[str(path)] = expected[str(target)] = f.digest(path)
        for name, sha in spec['diagnosticFiles'].items():
            src = ROOT/name; assert f.digest(src) == sha
            dst = here/'diagnostic-evidence'/name; dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src,dst); assert f.digest(dst) == sha
            expected[str(src)] = expected[str(dst)] = sha
        f.save(here/'continuation-join.json', join)

    def completed_finite(folder, signature, channels, maps):
        # The pinned original traceback lies after all three finite validators.
        # Restore their saved results and counters; do not repeat their arithmetic.
        owner = folder.parent.name
        assert owner in ('LEFT', 'RIGHT', OWNER)
        saved = f.unpickle(folder/'finite-material-boundary.pickle')
        current = f.unpickle(folder/'open-current-route.pickle')
        candidates = f.unpickle(folder/'transported-candidates.pickle')
        cc = ns['counts']; cc['finiteCandidates'] += len(candidates)
        cc['finiteBasisDirections'] += sum(int(v['INFO']['NULLITY']) for v in channels['CANDIDATES'])
        cc['finiteSourcePairs'] += len(current['pairs']); cc['finitePairEntries'] += int(current['coverage'].sum())
        cc['finiteControls'] += len(saved['controls'])
        ns['maxima']['finiteCurrent'] = max(ns['maxima']['finiteCurrent'], h.b.norm(current['differences']))
        # The prior in-memory source-pullback maxima were not checkpointed.
        # Their passing guards are retained; do not report invented zero maxima.
        ns['maxima']['sourcePullback'] = None; ns['maxima']['sourcePullbackScaled'] = None
        reuse.append({'owner': owner, 'resultSha256': f.digest(folder/'finite-material-boundary.pickle'),
                      'originalFunctionAST': h.ast_sha(next(n for n in ast.parse((OLD/'validate.py').read_text()).body
                                                          if getattr(n, 'name', None) == 'finite_validate')),
                      'completionWitness': 'Pinned original validator reached continuum unit assertion after this full finite call returned.'})
        f.save(here/'completed-finite-validation-reuse.json', reuse)
        return saved

    def restore_prefix(folder, end, source, maps, signature):
        assert folder == BASE/'material-families'/OWNER/'continuum' and end == 'RIGHT'
        saved = f.unpickle(folder/'right-material-boundary.pickle')
        material, common, proofs = saved['material'], saved['commonEulerian'], saved['proofs']; X, U, T = maps
        pending = f.unpickle(folder/'native-checkpoints/RIGHT-complete-before-guard.pickle')
        prepared = f.unpickle(folder/'native-checkpoints/RIGHT-prepared.pickle')
        whole = f.unpickle(folder/'material-boundary.pickle'); inp = f.unpickle(folder/'native-input.pickle')
        ordered = [v for direction in ('outgoing','incoming') for v in material['clusters'] if v['info']['direction'] == direction]
        offsets = source['offsets']; assert int(offsets[-1]) == 7
        f.save(here/'completed-continuum-prefix-reuse.json', {
            'originalValidatorSha256': OLD_SHA, 'failedLine': 182,
            'prefixAlreadyPassed': True, 'savedOffsetsSource': 'actual accepted own continuum offsets',
            'originalArraysRetained': True, 'newMapOrCurrentConstruction': False})
        return saved, material, common, proofs, X, U, T, pending, prepared, whole, inp, ordered, offsets

    def joined(folder, source, signature):
        target = here/'unit-overlay'; target.mkdir()
        raw_path = folder/'dimension-state.pickle'; raw = f.unpickle(raw_path)
        mcp_path = f.M/'S11c_d_remaining_case_modes_checkpoint.json'
        bcp_path = f.M/'S11c_d_remaining_case_boundary_checkpoint.json'
        checkpoints = [json.loads(p.read_text()) for p in (mcp_path, bcp_path)]
        assert [v['status'] for v in checkpoints] == ['ACCEPTED_FOUR_CASE_MODE_SUBSPACES','ACCEPTED_FOUR_CASE_BOUNDARY_MAPS']
        mr, br = [Path(v['runDirectory']) for v in checkpoints]
        mc, bc = [json.loads((r/'checks.json').read_text()) for r in (mr, br)]
        for r, cp in zip((mr, br), checkpoints): assert f.digest(r/'checks.json') == cp['checksSha256']
        modal_path = mr/'cases'/LABEL/'right/modal.pickle'
        boundary_path = br/'boundary-cases'/LABEL/'right/continuum-boundary.pickle'
        sha = f.digest(modal_path)
        assert sha == mc['artifacts'][str(modal_path.relative_to(mr))]['sha256'] == bc['inputPackets'][str(modal_path)] == ns['manifest']['inputPackets'][str(modal_path)]
        assert f.digest(boundary_path) == bc['artifacts'][str(boundary_path.relative_to(br))]['sha256']
        same(f.unpickle(boundary_path), source)
        modal, known = f.unpickle(modal_path)
        all_symbols = set().union(*(v.free_symbols for n, v in modal['SYMBOLIC_OPERANDS'].items()
                                  if n in ('CURRENT_SLAB','CURRENT_BULK')))
        tables = [f.unpickle(folder/('right-'+n+'-material-current-tables.pickle')) for n in ('slab','bulk')]
        arguments = tables[0]['originalArguments']; new_arguments = tables[0]['materialArguments']
        assert len(set(arguments)) == 4 and set(arguments) <= all_symbols
        kl, kr, ql, qr = arguments; ml, mm, _, _ = new_arguments
        assert tuple(str(v) for v in arguments) == ('s11cdCurrentLeftMomentum','s11cdCurrentRightMomentum','s11cdAcousticLeftNormalMomentum','s11cdAcousticRightNormalMomentum')
        for table, name in zip(tables, ('slab','bulk')):
            same(table['original'], source['currentTables'][name]); assert table['originalArguments'] == arguments and table['materialArguments'] == new_arguments
        def check_units(candidate_known, candidate_raw, candidate_arguments, candidate_material):
            assert candidate_arguments == arguments and candidate_material == new_arguments
            for atom in arguments:
                assert atom in candidate_known and tuple(candidate_known[atom]) == tuple(known[atom]) == (-1,0,0)
                assert atom not in candidate_raw['known']
            for old, new in ((kl,ml),(kr,mm)):
                assert tuple(candidate_raw['known'][new]) == tuple(candidate_raw['unknown'][old])
                assert new.assumptions0 == old.assumptions0
                assert all(v not in candidate_raw['solution'] for v in candidate_raw['unknown'][old])
        check_units(known, raw, arguments, new_arguments)
        controls = []
        for atom in arguments:
            changed = dict(known); changed[atom] = (known[atom][0]+1,*known[atom][1:])
            try: check_units(changed, raw, arguments, new_arguments)
            except AssertionError: controls.append({'kind':'source-unit','atom':atom,'before':known[atom],'after':changed[atom],'rejected':True})
            else: raise AssertionError('changed physical source unit accepted')
        for old, new in ((kl,ml),(kr,mm)):
            changed = dict(raw, known=dict(raw['known'])); changed['known'][new] = raw['known'][new][::-1]
            try: check_units(known, changed, arguments, new_arguments)
            except AssertionError: controls.append({'kind':'native-alias','atom':new,'before':raw['known'][new],'after':changed['known'][new],'rejected':True})
            else: raise AssertionError('changed native unit alias accepted')
        for args, mat in ((arguments[::-1],new_arguments),(arguments,new_arguments[::-1])):
            try: check_units(known,raw,args,mat)
            except AssertionError: controls.append({'kind':'argument-address','original':(arguments,new_arguments),'changed':(args,mat),'rejected':True})
            else: raise AssertionError('wrong current argument address accepted')
        overlay = {v:known[v] for v in arguments}; overlay.update({ml:known[kl],mm:known[kr]})
        completed = dict(raw, known=dict(raw['known'])); completed['known'].update(overlay)
        for k, value in raw.items():
            if k != 'known': same(completed[k],value)
        same({k:v for k,v in completed['known'].items() if k not in overlay},
             {k:v for k,v in raw['known'].items() if k not in overlay})
        f.atomic_pickle(target/'unit-source-pairs.pickle', {'rawKnown':{v:raw['known'].get(v) for v in overlay},
            'rawUnknown':{v:raw['unknown'].get(v) for v in overlay},'acceptedSourceUnits':{v:known[v] for v in arguments},
            'nativeAssignments':((ml,kl),(mm,kr)),'joinedKnown':overlay,'sourceArguments':arguments,'materialArguments':new_arguments,
            'acceptedModal':str(modal_path),'actualBoundarySource':str(boundary_path)})
        f.atomic_pickle(target/'dimension-state-view.pickle', completed)
        f.atomic_pickle(target/'mutation-controls.pickle', controls)
        # The overlay is a supplemental metadata view. No producer packet changes.
        metadata = {'status':'JOINED_ACCEPTED_CURRENT_UNIT_DECLARATIONS','originalDimensionState':str(raw_path),
                    'originalDimensionStateSha256':f.digest(raw_path),'acceptedModal':str(modal_path),'acceptedModalSha256':sha,
                    'acceptedBoundary':str(boundary_path),'acceptedBoundarySha256':f.digest(boundary_path),
                    'sourceArguments':4,'materialAliases':2,'respondingMutations':len(controls),
                    'nativeCurrentClassAST':native_current_ast,
                    'nativeMaterialEndAST':native_material_ast,
                    'originalUnknownHistoryRetained':True,'changedNumericalOrSymbolicCurrentOperands':0}
        f.save(target/'checks.json',metadata)
        for p in (mcp_path,bcp_path,mr/'checks.json',br/'checks.json',modal_path,boundary_path):
            digest=f.digest(p); assert str(p) not in ns['expected'] or ns['expected'][str(p)]==digest
            ns['expected'][str(p)]=digest
        for p in target.iterdir():
            unit_outputs[str(p.relative_to(here))]={'sha256':f.digest(p),'bytes':p.stat().st_size}
            ns['expected'][str(p)]=f.digest(p)
        return completed

    def retain_summary(path, value):
        old = OLD/path.name
        if old.exists():
            assert json.loads(old.read_text()) == value; shutil.copyfile(old,path); assert f.digest(path)==f.digest(old)
        else: f.save(path,value)

    def finish(answer, where):
        assert where == here and len(reuse)==3 and unit_outputs
        answer['metadataSupplement']={'runDirectory':str(here),'artifacts':unit_outputs,
            'scope':'Four accepted current momentum declarations and their two exact native material aliases; raw producer metadata and all physical operands remain immutable.'}
        answer['continuation']={'originalValidatorDirectory':str(OLD),'originalValidatorSha256':OLD_SHA,
            'reverseAST':join,'completedFiniteFamiliesReused':3,'completedBaselineFamilySummariesReused':2,
            'completedNewContinuumPrefixReused':True,'sourcePullbackMaximumNotCheckpointed':True}
        answer['sourceFiles'] = dict(answer['sourceFiles'])
        for path in (Path(__file__),PLAN,REPAIR):
            answer['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    ns.update(register_continuation=register,reuse_completed_finite=completed_finite,
              restore_continuum_prefix=restore_prefix,joined_dimensions=joined,
              retain_family_summary=retain_summary,finish_provenance=finish)
    exec(compile(tree,str(OLD/'validate.py'),'exec'),ns)
