#!/usr/bin/env python3
"""New bounded centre-load diagnostic on exact saved b/c2 operands.

No producer is imported. No kernel, projection, integral, derivative, mode or
response is reconstructed. Interpretation belongs in the companion report.
"""
import argparse
import ast
from fractions import Fraction
import gc
import hashlib
import json
import os
import pickle
from pathlib import Path
import time

import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

M = Path(__file__).resolve().parent
ROOT = M.parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
        if hasattr(os, 'posix_fadvise'):
            os.posix_fadvise(f.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    return h.hexdigest()


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def save_value(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as f:
        pickle.dump(value, f, protocol=4)
    if isinstance(value, sp.Basic) and path.stat().st_size < 65536:
        with path.with_suffix('.srepr').open('x') as f:
            f.write(sp.srepr(value) + '\n')
    return {'path': str(path), 'sha256': digest(path), 'bytes': path.stat().st_size}


def file_route(path):
    path = Path(path)
    links = []
    for part in (*reversed(path.parents), path):
        if part.is_symlink():
            links.append({'path': str(part), 'target': os.readlink(part)})
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'links': links, 'sha256': digest(path)}


class Inputs:
    def __init__(self, specification, out):
        self.out = out
        self.specification = specification
        self.pins = {r['path']: r for r in specification['files']}
        self.seen = {}
        self.restored = {}
        self.codec = None

    def pin(self, path):
        path = str(path)
        if path not in self.seen:
            expected = self.pins[path]
            actual = file_route(path)
            save_json(self.out/'consumed'/f'{len(self.seen):03d}.json', actual)
            if actual != expected:
                raise ValueError(('consumed file changed', path))
            self.seen[path] = actual
        return Path(path)

    def json(self, path):
        return json.loads(self.pin(path).read_text())

    def literal(self, route):
        path = self.pin(route['record']['path'])
        if str(path.resolve(strict=True)) != route['record']['canonicalPath']:
            raise ValueError('changed canonical saved record')
        with path.open('rb') as f:
            f.seek(route['byteOffset']); raw = f.read(route['bytes'])
        if hashlib.sha256(raw).hexdigest() != route['sha256']:
            raise ValueError(('saved field changed', route['selector']))
        return raw.decode()

    def restore(self, route):
        key = (route['record']['path'], route['byteOffset'], route['bytes'], route['sha256'])
        if key not in self.restored:
            text = self.literal(route)
            # Same native generated-export codec; no scientific producer main.
            value = self.codec(text)
            self.restored[key] = value
        return self.restored[key]

    def postcheck(self):
        actual = {p: file_route(Path(p)) for p in self.seen}
        save_json(self.out/'posthashes.json', actual)
        if actual != self.seen:
            raise ValueError('consumed path changed')
        return len(actual)


def source_codec(source):
    # Compile the literal native _RELATIONALS and _restore definitions only.
    selected = []
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '_RELATIONALS' for t in node.targets):
            selected.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name == '_restore':
            selected.append(node)
    if len(selected) != 2:
        raise ValueError('native codec source shape changed')
    module = ast.Module(body=selected, type_ignores=[])
    namespace = {'sp': sp, 'Str': Str, 'ExprCondPair': ExprCondPair}
    exec(compile(module, '<saved-native-codec>', 'exec'), namespace)
    return namespace['_restore']


def raw_coeff(expr, slot):
    return expr.coeff(slot)


def zero_property(value):
    return value.is_zero


def nonzero_property(value):
    return value.is_nonzero


def structural_equal(left, right):
    return left == right


def raw_divide(left, right):
    return left / right


def replace(expr, replacements):
    return expr.xreplace(replacements)


FUNCTIONS = {'coeff': raw_coeff, 'add': sp.Add, 'mul': sp.Mul,
             'divide': raw_divide, 'cancel': sp.cancel, 'replace': replace,
             'expand': sp.expand, 'fraction': sp.fraction,
             'zero_property': zero_property, 'nonzero_property': nonzero_property,
             'structural_equal': structural_equal}


class Journal:
    def __init__(self, out, context):
        self.out = out
        self.context = context
        self.cache = {}
        self.count = 0
        self.receipts = []

    def op(self, label, operation, *args):
        directory = self.out/'operations'/f'{self.count:04d}'
        self.count += 1
        inp = {'operation': operation, 'args': args, 'physicalContext': self.context}
        input_route = save_value(directory/'input.pickle', inp)
        key = (operation, input_route['sha256'])
        if key in self.cache:
            value, owner = self.cache[key]
            result_route = owner['result']
            decision = 'reuse-exact-preceding-new-call'
        else:
            value = FUNCTIONS[operation](*args)
            result_route = save_value(directory/'value.pickle', value)
            owner = None
            decision = 'new'
        receipt = {'label': label, 'operation': operation, 'decision': decision,
                   'input': input_route, 'result': result_route,
                   'ownerReceipt': owner['receipt'] if owner else None,
                   'receipt': str(directory/'completed.json')}
        save_json(directory/'completed.json', receipt)
        if owner is None:
            self.cache[key] = (value, receipt)
        self.receipts.append(receipt)
        return value


def scalar_units(expr, schema):
    """Dimension metadata of the newly computed scalar, using native schema."""
    if expr.is_Number:
        return (Fraction(0),) * 3
    if expr.is_Symbol:
        return tuple(Fraction(v) for v in schema[expr.name])
    if expr.is_Mul:
        children = [scalar_units(v, schema) for v in expr.args]
        return tuple(sum(v[i] for v in children) for i in range(3))
    if expr.is_Pow and expr.exp.is_Rational:
        factor = Fraction(int(expr.exp.p), int(expr.exp.q))
        return tuple(factor*v for v in scalar_units(expr.base, schema))
    if expr.is_Add:
        children = [scalar_units(v, schema) for v in expr.args]
        if not all(v == children[0] for v in children):
            raise ValueError('scalar units inconsistent')
        return children[0]
    raise ValueError(('unsupported scalar dimensional form', expr.func))


def split_slots(j, label, expr, slots):
    coefs = tuple(j.op(label+'/coefficient/'+str(i), 'coeff', expr, s) for i, s in enumerate(slots))
    terms = tuple(j.op(label+'/slot-term/'+str(i), 'mul', c, s) for i, (c, s) in enumerate(zip(coefs, slots)))
    content = j.op(label+'/slot-content', 'add', *terms)
    neg = j.op(label+'/negative-content', 'mul', sp.S.NegativeOne, content)
    remainder_raw = j.op(label+'/non-slot-remainder-input', 'add', expr, neg)
    remainder = j.op(label+'/non-slot-remainder', 'expand', remainder_raw)
    flags = [j.op(label+'/remainder-slot/'+str(i), 'coeff', remainder, s) for i, s in enumerate(slots)]
    statuses = [j.op(label+'/remainder-slot-zero/'+str(i), 'zero_property', value) for i, value in enumerate(flags)]
    reconstruction_input = j.op(label+'/reconstruction-input', 'add', content, remainder)
    reconstruction = j.op(label+'/reconstructed', 'expand', reconstruction_input)
    identity = j.op(label+'/reconstruction-identity', 'structural_equal', reconstruction, expr)
    save_json(j.out/(label+'-split.json'), {'slotCoefficientsRecorded': len(coefs),
              'remainderSlotZeros': statuses, 'reconstructionIdentity': identity})
    return coefs, remainder, all(s is True for s in statuses) and identity


def derive(j, label, centre, thickness, slots, centre_dim, face_dim, schema, scalar_roles):
    c, outside, c_ok = split_slots(j, label+'/centre', centre, slots)
    e, e_outside, e_ok = split_slots(j, label+'/thickness', thickness, slots)
    outside_zero = j.op(label+'/centre-non-slot-zero', 'zero_property', outside)
    ratios = []; face_records = []
    for face_index in range(2):
        ci = c[2*face_index:2*face_index+2]; ei = e[2*face_index:2*face_index+2]
        candidate = None; coefficient_slot = None
        for k in range(2):
            flag = j.op(label+f'/face{face_index}/denominator-zero{k}', 'zero_property', ei[k])
            if flag is not True:
                coefficient_slot = k
                raw = j.op(label+f'/face{face_index}/candidate-quotient', 'divide', ci[k], ei[k])
                candidate = j.op(label+f'/face{face_index}/reduced-candidate', 'cancel', raw)
                break
        if candidate is None:
            face_records.append({'candidateAvailable': False}); ratios.append(None); continue
        residuals = []
        for k in range(2):
            product = j.op(label+f'/face{face_index}/product{k}', 'mul', candidate, ei[k])
            neg = j.op(label+f'/face{face_index}/negative-product{k}', 'mul', sp.S.NegativeOne, product)
            difference = j.op(label+f'/face{face_index}/difference{k}', 'add', ci[k], neg)
            residual = j.op(label+f'/face{face_index}/residual{k}', 'cancel', difference)
            residuals.append(j.op(label+f'/face{face_index}/residual-zero{k}', 'zero_property', residual))
        names = sorted(a.name for a in candidate.free_symbols)
        scalar_ok = set(names) <= set(scalar_roles) and not candidate.atoms(sp.Function, sp.Integral, sp.Derivative)
        # Candidate division is discovery only: coefficient identities, including
        # degenerate numerator branches, determine validity of the reduced map.
        numerator, denominator = j.op(label+f'/face{face_index}/candidate-fraction', 'fraction', candidate)
        denominator_route = save_value(j.out/(label+f'/face{face_index}-candidate-denominator.pickle'), denominator)
        domain = j.op(label+f'/face{face_index}/candidate-domain', 'nonzero_property', denominator)
        try:
            r_units = scalar_units(candidate, schema)
            mapped_units = tuple(a+b for a, b in zip(r_units, face_dim))
            units_ok = mapped_units == centre_dim
            units = {'ratio': [str(v) for v in r_units], 'mapped': [str(v) for v in mapped_units]}
        except (KeyError, ValueError):
            units_ok = False; units = None
        ok = c_ok and e_ok and all(v is True for v in residuals) and scalar_ok and units_ok and domain is True
        face_records.append({'candidateAvailable': True, 'candidateCoefficientSlot': coefficient_slot,
            'slotResidualZeros': residuals, 'freeScalarSymbols': names, 'commutationRolesJoined': bool(scalar_ok),
            'units': units, 'unitsJoined': bool(units_ok), 'candidateDenominator': denominator_route,
            'candidateDenominatorNonzero': domain, 'coefficientIdentityExtendsCancelledFactors': all(v is True for v in residuals),
            'safeSlotMap': bool(ok)})
        ratios.append(candidate if ok else None)
    summary = {'centreOutsideSlotZero': outside_zero, 'faceRelations': face_records,
               'completeCentreSlotImage': outside_zero is True, 'safeSlotMap': all(r is not None for r in ratios)}
    save_json(j.out/(label+'-relation.json'), summary)
    return ratios, summary


def assemble(j, label, ratios, parity_sum, parity_difference):
    total = j.op(label+'/sum-weight', 'add', ratios[0], ratios[1])
    negative = j.op(label+'/negative-minus-weight', 'mul', sp.S.NegativeOne, ratios[1])
    difference = j.op(label+'/difference-weight', 'add', ratios[0], negative)
    even = j.op(label+'/sum-product', 'mul', total, parity_sum)
    odd = j.op(label+'/difference-product', 'mul', difference, parity_difference)
    result = j.op(label+'/result', 'add', even, odd)
    return result


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args(); out = args.run_directory.resolve()
    out.relative_to(ROOT/'_scratch/s11c'); out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    specification_path = M/'S11c_d_centre_compatibility_run_inputs.json'
    specification_route = file_route(specification_path)
    specification = json.loads(specification_path.read_text())
    inp = Inputs(specification, out)
    save_json(out/'input-specification.json', specification)
    save_json(out/'input-specification-route.json', specification_route)
    for item in specification['files']:
        inp.pin(item['path'])
    export = inp.pin(specification['bExport'])
    with export.open() as f:
        prefix = []
        for line in f:
            if line.startswith('BUILD_INPUT_DIGESTS'): break
            prefix.append(line)
    inp.codec = source_codec(''.join(prefix))
    source = inp.pin(specification['c2Source']).read_text()
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'DIMENSION_SCHEMA' for t in n.targets))
    schema = ast.literal_eval(node.value)
    saved = inp.json(specification['navigation']['saved-operands'])
    reuse = inp.json(specification['navigation']['reuse-inputs'])
    joins = inp.json(specification['navigation']['join-inputs'])
    closure = inp.json(specification['navigation']['closure-inputs'])
    case_results = []; all_receipts = []
    for case_index, joined in enumerate(joins['rows']):
        case = joined['case']; name = '__'.join(case); directory = out/'cases'/name
        directory.mkdir(parents=True)
        centre_route = next(r['value'] for r in saved['records'] if r['case'] == case and r['value']['selector'][-1] == 'CENTER_FACE_GENERALIZED_ROW')
        face_route = next(r['normalizedFaceWork'] for r in reuse['bFaceInputs'] if r['case'] == case)
        parity = next(r for r in reuse['c2SavedClosures'] if r['tag'].endswith(name.replace('__','_')) and 'PARITY_BLOCKS' in r['tag'])
        face_returns = next(r for r in reuse['c2SavedClosures'] if r['tag'].endswith(name.replace('__','_')) and 'TERM_ORIGINS' in r['tag'])
        context = {'case': case, 'centreRow': centre_route, 'normalizedFaceRow': face_route,
                   'openRowJoin': joined, 'closedParity': parity, 'closedFaces': face_returns,
                   'closureProvenance': [r for r in closure['closureFaces'] if r['tag'].endswith(name.replace('__','_'))],
                   'governingInputPins': specification['files'], 'sourceMeaning':'Unrestricted row after closure/profile replacement/retained_shape/physical_fields, before extract.'}
        save_json(directory/'physical-input.json', context)
        j = Journal(directory, context)
        centre = inp.restore(centre_route); thickness = inp.restore(face_route)
        # Restore the exact row and units/normalization before source guards.
        original_row = inp.restore(joined['openEw']); normalization = inp.restore(joined['normalization'])
        save_value(directory/'saved-normalization-view.pickle', normalization)
        centre_dim = tuple(Fraction(int(v)) for v in inp.restore(joined['centreDimension']))
        open_dim = tuple(Fraction(int(v)) for v in inp.restore(joined['openEwDimension']))
        face_dims = [tuple(Fraction(int(v)) for v in inp.restore(f['dimension'])) for f in face_returns['fields']]
        save_json(directory/'native-row-units.json', {'centre':list(map(str,centre_dim)),
            'openEw':list(map(str,open_dim)), 'closedFaces':[list(map(str,d)) for d in face_dims]})
        if not all(d == open_dim for d in face_dims):
            raise ValueError('saved closed face units do not match their actual open row')
        face_dim = open_dim
        available = {a.name: a for a in centre.atoms(sp.Symbol)}
        slots = tuple(available[n] for n in specification['slotOrder'])
        oc, or_, o_ok = split_slots(j, 'open-row-join', original_row, slots)
        fc, fr_, f_ok = split_slots(j, 'face-row-join', thickness, slots)
        identities = [j.op('open-row-slot-identity/'+str(k), 'structural_equal', a, b) for k,(a,b) in enumerate(zip(oc,fc))]
        save_json(directory/'open-row-input-join.json', {'slotIdentities':identities,'openLinear':o_ok,'faceLinear':f_ok,'literalPreflight':joined['literalSlotTermsEqual']})
        if not (o_ok and f_ok and all(identities)):
            raise ValueError('actual consumed open row does not join normalized face load')
        ratios, relation = derive(j, 'baseline', centre, thickness, slots, centre_dim, face_dim, schema, specification['constantParameterRoles'])
        result = None
        if relation['safeSlotMap']:
            parity_values = {r['faceOrParity']: inp.restore(r['value']) for r in parity['fields']}
            result = assemble(j,'baseline-assembly',ratios,parity_values['FACE_SUM'],parity_values['FACE_DIFFERENCE'])
            output = save_value(directory/'centre-slot-image.pickle', result)
            zero = j.op('baseline-result-zero', 'zero_property', result)
        else:
            output = None; zero = None
        controls = []
        if case_index == 0:
            neg_p = j.op('controls/negative-pressure-slot','mul',sp.S.NegativeOne,slots[0])
            neg_j = j.op('controls/negative-jet-slot','mul',sp.S.NegativeOne,slots[1])
            mutations = [('one-pressure-sign',{slots[0]:neg_p}),
                         ('exchange-face-jets',{slots[1]:slots[3],slots[3]:slots[1]}),
                         ('whole-face-sign',{slots[0]:neg_p,slots[1]:neg_j})]
            for control, replacements in mutations:
                label='controls/'+control
                changed = j.op(label+'/imported-centre-row','replace',centre,replacements)
                equal_input = j.op(label+'/same-input','structural_equal',changed,centre)
                changed_ratios, cr = derive(j,label,changed,thickness,slots,centre_dim,face_dim,schema,specification['constantParameterRoles'])
                cr.update(name=control,importedInputUnchanged=equal_input)
                if cr['safeSlotMap'] and result is not None:
                    cv=assemble(j,label+'/assembly',changed_ratios,parity_values['FACE_SUM'],parity_values['FACE_DIFFERENCE'])
                    cr['result']=save_value(directory/(label+'-result.pickle'),cv)
                    cr['assemblyEqualsBaseline']=j.op(label+'/assembly-equality','structural_equal',cv,result)
                    cr['resultZeroProperty']=j.op(label+'/assembly-zero','zero_property',cv)
                controls.append(cr)
                save_json(directory/(label+'-summary.json'),cr)
        summary={'case':case,'relation':relation,'result':output,'zeroProperty':zero,'controls':controls,
                 'savedLiteralZeroParityRoutes':[f['value'] for f in parity['fields'] if inp.literal(f['value']).strip()=='Integer(0)'],
                 'closedLoadIndependentlyDerived':False,
                 'newCalls':sum(r['decision']=='new' for r in j.receipts),'exactPrecedingUses':sum(r['decision']!='new' for r in j.receipts),
                 'receiptCount':len(j.receipts),'fullCentreDynamicsConstructed':False,'oldClosureCalls':0}
        save_json(directory/'summary.json',summary)
        case_results.append(summary); all_receipts.extend(j.receipts)
        if case_index==0 and relation['safeSlotMap']:
            assembly_control=next(c for c in controls if c['name']=='whole-face-sign')
            if not (assembly_control['safeSlotMap'] and assembly_control.get('assemblyEqualsBaseline') is False):
                raise ValueError('whole-face routing control did not exercise final assembly')
        inp.restored.clear(); j.cache.clear(); sp.core.cache.clear_cache(); gc.collect()
    save_json(out/'operation-catalogue.json',all_receipts)
    count=inp.postcheck()
    spec_post = file_route(specification_path)
    save_json(out/'input-specification-posthash.json',spec_post)
    if spec_post != specification_route:
        raise ValueError('run input specification changed')
    # A finite manifest of newly emitted files, not a recursive ancestor copy.
    artifacts=[{'path':str(p),'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(out.rglob('*')) if p.is_file()]
    save_json(out/'artifacts.json',artifacts)
    checks={'status':'BOUNDED_CENTRE_DIAGNOSTIC_RECORDED','cases':case_results,'consumedFiles':count,
            'newCalls':sum(r['decision']=='new' for r in all_receipts),'exactPrecedingUses':sum(r['decision']!='new' for r in all_receipts),
            'workerSeconds':time.monotonic()-started,'oldClosureCalls':0,'newResponseCalls':0,'centreDynamicsConstructed':False,
            'implementationCommit':specification['implementationCommit'],'helperSha256':digest(Path(__file__)),
            'inputSpecification':specification_route,'artifacts':file_route(out/'artifacts.json')}
    save_json(out/'checks.json',checks)
    print((out/'checks.json').read_text(),end='')


if __name__=='__main__':
    main()
