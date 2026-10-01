#!/usr/bin/env python3
"""Extract native consumer source text with stdlib only; never restore science.

Constructor strings are data. This module parses their syntax and literal
addresses but does not import a producer, call _restore, or evaluate a symbolic
constructor. build_records() is deterministic for the pinned source bytes.
"""
import ast
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DESTINATION = BASE / '_measurements/S11c_upstream_mixed_consumer_native.json'
CASE = ('LAB_HELD', 'RHO4_CONSTANT')
PRESSURE_NAMES = {
    'delta_p_plus', 'delta_p_minus',
    'd_w_delta_p_plus', 'd_w_delta_p_minus',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def text_sha(value):
    return hashlib.sha256(value.encode()).hexdigest()


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def literal_record(path, key):
    """Read one original _restore string via its outer Python AST Constant."""
    found = False
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith('    ' + repr(key) + ':'):
                found, key_line = True, number
            elif found and line.startswith("        'value': "):
                node = ast.parse('{' + line.strip().rstrip(',') + '}', mode='eval').body
                require(isinstance(node, ast.Dict) and len(node.values) == 1,
                        'one native value mapping')
                call = node.values[0]
                require(isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                        and call.func.id == '_restore' and len(call.args) == 1
                        and not call.keywords, 'native restore-string schema')
                argument = call.args[0]
                require(isinstance(argument, ast.Constant) and isinstance(argument.value, str),
                        'native constructor is a literal string')
                text = argument.value
                return text, dict(
                    source=str(path), key=key, keyLine=key_line, valueLine=number,
                    sourceLineSha256=text_sha(line), sourceLineBytes=len(line.encode()),
                    constructorSha256=text_sha(text), constructorBytes=len(text.encode()),
                    extraction='AST Constant argument of native _restore call; no restoration')
    raise ValueError('missing native literal ' + key)


def tuple_arguments(value):
    """Navigate one serialized Tuple, preserving each exact argument substring."""
    require(value.startswith('Tuple(') and value.endswith(')'), 'serialized Tuple')
    start, depth, quote, escaped, result = 6, 0, None, False, []
    for index, char in enumerate(value[6:-1], 6):
        if quote:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
        elif char in "'\"":
            quote = char
        elif char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            require(depth >= 0, 'balanced constructor parentheses')
        elif char == ',' and depth == 0:
            result.append(value[start:index].strip())
            start = index + 1
    require(depth == 0 and quote is None and not escaped, 'complete constructor substring')
    result.append(value[start:-1].strip())
    return result


def literal_key(value):
    """Only Str/Integer/Tuple labels; no mathematical constructor evaluation."""
    node = ast.parse(value, mode='eval').body

    def visit(item):
        require(isinstance(item, ast.Call) and isinstance(item.func, ast.Name)
                and not item.keywords, 'literal key schema')
        if item.func.id == 'Tuple':
            return tuple(visit(arg) for arg in item.args)
        require(item.func.id in ('Str', 'Integer') and len(item.args) == 1,
                'only string/integer key atoms')
        atom = ast.literal_eval(item.args[0])
        require(type(atom) is (str if item.func.id == 'Str' else int), 'literal key type')
        return atom

    return visit(node)


def named(value, name):
    result = []
    for row in tuple_arguments(value):
        pair = tuple_arguments(row)
        require(len(pair) == 2, 'named pair arity')
        if literal_key(pair[0]) == name:
            result.append(pair[1])
    require(len(result) == 1, 'unique named source ' + name)
    return result[0]


def selected_case(value, key):
    result = []
    for index, row in enumerate(tuple_arguments(value)):
        pair = tuple_arguments(row)
        require(len(pair) == 2, 'case pair arity')
        if literal_key(pair[0]) == key:
            result.append((index, pair[1]))
    require(len(result) == 1, 'unique native case ' + repr(key))
    return result[0]


def case_record(text, provenance, key):
    index, selected = selected_case(text, key)
    return dict(
        provenance=provenance, case=list(key), caseIndex=index,
        caseConstructorText=selected, caseConstructorSha256=text_sha(selected),
        valueConstructorText=named(selected, 'VALUE'),
        scope='Complete original selected case including original grades/dimensions; not restored')


def pressure_occurrences(node):
    """Inventory names rather than infer that only the expected two slots exist."""
    result = []
    for item in ast.walk(node):
        if (isinstance(item, ast.Call) and isinstance(item.func, ast.Name)
                and item.func.id in ('Symbol', 'Function') and item.args
                and isinstance(item.args[0], ast.Constant)
                and isinstance(item.args[0].value, str)):
            name = item.args[0].value
            if 'delta_p' in name:
                result.append(dict(constructor=item.func.id, name=name))
    return result


def row_census(row_name, row_text):
    expanded = named(row_text, 'EXPANDED')
    entries = tuple_arguments(expanded) if row_name == 'U_BODY_BALANCE' else [expanded]
    records = []
    for component, expression in enumerate(entries):
        tree = ast.parse(expression, mode='eval').body
        occurrences = pressure_occurrences(tree)
        require(all(item['constructor'] == 'Symbol' and item['name'] in PRESSURE_NAMES
                    for item in occurrences), 'unrecognized pressure derivative/function: ' + repr(occurrences))
        is_add = (isinstance(tree, ast.Call) and isinstance(tree.func, ast.Name)
                  and tree.func.id == 'Add')
        children = tree.args if is_add else [tree]
        index, selected, covered = [], [], []
        for child_index, child in enumerate(children):
            # Native constructors are single-line ASCII. AST byte columns then
            # slice the exact original text without rescanning the full row for
            # each of its hundreds of children.
            require(expression.isascii() and child.lineno == child.end_lineno == 1,
                    'single-line ASCII native source span')
            text = expression[child.col_offset:child.end_col_offset]
            require(bool(text), 'exact native child source segment')
            hits = pressure_occurrences(child)
            covered.extend(hits)
            item = dict(childIndex=child_index, constructorSha256=text_sha(text),
                        constructorBytes=len(text.encode()), pressureSymbolOccurrences=hits)
            index.append(item)
            if hits:
                selected.append(dict(item, constructorText=text))
        order = lambda values: sorted((x['constructor'], x['name']) for x in values)
        require(order(covered) == order(occurrences), 'complete pressure/jet occurrence coverage')
        address = ['slab_operator', list(CASE), 'VALUE', row_name, 'EXPANDED']
        if row_name == 'U_BODY_BALANCE':
            address.append(component)
        records.append(dict(
            address=address, rowConstructorSha256=text_sha(row_text),
            expandedConstructorSha256=text_sha(expression), expandedConstructorBytes=len(expression.encode()),
            outerConstructor=getattr(getattr(tree, 'func', None), 'id', None),
            totalChildren=len(children), pressureNameQuery='substring delta_p in every Symbol/Function name',
            pressureSymbolOccurrenceCounts={name: sum(x['name'] == name for x in occurrences)
                                           for name in sorted(PRESSURE_NAMES)},
            completeTopLevelChildIndex=index, selectedPressureJetChildren=selected,
            allPressureOccurrencesCovered=True))
    return records


def build_records(base=BASE):
    base = Path(base).resolve()
    names = ('S11c_a_exports.py', 'S11c_b_exports.py', 'S11c_c1_exports.py',
             'S11c_c2_selfenergy_fold_sympy_audit.py', 'S11c_b_brane_operator_sympy_audit.py')
    paths = {name: base / 'scripts' / name for name in names}
    pins = {str(path): dict(sha256=file_sha(path), bytes=path.stat().st_size)
            for path in paths.values()}
    data = dict(
        status='AUTHOR_SOURCE_PREPARATION_NOT_SCIENTIFIC_RESULT_OR_REVIEW_CLEARANCE',
        scope='Selected upper-face eta*sigma_W consumer preparation; both-face source metadata preserved. No scientific evaluation.',
        sourcePins=pins,
        extractionMethod=dict(
            outer='AST Constant of the exact quoted _restore argument; never call _restore.',
            cases='Balanced constructor-text Tuple navigation; literal labels checked by AST schema.',
            rows='AST census of all pressure-named Symbol/Function nodes; exact source text for every matching Add child, hash/index of every child. No algebra.'))
    text, meta = literal_record(paths['S11c_b_exports.py'], 'mu_theta_operator')
    data['chemicalSource'] = case_record(text, meta, CASE)
    text, meta = literal_record(paths['S11c_c1_exports.py'], 's11c_c1_face_response')
    response_cases = named(text, 'CASES')
    data['faceResponseSources'] = dict(
        provenance=meta, outerCasePath=['CASES'],
        cases=[case_record(response_cases, meta, ('LAB_HELD', face, 'RHO4_CONSTANT'))
               for face in (1, -1)],
        note='Both complete native response cases; source equality between faces is not assumed.')
    geometry = {}
    for key in ('face_velocity', 'face_normal', 'face_shift', 'background_density_map'):
        a_text, a_meta = literal_record(paths['S11c_a_exports.py'], key)
        b_text, b_meta = literal_record(paths['S11c_b_exports.py'], key)
        keys = [('RHO4_CONSTANT',)] if key == 'background_density_map' else [
            ('LAB_HELD', face, 'DELTA_W', 'RHO4_CONSTANT') if key == 'face_shift'
            else ('LAB_HELD', face, 'DELTA_W') for face in (1, -1)]
        geometry[key] = dict(
            inheritedSource='S11c_b_exports.py consumed by native c2 Inputs',
            aLiteralSource=a_meta, bLiteralSource=b_meta,
            aAndBConstructorByteIdentity=a_text == b_text,
            cases=[case_record(b_text, b_meta, key_tuple) for key_tuple in keys])
        if a_text != b_text:
            geometry[key]['originalACases'] = [case_record(a_text, a_meta, key_tuple)
                                               for key_tuple in keys]
            geometry[key]['comparisonLimit'] = (
                'Constructor text differs; no algebraic equality or scientific discrepancy inferred.')
    data['geometry'] = geometry
    slab, meta = literal_record(paths['S11c_b_exports.py'], 'slab_operator')
    index, case = selected_case(slab, CASE)
    value = named(case, 'VALUE')
    data['slabConsumers'] = dict(
        provenance=meta, case=list(CASE), caseIndex=index,
        caseConstructorSha256=text_sha(case), caseConstructorBytes=len(case.encode()),
        caseGradesConstructorText=named(case, 'MULTIGRADE'),
        caseDimensionsConstructorText=named(case, 'DIMENSION_L_T_M'),
        rows=[record for name in ('U_BODY_BALANCE', 'THETA_BALANCE', 'E_W_BALANCE')
              for record in row_census(name, named(value, name))])
    literal_u = named(named(value, 'FACE_GENERALIZED_FORCE_ROWS'), 'U')
    data['physicalGeneralizedForceULiteral'] = dict(
        address=['slab_operator', list(CASE), 'VALUE', 'FACE_GENERALIZED_FORCE_ROWS', 'U'],
        constructorText=literal_u, constructorSha256=text_sha(literal_u))
    data['nativeSourceRouting'] = [
        dict(source=str(paths['S11c_c2_selfenergy_fold_sympy_audit.py']), lines=[558, 576],
             role='Live density binding; divide epsilon; identify c1 chemical/velocity with complete selected native amplitudes.'),
        dict(source=str(paths['S11c_c2_selfenergy_fold_sympy_audit.py']), lines=[478, 548],
             role='Physical pressure to reference value and lab-normal jet.'),
        dict(source=str(paths['S11c_c2_selfenergy_fold_sympy_audit.py']), lines=[609, 655],
             role='Expanded slab rows, pressure-slot substitution, retained rectangle, physical fields and weak extraction.'),
        dict(source=str(paths['S11c_b_brane_operator_sympy_audit.py']), lines=[2877, 2893],
             role='Mass closure already folded; do not append a second flux equation.'),
        dict(source=str(paths['S11c_b_brane_operator_sympy_audit.py']), lines=[3021, 3098],
             role='Mechanical stored-row orientation already applied; do not apply it again.')]
    data['limitations'] = [
        'Constructor text and syntactic coverage are not runtime-restoration or algebraic evidence.',
        'The complete chemical source still needs an actual transverse-trial evaluation.',
        'This census supplies no b contraction, modal excitation, both-face cancellation or leakage result.',
        'Lower-face source and consumer entries are unchanged provenance; new correction scope is upper face only.',
        'No producer, root/mode, forced-field, power, physical input or equation was computed or changed.']
    for path, pin in pins.items():
        require(file_sha(path) == pin['sha256'], 'source changed during census: ' + path)
    return data


def main():
    records = build_records()
    DESTINATION.write_text(json.dumps(records, indent=2) + '\n')
    print(json.dumps(dict(path=str(DESTINATION), bytes=DESTINATION.stat().st_size,
                          sha256=file_sha(DESTINATION),
                          pressureRows=[dict(address=r['address'], matchedChildren=len(r['selectedPressureJetChildren']),
                                             counts=r['pressureSymbolOccurrenceCounts'])
                                        for r in records['slabConsumers']['rows']])))


if __name__ == '__main__':
    main()
