#!/usr/bin/env python3
"""Fresh unfinished metadata continuation; preserve the failed reader."""
import ast
import copy
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_remainder_inputs as old

ORIGINAL_SHA = '59983753e436d655f7ff48a6e7860462d06442a696aecdd7d3cf12741438765e'
FAILED = old.F/'remainder-inputs-recovery-01/failed-file-inventory.json'
FAILED_SHA = '258e89393014c983ab7f395c98dc50f1daefd1cc11a000f629b6f2f9075c3a24'


def printer_metadata(reader, printer_file, printer):
    namespace = Path(old.inspect.getsourcefile(old.sp.lambdify))
    module = ast.parse(namespace.read_text())
    assignment = next(n for n in module.body if isinstance(n, ast.Assign) and
                      any(isinstance(t, ast.Name) and t.id == 'NUMPY_TRANSLATIONS' for t in n.targets))
    translation = ast.literal_eval(assignment.value)['Heaviside']
    old.require(translation == 'heaviside', 'actual native numpy namespace translation')
    method = printer._print_Heaviside
    return {'file': reader.retain(printer_file), 'dedicatedMethodPresent': callable(method),
            'body': old.inspect.getsource(method) if callable(method) else None,
            'namespaceFile': reader.retain(namespace), 'namespaceAssignment': ast.unparse(assignment),
            'actualTranslation': translation, 'lambdifyInvoked': False}


def main():
    original = Path(old.__file__).resolve()
    old.require(old.saved.digest(original) == ORIGINAL_SHA, 'immutable failed reader')
    parsed = ast.parse(original.read_text()); fn = next(n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    before = copy.deepcopy(fn); matches = []
    for node in ast.walk(fn):
        if isinstance(node, ast.Dict):
            for i, key in enumerate(node.keys):
                if isinstance(key, ast.Constant) and key.value == 'nativeHeavisidePrinter':
                    matches.append((node, i, node.values[i])); node.values[i] = ast.parse('printer_metadata(reader, printer_file, NumPyPrinter)', mode='eval').body
    old.require(len(matches) == 1, 'exact single source-metadata repair')
    repaired = copy.deepcopy(fn); matches[0][0].values[matches[0][1]] = matches[0][2]
    old.require(ast.dump(fn) == ast.dump(before), 'whole main reverse AST except single metadata field')
    final_index = next(i for i, n in enumerate(repaired.body) if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call) and ast.unparse(n.value.func) == 'reader.postcheck')
    extra = ast.parse("""
reader.retain(original, ORIGINAL_SHA)
failed = reader.json(FAILED, FAILED_SHA)
for path, record in failed.items(): reader.retain(path, record['sha256'])
journal.json('reader-metadata-repair.json', {'original': reader.retain(original), 'failedInventory': reader.retain(FAILED), 'wholeMainReverseAST': True, 'newScientificCalls': 0, 'onlyChangedNativeSourceMetadata': True})
""").body
    repaired.body[final_index:final_index] = extra
    module = ast.fix_missing_locations(ast.Module(body=[repaired], type_ignores=[]))
    env = dict(vars(old), printer_metadata=printer_metadata, original=original, ORIGINAL_SHA=ORIGINAL_SHA, FAILED=FAILED, FAILED_SHA=FAILED_SHA,
               __file__=str(Path(__file__).resolve()), NAME='S11c_d_remaining_case_frequency_remainder_inputs_recover')
    exec(compile(module, '<saved-reader-metadata-continuation>', 'exec'), env); env['main']()


if __name__ == '__main__': main()
