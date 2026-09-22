#!/usr/bin/env python3
"""Finish only saved row metadata; retain all three no-science failures."""
import ast
import copy
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_remainder_inputs as original

NAME = 'S11c_d_remaining_case_frequency_remainder_inputs_finish'
SOURCE_SHA = '59983753e436d655f7ff48a6e7860462d06442a696aecdd7d3cf12741438765e'
FAILED = original.F/'remainder-inputs-recovery-03/failed-file-inventory.json'
FAILED_SHA = '1cbb9e191d0e166e904e115513a979b6d43bfd9739d92c6c70d8b34fc043a004'
CALLER = original.F/'remainder-inputs-recovery-02/complete/native-caller-and-unsaved-internal-history.json'
CALLER_SHA = '848cf9501bee0ff0c94d90cc2b4e06f975666a19a28fc915529eda575bdf7e70'


class MetadataJournal(original.saved.Journal):
    def json(self, name, value):
        def scalar(v):
            if isinstance(v, original.sp.Integer): return int(v)
            raise TypeError(('unsupported saved metadata type', type(v).__name__))
        # Serialize before opening a new file. This changes metadata only.
        native_json = json.loads(json.dumps(value, default=scalar))
        return super().json(name, native_json)


def reuse_caller(reader, journal):
    record = reader.retain(CALLER, CALLER_SHA)
    metadata = reader.json(CALLER, CALLER_SHA)
    for entry in metadata['native'].values(): reader.retain(entry['file']['logical'], entry['file']['sha256'])
    for rec in [metadata['complexCaller'], metadata['nativeHeavisidePrinter']['file'], metadata['nativeHeavisidePrinter']['namespaceFile']]: reader.retain(rec['logical'], rec['sha256'])
    target = journal.base/'native-caller-and-unsaved-internal-history.json'
    original.require(not target.exists() and not target.is_symlink(), 'fresh exact metadata reference')
    target.symlink_to(Path(record['canonical']))
    journal.artifacts[target.name] = {'path': str(target), 'sha256': CALLER_SHA, 'bytes': record['bytes']}
    reader.retain(target, CALLER_SHA)


def finish_provenance(reader, journal):
    for path, record in reader.json(FAILED, FAILED_SHA).items(): reader.retain(path, record['sha256'])
    journal.json('metadata-finish-joins.json', {'failedInventory': reader.retain(FAILED), 'completedCallerReused': reader.retain(CALLER),
        'original': reader.retain(original.__file__, SOURCE_SHA), 'wholeOriginalMainReverseAST': True,
        'JSONIntegerBoundaryOnly': True, 'nativePhysicalTypesUnchanged': True, 'newScientificCalls': 0})


def main():
    original.require(original.saved.digest(original.__file__) == SOURCE_SHA, 'immutable original reader')
    parsed = ast.parse(Path(original.__file__).read_text()); fn = next(n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    before = copy.deepcopy(fn)
    start = next(i for i,n in enumerate(fn.body) if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'native')
    end = next(i for i,n in enumerate(fn.body) if isinstance(n, ast.FunctionDef) and n.name == 'forbidden')
    removed = copy.deepcopy(fn.body[start:end]); fn.body[start:end] = ast.parse('reuse_caller(reader, journal)').body
    changed = copy.deepcopy(fn); changed.body[start:start+1] = removed
    original.require(ast.dump(changed) == ast.dump(before), 'exact completed source-metadata prefix replacement')
    for node in ast.walk(fn):
        if isinstance(node, ast.Call) and ast.unparse(node.func) == 'saved.Journal': node.func = ast.Name(id='MetadataJournal',ctx=ast.Load())
    at = next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func) == 'reader.postcheck')
    fn.body[at:at] = ast.parse('finish_provenance(reader, journal)').body
    env = dict(vars(original), __file__=str(Path(__file__).resolve()), NAME=NAME, MetadataJournal=MetadataJournal, reuse_caller=reuse_caller, finish_provenance=finish_provenance)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])), '<saved-row-metadata-finish>', 'exec'), env); env['main']()


if __name__ == '__main__': main()
