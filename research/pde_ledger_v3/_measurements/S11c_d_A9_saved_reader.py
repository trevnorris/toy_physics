#!/usr/bin/env python3
"""Bounded navigation of selected saved A9 operands; no new physics operation.

Native SymPy/NumPy pickle codecs are retained. Producer modules cannot be
imported by this unpickler. The output is a navigation summary, not an export
or independent validation, and never treats omitted display content as zero.
"""
import argparse
import collections
import hashlib
import json
import os
import pickle
from pathlib import Path
import time

import numpy as np
import sympy as sp

M = Path(__file__).resolve().parent


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def route(path):
    links = []
    for part in (*reversed(path.parents), path):
        if part.is_symlink():
            links.append({'path': str(part), 'target': os.readlink(part)})
    resolved = path.resolve(strict=True)
    stat = resolved.stat()
    return {'logical': str(path), 'canonical': str(resolved), 'links': links,
            'bytes': stat.st_size, 'device': stat.st_dev, 'inode': stat.st_ino}


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError(('unapproved saved codec module', module, name))
        return super().find_class(module, name)


def label(value):
    if isinstance(value, str):
        return value
    if type(value) in (int, float, bool, type(None)):
        return repr(value)
    if isinstance(value, tuple):
        return '(' + ', '.join(label(v) for v in value) + ')'
    if isinstance(value, sp.Symbol):
        return value.name
    return '<' + type(value).__module__ + '.' + type(value).__name__ + '>'


def describe(value, depth=2):
    kind = type(value).__module__ + '.' + type(value).__name__
    if value is None or type(value) in (str, int, bool):
        return {'type': kind, 'value': value}
    if type(value) in (float, complex):
        return {'type': kind, 'literal': repr(value)}
    if isinstance(value, np.ndarray):
        return {'type': kind, 'shape': list(value.shape), 'dtype': str(value.dtype),
                'valuesNotSummarized': True}
    if isinstance(value, dict):
        entries = [{'key': label(k), 'keyType': type(k).__name__,
                    **({'value': describe(v, depth-1)} if depth > 0 else {})}
                   for k, v in value.items()]
        return {'type': kind, 'count': len(value), 'entries': entries}
    if isinstance(value, (tuple, list)):
        return {'type': kind, 'count': len(value),
                **({'items': [describe(v, depth-1) for v in value]} if depth > 0 else {})}
    if isinstance(value, (sp.Basic, sp.MatrixBase)):
        pending = [value]
        seen = set()
        heads = collections.Counter()
        symbols = set()
        while pending and len(seen) < 512:
            node = pending.pop()
            if id(node) in seen:
                continue
            seen.add(id(node))
            heads[type(node).__name__] += 1
            if isinstance(node, sp.Symbol):
                symbols.add(node.name)
            if isinstance(node, sp.MatrixBase) and not isinstance(node, sp.Basic):
                pending.extend(list(node))
            elif isinstance(node, sp.Basic):
                pending.extend(node.args)
        out = {'type': kind, 'visitedNodeIdentities': len(seen),
               'structureDisplayComplete': not pending,
               'observedHeads': dict(heads), 'observedSymbolNames': sorted(symbols)}
        if isinstance(value, sp.MatrixBase):
            out['shape'] = [value.rows, value.cols]
        # Only format small saved expressions. Formatting is not computation of
        # a residual or reconstruction of a derivative/integral/response.
        if not pending and len(seen) <= 48:
            out['savedExpressionText'] = sp.srepr(value)
        return out
    if isinstance(value, np.generic):
        return {'type': kind, 'shape': list(value.shape), 'dtype': str(value.dtype),
                'literal': repr(value)}
    return {'type': kind, 'displayUnavailable': True}


def selections():
    initial = json.loads((M/'S11c_d_A9_input_map.json').read_text())
    later = json.loads((M/'S11c_d_A9_saved_dependency_routes.json').read_text())
    rows = []
    for cp in later['checkpoints']:
        name = Path(cp['checkpoint']['path']).name
        for item in cp['routes']:
            key = item['manifestKey']
            fields = None
            if name == 'S11c_d_uniform_source_checkpoint.json' and key == 'uniform-source.pickle':
                fields = [['records', end, 'coupling'] for end in ('REFERENCE', 'LEFT', 'RIGHT')]
                fields += [['differences'], ['units'], ['profileEndpoints']]
            elif name == 'S11c_d_profile_form_checkpoint.json' and key == 'binding-checks.pickle':
                fields = [['moments']]
            elif name == 'S11c_d_remaining_case_coordinate_inputs_checkpoint.json' and key.endswith('/density-operands.pickle'):
                fields = [[]]
            elif name == 'S11c_d_remaining_case_coordinate_sources_checkpoint.json' and key.endswith('/coordinate-source.pickle'):
                fields = [['records'], ['densityAdvection']]
            elif name == 'S11c_d_reduced_action_source_checkpoint.json' and key == 'actions.pickle':
                fields = [[]]
            if fields is not None:
                rows.append({'path': item['logicalPath'], 'sha256': item['manifestSha256'],
                             'bytes': item['bytes'], 'fields': fields,
                             'checkpoint': cp['checkpoint']['path'], 'role': name})
    for cp in initial['checkpoints']:
        if Path(cp['checkpoint']).name != 'S11c_d_remaining_case_flux_checkpoint.json':
            continue
        for item in cp['selectedArtifactRoutes']:
            if item['manifestKey'].endswith('/reduced-action.pickle'):
                rows.append({'path': item['logicalPath'], 'sha256': item['manifestSha256'],
                             'bytes': item['manifestBytes'], 'fields': [['rows'], ['payloads'], ['reductionState']],
                             'checkpoint': cp['checkpoint'], 'role': 'own reduced action'})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(M.parents[2] / '_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    inputs = selections()
    source_paths = [Path(__file__).resolve(), M/'S11c_d_A9_input_map.json',
                    M/'S11c_d_A9_saved_dependency_routes.json',
                    M/'S11c_d_A9_saved_reader_plan.md',
                    M/'S11c_d_scattering_form_acceptance.json',
                    M.parent/'directives/S11c_d_SCATTERING_FORM_AMENDMENT.md']
    source_hashes = {str(p): digest(p) for p in source_paths}
    states = {v['path']: route(Path(v['path'])) for v in inputs}
    for row in inputs:
        if states[row['path']]['bytes'] != row['bytes']:
            raise ValueError(('saved input size changed', row['path']))
    save(base/'inputs.json', {'sourceFiles': source_hashes, 'inputs': inputs,
                             'routes': states, 'scope': 'Saved operand navigation only.'})
    groups = {}
    for row in inputs:
        groups.setdefault(states[row['path']]['canonical'], []).append(row)
    summaries = []
    for index, (canonical, consumers) in enumerate(groups.items()):
        path = Path(canonical)
        before = digest(path)
        if any(v['sha256'] != before for v in consumers):
            raise ValueError(('saved manifest hash mismatch', canonical))
        with path.open('rb') as stream:
            saved = SavedCodec(stream).load()
        selections_out = []
        for consumer in consumers:
            for fields in consumer['fields']:
                value = saved
                present = True
                for key in fields:
                    if not isinstance(value, dict) or key not in value:
                        present = False
                        break
                    value = value[key]
                selections_out.append({'consumer': consumer['path'], 'fieldPath': fields,
                                       'fieldPresent': present,
                                       **({'summary': describe(value, 1 if fields in (['records'], ['payloads'], ['reductionState']) else 3)}
                                          if present else {'disposition': 'Requested field absent; consult root keys. No value inferred.'})})
        output = base / ('saved-%02d.json' % index)
        save(output, {'canonical': canonical, 'sha256': before,
                      'consumers': consumers, 'root': describe(saved, 0),
                      'selections': selections_out})
        summaries.append({'path': str(output), 'sha256': digest(output)})
        del saved, selections_out
        if digest(path) != before:
            raise ValueError(('input changed after read', canonical))
    for logical, state in states.items():
        if route(Path(logical)) != state:
            raise ValueError(('consumed route changed', logical))
    for canonical, consumers in groups.items():
        if digest(Path(canonical)) != consumers[0]['sha256']:
            raise ValueError(('final consumed hash changed', canonical))
    if {str(p): digest(p) for p in source_paths} != source_hashes:
        raise ValueError('reader source changed')
    checks = {'status': 'COMPLETED_SAVED_A9_INPUT_NAVIGATION_NOT_PHYSICS_ACCEPTANCE',
              'sourceFiles': source_hashes, 'logicalInputs': len(inputs),
              'nativeCodecRestorations': len(groups), 'summaries': summaries,
              'scientificConstructionCalls': 0, 'independentPhysicsValidation': False,
              'allConsumedHashesAndRoutesUnchanged': True,
              'wallSeconds': time.monotonic()-started}
    rendered = json.dumps(checks, indent=2) + '\n'
    with (base/'checks.json').open('x') as stream:
        stream.write(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
