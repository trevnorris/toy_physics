#!/usr/bin/env python3
"""Recover independent transcript validation from saved scheduling-test operands."""
import argparse
import ast
import copy
import json
from pathlib import Path
import shutil
import time

import numpy as np
import S11c_d_parallel_momentum_preflight as focused
from S11c_d_output_codec import PayloadDecoder, decoded_lines

p = focused.p


def exact(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(exact(a[k], b[k]) for k in a)
    if isinstance(a, (tuple, list)):
        return len(a) == len(b) and all(exact(x, y) for x, y in zip(a, b))
    if isinstance(a, np.ndarray):
        return np.array_equal(a, b)
    return a == b


def source_join(previous, pins):
    changed = []
    name = '_measurements/S11c_d_parallel_momentum_preflight.py'
    for key, sha in pins.items():
        frozen = previous / 'source' / key
        p.require(p.digest(frozen) == sha, ('frozen source', key))
        if p.digest(p.ROOT / key) == sha:
            continue
        p.require(key == name, ('unexpected source change', key))
        old = ast.parse(frozen.read_text())
        new = ast.parse((p.ROOT / key).read_text())
        reset = ast.dump(ast.parse('p.engine.PAYLOAD_ENCODER=p.engine.PayloadEncoder()').body[0])
        removed = []

        class Strip(ast.NodeTransformer):
            def visit_Assign(self, node):
                if ast.dump(node) == reset:
                    removed.append(node)
                    return None
                return self.generic_visit(node)

        restored = Strip().visit(copy.deepcopy(new))
        p.require(len(removed) == 2 and ast.dump(restored) == ast.dump(old),
                  'whole preflight AST join after two encoder resets')
        changed.append({'path': key, 'before': sha, 'after': p.digest(p.ROOT / key),
                        'encoderResetsOnly': len(removed)})
    p.require(len(changed) == 1, 'repair census')
    return changed


def decoded(path, decoder=None):
    decoder = decoder or PayloadDecoder()
    entries = {}
    for line in path.read_text().splitlines():
        tag, sep, body = line.partition(': ')
        p.require(sep and tag not in entries, 'original transcript syntax/census')
        entries[tag] = p.native._restore(decoder.decode(body))
    return entries


def worker_packets(directory, preflight):
    manifest = json.loads((directory / 'workers.json').read_text())
    p.require(manifest['status'] == 'completed' and len(manifest['outcomes']) == 4,
              'four completed children')
    p.require(all(v['exitCode'] == 0 and v['stderrBytes'] == 0 for v in manifest['outcomes']),
              'child outcomes')
    packets = {}
    for key, sha in manifest['resultFiles'].items():
        task = tuple(map(int, key.split('-')))
        d = directory / ('worker-' + key)
        checks = json.loads((d / 'checks.json').read_text())
        p.require(p.digest(d / 'result.pickle') == sha == checks['resultSha256'], 'child packet hash')
        p.require((d / 'stderr').stat().st_size == 0 and (d / 'stdout').stat().st_size == 0,
                  'child logs')
        packet = p.unpickle(d / 'result.pickle')
        p.require(packet['task'] == task and packet['sourceFiles'] == preflight['sourceFiles']
                  and packet['provenance'] == preflight['provenance']
                  and packet['methodJoins'] == p.METHOD_JOINS, 'child source/method join')
        for record in packet['partialArtifacts']:
            p.require(p.digest(d / record['path']) == record['sha256'], 'child partial hash')
        packets[task] = packet
    p.require(set(packets) == {(ti, i) for ti in range(2) for i in (1, 2)}, 'task coverage')
    return packets, manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--previous', type=Path, required=True)
    args = parser.parse_args()
    base, previous = args.directory.resolve(), args.previous.resolve()
    base.relative_to(p.STORE); previous.relative_to(p.STORE)
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    # Preserve a complete pre-recovery identity inventory; no original file is modified.
    before = {str(q.relative_to(previous)): p.digest(q) for q in previous.rglob('*') if q.is_file()}
    p.save(base / 'original-inventory.json', before)
    old = json.loads((previous / 'parallel-preflight.json').read_text())
    repair = source_join(previous, old['sourceFiles'])
    data = p.prepare(base)
    r, bound, rows, variables, finest, provenance, pins = data
    recovery_path = Path(__file__).resolve()
    shutil.copyfile(recovery_path, base / 'source' / '_measurements' / recovery_path.name)
    recovery_sha = p.digest(recovery_path)
    p.require(provenance == old['provenance'] and p.METHOD_JOINS == old['methodJoins'],
              'unchanged physical/execution provenance')
    for name in p.PACKETS:
        p.require(p.digest(base / name) == before[name], 'accepted input packet join')
    launch = json.loads((previous.parent / 'launch.json').read_text())
    p.require(set(launch['threadEnvironment'].values()) == {'1'}, 'recorded native thread limits')
    prefix, prefix_manifest = worker_packets(previous / 'parallel-prefix', old)
    complete, complete_manifest = worker_packets(previous / 'parallel-full', old)
    prefix_residuals = {}
    for task, packet in prefix.items():
        serial = p.unpickle(previous / f'serial-prefix-{task[0]}-{task[1]}' / 'prefix-256.pickle')
        prefix_residuals[str(task)] = focused.same_group(serial, packet['group'])
        p.require(serial['setting'] == packet['setting'] and serial['nodeCount'] == 65536,
                  'native prefix domain/count')
        # Replay only the literal frequency prefix on saved node positions; no integrals.
        tracer = p.ResumableMomentum(bound['rows'], bound['sources'], r)
        count = 0
        for index, (points, weights) in enumerate(tracer.batches(variables, serial['setting'],
                                                                bound['pairs'], serial['width']), 1):
            env = {v: points[:, i] for i, v in enumerate(variables)}
            env[r.regulator] = np.full(len(weights), serial['setting']['regulator'])
            for row in rows:
                for factor in row['factors']:
                    tracer.trace_frequency(bound['sources'][(task[0], factor['sourceIndex'])],
                                           env, serial['setting'], {})
            count += len(weights)
            if index == serial['batchCount']:
                break
        p.require(count == serial['nodeCount'] and
                  tracer.source_frequency_census == packet['telemetry']['sourceFrequencyCensus'],
                  'literal prefix frequency census')

    results = {}
    for kind in ('serial-full', 'parallel-full'):
        original = previous / kind
        packet = p.unpickle(original / 'three-momentum-source.pickle')
        p.require(packet['boundPacketSha256'] == provenance['BOUND_PACKET_SHA256'], 'full bound join')
        result = packet['result']; results[kind] = result
        for record in (*result['recordArtifacts'], *result['partialArtifacts']):
            p.require(p.digest(original / record['path']) == record['sha256'], 'full saved record hash')
        serial_peak = 0
        for item in result['records']:
            rebuilt = p.native.full_action(r, bound, finest[item['test']], item['group'])
            p.require(all(exact(item[k], rebuilt[k]) for k in rebuilt), 'saved action/native terms join')
            if item['index']:
                worker_group = complete[(item['test'], item['index'])]['group']
                serial_peak = max(serial_peak, worker_group['peakWorkspaceEstimateBytes'])
                expected_group = dict(worker_group)
                if kind == 'serial-full':
                    expected_group['peakWorkspaceEstimateBytes'] = serial_peak
                p.require(exact(item['group'], expected_group), 'full saved worker group/peak join')
            else:
                ref = next(g for g in finest[item['test']]['groups'] if g['variables'] == variables)
                p.require(exact(item['group'], ref), 'retained baseline group')
        d = base / kind; d.mkdir()
        shutil.copyfile(original / 'three-momentum-source.pickle', d / 'three-momentum-source.pickle')

    left_result, right_result = results['serial-full'], results['parallel-full']
    telemetry_differences = []
    for x, y in zip(left_result['records'], right_result['records'], strict=True):
        focused.same_group(x['group'], y['group'])
        lcopy, rcopy = copy.deepcopy(x), copy.deepcopy(y)
        lp = lcopy['group'].pop('peakWorkspaceEstimateBytes')
        rp = rcopy['group'].pop('peakWorkspaceEstimateBytes')
        if lp != rp:
            telemetry_differences.append({'test': x['test'], 'index': x['index'],
                                         'serialCumulativePeakBytes': lp, 'workerPeakBytes': rp})
        p.require(exact(lcopy, rcopy), 'complete serial/parallel physical record equality')
    for key in left_result:
        if key not in ('records', 'recordArtifacts'):
            p.require(exact(left_result[key], right_result[key]), ('complete result census', key))
    p.require([v['path'] for v in left_result['recordArtifacts']] ==
              [v['path'] for v in right_result['recordArtifacts']], 'separately verified record census')
    # The failed standalone transcript must decode only with the preceding stream.
    rejected = False
    try:
        list(decoded_lines(previous / 'parallel-full' / 'full.out'))
    except ValueError as exc:
        p.require(str(exc) == 'undefined shared-payload reference', 'unexpected codec failure')
        rejected = True
    p.require(rejected, 'missing isolated-stream negative control')
    shared_decoder = PayloadDecoder()
    original_left = decoded(previous / 'serial-full' / 'full.out', shared_decoder)
    original_right = decoded(previous / 'parallel-full' / 'full.out', shared_decoder)
    p.require(original_left == original_right, 'original concatenated payload identity')
    entries = {}; keysets = {}; paths = {}
    for kind in ('serial-full', 'parallel-full'):
        # Fresh encoder per standalone transcript; constructor/emitter/decoder unchanged.
        p.engine.PAYLOAD_ENCODER = p.engine.PayloadEncoder()
        p.engine.EMISSION_LINES.clear()
        entries[kind], keysets[kind], paths[kind] = p.native.emit_and_replay(
            base / kind, results[kind], r, bound, rows, finest, provenance)
    differences = {kind: [k for k in set(entries[kind]) | set(original_left)
                         if entries[kind].get(k) != original_left.get(k)] for kind in entries}
    p.save(base / 'emission-differences.json', differences)
    p.require(not any(differences.values()) and entries['serial-full'] == entries['parallel-full']
              and keysets['serial-full'] == keysets['parallel-full']
              and paths['serial-full'] == paths['parallel-full'], 'standalone emission/metadata identity')
    # The failed run did not persist its monotonic serial timer. Report only the
    # recorded file-time interval, explicitly distinguished from an exact speedup.
    first = previous / 'serial-prefix-0-1' / 'prefix-64.pickle'
    last = previous / 'serial-prefix-1-2' / 'prefix-256.pickle'
    interval = (last.stat().st_mtime_ns - first.stat().st_mtime_ns) / 1e9
    phase_interval = (last.stat().st_mtime_ns - (previous / 'parallel-preflight.json').stat().st_mtime_ns) / 1e9
    timing = {'partialSerialFileTimeIntervalSeconds': interval,
              'serialPhaseFileTimeIntervalSeconds': phase_interval,
              'parallelMonotonicWallSeconds': prefix_manifest['wallSeconds'],
              'partialSerialToParallelRatio': interval / prefix_manifest['wallSeconds'],
              'phaseFileTimeRatio': phase_interval / prefix_manifest['wallSeconds'],
              'scope': 'Serial monotonic timer was not saved before decoder failure. File-time intervals reconstruct approximate speedup only; the partial interval omits the first 64 serial batches.'}
    result = {'sourceFiles': pins, 'originalSourceFiles': old['sourceFiles'], 'provenance': provenance,
              'methodJoins': p.METHOD_JOINS, 'repair': repair, 'recoverySourceSha256': recovery_sha,
              'prefixResiduals': prefix_residuals, 'frequencyCensusExact': True,
              'fullNumericalResultsExactlyEqual': True, 'fullRecordCount': len(left_result['records']),
              'workspaceTelemetryDifferences': telemetry_differences,
              'workspaceScope': 'Each serial peak equals the running maximum of independent worker peaks; every other record operand agrees exactly. Pickle hashes are checked per stream, not required to match across object-sharing layouts.',
              'fullActionResiduals': [float(np.max(abs(x['action'] - y['action'])))
                                     for x, y in zip(left_result['records'], right_result['records'], strict=True)],
              'emissionDifferences': differences, 'tagCount': len(entries['serial-full']),
              'writeKeyCount': len(keysets['serial-full']), 'metadataPaths': paths['serial-full'],
              'originalIsolatedStreamRejected': rejected, 'originalCombinedStreamExact': True,
              'prefixWorkerManifest': prefix_manifest, 'fullWorkerManifest': complete_manifest,
              'workerPeakRssKiB': {str(t): v['peakRssKiB'] for t, v in prefix.items()},
              'threadEnvironment': launch['threadEnvironment'], 'timing': timing,
              'originalInventorySha256': p.digest(base / 'original-inventory.json'),
              'transcripts': {kind: {'sha256': p.digest(base / kind / 'full.out'),
                                    'bytes': (base / kind / 'full.out').stat().st_size} for kind in entries},
              'packetHashes': {kind: p.digest(base / kind / 'three-momentum-source.pickle') for kind in entries},
              'wallSeconds': time.monotonic() - started,
              'scope': 'Saved scheduling-test quadratures only. No numerical integration rerun, convergence, physical tails, Abel limit, scattering or poles.'}
    p.require(before == {n: p.digest(previous / n) for n in before}, 'original artifacts changed')
    p.require(pins == {n: p.digest(p.ROOT / n) for n in pins}
              and p.digest(recovery_path) == recovery_sha, 'current recovery sources changed')
    p.require(all(result['packetHashes'][kind] == before[kind + '/three-momentum-source.pickle']
                  for kind in entries) and not p.engine.PHYSICAL_METADATA.dimensions.constraints,
              'post-emission packet/dimension guard')
    p.save(base / 'checks.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
