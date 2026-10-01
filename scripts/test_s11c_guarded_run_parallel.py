"""Pooled guard tooling tests: no scientific imports, payloads or jobs."""
import copy
import fcntl
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import s11c_guarded_run as guard

GIB = 1024**3


def record(memory=4, cpu=7, pid=None):
    return {'pool': 'test', 'poolMemoryMax': guard.POOL_MEMORY_MAX,
            'memoryMax': memory * GIB, 'tasksMax': 32, 'cpu': cpu,
            'owner': guard.process_identity(os.getpid() if pid is None else pid),
            'logDirectory': '/synthetic'}


def live_record(saved, used=1):
    return {'memory.max': str(saved['memoryMax']), 'memory.swap.max': '0',
            'pids.max': str(saved['tasksMax']), 'memory.current': str(used * GIB),
            'systemd.CPUAffinity': str(saved['cpu'])}


def racing_reservation(directory, number, memory, ready, release, result):
    """Each subprocess invokes the real serialized admission implementation."""
    spec = {**guard.resource_limits(memory, 32, 'test'), 'pool': 'test',
            'unit': 's11c-guard-test-' + str(number), 'logDirectory': directory}
    ready.wait()
    try:
        with patch.object(guard, 'active_units', return_value=[]), \
             patch.object(guard, 'pool_cpus', return_value=[7, 6, 5]), \
             patch.object(guard, 'available_memory', return_value=32 * GIB):
            guard.reserve_pool(spec, Path(directory))
        result.put(('admitted', spec['cpu'], spec['memoryMax']))
        release.wait()
        with patch.object(guard, 'active_units', return_value=[]):
            guard.release_pool(spec, Path(directory))
    except RuntimeError as error:
        result.put(('refused', str(error)))


class ParallelAdmission(unittest.TestCase):
    def test_budget_and_outstanding_host_headroom(self):
        saved = record(8)
        reservations = {'s11c-guard-test.service': saved}
        live = {unit: live_record(value, used=3) for unit, value in reservations.items()}
        plan = guard.admission_plan(reservations, live, 8 * GIB, [7, 6], 17 * GIB)
        self.assertEqual(plan['cpu'], 6)
        self.assertEqual(plan['outstandingHeadroomBytes'], 5 * GIB)
        with self.assertRaisesRegex(RuntimeError, 'aggregate'):
            guard.admission_plan(reservations, live, 9 * GIB, [7, 6], 64 * GIB)
        with self.assertRaisesRegex(RuntimeError, 'host reserve'):
            guard.admission_plan(reservations, live, 8 * GIB, [7, 6], 17 * GIB - 1)

    def test_pending_and_orphan_services_still_count(self):
        saved = record(10)
        records = {'s11c-guard-test.service': saved}
        # A launch reserved before its service appears still reserves all memory.
        with self.assertRaisesRegex(RuntimeError, 'aggregate'):
            guard.admission_plan(records, {}, 7 * GIB, [7, 6], 64 * GIB)
        # A launcher gone but its verified service alive is still accountable.
        result = guard.admission_plan(records, {next(iter(records)): live_record(saved)},
                                      6 * GIB, [7, 6], 64 * GIB, identity=lambda _: None)
        self.assertEqual(result['cpu'], 6)
        # A dead launcher with no service must never be silently reclaimed.
        with self.assertRaisesRegex(RuntimeError, 'stale reservation'):
            guard.admission_plan(records, {}, GIB, [7, 6], 64 * GIB, identity=lambda _: None)

    def test_unknown_or_mismatched_service_refused(self):
        saved = record()
        with self.assertRaisesRegex(RuntimeError, 'unaccounted'):
            guard.admission_plan({}, {'orphan.service': live_record(saved)}, GIB, [7], 64 * GIB)
        for key, value in [('memory.max', 'max'), ('memory.swap.max', '1'),
                           ('pids.max', '64'), ('systemd.CPUAffinity', '6-7')]:
            actual = live_record(saved)
            actual[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(RuntimeError, 'does not match'):
                guard.admission_plan({'known': saved}, {'known': actual}, GIB, [7, 6], 64 * GIB)

    def test_cpu_and_registry_integrity_refused(self):
        saved = record()
        with self.assertRaisesRegex(RuntimeError, 'physical CPU'):
            guard.admission_plan({'known': saved}, {}, GIB, [7], 64 * GIB)
        duplicate = copy.deepcopy(saved)
        with self.assertRaisesRegex(RuntimeError, 'overlapping'):
            guard.admission_plan({'a': saved, 'b': duplicate}, {}, GIB, [7, 6], 64 * GIB)
        for key, value in [('memoryMax', 0), ('memoryMax', 17 * GIB),
                           ('poolMemoryMax', 32 * GIB), ('tasksMax', 128)]:
            changed = copy.deepcopy(saved)
            changed[key] = value
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                guard.admission_plan({'known': changed}, {}, GIB, [7, 6], 64 * GIB)

    def test_release_requires_terminal_service_and_owner(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            saved = record()
            spec = {**saved, 'unit': 's11c-guard-test'}
            guard.write_reservations(root, {spec['unit'] + '.service': saved})
            with patch.object(guard, 'active_units', return_value=[spec['unit'] + '.service']):
                with self.assertRaisesRegex(RuntimeError, 'still active'):
                    guard.release_pool(spec, root)
            with patch.object(guard, 'active_units', return_value=[]):
                changed = copy.deepcopy(spec)
                changed['owner']['startTicks'] = 'wrong'
                with self.assertRaisesRegex(RuntimeError, 'ownership'):
                    guard.release_pool(changed, root)
                guard.release_pool(spec, root)
            self.assertEqual(guard.read_reservations(root), {})

    def test_shared_pool_and_exclusive_legacy_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'guard.lock'
            shared = guard.acquire_execution_lock(path, pooled=True)
            shared2 = guard.acquire_execution_lock(path, pooled=True)
            with self.assertRaises(BlockingIOError):
                guard.acquire_execution_lock(path)
            shared.close()
            shared2.close()
            exclusive = guard.acquire_execution_lock(path)
            with self.assertRaises(BlockingIOError):
                guard.acquire_execution_lock(path, pooled=True)
            exclusive.close()

    def test_simultaneous_admissions_serialize_without_overbooking(self):
        context = multiprocessing.get_context('fork')
        with tempfile.TemporaryDirectory() as directory:
            ready, release, result = context.Event(), context.Event(), context.Queue()
            children = [context.Process(target=racing_reservation,
                args=(directory, index, 6, ready, release, result)) for index in range(3)]
            try:
                for child in children:
                    child.start()
                ready.set()
                results = [result.get(timeout=10) for _ in children]
                accepted = [item for item in results if item[0] == 'admitted']
                self.assertEqual(len(accepted), 2, results)
                self.assertEqual(len({item[1] for item in accepted}), 2)
                self.assertEqual(sum(item[2] for item in accepted), 12 * GIB)
                self.assertEqual(len(guard.read_reservations(Path(directory))), 2)
                release.set()
                for child in children:
                    child.join(timeout=10)
                    self.assertEqual(child.exitcode, 0)
                self.assertEqual(guard.read_reservations(Path(directory)), {})
            finally:
                release.set()
                for child in children:
                    if child.is_alive():
                        child.terminate()
                    child.join(timeout=10)


class ParallelEnforcement(unittest.TestCase):
    def profile(self):
        spec = {**guard.resource_limits(8, 32, 'test'), 'pool': 'test',
                'poolMemoryMax': guard.POOL_MEMORY_MAX, 'cpu': 7,
                'priorityPolicy': 'desktop-managed',
                'affinityEnforcement': 'systemd-and-inherited-process'}
        actual = {'memory.max': str(8 * GIB), 'memory.swap.max': '0',
                  'pids.max': '32', 'nice': 6, 'affinity': [7],
                  'threads': dict.fromkeys(guard.THREADS, '1'),
                  'nativeAddressSpaceBytes': [8 * GIB] * 2, 'systemd.CPUAffinity': '7'}
        return spec, actual

    def test_live_service_affinity_without_delegated_cpuset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            group = root / 'test-service'
            group.mkdir()
            for name, value in {'memory.max': str(8 * GIB), 'memory.swap.max': '0',
                                'memory.current': str(GIB), 'pids.max': '32'}.items():
                (group / name).write_text(value)
            self.assertFalse((group / 'cpuset.cpus.effective').exists())
            def resolve_path(value):
                return root if value == '/sys/fs/cgroup' else Path(value)
            with patch.object(guard, 'Path', side_effect=resolve_path), \
                 patch.object(guard.subprocess, 'check_output', return_value=
                    'ControlGroup=/test-service\nRuntimeMaxUSec=infinity\nRestart=no\nCPUAffinity=7\n'):
                actual = guard.service_resources('s11c-guard-test.service')
            self.assertEqual(actual['systemd.CPUAffinity'], '7')
            self.assertEqual(actual['memory.max'], str(8 * GIB))
            self.assertNotIn('cpuset.cpus.effective', actual)

    def test_native_cgroup_cpu_threads_and_desktop_priority(self):
        spec, actual = self.profile()
        self.assertTrue(guard.limits_match(actual, spec))
        for key, value in [('nativeAddressSpaceBytes', [-1, -1]),
                           ('memory.max', 'max'), ('memory.swap.max', '1'),
                           ('systemd.CPUAffinity', '6-7'), ('affinity', [6]),
                           ('pids.max', 'max'), ('threads', dict.fromkeys(guard.THREADS, '2'))]:
            changed = copy.deepcopy(actual)
            changed[key] = value
            with self.subTest(key=key):
                self.assertFalse(guard.limits_match(changed, spec))

    def test_no_elapsed_deadline_and_memory_protection(self):
        for low_memory in (False, True):
            with self.subTest(low_memory=low_memory), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                group = root / 'cgroup'
                group.mkdir()
                spec, actual = self.profile()
                for name in ('memory.max', 'memory.swap.max', 'pids.max'):
                    (group / name).write_text(actual[name])
                spec.update(unit='synthetic', command=['synthetic'], cwd=directory)
                self.assertFalse((group / 'cpuset.cpus.effective').exists())
                manifest = root / 'invocation.json'
                manifest.write_text(json.dumps(spec))
                child = MagicMock()
                child.pid = 123
                child.wait.side_effect = [subprocess.TimeoutExpired('synthetic', 2),
                                         subprocess.TimeoutExpired('synthetic', 2), 0]
                with patch.object(guard, 'group_path', return_value=group), \
                     patch.object(guard.os, 'sched_setaffinity'), \
                     patch.object(guard.os, 'sched_getaffinity', return_value={7}), \
                     patch.object(guard.os, 'getpriority', return_value=6), \
                     patch.object(guard.resource, 'setrlimit') as native, \
                     patch.object(guard.resource, 'getrlimit', return_value=(8 * GIB, 8 * GIB)), \
                     patch.dict(os.environ, dict.fromkeys(guard.THREADS, '1')), \
                     patch.object(guard, 'available_memory', return_value=24 * GIB), \
                     patch.object(guard, 'sample', return_value={'hostAvailableBytes': (2 if low_memory else 24) * GIB}), \
                     patch.object(guard.subprocess, 'check_output', return_value='RuntimeMaxUSec=infinity\nRestart=no\nCPUAffinity=7\n'), \
                     patch.object(guard.subprocess, 'Popen', return_value=child) as launch, \
                     patch.object(guard, 'stop'), \
                     patch.object(guard.time, 'monotonic', side_effect=[0, 100000000]):
                    code = guard.child_main(manifest)
                native.assert_called_once_with(guard.resource.RLIMIT_AS, (8 * GIB, 8 * GIB))
                self.assertEqual(launch.call_args.kwargs['env']['S11C_POOLED_GUARD_MANIFEST'], str(manifest))
                outcome = json.loads((root / 'child-outcome.json').read_text())
                self.assertEqual(code, 124 if low_memory else 0)
                self.assertEqual(outcome['wallSeconds'], 100000000)
                self.assertEqual(bool(outcome['guardReason']), low_memory)


class SupervisorSharedPrerequisite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1] / 'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'
        spec = importlib.util.spec_from_file_location('supervisor_tooling', path)
        cls.supervisor = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.supervisor)

    def test_default_exclusive_and_missing_pooled_evidence_fatal(self):
        self.assertEqual(self.supervisor.prerequisite_lock_mode(False), fcntl.LOCK_EX)
        with patch.dict(os.environ, {}, clear=True), self.assertRaises(KeyError):
            self.supervisor.prerequisite_lock_mode(True)

    def test_verified_shared_read_and_wrong_context_refusal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / 'invocation.json'
            spec = {'pool': 'test', 'priorityPolicy': 'desktop-managed', 'cpu': 7, 'memoryMax': 8 * GIB}
            manifest.write_text(json.dumps(spec))
            actual_group = next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines()
                                if line.startswith('0::'))
            validation = {'verified': True, 'actual': {'group': str(Path('/sys/fs/cgroup') / actual_group.lstrip('/'))}}
            validation_path = root / 'limit-validation.json'
            validation_path.write_text(json.dumps(validation))
            with patch.object(self.supervisor, 'STORE', root), \
                 patch.dict(os.environ, {'S11C_POOLED_GUARD_MANIFEST': str(manifest)}), \
                 patch.object(self.supervisor.os, 'sched_getaffinity', return_value={7}), \
                 patch.object(self.supervisor.resource, 'getrlimit', return_value=(8 * GIB, 8 * GIB)):
                self.assertEqual(self.supervisor.prerequisite_lock_mode(True), fcntl.LOCK_SH)
                validation['actual']['group'] = '/wrong/cgroup'
                validation_path.write_text(json.dumps(validation))
                with self.assertRaisesRegex(ValueError, 'verified pooled containment'):
                    self.supervisor.prerequisite_lock_mode(True)


if __name__ == '__main__':
    unittest.main()
