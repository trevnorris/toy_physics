"""Resource-profile and fail-closed checks; no scientific calculation."""
import copy
import unittest

import s11c_guarded_run as guard


class GuardLimits(unittest.TestCase):
    def profile(self, memory=2, tasks=32):
        spec = {**guard.resource_limits(memory, tasks), 'cpu': 3}
        actual = {'memory.max': str(memory * 1024**3), 'memory.swap.max': '0',
                  'pids.max': str(tasks), 'nice': 15, 'affinity': [3],
                  'threads': {name: '1' for name in guard.THREADS}}
        return spec, actual

    def test_existing_default_unchanged(self):
        self.assertEqual(guard.resource_limits(),
                         {'memoryMax': 2147483648, 'swapMax': 0, 'tasksMax': 32})
        spec, actual = self.profile()
        self.assertTrue(guard.limits_match(actual, spec))

    def test_authorized_profiles(self):
        for memory, tasks in [(2, 64), (4, 64), (4, 32)]:
            with self.subTest(memory=memory, tasks=tasks):
                spec, actual = self.profile(memory, tasks)
                self.assertTrue(guard.limits_match(actual, spec))

    def test_reject_requested_unbounded_or_unsupported_profile(self):
        for memory, tasks in [(0, 32), (8, 64), (2, 0), (2, 128)]:
            with self.subTest(memory=memory, tasks=tasks):
                with self.assertRaises(ValueError):
                    guard.resource_limits(memory, tasks)

    def test_reject_enforcement_mismatch_and_relaxed_safeguards(self):
        spec, actual = self.profile(2, 64)
        for name, value in [('memory.max', 'max'), ('memory.max', str(4 * 1024**3)),
                            ('memory.swap.max', '1'), ('pids.max', 'max'),
                            ('pids.max', '32'), ('nice', 0), ('affinity', [3, 4])]:
            with self.subTest(name=name, value=value):
                changed = copy.deepcopy(actual)
                changed[name] = value
                self.assertFalse(guard.limits_match(changed, spec))
        changed = copy.deepcopy(actual)
        changed['threads']['OMP_NUM_THREADS'] = '2'
        self.assertFalse(guard.limits_match(changed, spec))

    def test_reject_unsupported_manifest_even_if_runtime_matches(self):
        spec, actual = self.profile(2, 64)
        for key, value, runtime_key in [('memoryMax', 8 * 1024**3, 'memory.max'),
                                        ('tasksMax', 128, 'pids.max'),
                                        ('swapMax', 1, 'memory.swap.max')]:
            with self.subTest(key=key):
                changed_spec, changed_actual = copy.deepcopy(spec), copy.deepcopy(actual)
                changed_spec[key] = value
                changed_actual[runtime_key] = str(value)
                self.assertFalse(guard.limits_match(changed_actual, changed_spec))


if __name__ == '__main__':
    unittest.main()
