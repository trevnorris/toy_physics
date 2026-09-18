"""Portable runner safety regressions. Run: python3 -m unittest -v test_verify."""
import json
import os
from pathlib import Path
import select
import subprocess
import sys
import tempfile
import unittest

import verify


class EvidenceTests(unittest.TestCase):
    def test_all_completed_control_sources_and_historical_diagnostics(self):
        # Cross-version report formats must replay without silently omitting
        # an older source mutation or treating an instrument failure as one.
        total = 0
        for name, spec in verify.specifications(list(verify.CONTRACTS)).items():
            records = {r['name']: r for r in verify.load(verify.BASE / spec['report'])['checks']}
            for control in spec['controls']:
                total += 1
                if control['expected'] == 'REJECTED':
                    original = records[control['name']]
                    self.assertTrue(verify.adjudicate(original['exit_status'], original['output'],
                        control['source'], control['required'], control['source_mutation']),
                        (name, control['name']))
                    self.assertFalse(verify.adjudicate(124, original['output'], control['source'],
                        control['required'], control['source_mutation']))
        self.assertEqual(total, 319)

    def test_mutation_source_drift_is_not_reconstructed_silently(self):
        record = {'name': 'drift', 'source': 'theorem x : True := by trivial\n',
                  'source_sha256': '0' * 64}
        with self.assertRaises(ValueError):
            verify.control_source(record)

    def test_only_intended_false_diagnostic_passes(self):
        source = 'theorem contract_control : (1 : Nat) = 2 := by norm_num\n'
        output = 'control.lean:1:45: error: unsolved goals\n⊢ False\n'
        self.assertTrue(verify.adjudicate(1, output, source, 'contract_control'))
        for bad in ['unknown identifier x', 'unexpected token', 'maximum heartbeats exceeded',
                    'failed to synthesize instance', 'warning: unused tactic',
                    'control.lean:1:3: error: unsolved goals\n⊢ False\n']:
            self.assertFalse(verify.adjudicate(1, output + bad, source, 'contract_control'), bad)
        self.assertFalse(verify.adjudicate(1, output, source, 'different_theorem'))
        self.assertFalse(verify.adjudicate(0, output, source, 'contract_control'))
        self.assertFalse(verify.adjudicate(1, 'control.lean:1:3: error: Tactic `rewrite` failed\n⊢ 1 = 2',
                                          source, 'contract_control', True))

    def test_wrapped_empty_and_forbidden_axioms(self):
        source = '#print axioms a\n#print axioms b\n#print axioms c\n'
        output = ("a depends on axioms: [propext,\n Classical.choice, Quot.sound]\n"
                  "b depends on axioms: []\nc does not depend on any axioms\n")
        self.assertEqual(verify.audit_axioms(output, source), 3)
        for forbidden in ['sorryAx', 'Lean.trustCompiler', 'customPhysics']:
            with self.assertRaises(ValueError):
                verify.audit_axioms(output.replace('propext', forbidden), source)
        with self.assertRaises(ValueError):
            verify.audit_axioms(output.splitlines()[0], source)

    def test_snapshot_paths_cannot_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.assertEqual(verify.inside(root, 'lean/../test'), root / 'test')
            for path in ['../outside', '/etc/passwd']:
                with self.assertRaises(ValueError):
                    verify.inside(root, path)
            (root / 'escape').symlink_to('/etc')
            with self.assertRaises(ValueError):
                verify.inside(root, 'escape/passwd')

    def test_completed_targets_include_s11_but_not_unfinished_d5(self):
        roots = [r for _, rs in verify.CONTRACTS.values() for r in rs]
        order = verify.build_order(roots)
        self.assertIn('S11D4Bulk', order)
        self.assertIn('S11NonlinearPole', order)
        self.assertNotIn('S11D5Invariants', order)
        self.assertEqual(len(order), len(set(order)))

    @unittest.skipUnless(sys.platform == 'linux', 'checks Linux process state')
    def test_timeout_stops_child_even_when_launcher_exits_on_term(self):
        # The child ignores SIGTERM; killing only its launcher would leak it.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            child_code = 'import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)'
            launcher = ('import subprocess,sys,time; from pathlib import Path; '
                        f'p=subprocess.Popen([sys.executable,"-c",{child_code!r}]); '
                        'Path("child.pid").write_text(str(p.pid)); time.sleep(60)')
            with self.assertRaises(subprocess.TimeoutExpired):
                verify.run_process([sys.executable, '-c', launcher], root, root / 'job.log', 1, dict(os.environ))
            pid = int((root / 'child.pid').read_text())
            try:
                fd = os.pidfd_open(pid)
            except ProcessLookupError:
                return
            try:
                self.assertTrue(select.select([fd], [], [], 2)[0], 'child survived timeout cleanup')
            finally:
                os.close(fd)
            state = Path(f'/proc/{pid}/stat')
            # A zombie is already dead; the init process reaps it separately.
            if state.exists():
                self.assertEqual(state.read_text().split()[2], 'Z')


if __name__ == '__main__':
    unittest.main()
