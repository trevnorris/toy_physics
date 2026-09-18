"""Exercise bootstrap control flow without downloads or a user's elan directory."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class SetupTests(unittest.TestCase):
    def fixture(self, root, existing_elan=False):
        source = Path(__file__).resolve().parent
        for name in ['setup.sh', 'requirements.txt', 'lean-toolchain', 'lake-manifest.json', 'lakefile.toml']:
            shutil.copyfile(source / name, root / name)
        bindir = root / 'bin'
        bindir.mkdir()
        # Use executable stubs solely at download/install boundaries. All
        # shell flow, file checks, logging and error propagation run normally.
        program = f'''#!{sys.executable}
import json,os,sys
from pathlib import Path
name=Path(sys.argv[0]).name;args=sys.argv[1:]
with open(os.environ['PDE_TEST_CALLS'],'a') as f:f.write(json.dumps([name,*args])+'\\n')
if name=='curl':
 p=Path(args[args.index('-o')+1])
 p.write_text(''' + repr('''#!/bin/sh
printf '%s\n' "$*" > installer-args.txt
mkdir -p "$ELAN_HOME/bin"
cp "$PDE_TEST_ELAN_STUB" "$ELAN_HOME/bin/elan"
''') + ''')
elif name=='python-test' and args[:2]==['-m','venv']:
 p=Path('.venv/bin/python');p.parent.mkdir(parents=True,exist_ok=True)
 p.write_bytes(Path(sys.argv[0]).read_bytes());p.chmod(0o755)
elif name=='python-test':
 os.execv(sys.executable,[sys.executable,*args])
elif name=='lake' and os.environ.get('PDE_TEST_FAIL_CACHE') and args[:2]==['exe','cache']:
 sys.exit(23)
'''
        for name in ['git', 'curl', 'zstd', 'tar', 'cc', 'lake', 'python-test', 'elan-stub']:
            p = bindir / name
            p.write_text(program)
            p.chmod(0o755)
        if existing_elan:
            shutil.copyfile(bindir / 'elan-stub', bindir / 'elan')
            (bindir / 'elan').chmod(0o755)
        env = {**os.environ, 'PATH': str(bindir) + ':/usr/bin:/bin',
               'ELAN_HOME': str(root / 'elan'), 'PYTHON': str(bindir / 'python-test'),
               'PDE_TEST_ELAN_STUB': str(bindir / 'elan-stub'),
               'PDE_TEST_CALLS': str(root / 'calls.jsonl')}
        return env

    def run_setup(self, root, env, *args):
        return subprocess.run(['bash', 'setup.sh', *args], cwd=root, env=env,
                              text=True, capture_output=True, timeout=10)

    def test_missing_elan_bootstraps_without_default_or_shell_modification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env = self.fixture(root)
            pins = {p: (root / p).read_bytes() for p in ['lake-manifest.json', 'lean-toolchain', 'lakefile.toml']}
            result = self.run_setup(root, env)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((root / 'installer-args.txt').read_text().strip(),
                             '-y --default-toolchain none --no-modify-path')
            calls = [json.loads(line) for line in (root / 'calls.jsonl').read_text().splitlines()]
            self.assertIn(['elan', 'toolchain', 'install', (root / 'lean-toolchain').read_text().strip()], calls)
            self.assertIn(['python', '-m', 'pip', 'install', '--disable-pip-version-check', '-r', 'requirements.txt'], calls)
            self.assertIn(['python', 'verify.py', '--doctor'], calls)
            self.assertFalse(any(c[:2] == ['lake', 'update'] for c in calls))
            self.assertFalse(any(c == ['lake', 'build'] for c in calls))
            self.assertTrue(any(c[:2] == ['lake', 'build'] and '+Physlib.Units.Dimension:olean' in c for c in calls))
            self.assertEqual(pins, {p: (root / p).read_bytes() for p in pins})

    def test_missing_manifest_stops_before_downloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env = self.fixture(root)
            (root / 'lake-manifest.json').unlink()
            result = self.run_setup(root, env)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Missing lake-manifest.json', result.stderr)
            calls = (root / 'calls.jsonl').read_text()
            self.assertNotIn('curl', calls)
            self.assertNotIn('lake', calls)

    def test_existing_elan_skips_download_and_cache_failure_is_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            env = self.fixture(root, existing_elan=True)
            env['PDE_TEST_FAIL_CACHE'] = '1'
            result = self.run_setup(root, env)
            self.assertEqual(result.returncode, 23)
            self.assertIn('Setup failed', result.stdout)
            self.assertIn('setup.log', result.stdout)
            calls = (root / 'calls.jsonl').read_text()
            self.assertNotIn('curl', calls)
            self.assertNotIn('--doctor', calls)


if __name__ == '__main__':
    unittest.main()
