"""Review dispatcher path validation; standard library only."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import s11c_review_cas as runner


class ReviewProbePaths(unittest.TestCase):
    def test_assigned_workspace_python_only(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            work = base / 'leg'; work.mkdir()
            script = work / 'probe.py'; script.write_text('print(1)\n')
            with patch.object(runner, 'STORE', base):
                self.assertEqual(runner.validate_script(work, script), (work, script))
                other = work / 'data.json'; other.write_text('{}')
                with self.assertRaises(ValueError):
                    runner.validate_script(work, other)

    def test_external_or_symlink_escape_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            work = base / 'leg'; work.mkdir()
            outside = base / 'outside.py'; outside.write_text('print(1)\n')
            link = work / 'linked.py'; link.symlink_to(outside)
            with patch.object(runner, 'STORE', base):
                for script in (outside, link):
                    with self.assertRaises(ValueError):
                        runner.validate_script(work, script)


if __name__ == '__main__':
    unittest.main()
