"""Make a partial unittest-discovery run visibly fail as the wrong runner."""

from __future__ import annotations

import os
import unittest


class IntendedRunnerContract(unittest.TestCase):
    def test_full_reduction_suite_requires_pytest(self) -> None:
        self.assertTrue(
            "PYTEST_CURRENT_TEST" in os.environ,
            "WRONG_RUNNER: use `python -m pytest`; unittest discovery does not "
            "collect the full reduction test population",
        )
