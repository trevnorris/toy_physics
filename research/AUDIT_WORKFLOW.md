# Paper Audit Workflow

Use this for each `research/<paper_name>/` directory.

1. Inventory the paper bundle.
   - Confirm the directory has `paper/*.tex`, `mathematica/*.wl`, and any `scripts/*.py`.
   - Note which manuscript sections each script is supposed to support.

2. Run every verification script.
   - Mathematica: `math -script <path/to/file.wl>`
   - Python: `python3 <path/to/file.py>`
   - Record exit codes and keep the raw stdout for comparison.

3. Normalize Mathematica `Output:` blocks.
   - Every `.wl` must end with `Output:` followed by the file's stdout.
   - Compare the embedded block to the live run.
   - Ignore the environment-specific OpenMP SHM warning when diffing.

4. Cross-check scripts against the manuscript.
   - Verify that the script assumptions match the current preferred branch in the paper.
   - Fix stale local paths after repo moves.
   - If a script implements a rejected or superseded branch, either rewrite it to the current paper or clearly quarantine it.

5. Prefer fixing scripts over weakening claims.
   - If the paper already has the cleaner formulation and the code is stale, update the code.
   - If the code exposes a real derivation gap, tighten the paper wording.

6. Re-run changed files and rebuild the paper.
   - Re-run every modified `.wl` and `.py`.
   - Rebuild the PDF with `pdflatex` twice.
   - Check for LaTeX errors and note any non-blocking warnings.

7. Close out with a short summary.
   - Which scripts now pass.
   - Whether embedded `Output:` blocks match.
   - Any remaining caveats or places where the paper is still interpretive rather than derived.
