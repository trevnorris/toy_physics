Run from the repository root with the review packet still at the location shown in the scripts. The reporting replay reads existing checkpoints, native sources, archived guard receipts and the original four input spans. It writes to a temporary directory below the review scripts. Run each command once with a fresh guard directory.

```sh
python3 scripts/s11c_guarded_run.py --pool s15-r2-b --memory-gib 2 --log-directory _scratch/s11c/s15-r2-b-reference-reproduce -- /usr/bin/python3 -B _scratch/s15/reading_comparator_build_review/r2/b/reference.py
python3 scripts/s11c_guarded_run.py --pool s15-r2-b --memory-gib 2 --log-directory _scratch/s11c/s15-r2-b-publication-reproduce -- /usr/bin/python3 -B _scratch/s15/reading_comparator_build_review/r2/b/check_publication.py
```

The original reference ran under `_scratch/s11c/s15-r2-b-reference01`. The first report and mutation run used `_scratch/s11c/s15-r2-b-publication01`. It completed the report, source dimension mutation and plane sign mutation, then stopped on a reviewer-side assertion that compared integer dictionary keys from a Python function to their published JSON strings. The reviewer corrected that representation check and ran only the remaining controls with the command below. The separate stdout files are literal copies of the three guard outputs. The all-controls command above uses the corrected script.

```sh
python3 scripts/s11c_guarded_run.py --pool s15-r2-b --memory-gib 2 --log-directory _scratch/s11c/s15-r2-b-points01 -- /usr/bin/python3 -B _scratch/s15/reading_comparator_build_review/r2/b/check_publication.py --remaining
```
