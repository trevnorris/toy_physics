# Measurements: S15 comparator decision list, review pass 1 lookups

Generator: `_scratch/s15/gen/s15_cmpdl_p1_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output. The legs' scripts were not re-run here; their filed stdout is retrieved.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/cmp/cmp_dl_baseline_p1.sha256
research/pde_ledger_v3/directives/S15_comparator_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_comparator_setup_lookups.md: OK
```

```
$ sha256sum _scratch/s15/cmp/cmp_dl_review_prompt_p1.md
d4f1b5dc4f76160577133b305f49ec50ecc6817ad3bcab06c394896f368f83a8  _scratch/s15/cmp/cmp_dl_review_prompt_p1.md
```

```
$ cat /var/projects/_review_scratch/s15-cmpdl-p1-codex/DONE /var/projects/_review_scratch/s15-cmpdl-p1-grok/DONE; ls -la /var/projects/_review_scratch/s15-cmpdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-cmpdl-p1-grok/REVIEW.md | cut -c1-110
exit 0
exit 0
-rw-rw-r-- 1 trevnorris trevnorris  2486 Oct 24 20:19 /var/projects/_review_scratch/s15-cmpdl-p1-codex/REVIEW.
-rw-rw-r-- 1 trevnorris trevnorris 13443 Oct 24 20:37 /var/projects/_review_scratch/s15-cmpdl-p1-grok/REVIEW.m
```

```
$ git status --porcelain | /usr/bin/grep -v -F -e 'S15_comparator_decision_list_review_p1_lookups.md' -e 'S15_comparator_decision_list_review_p1_disposition.md' | diff - <(sed -n '1,3p' _scratch/s15/cmp/status_before_cmp_dl_p1.txt) && echo 'status unchanged apart from the p1 disposition and lookups files'
status unchanged apart from the p1 disposition and lookups files
```

```
$ sed -n '4p' _scratch/s15/cmp/status_before_cmp_dl_p1.txt; git rev-parse HEAD
3192123b3c9029ec949197efecf18e5e1f86fcdb
3192123b3c9029ec949197efecf18e5e1f86fcdb
```

## M1. The verdicts

```
$ /usr/bin/grep -n -E '\*\*CLEAR|\*\*FINDINGS|^## \*\*CLEAR' /var/projects/_review_scratch/s15-cmpdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s15-cmpdl-p1-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-140
s15-cmpdl-p1-codex/REVIEW.md:35:**CLEAR**: nothing outstanding changes what the comparator would measure or what may be claimed from it.
s15-cmpdl-p1-grok/REVIEW.md:3:**CLEAR.** Nothing outstanding in `directives/S15_comparator_decision_list.md` changes what the comparator wou
s15-cmpdl-p1-grok/REVIEW.md:97:## **CLEAR**
```

## M2. Selected leg stdout

```
$ /usr/bin/grep -n -F 'current_bytes_equal' /var/projects/_review_scratch/s15-cmpdl-p1-codex/12_stream_audit.stdout | cut -c1-160
43:current_bytes_equal_3192123b research/pde_ledger_v3/scripts/S15_two_drain_region.py: True
44:current_bytes_equal_3192123b research/pde_ledger_v3/mathematica/S15_two_drain_region.wl: True
45:current_bytes_equal_91248b67 research/pde_ledger_v3/directives/S15_decision_list.md: True
46:current_bytes_equal_f08db34a research/pde_ledger_v3/scripts/S16_cross_engine_comparator.py: True
47:current_bytes_equal_f08db34a research/pde_ledger_v3/scripts/test_S16_cross_engine_comparator.py: True
```

```
$ cat /var/projects/_review_scratch/s15-cmpdl-p1-grok/out/group_coverage.stdout | tail -3
absent 0
```

```
$ /usr/bin/grep -n -E '^(records|shared tags containing|dimension_argument/name count|partition/unlabelled_coefficient count)' /var/projects/_review_scratch/s15-cmpdl-p1-grok/out/adjudication_certificate.txt | cut -c1-160
1:records py=2231 wl=1370 unique py=1360 wl=1354 shared=401
2:shared tags containing labels_ : 0
69:dimension_argument/name count=208 unique=208
70:partition/unlabelled_coefficient count=27 distinct=['0']
```

## M3. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-cmpdl-p1-codex/EVIDENCE.md /var/projects/_review_scratch/s15-cmpdl-p1-codex/PRE_ARTIFACT_REQUIREMENTS.md /var/projects/_review_scratch/s15-cmpdl-p1-codex/12_stream_audit.py /var/projects/_review_scratch/s15-cmpdl-p1-codex/12_stream_audit.stdout /var/projects/_review_scratch/s15-cmpdl-p1-codex/13_contract_crosswalk.py /var/projects/_review_scratch/s15-cmpdl-p1-codex/13_contract_crosswalk.stdout /var/projects/_review_scratch/s15-cmpdl-p1-grok/scripts/*.py /var/projects/_review_scratch/s15-cmpdl-p1-grok/out/*.txt /var/projects/_review_scratch/s15-cmpdl-p1-grok/out/*.stdout | sed 's|/var/projects/_review_scratch/||'
42466f925022cc836c0cf9ad66b3869cafd478c60436a4dc6e6fbf515ea3240d  s15-cmpdl-p1-codex/EVIDENCE.md
8aa26d6b442d51f31d2cb4e8d29c3e4ca43c5d7bef6dab076baf3c20c3ed04ff  s15-cmpdl-p1-codex/PRE_ARTIFACT_REQUIREMENTS.md
40ec7d1aabc165e631c965033afb8c42c4cd47ad8c7e7bf67656d044a4d2be84  s15-cmpdl-p1-codex/12_stream_audit.py
1be8d90d766c7bf92bed4e88839afe350268546b3ce7402a88df714fe90d98bc  s15-cmpdl-p1-codex/12_stream_audit.stdout
e4f2af708f795ea0a53b7c0f28a9d0f6bdf0fec4a835791d4b1261a4815f73cf  s15-cmpdl-p1-codex/13_contract_crosswalk.py
2b4ec4a21b8c2e63219dbdd92ab5de250668ea3fe157abdf89438d87ca5c74ef  s15-cmpdl-p1-codex/13_contract_crosswalk.stdout
c3d7ed266f635beffe58cab1a9ab5aca598f9d58b152939c6c290448d680367e  s15-cmpdl-p1-grok/scripts/adjudication_certificate.py
c1f27383d63fecbb6f245df28c345cff1ad6551c5da66ee575b0da1ef1ac582d  s15-cmpdl-p1-grok/scripts/group_coverage.py
1dc38b46fb18c3af06cb305bb95a36949e12498e538fbcfcf1e872385753f0b0  s15-cmpdl-p1-grok/scripts/record_identity.py
1c6372c6215a42f703d3306babd03ba3e3743fd06ced280c486e2552dbef67e8  s15-cmpdl-p1-grok/scripts/repeated_records.py
f6d8ba1b00357bf5bf8c4a48c2164cc32d08962df156ef26fbb7387492a65f09  s15-cmpdl-p1-grok/scripts/targeted_prefixes.py
404c3040775238553d0668b238e2cf2b92a2dff1d4cec5a670871156fe1d126f  s15-cmpdl-p1-grok/scripts/transcript_structure.py
f4c11f2ed145dd70e41453093fdfb60ed01ddd282a1d5e0fd98d594749511bb8  s15-cmpdl-p1-grok/out/adjudication_certificate.txt
6e647c3e8bae9eb38ddaa5453c8a4c0298e5a8f482b7e8cc9909c6ed5bb2539e  s15-cmpdl-p1-grok/out/record_identity.txt
571b5e5a55b199f819c81dbf5bfb1de5ec734d6eb54106f3d2d48bc42e85ddd0  s15-cmpdl-p1-grok/out/repeated_records.txt
7169c1fe0349c9782acf6a47d16f952488f99e35b08e51f8b3685b398ae5c9a8  s15-cmpdl-p1-grok/out/targeted_prefixes.txt
381a12b537b7012a18d5f502be3a0185ae6af6fcba83176b095f7b37fd91d903  s15-cmpdl-p1-grok/out/transcript_structure.txt
46fa52b77dbd7c8a80ea0066239b13c1899d6dfeb473f6260b6681bd7ba582da  s15-cmpdl-p1-grok/out/adjudication_certificate.stdout
afc0a1a4f54f61f09477641d1bd519edd4fab0f8aab7770fe102412b102a7dad  s15-cmpdl-p1-grok/out/group_coverage.stdout
```
