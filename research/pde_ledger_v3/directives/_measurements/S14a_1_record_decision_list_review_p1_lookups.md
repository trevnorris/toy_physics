# Measurements: S14a part 1 record decision list review pass 1 lookups (generated 2026-10-21 15:04)

Generator: `_scratch/s14a/gen/recdl_p1_lookups.sh` (grep/sed/sha256sum/diff/git/cat only; this file is written
by it). Mechanical lookups only.

## D1. The pass: baseline, tree, prompt, reports

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/rec/recdl_baseline_p1.sha256
research/pde_ledger_v3/directives/S14a_1_record_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_record_setup_lookups.md: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/rec/status_before_recdl_p1.txt && echo 'status unchanged'
3d2
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_record_decision_list_review_p1_lookups.md
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/rec/recdl_review_prompt_p1.md
5cf115161db9d710a606cb72b4b3a0679ba265623ee36a30ef7f68b98bfa530e  /var/projects/toy_physics/_scratch/s14a/rec/recdl_review_prompt_p1.md
```

```
$ cat /var/projects/_review_scratch/s14a1-recdl-p1-codex/DONE; grep -A1 'tokens used' /var/projects/_review_scratch/s14a1-recdl-p1-codex/codex_review.txt | tail -1; cat /var/projects/_review_scratch/s14a1-recdl-p1-grok/DONE
exit 0
174,817
exit 0
```

```
$ grep -nE 'CLEAR|FINDINGS' /var/projects/_review_scratch/s14a1-recdl-p1-codex/REVIEW.md | cut -c1-160; grep -nE '^\*\*Verdict|^\*\*FINDINGS|^1\. \*\*Medium' /var/projects/_review_scratch/s14a1-recdl-p1-grok/REVIEW.md | cut -c1-200
32:**CLEAR**: nothing outstanding changes what the record may claim or what the register pass may add.
5:**Verdict.** **FINDINGS.** One item changes what the register will still say about the accepted S12 pass.
108:**FINDINGS**
110:1. **Medium.** Item 9 updates the S12 pass heading and then says nothing else in that pass changes, so `SUBSTRATE_REQUIREMENTS.md:1300–1302` still says the accepted pass is awaiting the narrow c
```

## D2. Grok 1: the S12 pass source sentence

```
$ grep -nE 'Stale status text|Nothing else in the S12' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_record_decision_list.md | cut -c1-200
108:     J2 are updated as the S12 record's round-3 disposition states them (L8). Nothing else in the S12 record or in
```

```
$ sed -n '/^## What must be true/,/^## Exclusions/p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_record_decision_list.md | grep -nE 'ac909aaf|J1|N1|STATUS' | cut -c1-200
58:   - The S12 record and register pass were accepted at `ac909aaf`. The stale status text and the known items J1 and
61:   - The run record's construction-order sentence (N1) states the authorship and order the build's round-2
63:   - `../../STATUS.md` states where S14a part 1 stands, as the record does.
```

```
$ sed -n '1298,1302p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
### S12 pass — local order-source compatibility and separate inventories · 2026-10-20; written for review

Source: `steps/S12_order_source.md`, §§1–7. Its constructions are accepted at `a3da43d9`, and the
comparator with its engine-row note at `fe485b73`; this is the round-3 repair of reviewed baseline
`fa5d1396`, awaiting the narrow check of H1/H2 (M19).
```

```
$ sed -n '/^\*\*Stale status text/,$p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S12_record_r3_review_disposition.md
**Stale status text.** The accepted files still describe S12's record and register pass as written for, or under,
review: the record header, the register's status and S12 pass heading, and `STATUS.md`. The S9b record's acceptance
left the same kind of text, and S12's item 9 then repaired it for the light closeout. The next step's record list
updates these lines to this acceptance, and picks up J1 and J2.
```

```
$ grep -nE 'written for review|under review|awaiting the narrow' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md /var/projects/toy_physics/STATUS.md | cut -c1-180
/var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md:3:**Status, 2026-10-20:** record and register pass written for review. The constructions are accepted at
/var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md:5:note below. This is the round-3 repair of reviewed baseline `fa5d1396`, awaiting the narrow check of H1
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:6:split, is accepted at `e070c7f3`. The S12 record and register pass (2026-10-20) are written for review;
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1298:### S12 pass — local order-source compatibility and separate inventories · 2026-10-20; written fo
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:1302:`fa5d1396`, awaiting the narrow check of H1/H2 (M19).
/var/projects/toy_physics/STATUS.md:1:# Project status — light sector closed for now; S12 record written for review
/var/projects/toy_physics/STATUS.md:65:| S10/S11 exports pin their pre-2026-09-11 producers. Regenerate them from the current producers under review, or revert the producer edits. 
/var/projects/toy_physics/STATUS.md:94:**S12's record is written for review:** its accepted constructions and comparator report the source
```

## D3. J1 and J2 as the S12 round-3 disposition locates them

```
$ sed -n '/^\*\*J1\./,/^## The leg disagreement/p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S12_record_r3_review_disposition.md
**J1. The record's account of its own edits no longer covers two status-only edits.** (Claude 1.) **CONFIRMED.**
- **What was there.** At `fa5d1396`, the S12 pass said "The light-closeout status lines now record its acceptance".
  The record said "The earlier light-closeout acceptance edits at `e070c7f3` are retained". Both sentences are gone
  now (counts 1 → 0, K3).
- **What is now unreported.** Two status edits are mentioned nowhere in the record's account (K3):
  - the S9b oracle record's status line, from "its own review is pending" to accepted at `e070c7f3`;
  - the register header's light-closing sentence.
- **Standing.** Both edits are correct, and the record decision list's item 9 orders them. Neither changes what S12
  claims or hands on.

**J2. A true count cites sections that do not hold it.** (Claude 2.) **CONFIRMED.**
- **The claim.** "The present repair of `fa5d1396` is also 18 → 18, with no additional source gain", cited to
  M19/M20 (record l.340–343, S12 pass l.1315–1317).
- **The citation.**
  - M19 and M20 hold no count.
  - The count is in M12: "before round-3 repair fa5d1396: entry headings cardinality 18" and "after working tree:
    entry headings cardinality 18".
  - No `**source**` field changed since `fa5d1396` (diff count 0) (K4).
- **Standing.** The claim is true, but its command is in an uncited section. It changes nothing claimed.

## The leg disagreement
```

```
$ sed -n '340,343p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S12_order_source.md
The present repair of `fa5d1396` is also 18 → 18, with no additional source gain. It carries the
real-fraction/A13 condition into `R-S12-02`, its by-target summary and the common handoff of the
S12 pass: the printed requirements rest on the real-fraction branch and A13 stays OPEN until
S1's gate decides. It also supplies the accepted-pass edit account below (M19/M20).
```

```
$ sed -n '1315,1317p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
The present repair of `fa5d1396` is also 18 → 18, with no additional source gain. It adds the
real-fraction/A13 qualification to the S12 handoffs and reports the accepted-pass edits below;
it changes no owner, computed object or comparison claim (M19/M20).
```

```
$ grep -nE 'carried as known items|J1 and J2' /var/projects/_review_scratch/s14a1-recdl-p1-grok/REVIEW.md | cut -c1-220
17:- the S12 acceptance at `ac909aaf`, with J1 and J2 carried as known items and N1 carried as the stale construction-order sentence;
43:**9, except the sentence in Finding 1.** The acceptance commit is `ac909aaf` ("Accept the S12 record and its register pass (user: accept, note both)"), which matches L1. J1 and J2 are the known items L8 quotes as carr
64:> The stale status text and the known items J1 and J2 are updated as the S12 record's round-3 disposition states them (L8). Nothing else in the S12 record or in the register's S12 pass changes.
68:> **Stale status text.** The accepted files still describe S12's record and register pass as written for, or under, review: the record header, the register's status and S12 pass heading, and `STATUS.md`. … The next 
```
