# Measurements: S14a part 1 comparator decision list, review pass 1 lookups (generated 2026-10-21 10:38)

Generator: `_scratch/s14a/gen/cmp_dl_p1_lookups.sh` (grep/sed/sha256sum/git only; this file is written by it).
Mechanical lookups only.

## D1. The pass: prompt, artifact baseline, tree, reports

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/cmp/cmp_dl_review_prompt_p1.md
e6e48e28a2a1a6b4927019529cad4d38eb7e455148b085ae2fca1342010c1572  /var/projects/toy_physics/_scratch/s14a/cmp/cmp_dl_review_prompt_p1.md
```

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/cmp/cmp_dl_baseline_p1.sha256
research/pde_ledger_v3/directives/S14a_1_comparator_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_setup_lookups.md: OK
/var/projects/toy_physics/_scratch/s14a/cmp/cmp_dl_review_prompt_p1.md: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/cmp/status_before_cmp_dl_p1.txt && echo 'status unchanged'
3d2
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_decision_list_review_p1_lookups.md
```

```
$ cat /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/DONE /var/projects/_review_scratch/s14a1-cmpdl-p1-grok/DONE
exit 0
exit 0
```

```
$ sed -n '1,12p' /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/codex_review.txt
OpenAI Codex v0.162.1
--------
workdir: /var/projects/_review_scratch/s14a1-cmpdl-p1-codex
model: gpt-6.1-sol
provider: openai
approval: never
sandbox: workspace-write [workdir, /tmp, $TMPDIR]
reasoning effort: xhigh
reasoning summaries: none
session id: 01a15f39-8ed0-73b0-bb8f-ba52beac5535
--------
user
```

```
$ grep -nE '^\*\*(FINDINGS|CLEAR)' /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/REVIEW.md /var/projects/_review_scratch/s14a1-cmpdl-p1-grok/REVIEW.md
/var/projects/_review_scratch/s14a1-cmpdl-p1-codex/REVIEW.md:13:**FINDINGS**
/var/projects/_review_scratch/s14a1-cmpdl-p1-grok/REVIEW.md:5:**CLEAR**: nothing outstanding changes what the comparator would measure or what may be claimed from it.
/var/projects/_review_scratch/s14a1-cmpdl-p1-grok/REVIEW.md:91:**CLEAR**: nothing outstanding changes what the comparator would measure or what may be claimed from it.
```

```
$ grep -nE '^[0-9]+\. \*\*' /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/REVIEW.md
15:1. **High — Whole-tag injectivity prevents complete comparison of objects split across tags on one side and grouped on the other.**
```

## D2. The list text Codex 1 quotes

```
$ grep -nF 'Each row pairs one tag of each stream' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_comparator_decision_list.md
46:   - Each row pairs one tag of each stream. It cites the object both denote, by the S14a part 1 decision list's item
```

```
$ grep -nF 'The table is injective' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_comparator_decision_list.md
48:   - The table is injective. No row is made by similarity of names, by position, or because two payloads are
```

```
$ grep -nF 'which components of the two payloads denote the same object' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_comparator_decision_list.md
57:     the comparator establishes which components of the two payloads denote the same object, citing the line in each
```

```
$ grep -nF 'listed as unjoined with' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_comparator_decision_list.md
51:   - Each tag occurs once per stream, and a repeat is reported. Every tag is on a row, or is listed as unjoined with
55:     listed as unjoined with that reason.
```

## D3. The governing list's objects Codex 1 cites

```
$ grep -nF 'every term that leaves at each face' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
72:     - every term that leaves at each face;
```

```
$ grep -nF 'Dimensional homogeneity' /var/projects/toy_physics/research/pde_ledger_v3/directives/S14a_1_decision_list.md
101:4. **Dimensional homogeneity** of the projected balances and of objects 2 and 3.
```

## D4. Split versus grouped packaging in the transcripts (tag names; payloads cut at 120 characters)

```
$ grep -nE '^General/material/face [01]/closed-projection leaving term: ' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out | cut -c1-120
88:General/material/face 0/closed-projection leaving term: General_below_n(x, y, z, General_face_lower(x, y, z, t), t)*G
95:General/material/face 1/closed-projection leaving term: -General_above_n(x, y, z, General_face_upper(x, y, z, t), t)*
```

```
$ grep -nE '^general projection/material/each closed-projection face term: ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | cut -c1-120
63:general projection/material/each closed-projection face term: {-(n1[x1, x2, x3, wLower[x1, x2, x3, t], t]*u1c4[x1, x2
```

```
$ grep -nE '^General/dimensions/material/flux [012]/dimension: ' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out
372:General/dimensions/material/flux 0/dimension: (1, -2, -1)
375:General/dimensions/material/flux 1/dimension: (1, -2, -1)
378:General/dimensions/material/flux 2/dimension: (1, -2, -1)
```

```
$ grep -nE '^general projection/dimensions/material/plane_flux: ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
294:general projection/dimensions/material/plane_flux: {{1, -2, -1}, {1, -2, -1}, {1, -2, -1}}
```

## D5. The engine lines that emit them

```
$ sed -n '193,196p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_projection.py
        for index, (edge, inside, atom) in enumerate(zip(leaving, interior_leaving, atoms)):
            object_(label + f"/{kind}/face {index}/interior leaving term", inside)
            object_(label + f"/{kind}/face {index}/closed-projection leaving term", edge)
            if kind == "material":
```

```
$ sed -n '322,323p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_projection.py
        for index, flux in enumerate(current["flux"]):
            audit_dimensions(prefix + f"/flux {index}", flux)
```

```
$ sed -n '168p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S14a_1_projection.wl
  alone[prefix <> "/" <> eq <> "/each closed-projection face term", closedFace];
```

```
$ sed -n '249,250p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S14a_1_projection.wl
  Do[alone[package[[1]] <> "/dimensions/" <> kind <> "/" <> key, dimension[obj[key]]],
    {key, {"density", "plane_flux", "faces", "atoms", "derivative", "left", "source", "supplied_source"}}];
```

## D6. Codex's witness script and its stdout

```
$ sha256sum /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/11_tag_grouping_witness.py /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/11_tag_grouping_witness.stdout
b297121c4fc1763901384b14a44c6545833c3c22f7ea28b84e7d2874b1873648  /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/11_tag_grouping_witness.py
343b8d0183c3f4f51bba8b36efdf4dc4deb0eeb9f2f08a06330a74dd469e31b3  /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/11_tag_grouping_witness.stdout
```

```
$ grep -nE 'Maximum directly bound|permits 3/3' /var/projects/_review_scratch/s14a1-cmpdl-p1-codex/11_tag_grouping_witness.stdout
62:Maximum directly bound counterpart components under tag injectivity: 1/3
66:An injective table on (tag, component-path) endpoints permits 3/3 role-correct bindings with this same stored packaging.
```
