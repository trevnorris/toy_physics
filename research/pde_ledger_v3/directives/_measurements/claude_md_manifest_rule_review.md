# Review: the run-manifest rule (CLAUDE.md, AGENTS.md, docs/development_pipeline.md, .gitignore)

User decision of 2026-11-02 ("Manifest only (Recommended)"): engine transcripts are regenerable, so a new one is not committed; each run commits a run manifest. The orchestrator wrote the amendment, so the legs are Codex (gpt-6.1-sol, xhigh) and Grok (grok-4.7, high), with one identical prompt (`directives/_legs/claude_md_manifest_rule_review_prompt_r0.md`; each leg got one leading line naming it). The leg reports are copied beside it in `directives/_legs/`.

## Round 0

Artifact: the working-tree diff of the four files against `eceafc57`, sha256 `7f0c7114…` (CLAUDE.md), `ce029cc5…` (AGENTS.md), `471a56c8…` (pipeline doc), `e15360bd…` (.gitignore).

| Leg | Report sha256 (recorded before either was read) | Verdict |
|---|---|---|
| Codex | `3974acdef9d30e696cd88226ebe3c7355998548ed66e6ad05c940fc09eacbd95` | CLEAR |
| Grok | `16251aa3f7ad32c261e4d5571de9d513f711b91354eeaa4bbb4b77b522318edd` | NOT CLEAR, one finding |

### Finding M1 (Grok): the staging ban covered every `.out`

Grok: "`CLAUDE.md:284` forbids staging every `.out`. The decision, `AGENTS.md:31`, and `.gitignore:28-31` cover new v3 engine transcripts only. Ninety-five other tracked `.out` files, including v2 audit transcripts and `research/pde_ledger_v3/_measurements/*.out`, stay addable under gitignore and become unstageable under the governing sentence."

Lookups (orchestrator, mechanical):

```
$ grep -n -F 'Never stage a `.out`' CLAUDE.md   # on the round-0 baseline
284:- ⛔ Never stage a `.out`, and ⛔ never `git add -f` one. When a run replaces a transcript that is still tracked,
$ git ls-files '*.out' | grep -v -c -E '^research/pde_ledger_v3/(scripts|mathematica)/out/[^/]+\.out$'
95
$ git ls-files 'research/pde_ledger_v3/_measurements/*.out'
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_EG.out
research/pde_ledger_v3/_measurements/S11c_c2_adjudication_verify_F.out
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_adjudication_verify.out
research/pde_ledger_v3/_measurements/S11c_c2_export_repair_rereview_claude_recompute_artifact.out
```

Disposition: **confirmed.** The sentence was meant for engine transcripts but read as a ban on every `.out`, including E1 stdout files in `_measurements/`.

What must be true: any `.out` other than a new engine transcript in the two v3 `out/` directories (or a preserved copy under them) can be staged and committed as ordinary evidence.

Repair (CLAUDE.md, Repository infrastructure): the ban now names engine transcripts only, "a `.out` in one of those two `out/` directories or a preserved copy under them", and says any other `.out` is committed as before.

Codex found no defect and stated conditions the rule must keep (manifests keep every relied-on line; comparators still read complete `.out` streams, retrieved or regenerated; ablation triples stay in `_measurements/`). The amended text already requires each of them.

## Round 1 (scoped to the M1 repair)

Artifact: sha256 `1c1ce0e6…` (CLAUDE.md); the other three files unchanged from round 0. Prompt: `directives/_legs/claude_md_manifest_rule_review_prompt_r1.md` (sha256 `4c53d3d3…`). The legs were handed this record as it stood after round 0 (sha256 `6184e7df…`).

| Leg | Report sha256 (recorded before either was read) | Verdict |
|---|---|---|
| Codex | `1dbae11179ff3b60cef502a35bca8befde0f5152fd7fae121c77e6e34322a77c` | CLEAR |
| Grok | `6c4a14def5b68bcb815ebed60b89dcfc3225c41a87eb9299156e4b7cdcad44f8` | CLEAR |

Both legs reconstructed the round-0 `CLAUDE.md` by putting the round-0 bullet back, and got its hash `7f0c7114…`. So the only change since round 0 is the repaired bullet.

Both legs ran `git check-ignore` against all 95 tracked `.out` files outside the two engine `out/` directories: none is ignored, with or without `--no-index`. Both ran it on new and preserved transcripts inside those directories: all four are ignored by `.gitignore:28-31`. Both confirmed the `git rm --cached` clause is unchanged.

Grok noted that `**/paper/*.out` (LaTeX aux, unchanged) also matches a `.out` that is not an engine transcript. No tracked path has a `paper` component, and the rule predates this amendment. Not a finding.

Disposition: **clear.** Nothing outstanding changes what is committed, what evidence exists, or what may be claimed.
