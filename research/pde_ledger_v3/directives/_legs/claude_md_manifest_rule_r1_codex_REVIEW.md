codex — Independent review, round 1 (scoped): run-manifest rule

Reviewed only the repair at `CLAUDE.md:284–288`, against M1 and the three questions in the review request. Findings: none within this scope.

The four artifact hashes and the round-0 review-record hash match the supplied baselines before and after review: `sha256sum CLAUDE.md AGENTS.md docs/development_pipeline.md .gitignore research/pde_ledger_v3/directives/_measurements/claude_md_manifest_rule_review.md`, with its literal output in [baseline.txt](baseline.txt) and [baseline_final.txt](baseline_final.txt). The requested diff against `eceafc57` is in [artifact_diff.txt](artifact_diff.txt). Reconstructing the prior paragraph in scratch produced the supplied round-0 SHA256 `7f0c71140a0046c0f79c382379de61c8af47152fd9302c982f4486a3c9fbd3a5` ([command and output](round0_reconstruction.txt)); `git diff --no-index` then shows that only this paragraph changed ([command and output](round1_only_diff.txt)).

1. **M1's required condition holds for the ordinary evidence checked.** M1's exact requirement is recorded in the `sed -n 1,280p` output in [round0.txt](round0.txt). The repaired `CLAUDE.md:284–288`, retrieved with `sed -n '284,288{=;p;}' CLAUDE.md`, says literally:

   ```text
   - ⛔ Never stage an engine transcript, meaning a `.out` in one of those two `out/` directories or a preserved copy
     under them, and ⛔ never `git add -f` one. Any other `.out` (a `_measurements/` stdout file, a leg's derivation
     output, a v2 audit transcript) is ordinary evidence and is committed as before. When a run replaces a
     transcript that is still tracked, the commit that adds the new manifest removes the old file from the index
     (`git rm --cached`); history keeps it.
   ```

   Source: [repair_lines.txt](repair_lines.txt). `git ls-files '*.out'` is saved in [tracked_out.txt](tracked_out.txt). Excluding the two v3 `out/` trees leaves **95** tracked paths ([exact command/list](outside_out.txt), [count command/output](outside_count.txt)). Feeding that entire list to `git check-ignore -v --no-index --stdin` gives **empty stdout, empty stderr, exit 1** ([exact pipeline/output](outside_ignore.txt)). Thus none has a matching ignore rule, including when index tracking is disregarded; the repair expressly permits their ordinary evidence treatment.

2. **The transcript and preserved-copy ban remains.** The quoted definition still covers both v3 engine-transcript directories and their preserved copies. `git check-ignore -v --no-index --` with one hypothetical new transcript and one preserved `.out` in each directory matches `.gitignore:28`, `:29`, `:30`, and `:31`, respectively ([four exact paths and literal matching output](ignore_probes.txt)). The unchanged surrounding paragraph still says new transcripts are not committed to git or git-annex ([sed command/output](claude_context.txt)).

3. **The tracked-transcript removal sentence keeps its meaning; no new scoped defect found.** The round-0/current diff shows the same condition and action: when a run replaces a still-tracked transcript, the manifest commit removes that old transcript from the index with `git rm --cached`, retaining history. Only its wrapping changed ([literal diff](round1_only_diff.txt)). The added ordinary-evidence permission does not alter that condition or the manifest's required evidence ([context](claude_context.txt)).

Commit IDs: none created; comparison base `eceafc57`. Files to read: this report and its linked command/output records. Open questions: none within the scoped repair. Repository changes: none made.

Verdict: CLEAR
