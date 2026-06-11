# Phase C Adversarial Audit Wrapper

Candidate: `{CANDIDATE_ID}`

```yaml
{DRY_RUN_BLOCK}
```

If `dry_run: true`, this prompt is non-binding and must not produce a campaign verdict.

## Frozen Directive

The following directive is included verbatim. Do not edit or paraphrase it.

```markdown
{FROZEN_DIRECTIVE}
```

## Candidate

```yaml
{CANDIDATE_YAML}
```

## Primary Sources

{PRIMARY_SOURCES}

## Provenance Slice

```yaml
{PROVENANCE_SLICE}
```

## Benchmarks

External-match checks must use these sourced benchmark entries. Do not adjudicate a match from model memory.

```yaml
{BENCHMARKS}
```

## Graph Context

```yaml
{GRAPH_CONTEXT}
```

Output a markdown adversarial report only when this is not a dry run. For dry runs, stop after confirming the prompt has enough context and record no verdict.
