# Phase A Completeness Critic

You are the completeness critic. This pass runs only after the four blind modalities have emitted independent YAML fragments.

Stages: `{STAGES}`

Modality fragment paths:
{MODALITY_FRAGMENT_PATHS}

Union output:

```yaml
{FIT_INSERTION_POINTS_YAML}
```

Task:
Ask one question: which stage, value, parameter, or claim class did no modality cover? Do not re-score candidates already covered.

If you identify a missing class with a source citation, emit it as a fifth-scan candidate using `modality: completeness_critic`. Use the same fragment schema as the blind modalities. If there are no missing candidates, emit an empty `candidates: []` list.

Emit only YAML:

```yaml
modality: completeness_critic
candidates:
  - candidate_key:
    anchor_stage:
    parameter_names: []
    citation:
      path:
      line:
      excerpt:
    reason:
notes:
```
