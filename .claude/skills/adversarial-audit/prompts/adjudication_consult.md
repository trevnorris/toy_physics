# Adjudication Consult Record

Use this template for the structured Claude plus Codex consult record.

```yaml
candidate_id: "{CANDIDATE_ID}"
parameter_names: []
adversarial_report:
defense_report:
consult_type: adversarial_audit_adjudication
dry_run: false
outcome: FIND_STANDS | FIND_FAILS | PARTIAL
requires_user_gate: false
conceptual_change: false
participants:
  claude:
    position:
    evidence: []
  codex:
    position:
    evidence: []
decision:
  summary:
  rationale:
  citations: []
```
