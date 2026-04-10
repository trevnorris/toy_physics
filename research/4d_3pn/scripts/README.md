# 4d_3pn SymPy Archive

This directory contains the referee-facing SymPy verification archive for the
`4d_3pn` paper.

## Layout

- `3pn_referee_master_sympy_audit.py`: top-level hash-checking runner
- `3pn_*_audit.py`: extracted stage audits, one file per theorem step

## Intended Use

Run the full paper-supportive suite with:

```bash
python3 research/4d_3pn/scripts/3pn_referee_master_sympy_audit.py
```

The master runner replays the stage files in order, checks each file against its pinned
SHA256 digest, and stops on the first failed symbolic identity.

## Scope

This archive verifies the symbolic theorem chain used by the paper within the declared
closure framework. It does not claim to solve a stronger problem than the manuscript states.
