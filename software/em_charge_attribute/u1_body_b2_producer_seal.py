#!/usr/bin/env python3
"""Record producer digests before a generated output becomes a stage input."""

from __future__ import annotations

import argparse
from pathlib import Path

from u1_body_b2_common import dump_yaml, require, sha256_file


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--producer", required=True)
    p.add_argument("--path", type=Path, action="append", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    try:
        rows = []
        for path in a.path:
            path = path.resolve()
            require(path.is_file(), "B2_A1_PRODUCER_DIGEST", f"missing:{path}")
            rows.append({"path": str(path), "sha256": sha256_file(path), "producer": a.producer})
        dump_yaml(a.output, {"schema_version": "U1_PHASE_B2_PRODUCER_SEAL_V1", "status": "PASS", "producer": a.producer, "outputs": rows})
        print(f"B2_PRODUCER_SEAL: PASS producer={a.producer} outputs={len(rows)}")
        return 0
    except Exception as exc:
        print(str(exc))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
