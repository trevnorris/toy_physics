#!/usr/bin/env python3
"""Check what the ontology's cited return/drain lines actually state."""

from pathlib import Path


path = Path("/var/projects/toy_physics/docs/toy_model_ontology_summary.md")
lines = path.read_text(encoding="utf-8").splitlines()

for number in (100, 1366):
    text = lines[number - 1]
    print("ONTOLOGY_LINE", number, text)
    lowered = text.lower()
    for word in ("material", "energy", "particle", "coherent"):
        print("ONTOLOGY_LINE_HAS_WORD", number, word, word in lowered)
