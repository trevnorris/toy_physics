#!/usr/bin/env python3
"""Detect drift in definitions intentionally repeated by the two ontology docs."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ONTOLOGY = ROOT / "docs" / "toy_model_ontology_summary.md"
COMPANION = ROOT / "docs" / "opposite_orientation_throat_coincidence.md"


def compact(path: Path) -> str:
    return "".join(path.read_text(encoding="utf-8").split())


def require(label: str, text: str, snippets: tuple[str, ...], failures: list[str]) -> None:
    missing = [snippet for snippet in snippets if snippet not in text]
    if missing:
        failures.append(f"{label}: missing {', '.join(missing)}")


def main() -> int:
    ontology = compact(ONTOLOGY)
    companion = compact(COMPANION)
    failures: list[str] = []

    shared = (
        r"h_m=\frac{h_++h_-}{2}",
        r"h_t=\frac{h_+-h_-}{2}",
        r"h_+\to-h_-",
        r"h_-\to-h_+",
        r"\mathcalC_{\rmth}=\mathcalI_{\rminternal}\circ\mathcalR_w",
    )
    require("canonical ontology shared definitions", ontology, shared, failures)
    require("opposite-orientation shared definitions", companion, shared, failures)

    require(
        "canonical flux identity",
        ontology,
        (r"\widetildeQ_{n,a}^{(w)}=s_aQ_{n,a}^{\rmnet}",),
        failures,
    )
    require(
        "companion flux identity",
        companion,
        (r"\widetildeQ_{n,\alpha}^{(w)}=sQ_{n,\alpha}^{\rmnet}",),
        failures,
    )
    require(
        "canonical graph-regime gap",
        ontology,
        (
            r"d_w^{\rmcore}=H_{\rmbr}+2h_t^{\rmtot}"
            r"-\ell_+^{\rmin}-\ell_-^{\rmin}",
        ),
        failures,
    )
    require(
        "companion graph-regime gap",
        companion,
        (
            r"d_w^{\rmcore}=H_0+2h_t^{\rmtot}"
            r"-\ell_+^{\rmin}-\ell_-^{\rmin}",
        ),
        failures,
    )

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1

    print(
        "PASS: shared h_m/h_t definitions, two-minus reflection map, "
        "flux identity, C_th map, and graph-regime core gap are synchronized."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
