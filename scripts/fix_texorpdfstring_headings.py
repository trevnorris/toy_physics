#!/usr/bin/env python3
"""Wrap raw math in section-like headings with \\texorpdfstring.

This fixes headings like:
    \\section{Something with $\\alpha^2$ and \\(w\\)}

into:
    \\section{Something with \\texorpdfstring{$\\alpha^2$}{alpha squared}
    and \\texorpdfstring{\\(w\\)}{w}}

Only inline math fragments inside section-like headings are wrapped. The whole
heading is not rewritten.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


HEADING_RE = re.compile(r"\\(?:section|subsection|subsubsection|paragraph)\*?\{")


MATH_COMMAND_MAP: dict[str, str] = {
    r"\alpha": "alpha",
    r"\beta": "beta",
    r"\gamma": "gamma",
    r"\delta": "delta",
    r"\epsilon": "epsilon",
    r"\varepsilon": "epsilon",
    r"\zeta": "zeta",
    r"\eta": "eta",
    r"\theta": "theta",
    r"\vartheta": "theta",
    r"\iota": "iota",
    r"\kappa": "kappa",
    r"\lambda": "lambda",
    r"\mu": "mu",
    r"\nu": "nu",
    r"\xi": "xi",
    r"\pi": "pi",
    r"\rho": "rho",
    r"\sigma": "sigma",
    r"\tau": "tau",
    r"\phi": "phi",
    r"\varphi": "phi",
    r"\chi": "chi",
    r"\psi": "psi",
    r"\omega": "omega",
    r"\Gamma": "Gamma",
    r"\Delta": "Delta",
    r"\Theta": "Theta",
    r"\Lambda": "Lambda",
    r"\Xi": "Xi",
    r"\Pi": "Pi",
    r"\Sigma": "Sigma",
    r"\Phi": "Phi",
    r"\Psi": "Psi",
    r"\Omega": "Omega",
    r"\ell": "ell",
    r"\to": " -> ",
    r"\rightarrow": " -> ",
    r"\leftarrow": " <- ",
    r"\leftrightarrow": " <-> ",
    r"\times": " x ",
    r"\cdot": " dot ",
    r"\sim": " ~ ",
    r"\simeq": " approx ",
    r"\approx": " approx ",
    r"\oplus": " + ",
    r"\pm": " +/- ",
    r"\mp": " -/+ ",
    r"\neq": " != ",
    r"\geq": " >= ",
    r"\leq": " <= ",
    r"\infty": "infinity",
    r"\partial": "partial",
    r"\nabla": "nabla",
    r"\boldsymbol": "",
    r"\mathbf": "",
    r"\mathrm": "",
    r"\mathcal": "",
    r"\mathsf": "",
    r"\mathfrak": "",
    r"\operatorname": "",
    r"\left": "",
    r"\right": "",
    r"\,": " ",
    r"\!": "",
}


def find_matching_brace(text: str, open_index: int) -> int:
    if open_index >= len(text) or text[open_index] != "{":
        raise ValueError("Expected opening brace")
    depth = 1
    i = open_index + 1
    while i < len(text):
        ch = text[i]
        if ch == "\\":
            i += 2
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("Unbalanced braces")


def escape_fallback(text: str) -> str:
    return (
        text.replace("\\", "")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def simplify_superscripts(text: str) -> str:
    text = re.sub(r"\^\{2\}|\^2\b", " squared", text)
    text = re.sub(r"\^\{3\}|\^3\b", " cubed", text)
    text = re.sub(r"\^\{([^}]*)\}", r" to the \1", text)
    text = re.sub(r"\^([A-Za-z0-9])", r" to the \1", text)
    return text


def simplify_subscripts(text: str) -> str:
    text = re.sub(r"_\{([^}]*)\}", r" \1", text)
    text = re.sub(r"_([A-Za-z0-9])", r" \1", text)
    return text


def math_to_plain(math: str) -> str:
    plain = math.replace("~", " ")
    plain = simplify_superscripts(plain)
    plain = simplify_subscripts(plain)
    for old, new in MATH_COMMAND_MAP.items():
        plain = plain.replace(old, new)
    plain = plain.replace("{", "").replace("}", "")
    plain = plain.replace("$", "")
    plain = re.sub(r"\s+", " ", plain).strip()
    return escape_fallback(plain)


def wrap_inline_math(fragment: str, open_delim: str, close_delim: str) -> str:
    inner = fragment[len(open_delim) : len(fragment) - len(close_delim)]
    return rf"\texorpdfstring{{{fragment}}}{{{math_to_plain(inner)}}}"


def extract_texorpdfstring(text: str, start: int) -> tuple[str, int]:
    token = r"\texorpdfstring"
    if not text.startswith(token, start):
        raise ValueError("Expected texorpdfstring")
    i = start + len(token)
    if i >= len(text) or text[i] != "{":
        raise ValueError("Malformed texorpdfstring")
    end_first = find_matching_brace(text, i)
    j = end_first + 1
    if j >= len(text) or text[j] != "{":
        raise ValueError("Malformed texorpdfstring")
    end_second = find_matching_brace(text, j)
    return text[start : end_second + 1], end_second + 1


def extract_dollar_math(text: str, start: int) -> tuple[str, int] | None:
    if text[start] != "$":
        return None
    double = start + 1 < len(text) and text[start + 1] == "$"
    delim = "$$" if double else "$"
    i = start + len(delim)
    while i < len(text):
        if text[i] == "\\":
            i += 2
            continue
        if text.startswith(delim, i):
            end = i + len(delim)
            return text[start:end], end
        i += 1
    return None


def extract_paren_math(text: str, start: int) -> tuple[str, int] | None:
    if not text.startswith(r"\(", start):
        return None
    i = start + 2
    while i < len(text):
        if text[i] == "\\":
            if text.startswith(r"\)", i):
                return text[start : i + 2], i + 2
            i += 2
            continue
        i += 1
    return None


def convert_heading_content(content: str) -> tuple[str, int]:
    out: list[str] = []
    changes = 0
    i = 0
    while i < len(content):
        if content.startswith(r"\texorpdfstring", i):
            macro, i = extract_texorpdfstring(content, i)
            out.append(macro)
            continue

        paren = extract_paren_math(content, i)
        if paren is not None:
            fragment, i = paren
            out.append(wrap_inline_math(fragment, r"\(", r"\)"))
            changes += 1
            continue

        if content[i] == "$":
            dollar = extract_dollar_math(content, i)
            if dollar is not None:
                fragment, i = dollar
                open_delim = "$$" if fragment.startswith("$$") else "$"
                out.append(wrap_inline_math(fragment, open_delim, open_delim))
                changes += 1
                continue

        out.append(content[i])
        i += 1

    return "".join(out), changes


def rewrite_tex(text: str) -> tuple[str, int]:
    out: list[str] = []
    changes = 0
    pos = 0

    while True:
        match = HEADING_RE.search(text, pos)
        if match is None:
            out.append(text[pos:])
            break

        out.append(text[pos : match.end()])
        open_brace = match.end() - 1
        close_brace = find_matching_brace(text, open_brace)
        content = text[open_brace + 1 : close_brace]
        rewritten, delta = convert_heading_content(content)
        out.append(rewritten)
        out.append("}")
        changes += delta
        pos = close_brace + 1

    return "".join(out), changes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Wrap raw math in section-like headings with \\texorpdfstring."
    )
    parser.add_argument("path", help="Path to the .tex file to update.")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report whether changes would be made without writing the file.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the rewritten file to stdout instead of writing in place.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    path = Path(args.path)
    if not path.exists():
        print(f"error: file not found: {path}", file=sys.stderr)
        return 2

    original = path.read_text()
    rewritten, changes = rewrite_tex(original)

    if args.check:
        print(f"{path}: heading-math replacements needed={changes}")
        return 1 if changes else 0

    if args.dry_run:
        sys.stdout.write(rewritten)
        return 0

    if rewritten != original:
        path.write_text(rewritten)
    print(f"{path}: applied {changes} heading-math replacements")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
