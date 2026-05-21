#!/usr/bin/env python3
"""
Render a prompt template by substituting {KEY} placeholders with values.

Usage:
  render_prompt.py <template-path> KEY=value [KEY=value ...]

Reads the template file, performs literal string substitution of each `{KEY}`
occurrence with the corresponding value, prints the result to stdout.

Uses simple string.replace (no regex), so values with special characters
(spaces, ampersands, backslashes, etc.) are handled literally.

Unknown placeholders in the template are left untouched.
"""
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    template_path = Path(sys.argv[1])
    if not template_path.exists():
        print(f"error: template not found: {template_path}", file=sys.stderr)
        sys.exit(2)
    text = template_path.read_text()
    for arg in sys.argv[2:]:
        if '=' not in arg:
            print(f"warning: ignoring malformed arg (no '='): {arg}", file=sys.stderr)
            continue
        key, value = arg.split('=', 1)
        text = text.replace('{' + key + '}', value)
    sys.stdout.write(text)


if __name__ == '__main__':
    main()
