#!/usr/bin/env python3
"""Re-run every appendix command in the survey and diff its stdout against the recorded literal output."""
import re, subprocess, sys
src = open('/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md').read()
app = src[src.index('## Appendix A.'):]
blocks = re.split(r'<a id="([a-z0-9]+)"></a>', app)
# blocks: [pre, id1, body1, id2, body2, ...]
for i in range(1, len(blocks), 2):
    aid, body = blocks[i], blocks[i+1]
    m = re.search(r'``````bash\n(.*?)\n``````', body, re.S)
    o = re.search(r'Literal output:\s*\n``````text\n(.*?)``````', body, re.S)
    if not m or not o:
        print(f'{aid}: could not parse command/output'); continue
    cmd, rec = m.group(1), o.group(1)
    r = subprocess.run(['bash','-c',cmd], cwd='/var/projects/toy_physics', capture_output=True, text=True)
    out = r.stdout
    same = out == rec or out.rstrip('\n') == rec.rstrip('\n')
    print(f'{aid}: exit={r.returncode} rec_len={len(rec)} new_len={len(out)} identical={same}')
    if not same:
        import difflib
        d = list(difflib.unified_diff(rec.splitlines(), out.splitlines(), 'recorded', 'rerun', lineterm='', n=0))
        print('\n'.join(d[:20]))
