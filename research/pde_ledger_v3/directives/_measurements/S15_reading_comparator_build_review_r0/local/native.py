"""Read frozen bounded native streams by line start, retaining original offset and hash."""
import hashlib
import json
import re
from pathlib import Path

LOCAL=Path(__file__).resolve().parent
PACKET=LOCAL.parent/'packet'/'research/pde_ledger_v3'
ROOT=next(p for p in LOCAL.parents if (p/'.git').exists())
SPEC=json.loads((PACKET/'directives/_measurements/S15_reading_comparator_run_manifest.json').read_text())['repository_inputs']
PREFIX=re.compile(rb'^([A-Za-z][a-zA-Z0-9_{},+ /-]*): ')

def records(needed):
    data={}
    manifests={}
    for stream,source in SPEC.items():
        want=needed.get(stream,set())
        with (ROOT/source['path']).open('rb') as f:
            h=hashlib.sha256()
            remaining=source['bytes']
            current=None
            while remaining:
                offset=f.tell()
                line=f.readline(remaining)
                if not line: raise ValueError((stream,'short input'))
                remaining-=len(line)
                h.update(line)
                match=PREFIX.match(line)
                if match:
                    name=match[1].decode()
                    current=(stream,name) if name in want else None
                    if current:
                        if current in data: raise ValueError(('duplicate',current))
                        data[current]=[offset+match.end(),line[match.end():]]
                elif current:
                    data[current][1]+=line
            assert h.hexdigest()==source['sha256']
            manifests[stream]=(f.tell(),h.hexdigest())
    return {key:(a,b.rstrip(b'\n').decode()) for key,(a,b) in data.items()},manifests
