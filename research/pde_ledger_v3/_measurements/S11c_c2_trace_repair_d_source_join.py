#!/usr/bin/env python3
"""Join scoped and full native cache operands without altering either producer."""
import argparse,hashlib,json,pickle,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'scripts'))
import sympy as sp
parser=argparse.ArgumentParser()
parser.add_argument('--run-root',type=Path,default=Path('/tmp/s11c-trace-repair-20260913'))
parser.add_argument('--report',type=Path,default=root/'_measurements/S11c_c2_trace_repair_d_source_joins.json')
args=parser.parse_args();base=args.run_root
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
manifest_paths={name:base/name/'manifest.json' for name in ('d_ends','d_full')}
manifest={name:json.loads(p.read_text()) for name,p in manifest_paths.items()}
report={'scope':'Structural joins of independently executed calls to the same native constructors; provenance checks, not an independent physics derivation.',
 'instrumentSha256':sha(Path(__file__)),
 'manifests':{n:{'path':str(p),'sha256':sha(p)} for n,p in manifest_paths.items()},'ends':{}}
for name,m in manifest.items():
 if m['exit_code']!=0 or m['source_hashes_before']!=m['source_hashes_after']:raise ValueError(('producer incomplete/changed',name))
 for path,pin in m['source_hashes_after'].items():
  if sha(root/path)!=pin or sha(base/name/'source'/path)!=pin:raise ValueError(('source pin',name,path))
common_sources=manifest['d_ends']['source_hashes_before'].keys()&manifest['d_full']['source_hashes_before'].keys()
if any(manifest['d_ends']['source_hashes_before'][p]!=manifest['d_full']['source_hashes_before'][p] for p in common_sources):raise ValueError('scoped/full source inputs differ')
for end in ('REFERENCE','LEFT','RIGHT'):
 payloads={};records={}
 for name,m in manifest.items():
  relative='symbols/'+end+'_LAB_HELD_RHO4_CONSTANT.pickle';p=base/name/relative
  if sha(p)!=m['artifacts'][relative]['sha256']:raise ValueError(('cache pin',name,end))
  payloads[name]=pickle.loads(p.read_bytes());records[name]={'path':str(p),'sha256':sha(p)}
 a,b=payloads['d_ends'],payloads['d_full']
 joins={label:a[i]==b[i] for i,label in ((0,'weakPencil'),(1,'curl'),(2,'weakUnits'),(4,'strongPencil'),(5,'energyBasis'))}
 common=a[3].keys()&b[3].keys();units={k:a[3][k]==b[3][k] for k in common}
 report['ends'][end]={'cachePins':records,'structuralJoins':joins,'sharedDimensionBindings':len(common),
  'differentSharedDimensionBindings':[sp.srepr(k) for k,v in units.items() if not v]}
path=args.report
path.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
if any(not all(v['structuralJoins'].values()) or v['differentSharedDimensionBindings'] for v in report['ends'].values()):raise SystemExit(1)
