#!/usr/bin/env python3
# MECHANICAL fact-lookup: parse the reconcile .out, tally the instrument's
# emitted nonzero_modular_numerator booleans per object. No CAS, no derived
# predicate — the instrument computed the flags via PIT; this only counts them.
import json, glob, collections

KEY = ['CARRIER_BRIDGE_RESIDUAL','SOURCE_BRIDGE_RESIDUAL',
       'CARRIER_CHANNEL','SOURCE_CHANNEL','CROSS_CHANNEL',
       'R_N','SPLIT_CHECK','SPLIT_SUM','EULERIAN_OPERAND','MATERIAL_OPERAND']

def blockkey(col):
    # col = [block, slot, face, [g1,g2]] (or similar); summarize by block only
    try:
        b = col[0] if isinstance(col,list) else col
        return str(b)
    except Exception:
        return str(col)

for path in sorted(glob.glob('/tmp/S11c_c2_N6_reconcile_sympy.*.out')):
    case = path.split('sympy.')[1].rsplit('.out',1)[0]
    per = {}
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or not line.startswith('{'): continue
            try: o=json.loads(line)
            except Exception: continue
            name=o.get('object','')
            data=o.get('data',{})
            if not isinstance(data,dict): continue
            nz=data.get('nonzero_modular_numerator')
            cols=data.get('columns')
            if nz is None or cols is None: continue
            short=name.replace('S11CC2_N6RC_','')
            probe=o.get('probe','')
            tag=short+('' if not probe else f'[{probe}]')
            total=len(nz); ncount=sum(1 for b in nz if b)
            # nonzero block/grade breakdown
            blocks=collections.Counter()
            for c,b in zip(cols,nz):
                if b: blocks[blockkey(c)]+=1
            per.setdefault(tag,[0,0,collections.Counter()])
            per[tag][0]+=total; per[tag][1]+=ncount
            per[tag][2].update(blocks)
    print(f'\n===== CASE {case} =====')
    for tag in sorted(per):
        total,ncount,blocks=per[tag]
        short=tag.split('[')[0]
        star = ' <<<' if short in KEY else ''
        bl=' '.join(f'{k}:{v}' for k,v in sorted(blocks.items())) if ncount else ''
        print(f'  {tag:42s} nonzero {ncount:3d}/{total:<4d} {bl}{star}')
