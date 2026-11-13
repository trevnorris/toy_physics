"""One-sided mutations at the normal per-object comparison dispatch."""
from native import records,PACKET
import gzip,json,io,hashlib,sys,tarfile
from collections import Counter
sys.path.insert(0,str(PACKET/'scripts'))
import S15_reading_comparator as c
import S15_reading_math as m
import S15_reading_refinement as ref
import sympy as s

archive=PACKET/'scripts/out/S15_reading_comparator_details.out'
with tarfile.open(archive) as tf:
    def load(n): return json.load(gzip.GzipFile(fileobj=tf.extractfile(n))) if n.endswith('.gz') else json.load(tf.extractfile(n))
    data=load('configuration.json'); coverage=load('coverage.json.gz')
    stored={r['id']:r for r in coverage}
    num=80
    obj=load(f'objects/{num:03d}.json.gz')
    print('component',obj['object'])
    needs={e+'_reading':set() for e in ('py','wl')}
    dim='E_geo/dimension_reference_residuals'
    selected=['velocity_part_1/operand_b','velocity_part_2/operand_b','velocity_part_2/residual',dim,'departure_existence_condition']
    info={r['field']:r for r in obj['rows'] if r['role']=='cross' and r['field'] in selected}
    for r in info.values():
        for refid in (*r['left'],*r['right']):
            e=refid.split(':')[0]; needs[e].add(stored[refid]['tag'])
raw,_=records(needs)
def fetch(rid):
    e=rid.split(':')[0];off,txt=raw[e,stored[rid]['tag']]
    assert int(rid.split(':')[1])==off
    return txt
def call(field,changed=None):
    row=info[field]; refs=row['left']+row['right']; buffers={}
    for i,rid in enumerate(refs):
        stream=rid.split(':')[0]; txt=fetch(rid)
        if i==1 and changed is not None:txt=changed
        body=(txt+'\n').encode(); buffers[stream]=io.BytesIO(body)
        c.LOCATOR[rid]=(stream,{'offset':0,'length':len(body)},hashlib.sha256(body[:-1]).hexdigest())
    runner=c.ObjectComparison(buffers)
    runner.compare(field,row['item'],row['left'],row['right'],'cross')
    return runner.rows[-1]
def comparefield(f,bad):
    old,new=call(f),call(f,bad)
    print(f,'native offsets',info[f]['left']+info[f]['right'],'statuses',old['status'],'->',new['status'])
    if new.get('samples'):
        sample=new['samples'][0]; print('new first-case values and residual',sample['values'][-1],sample['residual'],'case statuses',new['sample_counts'])
    else:print('residual',old['residual'],'->',new['residual'])
w1=fetch(info['velocity_part_1/operand_b']['right'][0]);w2=fetch(info['velocity_part_2/operand_b']['right'][0])
print('untouched WL selected-2 without native value',w2)
comparefield('velocity_part_2/operand_b',w1)
comparefield(dim,'{}')
state=call('departure_existence_condition')
print('sufficient-condition native text typed outcome',state['status'],'literal',state['residual'])
t=s.Symbol('t');x=s.Function('X_1_1')(t);y=s.Function('X_2_1')(t)
# Keep the original's stationary test for the reported residual=128. Exercise its
# jet selection when a Hessian/time acceleration is present in a changed operand.
base=s.diff(x,t)+s.diff(y,t)+s.diff(y,t,2)
after=ref.stationary_case(base,['', 'X2_and_radii_velocity_zero',0,''])
print('target/radii stationary exact replacement applied to changed jet sum',base,'->',after)
