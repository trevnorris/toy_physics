from native import PACKET, records
import sys
sys.path.insert(0,str(PACKET/'scripts'))
import sympy as s
import S15_reading_comparator as c
import S15_reading_context_checks as k
from collections import defaultdict

source_tags={
 'py':['family/changing_function_heads','family/held','family/jet_notation','family/centre_velocity_changes','family/centre_substitution','profile_family/separation_chart','profile_family/chart_condition','profile_family/substitution','family/response_dimension_signatures'],
 'wl':['family/field_change_notation','family/changing_heads','family/held','family/velocity_changes','family/separation_vector','profile_family/separation_chart','profile_family/domain','profile_family/definitions','family/dimension_signatures_including_changes']}
extra={
 'py':['squared_distance','change_of_variable','distance_endpoints','power_primitive_operand','power_primitive'],
 'wl':['distance_substitution','power_primitive_operand','power_primitive_branches']}
need={e+'_reading':set(tags)|{'item3_kernel/'+z for z in extra[e]} for e,tags in source_tags.items()}
raw,hashes=records(need)
base={(e,tag):v for (e,tag),(_,v) in raw.items()}
print('verified full/prefix bytes and sha256',hashes)
contexts=[{tag:[(e+'_reading',tag)] for tag in source_tags[e]} for e in ('py','wl')]
def family(data,field):
    out=[r for r in c.family_context(*contexts,data.__getitem__) if r['field']==field]
    assert len(out)==1
    return out[0]
def mutate(key,before,after,field):
    copy=base.copy()
    assert copy[key].count(before)==1
    copy[key]=copy[key].replace(before,after)
    a,b=family(base,field),family(copy,field)
    print('one-sided',key, 'offset',raw[key][0], 'edit',repr(before),'->',repr(after))
    print('result',a['status'],'->',b['status'],'same compared values?',(a['normalized_py'],a['normalized_wl'])==(b['normalized_py'],b['normalized_wl']))
    return copy
domain=mutate(('wl_reading','profile_family/domain'),'Reals','Integers','profile_family/chart_condition')
print('frozen and mutated Wolfram domains',base['wl_reading','profile_family/domain'],domain['wl_reading','profile_family/domain'])
print('admissible generic exponents',c.algebra.GENERIC_EXP[0],next(c.algebra.exponent_samples(0,()))[1])
binding=mutate(('py_reading','family/centre_substitution'),'X_1_1(t):','X_2_1(t):','family/centre_substitution')
print('frozen and corrupted SymPy target',base['py_reading','family/centre_substitution'].split(':')[0],binding['py_reading','family/centre_substitution'].split(':')[0])
def kernels(data):
    # Run production kernel_checks unmodified; stop at the first unused log record.
    class Cat(k.Catalogue):
        def __init__(self):
            self.rows=[]
        def get(self,e,tag,index=0):
            key=(e+'_reading',tag)
            if tag=='item3_kernel/log_primitive_operand': raise StopIteration
            return key,data[key]
    cat=Cat()
    try:k.kernel_checks(cat)
    except StopIteration: pass
    return cat.rows
endpoint=('py_reading','item3_kernel/distance_endpoints')
bad=base.copy()
bad[endpoint]=base[endpoint].replace('-integration_radius + kernel_distance','integration_radius - kernel_distance')
assert base[endpoint]!=bad[endpoint]
for title,data in [('accepted',base),('one-sided sign',bad)]:
    print(title,endpoint,raw[endpoint][0],'native',data[endpoint], 'r=2 d=5',s.sympify(data[endpoint])[0].subs({'integration_radius':2,'kernel_distance':5}))
    r=next(row for row in kernels(data) if row['field'].startswith('positive endpoint 0'))
    print(r['field'],r['translated_py'],r['translated_wl'], r['as_computed'],r['status'])
power=('wl_reading','item3_kernel/power_primitive_branches')
change=base.copy()
change[power]=change[power].replace('Log[squared]','-Log[squared]')
for title,data in [('accepted',base),('one-sided branch sign',change)]:
    for row in kernels(data):
        if row['field'].startswith('power primitive, exponent'):
            print(title,row['field'],row['translated_py'],row['translated_wl'],row['as_computed'],row['status'])
