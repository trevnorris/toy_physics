"""Compare the scoped S15 reading outputs. Run through s11c_guarded_run.py.

The operand catalogue preserves native data. Every comparison references its
operands by stable record ID; gzip files contain the results for one Q each.
"""
import argparse
import ast
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import re
import sys
import time

import sympy as s

import S15_reading_input_inventory as inputs
import S15_reading_math as algebra
from S15_cross_engine_comparator import GROUPS, PIECES, InputError


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def emit(stream,data):
    stream.write(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')


def save(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(tmp,'wt') as f:
        json.dump(data,f,separators=(',',':'),allow_nan=False)
    tmp.replace(path)


def required_names():
    result=set()
    groups=['+'.join(g) for g in GROUPS]
    for component in ('component_1','component_2','component_3'):
        for surface in ('inner','outer'):
            for piece in PIECES:
                root=f'{surface}/{piece}/{component}'
                result.add(root)
                if piece=='pressure':
                    result.update((root+'/first_order_choice',root+'/whole_minus_expansion'))
                    result.update(root+'/expansion_labels_'+g for g in groups)
                elif piece!='other':
                    result.update(root+'/labels_'+g for g in groups)
        region='region/'+component+'/'
        result.update(region+field for field in ('momentum','momentum_rate','integrated_force'))
        result.update(region+field+g for field in ('momentum_labels_','rate_labels_') for g in groups)
        result.update(f'{loc}/sum_{kind}/{component}' for loc in ('inner','outer','both_surfaces')
                      for kind in ('rest','transport','crossing'))
    return result


def object_key(name):
    return re.sub(r'labels_\{((?:drain1|drain2|remainder)(?:, (?:drain1|drain2|remainder))*)\}',
                  lambda m:'labels_'+'+'.join(m[1].split(', ')),name)


def read_inputs():
    streams={}
    digests={}
    counts={}
    manifest_names={}
    # Read from the final accepted manifests, including the verified base prefix.
    for engine,side in (('py','sympy'),('wl','wolfram')):
        manifest=json.loads(inputs.accepted(inputs.MEASUREMENTS/f'S15_separation_departure_{side}_run_manifest.json',inputs.READING))
        reading=inputs.ROOT/manifest['transcript']
        original=reading.with_name('S15_two_drain_region.out')
        record=inputs.accepted(inputs.MEASUREMENTS/'S15_build_run_record.md',inputs.BASE).decode()
        section=record.split('## '+manifest['engine']+'\n',1)[1].split('\n## ',1)[0]
        old_hash,=re.findall(r'- Transcript SHA256: `([0-9a-f]{64})`',section)
        for kind,path in [('original',original),('reading',reading)]:
            prefix=engine=='wl' and kind=='reading'
            tags,names,digest,size=inputs.record_index(path,'END_OF_BASE_PASS' if prefix else 'END_OF_TRANSCRIPT')
            if prefix:
                ref=manifest['final_transcript_production']['verification']
                assert (digest,size)==(ref['prefix_sha256'],ref['prefix_bytes']),path
            else:
                assert digest==(old_hash if kind=='original' else manifest['transcript_sha256']),path
            stream=engine+'_'+kind
            streams[stream]=(path,tags)
            counts[stream]=sum(map(len,tags.values()))
            digests[stream]={'path':str(path.relative_to(inputs.ROOT)),'sha256':digest,'bytes':size}
            if kind=='reading':
                manifest_names[engine]=names
    return streams,manifest_names,digests,counts


STANDARD={'Q':1,'separation_derivative':2,'E_geo':2,'E_resp':2,'departure_existence_condition':2,
          'departure_partition/operand_a':2,'departure_partition/operand_b':2,'departure_partition/residual':2,
          'item3/family_integral':3,'item3/existence_condition':3,'item3/object':3,
          'item3/Q_closed_operand':3,'item3/status':3,'item3/exponent_branch_inventory':3}
for name in ('velocity_part_1','velocity_part_2'):
    for role in ('operand_a','operand_b','residual','object'):
        STANDARD[name+'/'+role]=5
for name in ('E_geo','E_resp','velocity_part_1','velocity_part_2'):
    fields=['dimension_operands','dimension_pair_residuals','dimension_reference_residuals']
    if name in ('E_geo','E_resp'): fields.append('dimension_reference_operand')
    for field in fields:
        STANDARD[name+'/'+field]=6


def input_plan(streams,manifests,folder):
    roster=required_names()
    names={key:{} for key in sorted(roster)}
    observed={}
    gaps={}
    for engine in ('py','wl'):
        raw=manifests[engine]
        ptags=streams[engine+'_reading'][1]
        printed={tag[8:-2] for tag in ptags if tag.startswith('reading/') and tag.endswith('/Q')}
        raw_keys=[object_key(n) for n in raw]
        printed_keys={object_key(n) for n in printed}
        gaps[engine]={'manifest_missing':sorted(roster-set(raw_keys)),'manifest_extra':sorted(set(raw_keys)-roster),
                      'stream_missing':sorted(roster-printed_keys),'stream_extra':sorted(printed_keys-roster),
                      'duplicate_manifest_keys':[k for k,v in Counter(raw_keys).items() if v>1]}
        if gaps[engine]['duplicate_manifest_keys']:
            raise InputError(f'ambiguous native manifest: {gaps[engine]}')
        known={}
        for native in sorted(set(raw)|printed,key=lambda x:(-len(x),x)):
            key=object_key(native)
            if key in names:
                names[key][engine]={'native':native,'fields':{},'original':{}}
            known[native]=key
        observed[engine]=known

    inventory=[]
    blobs={}
    locator={}
    global_tags=defaultdict(dict)
    with gzip.open(folder/'operands.jsonl.gz','wt',compresslevel=2) as catalog:
        for stream,(path,tags) in streams.items():
            engine,kind=stream.split('_')
            original_tags=streams[engine+'_original'][1]
            with path.open('rb') as f:
                for tag,entries in tags.items():
                    native=None
                    if tag.startswith('reading/'):
                        parent=tag[len('reading/'):]
                        while '/' in parent:
                            parent=parent.rsplit('/',1)[0]
                            if parent in observed[engine]:
                                native=parent
                                break
                    suffix=tag[9+len(native):] if native is not None else None
                    key=observed[engine].get(native)
                    for ordinal,span in enumerate(entries):
                        rid=f'{stream}:{span["offset"]}'
                        role='source'
                        if ('/item4/' in tag or tag.startswith(('item4_kernel/','ITEM4_','END_OF_ITEM4'))):
                            role='excluded_item4'
                        elif tag=='WORKING_ON':
                            role='progress'
                        data=inputs.payload(f,span)
                        if tag=='EXACT_OBJECT_GAP' and engine=='py' and ast.literal_eval(data.decode())[1]=='item4':
                            role='excluded_item4'
                        if role=='source' and kind=='reading':
                            if native is not None:
                                role='scoped' if key in names else 'extra_object'
                                if key in names:
                                    names[key][engine]['fields'].setdefault(suffix,[]).append(rid)
                            elif tag in original_tags:
                                role='original_context'
                            elif tag.startswith(('item3_kernel/','family/','profile_family/')) or tag in (
                                    'reading/dimension_argument_audit','Q_manifest','END_OF_BASE_PASS','END_OF_TRANSCRIPT'):
                                role='reading_context'
                            else:
                                role='unresolved_context'
                            if native is None: global_tags[stream].setdefault(tag,[]).append(rid)
                        rec={'id':rid,'tag':tag,'occurrence':ordinal,'span':span,'role':role}
                        if role not in ('excluded_item4','progress'):
                            digest=hashlib.sha256(data).hexdigest()
                            rec['payload']=digest
                            if digest not in blobs:
                                emit(catalog,{'sha256':digest,'native_operand':data.decode()})
                                blobs[digest]=True
                            locator[rid]=(stream,span,digest)
                        inventory.append(rec)
            if kind=='original':
                for obj in names.values():
                    if engine not in obj:continue
                    row=obj[engine]
                    native=row['native']
                    for source in inputs.sum_operands(native) or [native]:
                        entries=tags.get(source,[])
                        row['original'][source]=[f'{stream}:{e["offset"]}' for e in entries]
    save(folder/'coverage.json.gz',inventory)
    return names,gaps,locator,global_tags


PATHS={}
LOCATOR={}
OBJECTS={}
ROOT=None
BUILD_ID=None


def fraction_rows(text,engine):
    if engine=='wl':
        text=text.replace('{','[').replace('}',']')
    shape=ast.literal_eval(text)
    if not isinstance(shape,(tuple,list)):
        raise InputError('dimension collection required')
    if any(not isinstance(v,(tuple,list)) or len(v)!=3 or any(type(x)!=int for x in v) for v in shape):
        raise InputError('invalid dimension vector')
    return sorted(tuple(v) for v in shape)


def decoded_label(payload,engine):
    return json.loads(payload) if engine=='wl' and payload.startswith('"') else payload


class ObjectComparison:
    def __init__(self,handles):
        self.handles=handles
        self.parsed={}
        self.pairs={}
        self.rows=[]
        self.loci=()
        self.branch_issue=None

    def text(self,rid):
        stream,span,_=LOCATOR[rid]
        return inputs.payload(self.handles[stream],span).decode()

    def operand(self,rid):
        stream,_,digest=LOCATOR[rid]
        engine=stream.split('_')[0]
        key=engine,digest
        if key not in self.parsed:
            self.parsed[key]=algebra.parse(self.text(rid),engine)
        return self.parsed[key]

    def arithmetic(self,left,right,field,item,loci):
        # All arithmetic publications refer to raw catalogue entries, including identical bytes.
        same_engine=LOCATOR[left][0][:2]==LOCATOR[right[0]][0][:2]
        if len(right)==1 and same_engine and LOCATOR[left][2]==LOCATOR[right[0]][2]:
            return {'status':'exact_identity','residual':'0','method':'identical native expression on same engine'}
        a=self.operand(left)
        b=s.Add(*(self.operand(r) for r in right)) if len(right)>1 else self.operand(right[0])
        res=algebra.exact_residual(a,b)
        if res==0:
            return {'status':'exact_identity','residual':'0','method':'exact algebra on translated operands'}
        if res.is_number and not res.free_symbols and not res.has(s.Integral,s.Piecewise) and res.is_zero is False:
            return {'status':'exact_difference','residual':str(res),'method':'nonzero constant residual'}
        result=algebra.sample_values(a,b,loci)
        # The catalogue and the ordered arithmetic tree preserve the exact residual
        # losslessly, including signs, without printing megabytes of expanded jets twice.
        result['residual']={'op':'subtract','left':left,'right':{'op':'sum','operands':right}}

        result['method']='native arithmetic difference and domain samples'
        return result

    def compare(self,key,item,left,right,role):
        row={'field':key,'item':item,'left':left,'right':right,'role':role}
        if len(left)!=1 or not right or (role=='cross' and len(right)!=1):
            row.update(status='unresolved',reason='absent, extra, or ambiguous occurrence')
            self.rows.append(row)
            return
        locus=self.loci if item==3 and key in ('item3/object','item3/Q_closed_operand','item3/family_integral') else ()
        mode=('dimension' if item==6 else 'branch' if key.endswith('exponent_branch_inventory')
              else 'text' if key.endswith(('existence_condition','/status','/comparison')) else 'math')
        identity=(mode,tuple((LOCATOR[r][0][:2],LOCATOR[r][2]) for r in left+right),locus)
        if identity in self.pairs:
            index=self.pairs[identity]
            row.update(status=self.rows[index]['status'],same_comparison_as=index)
            self.rows.append(row)
            return
        try:
            if mode=='dimension':
                a=fraction_rows(self.text(left[0]),'py')
                b=fraction_rows(self.text(right[0]),'wl')
                row.update(status='exact_identity' if a==b else 'dimension_difference',basis=['Density','Length','Time'],
                           residual={'py_only':list((Counter(a)-Counter(b)).elements()),'wl_only':list((Counter(b)-Counter(a)).elements())})
            elif mode=='text':
                a=decoded_label(self.text(left[0]),'py')
                b=decoded_label(self.text(right[0]),'wl')
                row.update(status='text_identical' if a==b else 'text_difference',residual=[] if a==b else [a,b])
            elif mode=='branch':
                a=algebra.branch_loci(self.text(left[0]),'py')
                b=algebra.branch_loci(self.text(right[0]),'wl')
                row.update(status='same_exponent_loci' if a==b else 'different_exponent_loci',
                           residual={'py_only':[str(q) for q in set(a)-set(b)],'wl_only':[str(q) for q in set(b)-set(a)]},
                           boundary_normalization='Equality and inequality boundaries only; native branch selection retained in value comparisons.')
            else:
                row.update(self.arithmetic(left[0],right,key,item,locus))
                if self.branch_issue and locus==() and item==3:
                    row['branch_coverage_issue']=self.branch_issue
        except (InputError,algebra.Unsupported,ValueError,TypeError,RecursionError,NotImplementedError) as exc:
            row.update(status='unresolved',reason=f'{type(exc).__name__}: {exc}')
        self.pairs[identity]=len(self.rows)
        self.rows.append(row)

    def run(self,obj):
        py=obj.get('py',{'fields':{}})
        wl=obj.get('wl',{'fields':{}})
        fields=(py['fields'],wl['fields'])
        branch='item3/exponent_branch_inventory'
        try:
            self.loci=tuple(sorted(set().union(*(algebra.branch_loci(self.text(f[branch][0]),e)
                                                for f,e in zip(fields,('py','wl'))))))
        except (KeyError,InputError,algebra.Unsupported,ValueError,IndexError) as exc:
            self.branch_issue=str(exc)
        for engine,native in obj.items():
            refs=native['original']
            linked=all(len(v)==1 for v in refs.values())
            self.compare('original_Q',1,native['fields'].get('Q',[]),
                         [v[0] for v in refs.values()] if linked else [],engine+'_same_engine')
        # Every native field is bound, explicitly qualified, or listed as unresolved.
        required_fields=set(STANDARD)-{'item3/Q_closed_operand','velocity_part_1/object','velocity_part_2/object'}
        if any('Closed-form' in self.text(r) for f in fields for r in f.get('item3/status',[])):
            required_fields.add('item3/Q_closed_operand')
        for key in sorted(set(fields[0])|set(fields[1])|required_fields):
            left,right=[f.get(key,[]) for f in fields]
            if key in ('velocity_part_1/object','velocity_part_2/object',
                       'velocity_part_1/object/comparison','velocity_part_2/object/comparison'):
                base=key.split('/')[0]
                left=fields[0].get(key.replace('/object','/residual'),[])
            if key.endswith('/comparison') and key.startswith('velocity_part') and not right:
                other=key.replace('/residual/','/object/')
                if other in fields[1]:continue
            base=key.removesuffix('/comparison')
            item=STANDARD.get(base)
            if item is None or not left or not right:
                reason=('Wolfram-only closed-family parameter or divisor audit' if
                        any(q in key for q in ('exponent_divisor','retained_Appell','exceptional_status')) else
                        'native per-engine field without established counterpart')
                if not left and not right:
                    reason='Required native field missing from both streams.'
                self.rows.append(dict(field=key,item=item or (3 if 'item3/' in key else None),left=left,right=right,
                                      role='cross',status='unpaired',reason=reason))
            else:
                self.compare(key,item,left,right,'cross')
        return self.rows


def object_worker(job):
    idx,name=job
    with ExitStack() as stack:
        handles={k:stack.enter_context(p.open('rb')) for k,p in PATHS.items()}
        computation=ObjectComparison(handles)
        rows=computation.run(OBJECTS[name])
    output={'object':name,'index':idx,'build_id':BUILD_ID,'exponent_loci':[list(map(str,q)) for q in computation.loci],
            'locus_intersections':'Samples target each listed boundary individually. Other simultaneous degeneracies are unverified.',
            'rows':rows}
    save(ROOT/'objects'/f'{idx:03d}.json.gz',output)
    algebra.translate.cache_clear()
    algebra.normal.cache_clear()
    algebra.outer_polynomial.cache_clear()
    s.core.cache.clear_cache()
    return idx,len(rows)


def family_context(py,wl,get):
    """Source-defined correspondences for the family and the real separation chart."""
    ns=algebra.original
    def tree(tag,engine):
        records=(py if engine=='py' else wl)[tag]
        if len(records)!=1: raise InputError('ambiguous family context '+tag)
        raw=get(records[0])
        if engine=='py': return ns.py_parse(raw)
        return ns.wl_parse(algebra.WL_SCALARS.sub(lambda m:algebra.CHART[m[0].replace(' ','')],raw))
    def val(node,engine,scope=()):
        return algebra.translate(node,engine,scope)
    def named(node,engine):
        if node.head=='Symbol' and node.value in algebra.CHART.values(): return node.value
        key=ns.key_text(node)
        return algebra.NAMES.name(key,engine)
    def held(node,engine):
        if ns.call(node,'Matrix'):
            return tuple(val(row.args[0],engine) for row in node.args[1].args)
        if node.head=='Sequence':
            return tuple(val(q,engine) for q in node.args)
        return val(node,engine)
    def pairs(node,key,value):
        return dict((key(a),value(b)) for a,b in (p.args for p in node.args))
    def cv(n,e):return val(n,e)
    a='family/changing_function_heads'
    b='family/field_change_notation'
    heads=pairs(tree(a,'py'),lambda q:named(q,'py'),lambda q:named(q,'py'))
    wlheads=pairs(tree(b,'wl'),lambda q:named(q,'wl'),lambda q:named(q,'wl'))
    def row(tag,other,left,right,**kw):
        return dict(field=tag,left=py[tag],right=wl[other],status='same_values' if left==right else 'value_difference',
                    normalized_py=str(left),normalized_wl=str(right),**kw)
    yield row(a,b,heads,wlheads)
    yield row(a,'family/changing_heads',sorted(heads),sorted(named(q,'wl') for q in tree('family/changing_heads','wl').args))
    tag='family/held'
    yield row(tag,tag,tuple(held(q,'py') for q in tree(tag,'py').args),tuple(held(q,'wl') for q in tree(tag,'wl').args))
    tag='family/centre_velocity_changes';other='family/velocity_changes'
    jets=pairs(tree('family/jet_notation','py'),ns.key_text,lambda q:cv(q,'py'))
    yield row(tag,other,pairs(tree(tag,'py'),lambda q:jets[ns.key_text(q)],lambda q:cv(q,'py')),
              pairs(tree(other,'wl'),lambda q:cv(q,'wl'),lambda q:cv(q,'wl')),jet_definition=py['family/jet_notation'])
    tag='family/centre_substitution';other='family/separation_vector'
    ep=s.Symbol('family_epsilon')
    shifts=[val(q.args[1],'py',(('family_epsilon',ep),)).diff(ep) for q in tree(tag,'py').args]
    yield row(tag,other,shifts,[cv(q,'wl') for q in tree(other,'wl').args])
    tag='profile_family/separation_chart'
    yield row(tag,tag,*[pairs(tree(tag,e),lambda q:cv(q,e),lambda q:cv(q,e)) for e in ('py','wl')])
    tag='profile_family/chart_condition';other='profile_family/domain'
    wlconds=tree(other,'wl').args
    yield row(tag,other,tuple(cv(q,'py') for q in tree(tag,'py').args),tuple(cv(q,'wl') for q in wlconds[:2]),
              additional_wolfram_clause='Explicit real domain for distance, direction and four exponents; retained in right operand.')
    tag='profile_family/substitution';other='profile_family/definitions'
    def body(node):
        if not ns.call(node,'Lambda'): raise InputError('expected profile lambda')
        scope=tuple((q.value,v) for q,v in zip(node.args[1].args,s.symbols('sample_radius sample_time')))
        return val(node.args[2],'py',scope)
    ws=tuple(zip(('profileRadius','profileTime'),s.symbols('sample_radius sample_time')))
    yield row(tag,other,pairs(tree(tag,'py'),lambda q:named(q,'py'),body),
              pairs(tree(other,'wl'),lambda q:named(q,'wl'),lambda q:val(q,'wl',ws)))
    tag='family/response_dimension_signatures';other='family/dimension_signatures_including_changes'
    sigpy=pairs(tree(tag,'py'),lambda q:named(q,'py'),lambda q:val(q,'py')[0])
    sigwl=pairs(tree(other,'wl'),lambda q:named(q,'wl'),lambda q:val(q,'wl'))
    yield row(tag,other,sigpy,{k:sigwl.get(k) for k in sigpy},
              note='Result dimensions of field changes. Python prints argument signatures alongside them.')


def compare_context(global_tags):
    py=global_tags['py_reading'];wl=global_tags['wl_reading']
    result=[]
    # Context printed at base start predates the reading objects. Bind the native
    # domain relation and the label texts; the older accepted comparator owns its other results.
    with ExitStack() as stack:
        handles={k:stack.enter_context(p.open('rb')) for k,p in PATHS.items()}
        get=ObjectComparison(handles).text
        for key in ('PRESSURE_DOMAIN','REGULARITY','PREMISES','CHOICES','BOOKKEEPING','dimension_convention',
                    'family/choice','profile_family/chart_note','DOMAIN_FOR_EVERY_OBJECT','dimension_basis'):
            a,b=py[key][0],wl[key][0]
            av,bv=get(a),get(b)
            row=dict(field=key,left=a,right=b)
            if key=='dimension_basis':
                vals=ast.literal_eval(av)
                other=[s.strip() for s in bv[1:-1].split(',')]
                row['status']='same_basis' if tuple(algebra.NAMES.name(q,'wl') for q in other)==tuple(vals) else 'basis_difference'
                if row['status']!='same_basis':raise InputError('unverified basis for dimension comparisons')
            elif key=='DOMAIN_FOR_EVERY_OBJECT':
                ns=algebra.original
                left=ns.py_parse(av).args[0]
                right=ns.wl_parse(bv).args[0]
                # Native selection by source role. In SymPy the relation is first; in Wolfram too.
                aa=ns.math_node(left,algebra.NAMES,'py')
                bb=ns.math_node(right,algebra.NAMES,'wl')
                row.update(status='same_relation' if aa==bb else 'relation_comparison_unresolved',py=str(aa),wl=str(bb))
            else:
                row.update(status='text_identical' if decoded_label(bv,'wl')==av else 'text_difference')
            result.append(row)
        result.extend(family_context(py,wl,get))
        for tag,records in py.items():
            if tag=='reading/dimension_argument_audit':
                expr=ast.parse(get(records[0]),mode='eval').body
                for index,(key,value) in enumerate(zip(expr.keys,expr.values)):
                    # These are each SymPy argument/bound's expected and observed sets.
                    def vectors(node):
                        if not (isinstance(node,ast.Call) and getattr(node.func,'id',None)=='frozenset'):
                            raise InputError('unexpected dimension-audit representation')
                        return sorted(ast.literal_eval(node.args[0])) if node.args else []
                    x,y=map(vectors,value.elts)
                    result.append(dict(field=tag,ordinal=index,left=records[0],argument=ast.literal_eval(key),
                                       expected=x,observed=y,status='unpaired_argument_check',
                                       right=wl.get('field_argument_dimension_signatures',[]),
                                       reason='Wolfram supplies signatures; SymPy checks each observed argument or integral bound.'))
        accounted=set()
        for row in result:
            for field in ('left','right','jet_definition'):
                ids=row.get(field,[])
                accounted.update(ids if isinstance(ids,list) else [ids])
        # Original inputs predate the reading family. Native construction diagnostics
        # have engine-specific operand order and multiplicity; do not pair by ordinal.
        for e,native in [('py',py),('wl',wl)]:
            for tag,ids in native.items():
                ids=[rid for rid in ids if rid not in accounted]
                if not ids:continue
                if tag.startswith('item3_kernel/'):
                    reason='Internal derivation or branch diagnostic; no verified cross-engine pairing for these occurrences.'
                    status='unpaired_kernel_diagnostic'
                elif tag.startswith(('family/','profile_family/','reading/')) or tag=='EXACT_OBJECT_GAP':
                    reason='Additional native reading record; counterpart not established.'
                    status='unpaired_reading_context'
                elif tag in ('Q_manifest','END_OF_TRANSCRIPT','END_OF_BASE_PASS'):
                    reason='Manifest or byte boundary independently verified in configuration and roster.'
                    status='verified_inventory'
                else:
                    reason='Source context from the original build; preserved by exact stream offset.'
                    status='source_context'
                result.append(dict(field=tag,engine=e,records=ids,status=status,reason=reason))
    return result


def main():
    global ROOT,BUILD_ID,OBJECTS,LOCATOR,PATHS
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--workers',type=int,default=1)
    args=ap.parse_args()
    sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
    ROOT=args.output
    ROOT.mkdir(parents=True,exist_ok=True)
    (ROOT/'objects').mkdir(exist_ok=True)
    if not 1<=args.workers<=len(os.sched_getaffinity(0)):
        ap.error('worker count exceeds inherited affinity or is not positive')
    sources={f.name:sha(f) for f in (Path(__file__),Path(algebra.__file__),Path(inputs.__file__),Path(algebra.original.__file__))}
    streams,manifests,digests,counts=read_inputs()
    BUILD_ID=hashlib.sha256(json.dumps([sources,digests],sort_keys=True).encode()).hexdigest()
    config=ROOT/'configuration.json'
    if config.exists() and json.loads(config.read_text())['build_id']!=BUILD_ID:
        raise InputError('source/input fingerprint changed; select a fresh output directory')
    PATHS={key:path for key,(path,_) in streams.items()}
    OBJECTS,gaps,LOCATOR,globals_=input_plan(streams,manifests,ROOT)
    save(config,dict(build_id=BUILD_ID,sources=sources,inputs=digests,occurrences=counts,roster=gaps,
                     required_objects=sorted(OBJECTS),forms=algebra.FORMS,strata=algebra.STRATA,
                     points=[list(map(str,q)) for q in algebra.POINTS],quadrature_orders=algebra.RULES,
                     precision='IEEE float64; tensor Gauss-Legendre for chart integrals; normalized Gauss-Jacobi Euler integral for Appell F1; Li2 via scipy.special.spence.',
                     evaluator='SymPy difference with original S15 profile jets, NumPy/SciPy numerical fallback. Geometric direction and distance computed from the two live centres.'))
    save(ROOT/'context.json.gz',compare_context(globals_))
    jobs=[]
    for index,name in enumerate(sorted(OBJECTS)):
        path=ROOT/'objects'/f'{index:03d}.json.gz'
        if path.exists():
            with gzip.open(path,'rt') as f: prior=json.load(f)
            if (prior['build_id'],prior['object'])!=(BUILD_ID,name):raise InputError('checkpoint mismatch')
        else:
            jobs.append((index,name))
    # Largest source blocks start early; each object is published atomically.
    jobs.sort(key=lambda job:sum(LOCATOR[r][1]['length'] for entry in OBJECTS[job[1]].values()
                                  for field in ('E_geo','E_resp','item3/object') for r in entry['fields'].get(field,[])),reverse=True)
    completed=len(OBJECTS)-len(jobs)
    with ProcessPoolExecutor(args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
        futures={pool.submit(object_worker,job):job for job in jobs}
        for future in as_completed(futures):
            index,n=future.result()
            completed+=1
            print(f'{completed}/{len(OBJECTS)} object {index} {futures[future][1]} records {n}',flush=True)
    summarize(ROOT)


def summarize(root):
    counts=defaultdict(Counter)
    byfield=defaultdict(Counter)
    exceptions=[]
    n=0
    for path in sorted((root/'objects').glob('*.json.gz')):
        with gzip.open(path,'rt') as f: obj=json.load(f)
        n+=1
        for row in obj['rows']:
            k=row['role']+' / '+row['field']
            counts[str(row['item'])][row['status']]+=1
            byfield[k][row['status']]+=1
            if row['status'] in ('exact_difference','sample_difference','unresolved','sample_coverage_incomplete','dimension_difference','different_exponent_loci'):
                exceptions.append({'object':obj['object'],'field':k,'status':row['status'],
                                   'reason':row.get('reason'),'sample_counts':row.get('sample_counts')})
    save(root/'summary.json',dict(objects=n,by_item=counts,by_field=byfield,exceptions=exceptions))


if __name__=='__main__':
    main()
