"""Independent reference for the two audited S15 publication fields.

Run under the repository guard. No comparator modules or candidate outputs are read.
"""
import ast
import gzip
import hashlib
import itertools
import json
import re
import tarfile
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent / 'packet'
P = Path('research/pde_ledger_v3')
ARCHIVE = PACKET / P / 'scripts/out/S15_reading_comparator_details.out'
COVERAGE = PACKET / P / 'directives/_measurements/S15_reading_comparator_coverage.json'
SUMMARY = PACKET / P / 'scripts/out/S15_reading_comparator.out'
MANIFEST = PACKET / P / 'directives/_measurements/S15_reading_comparator_run_manifest.json'
VARS = ['a_1', 'a_2', 'b_1', 'b_2']
ARCHIVED_BYTES = {
    'py': (192065031, '10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6'),
    'wl': (333550018, '75622128d04da2d934c2a52468505b4438fff9f17a902d7a6cd35993d769b6c4'),
}


def member(archive, key):
    data = archive.extractfile(key).read()
    return json.loads(gzip.decompress(data) if key.endswith('.gz') else data)


def expression(text, point):
    # The two native inventories emit affine arithmetic with rational constants.
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return Q(node.value)
        if isinstance(node, ast.Name):
            return point[VARS.index(node.id)]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Div):
                return a / b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
        raise ValueError(ast.dump(node))
    return visit(ast.parse(text.strip(), mode='eval').body)


@lru_cache(None)
def plane(lhs, rhs):
    zero = (Q(0),) * 4
    const = expression(lhs, zero) - expression(rhs, zero)
    coeff = []
    for axis in range(4):
        unit = tuple(Q(i == axis) for i in range(4))
        coeff.append(expression(lhs, unit) - expression(rhs, unit) - const)
    row = tuple(coeff + [const])
    first = next(x for x in coeff if x)
    return tuple(c / first for c in row)


@lru_cache(None)
def inventory(engine, raw):
    """Unpack all atomic predicates, whether single or part of a conjunction."""
    planes = set()
    if engine == 'py':
        for term in ast.literal_eval(raw):
            if term == 'True':
                continue
            relation = ast.parse(term, mode='eval').body
            assert relation.func.id == 'Ne' and len(relation.args) == 2
            planes.add(plane(*(ast.unparse(a) for a in relation.args)))
    else:
        assert raw[0] == '{' and raw[-1] == '}'
        if raw != '{}':
            for clause in raw[1:-1].split(', '):
                for eq in clause.split(' && '):
                    converted = re.sub(r'([ab])\1\[([12])\]', r'\1_\2', eq)
                    left, right = converted.split(' == ')
                    planes.add(plane(left, right))
    return tuple(sorted(planes))


def read_stream(engine):
    path = P / ('scripts' if engine == 'py' else 'mathematica') / 'out/S15_separation_departure.out'
    size, digest = ARCHIVED_BYTES[engine]
    h = hashlib.sha256()
    branches, operands, context = {}, {}, {}
    used = 0
    with path.open('rb') as handle:
        while used < size:
            line = handle.readline(size - used)
            assert line
            h.update(line)
            head, sep, buf = line.partition(b': ')
            off = used + len(head) + 2
            used += len(line)
            if not sep:
                continue
            if (head.endswith(b'/item3/exponent_branch_inventory') or head.endswith(b'/item3/object')) and head.startswith(b'reading/'):
                tag = head.decode().removeprefix('reading/')
                name, part = tag.rsplit('/item3/', 1)
                name = re.sub(r'\{([^}]+)\}', lambda m: m[1].replace(', ', '+'), name)
                if part == 'object':
                    assert name not in operands
                    operands[name] = f'{engine}_reading:{off}'
                else:
                    assert name not in branches
                    branches[name] = buf.decode().strip()
            elif head in {b'family/velocity_changes', b'family/centre_velocity_changes', b'reading/dimension_argument_audit', b'family/dimension_signatures_including_changes', b'field_argument_dimension_signatures'}:
                context[head.decode()] = (f'{engine}_reading:{off}', buf.decode().strip())
    assert h.hexdigest() == digest and used == size, (engine, h.hexdigest())
    print('native', engine, used, h.hexdigest(), 'inventories', len(branches), 'emitted objects', len(operands))
    return branches, operands, context


def solution_space(chosen):
    """Rational Gaussian elimination for an affine plane intersection."""
    mat = [list(row[:4]) + [-row[-1]] for row in chosen]
    pivots = []
    for col in range(4):
        level = len(pivots)
        pivot = next((i for i in range(level, len(mat)) if mat[i][col]), None)
        if pivot is None:
            continue
        mat[level], mat[pivot] = mat[pivot], mat[level]
        mat[level] = [v / mat[level][col] for v in mat[level]]
        for i, row in enumerate(mat):
            if i != level and row[col]:
                mat[i] = [a - row[col]*b for a, b in zip(row, mat[level])]
        pivots.append(col)
    if any(row[-1] for row in mat[len(pivots):]):
        return None
    return mat[:len(pivots)], pivots


def signature(rows, point):
    lifted = (*point, Q(1))
    return tuple(i for i, row in enumerate(rows) if not sum(a*b for a,b in zip(row, lifted)))


@lru_cache(None)
def arrangement(rows):
    # Every realizable zero set has a basis of at most four independent normals.
    # For each basis construct a rational point in the implied flat, avoiding
    # every proper cut. Each excluded plane restricts to a polynomial degree <=4
    # on the free-coordinate moment curve, hence at most 4 bad integer parameters.
    flats = {}
    for n in range(5):
        for chosen in itertools.combinations(rows, n):
            if len({row[:4] for row in chosen}) != n:
                continue
            solved = solution_space(chosen)
            if solved is None:
                continue
            echelon, pivots = solved
            if len(pivots) != n:
                continue
            free = [j for j in range(4) if j not in pivots]
            def witness(t):
                point = [Q(0)]*4
                for k, j in enumerate(free, start=1):
                    point[j] = Q(t)**k
                for row, j in zip(echelon, pivots):
                    point[j] = row[-1] - sum(row[i]*point[i] for i in free)
                return point
            def identically_zero(row):
                # Substitute pivot coordinates and inspect the affine coefficients.
                return (row[-1] + sum(row[j]*eq[-1] for eq,j in zip(echelon,pivots)) == 0
                        and all(row[i] == sum(row[j]*eq[i] for eq,j in zip(echelon,pivots)) for i in free))
            implied = tuple(i for i,row in enumerate(rows) if identically_zero(row))
            if implied in flats:
                assert flats[implied] == n
                continue
            for t in range(1, 4*len(rows)+2):
                if signature(rows, witness(t)) == implied:
                    flats[implied] = n
                    break
            else:
                raise AssertionError(('missing witness', implied))
    return flats


def native_slots(py, wl, source):
    src = source.splitlines()
    time = re.search(r'timeDim=\{([\d,]+)\}', source)[1]
    assert 'centres = Table[xx[k, j][t], {k, 2}, {j, 3}];' in src
    line, = [i for i,s in enumerate(src,1) if 'register[centreVelocityChange[j]' in s]
    decl = src[line-1]
    m = re.fullmatch(r'Do\[register\[centreVelocityChange\[j\],lenDim-timeDim,\{(\w+)\}\],\{j,(\d+)\}\];', decl)
    assert m and m[1]=='timeDim'
    timevec = list(map(int,time.split(',')))
    # The native source prints the argument table before registration, and only
    # the result-dimensional table thereafter. Recover the type from the declaration.
    assert 'centreVelocityChange' not in wl['field_argument_dimension_signatures'][1]
    oldrules = wl['family/velocity_changes']
    mapping = re.findall(r'Derivative\[1\]\[xx\[1, (\d)\]\]\[t\] -> centreVelocityChange\[(\d)\]\[t\]',oldrules[1])
    assert len(mapping) == int(m[2])
    native_audit = py['reading/dimension_argument_audit']
    audited = ast.parse(native_audit[1], mode='eval').body
    refs = {}
    for key, value in zip(audited.keys, audited.values):
        triple = ast.literal_eval(key)
        if not re.fullmatch(r'X_1_\d_velocity_change', str(triple[0])):
            continue
        j = triple[0].split('_')[2]
        assert triple[1:] == (0,'t')
        assert (j,j) in mapping and f'X_1_{j}_velocity(t): {triple[0]}(t)' in py['family/centre_velocity_changes'][1]
        assert f'"centreVelocityChange[{j}]" -> {{0, 1, -1}}' in wl['family/dimension_signatures_including_changes'][1]
        dims = []
        for dim in value.elts:
            assert dim.func.id=='frozenset'
            dims.append([list(row) for row in ast.literal_eval(dim.args[0])])
        assert dims == [[timevec],[timevec]]
        refs[triple[0]] = dict(argument=list(triple), expected=dims[0], observed=dims[1], py_record=native_audit[0],
                wl_velocity_rules=oldrules[0], wl_head=f'centreVelocityChange[{j}]',wl_slot=1, wl_expected=timevec,
                status='matches_source_declared_slot', wl_source=dict(line=line,declaration=decl,sha256=hashlib.sha256(source.encode()).hexdigest(),
                    path=str(P / 'mathematica/S15_separation_departure.wl'),revision='cc62dc9f'))
    assert len(refs) == 3
    return refs


def main():
    py, pyids, pycontext = read_stream('py')
    wl, wlids, wlcontext = read_stream('wl')
    assert py.keys()==wl.keys()==pyids.keys()==wlids.keys() and len(py)==336
    wl_source = (PACKET/P/'mathematica/S15_separation_departure.wl').read_text()
    slots = native_slots(pycontext, wlcontext, wl_source)
    products = []
    phases, aliases, forms, shapes = Counter(), Counter(), Counter(), Counter()
    with tarfile.open(ARCHIVE) as ar:
        cfg = member(ar,'configuration.json')
        archive_names = ar.getnames()
        for i,name in enumerate(cfg['required_objects']):
            pyplanes, wlplanes = inventory('py',py[name]), inventory('wl',wl[name])
            assert pyplanes==wlplanes, (name, pyplanes, wlplanes)
            obj = member(ar,f'objects/{i:03}.json.gz')
            assert obj['object']==name and obj['index']==i
            row, = [r for r in obj['rows'] if r['field']=='item3/object']
            operands = row['left']+row['right']
            assert operands==[pyids[name],wlids[name]], name
            chain = set()
            while 'same_comparison_as' in row:
                ix = row['same_comparison_as']
                assert ix not in chain
                chain.add(ix)
                row = obj['rows'][ix]
            aliases[(bool(pyplanes),bool(chain))]+=1
            if not pyplanes:
                continue
            rows = pyplanes
            sample_rows = {'base': row.get('samples',[])}
            for section in ('refinement','log_branches'):
                path = f'{section}/{i:03}.json.gz'
                if path not in archive_names:
                    continue
                extra=member(ar,path)
                assert extra['object']==name and extra['build_id']==obj['build_id']
                if section=='refinement':
                    field, = [f for f in extra['fields'] if f['field']=='item3/object']
                    assert field['left']+field['right']==operands
                    sample_rows['additional'] = field['additional_samples']
                    sample_rows['precision'] = field['rechecks']
                else:
                    sample_rows['log'] = [f for f in extra['rechecks'] if f['field']=='item3/object']
            points = set()
            recovered, needs_recheck = set(), set()
            for phase, samples in sample_rows.items():
                phases[phase] += len(samples)
                for s in samples:
                    key = (tuple(s['case']),tuple(s['exponents']))
                    forms[s['case'][0]]+=1
                    if s['status'] in ('unresolved','separated'):
                        needs_recheck.add(key)
                        continue
                    assert s['status'] in ('within_estimate','within_precision_estimate'), s
                    if phase in ('log','precision'):
                        recovered.add(key)
                    points.add(tuple(map(Q,s['exponents'])))
            assert needs_recheck == recovered, (name,len(needs_recheck),len(recovered))
            sampled = {signature(rows,p) for p in points}
            attainable = arrangement(rows)
            assert sampled <= attainable.keys()
            missing = sorted(attainable.keys()-sampled)
            shapes[(len(rows),len(sampled),len(missing))]+=1
            products.append(dict(object=name,index=i,operand_ids=operands,planes=[list(map(str,r)) for r in rows],
                    sampled=sorted(sampled),untested=missing,ranks=dict(Counter(attainable[s] for s in missing))))
    print('aliases (has planes, resolved via archive alias):',dict(sorted(aliases.items())))
    print('scoped sample rows by archive phase',dict(sorted(phases.items())),'by form',dict(sorted(forms.items())))
    print('(planes, sampled selections, untested selections) -> objects:',dict(sorted(shapes.items())))
    print('source-slot references:',[(k,s['wl_head'],s['wl_slot'],s['wl_expected']) for k,s in sorted(slots.items())])
    result={'objects':products,'slots':slots}
    (HERE/'reference.json').write_text(json.dumps(result,separators=(',',':')))
    print('native arrangement totals',len(products),sum(bool(r['untested']) for r in products),sum(len(r['sampled']) for r in products),sum(len(r['untested']) for r in products))
    print('independent native reference written',hashlib.sha256((HERE/'reference.json').read_bytes()).hexdigest())

if __name__=='__main__':
    main()
