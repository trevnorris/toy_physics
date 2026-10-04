"""Exact affine geometry only. Actual scientific plans run under the shared guard.

Tests use an unrelated synthetic quadratic field/domain. No quadrature nodes,
source coefficients, kernels or transform values are constructed here.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
from math import isqrt


def require(value, message):
    if value is not True:
        raise ValueError(message)


def rational(value):
    require(type(value) in (int, str, F), 'exact rational input only')
    return F(value)


@total_ordering
@dataclass(frozen=True, eq=False)
class Quad:
    a: F
    b: F
    d: F

    def __post_init__(self):
        for name in ('a', 'b', 'd'):
            object.__setattr__(self, name, rational(getattr(self, name)))
        require(self.d > 0, 'positive basis square')
        square = isqrt(self.d.numerator)**2 == self.d.numerator and isqrt(self.d.denominator)**2 == self.d.denominator
        require(not square, 'genuine quadratic field')

    def other(self, value):
        if isinstance(value, Quad):
            require(value.d == self.d, 'same exact field')
            return value
        return Quad(rational(value), 0, self.d)

    def __add__(self, value):
        z = self.other(value)
        return Quad(self.a+z.a, self.b+z.b, self.d)
    __radd__ = __add__

    def __neg__(self): return Quad(-self.a, -self.b, self.d)
    def __sub__(self, value): return self + (-self.other(value))
    def __rsub__(self, value): return self.other(value) - self

    def __mul__(self, value):
        z = self.other(value)
        return Quad(self.a*z.a+self.b*z.b*self.d, self.a*z.b+self.b*z.a, self.d)
    __rmul__ = __mul__

    def __truediv__(self, value):
        z = self.other(value)
        denominator = z.a*z.a-z.b*z.b*self.d
        require(denominator != 0, 'nonzero exact divisor')
        return self * Quad(z.a/denominator, -z.b/denominator, self.d)

    def sign(self):
        if not self.b: return (self.a > 0)-(self.a < 0)
        if not self.a: return (self.b > 0)-(self.b < 0)
        if self.a > 0 and self.b > 0: return 1
        if self.a < 0 and self.b < 0: return -1
        comparison = self.a*self.a-self.b*self.b*self.d
        require(comparison != 0, 'irrational basis cannot cancel nonzero rationals')
        return ((comparison > 0)-(comparison < 0)) * (1 if self.a > 0 else -1)

    def __eq__(self, value):
        if not isinstance(value, (Quad, int, F)): return False
        z = self.other(value)
        return self.a == z.a and self.b == z.b

    def __lt__(self, value): return (self-self.other(value)).sign() < 0
    def __hash__(self): return hash(self.a) if not self.b else hash((self.a, self.b, self.d))
    def packed(self): return {'a':str(self.a), 'b':str(self.b), 'basisSquare':str(self.d)}


@dataclass(frozen=True)
class Line:
    slope: F
    intercept: Quad
    labels: tuple

    def __post_init__(self):
        object.__setattr__(self, 'slope', rational(self.slope))
        require(isinstance(self.intercept, Quad) and isinstance(self.labels, tuple)
                and all(isinstance(v, str) for v in self.labels), 'exact labelled line types')

    def at(self, x): return self.slope*x+self.intercept
    def key(self): return (self.slope, self.intercept)
    def packed(self):
        return {'slope':str(self.slope), 'intercept':self.intercept.packed(), 'labels':list(self.labels)}


def canonical_lines(lines):
    grouped = {}
    for line in lines:
        require(isinstance(line, Line) and line.labels and len(set(line.labels)) == len(line.labels), 'labelled exact affine line')
        grouped.setdefault(line.key(), []).extend(line.labels)
    labels = [label for values in grouped.values() for label in values]
    require(len(labels) == len(set(labels)), 'unique original line label')
    return [Line(slope, intercept, tuple(sorted(labels))) for (slope, intercept), labels in sorted(grouped.items())]


def arrangement(lines, xlo, xhi, ylo, yhi, mandatory_cuts):
    """Complete vertical decomposition; exact cuts and coalescing labels retained."""
    require(xlo < xhi and ylo < yhi, 'positive rectangle')
    lines = canonical_lines(lines)
    require(any(l.slope == 0 and l.intercept == ylo for l in lines), 'lower box boundary supplied')
    require(any(l.slope == 0 and l.intercept == yhi for l in lines), 'upper box boundary supplied')
    cuts = {xlo:['box:xlo'], xhi:['box:xhi']}
    for x, label in mandatory_cuts:
        if xlo <= x <= xhi: cuts.setdefault(x, []).append(label)
    intersections = []
    for i, left in enumerate(lines):
        for j in range(i+1, len(lines)):
            right = lines[j]
            if left.slope == right.slope: continue
            x = (right.intercept-left.intercept)/(left.slope-right.slope)
            y = left.at(x)
            inside = xlo <= x <= xhi and ylo <= y <= yhi
            intersections.append({'lines':[i,j], 'x':x, 'y':y, 'inside':inside})
            if inside: cuts.setdefault(x, []).append('intersection:'+str(i)+':'+str(j))
    ordered = sorted(cuts)
    slabs = []
    for slab_id, (left, right) in enumerate(zip(ordered, ordered[1:])):
        middle = (left+right)/2
        ids = [i for i,l in enumerate(lines) if ylo <= l.at(middle) <= yhi]
        ids.sort(key=lambda i:lines[i].at(middle))
        cells = []
        for low, high in zip(ids, ids[1:]):
            wleft = lines[high].at(left)-lines[low].at(left)
            wright = lines[high].at(right)-lines[low].at(right)
            cells.append({'lowerLine':low, 'upperLine':high, 'leftWidth':wleft,
                'rightWidth':wright, 'area':(right-left)*(wleft+wright)/2})
        slabs.append({'id':slab_id, 'left':left, 'right':right, 'middle':middle,
            'boundaryLines':ids, 'cells':cells})
    return {'box':[xlo,xhi,ylo,yhi], 'lines':lines, 'cuts':ordered,
        'cutLabels':[{'x':x,'labels':sorted(cuts[x])} for x in ordered],
        'intersections':intersections, 'slabs':slabs}


def audit(plan, required_lines, required_cuts):
    """Independent incidence/coverage checks against the original line specification."""
    xlo,xhi,ylo,yhi = plan['box'];lines = plan['lines'];cuts = plan['cuts']
    require(cuts[0] == xlo and cuts[-1] == xhi and all(a < b for a,b in zip(cuts,cuts[1:])), 'ordered exact x coverage')
    require(len(plan['slabs']) == len(cuts)-1, 'complete slab census')
    require(len(plan['cutLabels']) == len(cuts) and [v['x'] for v in plan['cutLabels']] == cuts, 'complete cut provenance')
    actual_labels = {label:l.key() for l in lines for label in l.labels}
    require(len(actual_labels) == sum(len(l.labels) for l in lines), 'unique actual labels')
    required = canonical_lines(required_lines)
    require(actual_labels == {label:l.key() for l in required for label in l.labels}, 'required collision/resolution incidence')
    for x,label in required_cuts:
        if xlo <= x <= xhi:
            require(x in cuts and label in plan['cutLabels'][cuts.index(x)]['labels'], 'required vertical cut incidence')
    # Recheck ALL pairwise crossings independently, including boundary crossings.
    crossings=[]
    for i,l in enumerate(lines):
        for j in range(i+1,len(lines)):
            r=lines[j]
            if l.slope == r.slope: continue
            x = (r.intercept-l.intercept)/(l.slope-r.slope);y=l.at(x)
            inside=xlo <= x <= xhi and ylo <= y <= yhi
            crossings.append({'lines':[i,j],'x':x,'y':y,'inside':inside})
            if inside:
                require(x in cuts, 'all in-box line crossings split')
                require('intersection:'+str(i)+':'+str(j) in plan['cutLabels'][cuts.index(x)]['labels'], 'crossing label retained')
    require(plan['intersections']==crossings,'complete exact crossing provenance')
    area = xlo*0;cell_count=0
    for index,slab in enumerate(plan['slabs']):
        left,right=cuts[index:index+2];mid=(left+right)/2
        require(slab['id']==index and slab['left']==left and slab['right']==right and slab['middle']==mid, 'slab argument join')
        required_ids=sorted((i for i,l in enumerate(lines) if ylo <= l.at(mid) <= yhi),key=lambda i:lines[i].at(mid))
        ids=slab['boundaryLines'];require(ids==required_ids and len(ids)>=2, 'complete clipped line order')
        require(lines[ids[0]].at(mid)==ylo and lines[ids[-1]].at(mid)==yhi, 'complete y coverage')
        require(len(slab['cells'])==len(ids)-1,'cell census')
        for cell,low,high in zip(slab['cells'],ids,ids[1:]):
            require(cell['lowerLine']==low and cell['upperLine']==high,'oriented adjacent cell')
            a=lines[high].at(left)-lines[low].at(left);b=lines[high].at(right)-lines[low].at(right)
            require(a>=0 and b>=0 and lines[high].at(mid)>lines[low].at(mid), 'exact endpoint and interior orientation')
            expected=(right-left)*(a+b)/2
            require(cell['leftWidth']==a and cell['rightWidth']==b and cell['area']==expected and expected>0,'exact positive cell area')
            area+=expected;cell_count+=1
    expected=(xhi-xlo)*(yhi-ylo)
    require(area==expected,'exact union area equals full rectangle')
    return {'slabs':len(plan['slabs']),'cells':cell_count,'area':area,'expectedArea':expected,
        'coverage':True,'disjointInteriors':True,'allRequiredLines':True,
        'squareSubstitutionJacobians':'2*(half-width)*z > 0 for 0<z<1; both nested half-widths positive in each open cell',
        'numericalNodesEvaluated':0}


def specifications(kappa, K, U, carrier, width, length):
    """Caller supplies actual joined inputs; no hidden physical constants."""
    zero=kappa*0
    packet_offsets=[F(0)]+[sign*F(n)/width for n in [1,2,4,8] for sign in [-1,1]]
    profile_offsets=[F(0)]+[sign*F(n)/length for n in [1,2,4,8] for sign in [-1,1]]
    carrier_cuts=[(carrier+d,'carrier:'+str(d)) for d in packet_offsets]
    kcuts=[(-kappa,'external:k-'),(kappa,'external:k+')]+carrier_cuts
    square=[Line(F(0),zero-K,('box:l-',)),Line(F(0),zero+K,('box:l+',)),
        Line(F(0),-kappa,('external:l-',)),Line(F(0),kappa,('external:l+',))]
    square += [Line(F(-1),n*kappa,('collision:sum:'+str(n),)) for n in [-2,0,2]]
    square += [Line(F(1),zero+d,('profile:transfer:'+str(d),)) for d in profile_offsets]
    square += [Line(F(0),v,('resolution:l:'+label,)) for v,label in carrier_cuts]
    height=[Line(F(0),zero,('box:Q0',)),Line(F(0),zero+U,('box:QU',))]
    height += [Line(F(0),zero+d,('profile:Q:'+str(d),)) for d in profile_offsets if d>0]
    # Both Y(k+Q) and Y(k-Q) branches: retain affine labels before clipping Q>=0.
    y_targets=[(-kappa,'branch:l-'),(kappa,'branch:l+')]+[(v,'resolution:'+label) for v,label in carrier_cuts]
    height += [Line(F(-sign),sign*v,('height:'+str(sign)+':'+label,)) for v,label in y_targets for sign in [-1,1]]
    return {'square':(square,kcuts,[-K+zero,K+zero,-K+zero,K+zero]),
            'height':(height,kcuts,[-K+zero,K+zero,zero,U+zero])}


def static_counts(summary):
    # Occurrences, not certified distinct irrational numerical coordinates.
    cells=summary['cells'];slabs=summary['slabs']
    return {'A24':{'outerNodeOccurrences':cells*48**2,'kNodeDescriptors':slabs*48},
            'A48':{'outerNodeOccurrences':cells*96**2,'kNodeDescriptors':slabs*96},
            'B_initial':{'outerNodeOccurrencesLowerBound':cells*15**2,'kNodeOccurrencesLowerBound':slabs*15,
                'adaptiveRefinements':'UNKNOWN','sameGeometricCells':True},
            'distinctYArguments':'NOT_DEDUCED_WITHOUT_ACTUAL_NODE_DESCRIPTORS',
            'oldBankRequestsReused':0,'oldBankReuseClaim':'No actual outer numerical arguments evaluated or matched'}
