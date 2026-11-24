// Continuous intersection of a moving affine triangle with a fixed closed
// unit box. The polar fraction y in [0,1] sweeps exactly one native prism.
// A duplicated triangle vertex gives the swept line segment at the weld.
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {
using Real=long double;
using Vec=std::array<Real,3>;
using Triangle=std::array<Vec,3>;
using Poly=std::array<Real,4>;
using Roots=std::vector<Real>;
constexpr Vec constant{0,0,1};

Vec minus(const Vec &a,const Vec &b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
Real det(const Vec &a,const Vec &b,const Vec &c) {
    return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]);
}
Real value(const Poly &p,Real y) {return ((p[3]*y+p[2])*y+p[1])*y+p[0];}
void interior(Roots &out,Real v) {if(v>0 && v<1)out.push_back(v);}
void events(const Poly &p,Roots &out) {
    // Include stationary points even for a double/tangent contact. Between
    // consecutive extrema the sign of a nonconstant polynomial changes at
    // most once. Constant zero polynomials introduce no boundaries.
    if(!p[1]&&!p[2]&&!p[3])return;
    if(!p[2]&&!p[3]){interior(out,-p[0]/p[1]);return;}
    Roots knots{0,1};
    if(p[3]==0){interior(knots,-p[1]/(2*p[2]));}
    else {
        Real a=3*p[3],b=2*p[2],d=b*b-4*a*p[1];
        if(d>=0){
            Real q=-.5L*(b+std::copysign(std::sqrt(d),b));
            if(q==0)interior(knots,-b/(2*a));
            else{interior(knots,q/a);interior(knots,p[1]/q);}
        }
    }
    std::sort(knots.begin(),knots.end());
    for(size_t i=1;i<knots.size();i++){
        Real lo=knots[i-1],hi=knots[i],vl=value(p,lo),vh=value(p,hi);
        interior(out,lo);
        if(!vl||!vh||((vl<0)==(vh<0)))continue;
        for(int iteration=0;iteration<80;iteration++){
            Real mid=(lo+hi)/2;
            if(mid==lo||mid==hi)break;
            if((value(p,mid)<0)==(vl<0))lo=mid;else hi=mid;
        }
        interior(out,(lo+hi)/2);
    }
}
struct Line {Vec l,d;};
using Lines=std::array<Line,10>;
Lines constraints(const Triangle &bottom,const Triangle &top,Real pad) {
    // h(a,b)=A*a+B*b+C >=0 in the 2D material barycentric chart. The
    // three barycentric sides close this polygon independently of degeneracy
    // of the mapped spatial triangle. The last row yields pair determinants.
    Lines lines{};
    lines[0].l={1,0,0}; lines[1].l={0,1,0}; lines[2].l={-1,-1,1}; lines[9].l=constant;
    for(int c=0;c<3;c++)for(int side=0;side<2;side++){
        auto abc=[&](const Triangle &v) {
            Real s=side?-1:1;
            return Vec{s*(v[1][c]-v[0][c]),s*(v[2][c]-v[0][c]),s*v[0][c]+side+pad};
        };
        auto l=abc(bottom);
        lines[3+2*c+side]={l,minus(abc(top),l)};
    }
    return lines;
}
bool feasible(const Lines &lines,Real y) {
    using Point=std::array<Real,2>;
    std::array<Point,16> poly{},tmp{};
    poly[0]={0,0};poly[1]={1,0};poly[2]={0,1};int n=3;
    for(int i=3;i<9;i++){
        Vec h;
        for(int c=0;c<3;c++)h[c]=lines[i].l[c]+y*lines[i].d[c];
        int m=0;
        for(int j=0;j<n;j++){
            const auto &a=poly[j],&b=poly[(j+1)%n];
            Real va=h[0]*a[0]+h[1]*a[1]+h[2],vb=h[0]*b[0]+h[1]*b[1]+h[2];
            if(va>=0)tmp[m++]=a;
            if((va<0 && vb>0)||(va>0 && vb<0)){
                Real s=va/(va-vb);
                tmp[m++]={a[0]+s*(b[0]-a[0]),a[1]+s*(b[1]-a[1])};
            }
        }
        if(!m)return false;
        std::swap(tmp,poly);n=m;
    }
    return true;
}
bool overlaps(const Lines &lines) {
    if(feasible(lines,0)||feasible(lines,.5)||feasible(lines,1))return true;
    Roots cuts{0,1};
    // For a pair of nonparallel lines, feasibility of their intersection is
    // determined by the signs of the pair determinant and all three-line
    // determinants. A nonempty clipped reference triangle has such a vertex.
    for(int i=0;i<8;i++)for(int j=i+1;j<9;j++)for(int k=j+1;k<10;k++){
        auto a=lines[i],b=lines[j],c=lines[k];
        events({det(a.l,b.l,c.l),det(a.d,b.l,c.l)+det(a.l,b.d,c.l)+det(a.l,b.l,c.d),
                det(a.d,b.d,c.l)+det(a.d,b.l,c.d)+det(a.l,b.d,c.d),det(a.d,b.d,c.d)},cuts);
    }
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    for(size_t i=1;i<cuts.size();i++)
        if(feasible(lines,cuts[i])||feasible(lines,(cuts[i-1]+cuts[i])/2))return true;
    return false;
}
}
extern "C" int intersects(const double *corners,const double *lo,const double *hi,double pad) {
    Triangle bot{},top{};
    for(int p=0;p<3;p++)for(int k=0;k<3;k++){
        Real width=static_cast<Real>(hi[k])-lo[k];
        if(!(width>0))return -1;
        bot[p][k]=(static_cast<Real>(corners[p*3+k])-lo[k])/width;
        top[p][k]=(static_cast<Real>(corners[9+p*3+k])-lo[k])/width;
    }
    return overlaps(constraints(bot,top,pad));
}
