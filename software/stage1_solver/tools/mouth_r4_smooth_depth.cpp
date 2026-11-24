// Analytic log-depth response with moving quadrature and labels, checked against
// signed real energies from the production two-metric constitutive evaluator.
#include "mouth_r4_smooth_points.hpp"
namespace smooth {
struct RefTangent {Jet S,sx,st,zx,zt;};
RefTangent natural(const Data&d,Jet x,Jet t){
    const double cut=d.x.v[(d.x.n-1)*2/3];
    if(x<Jet(cut)){
        auto p=profile(d,t.v,true);
        Jet scale=d.core_ratio*Jet{p.b,p.bt*t.d};
        Jet slope=d.core_ratio*Jet{p.bt,p.btt*t.d};
        return {x*scale,scale,x*slope,0.,d.ref_length};
    }
    Jet e=-log_(x)/d.pressure;
    auto locate=[](Grid g,Jet v){return std::clamp(int(std::upper_bound(g.v,g.v+g.n,v,
                                  [](Jet a,double b){return a<Jet(b);})-g.v)-1,0,g.n-2);};
    const int i=locate(d.reg_e,e),j=locate(d.reg_t,t),ne=d.reg_e.n;
    const double de=d.reg_e.v[i+1]-d.reg_e.v[i],dt=d.reg_t.v[j+1]-d.reg_t.v[j];
    Jet u=(e-d.reg_e.v[i])/de,v=(t-d.reg_t.v[j])/dt;
    bool upper=Jet(1)<u+v;
    int base=j*ne+i;double sign=upper?-1.:1.;
    const double *p=d.labels+2*(base+(upper?ne+1:0)),*r=d.labels+2*(base+(upper?ne:1)),*z=d.labels+2*(base+(upper?1:ne));
    if(upper){u=1-u;v=1-v;}
    return {4.+p[0]+u*(r[0]-p[0])+v*(z[0]-p[0]),(r[0]-p[0])*sign/(-d.pressure*x*de),
            (z[0]-p[0])*sign/dt,(r[1]-p[1])*sign/(-d.pressure*x*de),(z[1]-p[1])*sign/dt};
}
// Use the production gradient and Hessian of the wall to differentiate its
// location and metric. Log length has the explicit variation dL=L dq.
struct MovingFrame {Jet radius,rx,rt,ra,length;double sine;};
MovingFrame frame(const Data&d,const Sample<Jet>& p,int direction){
    Vec v{},dotp{};for(int c=0;c<3;c++){v[c]=p.pos[c].v;dotp[c]=p.pos[c].d;}
    const Map r=radius(d,v);
    return {{r.R,dot(r.grad,dotp)},{r.grad[0],dot(r.hess[0],dotp)},{r.grad[1],dot(r.hess[1],dotp)},
            {r.grad[2],dot(r.hess[2],dotp)},{d.length,direction*d.length},std::sin(v[2])};
}
Jet unit_density(const Data&d,const Sample<Jet>&p,int bank,const MovingFrame&c){
    const auto ref=natural(d,p.xi[bank],p.pos[1]);
    const auto det=ref.sx*ref.zt-ref.zx*ref.st;
    const auto radial=c.rx/p.dx[bank],axial=c.rt-radial*p.dt[bank];
    const auto F00=(radial*ref.zt-axial*ref.zx)/det, F01=(-radial*ref.st+axial*ref.sx)/det;
    const auto F10=-c.length*ref.zx/det, F11=c.length*ref.sx/det;
    const auto p_arc=c.radius/ref.S, F02=c.ra/ref.S;
    const auto reference_volume=p.dx[bank]*abs_(det)*ref.S*ref.S;
    const auto J=c.rx*c.length*c.radius*c.radius/reference_volume;
    if(!(J.v>0 && std::isfinite(J.v+J.d)))throw std::runtime_error("nonpositive material orientation in smooth chart");
    const auto W=.5*(F00*F00+F01*F01+F10*F10+F11*F11+F02*F02+2*p_arc*p_arc-4.)-log_(J)+.25*(J-1.)*(J-1.);
    return 4*pi*c.sine*reference_volume*W;
}
// Separate real-valued energy path using the original 3x3 reference tensor and
// 4D angular constitutive law, independent of the compact tangent expressions.
double production(const Data&d,const Point&p,int b,const Map&r){
    const auto ref=reference(d,p.xi[b],p.pos[1],profile(d,p.pos[1],true),p.xi[b]<=d.x.v[(d.x.n-1)*2/3]);
    const auto f=forward(p,b);
    const Mat C{r.grad,Vec{0.,d.length,0.},Vec{0.,0.,r.R}};
    const Mat F=product(C,product(f,ref.inv));
    const double hoop=r.R/ref.S, J=determinant(C)*hoop/(p.dx[b]*ref.det2*ref.S);
    if(!(J>0. && std::isfinite(J)))throw std::runtime_error("invalid signed chart");
    const double W=.5*(norm2(F)+hoop*hoop-4)-std::log(J)+.25*(J-1.)*(J-1.);
    return 4*pi*std::sin(p.pos[2])*determinant(C)*r.R*W/J;
}
} // namespace smooth

extern "C" int smooth_depth_jet(const Data*d,int level,const double*order,const unsigned char*bucket,
           Grid rule,const double*weights,int direction,double*row,long long*count){
    using namespace smooth;
    long double sums[3][2][3]{}; // static energy, total depth force, explicit order work
    std::vector<Jet> t;for(double t0:axis(d->field_t))t.emplace_back(t0,direction*rate(t0,d->length));
    try{
        traverse(*d,level,t,d->field_t,order,bucket,rule,weights,[&](const Sample<Jet>&p){
            const auto c=frame(*d,p,direction);
            for(int b=0;b<2;b++){
                auto density=unit_density(*d,p,b,c),E=density*p.chi*p.weight;
                sums[p.bucket][b][0]+=E.v;
                sums[p.bucket][b][1]+=-direction*E.d;
                sums[p.bucket][b][2]+=density.v*p.minus_dchi.v*p.weight.v;
            }
            ++*count;
        });
    } catch(const std::exception &){return 1;}
    int k=0;for(auto &region:sums)for(auto &bank:region)for(auto value:bank)row[k++]=double(value);
    return 0;
}
extern "C" int smooth_depth_fd(const Data*d,Grid fixed_t,int level,const double*order,const unsigned char*bucket,
           Grid rule,const double*weights,double*row,long long*count){
    using namespace smooth;
    long double sums[3][2]{};
    try{
        traverse(*d,level,axis(d->field_t),fixed_t,order,bucket,rule,weights,[&](const Point&p){
            const auto map=radius(*d,p.pos);
            for(int b=0;b<2;b++)sums[p.bucket][b]+=(long double)p.weight*p.chi*production(*d,p,b,map);
            ++*count;
        });
    } catch(const std::exception &){return 1;}
    int k=0;for(auto &region:sums)for(auto value:region)row[k++]=double(value);
    return 0;
}
