// Common physical quadrature for the analytic and sampled two-metric inverse charts.
// Each radial integration interval is cut at both charts' welds and every active
// triangular registry boundary. Moving signed knots and cuts carry their own tangent.
#pragma once
#ifndef SMOOTH_HAS_SOLID
#include "mouth_free_shape_solid.cpp"
#endif
#include <stdexcept>

namespace smooth {
constexpr double amp=.12;
struct Jet {
    double v,d;
    Jet(double a=0.,double b=0.):v(a),d(b){}
};
inline Jet operator+(Jet a,Jet b){return {a.v+b.v,a.d+b.d};}
inline Jet operator-(Jet a,Jet b){return {a.v-b.v,a.d-b.d};}
inline Jet operator-(Jet a){return {-a.v,-a.d};}
inline Jet operator*(Jet a,Jet b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
inline Jet operator/(Jet a,Jet b){double r=a.v/b.v;return {r,(a.d-r*b.d)/b.v};}
inline bool operator<(Jet a,Jet b){return a.v<b.v || (a.v==b.v&&a.d<b.d);}
inline bool operator==(Jet a,Jet b){return a.v==b.v && a.d==b.d;}
inline double val(double x){return x;} inline double val(Jet x){return x.v;}
inline double sin_(double x){return std::sin(x);} inline Jet sin_(Jet x){return {std::sin(x.v),std::cos(x.v)*x.d};}
inline double cos_(double x){return std::cos(x);} inline Jet cos_(Jet x){return {std::cos(x.v),-std::sin(x.v)*x.d};}
inline double exp_(double x){return std::exp(x);} inline Jet exp_(Jet x){double e=std::exp(x.v);return {e,e*x.d};}
inline double log_(double x){return std::log(x);} inline Jet log_(Jet x){return {std::log(x.v),x.d/x.v};}
inline double abs_(double x){return std::abs(x);} inline Jet abs_(Jet x){return x.v<0?-x:x;}
template<typename S> S amplitude(S x){return amp*x*(1.-x*x);}
template<typename S> S gradient(S x){return amp*(1.-3.*x*x);}
template<typename S> S mix(S a,S b,S v){return a+(b-a)*v;}
template<typename S> void knots(std::vector<S>& a){std::sort(a.begin(),a.end());a.erase(std::unique(a.begin(),a.end()),a.end());}
template<typename S> int bin(const std::vector<S>& axis,S t){return std::clamp(int(std::upper_bound(axis.begin(),axis.end(),t)-axis.begin())-1,0,int(axis.size())-2);}
inline std::vector<double> axis(Grid grid){return {grid.v,grid.v+grid.n};}
template<typename S> std::vector<S> refine(const std::vector<S>& native,int n){
    std::vector<S> r;for(size_t i=0;i+1<native.size();i++)for(int k=0;k<n;k++)r.push_back(native[i]+(native[i+1]-native[i])*(double(k)/n));
    r.push_back(native.back());return r;
}
inline double rate(double t,double L){return (t*L<L-.12?std::pow(t*L/(L-.12),25):1.)-t;}
// The normed order mesh is moved on the saved physical law. The reference cells
// and the localization bands stay at their exact checkpoint-10 chart coordinates.
template<typename S> struct Sample {
    std::array<S,3> pos{}; S chi{},minus_dchi{},weight{},xi[2]{},dx[2]{},dt[2]{};
    int cell[3]{},cell_id=0,bucket=0;
};
using Point=Sample<double>;
inline Mat forward(const Point &p,int b){return {Vec{1/p.dx[b],-p.dt[b]/p.dx[b],0},Vec{0,1,0},Vec{0,0,1}};}
// The monotone cubic has a unique root in [0,1]. Newton from the label radius
// is contractive here; differentiate its implicit equation for the root tangent.
inline double cubic_root(double y,double sine){
    double r=y;for(int i=0;i<9;i++){
        double dx=(r+amplitude(r)*sine-y)/(1+gradient(r)*sine);r-=dx;
        if(std::abs(dx)<1e-16)break;
    }return r;
}
inline Jet cubic_root(Jet y,Jet sine){double r=cubic_root(y.v,sine.v);return {r,(y.d-amplitude(r)*sine.d)/(1+gradient(r)*sine.v)};}

template<typename S,typename Fn>
void traverse(const Data& d,int level,const std::vector<S>& native_t,Grid fixed_t,
        const double* order,const unsigned char* bucket,Grid qr,const double* w,Fn act){
    const auto x=refine(axis(d.field_x),level);
    const auto t=refine(native_t,level);
    const auto native_x=axis(d.field_x);
    const int nx=d.field_x.n,nt=d.field_t.n,na=d.field_a.n;
    auto node=[&](int i,int j,int k){return order[(i*nt+j)*na+k];};
    std::vector<bool> ordered(nt-1,false);
    for(int i=0;i<nx;i++)for(int j=0;j<nt;j++)for(int k=0;k<na;k++)if(node(i,j,k)>0){if(j)ordered[j-1]=true;if(j+1<nt)ordered[j]=true;}
    std::vector<S> tc=t;
    for(auto v:axis(fixed_t))tc.emplace_back(v);
    for(auto v:axis(d.reg_t))tc.emplace_back(v);
    for(int k=0;k<=d.spline_n;k++)tc.emplace_back(double(k)/d.spline_n);
    tc.emplace_back(.06);tc.emplace_back(.62);knots(tc);
    std::vector<S> fixed_axis;for(auto v:axis(fixed_t))fixed_axis.emplace_back(v);
    std::vector<S> registry;for(auto v:axis(d.reg_t))registry.emplace_back(v);
    const double weld=d.x.v[(d.x.n-1)*2/3];
    const auto energies=axis(d.reg_e);
    for(size_t j=0;j+1<tc.size();j++){
        const S width=tc[j+1]-tc[j], tm=tc[j]+.5*width;
        const int J=bin(native_t,tm),jz=bin(t,tm),home_j=bin(fixed_axis,tm),reg_j=bin(registry,tm);
        if(!ordered[J])continue;
        const S bottom=t[jz],top=t[jz+1], height=top-bottom;
        const S s0=sin_(pi*bottom),s1=sin_(pi*top), slope=(s1-s0)/height;
        const S hn=native_t[J+1]-native_t[J];
        for(int Z=0;Z<qr.n;Z++){
            const S z=tc[j]+width*qr.v[Z], sine=sin_(pi*z), ds=pi*cos_(pi*z);
            const S chord=s0+(z-bottom)*slope, cz=(z-native_t[J])/hn;
            const S rz=(z-registry[reg_j])/(registry[reg_j+1]-registry[reg_j]);
            std::vector<S> yc{S(weld)};
            for(size_t e=0;e<energies.size();e++){
                const double by=std::exp(-d.pressure*energies[e]);
                if(by>weld&&by<1.)yc.emplace_back(by);
                if(e+1<energies.size()){
                    const S diag=exp_(-d.pressure*(energies[e]+(energies[e+1]-energies[e])*(1.-rz)));
                    if(S(weld)<diag&&diag<S(1.))yc.push_back(diag);
                }
            }
            const S knot_velocity=mix(S(rate(val(native_t[J]),d.length)),S(rate(val(native_t[J+1]),d.length)),cz);
            std::vector<S> radial;for(double u:x)radial.emplace_back(u);
            std::vector<S> samples;for(double u:x)samples.push_back(u+amplitude(u)*chord);
            for(S y:yc){
                radial.push_back(cubic_root(y,sine));
                const int i=bin(samples,y);
                radial.push_back(x[i]+(x[i+1]-x[i])*(y-samples[i])/(samples[i+1]-samples[i]));
            }
            knots(radial);
            for(size_t i=0;i+1<radial.size();i++){
                const S wx=radial[i+1]-radial[i],xm=radial[i]+.5*wx;
                // Radial cell boundaries are stationary for the depth motion. At a zero-width
                // ribbon entries are selected by the limiting sign using the moving cut.
                int I=std::clamp(int(std::upper_bound(x.begin(),x.end(),xm,
                    [](S a,double b){return a<S(b);})-x.begin())-1,0,int(x.size())-2);
                const int i0=I/level;
                const double g0=amplitude(x[I]),g1=amplitude(x[I+1]);
                const double gp=(g1-g0)/(x[I+1]-x[I]);
                for(int X=0;X<qr.n;X++){
                    const S r=radial[i]+wx*qr.v[X],g=g0+gp*(r-x[I]);
                    Sample<S> p; p.pos[0]=r;p.pos[1]=z;
                    p.xi[0]=r+amplitude(r)*sine; p.dx[0]=1.+gradient(r)*sine;p.dt[0]=amplitude(r)*ds;
                    p.xi[1]=r+g*chord;p.dx[1]=1.+gp*chord;p.dt[1]=g*slope;
                    p.cell[0]=i0;p.cell[1]=home_j;
                    const S cr=(r-native_x[i0])/(native_x[i0+1]-native_x[i0]);
                    const S wt=width*wx*w[Z]*w[X];
                    for(int k=0;k<na-1;k++){
                        double left[2][2],right[2][2],largest=0.;
                        for(int ii=0;ii<2;ii++)for(int zz=0;zz<2;zz++){
                            left[ii][zz]=node(i0+ii,J+zz,k);right[ii][zz]=node(i0+ii,J+zz,k+1);
                            largest=std::max(largest,std::max(left[ii][zz],right[ii][zz]));
                        }
                        if(largest==0.)continue;
                        auto radial_chi=[&](double c[2][2],int zz){return c[0][zz]+cr*(c[1][zz]-c[0][zz]);};
                        const S lb=radial_chi(left,0),lt=radial_chi(left,1),rb=radial_chi(right,0),rt=radial_chi(right,1);
                        const S lc=mix(lb,lt,cz),rc=mix(rb,rt,cz),ls=(lt-lb)/hn,rs=(rt-rb)/hn;
                        p.cell[2]=k;p.cell_id=(i0*(nt-1)+home_j)*(na-1)+k;p.bucket=bucket[p.cell_id];
                        double angle=d.field_a.v[k],da=d.field_a.v[k+1]-angle;
                        for(int A=0;A<qr.n;A++){
                            double f=qr.v[A];p.pos[2]=angle+da*f;p.weight=wt*(da*w[A]);
                            p.chi=lc+(rc-lc)*f;p.minus_dchi=(ls+(rs-ls)*f)*knot_velocity;act(p);
                        }
                    }
                }
            }
        }
    }
}
template<typename Fn> void walk(const Data &d,int level,const double*order,const unsigned char*bucket,
               Grid rule,const double*weights,Fn action){
    traverse(d,level,axis(d.field_t),d.field_t,order,bucket,rule,weights,action);
}
} // namespace smooth
