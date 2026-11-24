// Five fixed-knot shape derivatives on the exact and trilinear inverse of one smooth map.
// The signed probe calls the same spatial() energy used by the production Eulerian evaluator.
#include "mouth_cp10_spatial_crystal.cpp"
#define SMOOTH_HAS_SOLID
#include "mouth_r4_smooth_points.hpp"

namespace {
constexpr int modes=5, probes=20, columns=2+modes+probes;
struct Direction {double R; Vec grad;};
std::array<Direction,modes> fixed_knot_motions(const Data &d,const Vec &pos){
    const double x=pos[0],t=pos[1],a=pos[2],u=t*(1-t),b=d.shape[0]+1.5;
    std::array<Direction,modes> v{};
    const double db[]={b*t,4*u,4*u*(2*t-1),16*u*u};
    const double dbt[]={b,4*(1-2*t),4*(-6*t*t+6*t-1),32*u*(1-2*t)};
    for(int i=0;i<4;i++)v[i]={x*db[i],{db[i],x*dbt[i],0}};
    const double y=(t-.06)/.56;
    if(y>0 && y<1){
        const double c=std::cos(a),p=(3*c*c-1)/2,h=16*y*y*(1-y)*(1-y),ht=32*y*(1-y)*(1-2*y)/.56;
        v[4]={x*x*x*h*p,{3*x*x*h*p,x*x*x*ht*p,-3*x*x*x*h*std::sin(a)*c}};
    }
    return v;
}
}

extern "C" int clearance_five(const Data *ds,int level,const double *order,const unsigned char *mask,
                     Grid rule,const double *weights,long double *out,double *jac,long long *count){
    const Data &d=ds[0];
    const size_t nc=(d.field_x.n-1)*(d.field_t.n-1)*(d.field_a.n-1);
    std::array<Geometry,probes+1> geom;
    for(int i=0;i<=probes;i++)geom[i]={ds[i],{0.,0.,0.},{1.,1.,1.}};
    jac[0]=jac[1]=std::numeric_limits<double>::infinity();
    const double cut=d.x.v[(d.x.n-1)*2/3];
    smooth::walk(d,level,order,mask,rule,weights,[&](const smooth::Point &p){
        const auto motion=fixed_knot_motions(d,p.pos);
        const double sn=std::sin(p.pos[2]);
        for(int bank=0;bank<2;bank++){
            const double xi=p.xi[bank];
            const auto ref=reference(d,xi,p.pos[1],profile(d,p.pos[1],true),xi<=cut);
            const Mat right=product(smooth::forward(p,bank),ref.inv);
            const double detright=1./(p.dx[bank]*ref.det2*ref.S);
            const auto base=spatial(geom[0],p.pos,right,detright,ref.S,sn,p.chi);
            if(!(base.J>0))throw std::runtime_error("Nonpositive active crystal orientation");
            jac[bank]=std::min(jac[bank],base.J);
            const int idx=(p.cell[0]*(d.field_t.n-1)+p.cell[1])*(d.field_a.n-1)+p.cell[2];
            if(idx<0 || size_t(idx)>=nc)throw std::runtime_error("Invalid native cell index");
            auto *row=out+(idx*2+bank)*columns;
            const double refvol=p.weight*base.volume/base.J, weight=refvol*p.chi;
            row[0]+=weight*base.W;
            row[1]+=refvol*p.minus_dchi*base.W;
            const double volstress=-1.+.5*base.J*(base.J-1.);
            const Mat F=product(base.C,right);
            Vec conjugate{};
            for(int k=0;k<3;k++)conjugate[k]=dot(F[0],right[k]);
            conjugate[0]+=volstress/base.wall.grad[0];
            const double hoop=(dot(F[2],F[2])+base.hoop*base.hoop+2*volstress)/base.R;
            for(int m=0;m<modes;m++)row[2+m]-=weight*(dot(conjugate,motion[m].grad)+hoop*motion[m].R);
            for(int i=0;i<probes;i++){
                const auto v=spatial(geom[1+i],p.pos,right,detright,ref.S,sn,p.chi);
                if(!(v.J>0))throw std::runtime_error("Invalid signed variation");
                row[2+modes+i]+=p.weight*p.chi*v.volume/v.J*v.W;
            }
        }
        ++*count;
    });
    return 0;
}
