// Forward derivative of the physical label gradient on the registered cubic
// collar. Geometry, order, quadrature rules and the physical wall are frozen.
#include "mouth_r4_eulerian_solid.cpp"
#define SMOOTH_HAS_SOLID
#include "mouth_r4_smooth_points.hpp"

namespace smooth {
inline double scalar_value(Jet x){return x.v;}
inline Jet scalar_log(Jet x){return log_(x);}
inline Jet scalar_sin(Jet x){return sin_(x);}
inline Jet scalar_cos(Jet x){return cos_(x);}
inline Jet scalar_abs(Jet x){return abs_(x);}
inline Jet scalar_sqrt(Jet x){double r=std::sqrt(x.v);return {r,x.d/(2*r)};}
inline Jet& operator+=(Jet& a,Jet b){return a=a+b;}
inline Jet& operator-=(Jet& a,Jet b){return a=a-b;}
inline Jet& operator*=(Jet& a,Jet b){return a=a*b;}
}

namespace {
using smooth::Jet;
struct WeldCut {Jet u; Vector3<Jet> normal;};

void split_weld(std::vector<WeldCut>& cuts,const Data&,const Vector3<Jet>& left,
                const Vector3<Jet>& right,double weld){
    cuts.clear();
    cuts.push_back({0.,{}});
    if(right[0].v!=left[0].v){
        const Jet u=(weld-left[0])/(right[0]-left[0]);
        if(0<u.v && u.v<1)cuts.push_back({u,{1.,0.,0.}});
    }
    cuts.push_back({1.,{}});
}
}

extern "C" int r4_physical_label_hvp(const Data* data,const double* labels,const double* order,
        const double* direction,Grid rule,const double* weights,double* product_out,int* problem){
    const Data& d=*data;
    if(!d.collar_coeff)return 1;
    const int nx=d.field_x.n,nt=d.field_t.n,na=d.field_a.n;
    const double weld=d.x.v[(d.x.n-1)*2/3];
    std::vector<WeldCut> cuts;
    std::vector<Jet> cut_work;
    cuts.reserve(3);
    cut_work.reserve(3);
    for(int i=0;i<nx-1;i++)for(int j=0;j<nt-1;j++)for(int k=0;k<na-1;k++){
        LabelCell<Jet> q{};
        Q1 base{};
        int ids[2][2][2];
        bool ordered=false,moving=false;
        for(int x=0;x<2;x++)for(int t=0;t<2;t++)for(int a=0;a<2;a++){
            const int n=ids[x][t][a]=((i+x)*nt+j+t)*na+k+a;
            q.order[x][t][a]=order[n];
            ordered |= order[n]>0;
            for(int c=0;c<3;c++){
                base.labels[x][t][a][c]=labels[3*n+c];
                q.labels[x][t][a][c]={labels[3*n+c],direction[3*n+c]};
                moving |= direction[3*n+c]!=0;
            }
        }
        Geometry g{d,{d.field_x.v[i],d.field_t.v[j],d.field_a.v[k]},
                     {d.field_x.v[i+1]-d.field_x.v[i],d.field_t.v[j+1]-d.field_t.v[j],
                      d.field_a.v[k+1]-d.field_a.v[k]}};
        euler::Orientation proof;
        int rc=euler::orientation(base,g.width,proof);
        if(!rc && ordered && moving){
            rc=euler::integrate_cell(g,q,ids,rule,weights,weld,true,cuts,cut_work,split_weld,
                [](const auto&,const auto&,Jet,Jet,double,double){},[](int,Jet){},
                [&](int n,Jet contribution){product_out[n]+=contribution.d;});
        }
        if(rc){problem[0]=i;problem[1]=j;problem[2]=k;return rc;}
    }
    return 0;
}
