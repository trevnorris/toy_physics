// Read-only cell energies for central differences of every free label at fixed geometry.
// Retains the evaluator's moving reference cuts. No analytic label gradient is used here.
#include "mouth_r4_eulerian_solid.cpp"

static long double crystal_energy(const Geometry &g,const Q1 &q,Grid rule,const double *weights){
    euler::Cuts cuts;
    const auto &d=g.d;
    const double weld=d.x.v[(d.x.n-1)*2/3];
    long double result=0;
    for(int t=0;t<rule.n;t++)for(int a=0;a<rule.n;a++){
        const double v=rule.v[t],theta=rule.v[a];
        const auto left=q.at(0,v,theta),right=q.at(1,v,theta);
        euler::radial_cuts(cuts,d,left.xi,right.xi,weld);
        for(size_t s=1;s<cuts.size();s++){
            const double width=cuts[s].u-cuts[s-1].u;
            for(int ix=0;ix<rule.n;ix++){
                const double u=cuts[s-1].u+width*rule.v[ix];
                auto p=left;
                p.chi+=u*(right.chi-left.chi);
                for(int c=0;c<3;c++){
                    p.xi[c]+=u*(right.xi[c]-left.xi[c]);
                    for(int k=0;k<3;k++)p.D[c][k]+=u*(right.D[c][k]-left.D[c][k]);
                }
                const auto state=euler::material(g,{u,v,theta},p,false);
                result+=weights[t]*weights[a]*weights[ix]*width*p.chi*state.unit;
            }
        }
    }
    return result;
}

static long double gauge_energy(const Geometry &g,const Q1 &q){
    long double harmonic=0,volume=0,barrier=0;
    const double half=.5/std::sqrt(3.);
    Q1 offsets=q;
    for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
        int bit[]={i,j,k};
        for(int c=0;c<3;c++) offsets.labels[i][j][k][c]-=g.lo[c]+bit[c]*g.width[c];
    }
    for(double i:{.5-half,.5+half})for(double j:{.5-half,.5+half})for(double k:{.5-half,.5+half}){
        Vec frac{i,j,k},pos;
        for(int c=0;c<3;c++)pos[c]=g.lo[c]+g.width[c]*frac[c];
        const Map r=radius(g.d,pos);
        const Mat C=physical_basis(r,g.d.length,g.width);
        const long double dv=.5L*pi*r.R*std::sin(pos[2])*determinant(C);
        volume+=dv;
        harmonic+=.5L*dv*norm2(product(offsets.at(i,j,k).D,invert(C)));
    }
    for(int i=0;i<2;i++)for(int j=0;j<2;j++)for(int k=0;k<2;k++){
        long double m[3][3];
        for(int c=0;c<3;c++){
            m[c][0]=(static_cast<long double>(q.labels[1][j][k][c])-q.labels[0][j][k][c])/g.width[0];
            m[c][1]=(static_cast<long double>(q.labels[i][1][k][c])-q.labels[i][0][k][c])/g.width[1];
            m[c][2]=(static_cast<long double>(q.labels[i][j][1][c])-q.labels[i][j][0][c])/g.width[2];
        }
        long double d=0;
        for(int c=0;c<3;c++)d+=m[0][c]*(m[1][(c+1)%3]*m[2][(c+2)%3]-m[2][(c+1)%3]*m[1][(c+2)%3]);
        if(!(d>0))throw std::domain_error("melt corner");
        barrier+=d-1-std::log(d);
    }
    return harmonic+volume*barrier/8;
}

extern "C" int r4_cell_work(const Data *data,const double *labels,const double *order,Grid rule,
    const double *weights,const int *enabled,long double *values,double *min_orientation){
    const auto &d=*data;const Grid axes[]={d.field_x,d.field_t,d.field_a};
    const int nt=d.field_t.n,na=d.field_a.n;
    *min_orientation=std::numeric_limits<double>::infinity();
    for(int i=0;i<axes[0].n-1;i++)for(int j=0;j<nt-1;j++)for(int k=0;k<na-1;k++){
        const int cell=(i*(nt-1)+j)*(na-1)+k;
        if(!enabled[cell])continue;
        int index[]={i,j,k};Geometry g{d,{},{}};
        for(int c=0;c<3;c++){g.lo[c]=axes[c].v[index[c]];g.width[c]=axes[c].v[index[c]+1]-g.lo[c];}
        Q1 q{};bool active=false;
        for(int x=0;x<2;x++)for(int t=0;t<2;t++)for(int a=0;a<2;a++){
            int id=((i+x)*nt+j+t)*na+k+a;
            q.order[x][t][a]=order[id];active|=order[id]>0;
            for(int c=0;c<3;c++)q.labels[x][t][a][c]=labels[3*id+c];
        }
        euler::Orientation proof;
        int rc=euler::orientation(q,g.width,proof);
        if(rc)return rc;
        *min_orientation=std::min(*min_orientation,proof.lower);
        try{values[cell]=active?crystal_energy(g,q,rule,weights):gauge_energy(g,q);}
        catch(std::domain_error&){return 4;}
    }
    return 0;
}
