// Independent attribution of core work to core-only and weld-crossing rays.
extern "C" int core_cutsum(const Data *source,const Real *labels,const Real *order,
                          Grid gauss,const Real *weights,Real *out) {
    const Data& d=*source;
    const Grid axes[]={d.field_x,d.field_t,d.field_a};
    const Real weld=d.x.v[(d.x.n-1)*2/3];
    Data legacy=d;
    legacy.collar_coeff=nullptr;
    for(int i=0;i<axes[0].n-1;i++)for(int j=0;j<axes[1].n-1;j++)for(int k=0;k<axes[2].n-1;k++){
        Geometry g{d,{},{}};int idx[]={i,j,k};
        for(int h=0;h<3;h++){
            g.lo[h]=axes[h].v[idx[h]];
            g.width[h]=axes[h].v[idx[h]+1]-g.lo[h];
        }
        Q1 q{};bool active=false;
        for(int X=0;X<2;X++)for(int T=0;T<2;T++)for(int A=0;A<2;A++){
            int n=((i+X)*axes[1].n+j+T)*axes[2].n+k+A;
            q.order[X][T][A]=order[n];active|=order[n]>0;
            for(int h=0;h<3;h++)q.labels[X][T][A][h]=labels[3*n+h];
        }
        if(!active)continue;
        for(int T=0;T<gauss.n;T++)for(int A=0;A<gauss.n;A++){
            const Real t=gauss.v[T],a=gauss.v[A];
            auto l=q.at(0,t,a),r=q.at(1,t,a);
            if(std::min(l.xi[0],r.xi[0])>=weld)continue;
            const int group=std::max(l.xi[0],r.xi[0])>weld;
            Real bank[2]={}, counts[2]={}; euler::Cuts cuts;
            for(int which=0;which<2;which++){
                euler::radial_cuts(cuts,which?d:legacy,l.xi,r.xi,weld);
                for(size_t n=1;n<cuts.size();n++){
                    Real left=cuts[n-1].u, width=cuts[n].u-left;
                    if(!(width>0))continue;
                    if(l.xi[0]+(left+width/2)*(r.xi[0]-l.xi[0])>weld)continue;
                    counts[which]++;
                    for(int X=0;X<gauss.n;X++){
                        const Real u=left+width*gauss.v[X];
                        auto p=q.at(u,t,a);
                        if(p.xi[0]>weld)return 2;
                        auto work=euler::material(g,{u,t,a},p,false);
                        bank[which]+=weights[X]*weights[T]*weights[A]*width*p.chi*work.unit;
                    }
                }
            }
            Real *row=out+group*7;
            row[0]+=1;row[1]+=counts[0];row[2]+=counts[1];row[3]+=bank[0];row[4]+=bank[1];
            row[5]+=bank[1]-bank[0];row[6]+=std::abs(bank[1]-bank[0]);
        }
    }
    return 0;
}
