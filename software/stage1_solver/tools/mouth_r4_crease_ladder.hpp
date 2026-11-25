// Signed work along a prescribed label direction at fixed physical mesh and order.
// Included after the production physical evaluator promoted to long double.
extern "C" int crease_ray_work(const Data *data, const Real *labels, const Real *direction,
        const Real *order, Grid rule, const Real *weights, Real alpha, Real *per_ray) {
    using namespace euler;
    const Data &d=*data;
    const int nx=d.field_x.n,nt=d.field_t.n,na=d.field_a.n;
    const Real weld=d.x.v[(d.x.n-1)*2/3];
    const Grid ax[]={d.field_x,d.field_t,d.field_a};
    Cuts cuts;
    for (int i=0;i<nx-1;++i) for (int j=0;j<nt-1;++j) for (int k=0;k<na-1;++k) {
        const int idx=(i*(nt-1)+j)*(na-1)+k;
        Q1 q{}, velocity{};
        bool active=false;
        for(int x=0;x<2;++x) for(int t=0;t<2;++t) for(int a=0;a<2;++a){
            const int n=((i+x)*nt+j+t)*na+k+a;
            active|=order[n]>0;
            q.order[x][t][a]=order[n];
            for(int c=0;c<3;++c){
                velocity.labels[x][t][a][c]=direction[3*n+c];
                q.labels[x][t][a][c]=labels[3*n+c]+alpha*direction[3*n+c];
            }
        }
        if(!active)continue;
        Geometry g{d,{ax[0].v[i],ax[1].v[j],ax[2].v[k]},
                     {ax[0].v[i+1]-ax[0].v[i],ax[1].v[j+1]-ax[1].v[j],ax[2].v[k+1]-ax[2].v[k]}};
        for(int T=0;T<rule.n;++T) for(int A=0;A<rule.n;++A){
            Real t=rule.v[T],a=rule.v[A];
            auto l=q.at(0,t,a),r=q.at(1,t,a);
            auto vl=velocity.at(0,t,a),vr=velocity.at(1,t,a);
            Vec delta{}; Mat dD{};
            for(int c=0;c<3;++c){
                delta[c]=r.xi[c]-l.xi[c];
                for(int j=0;j<3;++j)dD[c][j]=r.D[c][j]-l.D[c][j];
            }
            const Real dc=r.chi-l.chi;
            radial_cuts(cuts,d,l.xi,r.xi,weld);
            std::vector<Real> speed(cuts.size(),0);
            for(size_t h=1;h+1<cuts.size();++h){
                const auto cut=cuts[h];
                const Real denominator=dot(cut.normal,delta);
                const Real numerator=(1-cut.u)*dot(cut.normal,vl.xi)+cut.u*dot(cut.normal,vr.xi);
                if(denominator==0 && numerator!=0) return 1;
                if(denominator!=0) speed[h]=-numerator/denominator;
            }
            Real *row=per_ray+3*((idx*rule.n+T)*rule.n+A);
            for(size_t s=1;s<cuts.size();++s){
                const Real span=cuts[s].u-cuts[s-1].u,dspan=speed[s]-speed[s-1];
                for(int X=0;X<rule.n;++X){
                    Real f=rule.v[X],u=cuts[s-1].u+f*span,du=speed[s-1]+f*dspan;
                    Real w=weights[T]*weights[A]*weights[X];
                    auto p=l;
                    p.chi+=u*dc;
                    for(int c=0;c<3;++c){
                        p.xi[c]+=u*delta[c];
                        for(int j=0;j<3;++j)p.D[c][j]+=u*dD[c][j];
                    }
                    PointWork m;
                    try {m=material(g,{u,t,a},p,true);} catch(const std::domain_error&){return 2;}
                    Real ordinary=0,integrand_u=m.flux+m.unit*dc+dot(m.value,delta);
                    for(int c=0;c<3;++c){
                        ordinary+=m.value[c]*(vl.xi[c]+u*(vr.xi[c]-vl.xi[c]));
                        for(int j=0;j<3;++j)
                            ordinary+=m.tangent[c][j]*(vl.D[c][j]+u*(vr.D[c][j]-vl.D[c][j]));
                        integrand_u+=dot(m.tangent[c],dD[c]);
                    }
                    row[0]+=w*span*p.chi*m.unit;
                    row[1]+=w*span*ordinary;
                    row[2]+=w*(dspan*p.chi*m.unit+span*du*integrand_u);
                }
            }
        }
    }
    return 0;
}
