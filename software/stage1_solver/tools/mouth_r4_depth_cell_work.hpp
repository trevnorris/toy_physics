// Read-only log-depth tangent of the native radial-cut quadrature.
// Included after the promoted physical material evaluator. Entries per spatial
// cell are the moving physical frame, transported labels and transported cuts.
extern "C" int depth_cell_gradient(const Data *dp, const Real *labels, const Real *order,
                                  Grid rule, const Real *weights, const Real *rates, Real *terms) {
    using namespace euler;
    const Data &d=*dp;
    const Grid axes[]={d.field_x,d.field_t,d.field_a};
    const int nx=axes[0].n,nt=axes[1].n,na=axes[2].n;
    const Real weld=d.x.v[(d.x.n-1)*2/3];
    Cuts cuts;
    for(int i=0;i<nx-1;i++)for(int j=0;j<nt-1;j++)for(int k=0;k<na-1;k++){
        Q1 q{};
        bool active=false;
        for(int x=0;x<2;x++)for(int t=0;t<2;t++)for(int a=0;a<2;a++){
            int n=((i+x)*nt+j+t)*na+k+a;
            q.order[x][t][a]=order[n];
            active|=order[n]>0;
            for(int c=0;c<3;c++)q.labels[x][t][a][c]=labels[3*n+c];
        }
        if(!active)continue;
        const int cell=(i*(nt-1)+j)*(na-1)+k;
        Geometry g{d,{axes[0].v[i],axes[1].v[j],axes[2].v[k]},
                     {axes[0].v[i+1]-axes[0].v[i],axes[1].v[j+1]-axes[1].v[j],axes[2].v[k+1]-axes[2].v[k]}};
        const Real dwidth=rates[j+1]-rates[j];
        for(int T=0;T<rule.n;T++)for(int A=0;A<rule.n;A++){
            Real t=rule.v[T],a=rule.v[A],rate=rates[j]+t*dwidth;
            auto l=q.at(0,t,a),r=q.at(1,t,a);
            Vec delta{};Mat dD{};
            for(int c=0;c<3;c++){
                delta[c]=r.xi[c]-l.xi[c];
                for(int s=0;s<3;s++)dD[c][s]=r.D[c][s]-l.D[c][s];
            }
            const Real dc=r.chi-l.chi;
            radial_cuts(cuts,d,l.xi,r.xi,weld);
            std::vector<Real> velocities(cuts.size());
            for(size_t c=1;c+1<cuts.size();c++){
                const auto &n=cuts[c].normal;
                const Real denom=dot(n,delta);
                if(denom==0){if(n[1]*rate!=0)return 1;}
                else velocities[c]=-n[1]*rate/denom;
            }
            for(size_t s=1;s<cuts.size();s++){
                const Real width=cuts[s].u-cuts[s-1].u, dspan=velocities[s]-velocities[s-1];
                for(int X=0;X<rule.n;X++){
                    const Real x=rule.v[X], u=cuts[s-1].u+x*width;
                    const Real speed=velocities[s-1]+x*dspan,w=weights[T]*weights[A]*weights[X];
                    auto p=l;
                    p.chi+=u*dc;
                    for(int c=0;c<3;c++){
                        p.xi[c]+=u*delta[c];
                        for(int f=0;f<3;f++)p.D[c][f]+=u*dD[c][f];
                    }
                    auto m=material(g,{u,t,a},p,true);
                    // At stationary labels the reference-volume measure in chi W is fixed.
                    const Ref ref=reference(d,p.xi[0],p.xi[1],profile(d,p.xi[1],true),p.xi[0]<=weld);
                    const Real volume=4*pi*determinant(p.D)*ref.det2*ref.S*ref.S*std::sin(p.xi[2]);
                    Vec pos={g.lo[0]+u*g.width[0],g.lo[1]+t*g.width[1],g.lo[2]+a*g.width[2]};
                    const auto rwall=radius(d,pos);
                    const Real dr=rwall.grad[1]*rate;
                    const Real hoop=rwall.R*std::sin(pos[2])/(ref.S*std::sin(p.xi[2]));
                    Mat C=physical_basis(rwall,d.length,g.width),dC{};
                    for(int c=0;c<3;c++)dC[0][c]=rwall.hess[c][1]*rate*g.width[c];
                    dC[0][1]+=rwall.grad[1]*dwidth;
                    dC[1][1]=d.length*(g.width[1]+dwidth);
                    dC[2][2]=dr*g.width[2];
                    const Mat right=product(invert(p.D),ref.inv),F=product(C,right),dF=product(dC,right);
                    const Real J=1/m.Jinv;
                    const Real divergence=rwall.hess[0][1]*rate/rwall.grad[0]+1+2*dr/rwall.R+dwidth/g.width[1];
                    Real geom=hoop*hoop*dr/rwall.R+(-1+Real(.5)*J*(J-1))*divergence;
                    for(int c=0;c<3;c++)geom+=dot(F[c],dF[c]);
                    terms[3*cell]+=w*width*p.chi*volume*geom;
                    terms[3*cell+1]+=w*width*(m.value[1]*rate+m.tangent[1][1]*dwidth);
                    Real integrand_u=m.flux+dot(m.value,delta)+m.unit*dc;
                    for(int c=0;c<3;c++)integrand_u+=dot(m.tangent[c],dD[c]);
                    terms[3*cell+2]+=w*(dspan*p.chi*m.unit+width*speed*integrand_u);
                }
            }
        }
    }
    return 0;
}
