// Directional work on the star of the last free material ring.
// Reuse the held solver's chart, registration, basis and quadrature crossings.
#include <mouth_free_shape_solid.cpp>

struct Work {
    double stretch{}, log_volume{}, compression{}, sample_base{}, sample_correction{}, moving_rule{};
};

double contract(const Mat &a, const Mat &b) {
    return dot(a[0],b[0])+dot(a[1],b[1])+dot(a[2],b[2]);
}

double radial_lift(const TriSlice &s, const double *lift, double x) {
    double delta = 0.;
    for(int n=0;n<6;++n) delta += lift[s.id[n]]*(s.n[n]+(x-s.x0)*s.nx[n]);
    return delta;
}

// Differentiate at fixed reference material position. The quadrature-point and
// weight motion from the discrete cut rule is integrated separately below.
void local_work(const Data &d, const TriSlice &s, const Profile &initial, const double *lift,
                double t, double a, double x, double weight, bool core, Work &out) {
    const double dx=x-s.x0, v=radial_lift(s,lift,x);
    Vec X=s.X; Mat M=s.M, dM{};
    for(int c=0;c<3;++c) { X[c]+=M[c][0]*dx; M[c][2]+=s.Max[c]*dx; }
    for(int n=0;n<6;++n) {
        const double u=lift[s.id[n]];
        dM[0][0]+=u*s.nx[n]; dM[0][1]+=u*s.nt[n]; dM[0][2]+=u*(s.na[n]+dx*s.nax[n]);
    }
    const auto ref=reference(d,x,t,initial,core);
    const auto wall=radius(d,X);
    Mat C={wall.grad,Vec{0.,d.length,0.},Vec{0.,0.,wall.R}}, dC{};
    for(int c=0;c<3;++c) dC[0][c]=wall.hess[0][c]*v;
    const double dR=wall.grad[0]*v;
    dC[2][2]=dR;
    const Mat B=product(M,ref.inv), F=product(C,B);
    const double hoop=wall.R*std::sin(X[2])/(ref.S*std::sin(a)), dhoop=hoop*dR/wall.R;
    const double J=wall.grad[0]*d.length*wall.R*determinant(M)*hoop/(ref.det2*ref.S);
    const double W=.5*(norm2(F)+hoop*hoop-4.)-std::log(J)+.25*(J-1.)*(J-1.);
    const double dlogJ=dC[0][0]/wall.grad[0]+2*dR/wall.R+contract(cofactor(M),dM)/determinant(M);
    const double dStretch=contract(F,product(dC,B))+contract(F,product(C,product(dM,ref.inv)))+hoop*dhoop;
    const auto bx=basis(d.base_x,X[0]), bt=basis(d.base_t,X[1]);
    double chi=0., from_base=0., from_correction=0.;
    for(int j=0;j<2;++j) for(int i=0;i<2;++i) {
        const double u=d.base[(bt.i+j)*d.base_x.n+bx.i+i];
        chi+=bx.w[i]*bt.w[j]*u; from_base+=bx.v[i]*bt.w[j]*u;
    }
    const auto fx=basis(d.field_x,X[0]), ft=basis(d.field_t,X[1]), fa=basis(d.field_a,X[2]);
    for(int i=0;i<2;++i) for(int k=0;k<2;++k) for(int j=0;j<2;++j) {
        const double u=d.correction[((fx.i+i)*d.field_a.n+fa.i+k)*d.field_t.n+ft.i+j];
        chi+=fx.w[i]*fa.w[k]*ft.w[j]*u; from_correction+=fx.v[i]*fa.w[k]*ft.w[j]*u;
    }
    const double wt=weight*4*pi*ref.det2*ref.S*ref.S*std::sin(a);
    out.stretch+=wt*chi*dStretch;
    out.log_volume-=wt*chi*dlogJ;
    out.compression+=wt*chi*.5*J*(J-1.)*dlogJ;
    out.sample_base+=wt*W*from_base*v;
    out.sample_correction+=wt*W*from_correction*v;
}

// slots 0/1 are the actual chi*W energy and its mobile correction on the
// affected cells. Slots 2..7 are dE per unit physical outward ring motion.
extern "C" int radial_ring_work(const Data *data, const double *lift, int differentiate, double *values) {
    const Data &d=*data;
    std::array<double,30> stats{};
    std::vector<double> gradients(differentiate ? 3*d.x.n*d.t.n*d.a.n : 0);
    Results result{differentiate ? gradients.data() : nullptr,nullptr,stats.data()};
    Work work{}; const int cutnode=(d.x.n-1)*2/3;
    std::vector<Edge> edges; std::vector<double> fixed;
    for(int ti=0;ti<d.depth_edges.n-1;++ti) for(int tz=0;tz<d.gauss.n;++tz) {
        const double dt=d.depth_edges.v[ti+1]-d.depth_edges.v[ti], t=d.depth_edges.v[ti]+dt*d.gauss.v[tz];
        const int j=cell(d.t,t);
        if(j<d.t.n-3) continue;
        const double v=(t-d.t.v[j])/(d.t.v[j+1]-d.t.v[j]);
        const int jr=cell(d.reg_t,t); const double vr=(t-d.reg_t.v[jr])/(d.reg_t.v[jr+1]-d.reg_t.v[jr]);
        const auto initial=profile(d,t,true); fixed.clear();
        for(int l=0;l<d.reg_e.n-1;++l) {
            if(d.reg_e.v[l]<-1e-15) continue;
            fixed.push_back(std::exp(-d.pressure*d.reg_e.v[l]));
            fixed.push_back(std::exp(-d.pressure*(d.reg_e.v[l]+(1-vr)*(d.reg_e.v[l+1]-d.reg_e.v[l]))));
        }
        std::sort(fixed.begin(),fixed.end());
        for(int ai=0;ai<d.angle_edges.n-1;++ai) for(int az=0;az<d.gauss.n;++az) {
            const double da=d.angle_edges.v[ai+1]-d.angle_edges.v[ai], a=d.angle_edges.v[ai]+da*d.gauss.v[az];
            const int k=cell(d.a,a); const double w=dt*d.weights[tz]*da*d.weights[az];
            for(int i=d.x.n-3;i<d.x.n-1;++i) {
                const double lo=d.x.v[i], hi=d.x.v[i+1], diagonal=lo+(1-v)*(hi-lo);
                for(int half=0;half<2;++half) {
                    const double l=half?diagonal:lo, r=half?hi:diagonal;
                    if(r-l<1e-15) continue;
                    TriSlice s(d,i,j,k,t,a,(l+r)/2);
                    edges={{l,-1},{r,-1}};
                    if(i>=cutnode) for(auto it=std::upper_bound(fixed.begin(),fixed.end(),l);it!=fixed.end() && *it<r;++it)
                        edges.push_back({*it,-1});
                    cut(s,l,r,d.base_x,0,edges); cut(s,l,r,d.base_t,1,edges);
                    cut(s,l,r,d.cut_x,0,edges); cut(s,l,r,d.cut_t,1,edges); cut(s,l,r,d.cut_a,2,edges);
                    std::sort(edges.begin(),edges.end(),[](Edge p,Edge q){return p.x<q.x;});
                    for(size_t z=0;z+1<edges.size();++z) {
                        const double xl=edges[z].x, delta=edges[z+1].x-xl;
                        if(delta<1e-15) continue;
                        double edge_work[2]{};
                        for(int h=0;h<d.gauss.n;++h) {
                            const double g=d.gauss.v[h], wh=w*d.weights[h], x=xl+g*delta;
                            double density=0., dx_density=0.;
                            const int err=sample(d,result,initial,s,t,a,x,wh*delta,i<cutnode,density,dx_density);
                            if(err) return err;
                            if(differentiate) {
                                local_work(d,s,initial,lift,t,a,x,wh*delta,i<cutnode,work);
                                edge_work[0]+=wh*((1-g)*delta*dx_density-density);
                                edge_work[1]+=wh*(g*delta*dx_density+density);
                            }
                        }
                        if(differentiate) for(int e=0;e<2;++e) {
                            const auto edge=edges[z+e];
                            // Only radial coordinates have a nonzero displacement.
                            if(edge.component==0)
                                work.moving_rule-=edge_work[e]*radial_lift(s,lift,edge.x)/s.M[0][0];
                        }
                    }
                }
            }
        }
    }
    std::array<double,8> output={stats[0],stats[1],work.stretch,work.log_volume,work.compression,
                                 work.sample_base,work.sample_correction,work.moving_rule};
    std::copy(output.begin(),output.end(),values);
    return 0;
}
