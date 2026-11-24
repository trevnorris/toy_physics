// Locate positive-order witnesses in every crystal prism, using the production
// reference-volume factors. The Python caller independently certifies a selected witness.
#include <mouth_free_shape_solid.cpp>

extern "C" void ordered_j_candidates(const Data *data, double *out) {
    const Data &d=*data;
    // Return the best finite sample in each prism; columns are J, b1,b2,y.
    for (int j=0;j<d.t.n-1;++j) for(int i=0;i<d.x.n-1;++i)
    for (int k=0;k<d.a.n-1;++k) for(int tri=0;tri<2;++tri,out+=4) {
        out[0]=std::numeric_limits<double>::infinity();
        const double sign=tri ? -1. : 1.;
        auto probe=[&](double b1,double b2,double y) {
            if (b1<0. || b2<0. || b1+b2>1. || y<0. || y>1.) return;
            double dx=d.x.v[i+1]-d.x.v[i],dt=d.t.v[j+1]-d.t.v[j],da=d.a.v[k+1]-d.a.v[k];
            double x=d.x.v[i+tri]+sign*b1*dx, t=d.t.v[j+tri]+sign*b2*dt, a=d.a.v[k]+y*da;
            if (x<=0.) return; // the physical center uses a radial limit
            // Build firmly inside the requested triangle and extrapolate its
            // affine x polynomial to the one-sided measurement position.
            double anchor=d.x.v[i]+dx*(1.-(t-d.t.v[j])/dt+(tri?.125:-.125));
            const TriSlice slice(d,i,j,k,t,a,anchor);
            Vec X=slice.X; Mat M=slice.M;
            for(int c=0;c<3;++c) {
                X[c]+=M[c][0]*(x-anchor); M[c][2]+=slice.Max[c]*(x-anchor);
            }
            const Basis fx=basis(d.field_x,X[0]),ft=basis(d.field_t,X[1]),fa=basis(d.field_a,X[2]);
            double chi=0.;
            for(int xi=0;xi<2;++xi)for(int ti=0;ti<2;++ti)for(int ai=0;ai<2;++ai)
                chi+=fx.w[xi]*ft.w[ti]*fa.w[ai]*d.correction[field_id(d,fx.i+xi,fa.i+ai,ft.i+ti)];
            if (!(chi>0.)) return;
            const Ref ref=reference(d,x,t,profile(d,t,true),i<(d.x.n-1)*2/3);
            const Map wall=radius(d,X);
            const double angular=a==0. ? M[2][2] : std::sin(X[2])/std::sin(a);
            double J=wall.grad[0]*d.length*std::pow(wall.R/ref.S,2)*angular/ref.det2*determinant(M);
            if (J>0. && J<out[0]) {out[0]=J;out[1]=b1;out[2]=b2;out[3]=y;}
        };
        for(int n=0;n<=3;++n) for(int m=0;m<=3-n;++m) for(int a=0;a<=4;++a)
            probe(n/3.,m/3.,a/4.);
        if (!std::isfinite(out[0])) continue;
        // Small determinant cells receive a dyadic search along the ordered
        // side of the spatial support boundary, including diagonal approaches.
        if(out[0]<.25) {
            for(double h=.125; h>1.e-15; h*=.5) {
                double b1=out[1],b2=out[2],a=out[3];
                for(int b=-1;b<=1;++b)for(int c=-1;c<=1;++c)for(int e=-1;e<=1;++e)
                    probe(b1+h*b,b2+h*c,a+h*e);
            }
        }
    }
}
