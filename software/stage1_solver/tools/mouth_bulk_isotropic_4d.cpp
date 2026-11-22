// Four-periodic-dimensional Morse material; include all images inside a physical cutoff.
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using Vec = std::array<double,4>;
struct Bond { int i, j; Vec shift; };
struct Cell {
    int n;
    Vec length;
    double cutoff, depth, alpha, rest;
    double skin = .28;
    std::vector<Vec> old;
    std::vector<Bond> pairs;
    Cell(int n_, const double* box, double rc, const double* p) :
        n(n_),cutoff(rc),depth(p[0]),alpha(p[1]),rest(p[2]) {
        std::copy(box,box+4,length.begin());
    }
    void neighbors(const double* x) {
        bool rebuild = old.empty();
        for(int i=0; !rebuild && i<n; ++i) {
            double ds=0;
            for(int k=0;k<4;++k){ double d=x[4*i+k]-old[i][k]; ds+=d*d; }
            rebuild=ds>=.249*skin*skin;
        }
        if(!rebuild) return;
        old.resize(n);
        for(int i=0;i<n;++i) std::copy(x+4*i,x+4*i+4,old[i].begin());
        pairs.clear();
        double reach=cutoff+skin;
        for(int i=0;i<n;++i) for(int j=i;j<n;++j) {
            Vec raw; std::array<int,4> lo,hi;
            for(int k=0;k<4;++k) {
                raw[k]=x[4*j+k]-x[4*i+k];
                lo[k]=int(std::ceil((-reach-raw[k])/length[k]));
                hi[k]=int(std::floor((reach-raw[k])/length[k]));
            }
            for(int a=lo[0];a<=hi[0];++a) for(int b=lo[1];b<=hi[1];++b)
            for(int c=lo[2];c<=hi[2];++c) for(int d=lo[3];d<=hi[3];++d) {
                // For periodic copies of a single particle count one of the opposite pair.
                if(i==j && std::array<int,4>{a,b,c,d}<=std::array<int,4>{}) continue;
                Vec shift{a*length[0],b*length[1],c*length[2],d*length[3]};
                double r2=0;
                for(int k=0;k<4;++k){double u=raw[k]+shift[k];r2+=u*u;}
                if(r2<reach*reach) pairs.push_back({i,j,shift});
            }
        }
    }
    // Use a C2 switch in the last quarter bond length. The same potential,
    // including its derivative through the switch, supplies energy and force.
    double evaluate(const double* x,double* grad,double* diag) {
        neighbors(x);
        std::fill(grad,grad+4*n,0.);
        if(diag) std::fill(diag,diag+4,0.);
        double total=0.;
        for(const auto& bond:pairs) {
            Vec v; double r2=0.;
            for(int k=0;k<4;++k){v[k]=x[4*bond.j+k]-x[4*bond.i+k]+bond.shift[k];r2+=v[k]*v[k];}
            if(r2>=cutoff*cutoff) continue;
            double r=std::sqrt(r2), w=std::exp(-alpha*(r-rest));
            double energy=depth*w*(w-2.),du=2.*alpha*depth*w*(1.-w);
            if(r>cutoff-.25) {
                double t=(cutoff-r)/.25,s=t*t*t*(10.+t*(-15.+6.*t));
                du=du*s-energy*120.*t*t*(1.-t)*(1.-t);
                energy*=s;
            }
            total+=energy;
            double g=du/r;
            for(int k=0;k<4;++k) {
                double force=g*v[k];
                grad[4*bond.i+k]-=force;
                grad[4*bond.j+k]+=force;
                if(diag) diag[k]+=force*v[k];
            }
        }
        return total;
    }
};
extern "C" {
void* make_cell(int n,const double* box,double cutoff,const double* params) {
    return new Cell(n,box,cutoff,params);
}
void free_cell(void* cell){delete static_cast<Cell*>(cell);}
double evaluate(void* cell,const double* x,double* grad,double* diag){
    return static_cast<Cell*>(cell)->evaluate(x,grad,diag);
}
int integrate(void* ptr,int steps,double dt,const double* noise,double* x){
    Cell& cell=*static_cast<Cell*>(ptr);
    int count=4*cell.n;
    std::vector<double> grad(count);
    for(int j=0;j<steps;++j){
        cell.evaluate(x,grad.data(),nullptr);
        for(int k=0;k<count;++k){
            double dx=noise[j*count+k]-dt*grad[k];
            if(!std::isfinite(dx) || std::abs(dx)>1.) return 1;
            x[k]+=dx;
        }
    }
    return 0;
}
}
