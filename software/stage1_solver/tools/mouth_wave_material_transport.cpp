// Axial variable-section material wave with conserved momentum and displacement
// gradient. Read cell masses, axial stiffnesses, depths and the common kick.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

struct Pair { double p, g; };
Pair operator+(Pair a, Pair b) { return {a.p+b.p, a.g+b.g}; }
Pair operator-(Pair a, Pair b) { return {a.p-b.p, a.g-b.g}; }
Pair operator*(double k, Pair a) { return {k*a.p, k*a.g}; }
double limited(double a, double b) {
    if (a*b<=0) return 0;
    return std::copysign(std::min({2*std::abs(a),2*std::abs(b), std::abs(a+b)*.5}),a);
}

struct Tube {
    int n;
    double dx;
    std::vector<double> m,k,u,c;
    std::vector<Pair> left,right,slope,flux;
    explicit Tube(int count):n(count),m(n+4),k(n+4),u(n+4),c(n+4),
        left(n+4),right(n+4),slope(n+4),flux(n+4) {}
    Pair physical(Pair a, int i) const {
        return {u[i]*a.p-k[i]*a.g, u[i]*a.g-a.p/m[i]};
    }
    double energy(Pair a,int i) const {return .5*(a.p*a.p/m[i]+k[i]*a.g*a.g);}
    double power(Pair a,int i) const {return u[i]*energy(a,i)-k[i]*a.g*a.p/m[i];}
    struct Balance {double feed,bulk,work;};
    Balance derivative(std::vector<Pair>& y, std::vector<Pair>& dy) {
        // In the uniform entrance lead retain only the outgoing leftgoing
        // characteristic. In the shear-free right-hand lead both are outgoing.
        const double z=m[2]*c[2];
        const double vleft=.5*(y[2].p/m[2]+c[2]*y[2].g);
        y[0]=y[1]={m[2]*vleft, m[2]*vleft/z};
        y[n+2]=y[n+3]=y[n+1];
        for(int i=1;i<=n+2;i++) {
            Pair a=y[i]-y[i-1], b=y[i+1]-y[i];
            slope[i]={limited(a.p,b.p), limited(a.g,b.g)};
            left[i]=y[i]-.5*slope[i];
            right[i]=y[i]+.5*slope[i];
        }
        for(int i=1;i<=n+1;i++) {
            double speed=std::max(u[i]+c[i],u[i+1]+c[i+1]);
            flux[i]=.5*(physical(right[i],i)+physical(left[i+1],i+1)-
                          speed*(left[i+1]-right[i]));
        }
        double work=0;
        for(int i=2;i<=n+1;i++) {
            dy[i]=(1/dx)*(flux[i-1]-flux[i]);
            // mU = Q, so there is no background work on the kinetic term.
            double gradient=(u[i]*(k[i+1]-k[i-1])-k[i]*(u[i+1]-u[i-1]))/(2*dx);
            work+=.5*gradient*y[i].g*y[i].g*dx;
        }
        return {-power(y[1],1),power(y[n+1],n+1),work};
    }
};

int main() {
    int n; double dz,q;
    if (!(std::cin>>n>>dz>>q) || n<20) return 2;
    Tube tube(n); tube.dx=dz;
    std::vector<Pair> y(n+4), pred(n+4),dy(n+4);
    double e0=0,top_speed=0;
    for(int i=2;i<=n+1;i++) {
        double v;
        if(!(std::cin>>tube.m[i]>>tube.k[i]>>v)) return 2;
        tube.u[i]=q/tube.m[i]; tube.c[i]=std::sqrt(tube.k[i]/tube.m[i]);
        top_speed=std::max(top_speed,tube.u[i]+tube.c[i]);
        y[i]={tube.m[i]*v,0};
        e0+=dz*tube.energy(y[i],i);
    }
    for(auto field:{&tube.m,&tube.k,&tube.u,&tube.c}) {
        (*field)[0]=(*field)[1]=(*field)[2];
        (*field)[n+2]=(*field)[n+3]=(*field)[n+1];
    }
    double feed=0,bulk=0,work=0,time=0;
    std::cout<<std::setprecision(16);
    for(double target:{2.,4.,8.,16.,48.,96.,160.}) {
        int steps=std::ceil((target-time)*top_speed/(.62*dz));
        double dt=(target-time)/steps;
        for(int s=0;s<steps;s++) {
            auto a=tube.derivative(y,dy);
            for(int i=2;i<=n+1;i++) pred[i]=y[i]+dt*dy[i];
            auto b=tube.derivative(pred,dy);
            for(int i=2;i<=n+1;i++) y[i]=.5*(y[i]+pred[i]+dt*dy[i]);
            feed+=.5*dt*(a.feed+b.feed);
            bulk+=.5*dt*(a.bulk+b.bulk);
            work+=.5*dt*(a.work+b.work);
        }
        time=target;
        double remaining=0;
        for(int i=2;i<=n+1;i++) remaining+=dz*tube.energy(y[i],i);
        double defect=(remaining+feed+bulk-work)/e0-1;
        std::cout<<target<<" "<<remaining/e0<<" "<<feed/e0<<" "<<bulk/e0<<" "
                 <<work/e0<<" "<<defect<<"\n";
        if(!std::isfinite(defect)) return 1;
    }
}
