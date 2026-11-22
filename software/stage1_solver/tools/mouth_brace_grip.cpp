// Overdamped open-normal neutral lattice. Python supplies Wiener increments,
// conversion history and the periodic image list of the archived OpenCell.
#include <algorithm>
#include <cmath>
#include <vector>

struct BondImage { int i,j; double x,y; };

extern "C" int evolve(int n, int pairs, const int *ii, const int *jj,
    const double *shift, const double *box, double cutoff, int power,
    const double *within, const double *between, const double *zref,
    const double *normal_load, double brace, double mobility, int hold_x,
    int steps, double dt, const double *order, const double *noise, double *pos) {
    std::vector<BondImage> bonds;
    std::vector<double> old(3*n), force(3*n);
    const double skin=.35, search2=(cutoff+skin)*(cutoff+skin);
    for (int k=0;k<steps;++k) {
        bool rebuild=k==0;
        for (int i=0;!rebuild && i<n;++i) {
            double d2=0;
            for(int a=0;a<3;++a) {double d=pos[3*i+a]-old[3*i+a]; d2+=d*d;}
            rebuild=d2>=.249*skin*skin;
        }
        if (rebuild) {
            std::copy(pos,pos+3*n,old.begin());
            bonds.clear();
            for(int p=0;p<pairs;++p) {
                int i=ii[p],j=jj[p];
                if (i==j) continue;
                double x=pos[3*j]-pos[3*i], y=pos[3*j+1]-pos[3*i+1];
                double dx=shift[2*p]-box[0]*std::nearbyint(x/box[0]);
                double dy=shift[2*p+1]-box[1]*std::nearbyint(y/box[1]);
                double z=pos[3*j+2]-pos[3*i+2];
                if ((x+dx)*(x+dx)+(y+dy)*(y+dy)+z*z<search2)
                    bonds.push_back({i,j,dx,dy});
            }
        }
        std::fill(force.begin(),force.end(),0.);
        for(const auto &b:bonds) {
            int i=3*b.i,j=3*b.j;
            double x=pos[j]-pos[i]+b.x, y=pos[j+1]-pos[i+1]+b.y, z=pos[j+2]-pos[i+2];
            double rr=x*x+y*y+z*z;
            if(rr>=cutoff*cutoff) continue;
            if(!(rr>1.e-12)) return 1;
            double r=std::sqrt(rr), t=std::min(1.,std::abs(z));
            double f=t*t*t*(10.+t*(-15.+6.*t));
            double df=std::copysign(30.*t*t*(1.-t)*(1.-t),z);
            if(power==2) {df*=2.*f; f*=f;}
            double u=std::exp(-within[1]*(r-within[2])), v=std::exp(-between[1]*(r-between[2]));
            double ea=within[0]*u*(u-2.), eb=between[0]*v*(v-2.);
            double ca=2.*within[0]*within[1]*u*(1.-u), cb=2.*between[0]*between[1]*v*(1.-v);
            double radial=(ca+f*(cb-ca))/r;
            double gradient[]={radial*x,radial*y,radial*z+df*(eb-ea)};
            for(int a=0;a<3;++a) {force[i+a]+=gradient[a]; force[j+a]-=gradient[a];}
        }
        double mean_x[2]={0.,0.};
        for(int i=0;i<n;++i) {
            force[3*i+2]+=order[k]*(normal_load[i]-brace*(pos[3*i+2]-zref[i]));
            for(int a=0;a<3;++a) {
                const int l=3*i+a;
                double delta=dt*mobility*force[l]+noise[k*3*n+l];
                if(!std::isfinite(delta) || std::abs(delta)>1.) return 2;
                pos[l]+=delta;
                if(a==0) mean_x[i/(n/2)]+=delta/(n/2);
            }
        }
        if(hold_x) for(int i=0;i<n;++i) pos[3*i]-=mean_x[i/(n/2)];
    }
    return 0;
}
