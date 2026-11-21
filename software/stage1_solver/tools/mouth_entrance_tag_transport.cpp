// A steady stream of staggered material rows, fed at the fixed geometric waist.
// Arc positions, spacing, hoop stretch and order are samples of a measured tag path.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using Vec=std::vector<double>;

struct Site {double arc,gap,hoop,order,h,contact;};
struct Bond {int i,j; double dx2,base_y,period,ymax,weight;};
struct Phase {std::vector<Site> site; std::vector<Bond> bonds;};
struct Law {
    Vec energy, slope, radial;
    double spacing, range2=9., increment=.0001;
    Law(const char* file,double depth,double alpha,double rest) {
        std::ifstream in(file); double e,f;
        if(!(in>>spacing) || spacing<=0) throw std::runtime_error("Bad wall law");
        while(in>>e>>f){energy.push_back(e);slope.push_back(f);}
        if(energy.size()<2) throw std::runtime_error("No contacts");
        for(int n=0;n<=std::lround(range2/increment);++n) {
            double r=std::sqrt(std::max(n*increment, .01)), v=std::exp(-alpha*(r-rest));
            double pot=depth*v*(v-2.),tension=2*depth*alpha*v*(1-v);
            if(r>2.5) {
                double t=(r-2.5)*2., sw=1-t*t*t*(10+t*(-15+6*t));
                tension=tension*sw-pot*60*t*t*(1-t)*(1-t);
            }
            radial.push_back(tension/r);
        }
    }
    double pair(double rr) const {
        if(rr<.01 || !std::isfinite(rr)) throw std::runtime_error("Overlapping material sites");
        if(rr>=range2) return 0;
        const double z=rr/increment;
        const size_t k=std::min(static_cast<size_t>(z),radial.size()-2);
        return radial[k]+(z-k)*(radial[k+1]-radial[k]);
    }
    std::array<double,2> wall(double u) const {
        double z=std::abs(u)/spacing;
        if(!std::isfinite(z) || z>=energy.size()-1) throw std::runtime_error("Uncalibrated contact displacement");
        size_t k=static_cast<size_t>(z);
        double t=z-k,delta=energy[k+1]-energy[k];
        double b=3*delta-spacing*(2*slope[k]+slope[k+1]);
        double c=spacing*(slope[k]+slope[k+1])-2*delta;
        return {energy[k]+t*(spacing*slope[k]+t*(b+t*c)),std::copysign(1.,u)*
                (slope[k]+t*(2*b+3*t*c)/spacing)};
    }
};
struct Stream {
    std::vector<Site> path;
    std::vector<Phase> phases;
    int rows,lanes,sub,preparation,table_sub;
    double period,sites,entry,depth,alpha,rest,low,high,middle,h;
    explicit Stream(char** a):rows(std::stoi(a[3])),lanes(std::stoi(a[4])),sub(std::stoi(a[5])),
        preparation(std::stoi(a[6])),h(std::stod(a[7])) {
        std::ifstream in(a[1]);int num;
        if(!(in>>period>>sites>>entry>>depth>>alpha>>rest>>low>>high>>middle>>table_sub>>num))
            throw std::runtime_error("Bad path header");
        Site s;
        while(in>>s.arc>>s.gap>>s.hoop>>s.order>>s.h>>s.contact) path.push_back(s);
        if((int)path.size()!=num || table_sub%(2*sub) || rows*table_sub+2*table_sub>=num
            || lanes<2 || lanes%2 || sub%8 || h<=0 || path[(rows-5)*table_sub].arc<high)
            throw std::runtime_error("Insufficient spatial support, phases or lanes");
        double maxdiff=0;
        for(int j=0;j<=2*sub;++j) {
            Phase p;
            for(int k=0;k<rows*lanes;++k) {
                int lane=k%lanes, ix=(k/lanes)*table_sub+(lane%2)*(table_sub/2)+j*(table_sub/(2*sub));
                p.site.push_back(path[ix]);
            }
            for(int i=0;i<count();++i) for(int k=i+1;k<count();++k) {
                if((k/lanes-i/lanes)>1 && p.site[k].arc-p.site[i].arc>4.)break;
                const auto& u=p.site[i];const auto& v=p.site[k];
                const double x=v.arc-u.arc;
                if(x*x>=9.)continue;
                const double height=(u.hoop+v.hoop)*std::sqrt(3.)*.25;
                p.bonds.push_back({i,k,x*x,((k%lanes)-(i%lanes))*height,lanes*height,
                    std::sqrt(9.-x*x),.5*(u.order+v.order)/entry});
                maxdiff=std::max(maxdiff,std::abs(u.hoop-v.hoop));
            }
            phases.push_back(std::move(p));
        }
        std::cout<<"{\"test\":\"strip\",\"rows\":"<<rows<<",\"lanes\":"<<lanes<<",\"sub\":"<<sub
            <<",\"preparation\":"<<preparation<<",\"step_bond_lengths\":"<<h
            <<",\"outlet_arc\":"<<path[rows*table_sub].arc<<",\"largest_bond_hoop_difference\":"<<maxdiff<<"}\n";
    }
    int count() const {return rows*lanes;}
    double anchor(int i,const Phase& p) const {return i%2 ? 0.:2*p.site[i].contact/sites;}
    double inlet(int arrival,int k) const {
        // Opposite small displacements over each pair of arrivals give zero mean
        // feed tilt. Preparation 2 supplies exactly centered incoming sections.
        if(preparation==2)return 0;
        return 1e-5*((arrival+(preparation==1?k:0))%2==0?1.:-1.);
    }
    void inject(Vec& u,int arrival) const {
        std::move_backward(u.begin(),u.end()-lanes,u.end());
        for(int k=0;k<lanes;++k)u[k]=inlet(arrival,k);
    }
    void velocity(const Vec& u,double amplitude,const Phase& p,const Law& law,Vec& du) const {
        for(int i=0;i<count();++i){
            double slip=u[i]-amplitude*p.site[i].h;
            if(i%2==0 && (!std::isfinite(slip) || std::abs(slip)>=law.spacing*(law.energy.size()-1))) {
                std::cerr<<"contact outside table: arc="<<p.site[i].arc<<" row="<<i/lanes<<" lane="<<i%lanes<<" slip="<<slip<<'\n';
                throw std::runtime_error("Unresolved contact");
            }
            du[i]=i%2?0.:-anchor(i,p)*law.wall(slip)[1];
        }
        // Sum every transverse periodic image in range. The hoop width at the
        // midpoint sets the transverse metric for each short symmetric pair.
        for(const auto& b:p.bonds) {
            const double raw=b.base_y+u[b.j]-u[b.i];
            double y=raw-b.period*std::floor((raw+b.ymax)/b.period);
            double f=0;
            for(;y<=b.ymax;y+=b.period) f+=law.pair(b.dx2+y*y)*y;
            f*=b.weight;du[b.i]+=f;du[b.j]-=f;
        }
    }
    void advance(Vec& u,double amp,int j,const Law& law,Vec& du,Vec& mid) const {
        // Freeze the advected metric at the macrostep midpoint. Integrate its
        // gradient flow with an embedded Euler/midpoint error estimate. Shorter
        // internal steps resolve the stiff repulsion; no physical force is capped.
        const auto& geom=phases[2*j+1];
        const double interval=period/sub;
        double elapsed=0, trial=interval;
        Vec half_slope(count());
        constexpr double tolerance=2e-6, max_travel=.025;
        while(elapsed<interval) {
            velocity(u,amp,geom,law,du);
            double max_speed=0;
            for(double v:du) max_speed=std::max(max_speed,std::abs(v));
            double dt=std::min({trial, interval-elapsed, max_travel/std::max(max_speed,1e-16)});
            while(true) {
                if(elapsed+dt==elapsed)throw std::runtime_error("Material timestep lost precision");
                for(int i=0;i<count();++i)mid[i]=u[i]+dt*.5*du[i];
                velocity(mid,amp,geom,law,half_slope);
                double err=0, travel=0;
                for(int i=0;i<count();++i) {
                    err=std::max(err,std::abs(half_slope[i]-du[i])*dt);
                    travel=std::max(travel,std::abs(half_slope[i])*dt);
                }
                if(err<=tolerance && travel<=max_travel) {
                    for(int i=0;i<count();++i)u[i]+=dt*half_slope[i];
                    elapsed+=dt;
                    trial=dt*std::clamp(.9*std::sqrt(tolerance/std::max(err,1e-30)),.5,2.);
                    break;
                }
                dt*=std::clamp(.85*std::sqrt(tolerance/std::max(err,1e-30)),.1,.5);
            }
        }
    }
    double force(const Vec& u,const Vec& ref,double amp,const Phase& p,const Law& law,double& work) const {
        double sum=0, fd=0;
        constexpr double dh=1e-5;
        for(int i=0;i<count();i+=2) {
            const auto& s=p.site[i];
            if(s.arc<low || s.arc>=high)continue;
            double a=u[i]-amp*s.h,b=ref[i],w=anchor(i,p)*sites/lanes;
            sum+=w*(law.wall(a)[1]-law.wall(b)[1]);
            fd+=w*((law.wall(a+dh)[0]-law.wall(b+dh)[0])-(law.wall(a-dh)[0]-law.wall(b-dh)[0]))/(2*dh);
        }
        work=std::max(work,std::abs(fd-sum));return sum;
    }
    void run(const Law& law) const {
        Vec ref(count()),mid(count()),du(count());
        for(int i=0;i<count();++i)ref[i]=1e-5*std::cos(1.84*i+.1);
        // Actual material injection and exit are present during both stages.
        for(int n=0;n<2*rows;++n) {
            for(int j=0;j<sub;++j){
                try {advance(ref,0,j,law,du,mid);}
                catch(...) {std::cerr<<"Unforced initialization: arrival="<<n<<" step="<<j<<'\n';throw;}
            }
            inject(ref,n+1);
        }
        Vec plus=ref,minus=ref;
        const std::array<double,3> gate={middle-16*(high-low)/12.,low,middle};
        std::array<double,3> tags{},tag2{};std::array<int,3> ntag{};
        int samples=0;double sp=0,sm=0,work=0,alignment=0;
        for(int n=0;n<2*rows;++n) {
            bool measure=n>=2*rows-16;
            for(int j=0;j<sub;++j) {
                const auto& p=phases[2*j];
                if(measure && j%(sub/8)==0) {
                    sp+=force(plus,ref,h,p,law,work);
                    sm+=force(minus,ref,-h,p,law,work);++samples;
                }
                Vec before;
                if(measure) {
                    before.resize(count());
                    for(int i=0;i<count();++i)before[i]=(plus[i]-minus[i])/(2*h);
                }
                advance(ref,0,j,law,du,mid); advance(plus,h,j,law,du,mid);advance(minus,-h,j,law,du,mid);
                if(measure) {
                    const auto& next=phases[2*j+2];
                    for(int i=0;i<count();i+=2) for(size_t k=0;k<gate.size();++k) {
                        double a=p.site[i].arc,b=next.site[i].arc;
                        if(a<=gate[k] && b>gate[k]) {
                            if(i/lanes>=n)throw std::runtime_error("Read a pre-drive tag");
                            double t=(gate[k]-a)/(b-a);
                            double v=(1-t)*before[i]+t*(plus[i]-minus[i])/(2*h);
                            tags[k]+=v;tag2[k]+=v*v;++ntag[k];
                            alignment=std::max(alignment,std::abs(plus[i]+minus[i]-2*ref[i])/(2*h));
                        }
                    }
                }
            }
            inject(ref,n+1);inject(plus,n+1);inject(minus,n+1);
            if(measure && (n+1)%8==0) {
                std::cout<<"{\"test\":\"fresh_window\",\"last_arrival\":"<<n+1<<",\"positive_wall\":"<<sp/samples
                    <<",\"negative_wall\":"<<sm/samples<<",\"contacted_tag_mean_per_step\":[";
                for(int k=0;k<3;++k) {
                    if(ntag[k]!=8*lanes/2)throw std::runtime_error("Incomplete crossing cohort");
                    std::cout<<(k?",":"")<<tags[k]/ntag[k];
                }
                std::cout<<"],\"tag_rms_per_step\":[";
                for(int k=0;k<3;++k)std::cout<<(k?",":"")<<std::sqrt(tag2[k]/ntag[k]);
                std::cout<<"],\"cohort_per_gate\":"<<ntag[0]<<",\"energy_force_error\":"<<work
                    <<",\"max_even_tag_fraction\":"<<alignment<<"}\n";
                if(!std::isfinite(sp+sm+tags[0]+tags[1]+tags[2])||work>1e-7)
                    throw std::runtime_error("Energy derivative or material motion check failed");
                samples=0;sp=sm=0;work=alignment=0;tags={};tag2={};ntag={};
            }
        }
    }
};
int main(int argc,char** argv) {
    try {
        if(argc!=8)throw std::runtime_error("usage: transport path wall rows lanes sub preparation amplitude");
        std::cout<<std::setprecision(16)<<std::unitbuf;
        Stream s(argv); Law law(argv[2],s.depth,s.alpha,s.rest); s.run(law);return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
