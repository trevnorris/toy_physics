// Cold, open, laterally mobile patch. Axial compression, throughflow and hoop size
// are controlled independently. Lengths are in bare material bond lengths.
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using Copies = std::vector<Vec>;
constexpr double pi = 3.14159265358979323846;
double dot(const Vec& a,const Vec& b) { return std::inner_product(a.begin(),a.end(),b.begin(),0.); }
double rms(const Vec& a) { return std::sqrt(dot(a,a)/a.size()); }
void array(const Vec& v) {
    std::cout << "[";
    for (size_t i=0;i<v.size();++i) std::cout << (i?",":"") << v[i];
    std::cout << "]";
}
struct Pair {int a,b; double dx2,y0;};
struct Profile {Vec p, pd; double norm;};
struct Forces {Vec internal,wall; double energy;};
struct Reading {double native,wall,transport,total;};

struct Patch {
    int length,lanes,sub,seed;
    bool mobile_graph;
    double gap,hoop,period,order,section_wall,h,rest,alpha,depth,entry,sites,cutoff;
    double width,wavelength,seed_sign;
    std::vector<Pair> pairs;
    Vec wall_energy,wall_slope;
    double wall_spacing=0.;
    Patch(int argc,char** a):
        length(std::stoi(a[1])),lanes(std::stoi(a[2])),sub(std::stoi(a[3])),seed(std::stoi(a[4])),
        mobile_graph(std::string(a[5])=="all"),gap(std::stod(a[6])),hoop(std::stod(a[7])),
        period(std::stod(a[8])),order(std::stod(a[9])),section_wall(std::stod(a[10])),
        h(std::stod(a[11])),rest(std::stod(a[12])),alpha(std::stod(a[13])),depth(std::stod(a[14])),
        entry(std::stod(a[15])),sites(std::stod(a[16])),cutoff(std::stod(a[17])),
        width(std::stod(a[18])),wavelength(std::stod(a[19])),seed_sign(std::stod(a[20])) {
        if(lanes%2 || length<2*width+4 || cutoff>=circumference()/2 || rest!=1.)
            throw std::runtime_error("Invalid periodic strip geometry");
        if(argc==22) {
            std::ifstream in(a[21]);
            if(!(in>>wall_spacing) || wall_spacing<=0) throw std::runtime_error("Invalid wall table");
            double energy,slope;
            while(in>>energy>>slope) {wall_energy.push_back(energy);wall_slope.push_back(slope);}
            if(wall_energy.size()<2) throw std::runtime_error("Incomplete wall table");
        }
        for(int i=0;i<count();++i) for(int j=i+1;j<count();++j) {
            double x=reference_x(j)-reference_x(i);
            double y=wrap((j%lanes-i%lanes)*row_height());
            bool nearest=std::abs(x*x+y*y/(hoop*hoop)-1.)<1e-8;
            if(mobile_graph ? std::abs(gap*x)<cutoff : nearest)
                pairs.push_back({i,j,gap*gap*x*x,y});
        }
    }
    int count() const {return length*lanes;}
    double row_height() const {return std::sqrt(3.)*.5*hoop;}
    double circumference() const {return row_height()*lanes;}
    double wrap(double y) const {return y-circumference()*std::nearbyint(y/circumference());}
    double reference_x(int i) const {return i/lanes+.5*((i%lanes)%2);}
    double weight() const {return sites/lanes;}
    double anchor(int i) const {return (i%lanes)%2==0 ? section_wall*2/sites : 0.;}
    Profile profile(double phase) const {
        Profile q{Vec(count()),Vec(count()),0.};
        for(int i=0;i<count();++i) {
            double x=reference_x(i)+phase-length*.5;
            if(std::abs(x)>=width) continue;
            double c=.5*(1+std::cos(pi*x/width)), wave=2*pi*x/wavelength;
            q.p[i]=c*std::cos(wave);
            q.pd[i]=(-.5*pi/width*std::sin(pi*x/width)*std::cos(wave)
                     -2*pi/wavelength*c*std::sin(wave))/period;
        }
        q.norm=dot(q.p,q.p);
        return q;
    }
    // Cubic Hermite interpolation keeps the gradient conjugate to the energy
    // even at finite relaxed offsets. The two arrays are normalized to unit curvature at zero.
    std::array<double,2> wall_law(double u) const {
        if(wall_energy.empty()) return {.5*u*u,u};
        double r=std::abs(u)/wall_spacing;
        size_t j=static_cast<size_t>(r);
        if(j+1>=wall_energy.size()) throw std::runtime_error("Outside the calibrated wall potential");
        double t=r-j, change=wall_energy[j+1]-wall_energy[j];
        double b=3*change-wall_spacing*(2*wall_slope[j]+wall_slope[j+1]);
        double c=wall_spacing*(wall_slope[j]+wall_slope[j+1])-2*change;
        return {wall_energy[j]+t*(wall_spacing*wall_slope[j]+t*(b+t*c)),
                (u>=0?1.:-1.)*(wall_slope[j]+t*(2*b+3*t*c)/wall_spacing)};
    }
    Forces forces(const Vec& u, bool need_energy=false) const {
        Forces out{Vec(count()),Vec(count()),0.};
        for(const auto& e:pairs) {
            double y=e.y0+u[e.b]-u[e.a];
            // Persistent bonds retain their original hoop image; mobile pairs can exchange neighbors.
            if(mobile_graph) y=wrap(y);
            double r2=e.dx2+y*y;
            if(mobile_graph && r2>=cutoff*cutoff) continue;
            if(r2<1e-16) throw std::runtime_error("Coincident sites");
            double r=std::sqrt(r2), v=std::exp(-alpha*(r-rest));
            double energy=depth*v*(v-2), tension=2*depth*alpha*v*(1-v);
            if(mobile_graph && r>cutoff-.5) {
                double t=(r-cutoff+.5)*2;
                double sw=1-t*t*t*(10+t*(-15+6*t));
                tension=tension*sw-energy*60*t*t*(1-t)*(1-t);
                energy*=sw;
            }
            double f=(order/entry)*tension*y/r;
            out.internal[e.a]+=f; out.internal[e.b]-=f;
            if(need_energy) out.energy+=weight()*(order/entry)*energy;
        }
        for(int i=0;i<count();++i) {
            if(anchor(i)==0) continue;
            auto law=wall_law(u[i]);
            out.wall[i]=-anchor(i)*law[1];
            if(need_energy) out.energy+=weight()*anchor(i)*law[0];
        }
        return out;
    }
    // Equality-constrained perturbation relative to the evolving unforced copy.
    // Each mobile site has unit drag. Only the prescribed displacement mode is held.
    Copies slopes(const Copies& us,const Vec& amps,const Profile& q,std::vector<Reading>* reads=nullptr) const {
        Copies v;
        Forces base;
        for(size_t c=0;c<us.size();++c) {
            Forces f=forces(us[c]);
            if(c==0) base=f;
            double ni=0,wa=0,transport=0;
            Vec out(count());
            for(int i=0;i<count();++i) {
                out[i]=f.internal[i]+f.wall[i];
                if(c) {
                    ni+=q.p[i]*(base.internal[i]-f.internal[i]);
                    wa+=q.p[i]*(base.wall[i]-f.wall[i]);
                    transport+=q.pd[i]*(2*amps[c]*q.p[i]-us[c][i]+us[0][i]);
                }
            }
            double load=ni+wa+transport;
            if(c) for(int i=0;i<count();++i) out[i]+=load/q.norm*q.p[i];
            if(reads) reads->push_back({weight()*ni,weight()*wa,weight()*transport,weight()*load});
            v.push_back(std::move(out));
        }
        return v;
    }
    void project(Copies& us,const Vec& amps,const Profile& q) const {
        for(size_t c=1;c<us.size();++c) {
            double err=amps[c]-(dot(q.p,us[c])-dot(q.p,us[0]))/q.norm;
            for(int i=0;i<count();++i) us[c][i]+=err*q.p[i];
        }
    }
    void inject(Copies& us) const {
        for(auto& u:us) {
            std::move_backward(u.begin(),u.end()-lanes,u.end());
            std::fill(u.begin(),u.begin()+lanes,0.);
        }
    }
    Vec center(const Vec& u) const {
        return Vec(u.begin()+(length/2-6)*lanes,u.begin()+(length/2+6)*lanes);
    }
    void run() const {
        const double dt=period/sub;
        std::vector<Profile> phases;
        for(int i=0;i<=2*sub;++i) phases.push_back(profile(i/(2.*sub)));
        Vec initial(count());
        std::mt19937 gen(std::abs(seed)+1234);
        std::uniform_real_distribution<double> rand(-1.,1.);
        for(int i=0;i<count();++i) initial[i]=seed_sign*1e-5*rand(gen);
        Copies u{initial}; Vec amps{0.};
        auto advance=[&](int j) {
            Copies m=u;
            auto k=slopes(u,amps,phases[2*j]);
            for(size_t c=0;c<u.size();++c) for(int i=0;i<count();++i) m[c][i]+=.5*dt*k[c][i];
            project(m,amps,phases[2*j+1]);
            k=slopes(m,amps,phases[2*j+1]);
            for(size_t c=0;c<u.size();++c) for(int i=0;i<count();++i) u[c][i]+=dt*k[c][i];
            project(u,amps,phases[2*j+2]);
        };
        std::vector<Vec> history;
        auto once=[&]() {
            for(int j=0;j<sub;++j) advance(j);
            inject(u);
            for(auto &v:u) for(double x:v) if(!std::isfinite(x)) throw std::runtime_error("Invalid offset");
        };
        for(int n=1;n<=6*length;++n) {
            once(); history.push_back(center(u[0]));
            if(n%length==0) {
                std::cout<<"{\"test\":\"passage\",\"sections_replaced\":"<<n<<",\"rms_offset\":"<<rms(history.back());
                std::cout<<",\"phase_repeat_rms\":[";
                for(int lag:{1,2,4,8}) {
                    double largest=0;
                    for(int off=0;off<8;++off) {
                        Vec diff=history[n-1-off];
                        for(size_t i=0;i<diff.size();++i) diff[i]-=history[n-1-off-lag][i];
                        largest=std::max(largest,rms(diff));
                    }
                    std::cout<<(lag==1?"":",")<<largest;
                }
                std::cout<<"]}\n";
            }
        }
        auto f0=forces(u[0],true);
        std::cout<<"{\"test\":\"released_state\",\"center\":";array(center(u[0]));
        std::cout<<",\"energy\":"<<f0.energy<<",\"elapsed\":"<<6*length*period<<"}\n";
        amps={0,h,-h,h/2,-h/2};
        Vec baseline=u[0];u.assign(amps.size(),baseline);
        const auto& q=phases[0];
        for(size_t c=1;c<u.size();++c) for(int i=0;i<count();++i) u[c][i]+=amps[c]*q.p[i];
        std::vector<Reading> first;
        slopes(u,amps,q,&first);
        auto odd=[&](const std::vector<Reading>& r,int c, double size) {
            return std::array<double,4>{(r[c].native-r[c+1].native)/(2*size),
                (r[c].wall-r[c+1].wall)/(2*size),(r[c].transport-r[c+1].transport)/(2*size),
                (r[c].total-r[c+1].total)/(2*size)};
        };
        auto inst=odd(first,1,h);
        double d2=(forces(u[1],true).energy+forces(u[2],true).energy-2*f0.energy)/(h*h);
        std::cout<<"{\"test\":\"instant\",\"energy_curvature\":"<<d2<<",\"native\":"<<inst[0]
                 <<",\"wall\":"<<inst[1]<<",\"transport\":"<<inst[2]<<",\"total\":"<<inst[3]<<"}\n";
        if(std::abs(d2-inst[0]-inst[1])>2e-4*std::max(1.,std::abs(d2)))
            throw std::runtime_error("Direct energy and force differ");
        struct Window {
            std::array<double,4> sum{};
            int samples=0;
            double half_err=0, min=1e300,max=-1e300, peak_offset=0, constraint=0, plus=0,minus=0;
        } window;
        for(int n=0;n<2*length;++n) {
            for(int j=0;j<sub;++j) {
                if(j%(sub/8)==0 && (n%length)>=length-8) {
                    const Profile& p=phases[2*j];std::vector<Reading> out;
                    slopes(u,amps,p,&out);
                    auto main=odd(out,1,h),half=odd(out,3,h/2);
                    for(int k=0;k<4;++k) window.sum[k]+=main[k];
                    window.min=std::min(window.min,main[3]);window.max=std::max(window.max,main[3]);
                    window.half_err=std::max(window.half_err,std::abs(main[3]-half[3]));
                    window.plus+=out[1].total;window.minus+=out[2].total;
                    for(size_t c=1;c<u.size();++c) {
                        window.constraint=std::max(window.constraint,std::abs((dot(p.p,u[c])-dot(p.p,u[0]))/p.norm-amps[c]));
                        for(int i=0;i<count();++i) window.peak_offset=std::max(window.peak_offset,std::abs(u[c][i]-u[0][i]));
                    }
                    ++window.samples;
                }
                advance(j);
            }
            inject(u);
            if((n+1)%length==0) {
                double norm=window.samples;
                std::cout<<"{\"test\":\"held\",\"sections_since_step\":"<<n+1<<",\"samples\":"<<norm
                    <<",\"native\":"<<window.sum[0]/norm<<",\"wall\":"<<window.sum[1]/norm
                    <<",\"transport\":"<<window.sum[2]/norm<<",\"total\":"<<window.sum[3]/norm
                    <<",\"min\":"<<window.min<<",\"max\":"<<window.max<<",\"half_error_max\":"<<window.half_err
                    <<",\"positive_force\":"<<window.plus/norm<<",\"negative_force\":"<<window.minus/norm
                    <<",\"constraint_error\":"<<window.constraint<<",\"max_copy_offset\":"<<window.peak_offset<<"}\n";
                if(!std::isfinite(window.sum[3]) || window.constraint>1e-10) throw std::runtime_error("Unresolved hold");
                window=Window{};
            }
        }
    }
};
int main(int argc,char** argv) {
    try {
        if(argc!=21 && argc!=22) throw std::runtime_error("Pass local patch preparation via Python driver");
        std::cout<<std::setprecision(16);
        Patch(argc,argv).run();
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n"; return 1;}
}
