// Homogeneous four-dimensional material under a fixed-shape, fluctuating-volume hold.
// Row-vector lattice, unwrapped coordinates, unit constituent mass. Potential and full
// periodic image enumeration are exactly the same implementations as the archive.
#include "mouth_bulk_thermal_cell.cpp"
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

namespace {
constexpr double params[] = {1., 7., 1.};
constexpr int nblocks = 16;
using Sites = std::vector<double>;
static_assert(sizeof(Matrix)==16*sizeof(double));
Matrix identity() {
    Matrix m{};
    for(int k=0;k<4;++k) m[k][k]=1.;
    return m;
}
Matrix product(const Matrix& a,const Matrix& b) {
    Matrix c{};
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
        for(int k=0;k<4;++k) c[i][j]+=a[i][k]*b[k][j];
    return c;
}
Matrix deformation(const Matrix& axis,double eps) {
    Matrix s=identity(),t=s;
    for(int l=1;l<=20;++l) {
        t=product(t,axis);
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) {t[i][j]*=eps/l;s[i][j]+=t[i][j];}
    }
    return s;
}
void transform(Sites& x,const Matrix& a) {
    for(size_t i=0;i<x.size();i+=4) {
        Vec tmp{};
        for(int j=0;j<4;++j) for(int k=0;k<4;++k) tmp[j]+=x[i+k]*a[k][j];
        std::copy(tmp.begin(),tmp.end(),x.begin()+i);
    }
}
double determinant(Matrix a) {
    double result=1.;
    for(int i=0;i<4;++i) {
        int p=i;
        for(int j=i+1;j<4;++j) if(std::abs(a[j][i])>std::abs(a[p][i])) p=j;
        if(p!=i) {std::swap(a[p],a[i]);result=-result;}
        result*=a[i][i];
        if(a[i][i]==0.) return 0.;
        for(int j=i+1;j<4;++j) {
            double f=a[j][i]/a[i][i];
            for(int k=i+1;k<4;++k) a[j][k]-=f*a[i][k];
        }
    }
    return result;
}
template<class T> void vec(std::ostream& os,const T& values) {
    os<<'[';bool first=true;
    for(const auto& v:values) {if(!first) os<<',';os<<v;first=false;}
    os<<']';
}
void matrix(std::ostream& os,const Matrix& a) {
    os<<'[';for(int i=0;i<4;++i) {if(i) os<<',';vec(os,a[i]);} os<<']';
}
void sites(std::ostream& os,const Sites& a,double scale) {
    os<<'[';
    for(size_t i=0;i<a.size();i+=4) {
        if(i) os<<',';
        Vec v;
        for(int k=0;k<4;++k) v[k]=scale*a[i+k];
        vec(os,v);
    }
    os<<']';
}
struct Input {
    int n{},history{};
    double temperature{},mobility{},rc{};
    std::string sha;
    Matrix box{},inv{},axis{};
    Sites x;
    explicit Input(const std::string& path) {
        std::ifstream is(path);
        if(!(is>>n>>history>>temperature>>mobility>>rc>>sha) || n!=128 || !(temperature>=0.) || !std::isfinite(temperature) || rc!=2.5 || !(mobility>0.) || !std::isfinite(mobility) || sha.size()!=64)
            throw std::runtime_error("invalid archived preparation header");
        for(auto* m:{&box,&inv,&axis}) for(auto& row:*m) for(auto& a:row)
            if(!(is>>a) || !std::isfinite(a)) throw std::runtime_error("invalid input matrix");
        x.resize(4*n);
        for(auto& a:x) if(!(is>>a)||!std::isfinite(a)) throw std::runtime_error("invalid input site");
        std::string extra;
        if(is>>extra) throw std::runtime_error("extra input records");
        const Matrix unit=product(box,inv);
        double trace=0.,norm=0.;
        for(int i=0;i<4;++i) for(int j=0;j<4;++j) {
            if(std::abs(unit[i][j]-(i==j))>1.e-10||std::abs(axis[i][j]-axis[j][i])>1.e-10) throw std::runtime_error("invalid inverse or axis");
            if(i==j) trace+=axis[i][j];
            norm+=axis[i][j]*axis[i][j];
        }
        if(std::abs(trace)>1.e-10||std::abs(norm-.5)>1.e-10||determinant(box)<=0.) throw std::runtime_error("invalid archived geometry");
    }
};
struct Forces {
    Sites grad;
    Matrix virial{};
    double energy=0.;
    explicit Forces(int n):grad(4*n) {}
};
// Fixed reference shape; physical x and box are the reference variables times scale.
// Retaining this list across log-volume changes is safe while it covers cutoff/scale.
struct Interactions {
    Cell cell;
    double reach=0., cutoff;
    std::vector<Vec> self;
    Interactions(int n,const Matrix& box,const Matrix& inv,double rc):
        cell(n,&box[0][0],&inv[0][0],rc,params),cutoff(rc) {}
    void ensure(const Sites& x,double scale) {
        const double needed=cutoff/scale;
        if(needed>reach || needed<.85*reach) {
            reach=1.03*needed;
            cell.reference.clear();
            self.clear();
            cell.each_image(Vec{},reach,[&](const Tile& t,const Vec& v,double){if(t>Tile{}) self.push_back(v);});
        }
        // Only the image-search radius is changed. The actual pair evaluation
        // below uses the physical cutoff, including its original analytic switch.
        cell.potential.cutoff=reach;
        cell.rebuild_if_needed(x.data());
        cell.potential.cutoff=cutoff;
    }
    template<bool details> void read(const Sites& x,double scale,Forces& f) {
        ensure(x,scale);
        f.energy=0.;f.virial={};std::fill(f.grad.begin(),f.grad.end(),0.);
        const int n=cell.n;
        const double max_r2=cutoff*cutoff, ss=scale*scale, ref_cut2=max_r2/ss;
        for(int i=0;i<n;++i) {
            Vec gi{};
            for(int h=cell.first[i];h<cell.first[i+1];++h) {
                const auto& gp=cell.groups[h];
                Vec org;
                for(int k=0;k<4;++k) org[k]=gp.shift[k]-x[4*i+k];
                for(int b=gp.begin;b<gp.end;++b) {
                    const int j=cell.neighbors[b];Vec u;double rr=0.;
                    for(int k=0;k<4;++k) {u[k]=x[4*j+k]+org[k];rr+=u[k]*u[k];}
                    if(rr>=ref_cut2 || rr*ss>=max_r2) continue;
                    if(!(rr>0.)) throw std::runtime_error("coincident constituent sites");
                    double e=0.;
                    const double g=cell.potential.force_prefactor<details>(rr*ss,e)*scale;
                    if constexpr(details) f.energy+=e;
                    for(int k=0;k<4;++k) {
                        const double du=g*u[k];gi[k]-=du;f.grad[4*j+k]+=du;
                        if constexpr(details) for(int l=0;l<4;++l) f.virial[k][l]+=du*u[l]*scale;
                    }
                }
            }
            for(int k=0;k<4;++k) f.grad[4*i+k]+=gi[k];
        }
        if constexpr(details) for(const auto& u:self) {
            double rr=std::inner_product(u.begin(),u.end(),u.begin(),0.)*ss;
            if(rr>=max_r2) continue;
            double e;
            const double g=cell.n*ss*cell.potential.force_prefactor<true>(rr,e);
            f.energy+=cell.n*e;
            for(int k=0;k<4;++k) for(int l=0;l<4;++l) f.virial[k][l]+=g*u[k]*u[l];
        }
    }
    std::vector<uint8_t> neighbors(const Sites& x,double scale) {
        ensure(x,scale);
        const int n=cell.n;
        std::vector<uint8_t> out(n*n,0);
        for(int i=0;i<n;++i) for(int h=cell.first[i];h<cell.first[i+1];++h) {
            const auto& gp=cell.groups[h];
            for(int b=gp.begin;b<gp.end;++b) {
                const int j=cell.neighbors[b];double rr=0.;
                for(int k=0;k<4;++k) {double d=scale*(x[4*j+k]+gp.shift[k]-x[4*i+k]);rr+=d*d;}
                const uint8_t level=rr<1.15*1.15 ? 2 : rr<1.30*1.30 ? 1 : 0;
                out[n*i+j]=std::max(out[n*i+j],level);
            }
        }
        return out;
    }
};
struct Tracking {
    Sites ref;
    double reference_scale;
    std::vector<uint8_t> original,current;
    int norig=0;
    Tracking(const Sites& x,double scale,Interactions& p):ref(x),reference_scale(scale),original(p.neighbors(x,scale)) {
        for(auto& flag:original) {flag=(flag==2);norig+=flag;}
        current=original;
    }
    void write(std::ostream& os,const Sites& x,double scale,Interactions& p) {
        int n=p.cell.n,lost=0,gain=0,largest=1;
        auto edges=p.neighbors(x,scale);
        std::vector<int> parent(n),size(n,1);std::iota(parent.begin(),parent.end(),0);
        auto root=[&](int i){while(i!=parent[i]) i=parent[i];return i;};
        Vec shift{};double square=0.;
        for(size_t i=0;i<x.size();++i) shift[i%4]+=(x[i]-ref[i])/n;
        for(size_t i=0;i<x.size();++i) {const double dx=reference_scale*(x[i]-ref[i]-shift[i%4]);square+=dx*dx;}
        for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) {
            int id=n*i+j;
            current[id]=(current[id]&&edges[id]>=1)||edges[id]==2;
            lost+=original[id]&&!current[id];gain+=!original[id]&&current[id];
            if(current[id]) {int a=root(i),b=root(j);if(a!=b){parent[b]=a;size[a]+=size[b];largest=std::max(largest,size[a]);}}
        }
        os<<"{\"original\":"<<norig<<",\"absent\":"<<lost<<",\"new\":"<<gain<<",\"max_cluster\":"<<largest
          <<",\"nonaffine_rms\":"<<std::sqrt(square/n)<<'}';
    }
};
struct Reading {
    Matrix stress{}; // tensile, including kinetic momentum flux
    double pressure=0.,shear=0.,kinetic_T=0.,u_site=0.,rho=0.,volume=0.,logv=0.;
    void add(const Reading& a) {
        for(int k=0;k<4;++k) for(int l=0;l<4;++l) stress[k][l]+=a.stress[k][l];
        pressure+=a.pressure;shear+=a.shear;kinetic_T+=a.kinetic_T;u_site+=a.u_site;rho+=a.rho;volume+=a.volume;logv+=a.logv;
    }
    void write(std::ostream& os,int samples=1) const {
        Matrix w=stress;for(auto& row:w) for(double& value:row) value/=samples;
        os<<"{\"shear\":"<<shear/samples<<",\"pressure\":"<<pressure/samples<<",\"temperature\":"<<kinetic_T/samples
          <<",\"potential_per_site\":"<<u_site/samples<<",\"density\":"<<rho/samples<<",\"volume\":"<<volume/samples
          <<",\"log_volume\":"<<logv/samples<<",\"stress_tensile\":";matrix(os,w);os<<'}';
    }
};
Reading observe(const Forces& f,const Sites& v,const Matrix& axis,double volume) {
    Reading r;r.volume=volume;r.logv=std::log(volume);int n=int(v.size())/4;r.rho=n/volume;r.u_site=f.energy/n;
    for(int k=0;k<4;++k) for(int l=0;l<4;++l) {
        double kinetic=0.;for(int i=0;i<n;++i) kinetic+=v[4*i+k]*v[4*i+l];
        r.stress[k][l]=(f.virial[k][l]-kinetic)/volume;
        r.shear+=r.stress[k][l]*axis[k][l];
        if(k==l) {r.pressure-=r.stress[k][k]/4;r.kinetic_T+=kinetic/(4*n);}
    }
    return r;
}
struct Evolution {
    Sites x,v;
    Matrix box,inv,axis;
    double scale=1.,base_volume,temperature,mobility,cutoff;
    int count;
    Forces f,trial;
    std::mt19937_64 thermal_rng,volume_rng;
    std::normal_distribution<double> normal_site,normal_volume;
    std::uniform_real_distribution<double> uniform;
    int attempts=0,accepts=0;
    Evolution(const Input& a,uint64_t seed,int direction,double trial_T):x(a.x),v(a.x.size()),box(a.box),inv(a.inv),axis(a.axis),
        base_volume(determinant(a.box)),temperature(trial_T),mobility(a.mobility),cutoff(a.rc),count(a.n),
        f(a.n),trial(a.n),thermal_rng(seed),volume_rng(seed^0x2d693923439c104ULL) {
        if(direction==1) {axis={};axis[0][3]=axis[3][0]=.5;}
        for(auto& speed:v) speed=std::sqrt(temperature)*normal_site(thermal_rng);
    }
    double volume() const {const double ss=scale*scale;return base_volume*ss*ss;}
    void strain(double eps) {
        auto a=deformation(axis,eps), b=deformation(axis,-eps);
        transform(x,a);box=product(box,a);inv=product(b,inv);
        if(std::abs(determinant(box)/base_volume-1.)>1.e-11) throw std::runtime_error("volume-changing shear");
    }
    void volume_move(Interactions& p,double h,double multiplier,double imposed) {
        // At fixed fractional positions and physical velocities the log-volume measure
        // is exp((N+1)q - (U+Pext*exp(q))/T). Momenta contribute no extra factor.
        // Evaluate even rejected branches with independent volume RNG consumption.
        const double dq=std::sqrt(2*mobility*multiplier*temperature*h)*normal_volume(volume_rng);
        const double u=uniform(volume_rng), proposed=scale*std::exp(.25*dq);
        if(!(proposed>0.) || !std::isfinite(proposed)) throw std::runtime_error("nonfinite trial scale");
        p.read<true>(x,proposed,trial);
        const double log_ratio=(count+1)*dq - (trial.energy-f.energy+imposed*volume()*std::expm1(dq))/temperature;
        if(!std::isfinite(log_ratio)) throw std::runtime_error("nonfinite pressure move");
        ++attempts;
        if(std::log(u)<log_ratio) {scale=proposed;std::swap(f,trial);++accepts;}
    }
    void advance(Interactions& p,double dt,bool details,bool volume_kick,double baro_dt,double multiplier,double imposed) {
        const double h=.5*dt, drift=h/scale, c=std::exp(-dt),noise=std::sqrt(-std::expm1(-2*dt)*temperature);
        for(size_t k=0;k<x.size();++k) {
            v[k]-=h*f.grad[k];x[k]+=drift*v[k];
            v[k]=c*v[k]+noise*normal_site(thermal_rng);x[k]+=drift*v[k];
            if(!std::isfinite(x[k])||!std::isfinite(v[k])) throw std::runtime_error("nonfinite receiver state");
        }
        if(details||volume_kick) p.read<true>(x,scale,f);else p.read<false>(x,scale,f);
        for(size_t k=0;k<v.size();++k) v[k]-=h*f.grad[k];
        if(volume_kick) volume_move(p,baro_dt,multiplier,imposed);
    }
    void write_state(std::ostream& os) const {
        Matrix b=box;for(auto& row:b)for(auto& a:row)a*=scale;
        os<<"{\"box\":";matrix(os,b);os<<",\"pos\":";sites(os,x,scale);os<<",\"vel\":";sites(os,v,1.);os<<'}';
    }
};
struct Mesh {
    int steps,baro_steps;
    double dt,baro_dt;
    Mesh(double duration,double target_dt,double baro_interval):steps(int(std::ceil(duration/target_dt))),
        dt(steps?duration/steps:target_dt) {
        baro_steps=std::max(1,int(std::lround(baro_interval/dt)));
        baro_dt=dt*baro_steps;
    }
};
int execute(int argc,char** argv) {
    if(argc<8||argc>14) throw std::invalid_argument("INPUT OUT SEED DT WARM MEASURE STRAIN [BARO_INTERVAL=.08] [MOBILITY_MULT=1] [PEXT=0] [AXIS=0] [TARGET_T=source] [SAVE_BLOCK_STATES=0]");
    Input input(argv[1]);const auto seed=std::stoull(argv[3]);
    const double dt=std::stod(argv[4]),warm=std::stod(argv[5]),measure=std::stod(argv[6]),eps=std::stod(argv[7]);
    const double baro_interval=argc>8?std::stod(argv[8]):.08,mult=argc>9?std::stod(argv[9]):1.,pext=argc>10?std::stod(argv[10]):0.;
    const int axis=argc>11?std::stoi(argv[11]):0;
    const double trial_T=argc>12?std::stod(argv[12]):input.temperature;
    const int snapshots=argc>13?std::stoi(argv[13]):0;
    if(snapshots!=0 && snapshots!=1) throw std::invalid_argument("invalid block-state output flag");
    if(!(trial_T>0.) || !std::isfinite(trial_T)) throw std::invalid_argument("invalid target temperature");
    if(!(dt>0&&warm>=0&&measure>0&&warm/dt<2.e9&&measure/dt<2.e9&&std::abs(eps)<=.3&&baro_interval>=dt&&baro_interval/dt<2.e9&&mult>0&&std::isfinite(mult)&&std::isfinite(pext)&&(axis==0||axis==1)))
        throw std::invalid_argument("invalid receiver parameters");
    Mesh spin(warm,dt,baro_interval),run(measure,dt,baro_interval);
    if(run.steps<nblocks) throw std::invalid_argument("need at least 16 held steps");
    Evolution state(input,seed,axis,trial_T);
    {
        Interactions p(input.n,state.box,state.inv,input.rc);
        p.read<true>(state.x,state.scale,state.f);
        for(int k=1;k<=spin.steps;++k)
            state.advance(p,spin.dt,k==spin.steps,k%spin.baro_steps==0,spin.baro_dt,mult,0.);
    }
    std::ofstream os(argv[2]);os.exceptions(std::ios::badbit|std::ios::failbit);os<<std::setprecision(17);
    os<<"{\"source_sha256\":\""<<input.sha<<"\",\"history\":"<<input.history<<",\"seed\":"<<seed<<",\"source_temperature\":"<<input.temperature
      <<",\"temperature\":"<<trial_T
      <<",\"n\":"<<input.n<<",\"mass\":1,\"gamma\":1,\"requested_dt\":"<<dt<<",\"dt\":"<<run.dt<<",\"warm_dt\":"<<spin.dt
      <<",\"warm\":"<<warm<<",\"measure\":"<<measure<<",\"strain\":"<<eps<<",\"axis_index\":"<<axis<<",\"axis\":";matrix(os,state.axis);
    os<<",\"pressure_external\":"<<pext<<",\"barostat_interval\":"<<run.baro_dt<<",\"warm_barostat_interval\":"<<spin.baro_dt
      <<",\"volume_mobility\":"<<input.mobility*mult<<",\"warm_attempts\":"<<state.attempts<<",\"warm_accepted\":"<<state.accepts
      <<",\"prepared_reading\":";observe(state.f,state.v,state.axis,state.volume()).write(os);
    os<<",\"prepared\":";state.write_state(os);
    state.strain(eps);
    Interactions p(input.n,state.box,state.inv,input.rc);
    p.read<true>(state.x,state.scale,state.f);
    os<<",\"initial_reading\":";observe(state.f,state.v,state.axis,state.volume()).write(os);
    Tracking tracking(state.x,state.scale,p);
    state.attempts=state.accepts=0;
    os<<",\"blocks\":[";
    for(int b=0;b<nblocks;++b) {
        int start=b*int64_t(run.steps)/nblocks, end=(b+1)*int64_t(run.steps)/nblocks;
        int attempts=state.attempts, accepts=state.accepts;
        Reading sum,last;
        for(int k=start+1;k<=end;++k) {
            state.advance(p,run.dt,true,k%run.baro_steps==0,run.baro_dt,mult,pext);
            last=observe(state.f,state.v,state.axis,state.volume());sum.add(last);
        }
        if(b) os<<',';
        os<<"{\"steps\":"<<end-start<<",\"initial_t\":"<<start*run.dt<<",\"end_t\":"<<end*run.dt
          <<",\"attempts\":"<<state.attempts-attempts<<",\"accepted\":"<<state.accepts-accepts<<",\"mean\":";sum.write(os,end-start);
        os<<",\"end_reading\":";last.write(os);os<<",\"structure\":";tracking.write(os,state.x,state.scale,p);
        if(snapshots) {os<<",\"state\":";state.write_state(os);}
        os<<'}';
    }
    os<<"],\"final\":";state.write_state(os);os<<"}\n";os.close();
    std::cout<<"completed "<<spin.steps<<" zero-pressure equilibration and "<<run.steps<<" measured inertial steps; "
             <<state.accepts<<" of "<<state.attempts<<" log-volume proposals accepted\n";
    return 0;
}
} // namespace
int main(int argc,char** argv) {
    try {return execute(argc,argv);} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
