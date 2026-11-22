// A 4D periodic neutral material with a common tangent lattice and a free normal cell row.
// Uses the archived all-image Morse potential and the receiver's physical kinetic stress.
#define main original_isotropic_receiver_main
#include "mouth_receiver_relax.cpp"
#undef main

namespace junction {
constexpr int nb=32;
struct Preparation {
    int n;double T,mu;std::string hash;Matrix box{},inv{};Sites x;std::vector<int> tags;
    explicit Preparation(const std::string& path) {
        std::ifstream in(path);
        if(!(in>>n>>T>>mu>>hash)||n<2||n>10000||T!=1.088||!(mu>0)||hash.size()!=64) throw std::runtime_error("invalid junction input header");
        for(auto* a:{&box,&inv}) for(auto& row:*a) for(auto& z:row) if(!(in>>z)||!std::isfinite(z)) throw std::runtime_error("bad cell");
        if(!(box[3][3]>0.&&determinant(box)>0.)) throw std::runtime_error("inverted junction cell");
        for(int j=0;j<3;++j) if(std::abs(box[j][3])>1e-12) throw std::runtime_error("face does not have declared w normal");
        auto check=product(box,inv);for(int i=0;i<4;++i)for(int j=0;j<4;++j) if(std::abs(check[i][j]-(i==j))>2e-11) throw std::runtime_error("inconsistent inverse");
        x.resize(4*n);tags.resize(n);
        for(auto& z:x) if(!(in>>z)||!std::isfinite(z)) throw std::runtime_error("bad sites");
        for(auto& t:tags) if(!(in>>t)||(t!=0&&t!=1)) throw std::runtime_error("bad parent IDs");
        std::string extra;if(in>>extra)throw std::runtime_error("extra fields");
    }
};
struct Shape {
    // At fixed fractional position: r=(y_x+b_x y_w, y_y+b_y y_w, y_z+b_z y_w, l y_w).
    Vec a{0,0,0,1};
    Vec map(const Vec& y) const {return {y[0]+a[0]*y[3],y[1]+a[1]*y[3],y[2]+a[2]*y[3],a[3]*y[3]};}
    double lower_stretch() const {
        double l=a[3],tr=1.+std::inner_product(a.begin(),a.end(),a.begin(),0.);
        return std::min(1.,std::sqrt(2*l*l/(tr+std::sqrt(std::max(0.,tr*tr-4*l*l)))));
    }
    Matrix box(const Matrix& reference) const {Matrix r=reference;for(auto& u:r)u=map(u);return r;}
};
int bin(double frac){int i=int(std::floor(nb*frac));return (i%nb+nb)%nb;}
template<class F> void segment(double a,double b,F apply) {
    double u=nb*std::min(a,b), v=nb*std::max(a,b), len=v-u;
    if(len<=1e-12){apply(bin((a+b)/2.),1.);return;}
    for(int i=int(std::floor(u));i<int(std::ceil(v));++i){double w=std::min(v,i+1.)-std::max(u,double(i));if(w>0.)apply((i%nb+nb)%nb,w/len);}
}
using Bins=std::array<Vec,nb>;
struct Kernel {
    Cell list;double reach=0.,height;std::vector<Vec> self;
    Kernel(const Preparation& p):list(p.n,&p.box[0][0],&p.inv[0][0],2.5,params),height(p.box[3][3]){}
    void ensure(const Sites& x,const Shape& shape){
        const double need=2.5/shape.lower_stretch();
        if(!(need>0) || !std::isfinite(need))throw std::runtime_error("degenerate trial geometry");
        if(need>reach || need<.83*reach){reach=1.05*need;list.reference.clear();self.clear();list.each_image(Vec{},reach,[&](const Tile& t,const Vec& s,double){if(t>Tile{})self.push_back(s);});}
        list.potential.cutoff=reach;list.rebuild_if_needed(x.data());list.potential.cutoff=2.5;
    }
    template<bool details,bool profile=false> void read(const Sites& x,const Shape& shape,Forces& f,Bins* assigned=nullptr){
        ensure(x,shape);std::fill(f.grad.begin(),f.grad.end(),0.);f.energy=0.;f.virial={};if constexpr(profile) *assigned={};
        for(int i=0;i<list.n;++i){
            Vec gi{};const double wi=x[4*i+3]/height;
            for(int g=list.first[i];g<list.first[i+1];++g){const auto& gr=list.groups[g];Vec o;for(int k=0;k<4;++k)o[k]=gr.shift[k]-x[4*i+k];
                for(int b=gr.begin;b<gr.end;++b){const int j=list.neighbors[b];Vec y;for(int k=0;k<4;++k)y[k]=x[4*j+k]+o[k];
                    Vec r=shape.map(y);double rr=std::inner_product(r.begin(),r.end(),r.begin(),0.);
                    if(rr>=6.25) continue;
                    if(!(rr>0.))throw std::runtime_error("coincident physical sites");
                    double e=0.;const double dr=list.potential.force_prefactor<details>(rr,e);
                    if constexpr(details)f.energy+=e;
                    for(int k=0;k<4;++k){const double v=dr*r[k];f.grad[4*j+k]+=v;gi[k]-=v;if constexpr(details)for(int l=0;l<4;++l)f.virial[k][l]+=v*r[l];}
                    if constexpr(profile)segment(wi,wi+y[3]/height,[&](int h,double portion){for(int k=0;k<4;++k)(*assigned)[h][k]+=portion*dr*r[k]*r[3];});
                }}
            for(int k=0;k<4;++k)f.grad[4*i+k]+=gi[k];
        }
        if constexpr(details)for(const Vec& y:self){Vec r=shape.map(y);double rr=std::inner_product(r.begin(),r.end(),r.begin(),0.);if(rr>=6.25)continue;double e;double g=list.potential.force_prefactor<true>(rr,e);f.energy+=e*list.n;
            for(int k=0;k<4;++k)for(int l=0;l<4;++l)f.virial[k][l]+=g*list.n*r[k]*r[l];
            if constexpr(profile)for(int i=0;i<list.n;++i){const double wi=x[4*i+3]/height;segment(wi,wi+y[3]/height,[&](int h,double portion){for(int k=0;k<4;++k)(*assigned)[h][k]+=portion*g*r[k]*r[3];});}}
    }
};
struct Control {
    Bins initial{},final{},pair{},bath{},forward{},reverse{};
    std::array<int,nb> n0{},n1{};std::array<std::array<long long,2>,nb> npos{},nneg{};
    static int wrap(int a){return (a%nb+nb)%nb;}
    void record(double a,double b,const double* v,int tag){
        int start=int(std::floor(nb*a)),end=int(std::floor(nb*b));
        for(int k=start+1;k<=end;++k){int h=wrap(k);++npos[h][tag];for(int j=0;j<4;++j)forward[h][j]+=v[j];}
        for(int k=start;k>end;--k){int h=wrap(k);++nneg[h][tag];for(int j=0;j<4;++j)reverse[h][j]+=v[j];}
    }
    void endpoint(const Sites& y,const Sites& v,double H,bool first){
        auto& p=first?initial:final;auto& n=first?n0:n1;
        for(size_t i=0;i<y.size();i+=4){int h=bin(y[i+3]/H);++n[h];for(int j=0;j<4;++j)p[h][j]+=v[i+j];}
    }
    void output(std::ostream& os)const{
        os<<'[';for(int h=0;h<nb;++h){int next=(h+1)%nb;Vec residual{};long long adv=0;
            for(int k=0;k<4;++k)residual[k]=final[h][k]-initial[h][k]-pair[h][k]-bath[h][k]-forward[h][k]+reverse[h][k]+forward[next][k]-reverse[next][k];
            for(int tag=0;tag<2;++tag)adv+=npos[h][tag]-nneg[h][tag]-npos[next][tag]+nneg[next][tag];
            if(n1[h]-n0[h]!=adv)throw std::runtime_error("slab mass accounting failed");
            for(auto e:residual)if(std::abs(e)>1e-6)throw std::runtime_error("slab momentum accounting failed");
            if(h)os<<',';os<<"{\"n0\":"<<n0[h]<<",\"n1\":"<<n1[h]<<",\"p0\":";vec(os,initial[h]);os<<",\"p1\":";vec(os,final[h]);os<<",\"pair\":";vec(os,pair[h]);os<<",\"bath\":";vec(os,bath[h]);os<<",\"cross_plus\":";vec(os,npos[h]);os<<",\"cross_minus\":";vec(os,nneg[h]);os<<",\"p_plus\":";vec(os,forward[h]);os<<",\"p_minus\":";vec(os,reverse[h]);os<<",\"residual\":";vec(os,residual);os<<'}';}os<<']';
    }
};
struct Profile {Bins traction{};std::array<double,nb> rho{},kin{},pop0{},pop1{};double weight=0.;
    void add(const Preparation& p,const Sites& y,const Sites& v,const Bins& vir,double volume,double time){
        const double multiplier=time*nb/volume;weight+=time;
        for(int h=0;h<nb;++h)for(int k=0;k<4;++k)traction[h][k]-=multiplier*vir[h][k];
        for(int i=0;i<p.n;++i){int h=bin(y[4*i+3]/p.box[3][3]);rho[h]+=multiplier;(p.tags[i]?pop1:pop0)[h]+=time;
            for(int j=0;j<4;++j){double vk=v[4*i+j];kin[h]+=time*vk*vk/4;traction[h][j]+=multiplier*vk*v[4*i+3];}}
    }
    void output(std::ostream& os)const{os<<'[';for(int h=0;h<nb;++h){if(h)os<<',';os<<"{\"rho\":"<<rho[h]/weight<<",\"n_ordered_parent\":"<<pop0[h]/weight<<",\"n_receiver_parent\":"<<pop1[h]/weight<<",\"kinetic_T\":";
        if(pop0[h]+pop1[h]>0)os<<kin[h]/(pop0[h]+pop1[h]);else os<<"null";os<<",\"pressure_w\":";auto t=traction[h];for(auto& z:t)z/=weight;vec(os,t);os<<'}';}os<<']';}
};
struct Dynamics {
    const Preparation& p;Sites y,v;Shape sh;Kernel mat;Forces f,trial;double v0;
    std::mt19937_64 thermal_rng,box_rng;std::normal_distribution<double> site_noise,box_noise;std::uniform_real_distribution<double> uni;
    int attempts=0,accepted=0;double max_abs_logH=0.,max_tilt=0.,box_energy=0.,bath_heat=0.;
    explicit Dynamics(const Preparation& input,uint64_t seed):p(input),y(p.x),v(p.x.size()),mat(p),f(p.n),trial(p.n),v0(determinant(p.box)),thermal_rng(seed),box_rng(seed^0xb3248f63ULL){for(auto& u:v)u=std::sqrt(p.T)*site_noise(thermal_rng);mat.read<true>(y,sh,f);}
    double volume() const{return v0*sh.a[3];}
    void kick(double h,Control* c){for(int i=0;i<p.n;++i){int m=bin(y[4*i+3]/mat.height);for(int k=0;k<4;++k){double impulse=-h*f.grad[4*i+k];v[4*i+k]+=impulse;if(c)c->pair[m][k]+=impulse;}}}
    void drift(double h,Control* c){for(int i=0;i<p.n;++i){int n=4*i;double old=y[n+3]/mat.height;double dw=h*v[n+3]/sh.a[3];for(int j=0;j<3;++j)y[n+j]+=h*v[n+j]-dw*sh.a[j];y[n+3]+=dw;if(c)c->record(old,y[n+3]/mat.height,v.data()+n,p.tags[i]);}}
    double trial_ratio(const Shape& proposed,double pext){
        mat.read<true>(y,proposed,trial);const double dq=std::log(proposed.a[3]/sh.a[3]);
        double result=(p.n+1)*dq-(trial.energy-f.energy+pext*volume()*std::expm1(dq))/p.T;
        if(!std::isfinite(result))throw std::runtime_error("nonfinite normal-cell step");return result;
    }
    bool decide(const Shape& proposed,double log_ratio,double logu){
        ++attempts;bool take=logu<log_ratio;
        if(take){box_energy+=trial.energy-f.energy;sh=proposed;std::swap(f,trial);++accepted;}
        max_abs_logH=std::max(max_abs_logH,std::abs(std::log(sh.a[3])));
        for(int i=0;i<3;++i)max_tilt=std::max(max_tilt,std::abs(sh.a[i]));return take;
    }
    void attempt_box(double delay,double mult,double pext){
        const double sigma=std::sqrt(2*p.mu*(128./p.n)*mult*p.T*delay),dq=sigma*box_noise(box_rng);
        Shape proposed=sh;proposed.a[3]*=std::exp(dq);for(int j=0;j<3;++j)proposed.a[j]+=sigma*box_noise(box_rng);
        const double ratio=trial_ratio(proposed,pext);decide(proposed,ratio,std::log(uni(box_rng)));
    }
    void advance(double dt,bool box_step,double baro_dt,double mult,double pext,Control* ctrl){
        kick(.5*dt,ctrl);drift(.5*dt,ctrl);double c=std::exp(-dt),sigma=std::sqrt(-std::expm1(-2*dt)*p.T);
        for(int i=0;i<p.n;++i){int h=bin(y[4*i+3]/mat.height);for(int j=0;j<4;++j){int k=4*i+j;double previous=v[k];v[k]=c*previous+sigma*site_noise(thermal_rng);bath_heat+=.5*(v[k]-previous)*(v[k]+previous);if(ctrl)ctrl->bath[h][j]+=v[k]-previous;}}
        drift(.5*dt,ctrl);if(box_step)mat.read<true>(y,sh,f);else mat.read<false>(y,sh,f);
        kick(.5*dt,ctrl);if(box_step)attempt_box(baro_dt,mult,pext);
        for(size_t k=0;k<y.size();++k)if(!std::isfinite(y[k])||!std::isfinite(v[k]))throw std::runtime_error("nonfinite material dynamics");
    }
    void state(std::ostream& os,bool velocities=true)const{
        os<<"{\"box\":";matrix(os,sh.box(p.box));os<<",\"pos\":[";
        for(int i=0;i<p.n;++i){if(i)os<<',';vec(os,sh.map({y[4*i],y[4*i+1],y[4*i+2],y[4*i+3]}));}os<<']';if(velocities){os<<",\"vel\":";sites(os,v,1.);}os<<'}';
    }
};

void check(const Preparation& p,const std::string& dest){Dynamics d(p,137);std::ofstream os(dest);os.exceptions(std::ios::badbit|std::ios::failbit);os<<std::setprecision(17);
    os<<"{\"n\":"<<p.n<<",\"input\":\""<<p.hash<<"\",\"U\":"<<d.f.energy<<",\"virial\":";matrix(os,d.f.virial);os<<",\"v\":";sites(os,d.v,1.);os<<",\"grad\":";sites(os,d.f.grad,1.);os<<",\"finite\":[";
    for(int k=0;k<4;++k){const double h=2.e-6;Shape a,b;
        if(k==3){a.a[3]=std::exp(h);b.a[3]=std::exp(-h);}else {a.a[k]=h/d.mat.height;b.a[k]=-h/d.mat.height;}
        Forces f1(p.n),f2(p.n);d.mat.read<true>(d.y,a,f1);d.mat.read<true>(d.y,b,f2);
        const double fd=(f1.energy-f2.energy)/(2*h),analytic=d.f.virial[k][3]/(k==3?1.:d.mat.height);
        if(k)os<<',';os<<"{\"k\":"<<k<<",\"slope\":"<<fd<<",\"virial\":"<<analytic<<'}';
        if(std::abs(fd-analytic)>4.e-6*std::max(1.,std::abs(analytic))+3.e-4)throw std::runtime_error("normal generalized-force derivative mismatch");}
    // Repeat in a tilted and expanded chart; add integer lattice vectors per site
    // to compare energy, local forces and virial, independently of their wrapping.
    Shape z;z.a={.043,-.071,.025,1.029};Forces changed(p.n),wound(p.n);d.mat.read<true>(d.y,z,changed);
    Sites moved=d.y;for(int i=0;i<p.n;++i){int tile=i%3-1;for(int j=0;j<4;++j)moved[4*i+j]+=tile*p.box[(i*3)%4][j];}
    d.mat.read<true>(moved,z,wound);
    double worst=0.;for(size_t i=0;i<moved.size();++i)worst=std::max(worst,std::abs(changed.grad[i]-wound.grad[i]));
    double ratio=(p.n+1)*std::log(z.a[3])-(changed.energy-d.f.energy)/p.T;
    double reverse=(p.n+1)*std::log(1/z.a[3])-(d.f.energy-changed.energy)/p.T;
    // Exercise the actual accepted/rejected branches using legal uniforms. No MD step is taken.
    const Sites original_y=d.y,original_v=d.v;const double firstU=d.f.energy;
    const double actual=d.trial_ratio(z,0.);
    if(!d.decide(z,actual,std::min(-1.,actual-1.))||d.y!=original_y||d.v!=original_v||d.sh.a!=z.a||std::abs(d.f.energy-changed.energy)>1.e-8)
        throw std::runtime_error("accepted change did not preserve fractions and physical momenta");
    bool rejected=false;
    for(double delta:{.01,-.01,.1,-.1,.25,-.25}){
        Shape try_z=z;try_z.a[3]*=std::exp(delta);double test=d.trial_ratio(try_z,0.);
        if(test<-.01){const auto gradient=d.f.grad;
            rejected=!d.decide(try_z,test,.5*test)&&d.sh.a==z.a&&d.y==original_y&&d.v==original_v&&d.f.grad==gradient;
            break;}}
    if(!rejected||std::abs(d.box_energy-(changed.energy-firstU))>1e-10)throw std::runtime_error("rejected proposal altered old state");
    os<<"],\"tilted_U\":"<<changed.energy<<",\"tilted_virial\":";matrix(os,changed.virial);
    os<<",\"accepted_preserved_velocities_and_fractions\":true,\"rejected_state_unchanged\":true,\"tilted_log_ratio\":"<<ratio<<",\"reverse_log_ratio\":"<<reverse<<",\"winding_dU\":"<<wound.energy-changed.energy<<",\"winding_max_dgrad\":"<<worst<<'}';
    if(std::abs(ratio+reverse)>1.e-10||std::abs(wound.energy-changed.energy)>1.e-7||worst>1.e-7)throw std::runtime_error("cell chart consistency failed");
    std::cout<<"validated four virtual force directions and tilted-cell winding\n";
}

int execute(int argc,char** argv){
    if(argc==4&&std::string(argv[1])=="check"){Preparation p(argv[2]);check(p,argv[3]);return 0;}
    if(argc<7||argc>11)throw std::invalid_argument("INPUT OUT SEED DT WARM MEASURE [BARO_INTERVAL=.08] [MOB_MULT=1] [Pext=0] [SNAPSHOT_CADENCE=0]");
    Preparation p(argv[1]);auto seed=std::stoull(argv[3]);const double dt=std::stod(argv[4]),warm=std::stod(argv[5]),time=std::stod(argv[6]);double baro=argc>7?std::stod(argv[7]):.08,mult=argc>8?std::stod(argv[8]):1.,pext=argc>9?std::stod(argv[9]):0.,stride=argc>10?std::stod(argv[10]):0.;
    if(!(dt>0&&warm>=0&&time>=16*dt&&warm/dt<2.e9&&time/dt<2.e9&&baro>=dt&&baro/dt<2.e9&&mult>0&&std::isfinite(mult)&&std::isfinite(pext)&&stride>=0.&&std::isfinite(stride)))throw std::invalid_argument("invalid junction parameters");
    Mesh spin(warm,dt,baro),run(time,dt,baro);Dynamics d(p,seed);Matrix xw{};xw[0][3]=xw[3][0]=.5;
    Reading input_reading=observe(d.f,d.v,xw,d.volume());
    for(int k=1;k<=spin.steps;++k)d.advance(spin.dt,k%spin.baro_steps==0,spin.baro_dt,mult,pext,nullptr);
    d.mat.read<true>(d.y,d.sh,d.f);Reading warm_reading=observe(d.f,d.v,xw,d.volume());
    std::ofstream os(argv[2]);os.exceptions(std::ios::badbit|std::ios::failbit);os<<std::setprecision(17);
    os<<"{\"input_hash\":\""<<p.hash<<"\",\"n\":"<<p.n<<",\"T\":"<<p.T<<",\"seed\":"<<seed<<",\"dt\":"<<run.dt<<",\"warm\":"<<warm<<",\"measure\":"<<time<<",\"baro_dt\":"<<run.baro_dt<<",\"mob_multiplier\":"<<mult<<",\"logH_mobility\":"<<p.mu*128./p.n*mult<<",\"Pext\":"<<pext<<",\"bins\":"<<nb<<",\"origin\":";vec(os,p.tags);os<<",\"input_reading\":";input_reading.write(os);os<<",\"prepared_reading\":";warm_reading.write(os);os<<",\"warm_bath_heat\":"<<d.bath_heat<<",\"warm_cell_dU\":"<<d.box_energy<<",\"prepared\":";d.state(os);os<<",\"blocks\":[";
    int sample_stride=std::max(1,int(std::lround(.125/run.dt)));int pos_stride=stride>0?std::max(1,int(std::lround(stride/run.dt))):0;
    for(int b=0;b<16;++b){int start=b*int64_t(run.steps)/16,end=(b+1)*int64_t(run.steps)/16;
        Control ctrl;ctrl.endpoint(d.y,d.v,d.mat.height,true);Profile pr;Reading sum,last;int represented=0,samples=0,a0=d.attempts,ac0=d.accepted;Bins assign; double bath0=d.bath_heat,cell0=d.box_energy;
        if(b)os<<',';os<<"{\"steps\":"<<end-start<<",\"t_start\":"<<start*run.dt<<",\"t_end\":"<<end*run.dt<<",\"snapshots\":[";bool first=true;
        for(int k=start+1;k<=end;++k){d.advance(run.dt,k%run.baro_steps==0,run.baro_dt,mult,pext,&ctrl);
            if((k-start)%sample_stride==0||k==end){int weight=k-start-represented;d.mat.read<true,true>(d.y,d.sh,d.f,&assign);last=observe(d.f,d.v,xw,d.volume());
                // Global and local tractions use exactly the same right-endpoint sampling times.
                for(int rep=0;rep<weight;++rep)sum.add(last);
                pr.add(p,d.y,d.v,assign,d.volume(),weight*run.dt);represented+=weight;++samples;}
            if(pos_stride && k%pos_stride==0 && k!=end){if(!first)os<<',';first=false;os<<"{\"t\":"<<k*run.dt<<",\"state\":";d.state(os,false);os<<'}';}}
        ctrl.endpoint(d.y,d.v,d.mat.height,false);os<<"],\"attempts\":"<<d.attempts-a0<<",\"accepted\":"<<d.accepted-ac0<<",\"bath_heat\":"<<d.bath_heat-bath0<<",\"cell_dU\":"<<d.box_energy-cell0<<",\"stress_samples\":"<<samples<<",\"mean\":";sum.write(os,represented);os<<",\"end_reading\":";last.write(os);os<<",\"state\":";d.state(os);os<<",\"profile\":";pr.output(os);os<<",\"balance\":";ctrl.output(os);os<<'}';}
    os<<"],\"max_abs_log_height_ratio\":"<<d.max_abs_logH<<",\"max_unwrapped_tilt_over_initial_H\":"<<d.max_tilt<<",\"final\":";d.state(os);os<<"}\n";os.close();
    std::cout<<"completed "<<spin.steps<<" warm, "<<run.steps<<" measured; "<<d.accepted<<" / "<<d.attempts<<" cell moves accepted\n";return 0;
}
} // junction
int main(int argc,char** argv){try{return junction::execute(argc,argv);}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
