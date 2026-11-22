// Inertial response of a four-dimensional Morse material meeting an ordered layer.
// The neighbor list, complete periodic-image sum and switched pair potential
// are shared with the archived thermal-cell calculation.
#include "mouth_bulk_thermal_cell.cpp"
#include "mouth_dense_budget.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

namespace {
constexpr int count = 192, fixed_count = 32, mobile_count = count - fixed_count;
constexpr double density = 1.65, hot = 1.024, gamma_bath = 1.;
const double spacing = std::pow(.5/density, .25);
const Vec box = {4*spacing, 4*spacing, 4*spacing, 6*spacing};
const double volume = box[0]*box[1]*box[2]*box[3];
const double morse[] = {1., 7., 1.};
using Sites = std::array<Vec,count>;
using Positions = std::array<double,4*count>;

void write_vec(std::ostream& os, const Vec& v) {
    os << '[' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << ']';
}
void add_scaled(Vec& dest, const Vec& source, double scale) {
    for (int k=0; k<4; ++k) dest[k] += scale*source[k];
}
Vec scaled(const Vec& vector, double scale) {
    Vec result{};
    add_scaled(result,vector,scale);
    return result;
}
double periodic(double w) {return w-box[3]*std::floor(w/box[3]);}
bool plus(double w) {return periodic(w)<box[3]/2;}
int zone(double w) {
    // Four zones: -far, -near, +near, +far. Near means <1.5 lattice steps from the held plane.
    const double u=w-box[3]*std::floor(w/box[3]+.5);
    return u< -1.5*spacing ? 0 : u<0 ? 1 : u<1.5*spacing ? 2 : 3;
}
std::array<double,16> diagonal(bool inverse) {
    std::array<double,16> a{};
    for (int k=0;k<4;++k) a[k*5]=inverse ? 1/box[k] : box[k];
    return a;
}

struct State {
    Positions x{};
    Sites v{};
    explicit State(std::mt19937_64& rng, std::normal_distribution<double>& normal, double cold) {
        int i=0;
        for (int w=0;w<6;++w) for (int x0=0;x0<4;++x0)
        for (int y=0;y<4;++y) for (int z=0;z<4;++z) {
            if ((x0+y+z+w)%2) continue;
            x[4*i]=x0*spacing; x[4*i+1]=y*spacing; x[4*i+2]=z*spacing; x[4*i+3]=w*spacing;
            if (i>=fixed_count)
                for (double& u:v[i]) u=normal(rng)*std::sqrt(plus(w*spacing)?hot:cold);
            ++i;
        }
        if (i!=count) throw std::runtime_error("lattice count mismatch");
    }
};

struct Forces {
    Sites mobile{}, from_sheet{};
    Vec wall{}, tagged{}, virial_diag{}, pair_at_zero{}, pair_at_half{};
    double energy=0.;
};
// F_on_i = U'(r) (r_j + shift - r_i)/r for each unique bond/image.
void force_read(Cell& cell, const State& s, Forces& result, bool measure) {
    result={};
    cell.rebuild_if_needed(s.x.data());
    std::array<int,count> at_zero{}, at_half{};
    for (int i=0; i<count; ++i) {
        const double w=s.x[4*i+3];
        at_zero[i]=int(std::floor(w/box[3]));
        at_half[i]=int(std::floor(w/box[3]-.5));
    }
    // All periods exceed the cutoff; no constituent interacts with its own image.
    for (int i=0;i<count;++i) {
        for (int h=cell.first[i];h<cell.first[i+1];++h) {
            const auto& group=cell.groups[h];
            const int tile_w=int(std::lround(group.shift[3]/box[3]));
            Vec origin{};
            for (int k=0;k<4;++k) origin[k]=group.shift[k]-s.x[4*i+k];
            for (int p=group.begin;p<group.end;++p) {
                const int j=cell.neighbors[p];
                Vec d{}, f{};
                double r2=0.;
                for (int k=0;k<4;++k) {d[k]=s.x[4*j+k]+origin[k]; r2+=d[k]*d[k];}
                if (r2>=6.25) continue;
                if (!(r2>0.)) throw std::runtime_error("coincident sites");
                double e=0.;
                const double g=measure ? cell.potential.force_prefactor<true>(r2,e)
                                       : cell.potential.force_prefactor<false>(r2,e);
                for (int k=0;k<4;++k) f[k]=g*d[k];
                if (i>=fixed_count) {
                    add_scaled(result.mobile[i],f,1.);
                    add_scaled(result.mobile[j],f,-1.);
                    if (measure) {
                        const int n0=at_zero[j]+tile_w-at_zero[i];
                        const int nh=at_half[j]+tile_w-at_half[i];
                        if (n0) add_scaled(result.pair_at_zero,f,n0);
                        if (nh) add_scaled(result.pair_at_half,f,nh);
                    }
                } else if (j>=fixed_count) {
                    add_scaled(result.mobile[j],f,-1.);
                    if (measure) {
                        add_scaled(result.from_sheet[j],f,-1.);
                        add_scaled(result.wall,f,1.);
                        if (i==0) add_scaled(result.tagged,f,1.);
                    }
                }
                if (measure) {
                    result.energy+=e;
                    for (int k=0;k<4;++k) result.virial_diag[k]-=f[k]*d[k];
                }
            }
        }
    }
}
struct Union {
    std::array<int,count> parent{}, size{};
    Union() {std::iota(parent.begin(),parent.end(),0); size.fill(1);}
    int root(int i) {while(i!=parent[i]) i=parent[i]; return i;}
    void link(int i,int j) {
        i=root(i); j=root(j);
        if(i!=j) {parent[j]=i; size[i]+=size[j];}
    }
    int largest(int first) {
        int m=0;
        for(int i=first;i<count;++i) m=std::max(m,size[root(i)]);
        return m;
    }
};
struct SiteInfo {
    std::set<int> short_bonds, long_bonds;
    std::array<int,4> nz{};
    std::array<double,4> kinetic{}, tang_bragg{}, norm_bragg{};
    int cluster_all=0, cluster_mobile=0;
};
SiteInfo sites_read(const State& s) {
    SiteInfo data;
    std::array<Vec,4> co{}, si{};
    Union connected, free_connected;
    for (int i=0;i<count;++i) {
        if(i>=fixed_count) {
            int z=zone(s.x[4*i+3]);
            ++data.nz[z];
            for(int k=0;k<4;++k) {
                double phase=(2*std::acos(-1.)/spacing)*s.x[4*i+k];
                co[z][k]+=std::cos(phase); si[z][k]+=std::sin(phase);
                data.kinetic[z]+=s.v[i][k]*s.v[i][k];
            }
        }
        for (int j=i+1;j<count;++j) {
            double d2=0.;
            for (int k=0;k<4;++k) {
                double d=s.x[4*j+k]-s.x[4*i+k];
                d-=box[k]*std::round(d/box[k]);
                d2+=d*d;
            }
            if (d2>=1.30*1.30) continue;
            connected.link(i,j);
            if (i>=fixed_count) free_connected.link(i,j);
            if(j>=fixed_count) {
                const int id=count*i+j;
                data.long_bonds.insert(id);
                if(d2<1.15*1.15) data.short_bonds.insert(id);
            }
        }
    }
    for(int z=0;z<4;++z) if(data.nz[z]) {
        data.kinetic[z]/=4.*data.nz[z];
        for(int k=0;k<4;++k) {
            const double b=(co[z][k]*co[z][k]+si[z][k]*si[z][k])/(data.nz[z]*data.nz[z]);
            if(k==3) data.norm_bragg[z]=b; else data.tang_bragg[z]+=b/3.;
        }
    }
    data.cluster_all=connected.largest(0);
    data.cluster_mobile=free_connected.largest(fixed_count);
    return data;
}
void write_sites(std::ostream& os,const SiteInfo& s,const SiteInfo& ref) {
    std::array<int,2> orig{},lost{},gained{};
    for(int id:ref.short_bonds) {int k=id/count>=fixed_count; ++orig[k]; if(!s.long_bonds.count(id)) ++lost[k];}
    for(int id:s.short_bonds) if(!ref.long_bonds.count(id)) ++gained[id/count>=fixed_count];
    os<<"{\"largest_connected_all\":"<<s.cluster_all<<",\"largest_connected_mobile\":"<<s.cluster_mobile
      <<",\"contact_columns\":[\"held_mobile\",\"mobile_mobile\"],\"starting_contacts\":["<<orig[0]<<','<<orig[1]
      <<"],\"starting_contacts_lost\":["<<lost[0]<<','<<lost[1]<<"],\"new_contacts\":["<<gained[0]<<','<<gained[1]
      <<"],\"zones\":[";
    for(int z=0;z<4;++z) {
        if(z) os<<',';
        os<<"{\"n\":"<<s.nz[z]<<",\"kinetic_T\":"<<s.kinetic[z]<<",\"tangent_bragg_intensity\":"<<s.tang_bragg[z]
          <<",\"normal_bragg_intensity\":"<<s.norm_bragg[z]<<'}';
    }
    os<<"]}";
}
void write_state(std::ostream& os,const State& state) {
    os<<"{\"x\":[";
    for(int i=0;i<count;++i) {
        if(i) os<<',';
        Vec x; std::copy_n(state.x.data()+4*i,4,x.begin()); write_vec(os,x);
    }
    os<<"],\"v\":[";
    for(int i=0;i<count;++i) {if(i) os<<',';write_vec(os,state.v[i]);}
    os<<"]}";
}

struct Block {
    dense_budget::Budget ledger;
    Vec wall{},tagged{},virial{};
    std::array<double,2> bath_work{};
    std::array<Vec,4> travel{};
    std::array<double,4> residence{};
    std::vector<std::array<double,7>> crossings;
};
void boundary(Block& b,const State& s,bool start) {
    for(int i=fixed_count;i<count;++i) {
        if(start) b.ledger.start(i,s.x[4*i+3],s.v[i],box[3]);
        else b.ledger.finish(i,s.x[4*i+3],s.v[i],box[3]);
    }
}
void kick(State& s,const Forces& f,double h,const Vec& drive,Block* b) {
    if(b) {
        b->ledger.pair_kick(1,0,scaled(f.pair_at_zero,h));
        b->ledger.pair_kick(0,1,scaled(f.pair_at_half,h));
        add_scaled(b->wall,f.wall,h);
        add_scaled(b->tagged,f.tagged,h);
        add_scaled(b->virial,f.virial_diag,h);
    }
    const Vec bd=scaled(drive,h);
    for(int i=fixed_count;i<count;++i) {
        add_scaled(s.v[i],f.mobile[i],h);
        add_scaled(s.v[i],drive,h);
        if(b) b->ledger.kick_mobile(s.x[4*i+3],bd,scaled(f.from_sheet[i],h),box[3]);
    }
}
void drift(State& s,double h,Block* b,int step,int half) {
    for(int i=fixed_count;i<count;++i) {
        const double a=s.x[4*i+3];
        for(int k=0;k<4;++k) s.x[4*i+k]+=h*s.v[i][k];
        if(!b) continue;
        const double w=s.x[4*i+3];
        b->ledger.drift(a,w,s.v[i],box[3]);
        for(int c=0;c<2;++c) {
            int n=dense_budget::crossings(a,w,.5*c*box[3],box[3]);
            if(n) b->crossings.push_back({double(step),double(half),double(i),double(c),double(n),s.v[i][0],s.v[i][3]});
        }
        // Midpoint region quadrature for drift diagnostics; the main crossing ledger is exact.
        const int z=zone(.5*(a+w));
        b->residence[z]+=h;
        add_scaled(b->travel[z],s.v[i],h);
        for(int k=0;k<4;++k) b->virial[k]+=h*s.v[i][k]*s.v[i][k];
    }
}
void evolve(State& s,Cell& cell,Forces& f,std::mt19937_64& rng,std::normal_distribution<double>& noise,
            double dt,double cold,const Vec& drive,Block* b,int step) {
    double h=.5*dt;
    kick(s,f,h,drive,b);
    drift(s,h,b,step,0);
    const double c=std::exp(-gamma_bath*dt), variance=-std::expm1(-2*gamma_bath*dt);
    for(int i=fixed_count;i<count;++i) {
        double w=s.x[4*i+3];
        double sigma=std::sqrt(variance*(plus(w)?hot:cold));
        Vec impulse{};
        double bath_work=0.;
        for(int k=0;k<4;++k) {
            double next=c*s.v[i][k]+sigma*noise(rng);
            impulse[k]=next-s.v[i][k];
            bath_work+=.5*(next+s.v[i][k])*impulse[k];
            s.v[i][k]=next;
        }
        if(b) {
            b->ledger.bath(w,impulse,box[3]);
            b->bath_work[plus(w)]+=bath_work;
        }
    }
    drift(s,h,b,step,1);
    force_read(cell,s,f,b!=nullptr);
    kick(s,f,h,drive,b);
    for(int i=fixed_count;i<count;++i) for(int k=0;k<4;++k)
        if(!std::isfinite(s.v[i][k])||!std::isfinite(s.x[4*i+k]))
            throw std::runtime_error("nonfinite material state");
}
void write_block(std::ostream& os,const Block& b,int steps,double dt) {
    os<<"{\"steps\":"<<steps<<",\"ledger\":"; b.ledger.json(os);
    os<<",\"sheet_mobile_force_impulse\":"; write_vec(os,b.wall);
    os<<",\"site0_mobile_force_impulse\":"; write_vec(os,b.tagged);
    os<<",\"pressure_diagonal\":"; write_vec(os,scaled(b.virial,1/(volume*steps*dt)));
    os<<",\"bath_kinetic_energy_delivered_minus_plus\":["<<b.bath_work[0]<<','<<b.bath_work[1]<<']';
    os<<",\"drift_zones\":[";
    for(int z=0;z<4;++z) {
        if(z) os<<',';
        os<<"{\"residence_site_time\":"<<b.residence[z]<<",\"integrated_displacement\":"; write_vec(os,b.travel[z]); os<<'}';
    }
    os<<"],\"events\":[";
    bool first=true;
    for(const auto& event:b.crossings) {
        if(!first) os<<',';
        first=false; os<<'[';
        for(int j=0;j<7;++j) {if(j) os<<',';os<<event[j];} os<<']';
    }
    os<<"]}";
}
int step_count(double t,double dt) {
    const auto n=std::lround(t/dt);
    if(n<=0 || n>2000000000L || std::abs(n*dt-t)>1.e-9*std::max(1.,t))
        throw std::invalid_argument("durations must be positive integer multiples of dt");
    return int(n);
}
} // namespace

int main(int argc,char** argv) {
    try {
        if(argc<6||argc>9) throw std::invalid_argument("OUT.json SEED DT AX AW [WARM=3] [MEASURE=3] [T_MINUS=1.024]");
        const auto seed=std::stoull(argv[2]);
        const double dt=std::stod(argv[3]), ax=std::stod(argv[4]), aw=std::stod(argv[5]);
        const double warm=argc>6?std::stod(argv[6]):3., measure=argc>7?std::stod(argv[7]):3.;
        const double cold=argc>8?std::stod(argv[8]):hot;
        if(!(dt>0) || !(cold>0) || !std::isfinite(ax)||!std::isfinite(aw)) throw std::invalid_argument("invalid physical input");
        const int nw=step_count(warm,dt), nm=step_count(measure,dt);
        if(nm<4) throw std::invalid_argument("need four nonempty measurement blocks");
        std::mt19937_64 rng(seed); std::normal_distribution<double> noise;
        State s(rng,noise,cold); Cell cell(count,diagonal(false).data(),diagonal(true).data(),2.5,morse);
        Forces f;
        force_read(cell,s,f,false);
        for(int step=0;step<nw;++step) evolve(s,cell,f,rng,noise,dt,cold,Vec{},nullptr,step);
        const State prep=s;
        const SiteInfo ref=sites_read(prep);
        force_read(cell,s,f,true);
        const double estart=f.energy; const Vec fstart=f.wall;
        std::ofstream os(argv[1]);
        os.exceptions(std::ios::failbit|std::ios::badbit);
        os<<std::setprecision(17)<<"{\"version\":1,\"n_sites\":"<<count<<",\"n_held\":"<<fixed_count<<",\"lattice_spacing\":"<<spacing
          <<",\"box\":"; write_vec(os,box); os<<",\"density\":"<<density<<",\"seed\":"<<seed<<",\"dt\":"<<dt
          <<",\"warm_time\":"<<warm<<",\"measurement_time\":"<<measure<<",\"drive\":["<<ax<<",0,0,"<<aw
          <<"],\"T_plus\":"<<hot<<",\"T_minus\":"<<cold<<",\"gamma\":"<<gamma_bath<<",\"mass\":1"
          <<",\"zone_names\":[\"minus_far\",\"minus_near\",\"plus_near\",\"plus_far\"]"
          <<",\"event_columns\":[\"measurement_step_zero_based\",\"half_drift_zero_based\",\"site\",\"plane_zero_or_mid\",\"signed_crossings\",\"vx\",\"vw\"]"
          <<",\"initial\":"; write_state(os,prep); os<<",\"prepared_structure\":"; write_sites(os,ref,ref);
        os<<",\"initial_pair_energy\":"<<estart<<",\"initial_sheet_force\":"; write_vec(os,fstart); os<<",\"blocks\":[";
        for(int b=0;b<4;++b) {
            const int begin=b*nm/4, end=(b+1)*nm/4;
            Block block; boundary(block,s,true);
            for(int step=begin;step<end;++step) evolve(s,cell,f,rng,noise,dt,cold,{ax,0,0,aw},&block,step);
            boundary(block,s,false);
            if(b) os<<',';
            write_block(os,block,end-begin,dt);
        }
        os<<"],\"final\":"; write_state(os,s); os<<",\"final_structure\":"; write_sites(os,sites_read(s),ref);
        os<<",\"final_pair_energy\":"<<f.energy<<",\"final_sheet_force\":"; write_vec(os,f.wall); os<<"}\n"; os.close();
        std::cout<<"completed "<<nw<<" undriven and "<<nm<<" measured steps for "<<mobile_count<<" mobile sites and "<<fixed_count<<" held sites\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
