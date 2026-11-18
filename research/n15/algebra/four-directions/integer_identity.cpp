// 在互异整数参数上验证六个整系数恒等式，并累计商群阶数的严格下界。
#include <gmp.h>
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class Integer {
public:
    mpz_t z;
    Integer() { mpz_init(z); }
    Integer(const Integer& other) { mpz_init_set(z,other.z); }
    Integer& operator=(const Integer& other) { mpz_set(z,other.z); return *this; }
    ~Integer() { mpz_clear(z); }
};

struct Gate { int op,a,b; };
struct Input {
    std::vector<Integer> p, dp, denom, scaling;
    Integer common;
    std::vector<std::vector<Integer>> numerators;
    std::vector<Gate> circuit;
    std::vector<int> outputs;
};

Integer read_integer(std::istream& in) {
    std::string s;
    Integer x;
    if(!(in>>s) || mpz_set_str(x.z,s.c_str(),10)) throw std::runtime_error("无效整数");
    return x;
}
std::vector<Integer> polynomial(std::istream& in) {
    size_t n;
    if(!(in>>n) || n==0 || n>10000) throw std::runtime_error("无效系数个数");
    std::vector<Integer> v;
    v.reserve(n);
    for(size_t i=0;i<n;++i) v.push_back(read_integer(in));
    return v;
}
Integer evaluate(const std::vector<Integer>& v,unsigned long k,const Integer& stride) {
    Integer a,part;
    size_t n=v.size(),length=n%24;
    if(!length) length=24;
    while(n) {
        mpz_set_ui(part.z,0);
        for(size_t j=0;j<length;++j) {
            mpz_mul_ui(part.z,part.z,k);
            mpz_add(part.z,part.z,v[--n].z);
        }
        mpz_mul(a.z,a.z,stride.z);
        mpz_add(a.z,a.z,part.z);
        length=24;
    }
    return a;
}
void strip(Integer& modulus,const Integer& forbidden) {
    Integer g;
    mpz_gcd(g.z,modulus.z,forbidden.z);
    while(mpz_cmp_ui(g.z,1)>0) {
        mpz_divexact(modulus.z,modulus.z,g.z);
        mpz_mul(g.z,g.z,g.z);
        mpz_gcd(g.z,modulus.z,g.z);
    }
}

long verify(const Input& in,unsigned long k,const Integer& forbidden) {
    Integer stride;
    mpz_ui_pow_ui(stride.z,k,24);
    Integer mod=evaluate(in.p,k,stride), derivative=evaluate(in.dp,k,stride),temp;
    mpz_abs(mod.z,mod.z);
    if(mpz_cmp_ui(mod.z,2)<0) return -20;
    strip(mod,forbidden);
    strip(mod,derivative);
    if(mpz_cmp_ui(mod.z,2)<0) return 0;
    std::vector<Integer> x(6);
    mpz_set_ui(x[0].z,k);
    mpz_mod(x[0].z,x[0].z,mod.z);
    mpz_mul(temp.z,in.common.z,derivative.z);
    if(!mpz_invert(temp.z,temp.z,mod.z)) return -21;
    for(size_t i=0;i<5;++i) {
        x[i+1]=evaluate(in.numerators[i],k,stride);
        mpz_mul(x[i+1].z,x[i+1].z,in.scaling[i].z);
        mpz_mod(x[i+1].z,x[i+1].z,mod.z);
        mpz_mul(x[i+1].z,x[i+1].z,temp.z);
        mpz_mod(x[i+1].z,x[i+1].z,mod.z);
    }
    std::vector<Integer> v(in.circuit.size());
    for(size_t i=0;i<v.size();++i) {
        const auto& t=in.circuit[i];
        auto z=v[i].z;
        if(t.op==0) mpz_set(z,x[t.a].z);
        else if(t.op==1) { mpz_set_si(z,t.a); mpz_mod(z,z,mod.z); }
        else if(t.op==2) {
            if(mpz_sgn(v[t.a].z)) mpz_sub(z,mod.z,v[t.a].z);
        } else if(t.op==3) {
            mpz_add(z,v[t.a].z,v[t.b].z);
            if(mpz_cmp(z,mod.z)>=0) mpz_sub(z,z,mod.z);
        } else {
            mpz_mul(z,v[t.a].z,v[t.b].z);
            mpz_mod(z,z,mod.z);
        }
    }
    for(size_t i=0;i<in.outputs.size();++i)
        if(mpz_sgn(v[in.outputs[i]].z)) return -1-long(i);
    return long(mpz_sizeinbase(mod.z,2))-1;
}

int main(int argc,char** argv) {
    if(argc!=7) return 1;
    std::ifstream file(argv[1]);
    Input data;
    data.p=polynomial(file);
    for(size_t i=1;i<data.p.size();++i) {
        data.dp.emplace_back();
        mpz_mul_ui(data.dp.back().z,data.p[i].z,i);
    }
    unsigned long start=std::stoul(argv[2]),cutoff=std::stoul(argv[3]),count=std::stoul(argv[4]);
    int threads=std::stoi(argv[5]);
    long long bound=std::stoll(argv[6]);
    if(!count || count>cutoff || start>2000000000UL || cutoff>1000000 || threads<1 || threads>6) return 2;
    Integer forbidden;
    mpz_set_ui(data.common.z,1);
    mpz_primorial_ui(forbidden.z,cutoff);
    mpz_mul(forbidden.z,forbidden.z,data.p.back().z);
    for(int j=0;j<5;++j) {
        data.denom.push_back(read_integer(file));
        mpz_lcm(data.common.z,data.common.z,data.denom.back().z);
        mpz_mul(forbidden.z,forbidden.z,data.denom.back().z);
        data.numerators.push_back(polynomial(file));
    }
    for(const auto& v:data.denom) {
        data.scaling.emplace_back();
        mpz_divexact(data.scaling.back().z,data.common.z,v.z);
    }
    size_t n;
    if(!(file>>n) || n>10000) return 3;
    for(size_t i=0;i<n;++i) {
        Gate t;
        if(!(file>>t.op>>t.a>>t.b) || t.op<0 || t.op>4) return 3;
        if(t.op==0 && (t.a<0 || t.a>=6)) return 3;
        if(t.op>=2 && (t.a<0 || size_t(t.a)>=i)) return 3;
        if(t.op>=3 && (t.b<0 || size_t(t.b)>=i)) return 3;
        data.circuit.push_back(t);
    }
    for(int j=0;j<6;++j) {
        int o;
        if(!(file>>o) || o<0 || size_t(o)>=n) return 3;
        data.outputs.push_back(o);
    }
    long long total=0;
    for(unsigned long first=0;first<count;first+=64) {
        unsigned long batch=std::min(64UL,count-first);
        std::vector<long> bits(batch);
        #pragma omp parallel for num_threads(threads) schedule(dynamic)
        for(int i=0;i<int(batch);++i) bits[i]=verify(data,start+first+i,forbidden);
        for(unsigned long i=0;i<batch;++i) {
            if(bits[i]<0) {
                std::cerr<<"failed "<<start+first+i<<" "<<bits[i]<<std::endl;
                return 4;
            }
            total+=bits[i];
            std::cout<<bits[i]<<'\n';
        }
        std::cout.flush();
        if(total>bound) return 0;
    }
    return 0;
}
