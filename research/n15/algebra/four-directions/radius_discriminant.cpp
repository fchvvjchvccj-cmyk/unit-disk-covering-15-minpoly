// 对给定的双变量模多项式，在覆盖次数上界的完整单位根网格上消去方向。
#include <flint/nmod_poly.h>
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using U=uint32_t;
constexpr U prime=998244353;
namespace fs=std::filesystem;
U add(U a,U b) { U c=a+b; return c>=prime ? c-prime:c; }
U sub(U a,U b) { return a>=b ? a-b:a+prime-b; }
U mul(U a,U b) { return uint64_t(a)*b%prime; }
U power(U a,U n) {
    U b=1;
    for(;n;n>>=1) { if(n&1) b=mul(a,b); a=mul(a,a); }
    return b;
}
std::vector<U> reversal(U n) {
    std::vector<U> v(n);
    for(U i=1;i<n;++i) v[i]=v[i/2]/2+(i%2)*(n/2);
    return v;
}
void butterflies(U* v,U n,const std::vector<U>& roots,bool inverse=false) {
    U grid=roots.size();
    for(U length=2;length<=n;length*=2) {
        for(U i=0;i<n;i+=length) for(U j=0;j<length/2;++j) {
            U k=j*(grid/length);
            U a=v[i+j],b=mul(v[i+j+length/2],roots[inverse ? (-k)&(grid-1):k]);
            v[i+j]=add(a,b); v[i+j+length/2]=sub(a,b);
        }
    }
}

int main(int argc,char** argv) {
    if(argc!=7) return 1;
    U width=std::stoul(argv[2]), n=std::stoul(argv[3]);
    int threads=std::stoi(argv[6]);
    if(width<2 || n<1 || n>=prime || threads<1 || threads>6) return 2;
    uint64_t bound=uint64_t(2*n-1)*(width-1);
    if(bound>=(1u<<23)) return 3;
    U grid=1,block=1;
    while(grid<=bound) grid*=2;
    while(block<width) block*=2;
    U groups=grid/block;
    std::vector<U> coeffs(size_t(n+1)*width),roots(grid);
    std::ifstream matrix(argv[1],std::ios::binary);
    if(!matrix.read(reinterpret_cast<char*>(coeffs.data()),4*coeffs.size()) || matrix.peek()!=EOF) return 4;
    U omega=power(3,(prime-1)/grid);
    roots[0]=1;
    for(U i=1;i<grid;++i) roots[i]=mul(roots[i-1],omega);
    auto rev=reversal(block);
    fs::path directory(argv[4]);
    fs::create_directories(directory);
    int complete=0;
    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for(int group=0;group<int(groups);++group) {
        auto path=directory/(std::to_string(group)+".bin");
        if(!fs::exists(path)) {
            std::vector<U> table(size_t(n+1)*block), result(block);
            std::vector<U> scales(width);
            for(U i=0;i<width;++i) scales[i]=roots[(uint64_t(i)*group)&(grid-1)];
            for(U j=0;j<=n;++j) {
                U* buf=table.data()+size_t(j)*block;
                for(U i=0;i<width;++i) buf[rev[i]]=mul(coeffs[size_t(j)*width+i],scales[i]);
                butterflies(buf,block,roots);
            }
            nmod_poly_t f,df;
            nmod_poly_init2(f,prime,n+1); nmod_poly_init2(df,prime,n+1);
            for(U k=0;k<block;++k) {
                if(!table[size_t(n)*block+k]) continue;
                for(U j=0;j<=n;++j) f->coeffs[j]=table[size_t(j)*block+k];
                _nmod_poly_set_length(f,n+1);
                nmod_poly_derivative(df,f);
                result[k]=nmod_poly_resultant(f,df);
            }
            nmod_poly_clear(f); nmod_poly_clear(df);
            auto temp=path; temp+=".partial";
            std::ofstream out(temp,std::ios::binary);
            out.write(reinterpret_cast<char*>(result.data()),4*block);
            out.close();
            if(!out) std::terminate();
            fs::rename(temp,path);
        }
        #pragma omp critical
        {
            if(++complete%128==0) std::cout<<"blocks "<<complete<<"/"<<groups<<std::endl;
        }
    }
    std::vector<U> samples(grid),tmp(block);
    auto reverse_all=reversal(grid);
    for(U j=0;j<groups;++j) {
        std::ifstream in(directory/(std::to_string(j)+".bin"),std::ios::binary);
        if(!in.read(reinterpret_cast<char*>(tmp.data()),4*block) || in.peek()!=EOF) return 5;
        for(U k=0;k<block;++k) samples[reverse_all[j+k*groups]]=tmp[k];
    }
    butterflies(samples.data(),grid,roots,true);
    U inverse=power(grid,prime-2);
    for(auto& v:samples) v=mul(v,inverse);
    while(!samples.empty() && !samples.back()) samples.pop_back();
    if(samples.empty() || samples.size()>bound+1) return 6;
    std::ofstream out(argv[5],std::ios::binary);
    out.write(reinterpret_cast<char*>(samples.data()),4*samples.size());
    std::cout<<"degree "<<samples.size()-1<<std::endl;
    return out ? 0:7;
}
