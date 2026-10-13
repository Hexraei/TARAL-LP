#pragma once
#include "taral.hpp"
#include <array>
#include <cmath>
#include <cstring>
#include <chrono>
// No dynamic bigint dependency. 4352-bit fixed signed magnitude, base exponent -2148.
// Finite binary64 products highest exponent<=2047; <=500000 terms add<=19bits.
namespace exact_dyadic {
__extension__ typedef unsigned __int128 u128; // GCC/Clang extension; __extension__ keeps -Wpedantic quiet. Same type as before.
struct Wide{std::array<uint64_t,68>a{};bool negative=false;};
struct Part{uint64_t mant=0;int exponent=0;bool negative=false;};
inline Part decode(double x){uint64_t b;std::memcpy(&b,&x,8);unsigned e=(b>>52)&2047;Part p;p.mant=b&((1ULL<<52)-1);p.negative=b>>63;if(e){p.mant|=1ULL<<52;p.exponent=int(e)-1075;}else p.exponent=-1074;return p;}
inline bool zero(const Wide&w){for(auto v:w.a)if(v)return false;return true;}
inline int magnitude(const Wide&a,const Wide&b){for(int i=67;i>=0;--i)if(a.a[i]!=b.a[i])return a.a[i]<b.a[i]?-1:1;return 0;}
inline bool add(Wide&a,const Wide&b){if(a.negative==b.negative){u128 carry=0;for(int i=0;i<68;++i){carry+=(u128)a.a[i]+b.a[i];a.a[i]=uint64_t(carry);carry>>=64;}if(carry)return false;}else{int cmp=magnitude(a,b);if(!cmp){a=Wide{};return true;}const Wide&big=cmp>0?a:b;const Wide&small=cmp>0?b:a;Wide out;out.negative=big.negative;uint64_t borrow=0;for(int i=0;i<68;++i){uint64_t v=big.a[i]-small.a[i],b1=big.a[i]<small.a[i],r=v-borrow,b2=v<borrow;out.a[i]=r;borrow=b1|b2;}if(borrow)return false;a=out;}if(zero(a))a.negative=false;return true;}
inline bool term(Wide&out,double x,double y){if(!std::isfinite(x)||!std::isfinite(y))return false;Part a=decode(x),b=decode(y);u128 value=(u128)a.mant*b.mant;int shift=a.exponent+b.exponent+2148;if(shift<0||shift+106>4352)return false;out=Wide{};out.negative=a.negative!=b.negative;int limb=shift/64,offset=shift%64;uint64_t lo=uint64_t(value),hi=uint64_t(value>>64);out.a[limb]|=lo<<offset;if(offset){out.a[limb+1]|=lo>>(64-offset);out.a[limb+1]|=hi<<offset;if(hi>>(64-offset))out.a[limb+2]|=hi>>(64-offset);}else out.a[limb+1]=hi;if(zero(out))out.negative=false;return true;}
inline int compare(const Wide&a,const Wide&b){if(a.negative!=b.negative)return a.negative?-1:1;int c=magnitude(a,b);return a.negative?-c:c;}
inline bool within(const Wide&activity,double lo,double up){Wide bound;if(std::isfinite(lo)){if(!term(bound,lo,1.)||compare(activity,bound)<0)return false;}else if(lo!= -kInf)return false;if(std::isfinite(up)){if(!term(bound,up,1.)||compare(activity,bound)>0)return false;}else if(up!=kInf)return false;return true;}
inline bool feasible(const Model&md,const std::vector<double>&x,std::chrono::steady_clock::time_point deadline){if(std::chrono::steady_clock::now()>=deadline)return false;size_t n=md.cols.size(),m=md.row_lo.size(),nnz=0;if(n>100000||m>100000||x.size()!=n||md.is_int.size()!=n||md.col_lo.size()!=n||md.col_up.size()!=n||md.row_up.size()!=m)return false;size_t scanned=0;for(auto&col:md.cols){if((scanned++%256)==0&&std::chrono::steady_clock::now()>=deadline)return false;if(col.size()>500000-nnz)return false;nnz+=col.size();}if(std::chrono::steady_clock::now()>=deadline)return false;std::vector<Wide> activity(m);if(std::chrono::steady_clock::now()>=deadline)return false;Wide w;size_t visited=0;for(size_t j=0;j<n;++j){if((j%256)==0&&std::chrono::steady_clock::now()>=deadline)return false;if(!std::isfinite(x[j])||(md.is_int[j]&&x[j]!=std::trunc(x[j]))||!term(w,x[j],1.)||!within(w,md.col_lo[j],md.col_up[j]))return false;for(auto&e:md.cols[j]){if(e.index<0||size_t(e.index)>=m||!term(w,e.value,x[j])||!add(activity[e.index],w))return false;if((++visited%1024)==0&&std::chrono::steady_clock::now()>=deadline)return false;}}for(size_t i=0;i<m;++i){if(i%256==0&&std::chrono::steady_clock::now()>=deadline)return false;if(!within(activity[i],md.row_lo[i],md.row_up[i]))return false;}return std::chrono::steady_clock::now()<deadline;}
}
